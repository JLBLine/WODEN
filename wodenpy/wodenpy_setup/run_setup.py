import wodenpy
from wodenpy.wodenpy_setup.git_helper import retrieve_gitdict
import numpy as np
import sys
import os
from astropy.io import fits

def get_parser():
    """
    Runs the argument parser to get command line inputs - used by sphinx and
    argparse extension to unpack the help below into the online readthedocs
    documentation.

    Returns
    -------
    parser : `argparse.ArgumentParser`
        The populated argument parser used by `run_woden.py`

    """
    import argparse
    from argparse import RawTextHelpFormatter

    class SmartFormatter(argparse.HelpFormatter):
        """Argparse by default ignores all \n and \t formatters. If you start
        a help class with R| the formatters will be respected."""
        def _split_lines(self, text, width):
            if text.startswith('R|'):
                return text[2:].splitlines()
            # this is the RawTextHelpFormatter._split_lines
            return argparse.HelpFormatter._split_lines(self, text, width)

    parser = argparse.ArgumentParser(description="Run the WODEN simulator and profit. "
                "WODEN is setup to simulate MWA-style observations, where the "
                "full frequency bandwidth is split into 24 'coarse' bands, each "
                "of which is split into fine channels. This naturally allows "
                "any simulation to be split across multiple GPUs as separate "
                "processes.", formatter_class=SmartFormatter)

    freq_group = parser.add_argument_group('FREQUENCY OPTIONS')
    freq_group.add_argument('--band_nums', default='all',
        help='Defaults to running 24 coarse bands. Alternatively, enter required'
             ' numbers delineated by commas, e.g. --band_nums=1,7,9')
    freq_group.add_argument('--lowest_channel_freq', default=False,
        help='Set the frequency (Hz) of the lowest channel for band 1. '
             'If using a metafits file, this will override the frequency in'
             ' the metafits')
    freq_group.add_argument('--coarse_band_width', type=float, default=1.28e+6,
        help='Set the width of each coarse band \
              If using a metafits file, this will override the frequency in '
              'the metafits')
    freq_group.add_argument('--num_freq_channels', default='obs',
        help='Number of fine frequency channels to simulate - defaults to '
             '--coarse_band_width / --freq_res')
    freq_group.add_argument('--freq_res', type=float, default=False,
        help='Fine channel frequnecy resolution (Hz) - will default to what'
             ' is in the metafits')

    time_group = parser.add_argument_group('TIME OPTIONS')
    time_group.add_argument('--num_time_steps', default=False,
        help='The number of time steps to simualte. Defaults to how many are in'
             'the metafits if using metafits')
    time_group.add_argument('--time_res', type=float,default=False,
        help='Time resolution (s) - will default to what is in the metafits '
              'if the metafits if using metafits')


    obs_group = parser.add_argument_group('OBSERVATION OPTIONS')
    obs_group.add_argument('--ra0', type=float, required=True,
        help='RA of the desired phase centre (deg)')
    obs_group.add_argument('--dec0', type=float, required=True,
        help='Dec of the desired phase centre (deg)')
    obs_group.add_argument('--date', default=False,
        help='Initial UTC date of the observatio in format YYYY-MM-DDThh:mm:ss '
             'This is used to set the LST and array precession. This is set '
             'automatically when reading a metafits but including this will '
             'override the date in the metafits')
    obs_group.add_argument('--no_precession', default=False, action='store_true',
        help='By default, WODEN rotates the array back to J2000 to match '
             'the input sky catalogue. Add this to switch off precession')

    tel_group = parser.add_argument_group('TELESCOPE OPTIONS')
    tel_group.add_argument('--latitude', default=-26.703319405555554, type=float,
        help='Latitude (deg) of the array - defaults to MWA at -26.703319405555554')
    tel_group.add_argument('--longitude', default=116.67081523611111, type=float,
        help='Longitude (deg) of the array - defaults to MWA at 116.67081523611111')
    tel_group.add_argument('--array_height', default=377.827, type=float,
        help='Height (m) of the array above sea level - defaults to MWA at 377.827')
    tel_group.add_argument('--array_layout', default=False,
        help='Instead of reading the array layout from the metafits file, read'
             ' from a text file. Store antenna positions as offset from array '
             'centre, in east, north, height coords (metres)')
    tel_group.add_argument('--primary_beam', default="none",
        help="R|Which primary beam to use in the simulation.\nOptions are:\n"
            "\t - MWA_FEE (MWA fully embedded element model)\n"
            "\t - MWA_FEE_interp (MWA fully embedded element model that has had)\n"
            "\t\t spherical harmonics interpolated over frequency\n"
            "\t - Gaussian (Analytic symmetric Gaussian)\n"
            "\t\t see --gauss_beam_FWHM and\n"
            "\t\t and --gauss_beam_ref_freq for\nfine control)\n"
            "\t - EDA2 (Analytic dipole with a ground mesh) \n"
            "\t - MWA_analy (MWA analytic model)\n"
            "\t - none (Don't use a primary beam at all)\n"
            "Defaults to --primary_beam=none")


    tel_group.add_argument('--gauss_beam_FWHM', default=False,
        help='The FWHM of the Gaussian beam in deg - WODEN defaults to using'
             ' 20 deg if this is not set')
    tel_group.add_argument('--gauss_beam_ref_freq', default=False,
        help='The frequency at which the gauss beam FWHM is set at. If not set,'
             ' WODEN will default to 150MHz.')
    tel_group.add_argument('--gauss_ra_point', default=False,
        help='The initial RA (deg) to point the Gaussian beam at. This will be '
              'used to calculate an hour angle at which the beam will remain '
              'pointed at for the duration of the observation. Defaults to the '
              'RA of the metafits if available, or the RA of the phase centre '
              'if not')
    tel_group.add_argument('--gauss_dec_point', default=False,
        help='The initial Dec (deg) to point the Gaussian beam at. Defaults '
        'to the Dec of the metafits if available, or the Dec of the phase centre'
        ' if not')

    tel_group.add_argument('--hdf5_beam_path', default=False,
        help='Location of the hdf5 file holding the FEE beam coefficients')
    tel_group.add_argument('--MWA_FEE_delays', default=False,
        help='R|A list of 16 delays to point the MWA FEE primary beam \n'
              'model enter as as list like: \n'
              '--MWA_FEE_delays=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]\n'
              'for a zenith pointing. This is read directly from\n'
              'the metafits if using a metafits file')
    tel_group.add_argument('--telescope_name', default='MWA',
        help='Name of telescope written out to the uvfits file, defaults to MWA')

    input_group = parser.add_argument_group('INPUT/OUTPUT OPTIONS')
    input_group.add_argument('--IAU_order', default=False,
        help='NEW IN WODEN versions >= 1.4.0. By default, the XX pol is now '
             'from the East-West aligned dipoles. Add --IAU_order to define '
             'the XX pol as North-South. Background: '
             'the first polaristaion output in the uvfits is '
             'called "XX". The IAU defines "XX" as aligned to North-South, '
             'however typically it is assumed that "XX" in a uvfits is from '
             'the East-West dipoles. Le sigh. For WODEN versions < 1.4.0, '
             'XX was always the N-S dipoles. ')
    input_group.add_argument('--cat_filename', required=True,
        help='Path to WODEN style sky model')
    input_group.add_argument('--metafits_filename',default=False,
        help='MWA style metafits file to base the simulation on. Array layout,'
             ' frequency and time parameters are all set by this option, but '
             'can be overridden using other arguments')
    input_group.add_argument('--output_uvfits_prepend',default='output',
        help='Prepend name for uvfits - will append band%%02d.uvfits %%band_num '
             'at the end. Defaults to "output".')
    input_group.add_argument('--sky_crop_components', default=True, action='store_true',
        help='This option is deprecated, but held here for compatibility.')
    input_group.add_argument('--sky_crop_sources', default=False, action='store_true',
        help='WODEN will crop out sky model information that is below the '
             'horizon for the given LST. By default, '
             'WODEN will include any COMPONENT above the horizon, regardless '
             'of which SOURCE it belongs to. If --sky_crop_source is included '
             'for each SOURCE in the sky model, if any COMPONENT is below the ' 'horizon, the entire source will be flagged')
    input_group.add_argument('--do_autos', default=False, action='store_true',
        help='By default, WODEN only calculates cross-correlations. Add this'
             ' flag to include auto-correlations.')

    sim_group = parser.add_argument_group('SIMULATOR OPTIONS')
    sim_group.add_argument('--precision', default='double',
        help='What precision to run WODEN at. Options are "double" or "float". '
             'Defaults to "double"')
    sim_group.add_argument('--remove_phase_tracking', default=False, action='store_true',
        help='By adding this flag, remove the phase tracking of the '
             'visibilities - use this to feed uvfits into the RTS')
    sim_group.add_argument('--no_tidy', default=False, action='store_true',
        help='Defaults to deleting output binary files from woden and json '
             'files. Add this flag to not delete those files')
    sim_group.add_argument('--chunking_size', type=float, default=False,
        help='The chunk size to break up the point sources into for processing '
             '- defaults to 0 (use default chunking in WODEN)')
    sim_group.add_argument('--dry_run', default=False, action='store_true',
        help='Add this to NOT call the WODEN executable - this will just write '
             'out the .json file and do nothing else')


    ##Add a number of hidden arguments. This means we can add attributes to
    ##the args object to conveniently pass things into functions, but without
    ##them showing up in --help
    parser.add_argument('--east', help=argparse.SUPPRESS)
    parser.add_argument('--north', help=argparse.SUPPRESS)
    parser.add_argument('--height', help=argparse.SUPPRESS)
    parser.add_argument('--num_antennas', help=argparse.SUPPRESS)
    parser.add_argument('--array_layout_name', help=argparse.SUPPRESS)

    return parser

def select_argument_and_check(parser_arg, parser_value,
                              metafits_arg, parser_string,
                              do_exit=True):
    """Some arguments taken from the argparse.parser should override settings
    from the metafits if present. If the parser argument `parser_arg` is
    defined (i.e. not False), update it to equal `parser_value`. If not defined,
    update `parser_arg` to `metafits_arg`, which is the value read in from
    the metafits file. If both `parser_arg` and `metafits_arg` are False,
    WODEN will fail, so exit with a message. Use `parser_string` to define
    which parser arguement has failed; this will be included in the error
    message.

    Parameters
    ----------
    parser_arg : attribute of `argparse.Namespace`
        The option in `args` to update
    parser_value : Expected type for `parser_arg`
        The value to set `parser_arg` to (e.g. float(parser_arg))
    metafits_arg : Expected type for `parser_arg`
        The value read in from the metafits if using metafits; False if not
    parser_string : string
        The parser option under test to be written out in the error message,
        e.g. "--MWA_FEE_delays"
    do_exit : Boolean
        Whether to call `sys.exit` upon both `parser_arg` and `metafits_arg`
        being False. Defaults to True

    Returns
    -------
    parser_arg : attribute of `argparse.Namespace`
        The update option in `args`
    """

    ##If a parser arg is there, reset it to the parser_value give
    if parser_arg:
        parser_arg = parser_value
    ##If not there, check if a metafits equivalent has been found
    else:
        if metafits_arg:
            parser_arg = metafits_arg
        else:

            error_message = ("ARGS ERROR: args.{:s} has not been set. \n"
            "Either specify using --{:s} or get from a metafits using "
            "--metafits_filename\nExiting now as WODEN cannot run").format(parser_string, parser_string)
            if do_exit:
                sys.exit(error_message)

    return parser_arg

def select_correct_enh(args):
    """Depending on whether we are reading the array layout from the metafits
    file or a text file, read in the correct amount of east,north,height coords.
    Sets `args.east`, `args.north`, `args.height`, `args.num_antennas`, and
    `args.array_layout_name`.

    Parameters
    ----------
    args : `argparse.Namespace`
        The populated arguments `args = parser.parse_args()`` as returned from
        the parser given by :func:`~run_woden.get_parser`
    """

    if args.array_layout == "from_the_metafits":
        ##Using metafits for array layout. Have previously read in e,n,h
        ##In the metafits it lists XX,YY for each antenna so we select every second one
        selection = np.arange(0,len(args.east),2)
        args.num_antennas = int(len(selection))

        args.east = args.east[selection]
        args.north = args.north[selection]
        args.height = args.height[selection]
        args.ant_names = args.ant_names[selection]

        array_layout = np.zeros((args.num_antennas,3))

        array_layout[:,0] = args.east
        array_layout[:,1] = args.north
        array_layout[:,2] = args.height

        args.array_layout_name = 'WODEN_array_layout.txt'

        np.savetxt(args.array_layout_name, array_layout)
    else:
        try:
            array_layout = np.loadtxt(args.array_layout)
            args.num_antennas,_ = array_layout.shape

            args.east = array_layout[:,0]
            args.north = array_layout[:,1]
            args.height = array_layout[:,2]
            
            args.ant_names = np.array(["%05d" %ant for ant in range(1,args.num_antennas + 1)])

        except:
            exit("Could not read array layout file:\n"
                 "\t{:s}\nExiting before woe beings".format(args.array_layout))

        args.array_layout_name = args.array_layout

def check_args(args):
    """Check that the combination of arguments parsed will work with the
    WODEN executable. Attempts to grab information from a metafits file if
    possible. Should error with helpful messages if a combination that won't
    work is attempted by the user

    Parameters
    ----------
    args : `argparse.Namespace`
        The populated arguments `args = parser.parse_args()`` as returned from
        the parser given by :func:`~run_woden.get_parser`

    Returns
    -------
    args : `argparse.Namespacer`
        The populated arguments which will now have been checked and had
        information from metafits incorporated if requested
    """

    if args.primary_beam not in ['MWA_FEE', 'Gaussian', 'EDA2', 'none', 'None',
                                 'MWA_FEE_interp', 'MWA_analy']:
        exit('Primary beam option --primary_beam must be one of:\n'
             '\t MWA_FEE, MWA_FEE_interp, Gaussian, EDA2, none\n'
             'User has entered --primary_beam={:s}\n'
             'Please fix and try again. Exiting now'.format(args.primary_beam))

    ##Be a little flexible in how people specify 'none'
    if args.primary_beam in ['None', 'none']:
        args.primary_beam = 'none'

    ##If we're using the MWA FEE beam, make sure we can find the stored
    ##spherical harmonics file
    if args.primary_beam == 'MWA_FEE':
        if args.hdf5_beam_path:
            if not os.path.isfile(args.hdf5_beam_path):
                exit('Could not open hdf5 MWA FEE path as specified by user as:\n'
                     '\t--hdf5_beam_path={:s}.\n'
                     'This will cause WODEN to fail, exiting now'.format(args.hdf5_beam_path))
        else:
            try:
                MWA_FEE_HDF5 = os.environ['MWA_FEE_HDF5']
                args.hdf5_beam_path = MWA_FEE_HDF5
                if not os.path.isfile(args.hdf5_beam_path):
                    exit('Could not open hdf5 MWA FEE path as specified by user as:\n'
                         '\t--environ["MWA_FEE_HDF5"]={:s}.\n'
                         'This will cause WODEN to fail, exiting now'.format(args.hdf5_beam_path))
            except KeyError:
                exit('To use MWA FEE beam, either --hdf5_beam_path or environment\n'
                     'variable MWA_FEE_HDF5 must point towards the file\n'
                     'mwa_full_embedded_element_pattern.h5. Exiting now as WODEN will fail.')

    ##If we're using the MWA FEE beam, make sure we can find the stored
    ##spherical harmonics file
    elif args.primary_beam == 'MWA_FEE_interp':
        if args.hdf5_beam_path:
            if not os.path.isfile(args.hdf5_beam_path):
                exit('Could not open hdf5 MWA FEE path as specified by user as:\n'
                     '\t--hdf5_beam_path={:s}.\n'
                     'This will cause WODEN to fail, exiting now'.format(args.hdf5_beam_path))
        else:
            try:
                MWA_FEE_HDF5_INTERP = os.environ['MWA_FEE_HDF5_INTERP']
                args.hdf5_beam_path = MWA_FEE_HDF5_INTERP
                if not os.path.isfile(args.hdf5_beam_path):
                    exit('Could not open hdf5 MWA FEE path as specified by user as:\n'
                         '\t--environ["MWA_FEE_HDF5_INTERP"]={:s}.\n'
                         'This will cause WODEN to fail, exiting now'.format(args.hdf5_beam_path))
            except KeyError:
                exit('To use MWA FEE intrep beam, either --hdf5_beam_path or environment\n'
                     'variable MWA_FEE_HDF5_INTERP must point towards the file\n'
                     'MWA_embedded_element_pattern_rev2_interp_167_197MHz.h5. Exiting now as WODEN will fail.')

    ##variables that will be filled by metafits if reading a metafits
    ##set them as False here for testing later on
    MWA_FEE_delays = False
    time_res = False
    freq_res = False
    freqcent = False
    lowest_channel_freq = False
    num_time_steps = False
    date = False
    array_layout = False

    ##read in args from the metafits if requested
    if args.metafits_filename:

        if not os.path.isfile(args.metafits_filename):
            exit('Could not open metafits specified by user as:\n'
                 '\t--metafits_filename={:s}.\n'
                 'Cannot get required observation settings, exiting now'.format(args.metafits_filename))

        with fits.open(args.metafits_filename) as f:
            date = f[0].header['DATE-OBS']

            ##Get the east, north, height antenna positions from the metafits
            ##Need to order the antennas via the Tile
            ##column to be consistent with hyperdrive
            antenna_order = np.argsort(f[1].data['Tile'])

            east = f[1].data['East'][antenna_order]
            north = f[1].data['North'][antenna_order]
            height = f[1].data['Height'][antenna_order]
            tilenames = f[1].data['Tilename'][antenna_order]
            
            args.east = east
            args.north = north
            args.height = height
            args.ant_names = tilenames

            ##Use this to signal that reading in array layout from metafits
            ##was successful
            array_layout = "from_the_metafits"

            ##Read observation parameters from the metafits file
            time_res = float(f[0].header['INTTIME'])
            freq_res = float(f[0].header['FINECHAN'])*1e+3
            freqcent = float(f[0].header['FREQCENT'])*1e+6
            b_width = float(f[0].header['BANDWDTH'])*1e+6
            lowest_channel_freq = freqcent - (b_width/2) - (freq_res/2)

            num_time_steps = int(f[0].header['NSCANS'])

            delays = np.array(f[0].header['DELAYS'].split(','),dtype=int)
            delays[np.where(delays == 32)] = 0
            MWA_FEE_delays = str(list(delays))

            ##If user hasn't specified a pointing for a Gaussian beam,
            ##fill in using the metafits file
            if not args.gauss_ra_point:
                args.gauss_ra_point = float(f[0].header['RA'])
            if not args.gauss_dec_point:
                args.gauss_dec_point = float(f[0].header['DEC'])

            f.close()

    ##Override metafits and/or load arguements
    args.lowest_channel_freq = select_argument_and_check(args.lowest_channel_freq,
                                  float(args.lowest_channel_freq),
                                  lowest_channel_freq, "lowest_channel_freq")

    args.num_time_steps = select_argument_and_check(args.num_time_steps,
                                  int(args.num_time_steps),
                                  num_time_steps, "num_time_steps")

    args.freq_res = select_argument_and_check(args.freq_res, args.freq_res,
                                  freq_res, "freq_res")

    args.time_res = select_argument_and_check(args.time_res, args.time_res,
                                  time_res, "time_res")

    args.date = select_argument_and_check(args.date, args.date,
                                  date, "date")

    args.array_layout = select_argument_and_check(args.array_layout, args.array_layout,
                                  array_layout, "array_layout")

    ##TODO change this from MWA_FEE_delays to MWA_delays (or allow both via
    ##some argparse magic)
    ##If the user has manually specified some MWA FEE delays, ensure they
    ##can be made into an array of 16 floats
    if args.MWA_FEE_delays:
        message = ("ERROR - failed to convert --MWA_FEE_delays into a list"
                   " of 16 floats correctly. You have entered:\n"
                   "    --MWA_FEE_delays={:s}\n"
                   "Exiting now.".format(args.MWA_FEE_delays))
        try:
            test_list = list(np.array(args.MWA_FEE_delays.strip('[]').split(','),dtype=float))
            if len(test_list) != 16:
                exit(message)
        except:
            exit(message)

    ##Do the test on MWA_FEE_delays only if this is an MWA_FEE simulation
    if args.primary_beam == 'MWA_FEE' or args.primary_beam == 'MWA_FEE_interp' or args.primary_beam == 'MWA_analy':
        args.MWA_FEE_delays = select_argument_and_check(args.MWA_FEE_delays,
                                      args.MWA_FEE_delays,
                                      MWA_FEE_delays, "MWA_FEE_delays")

    ##Set the band numbers we are simulating in this run
    if args.num_freq_channels == 'obs':
        args.num_freq_channels = int(np.floor(args.coarse_band_width / args.freq_res))
    else:
        args.num_freq_channels = int(args.num_freq_channels)

    if args.band_nums == 'all':
        args.band_nums = range(1,25)
    else:
        try:
            args.band_nums = list(np.array(args.band_nums.split(','),dtype=int))
        except:
            message = ("ERROR - failed to convert --band_nums into a list of ints"
                       " correctly. You entered:\n"
                       "    --band_nums={:s}\n"
                       "Exiting now.")
            exit(message)

    ##If pointing for Gaussian beam is not set, point it at the phase centre
    if args.primary_beam == 'Gaussian':
        ##Explicitly test if False here as it could point at 0.0 deg, which
        ##get interpreted as False in a simple if statement doh
        if args.gauss_ra_point is False:
            args.gauss_ra_point = args.ra0
        else:
            args.gauss_ra_point = float(args.gauss_ra_point)

        if args.gauss_dec_point is False:
            args.gauss_dec_point = args.dec0
        else:
            args.gauss_dec_point = float(args.gauss_dec_point)

        if args.gauss_beam_ref_freq: args.gauss_beam_ref_freq = float(args.gauss_beam_ref_freq)
        if args.gauss_beam_FWHM: args.gauss_beam_FWHM = float(args.gauss_beam_FWHM)

    if args.precision not in ['double', 'float']:
        print(f"Arg --precision={args.precision} is not valid. Should be either"
              "'double' or 'float'. Setting to 'double'")
        args.precision='double'

    ##Either read the array layout from a file or use what was in the metafits
    select_correct_enh(args)

    ##If user asks to crop by source, change cropping by component to False
    if args.sky_crop_sources:
        args.sky_crop_components = False

    return args

def get_code_version():
    """
    Returns either the git hash if installed via a git repo, or the __version__
    if installed from a release

    Returns
    -------
    version : string
        Either the git commit or release version
    """

    git_dict = retrieve_gitdict()
    if git_dict:
        version = git_dict['describe']

    else:
        version = wodenpy.__version__
    
    print(f"You are using WODEN commit: {version}")

    return version
