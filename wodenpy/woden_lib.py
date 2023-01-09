# from __future__ import print_function
# from astropy.io import fits
import numpy as np
# from struct import unpack
import subprocess
# import os

##version of release, a fall back if this isn't in a git repo
VERSION = "2.0.0"

def command(cmd):
    """
    Runs the command string `cmd` using `subprocess.call`

    Parameters
    ----------
    cmd : string
         The command to run on the command line
    """
    subprocess.call(cmd,shell=True)
    # print(cmd)

def write_json(json_name=None, jd_date=None, lst=None, args=None):
    """
    Populate and write out a .json parameter file used to run WODEN.
    Is later used on the command line to run WODEN.

    Parameters
    ----------
    json_name : string
        Name out .json file to save outputs to
    jd_date : float
        Initial Julian date of simulation (days)
    lst : float
        Local sidereal time of the simulate (degrees)
    args : `argparse.Namespace`
        The args as returned by :func:`~run_woden.check_args`, which takes
        the args as return by `args=parser.parse_args()`, from the `parser`
        returned by :func:`~run_woden.get_parser`.

    """

    with open(json_name,'w+') as outfile:

        outfile.write('{\n')
        outfile.write('  "ra0": {:.16f},\n'.format(args.ra0))
        outfile.write('  "dec0": {:.16f},\n'.format(args.dec0))
        outfile.write('  "num_freqs": {:d},\n'.format(args.num_freq_channels))
        outfile.write('  "num_time_steps": {:d},\n'.format(args.num_time_steps))
        outfile.write('  "cat_filename": "{:s}",\n'.format(args.cat_filename))
        outfile.write('  "time_res": {:.5f},\n'.format(args.time_res))
        outfile.write('  "frequency_resolution": {:.3f},\n'.format(args.freq_res))
        outfile.write('  "chunking_size": {:d},\n'.format(int(args.chunking_size)))
        outfile.write('  "jd_date": {:.16f},\n'.format(jd_date))
        outfile.write('  "LST": {:.16f},\n'.format(lst))

        ##sometimes some bands finish before others start on a super cluster,
        ##and delete this array layout if generated from metafits. So append
        ##a band number to prevent this from happening

        if args.array_layout == 'from_the_metafits':

            band_num = json_name.split('_')[-1].split('.')[0]
            band_array_layout = f"WODEN_array_layout_band{band_num}.txt"
            command(f"cp WODEN_array_layout.txt {band_array_layout}")
        else:
            band_array_layout = args.array_layout_name

        outfile.write('  "array_layout": "{:s}",\n'.format(band_array_layout))
        outfile.write('  "lowest_channel_freq": {:.10e},\n'.format(args.lowest_channel_freq))
        outfile.write('  "latitude": {:.16f},\n'.format(args.latitude))
        outfile.write('  "coarse_band_width": {:.10e},\n'.format(args.coarse_band_width))

        if args.sky_crop_components:
            outfile.write('  "sky_crop_components": "True",\n')

        if args.no_precession:
            outfile.write('  "no_precession": "True",\n')

        if args.primary_beam == 'Gaussian':
            outfile.write('  "use_gaussian_beam": "True",\n')
            if args.gauss_beam_FWHM:
                outfile.write('  "gauss_beam_FWHM": %.10f,\n' %float(args.gauss_beam_FWHM))

            if args.gauss_beam_ref_freq:
                outfile.write('  "gauss_beam_ref_freq": %.10f,\n' %float(args.gauss_beam_ref_freq))

            outfile.write('  "gauss_ra_point": %.8f,\n' %float(args.gauss_ra_point))
            outfile.write('  "gauss_dec_point": %.8f,\n' %float(args.gauss_dec_point))

        elif args.primary_beam == 'MWA_FEE':
            outfile.write('  "use_FEE_beam": "True",\n')
            outfile.write('  "hdf5_beam_path": "%s",\n' %args.hdf5_beam_path)
            outfile.write('  "FEE_delays": %s,\n ' %args.MWA_FEE_delays)

        elif args.primary_beam == 'MWA_FEE_interp':
            outfile.write('  "use_FEE_interp_beam": "True",\n')
            outfile.write('  "hdf5_beam_path": "%s",\n' %args.hdf5_beam_path)
            outfile.write('  "FEE_delays": %s,\n ' %args.MWA_FEE_delays)

        elif args.primary_beam == 'MWA_analy':
            outfile.write('  "use_MWA_analy_beam": "True",\n')
            outfile.write('  "FEE_delays": %s,\n ' %args.MWA_FEE_delays)

        elif args.primary_beam == 'EDA2':
            outfile.write('  "use_EDA2_beam": "True",\n')

        if args.do_autos:
            outfile.write('  "do_autos": "True",\n')

        if len(args.band_nums) == 1:
            band_str = '[%d]' %args.band_nums[0]
        else:

            band_str = '[%d' %args.band_nums[0]
            for band in args.band_nums[1:-1]:
                band_str += ',%d' %band
            band_str += ',%d]' %args.band_nums[-1]
        outfile.write('  "band_nums": %s\n' %band_str)
        outfile.write('}\n')