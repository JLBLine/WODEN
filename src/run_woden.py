#!/usr/bin/env python
'''Wrapper script to generate json input files for, and to run,
the GPU WODEN code.
'''
from __future__ import print_function
from astropy.io import fits
from astropy.time import Time
from numpy import *
from struct import unpack
from subprocess import call, check_output
# from jdcal import gcal2jd
from os import environ
import os
import warnings

R2D = 180.0 / pi
D2R = pi / 180.0
MWA_LAT = -26.7033194444
VELC = 299792458.0
SOLAR2SIDEREAL = 1.00274

WODEN_DIR = environ['WODEN_DIR']

def command(cmd):
    call(cmd,shell=True)
    # print(cmd)

def calc_jdcal(date):
    '''Takes a string format date and returns julian date in
    the two chunks a uvfits file likes'''

    t = Time(date)
    jd = t.jd

    jd_day = floor(jd)
    jd_fraction = (jd - jd_day)

    ##The header of the uvdata file takes the integer, and
    ##then the fraction goes into the data array for PTYPE5
    return jd_day, jd_fraction

def RTS_encode_baseline(b1, b2):
    '''The ancient aips/miriad extended way of encoding a baseline.
    Needed for populating the uvfits file'''
    if b2 > 255:
        return b1*2048 + b2 + 65536
    else:
        return b1*256 + b2

def RTS_decode_baseline(blcode):
    blcode = int(blcode)

    if blcode > 65536:
        blcode -= 65536
        b2 = int(blcode % 2048)
        b1 = int((blcode - b2) / 2048)
    else:
        b2 = int(blcode % 256)
        b1 = int((blcode - b2) / 256)

    return b1,b2


def make_antenna_table(XYZ_array=None,telescope_name=None,num_antennas=None,
                       freq_cent=None,date=None):
    """Write an antenna table for a uvfits file"""

    ##Make some values for certain columns
    annnames = array(["%05d" %ant for ant in range(1,num_antennas+1)])
    xlabels = array(['X']*num_antennas)
    ylabels = array(['Y']*num_antennas)

    ##Make a number of FITS columns to create the antenna table from
    col1 = fits.Column(array=annnames,name='ANNAME',format='5A')
    col2 = fits.Column(array=XYZ_array,name='STABXYZ',format='3D')
    ##col3 makes an empty array, and the format states skip reading this column
    ##Just replicating the example uvfits I've been using
    col3 = fits.Column(array=array([]),name='ORBPARM',format='0D')
    col4 = fits.Column(array=arange(1,num_antennas+1),name='NOSTA',format='1J')
    col5 = fits.Column(array=zeros(num_antennas),name='MNTSTA',format='1J')
    col6 = fits.Column(array=zeros(num_antennas),name='STAXOF',format='1E')
    col7 = fits.Column(array=xlabels,name='POLTYA',format='1A')
    col8 = fits.Column(array=zeros(num_antennas),name='POLAA',format='1E')
    col9 = fits.Column(array=zeros(num_antennas),name='POLCALA',format='1E')
    col10 = fits.Column(array=ylabels,name='POLTYB',format='1A')
    col11 = fits.Column(array=zeros(num_antennas),name='POLAB',format='1E')
    col12 = fits.Column(array=zeros(num_antennas),name='POLCALB',format='1E')

    ##Stick the columns into a ColDefs
    coldefs = fits.ColDefs([col1,col2,col3,col4,col5,col6, \
                            col7,col8,col9,col10,col11,col12])

    ##Use the columns to for a BinTableHDU object. This is shoved into the
    ##uvfits file later
    ##Astropy doesn't like the fact we have a zero sized column (col3 see above)
    ##so supress the warning when making the BinTableHDU
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        hdu_ant = fits.BinTableHDU.from_columns(coldefs, name="AIPS AN")

    ##-----Add some header values that seem to be needed by casa/RTS/WSClean

    ##Absolute reference point of the centre of the array - this is hardcoded to
    ##the MWA at the moment; X,Y,Z are relative to this
    ##TODO get the equation to work this out from lat/lon
    hdu_ant.header['ARRAYX']  = -2559453.2906
    hdu_ant.header['ARRAYY']  = 5095371.73544
    hdu_ant.header['ARRAYZ']  = -2849056.77357

    hdu_ant.header['FREQ']    = freq_cent
    hdu_ant.header['RDATE']   = date
    hdu_ant.header['TIMSYS']  = 'UTC     '
    hdu_ant.header['ARRNAM']  = telescope_name
    hdu_ant.header['NUMORB']  = 0
    hdu_ant.header['NOPCAL']  = 0
    hdu_ant.header['POLTYPE'] = '        '
    hdu_ant.header['CREATOR']   = 'WODEN_uvfits_writer'

    return hdu_ant

def create_uvfits(v_container=None,freq_cent=None, ra_point=None, dec_point=None,
                  output_uvfits_name=None,uu=None,vv=None,ww=None,
                  baselines_array=None, date_array=None, date=None,
                  central_freq_chan=None,ch_width=None,
                  int_jd=None, hdu_ant=None, gitlabel=False):
    '''Takes visibility data and writes out a uvfits format file'''
    ##TODO replace all of this with an interface with pyuvdata

    uvparnames = ['UU','VV','WW','BASELINE','DATE']
    parvals = [uu,vv,ww,baselines_array,date_array]

    uvhdu = fits.GroupData(v_container,parnames=uvparnames,pardata=parvals,bitpix=-32)
    uvhdu = fits.GroupsHDU(uvhdu)

    ###uvfits standards
    uvhdu.header['CTYPE2'] = 'COMPLEX '
    uvhdu.header['CRVAL2'] = 1.0
    uvhdu.header['CRPIX2'] = 1.0
    uvhdu.header['CDELT2'] = 1.0

    ##This means it's linearly polarised
    uvhdu.header['CTYPE3'] = 'STOKES '
    uvhdu.header['CRVAL3'] = -5.0
    uvhdu.header['CRPIX3'] =  1.0
    uvhdu.header['CDELT3'] = -1.0

    uvhdu.header['CTYPE4'] = 'FREQ'
    uvhdu.header['CRVAL4'] = freq_cent  ##Middle pixel value in Hz
    uvhdu.header['CRPIX4'] = int(central_freq_chan) + 1 ##Middle pixel number
    uvhdu.header['CDELT4'] = ch_width

    uvhdu.header['CTYPE5'] = 'RA'
    uvhdu.header['CRVAL5'] = ra_point
    uvhdu.header['CRPIX5'] = 1.0
    uvhdu.header['CDELT5'] = 1.0

    uvhdu.header['CTYPE6'] = 'DEC'
    uvhdu.header['CRVAL6'] = dec_point
    uvhdu.header['CRPIX6'] = 1.0
    uvhdu.header['CDELT6'] = 1.0

    ## Write the parameters scaling explictly because they are omitted if default 1/0
    uvhdu.header['PSCAL1'] = 1.0
    uvhdu.header['PZERO1'] = 0.0
    uvhdu.header['PSCAL2'] = 1.0
    uvhdu.header['PZERO2'] = 0.0
    uvhdu.header['PSCAL3'] = 1.0
    uvhdu.header['PZERO3'] = 0.0
    uvhdu.header['PSCAL4'] = 1.0
    uvhdu.header['PZERO4'] = 0.0
    uvhdu.header['PSCAL5'] = 1.0

    int_jd, float_jd = calc_jdcal(date)
    uvhdu.header['PZERO5'] = float(int_jd)

    ##Old observation parameters that were/are needed in CHIPS
    uvhdu.header['OBJECT']  = 'Undefined'
    uvhdu.header['OBSRA']   = ra_point
    uvhdu.header['OBSDEC']  = dec_point
    # uvhdu.header['TELESCOP'] = 'MWA'

    ##Add in the gitlabel so we know what version generated the file
    if gitlabel: uvhdu.header['GITLABEL'] = gitlabel

    ## Create hdulist and write out file
    hdulist = fits.HDUList(hdus=[uvhdu,hdu_ant])
    hdulist.writeto(output_uvfits_name,overwrite=True)
    hdulist.close()

def enh2xyz(east,north,height,latitude=MWA_LAT*D2R):
    '''Takes east, north, height coords out of the metafits file,
    and returns local X,Y,Z coords to put in the uvfits file'''
    sl = sin(latitude)
    cl = cos(latitude)
    X = -north*sl + height*cl
    Y = east
    Z = north*cl + height*sl
    return X,Y,Z

def load_data(filename=None,num_baselines=None,num_freq_channels=None,num_time_steps=None):
    '''Read the WODEN binary output and shove into a numpy arrays, ready to be put
    into a uvfits file'''

    with open(filename,'rb') as f:
        read_data = f.read()
    f.close()

    data = frombuffer(read_data,dtype=float32)

    n_data = num_time_steps * num_baselines
    v_container = zeros((n_data,1,1,num_freq_channels,4,3))
    uus = zeros(n_data)
    vvs = zeros(n_data)
    wws = zeros(n_data)

    num_visi = num_time_steps * num_freq_channels * num_baselines

    u_base = 0
    v_base = num_visi
    w_base = 2*num_visi
    re_XX_base = 3*num_visi
    im_XX_base = 4*num_visi
    re_XY_base = 5*num_visi
    im_XY_base = 6*num_visi
    re_YX_base = 7*num_visi
    im_YX_base = 8*num_visi
    re_YY_base = 9*num_visi
    im_YY_base = 10*num_visi

    num_cols = 11
    for time_ind in arange(num_time_steps):

        time_step = num_baselines * time_ind * num_freq_channels
        u_ind = u_base + time_step
        v_ind = v_base + time_step
        w_ind = w_base + time_step

        uus[time_ind*num_baselines:(time_ind + 1)*num_baselines] = data[u_ind:u_ind+num_baselines] / VELC
        vvs[time_ind*num_baselines:(time_ind + 1)*num_baselines] = data[v_ind:v_ind+num_baselines] / VELC
        wws[time_ind*num_baselines:(time_ind + 1)*num_baselines] = data[w_ind:w_ind+num_baselines] / VELC

        for freq_ind in arange(num_freq_channels):

            freq_step = num_baselines * (time_ind * num_freq_channels + freq_ind)

            real_XX_ind = re_XX_base + freq_step
            imag_XX_ind = im_XX_base + freq_step
            real_YY_ind = re_YY_base + freq_step
            imag_YY_ind = im_YY_base + freq_step
            real_XY_ind = re_XY_base + freq_step
            imag_XY_ind = im_XY_base + freq_step
            real_YX_ind = re_YX_base + freq_step
            imag_YX_ind = im_YX_base + freq_step

            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,0,0] = data[real_XX_ind:real_XX_ind+num_baselines]
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,0,1] = data[imag_XX_ind:imag_XX_ind+num_baselines]
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,1,0] = data[real_YY_ind:real_YY_ind+num_baselines]
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,1,1] = data[imag_YY_ind:imag_YY_ind+num_baselines]

            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,2,0] = data[real_XY_ind:real_XY_ind+num_baselines]
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,2,1] = data[imag_XY_ind:imag_XY_ind+num_baselines]
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,3,0] = data[real_YX_ind:real_YX_ind+num_baselines]
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,3,1] = data[imag_YX_ind:imag_YX_ind+num_baselines]

            ##Set the weight for everything to one
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,0,2] = ones(num_baselines)
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,1,2] = ones(num_baselines)
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,2,2] = ones(num_baselines)
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,3,2] = ones(num_baselines)

    return uus, vvs, wws, v_container

def write_json(num_time_steps=None, num_freqs=None,
               band_nums=None, json_name=None, freq_res=None,
               time_res=None, jd_date=None, args=None):
    '''Populate a json parameter file used to run WODEN'''

    outfile = open(json_name,'w+')

    outfile.write('{\n')
    outfile.write('  "ra0": %.10f,\n' %args.ra0)
    outfile.write('  "dec0": %.10f,\n' %args.dec0)
    outfile.write('  "num_freqs": %d,\n' %num_freqs)
    outfile.write('  "num_time_steps": %d,\n' %num_time_steps)
    outfile.write('  "cat_filename": "%s",\n' %args.cat_filename)
    outfile.write('  "metafits_filename": "%s",\n' %args.metafits_filename)
    outfile.write('  "time_res": %.5f,\n' %time_res)
    outfile.write('  "frequency_resolution": %.3f,\n' %freq_res)
    outfile.write('  "chunking_size": %d,\n' %args.chunking_size)
    outfile.write('  "jd_date": %.16f,\n' %jd_date)

    if args.sky_crop_components:
        outfile.write('  "sky_crop_components": True,\n')

    if args.use_gaussian_beam:
        outfile.write('  "use_gaussian_beam": True,\n')

    if args.gauss_beam_FWHM:
        outfile.write('  "gauss_beam_FWHM": %.10f,\n' %float(args.gauss_beam_FWHM))

    if args.use_gaussian_beam:
        outfile.write('  "gauss_beam_ref_freq": %.10f,\n' %float(args.gauss_beam_ref_freq))

    if args.use_FEE_beam:
        outfile.write('  "use_FEE_beam": True,\n')
        outfile.write('  "hdf5_beam_path": "%s",\n' %args.hdf5_beam_path)

    if args.EDA2_sim:
        outfile.write('  "EDA2_sim": True,\n')

    if args.array_layout:
        outfile.write('  "array_layout": "%s",\n' %args.array_layout)

    if len(band_nums) == 1:
        band_str = '[%d]' %band_nums[0]
    else:

        band_str = '[%d' %band_nums[0]
        for band in band_nums[1:-1]:
            band_str += ',%d' %band
        band_str += ',%d]' %band_nums[-1]

    outfile.write('  "band_nums": %s\n' %band_str)
    outfile.write('}\n')
    outfile.close()


def make_baseline_date_arrays(num_antennas, date, num_time_steps, time_res):
    """Makes the BASELINE and DATE arrays needed in the uvfits file
    The BASELINE array encode which two antennas formed the baseline
    The DATE array contains the fractional jd date, that is added to the
    header value PZERO5, to specify the time each visibility was recorded at

    Arguements:
        num_antennas - the number of antennas in the antenna table
        date - initial UTC date in format YYYY-MM-DDThh:mm:ss
        num_time_steps - number of time steps in the data
        time_res - integration time of the data (seconds)
    """

    num_baselines = int(((num_antennas - 1)*num_antennas) / 2)
    template_baselines = empty(num_baselines)

    ##Loop over all antenna combinations and encode the baseline pair
    baseline_ind = 0
    for b1 in arange(num_antennas - 1):
        for b2 in arange(b1+1,num_antennas):
            template_baselines[baseline_ind] = RTS_encode_baseline(b1+1, b2+1)
            baseline_ind += 1

    ##Calculate the Julian date, which get's split up into the header (int_jd)
    ##and DATE array (float_jd)
    ##array in the
    int_jd, float_jd = calc_jdcal(initial_date)
    jd_date = int_jd + float_jd

    ##Need an array the length of number of baselines worth of the fractional jd date
    float_jd_array = ones(num_baselines)*float_jd

    ##Create empty data structures for final uvfits file
    n_data = num_time_steps * num_baselines
    baselines_array = zeros(n_data)
    date_array = zeros(n_data)

    for time_ind,time in enumerate(arange(0,num_time_steps*time_res,time_res)):
        time_ind_lower = time_ind*num_baselines
        baselines_array[time_ind_lower:time_ind_lower+num_baselines] = template_baselines

        ##Fill in the fractional julian date, after adding on the appropriate amount of
        ##time - /(24*60*60) because julian number is a fraction of a whole day
        adjust_float_jd_array = float_jd_array + (float(time) / (24.0*60.0*60.0))
        date_array[time_ind_lower:time_ind_lower+num_baselines] = adjust_float_jd_array

    return baselines_array, date_array

if __name__ == "__main__":
    import argparse

    ##Find out where the git repo is, cd in and grab the git label
    ##TODO do this in a better way
    fileloc = os.path.realpath(__file__)
    cwd = os.getcwd()
    os.chdir(('/').join(fileloc.split('/')[:-1]))
    gitlabel = check_output(["git", "describe", "--always"],universal_newlines=True).strip()
    ##Get back to where we were before
    os.chdir(cwd)

    print("You are using WODEN commit %s" %gitlabel)

    parser = argparse.ArgumentParser(description='Run the woden simulator')
    parser.add_argument('--ra0', type=float,
        help='RA of the desired phase centre (deg)')
    parser.add_argument('--dec0', type=float,
        help='Dec of the desired phase centre (deg)')
    parser.add_argument('--cat_filename',
        help='RTS-v2-like srclist to simulate')
    parser.add_argument('--metafits_filename',
        help='Metafits file to base the simulation on')
    parser.add_argument('--band_nums', default='all',
        help='Defaults to running all 24 course bands. Alternatively, enter required numbers delineated by commas, e.g. --band_nums=1,7,9')
    parser.add_argument('--output_uvfits_prepend',default='output',
        help='Prepend name for uvfits - will append band%%02d.uvfits %%band_num at the end. Defaults to "output".')
    parser.add_argument('--num_freq_channels', default='obs',
        help='Number of fine frequency channels to simulate - defaults to 1.28MHz / --freq_res')
    parser.add_argument('--num_time_steps', default='obs',
        help='The number of time steps to simualte - defaults to how many are in the metafits')
    parser.add_argument('--freq_res', type=float,default=False,
        help='Fine channel frequnecy resolution (Hz) - will default to what is in the metafits')
    parser.add_argument('--time_res', type=float,default=False,
        help='Time resolution (s) - will default to what is in the metafits')
    parser.add_argument('--no_tidy', default=False, action='store_true',
        help='Defaults to deleting output binary files from woden and json files. Add this flag to not delete those files')
    parser.add_argument('--nvprof', default=False, action='store_true',
        help='Add to switch on the nvidia profiler when running woden')
    parser.add_argument('--sky_crop_components', default=False, action='store_true',
        help='WODEN will crop out sky model information that is below the horizon for the given LST. By default, for each SOURCE in the sky model, if any COMPONENT is below the horizon, the entire source will be flagged. If --sky_crop_components is included WODEN will include any COMPONENT above the horizon, regardless of which SOURCE it belongs to.')
    parser.add_argument('--use_gaussian_beam', default=False, action='store_true',
        help='Apply a gaussian beam centred on the pointing centre in the metafits file')
    parser.add_argument('--gauss_beam_FWHM', default=False,
        help='The FWHM of the Gaussian beam in deg - WODEN defaults to using 20 deg if this is not set')
    parser.add_argument('--gauss_beam_ref_freq', default=False,
        help='The frequency at which the gauss beam FWHM is set at. If not set, WODEN will default to 150MHz.')
    parser.add_argument('--use_FEE_beam', default=False, action='store_true',
        help='Use the FEE MWA beam model, based on the settings in the metafits file')
    parser.add_argument('--hdf5_beam_path', default='/home/jline/software/useful/MWA_embedded_element_pattern_V02.h5',
        help='Location of the hdf5 file holding the FEE beam coefficients')
    parser.add_argument('--chunking_size', type=int, default=0, help='The chunk size to break up the point sources into for processing - defaults to 0 (do not perform chunking)')

    parser.add_argument('--EDA2_sim', default=False, action='store_true',
        help='If doing an EDA2 simulation, need this flag to read all the antenna positions from a modified metafits file')

    parser.add_argument('--telescope_name', default='MWA',
        help='Name of telescope written out to the uvfits file, defaults to MWA')

    parser.add_argument('--array_layout', default=False,
        help='Instead of reading the array layout from the metafits file, read from a text file. \
              Store antenna positions as offset from array centre, in east, north, height coords (metres)')

    args = parser.parse_args()

    # def false_test(option,name):
    #     if option == False:
    #         print '-----------------------------------------------------'
    #         print '%s not entered but is required. Exiting now!!' %name
    #         print '-----------------------------------------------------'
    #         exit()
    #     else:
    #         return option

    if args.band_nums == 'all':
        band_nums = range(1,25)
    else:
        try:
            band_nums = list(array(args.band_nums.split(','),dtype=int))
        except:
            print('-----------------------------------------------------')
            print('Failed to convert --band_nums into something sensible. Exiting now!!')
            print('-----------------------------------------------------')
            exit()

    ##Setup simulation parameters
    num_time_steps = args.num_time_steps
    num_freq_channels = args.num_freq_channels

    with fits.open(args.metafits_filename) as f:
        initial_date = f[0].header['DATE-OBS']

        ##Get the east, north, height antenna positions from the metafits
        east = f[1].data['East']
        north = f[1].data['North']
        height = f[1].data['Height']

        ##Read observation parameters from the metafits file
        time_res = float(f[0].header['INTTIME'])
        ch_width = float(f[0].header['FINECHAN'])*1e+3
        freqcent = float(f[0].header['FREQCENT'])*1e+6
        b_width = float(f[0].header['BANDWDTH'])*1e+6
        base_low_freq = freqcent - (b_width/2) - (ch_width/2)

        ##Replace default obs settings if specified in arguments
        if args.time_res:
            time_res = args.time_res

        if args.freq_res:
            ch_width = args.freq_res

        if args.num_freq_channels == 'obs':
            num_freq_channels = int(floor(1.28e+6) / ch_width)
        else:
            num_freq_channels = int(args.num_freq_channels)

        if args.num_time_steps == 'obs':
            num_time_steps = int(f[0].header['NSCANS'])
        else:
            num_time_steps = int(args.num_time_steps)

        f.close()

    int_jd, float_jd = calc_jdcal(initial_date)
    jd_date = int_jd + float_jd

    ##Write json file
    json_name = 'run_woden_%s.json' %args.band_nums
    write_json(num_time_steps=num_time_steps, num_freqs=num_freq_channels,
               band_nums=band_nums, json_name=json_name, freq_res=ch_width,
               time_res=time_res, jd_date=jd_date, args=args)

    ##Check the uvfits prepend to make sure we end in .uvfits
    output_uvfits_prepend = args.output_uvfits_prepend
    if output_uvfits_prepend[-7:] == '.uvfits': output_uvfits_prepend = output_uvfits_prepend[-7:]

    if args.nvprof:
        command('nvprof --log-file %s_%02d_nvprof.txt  %s/woden %s' %(output_uvfits_prepend,band_nums[0],WODEN_DIR,json_name))
    else:
        command('%s/woden %s' %(WODEN_DIR,json_name))

    # if args.EDA2_sim:
    #     ##If doing EDA2, we have 256 stations, so select them all but one
    #     selection = arange(len(east) - 1)
    # else:
    #
    #
    # ##TODO if we work out a different way to input east,north,height layout
    # ##this will need to change
    # num_antennas = int(len(selection))
    #
    # ##Prepare the uvfits information
    # ##Create and fill a layout array
    # array_layout = zeros((num_antennas,3))
    # ##Tiles are listed as YY,XX,YY,XX,etc so only use half positions
    #
    # array_layout[:,0] = east[selection]
    # array_layout[:,1] = north[selection]
    # array_layout[:,2] = height[selection]


    if args.array_layout:
        try:
            array_layout = loadtxt(args.array_layout)
            num_antennas,_ = array_layout.shape

            east = array_layout[:,0]
            north = array_layout[:,1]
            height = array_layout[:,2]

        except:
            exit("Could not open array layout file:\n%s\nexiting before woe beings")

    else:
        ##This is an MWA simulation, in the metafits it lists XX,YY for each
        ##antenna so we select every second one
        selection = arange(0,len(east),2)
        num_antennas = int(len(selection))

        east = east[selection]
        north = north[selection]
        height = height[selection]

        array_layout = zeros((num_antennas,3))

        array_layout[:,0] = east
        array_layout[:,1] = north
        array_layout[:,2] = height


    X,Y,Z = enh2xyz(east, north, height,MWA_LAT*D2R)

    ##Get the central frequency channels, used in the uvfits header
    central_freq_chan = int(floor(num_freq_channels / 2.0))

    ##Loop over coarse frequency band
    for band in band_nums:
        print('Converting binary to uvfits band',band)

        output_uvfits_name = output_uvfits_prepend + '_band%02d.uvfits' %band

        band_low_freq = base_low_freq + (band - 1)*1.28e+6
        central_freq_chan_value = band_low_freq + (1.28e+6 / 2.0)
        num_baselines = int(((num_antennas - 1)*num_antennas) / 2)

        filename = "output_visi_band%02d.dat" %band
        uus,vvs,wws,v_container = load_data(filename=filename,num_baselines=num_baselines,
                                            num_freq_channels=num_freq_channels,num_time_steps=num_time_steps)

        ##X,Y,Z are stored in a 2D array in units of seconds in the uvfits file
        XYZ_array = empty((num_antennas,3))
        XYZ_array[:,0] = X / VELC
        XYZ_array[:,1] = Y / VELC
        XYZ_array[:,2] = Z / VELC

        hdu_ant = make_antenna_table(XYZ_array=XYZ_array,telescope_name=args.telescope_name,
                      num_antennas=num_antennas, freq_cent=central_freq_chan_value,
                      date=initial_date)

        baselines_array, date_array = make_baseline_date_arrays(num_antennas,
                                      initial_date, num_time_steps, time_res)

        create_uvfits(v_container=v_container,freq_cent=central_freq_chan_value,ra_point=args.ra0,dec_point=args.dec0,
                    output_uvfits_name=output_uvfits_name,uu=uus,vv=vvs,ww=wws,baselines_array=baselines_array,
                    date_array=date_array,date=initial_date,central_freq_chan=central_freq_chan,ch_width=ch_width,
                    int_jd=int_jd, hdu_ant=hdu_ant, gitlabel=gitlabel)


        ##Tidy up or not
        if args.no_tidy:
            pass
        else:
            command("rm %s" %filename)

    ##Tidy up or not
    if args.no_tidy:
        pass
    else:
        command("rm %s" %json_name)
