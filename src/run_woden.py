#!/usr/bin/env python
"""Wrapper script to generate json input files for, and to run,
the GPU WODEN code. Author: J.L.B. Line
"""
from __future__ import print_function
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import EarthLocation
from astropy import units as u

from erfa import gd2gc

import numpy as np
from struct import unpack
from subprocess import call, check_output
import os
import warnings
import sys

##Constants
R2D = 180.0 / np.pi
D2R = np.pi / 180.0
MWA_LAT = -26.7033194444
MWA_LONG = 116.670813889
MWA_HEIGHT = 377.0
VELC = 299792458.0
SOLAR2SIDEREAL = 1.00274

read_the_docs_build = os.environ.get('READTHEDOCS', None) == 'True'

if read_the_docs_build:
    WODEN_DIR = "not-needed"
else:
    WODEN_DIR = os.environ['WODEN_DIR']

def command(cmd):
    """
    Runs the command string `cmd` using `subprocess.call`

    Parameters
    ----------
    cmd : string
         The command to run on the command line
    """
    call(cmd,shell=True)
    # print(cmd)

def calc_jdcal(date):
    """Takes a string calendar date-time and returns julian date by using
    `astropy.time.Time`_, so date can be formatted anyway `astropy.time.Time`_
    accepts. Returns the date in two parts, an integer day, and a fractional day.
    The header of a uvfits file takes the integer part, and the fractional
    part goes into the DATE array

    .. _astropy.time.Time: https://docs.astropy.org/en/stable/time/

    Parameters
    ----------
    date : string
        UTC date/time (e.g. in format YYYY-MM-DDThh:mm:ss or similar)

    Returns
    -------
    jd_day : float
        floor of the julian date
    jd_fraction : float
        remaining fraction of the julian date

    """

    t = Time(date)
    jd = t.jd

    jd_day = np.floor(jd)
    jd_fraction = (jd - jd_day)

    ##The header of the uvdata file takes the integer, and
    ##then the fraction goes into the data array for PTYPE5
    return jd_day, jd_fraction

def get_LST(latitude=None,longitude=None,date=None,height=None):
    """
    For the given Earth location and UTC date return the local sidereal time
    in degrees. Uses `astropy.time.Time`_ and `astropy.coordinates.EarthLocation`_


    .. _astropy.time.Time: https://docs.astropy.org/en/stable/time/
    .. _astropy.coordinates.EarthLocation: https://docs.astropy.org/en/stable/api/astropy.coordinates.EarthLocation.html?highlight=EarthLocation


    Parameters
    ----------
    latitude : float
        Latitude of location on Earth (deg)
    longitude : float
        Longitude of location on Earth (deg)
    date : string
            UTC date/time (e.g. in format YYYY-MM-DDThh:mm:ss or similar)

    Returns
    -------

    """
    ##Setup location
    observing_location = EarthLocation(lat=latitude*u.deg, lon=longitude*u.deg, height=height)
    ##Setup time at that locatoin
    observing_time = Time(date, scale='utc', location=observing_location)
    ##Grab the LST
    LST = observing_time.sidereal_time('apparent')
    return LST.value*15.0

def RTS_encode_baseline(b1, b2):
    """The ancient aips/miriad extended way of encoding a baseline by antenna
    number, by multiplying antenna number 1 by 256, and adding the second
    antenna number. (e.g. `b1*256 + b2`). Needed for populating the uvfits files.

    Uses the RTS extension for antennas higher that 255, where encoding happens
    as `b1*2048 + b2 + 65536`

    Parameters
    ----------
    b1 : int
        Index of first antenna
    b2 : int
        Index of second antenna

    Returns
    -------
    baseline_number : int
        Encdoded baseline number for uvfits 'BASELINE' array

    """
    if b2 > 255:
        return b1*2048 + b2 + 65536
    else:
        return b1*256 + b2

def RTS_decode_baseline(blcode):
    """The ancient aips/miriad extended way of decoding a baseline. Takes
    the baseline code from the 'BASELINE' array of a uvfits, and returns the
    index of antennas 1 and 2 that form the baseline.

    Parameters
    ----------
    blcode : int
        Baseline code from a uvfits file encoded the two antennas

    Returns
    -------
    b1 : int
        Index of first antenna
    b2 : int
        Index of second antenna

    """
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
                       freq_cent=None,date=None,
                       longitude=None, latitude=None, array_height=None):
    """Write an antenna table for a uvfits file. This is the first table in
    the uvfits file that encodes antenna positions, with some header keywords.
    Uses `astropy.io.fits.BinTableHDU`_ to create the table.

    .. _astropy.io.fits.BinTableHDU: https://docs.astropy.org/en/stable/io/fits/api/tables.html

    Parameters
    ----------
    XYZ_array : float array
        Array with shape = (num_antennas, 3), containg the local :math:`X,Y,Z`
        coorindates of the antenna locations (meters)
    telescope_name : string
        Name to assign to the 'ARRNAM' header keyword
    num_antennas : int
        Number of antennas in the array
    freq_cent : float
        Central frequency to assign to the 'FREQ' header keyword (Hz)
    date : string
        UTC date/time in format YYYY-MM-DDThh:mm:ss to give to the 'RDATE'
        header keyword
    longitude : float
        Longitude of the array (deg)
    latitude : float
        Latitude of the array (deg)
    array_height : float
        Height of the array above sea level (metres)

    Returns
    -------
    hdu_ant : `astropy.io.fits.hdu.table.BinTableHDU`
        Populated uvfits antenna table
    """

    ##Make some values for certain columns
    annnames = np.array(["%05d" %ant for ant in range(1,num_antennas+1)])
    xlabels = np.array(['X']*num_antennas)
    ylabels = np.array(['Y']*num_antennas)

    ##Make a number of FITS columns to create the antenna table from
    col1 = fits.Column(array=annnames,name='ANNAME',format='5A')
    col2 = fits.Column(array=XYZ_array,name='STABXYZ',format='3D')
    ##col3 makes an empty array, and the format states skip reading this column
    ##Just replicating the example uvfits I've been using
    col3 = fits.Column(array=np.array([]),name='ORBPARM',format='0D')
    col4 = fits.Column(array=np.arange(1,num_antennas+1),name='NOSTA',format='1J')
    col5 = fits.Column(array=np.zeros(num_antennas),name='MNTSTA',format='1J')
    col6 = fits.Column(array=np.zeros(num_antennas),name='STAXOF',format='1E')
    col7 = fits.Column(array=xlabels,name='POLTYA',format='1A')
    col8 = fits.Column(array=np.zeros(num_antennas),name='POLAA',format='1E')
    col9 = fits.Column(array=np.zeros(num_antennas),name='POLCALA',format='1E')
    col10 = fits.Column(array=ylabels,name='POLTYB',format='1A')
    col11 = fits.Column(array=np.zeros(num_antennas),name='POLAB',format='1E')
    col12 = fits.Column(array=np.zeros(num_antennas),name='POLCALB',format='1E')

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
    ##Absolute reference point of the centre of the array
    ##Use erfa to calculate this
    arrX, arrY, arrZ = gd2gc(1, longitude*D2R, latitude*D2R, array_height)

    hdu_ant.header['ARRAYX']  = arrX
    hdu_ant.header['ARRAYY']  = arrY
    hdu_ant.header['ARRAYZ']  = arrZ

    hdu_ant.header['FREQ']    = freq_cent
    hdu_ant.header['RDATE']   = date
    hdu_ant.header['TIMSYS']  = 'UTC     '
    hdu_ant.header['ARRNAM']  = telescope_name
    hdu_ant.header['NUMORB']  = 0
    hdu_ant.header['NOPCAL']  = 0
    hdu_ant.header['POLTYPE'] = '        '
    hdu_ant.header['CREATOR'] = 'WODEN_uvfits_writer'

    return hdu_ant

def create_uvfits(v_container=None,freq_cent=None,
                  central_freq_chan=None,ch_width=None,
                  ra_point=None, dec_point=None,
                  output_uvfits_name=None,uu=None,vv=None,ww=None,
                  baselines_array=None, date_array=None,
                  int_jd=None, hdu_ant=None, gitlabel=False):
    """
    Takes visibility data read in from WODEN binary files, predefined
    BASELINE and DATE arrays and an antenna table, and writes them out
    together into a `uvfits` file. Uses multiple `astropy.io.fits` routines to
    create the final uvfits file. Uses `GroupData`_ and `GroupsHDU`_ to create
    the data table, and then combines with the antenna table in a uvfits
    via `HDUList`_.

    .. _GroupData: https://docs.astropy.org/en/stable/io/fits/api/hdus.html?highlight=GroupsHDU#groupdata
    .. _GroupsHDU: https://docs.astropy.org/en/stable/io/fits/api/hdus.html?highlight=GroupsHDU#groupshdu
    .. _HDUList: https://docs.astropy.org/en/stable/io/fits/api/hdulists.html?highlight=HDUList#hdulist

    Parameters
    ----------
    v_container : float array
        Data container for the visibility data. Should have
        `shape = (num_time_steps*num_baselines, 1, 1, num_freq_channels, 4, 3)`
        and should contain instrumental linear polarisation visibilities.
        The axes should change as:

        - 1st axis: ordered by baseline (fastest changing) and then time step
          (slowest changing).
        - 4th axis: ordered with increasing frequency
        - 5th axis: ordered by polarisation in the order of XX,YY,XY,YX
        - 6th axis: ordered by real visi part, imaginary visi part, weighting
    freq_cent : float
        Frequency of the central frequency channel (Hz)
    central_freq_chan : int
        Index of the central frequency channel (zero indexed)
    ch_width : float
        Resolutiom of frequency channels (Hz)
    ra_point : float
        Right ascension of the observation, to set header keyword 'OBSRA' with (deg)
    dec_point : float
        Declination of the observation, to set header keyword 'OBSDEC' with (deg)
    output_uvfits_name : string
        Name for the output uvfits file
    uu : float array
        :math`u` coordinates (seconds).
    vv : float array
        :math`v` coordinates (seconds).
    ww : float array
        :math`w` coordinates (seconds).
    baselines_array : int/float array
        Baseline coding needed for 'BASELINE' array
    date_array : float array
        Fractional julian date array to put in 'DATE' array (days)
    int_jd : int
        Integer julian date to put in the header as 'PZERO5'
    hdu_ant : `astropy.io.fits.hdu.table.BinTableHDU`
        Populated uvfits antenna table
    gitlabel : string
        Optional string to add as 'GITLABEL' in the header. Used by WODEN to
        add the git commit of the code for this run
    """

    ##TODO replace all of this with an interface with pyuvdata

    if not uu.shape[0]==vv.shape[0]==ww.shape[0]==baselines_array.shape[0]==date_array.shape[0]==v_container.shape[0]:
        sys.exit("run_woden.create_uvfits: The first dimension of the arrays:\n"
                 "v_container, uu, vv, ww, baselines_array, date_array\n"
                 "must be equal to make a uvfits file. Exiting now.")

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

    # int_jd, float_jd = calc_jdcal(date)
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

def enh2xyz(east, north, height, latitude=MWA_LAT*D2R):
    """
    Takes local east, north, height coords for a given latitude (radians)
    and returns local X,Y,Z coords to put in the uvfits antenna table

    Parameters
    ----------
    east : float
        Local east coorindate (metres)
    north : float
        Local north coorindate (metres)
    height : float
        Local height coorindate (metres)
    latitude : float
        Latitude of the array - defaults to MWA location (radians)

    Returns
    -------
    X : float
        Local X antenna location
    Y : float
        Local Y antenna location
    Z : float
        Local Z antenna location
    """

    sl = np.sin(latitude)
    cl = np.cos(latitude)
    X = -north*sl + height*cl
    Y = east
    Z = north*cl + height*sl
    return X,Y,Z

def load_data(filename=None,num_baselines=None,num_freq_channels=None,num_time_steps=None):
    """
    Read the WODEN binary output and shove into a numpy arrays, ready to be put
    into a uvfits file. Data in WODEN binaries is ordered by baseline (fastest
    changing), frequency, and time (slowest changing). Visibility coords and
    data are read in, with the visi data output into an array of
    `shape=(num_time_steps*num_baselines,1,1,num_freq_channels,4,3))`, which is
    appropriate for a uvfits file.

    Parameters
    ----------
    filename : string
        Name of WODEN binary file to read from
    num_baselines : int
        Number of baselines in the binary file
    num_freq_channels : int
        Number of frequencies in the binary file
    num_time_steps : int
        Number of time steps in the binary file

    Returns
    -------
    uus : float array
        The :math:`u` coordinates (seconds)
    vvs : float array
        The :math:`v` coordinates (seconds)
    wws : float array
        The :math:`w` coordinates (seconds)
    v_container : float array
        Visibility data with
        `shape=(num_time_steps*num_baselines,1,1,num_freq_channels,4,3))`
    """

    with open(filename,'rb') as f:
        read_data = f.read()
    f.close()

    data = np.frombuffer(read_data,dtype=np.float32)

    n_data = num_time_steps * num_baselines
    v_container = np.zeros((n_data,1,1,num_freq_channels,4,3))
    uus = np.zeros(n_data)
    vvs = np.zeros(n_data)
    wws = np.zeros(n_data)

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
    for time_ind in np.arange(num_time_steps):

        time_step = num_baselines * time_ind * num_freq_channels
        u_ind = u_base + time_step
        v_ind = v_base + time_step
        w_ind = w_base + time_step

        uus[time_ind*num_baselines:(time_ind + 1)*num_baselines] = data[u_ind:u_ind+num_baselines] / VELC
        vvs[time_ind*num_baselines:(time_ind + 1)*num_baselines] = data[v_ind:v_ind+num_baselines] / VELC
        wws[time_ind*num_baselines:(time_ind + 1)*num_baselines] = data[w_ind:w_ind+num_baselines] / VELC

        for freq_ind in np.arange(num_freq_channels):

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
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,0,2] = np.ones(num_baselines)
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,1,2] = np.ones(num_baselines)
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,2,2] = np.ones(num_baselines)
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,3,2] = np.ones(num_baselines)

    return uus, vvs, wws, v_container

def write_json(num_time_steps=None, num_freqs=None,
               band_nums=None, json_name=None, freq_res=None,
               time_res=None, jd_date=None,
               lst=None, lowest_channel_freq=None, FEE_delays=None,
               latitude=None, array_layout_name=None,
               args=None):
    """
    Populate and write out a .json parameter file used to run WODEN.
    Is later used on the command line to run WODEN.

    Parameters
    ----------
    num_time_steps : int
        Number of time steps to simulate
    num_freqs : int
        Number of fine frequencies per coarse band to simulate
    band_nums : int array/list
        Which coarse bands to simulte
    json_name : string
        Name out .json file to save outputs to
    freq_res : float
        Frequency resolution of the fine channels (Hz)
    time_res : float
        Time resolution of the simulation (s)
    jd_date : float
        Initial Julian date of simulation (days)
    lst : float
        Local sidereal time of the simulate (degrees)
    lowest_channel_freq : float
        Frequency of the lowest fine channel (of coarse band 1)
    FEE_delays : int/float array
        Delays to use with MWA FEE if simulating with FEE beam
    latitude : float
        Latitude of the array centre (degrees)
    array_layout_name : string
        Path to array text file containing e,n,h array layout
    args : `argparse.Namespace`
        The args as return by `args=parser.parse_args()`, from the `parser`
        returned by  :func:`~run_woden.get_parser`.

    """

    with open(json_name,'w+') as outfile:

        outfile.write('{\n')
        outfile.write('  "ra0": %.10f,\n' %args.ra0)
        outfile.write('  "dec0": %.10f,\n' %args.dec0)
        outfile.write('  "num_freqs": %d,\n' %num_freqs)
        outfile.write('  "num_time_steps": %d,\n' %num_time_steps)
        outfile.write('  "cat_filename": "%s",\n' %args.cat_filename)
        outfile.write('  "time_res": %.5f,\n' %time_res)
        outfile.write('  "frequency_resolution": %.3f,\n' %freq_res)
        outfile.write('  "chunking_size": %d,\n' %int(args.chunking_size))
        outfile.write('  "jd_date": %.16f,\n' %jd_date)
        outfile.write('  "LST": %.8f,\n' %lst)
        outfile.write('  "array_layout": "%s",\n' %array_layout_name)
        outfile.write('  "lowest_channel_freq": %.10e,\n' %lowest_channel_freq)
        outfile.write('  "latitude": %.8f,\n' %latitude)
        outfile.write('  "coarse_band_width": %.10e,\n' %float(args.coarse_band_width))

        if args.sky_crop_components:
            outfile.write('  "sky_crop_components": "True",\n')

        if args.primary_beam == 'Gaussian':
            outfile.write('  "use_gaussian_beam": "True",\n')
            if args.gauss_beam_FWHM:
                outfile.write('  "gauss_beam_FWHM": %.10f,\n' %float(args.gauss_beam_FWHM))

            if args.gauss_beam_ref_freq:
                outfile.write('  "gauss_beam_ref_freq": %.10f,\n' %float(args.gauss_beam_ref_freq))

            if args.gauss_ra_point:
                outfile.write('  "gauss_ra_point": %.8f,\n' %float(args.gauss_ra_point))
            else:
                outfile.write('  "gauss_ra_point": %.8f,\n' %float(args.ra0))
            if args.gauss_dec_point:
                outfile.write('  "gauss_dec_point": %.8f,\n' %float(args.gauss_dec_point))
            else:
                outfile.write('  "gauss_dec_point": %.8f,\n' %float(args.dec0))

        elif args.primary_beam == 'MWA_FEE':
            outfile.write('  "use_FEE_beam": "True",\n')
            outfile.write('  "hdf5_beam_path": "%s",\n' %args.hdf5_beam_path)
            outfile.write('  "FEE_delays": %s,\n ' %FEE_delays)

        elif args.primary_beam == 'EDA2':
            outfile.write('  "use_EDA2_beam": "True",\n')

        if len(band_nums) == 1:
            band_str = '[%d]' %band_nums[0]
        else:

            band_str = '[%d' %band_nums[0]
            for band in band_nums[1:-1]:
                band_str += ',%d' %band
            band_str += ',%d]' %band_nums[-1]

        outfile.write('  "band_nums": %s\n' %band_str)
        outfile.write('}\n')


def make_baseline_date_arrays(num_antennas, date, num_time_steps, time_res):
    """Makes the BASELINE and DATE arrays needed in the uvfits file
    The BASELINE array encode which two antennas formed the baseline
    The DATE array contains the fractional jd date, that is added to the
    header value PZERO5, to specify the time each visibility was recorded at

    Parameters
    -----------
        num_antennas : int
            The number of antennas in the antenna table
        date : string
            Initial UTC date/time (e.g. in format YYYY-MM-DDThh:mm:ss or similar)
        num_time_steps : int
            Number of time steps in the data
        time_res : float
            Integration time of the data (seconds)

    Returns
    -------
    baselines_array : float array
        The uvfits antenna encoded array needed for BASELINE
    date_array : float array
        The fractional part of the julian date needed for DATE (days)
    """

    num_baselines = int(((num_antennas - 1)*num_antennas) / 2)
    template_baselines = np.empty(num_baselines)

    ##Loop over all antenna combinations and encode the baseline pair
    baseline_ind = 0
    for b1 in np.arange(num_antennas - 1):
        for b2 in np.arange(b1+1,num_antennas):
            template_baselines[baseline_ind] = RTS_encode_baseline(b1+1, b2+1)
            baseline_ind += 1

    ##Calculate the Julian date, which get's split up into the header (int_jd)
    ##and DATE array (float_jd)
    ##array in the
    int_jd, float_jd = calc_jdcal(date)

    ##Need an array the length of number of baselines worth of the fractional jd date
    float_jd_array = np.ones(num_baselines)*float_jd

    ##Create empty data structures for final uvfits file
    n_data = num_time_steps * num_baselines
    baselines_array = np.zeros(n_data)
    date_array = np.zeros(n_data)

    for time_ind,time in enumerate(np.arange(0,num_time_steps*time_res,time_res)):
        time_ind_lower = time_ind*num_baselines
        baselines_array[time_ind_lower:time_ind_lower+num_baselines] = template_baselines

        ##Fill in the fractional julian date, after adding on the appropriate amount of
        ##time - /(24*60*60) because julian number is a fraction of a whole day
        adjust_float_jd_array = float_jd_array + (float(time) / (24.0*60.0*60.0))
        date_array[time_ind_lower:time_ind_lower+num_baselines] = adjust_float_jd_array

    return baselines_array, date_array

def remove_phase_tracking(frequencies=None, wws_seconds=None,
                          num_time_steps=None, v_container=None,
                          num_baselines=None):
    """
    WARNING - currently does not change the :math:`u,v,w` coordinates, so they
    are still defined via the original phase centre. This function really is
    just to feed uvfits into the RTS (which generates it's own u,v,w using the
    antenna table)

    Undoes phase tracking applied by WODEN - to phase track, a phase was applied
    to counter the delay term caused by :math:`w` term of baseline - so just
    apply the opposite effect of the w term, i.e.

    .. math::
        V^\\prime = V \\exp(2\pi i w)

    where :math:`V` is the phase tracked visibility and :math:`V^\\prime` is
    the visibility after removing phase tracking.

    Parameters
    ----------

    frequencies : float array
        Frequencies of all fine channels (Hz)
    wws_seconds : float array
        The :math:`w` coordinates (seconds)
    num_baselines : int
        Number of baselines
    v_container : float array
        Complex visibility data out of WODEN with phase tracking, with
        `shape=(num_time_steps*num_baselines,1,1,num_freq_channels,4,3))`

    Returns
    -------
    v_container : float array
        Same visibility data as before, with phase tracking returned.

    """

    sign = 1
    PhaseConst = 2j * pi * sign

    num_freqs = len(frequencies)

    for time_ind in np.arange(num_time_steps):

        these_wws_secs = wws_seconds[time_ind*num_baselines:(time_ind + 1)*num_baselines]

        for freq_ind, freq in enumerate(frequencies):

            xx_re = v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,0,0]
            xx_im = v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,0,1]

            yy_re = v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,1,0]
            yy_im = v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,1,1]

            xy_re = v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,2,0]
            xy_im = v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,2,1]

            yx_re = v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,3,0]
            yx_im = v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,3,1]

            xx_comp = xx_re + 1j*xx_im
            yy_comp = yy_re + 1j*yy_im
            xy_comp = xy_re + 1j*xy_im
            yx_comp = yx_re + 1j*yx_im

            ##theory - so normal phase delay is caused by path difference across
            ##a base line, which is u*l + v*m + w*n
            ##To phase track, you insert a phase to make sure there is no w contribution at
            ##phase centre; this is when n = 1, so you insert a phase thus:
            ##a base line, which is u*l + v*m + w*(n - 1)
            ##So we just need to remove the effect of the -w term

            wws = these_wws_secs * freq
            phase_rotate = np.exp( PhaseConst * wws)
            xx_comp = xx_comp * phase_rotate
            yy_comp = yy_comp * phase_rotate
            xy_comp = xy_comp * phase_rotate
            yx_comp = yx_comp * phase_rotate

            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,0,0] = real(xx_comp)
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,0,1] = imag(xx_comp)
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,1,0] = real(yy_comp)
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,1,1] = imag(yy_comp)
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,2,0] = real(xy_comp)
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,2,1] = imag(xy_comp)
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,3,0] = real(yx_comp)
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,3,1] = imag(yx_comp)

    return v_container

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
    freq_group.add_argument('--coarse_band_width', default=1.28e+6,
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
    time_group.add_argument('--num_time_steps', default='False',
        help='The number of time steps to simualte. Defaults to how many are in'
             'the metafits if using metafits')
    time_group.add_argument('--time_res', type=float,default=False,
        help='Time resolution (s) - will default to what is in the metafits '
              'if the metafits if using metafits')


    obs_group = parser.add_argument_group('OBSERVATION OPTIONS')
    obs_group.add_argument('--ra0', type=float,
        help='RA of the desired phase centre (deg)')
    obs_group.add_argument('--dec0', type=float,
        help='Dec of the desired phase centre (deg)')
    obs_group.add_argument('--date', default=False,
        help='Initial UTC date of the observatio in format YYYY-MM-DDThh:mm:ss '
             'This is used to set the LST and array precession. This is set '
             'automatically when reading a metafits but including this will '
             'override the date in the metafits')

    tel_group = parser.add_argument_group('TELESCOPE OPTIONS')
    tel_group.add_argument('--latitude', default=MWA_LAT, type=float,
        help='Latitude (deg) of the array - defaults to MWA at -26.7033194444')
    tel_group.add_argument('--longitude', default=MWA_LONG, type=float,
        help='Longitude (deg) of the array - defaults to MWA at 116.670813889')
    tel_group.add_argument('--array_height', default=MWA_HEIGHT, type=float,
        help='Height (m) of the array above sea level - defaults to MWA at 377.0')
    tel_group.add_argument('--array_layout', default=False,
        help='Instead of reading the array layout from the metafits file, read'
             ' from a text file. Store antenna positions as offset from array '
             'centre, in east, north, height coords (metres)')
    tel_group.add_argument('--primary_beam', default="MWA_FEE",
        help="R|Which primary beam to use in the simulation.\nOptions are:\n"
            "\t MWA_FEE (MWA fully embedded element model)\n"
            "\t Gaussian (Analytic symmetric Gaussian)\n"
            "\t\t see --gauss_beam_FWHM and\n"
            "\t\t and --gauss_beam_ref_freq for\nfine control)\n"
            "\t EDA2 (Analytic dipole with a ground mesh) \n"
            "\t none (Don't use a primary beam at all)")

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
    tel_group.add_argument('--MWA_FEE_delays',
        help='R|A list of 16 delays to point the MWA FEE primary beam \n'
              'model enter as as list like: \n'
              '--MWA_FEE_delays=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]\n'
              'for a zenith pointing. This is read directly from\n'
              'the metafits if using a metafits file')
    tel_group.add_argument('--telescope_name', default='MWA',
        help='Name of telescope written out to the uvfits file, defaults to MWA')


    input_group = parser.add_argument_group('INPUT/OUTPUT OPTIONS')
    input_group.add_argument('--cat_filename',
        help='Path to WODEN style sky model')
    input_group.add_argument('--metafits_filename',default=False,
        help='MWA style metafits file to base the simulation on. Array layout,'
             ' frequency and time parameters are all set by this option, but '
             'can be overridden using other arguments')
    input_group.add_argument('--output_uvfits_prepend',default='output',
        help='Prepend name for uvfits - will append band%%02d.uvfits %%band_num '
             'at the end. Defaults to "output".')
    input_group.add_argument('--sky_crop_components', default=False, action='store_true',
        help='WODEN will crop out sky model information that is below the '
             'horizon for the given LST. By default, for each SOURCE in the '
             'sky model, if any COMPONENT is below the horizon, the entire '
             'source will be flagged. If --sky_crop_components is included '
             'WODEN will include any COMPONENT above the horizon, regardless '
             'of which SOURCE it belongs to.')

    sim_group = parser.add_argument_group('SIMULATOR OPTIONS')
    sim_group.add_argument('--remove_phase_tracking', default=False, action='store_true',
        help='By adding this flag, remove the phase tracking of the '
             'visibilities - use this to feed uvfits into the RTS')
    sim_group.add_argument('--no_tidy', default=False, action='store_true',
        help='Defaults to deleting output binary files from woden and json '
             'files. Add this flag to not delete those files')
    sim_group.add_argument('--nvprof', default=False, action='store_true',
        help='Add to switch on the nvidia profiler when running woden')
    sim_group.add_argument('--chunking_size', type=float, default=0,
        help='The chunk size to break up the point sources into for processing '
             '- defaults to 0 (do not perform chunking)')

    return parser

if __name__ == "__main__":

    ##Find out where the git repo is, cd in and grab the git label
    ##TODO do this in a better way
    fileloc = os.path.realpath(__file__)
    cwd = os.getcwd()
    os.chdir(('/').join(fileloc.split('/')[:-1]))
    gitlabel = check_output(["git", "describe", "--always"],universal_newlines=True).strip()
    ##Get back to where we were before
    os.chdir(cwd)

    print("You are using WODEN commit %s" %gitlabel)

    ##Grab the parser and parse some args
    parser = get_parser()
    args = parser.parse_args()

    # def false_test(option,name):
    #     if option == False:
    #         print '-----------------------------------------------------'
    #         print '%s not entered but is required. Exiting now!!' %name
    #         print '-----------------------------------------------------'
    #         exit()
    #     else:
    #         return option

    if args.primary_beam not in ['MWA_FEE', 'Gaussian', 'EDA2', 'none', 'None']:
        exit('Primary beam option --primary_beam must be one of:\n'
             '\t MWA_FEE, Gaussian, EDA2, none\n'
             'User has entered --primary_beam={:s}\n'
             'Please fix and try again. Exiting now'.format(args.primary_beam))

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
            except KeyError:
                exit('To use MWA FEE beam, either --hdf5_beam_path or environment\n'
                     'variable MWA_FEE_HDF5 must point towards the file\n'
                     'mwa_full_embedded_element_pattern.h5. Exiting now as WODEN will fail.')

    FEE_delays = False

    if args.metafits_filename:

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
            lowest_channel_freq = freqcent - (b_width/2) - (ch_width/2)

            num_time_steps = int(f[0].header['NSCANS'])

            delays = np.array(f[0].header['DELAYS'].split(','),dtype=int)
            delays[np.where(delays == 32)] = 0
            FEE_delays = str(list(delays))

            ##If user hasn't specified a pointing for a Gaussian beam,
            ##fill in using the metafits file
            if not args.gauss_ra_point:
                args.gauss_ra_point = float(f[0].header['RA'])
            if not args.gauss_dec_point:
                args.gauss_dec_point = float(f[0].header['DEC'])

            f.close()

    else:
        ##TODO put in a check that all other agruments needed are included
        pass

    ##Override metafits/load arguements

    coarse_band_width = float(args.coarse_band_width)

    if args.lowest_channel_freq: lowest_channel_freq = float(args.lowest_channel_freq)
    if args.num_time_steps: num_time_steps = int(args.num_time_steps)

    if args.num_freq_channels == 'obs':
        num_freq_channels = int(np.floor(coarse_band_width / ch_width))
    else:
        num_freq_channels = int(args.num_freq_channels)

    if args.time_res: time_res = args.time_res
    if args.freq_res: ch_width = args.freq_res
    if args.date: initial_date = args.date

    latitude = float(args.latitude)
    longitude = float(args.longitude)
    array_height = float(args.array_height)
    lst_deg = get_LST(latitude=latitude,longitude=longitude,height=array_height,
                      date=initial_date)

    if args.MWA_FEE_delays: FEE_delays = args.MWA_FEE_delays

    if args.band_nums == 'all':
        band_nums = range(1,25)
    else:
        try:
            band_nums = list(np.array(args.band_nums.split(','),dtype=int))
        except:
            print('-----------------------------------------------------')
            print('Failed to convert --band_nums into something sensible. Exiting now!!')
            print('-----------------------------------------------------')
            exit()

    int_jd, float_jd = calc_jdcal(initial_date)
    jd_date = int_jd + float_jd


    if args.array_layout:
        try:
            array_layout = np.loadtxt(args.array_layout)
            num_antennas,_ = array_layout.shape

            east = array_layout[:,0]
            north = array_layout[:,1]
            height = array_layout[:,2]

        except:
            exit("Could not open array layout file:\n"
                 "\t{:s}\nExiting before woe beings".format(args.array_layout))

        array_layout_name = args.array_layout

    else:
        ##Using metafits for array layout. In the metafits it lists XX,YY for each
        ##antenna so we select every second one
        selection = np.arange(0,len(east),2)
        num_antennas = int(len(selection))

        east = east[selection]
        north = north[selection]
        height = height[selection]

        array_layout = np.zeros((num_antennas,3))

        array_layout[:,0] = east
        array_layout[:,1] = north
        array_layout[:,2] = height

        array_layout_name = 'WODEN_array_layout.txt'

        np.savetxt(array_layout_name,array_layout)

    ##Write json file
    json_name = 'run_woden_%s.json' %args.band_nums

    write_json(num_time_steps=num_time_steps, num_freqs=num_freq_channels,
                   band_nums=band_nums, json_name=json_name, freq_res=ch_width,
                   time_res=time_res, jd_date=jd_date,
                   lst=lst_deg, lowest_channel_freq=lowest_channel_freq,
                   FEE_delays=FEE_delays,
                   latitude=latitude,
                   array_layout_name=array_layout_name,args=args)

    ##Check the uvfits prepend to make sure we end in .uvfits
    output_uvfits_prepend = args.output_uvfits_prepend
    if output_uvfits_prepend[-7:] == '.uvfits': output_uvfits_prepend = output_uvfits_prepend[-7:]

    if args.nvprof:
        command('nvprof --log-file %s_%02d_nvprof.txt  %s/woden %s' %(output_uvfits_prepend,band_nums[0],WODEN_DIR,json_name))
    else:
        command('%s/woden %s' %(WODEN_DIR,json_name))



    X,Y,Z = enh2xyz(east, north, height, latitude*D2R)

    ##Get the central frequency channels, used in the uvfits header
    central_freq_chan = int(np.floor(num_freq_channels / 2.0))

    ##Loop over coarse frequency band
    for band in band_nums:
        print('Converting binary to uvfits band',band)

        output_uvfits_name = output_uvfits_prepend + '_band%02d.uvfits' %band

        band_low_freq = lowest_channel_freq + (band - 1)*coarse_band_width
        central_freq_chan_value = band_low_freq + central_freq_chan*ch_width
        num_baselines = int(((num_antennas - 1)*num_antennas) / 2)

        filename = "output_visi_band%02d.dat" %band
        uus,vvs,wws,v_container = load_data(filename=filename,num_baselines=num_baselines,
                                            num_freq_channels=num_freq_channels,num_time_steps=num_time_steps)


        if args.remove_phase_tracking:
            frequencies = band_low_freq + np.arange(num_freq_channels)*ch_width

            v_container = remove_phase_tracking(frequencies=frequencies,
                                      wws_seconds=wws,
                                      num_time_steps=num_time_steps,
                                      v_container=v_container,
                                      num_baselines=num_baselines)

        ##X,Y,Z are stored in a 2D array in units of seconds in the uvfits file
        XYZ_array = np.empty((num_antennas,3))
        XYZ_array[:,0] = X
        XYZ_array[:,1] = Y
        XYZ_array[:,2] = Z

        hdu_ant = make_antenna_table(XYZ_array=XYZ_array,telescope_name=args.telescope_name,
                      num_antennas=num_antennas, freq_cent=central_freq_chan_value,
                      date=initial_date,
                      longitude=longitude, latitude=latitude, array_height=array_height)

        baselines_array, date_array = make_baseline_date_arrays(num_antennas,
                                      initial_date, num_time_steps, time_res)

        create_uvfits(v_container=v_container,freq_cent=central_freq_chan_value,ra_point=args.ra0,dec_point=args.dec0,
                    output_uvfits_name=output_uvfits_name,uu=uus,vv=vvs,ww=wws,baselines_array=baselines_array,
                    date_array=date_array,central_freq_chan=central_freq_chan,ch_width=ch_width,
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
