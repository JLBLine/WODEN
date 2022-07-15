#!/usr/bin/env python3
"""Wrapper script to generate json input files for, and to run,
the GPU WODEN code. Author: J.L.B. Line
"""
from __future__ import print_function
from astropy.io import fits
from astropy.time import Time, TimeDelta
from astropy.coordinates import EarthLocation
from astropy import units as u
from erfa import gd2gc
import numpy as np
from struct import unpack
import subprocess
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

##version of release, a fall back if this isn't in a git repo
VERSION = "1.4.0"

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

def get_uvfits_date_and_position_constants(latitude=None,longitude=None,
                                           date=None,height=None):
    """
    Returns a number of date and time based values that are needed for uvfits
    headers. For the given Earth location and UTC date return the local sidereal
    time (deg), the Greenwich sidereal time at 0 hours on the given date (deg),
    the rotational speed of Earth on the given date (in degrees per day), and
    the difference between UT1 and UTC.
    Uses `astropy.time.Time`_ and `astropy.coordinates.EarthLocation`_ to make
    the calculations.

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
    LST_deg : float
        Local sidereal time (degrees)
    GST0_deg : float
        Greenwich sidereal time at 0 hours on the given date (degrees)
    DEGPDY : float
        Rotational speed of Earth on the given date (degrees per day)
    ut1utc : float
        Difference between UT1 and UTC (secs)
    """

    ##Setup location
    observing_location = EarthLocation(lat=latitude*u.deg, lon=longitude*u.deg, height=height)
    ##Setup time at that location

    observing_time = Time(date, scale='utc', location=observing_location)
    ##Grab the LST
    LST = observing_time.sidereal_time('apparent')
    LST_deg = LST.value*15.0

    ##uvfits file needs to know the greenwich sidereal time at 0 hours
    ##on the date in question
    zero_date = date.split('T')[0] + "T00:00:00"
    zero_time = Time(zero_date, scale='utc', location=observing_location)
    GST0 = zero_time.sidereal_time('apparent', 'greenwich')
    GST0_deg = GST0.value*15.0

    ##It also needs to know the rotational rate of the Earth on that day, in
    ##units of degrees per day
    ##Do this by measuring the LST exactly a day later

    date_plus_one_day =  observing_time + TimeDelta(1*u.day)
    LST_plusone = date_plus_one_day.sidereal_time('apparent')

    LST_plusone_deg = LST_plusone.value*15.0

    DEGPDY = 360.0 + (LST_plusone_deg - LST_deg)

    ut1utc = float(observing_time.delta_ut1_utc)

    return LST_deg, GST0_deg, DEGPDY, ut1utc

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


def make_antenna_table(XYZ_array=None, telescope_name=None,num_antennas=None,
                       freq_cent=None, date=None, gst0_deg=None, degpdy=None,
                       ut1utc=None, longitude=None, latitude=None, array_height=None):
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
    gst0_deg : float
        Greenwich sidereal time at 0 hours on the given date (degrees)
    degpdy : float
        Rotational speed of Earth on the given date (degrees per day)
    ut1utc : float
        Difference between UT1 and UTC (secs)
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
    col1 = fits.Column(array=annnames,name='ANNAME',format='8A')
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

    hdu_ant.header['FREQ'] = freq_cent
    hdu_ant.header['RDATE'] = date
    hdu_ant.header['GSTIA0'] = gst0_deg
    hdu_ant.header['DEGPDY'] = degpdy

    hdu_ant.header['UT1UTC'] = ut1utc
    hdu_ant.header['XYZHAND'] = 'RIGHT'
    hdu_ant.header['FRAME'] = '????'

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
                  longitude=None, latitude=None, array_height=None,
                  telescope_name=None,
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
        - 2nd, 3rd axes: essentially do nothing in these uvfits, are placeholders
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
        Integer julian date to put in the header as 'PZERO4'
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

    antenna1_array = np.empty(len(baselines_array))
    antenna2_array = np.empty(len(baselines_array))

    for antind, baseline in enumerate(baselines_array):
        ant1, ant2 = RTS_decode_baseline(baseline)

        antenna1_array[antind] = ant1
        antenna2_array[antind] = ant2

    ##stick a bunch of ones in why not
    subarray = np.ones(len(baselines_array))

    uvparnames = ['UU','VV','WW','DATE','BASELINE', 'ANTENNA1', 'ANTENNA2', 'SUBARRAY']
    parvals = [uu,vv,ww,date_array,baselines_array, antenna1_array, antenna2_array, subarray]

    # Optional INTTIM length of time data were integrated over (seconds)

    uvhdu = fits.GroupData(v_container,parnames=uvparnames,pardata=parvals,bitpix=-32)
    uvhdu = fits.GroupsHDU(uvhdu)

    ## Write the parameters scaling explictly because they are omitted if default 1/0
    uvhdu.header['PSCAL1'] = 1.0
    uvhdu.header['PZERO1'] = 0.0
    uvhdu.header['PSCAL2'] = 1.0
    uvhdu.header['PZERO2'] = 0.0
    uvhdu.header['PSCAL3'] = 1.0
    uvhdu.header['PZERO3'] = 0.0
    uvhdu.header['PSCAL4'] = 1.0
    uvhdu.header['PZERO4'] = float(int_jd)
    uvhdu.header['PSCAL5'] = 1.0
    uvhdu.header['PZERO5'] = 0.0
    uvhdu.header['PSCAL6'] = 1.0
    uvhdu.header['PZERO6'] = 0.0
    uvhdu.header['PSCAL7'] = 1.0
    uvhdu.header['PZERO7'] = 0.0
    uvhdu.header['PSCAL8'] = 1.0
    uvhdu.header['PZERO8'] = 0.0

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

    ##We're outputting into J2000
    uvhdu.header['EPOCH'] = 2000.0

    ##Old observation parameters that were/are needed in CHIPS
    uvhdu.header['OBJECT']  = 'Undefined'
    uvhdu.header['OBSRA']   = ra_point
    uvhdu.header['OBSDEC']  = dec_point

    uvhdu.header['TELESCOP'] = telescope_name
    uvhdu.header['LAT']      = latitude
    uvhdu.header['LON']      = longitude
    uvhdu.header['ALT']      = array_height
    ## For everything WODEN can simulate, there are no extra instruments on
    ## the telescope (I guess this is for different feed horns on a dish and
    ## similar things to that?) so just set it to the telescope name
    uvhdu.header['INSTRUME'] = telescope_name

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

def load_data(filename=None,num_baselines=None,num_freq_channels=None,num_time_steps=None,
              precision=None):
    """
    Read the WODEN binary output and shove into a numpy arrays, ready to be put
    into a uvfits file. Data in WODEN binaries is ordered by baseline (fastest
    changing), frequency, and time (slowest changing). Visibility coords and
    data are read in, with the visi data output into an array of
    `shape=(num_time_steps*num_baselines,1,1,num_freq_channels,4,3))`, which is
    appropriate for a uvfits file. Needs to know whether WODEN was run with
    'float' (32 bit) or 'double' (64 bit) precision to read in the data
    correctly.

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
    precision : string
        Precision WODEN was run with - either 'float' or 'double'

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

    ##numpy needs to know if we have 32 (float) or 64 (double) bit precision
    if precision == 'float':
        data = np.frombuffer(read_data,dtype=np.float32)
    elif precision == 'double':
        data = np.frombuffer(read_data,dtype=np.float64)

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

        if len(args.band_nums) == 1:
            band_str = '[%d]' %args.band_nums[0]
        else:

            band_str = '[%d' %args.band_nums[0]
            for band in args.band_nums[1:-1]:
                band_str += ',%d' %band
            band_str += ',%d]' %args.band_nums[-1]

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
        ##Should be the central time of the integration so add half a time resolution
        adjust_float_jd_array = float_jd_array + ((float(time) + time_res/2.0) / (24.0*60.0*60.0))
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
    PhaseConst = 2j * np.pi * sign

    num_freqs = len(frequencies)

    # print("FREQS", frequencies)

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

            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,0,0] = np.real(xx_comp)
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,0,1] = np.imag(xx_comp)
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,1,0] = np.real(yy_comp)
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,1,1] = np.imag(yy_comp)
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,2,0] = np.real(xy_comp)
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,2,1] = np.imag(xy_comp)
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,3,0] = np.real(yx_comp)
            v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,3,1] = np.imag(yx_comp)

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
    input_group.add_argument('--cat_filename', required=True,
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
            args.east = f[1].data['East']
            args.north = f[1].data['North']
            args.height = f[1].data['Height']

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

    return args

def get_code_version(woden_exe):
    """
    Checks the environment variables and absolute path of `run_woden.py` to
    find out where the script lives on the current machine. If the code
    resides within a git repo, query to the git repo to find out the current
    commit tag. If not a git repo (i.e. a release version) uses the stored
    VERSION string to at least supply a version.

    Parameters
    ----------
    woden_exe : string
        The name of the woden executable to be used. Either 'woden_float' or
        'woden_double'

    Returns
    -------
    gitlabel : string
        Either the git commit, release version, or a note saying unknown (the
        version could not be found)
    WODEN_DIR : string
        The directory in which `woden_exe` lives in
    """

    WODEN_DIR = 'unset'

    if read_the_docs_build:
        WODEN_DIR = "not-needed"
    else:
        ##If the user is using 'init_WODEN.sh' in their bashrc, look for it
        try:
            WODEN_DIR = os.environ['WODEN_DIR']
        ##If it doesn't exist, try and find it in the path
        except KeyError:
            for path in os.environ["PATH"].split(os.pathsep):
                woden_exe_path = os.path.join(path, woden_exe)
                if os.path.isfile(woden_exe_path):
                    print(os.access(woden_exe_path, os.X_OK))

                    WODEN_DIR = '/'.join(woden_exe_path.split('/')[:-1])

    ##Print where we found woden executable
    print(f'Using the WODEN living here: {WODEN_DIR}/{woden_exe}')

    ##Find out where the git repo is, cd in and grab the git label
    ##TODO do this in a better way
    fileloc = os.path.realpath(__file__)
    cwd = os.getcwd()
    os.chdir(('/').join(fileloc.split('/')[:-1]))

    try:
        gitlabel = subprocess.check_output(["git", "describe", "--always"],
                                       stderr=subprocess.STDOUT,
                                       universal_newlines=True).strip()

    except subprocess.CalledProcessError as err:
        gitlabel = f"Release {VERSION}"


    ##Get back to where we were before
    os.chdir(cwd)

    print(f"You are using WODEN commit: {gitlabel}")

    return gitlabel, WODEN_DIR

if __name__ == "__main__":

    ##If we're at readthe docs, we don't need to know where woden exe lives
    read_the_docs_build = os.environ.get('READTHEDOCS', None) == 'True'

    ##Grab the parser and parse some args
    parser = get_parser()
    args = parser.parse_args()

    args = check_args(args)

    if args.precision == 'float':
        woden_exe = "woden_float"
    elif args.precision == 'double':
        woden_exe = "woden_double"

    ##Find out what git/release version we are using, and where the code lives
    gitlabel, WODEN_DIR = get_code_version(woden_exe)

    lst_deg, gst0_deg, degpdy, ut1utc = get_uvfits_date_and_position_constants(latitude=args.latitude, longitude=args.longitude,
                        height=args.array_height, date=args.date)

    int_jd, float_jd = calc_jdcal(args.date)
    jd_date = int_jd + float_jd

    ##Write json file
    json_band_str =  '-'.join(map(str, args.band_nums))
    json_name = 'run_woden_{:s}.json'.format(json_band_str)

    write_json(json_name=json_name, jd_date=jd_date, lst=lst_deg, args=args)

    ##Check the uvfits prepend to make sure we end in .uvfits
    output_uvfits_prepend = args.output_uvfits_prepend
    if output_uvfits_prepend[-7:] == '.uvfits': output_uvfits_prepend = output_uvfits_prepend[-7:]

    if args.dry_run:
        pass
    else:

        command(f'{WODEN_DIR}/{woden_exe} {json_name}')

        X,Y,Z = enh2xyz(args.east, args.north, args.height, args.latitude*D2R)
        ##X,Y,Z are stored in a 2D array in units of seconds in the uvfits file
        XYZ_array = np.empty((args.num_antennas,3))
        XYZ_array[:,0] = X
        XYZ_array[:,1] = Y
        XYZ_array[:,2] = Z

        ##Get the central frequency channels, used in the uvfits header
        central_freq_chan = int(np.floor(args.num_freq_channels / 2.0))
        ##Useful number
        num_baselines = int(((args.num_antennas - 1)*args.num_antennas) / 2)

        ##Loop over coarse frequency band
        for band in args.band_nums:
            print('Converting binary to uvfits band',band)

            output_uvfits_name = output_uvfits_prepend + '_band%02d.uvfits' %band

            band_low_freq = args.lowest_channel_freq + (band - 1)*args.coarse_band_width
            central_freq_chan_value = band_low_freq + central_freq_chan*args.freq_res


            filename = "output_visi_band{:02d}.dat".format(band)
            uus,vvs,wws,v_container = load_data(filename=filename,num_baselines=num_baselines,
                                                num_freq_channels=args.num_freq_channels,
                                                num_time_steps=args.num_time_steps,
                                                precision=args.precision)


            if args.remove_phase_tracking:
                frequencies = band_low_freq + np.arange(args.num_freq_channels)*args.freq_res

                v_container = remove_phase_tracking(frequencies=frequencies,
                                          wws_seconds=wws,
                                          num_time_steps=args.num_time_steps,
                                          v_container=v_container,
                                          num_baselines=num_baselines)

            hdu_ant = make_antenna_table(XYZ_array=XYZ_array,telescope_name=args.telescope_name,
                          num_antennas=args.num_antennas, freq_cent=central_freq_chan_value,
                          date=args.date, gst0_deg=gst0_deg, degpdy=degpdy,
                          ut1utc=ut1utc, longitude=args.longitude, latitude=args.latitude,
                          array_height=args.array_height)

            baselines_array, date_array = make_baseline_date_arrays(args.num_antennas,
                                          args.date, args.num_time_steps, args.time_res)

            create_uvfits(v_container=v_container, freq_cent=central_freq_chan_value,
                          ra_point=args.ra0, dec_point=args.dec0,
                          output_uvfits_name=output_uvfits_name,
                          uu=uus, vv=vvs, ww=wws, baselines_array=baselines_array,
                          date_array=date_array,
                          central_freq_chan=central_freq_chan,
                          ch_width=args.freq_res, int_jd=int_jd,
                          hdu_ant=hdu_ant, gitlabel=gitlabel,
                          longitude=args.longitude, latitude=args.latitude,
                          array_height=args.array_height,
                          telescope_name=args.telescope_name,)


            ##Tidy up or not
            if args.no_tidy:
                pass
            else:
                command("rm {:s}".format(filename))
                # if args.array_layout == 'from_the_metafits':
                #     command("rm WODEN_array_layout_band{:d}.txt".format(band))

        ##Tidy up or not
        if args.no_tidy:
            pass
        else:
            command("rm {:s}".format(json_name))
            ##if we generated a text file containing the array layout
            ##from the metafits, delete it now
            # if args.array_layout == 'from_the_metafits':
            command("rm WODEN_array_layout_band{:s}.txt".format(json_band_str))
            #     command("rm {:s}".format("WODEN_array_layout.txt"))
