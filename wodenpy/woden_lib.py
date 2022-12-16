from __future__ import print_function
from copy import deepcopy
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
import ctypes
import importlib_resources
import wodenpy

##Constants
R2D = 180.0 / np.pi
D2R = np.pi / 180.0
# MWA_LAT = -26.703319405555554
# MWA_LONG = 116.67081523611111
# MWA_HEIGHT = 377.827
VELC = 299792458.0
SOLAR2SIDEREAL = 1.00274

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

    lst_type = 'mean'

    observing_time = Time(date, scale='utc', location=observing_location)
    ##Grab the LST
    LST = observing_time.sidereal_time(lst_type)
    LST_deg = LST.value*15.0

    ##uvfits file needs to know the greenwich sidereal time at 0 hours
    ##on the date in question
    zero_date = date.split('T')[0] + "T00:00:00"
    zero_time = Time(zero_date, scale='utc', location=observing_location)
    GST0 = zero_time.sidereal_time(lst_type, 'greenwich')
    GST0_deg = GST0.value*15.0

    ##It also needs to know the rotational rate of the Earth on that day, in
    ##units of degrees per day
    ##Do this by measuring the LST exactly a day later

    date_plus_one_day =  observing_time + TimeDelta(1*u.day)
    LST_plusone = date_plus_one_day.sidereal_time(lst_type)

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
                       ut1utc=None, longitude=None, latitude=None, array_height=None, ant_names=False):
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
    if type(ant_names) == bool:
    
        ant_names = np.array(["%05d" %ant for ant in range(1,num_antennas+1)])
    xlabels = np.array(['X']*num_antennas)
    ylabels = np.array(['Y']*num_antennas)

    ##Make a number of FITS columns to create the antenna table from
    col1 = fits.Column(array=ant_names,name='ANNAME',format='8A')
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
                  int_jd=None, hdu_ant=None, gitlabel=False,
                  IAU_order=False):
    """
    Takes visibility data read in from WODEN binary files, predefined
    BASELINE and DATE arrays and an antenna table, and writes them out
    together into a `uvfits` file. Uses multiple `astropy.io.fits` routines to
    create the final uvfits file. Uses `GroupData`_ and `GroupsHDU`_ to create
    the data table, and then combines with the antenna table in a uvfits
    via `HDUList`_.

    Will only work for data as ordered as coming out of the WODEN C/CUDA code
    (where XX = NS). See `--IAU_order` for more explanation.

    .. _GroupData: https://docs.astropy.org/en/stable/io/fits/api/hdus.html?highlight=GroupsHDU#groupdata
    .. _GroupsHDU: https://docs.astropy.org/en/stable/io/fits/api/hdus.html?highlight=GroupsHDU#groupshdu
    .. _HDUList: https://docs.astropy.org/en/stable/io/fits/api/hdulists.html?highlight=HDUList#hdulist

    Parameters
    ----------
    v_container : float array
        Data container for the visibility data. Should have
        `shape = (num_time_steps*num_baselines, 1, 1, num_freq_channels, 4, 3)`
        and should contain instrumental linear polarisation visibilities. The axes should change as:

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
    IAU_order : Boolean
        By default, the visibilities out of the CUDA/C code have 
        XX = North-South, which is the the IAU ordering. Turns out most people
        want `uvfits` with XX = East-West. So when I`AU_order == True`, do
        not reorder the input data, and add a header value of `IAUORDER` = True.
        If `IAU_order == False`, then the XX is flipped to be East-West by
        reordering the data in 
    """

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

    ##stick a bunch of ones in why not, keeps pyuvdata happy
    subarray = np.ones(len(baselines_array))

    uvparnames = ['UU','VV','WW','DATE','BASELINE', 'ANTENNA1', 'ANTENNA2', 'SUBARRAY']
    parvals = [uu,vv,ww,date_array,baselines_array, antenna1_array, antenna2_array, subarray]

    ##Data out of WODEN C/CUDA code is in IAU pol order, so do nothing
    if IAU_order:
        pass

    ##Swap polarisations from NS-NS, EW-EW, NS-EW, EW-NS
    ##                   to   EW-EW, NS-NS, EW-NS, NS-EW
    else:
        old_data = deepcopy(v_container)

        v_container[:, :, :, :, 0, :] = old_data[:, :, :, :, 1, :]
        v_container[:, :, :, :, 1, :] = old_data[:, :, :, :, 0, :]
        v_container[:, :, :, :, 2, :] = old_data[:, :, :, :, 3, :]
        v_container[:, :, :, :, 3, :] = old_data[:, :, :, :, 2, :]

        ##Get rid of old array, no longer needed
        del old_data

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
    if gitlabel: uvhdu.header['GITLABEL'] = str(gitlabel)

    if IAU_order:
        uvhdu.header['IAUORDER'] = True
    else:
        uvhdu.header['IAUORDER'] = False

    ## Create hdulist and write out file
    hdulist = fits.HDUList(hdus=[uvhdu,hdu_ant])
    hdulist.writeto(output_uvfits_name,overwrite=True)
    hdulist.close()

def enh2xyz(east, north, height, latitude):
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

def load_data(visibility_set=None,num_baselines=None,num_freq_channels=None,
              num_time_steps=None, precision=None, do_autos=False, num_ants=0):
    """
    Read the WODEN binary output and shove into a numpy arrays, ready to be put
    into a uvfits file. By default, WODEN only outputs cross-correlations.
    In this case, the output binary is ordered by baseline (fastest
    changing), frequency, and time (slowest changing). Visibility coords and
    data are read in, with the visi data output into an array of
    `shape=(num_time_steps*num_baselines,1,1,num_freq_channels,4,3))`, which is
    appropriate for a uvfits file. Needs to know whether WODEN was run with
    'float' (32 bit) or 'double' (64 bit) precision to read in the data
    correctly.

    If WODEN was run with `do_autos=True`, then auto correlations are
    also in the output binary. These are stored AFTER the cross-correlations,
    and orders by antenna (fastest changing), frequency, and time (slowest changing).
    In this case the data are output into an array of
    `shape=(num_time_steps*(num_baselines+num_ants),1,1,num_freq_channels,4,3))`,
    where we say a baseline is only defined between to different antennas.
    The visibilities are output to match the BASELINE array, which orders
    the autos and crosses via antenna pairs as (1,1), (1,2), (1,3) .. (2,2),
    (2,3) etc etc meaning the autos and crosses are mixed.

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
    do_autos : Boolean
        if True, data has auto-correlations in
    do_autos : int
        Number of antennas in the array

    Returns
    -------
    uus : float array
        The :math:`u` coordinates (seconds). These are zero for auto-correlations.
    vvs : float array
        The :math:`v` coordinates (seconds). These are zero for auto-correlations.
    wws : float array
        The :math:`w` coordinates (seconds). These are zero for auto-correlations.
    v_container : float array
        Visibility data with
        `shape=(num_time_steps*num_baselines,1,1,num_freq_channels,4,3))`
    """

    ##If not doing autos, ensure this number is zero
    if do_autos == False:
        num_ants = 0

    n_data = num_time_steps * (num_baselines + num_ants)
    uus = np.zeros(n_data)
    vvs = np.zeros(n_data)
    wws = np.zeros(n_data)
    v_container = np.zeros((n_data,1,1,num_freq_channels,4,3))

    num_visi = num_time_steps * num_freq_channels * (num_baselines + num_ants)
    num_cross = num_time_steps * num_freq_channels * num_baselines

    ##Righto, this converts from the ctype POINTER into a numpy array
    ##This is grabbing all the lovely things calculated by the GPU
    us_metres = np.ctypeslib.as_array(visibility_set.us_metres, shape=(num_visi,))
    vs_metres = np.ctypeslib.as_array(visibility_set.vs_metres, shape=(num_visi,))
    ws_metres = np.ctypeslib.as_array(visibility_set.ws_metres, shape=(num_visi,))
    visi_XX_real = np.ctypeslib.as_array(visibility_set.sum_visi_XX_real, shape=(num_visi,))
    visi_XX_imag = np.ctypeslib.as_array(visibility_set.sum_visi_XX_imag, shape=(num_visi,))
    visi_XY_real = np.ctypeslib.as_array(visibility_set.sum_visi_XY_real, shape=(num_visi,))
    visi_XY_imag = np.ctypeslib.as_array(visibility_set.sum_visi_XY_imag, shape=(num_visi,))
    visi_YX_real = np.ctypeslib.as_array(visibility_set.sum_visi_YX_real, shape=(num_visi,))
    visi_YX_imag = np.ctypeslib.as_array(visibility_set.sum_visi_YX_imag, shape=(num_visi,))
    visi_YY_real = np.ctypeslib.as_array(visibility_set.sum_visi_YY_real, shape=(num_visi,))
    visi_YY_imag = np.ctypeslib.as_array(visibility_set.sum_visi_YY_imag, shape=(num_visi,))

    ##If doing auto-correlations, need some mapping arrays so we can
    ##shove the correct data into the correct spots
    if do_autos:
        cross_map = []
        auto_map = []

        visi_map = 0
        for b1 in np.arange(num_ants):
            for b2 in np.arange(b1,num_ants):
                if b1 == b2:
                    auto_map.append(visi_map)
                else:
                    cross_map.append(visi_map)
                visi_map += 1

        cross_map = np.array(cross_map, dtype=int)
        auto_map = np.array(auto_map, dtype=int)

    ##We only care about the u,v,w for the cross-correlations, so fill
    ##them only
    for time_ind in np.arange(num_time_steps):

        time_step = num_baselines * time_ind * num_freq_channels

        ##TODO now that we are reading the u,v,w directly from memory out of
        ##GPU, we don't have to store copies of u,v,w for every frequency like
        ##we did when writing to file. Will be a moderate save on GPU memory
        if do_autos:
            time_base = int(time_ind*(num_baselines + num_ants))
            this_cross_map = cross_map + time_base

            ##The baseline length
            uus[this_cross_map] = us_metres[time_step:time_step+num_baselines] / VELC
            vvs[this_cross_map] = vs_metres[time_step:time_step+num_baselines] / VELC
            wws[this_cross_map] = ws_metres[time_step:time_step+num_baselines] / VELC

        else:

            uus[time_ind*num_baselines:(time_ind + 1)*num_baselines] = us_metres[time_step:time_step+num_baselines] / VELC
            vvs[time_ind*num_baselines:(time_ind + 1)*num_baselines] = vs_metres[time_step:time_step+num_baselines] / VELC
            wws[time_ind*num_baselines:(time_ind + 1)*num_baselines] = ws_metres[time_step:time_step+num_baselines] / VELC

    for time_ind in np.arange(num_time_steps):
        for freq_ind in np.arange(num_freq_channels):

            freq_step = num_baselines * (time_ind * num_freq_channels + freq_ind)

            cross_XX_re = visi_XX_real[freq_step:freq_step+num_baselines]
            cross_XX_im = visi_XX_imag[freq_step:freq_step+num_baselines]
            cross_YY_re = visi_YY_real[freq_step:freq_step+num_baselines]
            cross_YY_im = visi_YY_imag[freq_step:freq_step+num_baselines]
            cross_XY_re = visi_XY_real[freq_step:freq_step+num_baselines]
            cross_XY_im = visi_XY_imag[freq_step:freq_step+num_baselines]
            cross_YX_re = visi_YX_real[freq_step:freq_step+num_baselines]
            cross_YX_im = visi_YX_imag[freq_step:freq_step+num_baselines]

            ##If doing auto-correlations, load up the autos and do some fancy
            ##mapping
            if do_autos:

                time_base = int(time_ind*(num_baselines + num_ants))
                this_cross_map = cross_map + time_base

                v_container[this_cross_map,0,0,freq_ind,0,0] = cross_XX_re
                v_container[this_cross_map,0,0,freq_ind,0,1] = cross_XX_im
                v_container[this_cross_map,0,0,freq_ind,1,0] = cross_YY_re
                v_container[this_cross_map,0,0,freq_ind,1,1] = cross_YY_im
                v_container[this_cross_map,0,0,freq_ind,2,0] = cross_XY_re
                v_container[this_cross_map,0,0,freq_ind,2,1] = cross_XY_im
                v_container[this_cross_map,0,0,freq_ind,3,0] = cross_YX_re
                v_container[this_cross_map,0,0,freq_ind,3,1] = cross_YX_im

                freq_step = num_ants * (time_ind * num_freq_channels + freq_ind)

                real_XX_ind = num_cross + freq_step
                imag_XX_ind = num_cross + freq_step
                real_YY_ind = num_cross + freq_step
                imag_YY_ind = num_cross + freq_step
                real_XY_ind = num_cross + freq_step
                imag_XY_ind = num_cross + freq_step
                real_YX_ind = num_cross + freq_step
                imag_YX_ind = num_cross + freq_step

                auto_XX_re = visi_XX_real[real_XX_ind:real_XX_ind+num_ants]
                auto_XX_im = visi_XX_imag[imag_XX_ind:imag_XX_ind+num_ants]
                auto_YY_re = visi_YY_real[real_YY_ind:real_YY_ind+num_ants]
                auto_YY_im = visi_YY_imag[imag_YY_ind:imag_YY_ind+num_ants]
                auto_XY_re = visi_XY_real[real_XY_ind:real_XY_ind+num_ants]
                auto_XY_im = visi_XY_imag[imag_XY_ind:imag_XY_ind+num_ants]
                auto_YX_re = visi_YX_real[real_YX_ind:real_YX_ind+num_ants]
                auto_YX_im = visi_YX_imag[imag_YX_ind:imag_YX_ind+num_ants]

                this_auto_map = auto_map + time_base

                v_container[this_auto_map,0,0,freq_ind,0,0] = auto_XX_re
                v_container[this_auto_map,0,0,freq_ind,0,1] = auto_XX_im
                v_container[this_auto_map,0,0,freq_ind,1,0] = auto_YY_re
                v_container[this_auto_map,0,0,freq_ind,1,1] = auto_YY_im
                v_container[this_auto_map,0,0,freq_ind,2,0] = auto_XY_re
                v_container[this_auto_map,0,0,freq_ind,2,1] = auto_XY_im
                v_container[this_auto_map,0,0,freq_ind,3,0] = auto_YX_re
                v_container[this_auto_map,0,0,freq_ind,3,1] = auto_YX_im

            ##Otherwise, everything is a cross-correlation so just bung em in
            else:
                v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,0,0] = cross_XX_re
                v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,0,1] = cross_XX_im
                v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,1,0] = cross_YY_re
                v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,1,1] = cross_YY_im
                v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,2,0] = cross_XY_re
                v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,2,1] = cross_XY_im
                v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,3,0] = cross_YX_re
                v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,3,1] = cross_YX_im

    ##Set the weights for everything to one
    v_container[:,0,0,:,:,2] = 1.0

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

def make_baseline_date_arrays(num_antennas, date, num_time_steps, time_res,
                              do_autos=False):
    """Makes the BASELINE and DATE arrays needed in the uvfits file
    The BASELINE array encode which two antennas formed the baseline
    The DATE array contains the fractional jd date, that is added to the
    header value PZERO5, to specify the time each visibility was recorded at.

    In WODEN, by default, only cross-correlations are simulated. To include
    auto-correlations in the BASELINE array, use `do_autos=True`.

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
        time_res : Boolean
            Whether to include auto-correlations (same antenna to antenna)

    Returns
    -------
    baselines_array : float array
        The uvfits antenna encoded array needed for BASELINE
    date_array : float array
        The fractional part of the julian date needed for DATE (days)
    """

    num_baselines = int(((num_antennas - 1)*num_antennas) / 2)

    if do_autos:
        num_baselines += num_antennas

    template_baselines = np.empty(num_baselines)

    ##Loop over all antenna combinations and encode the baseline pair
    baseline_ind = 0

    if do_autos:
        for b1 in np.arange(num_antennas):
            for b2 in np.arange(b1,num_antennas):
                template_baselines[baseline_ind] = RTS_encode_baseline(b1+1, b2+1)
                baseline_ind += 1
    else:
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

# ## define the _components_t struct and it's fields. This thing is inside
# ## simple_C_lib.so
class Visi_Set(ctypes.Structure):
    _fields_ = [("us_metres", ctypes.POINTER(ctypes.c_double)),
                ("vs_metres", ctypes.POINTER(ctypes.c_double)),
                ("ws_metres", ctypes.POINTER(ctypes.c_double)),
                ("allsteps_sha0s", ctypes.POINTER(ctypes.c_double)),
                ("allsteps_cha0s", ctypes.POINTER(ctypes.c_double)),
                ("allsteps_lsts", ctypes.POINTER(ctypes.c_double)),
                ("allsteps_wavelengths", ctypes.POINTER(ctypes.c_double)),
                ("channel_frequencies", ctypes.POINTER(ctypes.c_double)),
                ("sum_visi_XX_real", ctypes.POINTER(ctypes.c_double)),
                ("sum_visi_XX_imag", ctypes.POINTER(ctypes.c_double)),
                ("sum_visi_XY_real", ctypes.POINTER(ctypes.c_double)),
                ("sum_visi_XY_imag", ctypes.POINTER(ctypes.c_double)),
                ("sum_visi_YX_real", ctypes.POINTER(ctypes.c_double)),
                ("sum_visi_YX_imag", ctypes.POINTER(ctypes.c_double)),
                ("sum_visi_YY_real", ctypes.POINTER(ctypes.c_double)),
                ("sum_visi_YY_imag", ctypes.POINTER(ctypes.c_double))]
    
def setup_visi_set(num_visis):
    visibility_set = Visi_Set()
    
    LP_c_double = ctypes.c_double*num_visis
    
    visibility_set.us_metres = LP_c_double()
    visibility_set.vs_metres = LP_c_double()
    visibility_set.ws_metres = LP_c_double()
    visibility_set.sum_visi_XX_real = LP_c_double()
    visibility_set.sum_visi_XX_imag = LP_c_double()
    visibility_set.sum_visi_XY_real = LP_c_double()
    visibility_set.sum_visi_XY_imag = LP_c_double()
    visibility_set.sum_visi_YX_real = LP_c_double()
    visibility_set.sum_visi_YX_imag = LP_c_double()
    visibility_set.sum_visi_YY_real = LP_c_double()
    visibility_set.sum_visi_YY_imag = LP_c_double()
    
    return visibility_set

def load_in_woden_library(precision: str):
    """Load in the WODEN C and CUDA code"""

    woden_lib = importlib_resources.files(wodenpy).joinpath(f"libwoden_{precision}.so")

    ## Read in the C library
    libwoden = ctypes.cdll.LoadLibrary(woden_lib)

    # ##Define the input and return types for the `test_RTS_calculate_MWA_analytic_beam` function
    run_woden = libwoden.run_woden

    run_woden.restype = ctypes.c_int
    run_woden.argtypes = [ctypes.c_char_p, ctypes.POINTER(Visi_Set)]

    return run_woden