from astropy.io import fits
import numpy as np
from wodenpy.observational.calc_obs import calc_jdcal
from erfa import gd2gc
import warnings
import sys
from copy import deepcopy

##Constants
R2D = 180.0 / np.pi
D2R = np.pi / 180.0
SOLAR2SIDEREAL = 1.00274


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
                  IAU_order=False, comment=False):
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
        By default, the visibilities out of the GPU/C code have 
        XX = North-South, which is the the IAU ordering. Turns out most people
        want `uvfits` with XX = East-West. So when `IAU_order == True`, do
        not reorder the input data, and add a header value of `IAUORDER` = True.
        If `IAU_order == False`, then the XX is flipped to be East-West by
        reordering the data in v_container
    comment : string
        Optional COMMENT to add to the uvfits header
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
        
    if comment:
        uvhdu.header['COMMENT'] = comment

    ## Create hdulist and write out file
    hdulist = fits.HDUList(hdus=[uvhdu,hdu_ant])
    hdulist.writeto(output_uvfits_name,overwrite=True)
    hdulist.close()
    
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