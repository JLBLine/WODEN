"""Functions/Classes to create and populate a `Woden_Settings` class structured equivalently
to a `woden_settings_t` struct in the C/C++/GPU code. Created dynamically based on the
required precision."""
import ctypes 
import sys
import os
import subprocess
from typing import Union
from wodenpy.array_layout.precession import RTS_Precess_LST_Lat_to_J2000
from wodenpy.use_libwoden.beam_settings import BeamTypes, BeamGroups
import numpy as np
import argparse
from logging import Logger
from wodenpy.wodenpy_setup.woden_logger import simple_logger
from ctypes import c_double, c_float, c_int, POINTER, c_ulong, c_char, create_string_buffer, c_uint, c_long, c_uint64

D2R = np.pi/180.0
SOLAR2SIDEREAL = 1.00274
DS2R  = 7.2722052166430399038487115353692196393452995355905e-5
DD2R = 0.017453292519943295769236907684886127134428718885417

def command(cmd):
    """
    Runs the command string `cmd` using `subprocess.call`

    Parameters
    ----------
    cmd : string
         The command to run on the command line
    """
    subprocess.call(cmd,shell=True)
    
class Woden_Settings_Python(object):
    """A class structured equivalently to a `woden_settings_t` struct, used by
    the C/C++/GPU code in libwoden_float.so or libwoden_double.so. Retain a
    copy of the settings in Python for easy access and manipulation."""
    def __init__(self):
        self.lst_base = None
        self.lst_obs_epoch_base = None
        self.ra0 = None
        self.dec0 = None
        self.sdec0 = None
        self.cdec0 = None
        self.num_baselines = None
        self.num_ants = None
        self.num_freqs = None
        self.frequency_resolution = None
        self.base_low_freq = None
        self.num_time_steps = None
        self.time_res = None
        self.cat_filename = None
        self.num_bands = None
        self.band_nums = None
        self.sky_crop_type = None
        self.beamtype = None
        self.gauss_beam_FWHM = None
        self.gauss_beam_ref_freq = None
        self.chunking_size = None
        self.hdf5_beam_path = None
        self.jd_date = None
        self.array_layout_file = None
        self.array_layout_file_path = None
        self.latitude = None
        self.latitude_obs_epoch_base = None
        self.longitude = None
        self.FEE_ideal_delays = None
        self.coarse_band_width = None
        self.gauss_ra_point = None
        self.gauss_dec_point = None
        self.num_cross = None
        self.num_autos = None
        self.num_visis = None
        self.base_band_freq = None
        self.do_precession = None
        self.lsts = None
        self.latitudes = None
        self.mjds = None
        self.do_autos = None
        self.use_dipamps = None
        self.mwa_dipole_amps = None
        self.single_everybeam_station = None
        self.off_cardinal_dipoles = None
        self.do_gpu = None
        self.normalise_primary_beam = None

def create_woden_settings_struct(precision : str = "double"):
    """Creates a `Woden_Settings` class structured equivalently to a `woden_settings_t`
    struct in the C/CUDA code. Created dynamically based on the `precision`,
    to match the compile time precision flag `-DUSE_DOUBLE` in the C code.

    Parameters
    ----------
    precision : str, optional
        Either "float" or "double:, by default "double"

    Returns
    -------
    Woden_Settings
        The Woden_Settings class structured equivalently to a `woden_settings_t` struct
    """
    
    if precision == "float":
        c_user_precision = c_float
    else:
        c_user_precision = c_double

    ##TODO gotta be a way to set the float or double fields via some kind of
    ##variable instead of making two different classes
    class Woden_Settings(ctypes.Structure):
        """A class structured equivalently to a `woden_settings_t` struct, used by 
        the C and CUDA code in libwoden_float.so or libwoden_double.so.
        
        Created by the function `create_woden_settings_struct`, which sets
        `user_precision_t` to either `c_float` or `c_double`.
        
        :cvar c_double lst_base:  Local sidereal time for first time step (radians)
        :cvar c_double lst_obs_epoch_base:  Local sidereal time for first time step (radians) for the observation epoch (e.g. in 2020 for a 2020 obs)
        :cvar c_double ra0:  Right ascension of phase centre (radians)
        :cvar c_double dec0:  Declination of phase centre (radians)
        :cvar c_double sdec0:  Sine of Declination of phase centre (radians)
        :cvar c_double cdec0:  Cosine of Declination of phase centre (radians)
        :cvar c_int num_baselines:  Number of baselines this array layout has
        :cvar c_int num_ants:  Number of antennas this array layout has (MWA calls this number of tiles)
        :cvar c_int num_freqs:  Number of frequencies per coarse band
        :cvar c_double frequency_resolution:  Frequency resolution of a fine channel (Hz)
        :cvar c_double base_low_freq:  The lowest fine channel frequency of band 1
        :cvar c_int num_time_steps:  Number of time steps to simulate
        :cvar c_double time_res:  Time resolution of simulation (seconds)
        :cvar c_int num_bands:  Number of coarse frequency bands to simulate 
        :cvar POINTER(c_int) band_nums:  Which number coarse bands to simulate (e.g 1,4,6) 
        :cvar c_int sky_crop_type:  Whether to crop sky models by SOURCE or COMPONENT 
        :cvar c_int beamtype:  What type of primary beam to simulate with 
        :cvar c_user_precision gauss_beam_FWHM:  FWHM of Gaussian primary beam (degrees)
        :cvar c_double gauss_beam_ref_freq:  Reference frequency for given Gaussian primary beam FWHM
        :cvar c_uint64 chunking_size:  Maximum number of COMPONENTs to include in a single chunk
        :cvar POINTER(c_char) hdf5_beam_path:  Path to *.hf file containing MWA FEE beam spherical harmonic information
        :cvar c_double jd_date:  Julian date at beginning of simulation
        :cvar c_int array_layout_file:  Do we have a path to the array layout or not 
        :cvar POINTER(c_char) array_layout_file_path:  ath to file containing E,N,H coords of array layout 
        :cvar c_double latitude:  Latitude of the array to simulate (radians) 
        :cvar c_double latitude_obs_epoch_base:  Latitude of the array at the observation epoch (radians) 
        :cvar c_user_precision longitude:  Longitude of the array to simulate (radians) 
        :cvar POINTER(c_int) FEE_ideal_delays:  Delay values specifying the pointing for the MWA FEE beam model 
        :cvar c_double coarse_band_width:  Frequency bandwidth of a single coarse band (Hz)
        :cvar c_double gauss_ra_point:  The initial Right Ascension to point the Gaussian beam at (radians)
        :cvar c_double gauss_dec_point:  The initial Declination to point the Gaussian beam at (radians)
        :cvar c_int num_cross:  Total number of cross-correlations to simulate, so freqs*times*baselines 
        :cvar c_int num_autos:  Total number of auto-correlations to simulate, so freqs*times*baselines 
        :cvar c_int num_visis:  Total number of visibilities to simulate, so num_cross + num_autos 
        :cvar c_double base_band_freq:  The lowest fine channel frequency in the current band being simulated
        :cvar c_int do_precession:  Boolean of whether to apply precession to the array layout or not
        :cvar POINTER(c_double) lsts:  Array to hold LSTs for all time centroids (these are different when precession is happening)
        :cvar POINTER(c_double) latitudes:  Array to hold latitudes for all time centroids (these are different when precession is happening)
        :cvar POINTER(c_double) mjds:  Array to hold modified julian dates for all time centroids
        :cvar c_int do_autos:  Boolean of whether to simulate autos or not (0 False, 1 True)
        :cvar c_int use_dipamps:  Boolean of whether to use dipole amplitudes, so have an individual beam per tile  (0 False, 1 True)
        :cvar POINTER(c_double) mwa_dipole_amps: Bespoke MWA dipole amplitudes for each antenna(tile). Should be 2*num_ants*16 long
        :cvar c_int single_everybeam_station: If using everybeam, add this to say we are only using a single station
        :cvar c_int off_cardinal_dipoles: Boolean of whether to use off-cardinal dipole equations to apply the beams gains to the Stokes IQUV parameters
        :cvar c_int use_gpu: Boolean of whether to use the GPU or not (0 False, 1 True)
        :cvar c_int verbose: Boolean of whether to do verbose logging or not (0 False, 1 True)
        :cvar c_int normalise_primary_beam: Boolean of whether to normalise the primary beam (0 False, 1 True)
        :cvar POINTER(c_char) beam_ms_path:  Path to the beam model MS file for everybeam simulations
        :cvar c_double eb_beam_ra0:  Right ascension to lock the EveryBeam primary beam centre to (radians)
        :cvar c_double eb_beam_dec0:  Declination to lock the EveryBeam primary beam centre to (radians)
        """
        
        _fields_ = [("lst_base", c_double),
                    ("lst_obs_epoch_base", c_double),
                    ("ra0", c_double),
                    ("dec0", c_double),
                    ("sdec0", c_double),
                    ("cdec0", c_double),
                    ("num_baselines", c_int),
                    ("num_ants", c_int),
                    ("num_freqs", c_int),
                    ("frequency_resolution", c_double),
                    ("base_low_freq", c_double),
                    ("num_time_steps", c_int),
                    ("time_res", c_double),
                    ("num_bands", c_int),
                    ("band_nums", POINTER(c_int)),
                    ("sky_crop_type", c_int),
                    ("beamtype", c_int),
                    ("gauss_beam_FWHM", c_user_precision),
                    ("gauss_beam_ref_freq", c_double),
                    ("chunking_size", c_ulong),
                    ("hdf5_beam_path", POINTER(c_char)),
                    ("jd_date", c_double),
                    ("array_layout_file", c_int),
                    ("array_layout_file_path", POINTER(c_char)),
                    ("latitude", c_double),
                    ("latitude_obs_epoch_base", c_double),
                    ("longitude", c_user_precision),
                    ("FEE_ideal_delays", POINTER(c_int)),
                    ("coarse_band_width", c_double),
                    ("gauss_ra_point", c_double),
                    ("gauss_dec_point", c_double),
                    ("num_cross", c_int),
                    ("num_autos", c_int),
                    ("num_visis", c_int),
                    ("base_band_freq", c_double),
                    ("do_precession", c_int),
                    ("lsts", POINTER(c_double)),
                    ("latitudes", POINTER(c_double)),
                    ("mjds", POINTER(c_double)),
                    ("do_autos", c_int),
                    ("use_dipamps", c_int),
                    ("mwa_dipole_amps", POINTER(c_double)),
                    ("single_everybeam_station", c_int),
                    ("off_cardinal_dipoles", c_int),
                    ("do_gpu", c_int),
                    ("verbose", c_int),
                    ("normalise_primary_beam", c_int),
                    ("beam_ms_path", POINTER(c_char)),
                    ("eb_beam_ra0", c_double),
                    ("eb_beam_dec0", c_double),]
        
    return Woden_Settings

##This call is so we can use it as a type annotation, and so sphinx can document the class
Woden_Settings = create_woden_settings_struct()
    
def fill_woden_settings_python(args : argparse.Namespace,
                               jd_date : float, lst : float):
    """Given the parsed and checked arguments in `args`, populate a
    `woden_settings` ctypes.Structure that can be passed into
    libwoden_float.so or libwoden_double.so, depending on the desired
    `precision`.

    Parameters
    ----------
    woden_settings : Woden_Settings
        An initialised Woden_Settings class struct
    args : argparse.Namespace
        The populated arguments `args = parser.parse_args()`` as returned from
        the parser given by :func:`~wodenpy.run_setup.check_args`
    jd_date : float
        The julian date of the first time step of the observation
    lst : float
        The LST of the first time step of the observation (degrees)

    Returns
    -------
    woden_settings : Woden_Settings_Python
       Populated Woden_Settings_Python class
    """
    woden_settings = Woden_Settings_Python()
    
    woden_settings.ra0 = args.ra0 * D2R
    woden_settings.dec0 = args.dec0 * D2R
    
    woden_settings.latitude = args.latitude * D2R
    woden_settings.latitude_obs_epoch_base = args.latitude * D2R
    
    woden_settings.lst_base = np.float64(lst * D2R)
    woden_settings.lst_obs_epoch_base = lst * D2R
    
    woden_settings.num_freqs = args.num_freq_channels
    woden_settings.frequency_resolution = args.freq_res
    
    woden_settings.base_low_freq = args.lowest_channel_freq
    
    woden_settings.coarse_band_width = float(args.coarse_band_width)
    
    woden_settings.num_time_steps = args.num_time_steps
    woden_settings.time_res = args.time_res
    
    woden_settings.jd_date = jd_date
    
    woden_settings.sky_crop_type = int(args.sky_crop_components)
    
    ##If MWA_FEE_delays is set, convert into a array and populate the
    ##woden_settings equivalent
    if args.MWA_FEE_delays:
        
        ##If using a a different set of dipole amplitudes for each tile,
        ##need to repeat the delays for each tile for hyperdrive
        if args.use_MWA_dipamps:
            num_beams = args.num_antennas
        else:
            num_beams = 1
            
        delays = np.array(args.MWA_FEE_delays.strip('[]').split(','))
        num_delays = len(delays)*num_beams
        woden_settings.FEE_ideal_delays = np.empty(num_delays, dtype=np.int32)
        
        for delay_ind, delay in enumerate(delays):
            for beam in range(num_beams):
                woden_settings.FEE_ideal_delays[beam*len(delays) + delay_ind] = int(delay)
        
    if args.primary_beam == 'none':
        woden_settings.beamtype = BeamTypes.NO_BEAM.value
    
    elif args.primary_beam == 'Gaussian':
        woden_settings.beamtype = BeamTypes.GAUSS_BEAM.value
        woden_settings.gauss_beam_FWHM = float(args.gauss_beam_FWHM)
        woden_settings.gauss_beam_ref_freq = float(args.gauss_beam_ref_freq)
        woden_settings.gauss_ra_point = float(args.gauss_ra_point)*D2R
        woden_settings.gauss_dec_point = float(args.gauss_dec_point)*D2R

    elif args.primary_beam == 'MWA_FEE':
        woden_settings.beamtype = BeamTypes.FEE_BEAM.value
        woden_settings.hdf5_beam_path = args.hdf5_beam_path
        
    elif args.primary_beam == 'EDA2':
        woden_settings.beamtype = BeamTypes.ANALY_DIPOLE.value
        
    elif args.primary_beam == 'MWA_FEE_interp':
        woden_settings.beamtype = BeamTypes.FEE_BEAM_INTERP.value
        woden_settings.hdf5_beam_path = args.hdf5_beam_path
        
    elif args.primary_beam == 'MWA_analy':
        woden_settings.beamtype = BeamTypes.MWA_ANALY.value
        
    elif args.primary_beam == 'everybeam_OSKAR':
        woden_settings.beamtype = BeamTypes.EB_OSKAR.value
        
    elif args.primary_beam == 'everybeam_LOFAR':
        woden_settings.beamtype = BeamTypes.EB_LOFAR.value
        
    elif args.primary_beam == 'everybeam_MWA':
        woden_settings.hdf5_beam_path = args.hdf5_beam_path
        woden_settings.beamtype = BeamTypes.EB_MWA.value
        
    elif args.primary_beam == 'uvbeam_MWA':
        woden_settings.beamtype = BeamTypes.UVB_MWA.value
        woden_settings.hdf5_beam_path = args.hdf5_beam_path
        
    elif args.primary_beam == 'uvbeam_HERA':
        woden_settings.beamtype = BeamTypes.UVB_HERA.value
        
    if woden_settings.beamtype in BeamGroups.eb_beam_values:
        woden_settings.eb_beam_ra0 = float(args.eb_ra_point)*D2R
        woden_settings.eb_beam_dec0 = float(args.eb_dec_point)*D2R
        
    if args.no_precession:
        woden_settings.do_precession = 0
    else:
        woden_settings.do_precession = 1
        
    woden_settings.do_autos = args.do_autos
    
    woden_settings.chunking_size = int(args.chunking_size)
    
    woden_settings.num_bands = len(args.band_nums)
    woden_settings.band_nums = np.array(args.band_nums, dtype=np.int32)
    
    ##Are we using dipole amplitudes?
    woden_settings.use_dipamps = args.use_MWA_dipamps
    
    ##If so, populate the array
    if args.use_MWA_dipamps:
        woden_settings.mwa_dipole_amps = args.dipamps.astype(np.float64)
    else:
        woden_settings.mwa_dipole_amps = np.ones(16, dtype=np.float64)
    
    ##Always set identical primary beam for EveryBeam MWA
    ##Need to implement mapping of amplitudes to the c++ calls.
    ##The everybeam c++ code maps the same dipole amplitudes to both X and Y
    ##pols so I think you'd need to run everything twice and be careful mapping
    ##the correct pol to output jones array index
    if args.primary_beam == 'everybeam_MWA':
        woden_settings.single_everybeam_station = 1
        woden_settings.hdf5_beam_path = args.hdf5_beam_path
        woden_settings.mwa_dipole_amps = np.ones(16, dtype=np.float64)
    elif np.isnan(args.station_id):
        woden_settings.single_everybeam_station = 0
    else:
        woden_settings.single_everybeam_station = 1
        
    if args.off_cardinal_dipoles or woden_settings.beamtype in BeamGroups.off_cardinal_beam_values:
        woden_settings.off_cardinal_dipoles = 1
    else:
        woden_settings.off_cardinal_dipoles = 0
        
    if args.cpu_mode:
        woden_settings.do_gpu = 0
    else:
        woden_settings.do_gpu = 1
        
    woden_settings.verbose = args.verbose
    
    if args.no_beam_normalisation:
        woden_settings.normalise_primary_beam = 0
    else:
        woden_settings.normalise_primary_beam = 1
        
    if args.pointed_ms_file_name:
        woden_settings.beam_ms_path = args.pointed_ms_file_name.as_posix()
    else:
        woden_settings.beam_ms_path = args.beam_ms_path
    
    return woden_settings
    
def setup_lsts_and_phase_centre(woden_settings_python : Woden_Settings_Python,
                                logger : Logger = False) -> Union[np.ndarray, np.ndarray]: # type: ignore
    """
    Calculate the Local Sidereal Time (LST) for each time step of an observation,
    and set the phase centre coordinates. If `woden_settings.do_precession == True`,
    the returned LSTs are precessed back to J2000.

    Parameters
    ----------
    woden_settings_python : Woden_Settings_Python
        A populated Woden_Settings object containing the observation parameters.
    logger : Logger, optional
        A logger object for logging messages, by default False. If False, a
        simple logger is created.

    Returns
    -------
    lsts : np.ndarray:
        An array of LST values, one for each time step of the observation.
    latitudes : np.ndarray:
        An array of latitude values, one for each time step of the observation
        (these should be the same if precession is not applied, and even if applied
        there should be tiny differences).

    """
    
    if not logger:
        logger = simple_logger()

    ##Used for calculating l,m,n for components
    woden_settings_python.sdec0 = np.sin(woden_settings_python.dec0)
    woden_settings_python.cdec0 = np.cos(woden_settings_python.dec0)

    # logger.info("Setting phase centre RA,DEC {:.5f}deg {:.5f}deg".format(woden_settings_python.ra0/DD2R, woden_settings_python.dec0/DD2R))
    
    ##Calculate all lsts for this observation
    ##Used in some python calcs later, and by the C code, so store a ctypes
    ##array as well as a numpy one
    lsts = np.empty(woden_settings_python.num_time_steps)

    # num_time_array = woden_settings_python.num_time_steps*c_double


    woden_settings_python.lsts = np.empty(woden_settings_python.num_time_steps, dtype=np.float64)
    woden_settings_python.latitudes = np.empty(woden_settings_python.num_time_steps, dtype=np.float64)
    woden_settings_python.mjds = np.empty(woden_settings_python.num_time_steps, dtype=np.float64)
    
    mjd = woden_settings_python.jd_date - 2400000.5
    
    for time_step in range(woden_settings_python.num_time_steps):
        
        ##Add on the angle accrued by current time step to the base LST
        lst_current = woden_settings_python.lst_obs_epoch_base + time_step*woden_settings_python.time_res*SOLAR2SIDEREAL*DS2R

        ##Add half a time_res so we are sampling centre of each time step
        lst_current += 0.5*woden_settings_python.time_res*SOLAR2SIDEREAL*DS2R

        if (woden_settings_python.do_precession):
            #
            ##Move the mjd to the time of the current step
            ##Want the LST at centre of time step so 0.5 adds half a time step extra
            mjd_current = mjd + ((time_step + 0.5)*woden_settings_python.time_res)/(24.0*60.0*60.0)
            woden_settings_python.mjds[time_step] = mjd_current
            
            lst_J2000, latitude_J2000 = RTS_Precess_LST_Lat_to_J2000(
                                 lst_current,
                                 woden_settings_python.latitude_obs_epoch_base,
                                 mjd_current)

            lsts[time_step] = lst_J2000
            woden_settings_python.lsts[time_step] = lst_J2000
            woden_settings_python.latitudes[time_step] = latitude_J2000

            if (time_step == 0):
                logger.info("Obs epoch initial LST was {:.10f} deg".format(lst_current/D2R) )
                logger.info("Setting initial J2000 LST to {:.10f} deg".format(lst_J2000/D2R) )
                logger.info("Setting initial mjd to {:.10f}".format(woden_settings_python.mjds[time_step]) )
                logger.info("After precession initial latitude of the array is {:.10f} deg".format(latitude_J2000/D2R) )
                woden_settings_python.lst_base = lst_J2000
                woden_settings_python.latitude = latitude_J2000
    
        else:
            
            lsts[time_step] = lst_current
            woden_settings_python.lsts[time_step] = lst_current
            woden_settings_python.latitudes[time_step] = woden_settings_python.latitude_obs_epoch_base
            if time_step == 0:
                logger.info("Obs epoch initial LST was {:.10f} deg".format(lst_current/D2R) )
        
    latitudes = np.ctypeslib.as_array(woden_settings_python.latitudes, shape=(woden_settings_python.num_time_steps, ))
    
    return lsts, latitudes

def convert_woden_settings_to_ctypes(woden_settings_python : Woden_Settings_Python,
                                     woden_settings_ctypes: Woden_Settings) -> Woden_Settings: #type: ignore
    """
    Converts Woden settings from a Python object to a ctypes object.
    Parameters
    ----------
    woden_settings_python : Woden_Settings_Python
        The Woden settings in Python object format.
    woden_settings_ctypes : Woden_Settings
        The Woden settings in ctypes object format to be populated.
    Returns
    -------
    Woden_Settings
        The populated Woden settings in ctypes object format.
    """
    
    woden_settings_ctypes.lst_base = woden_settings_python.lst_base
    woden_settings_ctypes.lst_obs_epoch_base = woden_settings_python.lst_obs_epoch_base
        
    woden_settings_ctypes.ra0 = woden_settings_python.ra0
    woden_settings_ctypes.dec0 = woden_settings_python.dec0
    woden_settings_ctypes.sdec0 = woden_settings_python.sdec0
    woden_settings_ctypes.cdec0 = woden_settings_python.cdec0
    woden_settings_ctypes.num_baselines = woden_settings_python.num_baselines
    woden_settings_ctypes.num_ants = woden_settings_python.num_ants
    woden_settings_ctypes.num_freqs = woden_settings_python.num_freqs
    woden_settings_ctypes.frequency_resolution = woden_settings_python.frequency_resolution
    woden_settings_ctypes.base_low_freq = woden_settings_python.base_low_freq
    woden_settings_ctypes.num_time_steps = woden_settings_python.num_time_steps
    woden_settings_ctypes.time_res = woden_settings_python.time_res
    woden_settings_ctypes.num_bands = woden_settings_python.num_bands
    woden_settings_ctypes.beamtype = woden_settings_python.beamtype
    
    if woden_settings_python.beamtype == BeamTypes.GAUSS_BEAM.value:
        woden_settings_ctypes.gauss_beam_FWHM = woden_settings_python.gauss_beam_FWHM
        woden_settings_ctypes.gauss_beam_ref_freq = woden_settings_python.gauss_beam_ref_freq
        woden_settings_ctypes.gauss_ra_point = woden_settings_python.gauss_ra_point
        woden_settings_ctypes.gauss_dec_point = woden_settings_python.gauss_dec_point
        
    woden_settings_ctypes.chunking_size = woden_settings_python.chunking_size
    woden_settings_ctypes.jd_date = woden_settings_python.jd_date
    woden_settings_ctypes.latitude = woden_settings_python.latitude
    woden_settings_ctypes.latitude_obs_epoch_base = woden_settings_python.latitude_obs_epoch_base
    # woden_settings_ctypes.longitude = woden_settings_python.longitude
    woden_settings_ctypes.coarse_band_width = woden_settings_python.coarse_band_width

    woden_settings_ctypes.num_cross = woden_settings_python.num_cross
    woden_settings_ctypes.num_autos = woden_settings_python.num_autos
    woden_settings_ctypes.num_visis = woden_settings_python.num_visis
    # woden_settings_ctypes.base_band_freq = woden_settings_python.base_band_freq
    woden_settings_ctypes.do_precession = woden_settings_python.do_precession
    woden_settings_ctypes.do_autos = woden_settings_python.do_autos
    woden_settings_ctypes.use_dipamps = woden_settings_python.use_dipamps
    woden_settings_ctypes.single_everybeam_station = woden_settings_python.single_everybeam_station
    woden_settings_ctypes.off_cardinal_dipoles = woden_settings_python.off_cardinal_dipoles
    woden_settings_ctypes.do_gpu = woden_settings_python.do_gpu
    
    woden_settings_ctypes.mjds = woden_settings_python.mjds.ctypes.data_as(POINTER(c_double))
    woden_settings_ctypes.latitudes = woden_settings_python.latitudes.ctypes.data_as(POINTER(c_double))
    woden_settings_ctypes.lsts = woden_settings_python.lsts.ctypes.data_as(POINTER(c_double))
    woden_settings_ctypes.band_nums = woden_settings_python.band_nums.ctypes.data_as(POINTER(c_int))
    
    if woden_settings_ctypes.beamtype in BeamGroups.needs_MWA_delays:
        woden_settings_ctypes.FEE_ideal_delays = woden_settings_python.FEE_ideal_delays.ctypes.data_as(POINTER(c_int))
        
    if woden_settings_ctypes.beamtype in BeamGroups.needs_MWA_hdf5_path:
        woden_settings_ctypes.hdf5_beam_path = create_string_buffer(woden_settings_python.hdf5_beam_path.encode('utf-8'))
        
    if woden_settings_ctypes.beamtype in BeamGroups.eb_beam_values:
        woden_settings_ctypes.beam_ms_path = create_string_buffer(woden_settings_python.beam_ms_path.encode('utf-8'))
        woden_settings_ctypes.eb_beam_ra0 = woden_settings_python.eb_beam_ra0
        woden_settings_ctypes.eb_beam_dec0 = woden_settings_python.eb_beam_dec0
    
    # if woden_settings_ctypes.use_dipamps:
    woden_settings_ctypes.mwa_dipole_amps = woden_settings_python.mwa_dipole_amps.ctypes.data_as(POINTER(c_double))
    
    woden_settings_ctypes.verbose = woden_settings_python.verbose
    woden_settings_ctypes.normalise_primary_beam = woden_settings_python.normalise_primary_beam
    
    return woden_settings_ctypes