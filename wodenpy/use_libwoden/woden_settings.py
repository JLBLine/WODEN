import ctypes 
import sys
import os
import subprocess
from typing import Union
from wodenpy.array_layout.precession import RTS_Precess_LST_Lat_to_J2000
from wodenpy.use_libwoden.beam_settings import BeamTypes
import numpy as np
import argparse

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
        :cvar POINTER(c_char) cat_filename:  Path to WODEN-style sky model
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
                    ("cat_filename", POINTER(c_char)),
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
                    ("mwa_dipole_amps", POINTER(c_double))]
        
    return Woden_Settings

##This call is so we can use it as a type annotation, and so sphinx can document the class
Woden_Settings = create_woden_settings_struct()
    
def create_woden_settings(woden_settings : Woden_Settings, # type: ignore
                          args : argparse.Namespace,
                          jd_date : float, lst : float) -> Woden_Settings:  ## type: ignore
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
    woden_settings : Woden_Settings
       Populated ctype struct that can be passed into the C/CUDA code
    """
    
    woden_settings.ra0 = args.ra0 * D2R
    woden_settings.dec0 = args.dec0 * D2R
    
    woden_settings.latitude = args.latitude * D2R
    woden_settings.latitude_obs_epoch_base = args.latitude * D2R
    
    woden_settings.lst_base = lst * D2R
    woden_settings.lst_obs_epoch_base = lst * D2R
    
    woden_settings.num_freqs = args.num_freq_channels
    woden_settings.frequency_resolution = args.freq_res
    
    woden_settings.base_low_freq = args.lowest_channel_freq
    
    woden_settings.coarse_band_width = float(args.coarse_band_width)
    
    woden_settings.num_time_steps = args.num_time_steps
    woden_settings.time_res = args.time_res
    
    woden_settings.cat_filename = create_string_buffer(args.cat_filename.encode('utf-8'))
    
    woden_settings.jd_date = jd_date
    
    woden_settings.sky_crop_type = int(args.sky_crop_components)
    
    ##If MWA_FEE_delays is set, convert into a array and populate the
    ##woden_settings equivalent
    if args.MWA_FEE_delays:
        
        ##If using a a different set of dipole amplitudes for each tile,
        ##need to repeat the delays for each tile
        if args.use_MWA_dipamps:
            num_beams = args.num_antennas
        else:
            num_beams = 1
        
        delays = np.array(args.MWA_FEE_delays.strip('[]').split(','))
        num_delays = len(delays)*num_beams
        woden_settings.FEE_ideal_delays = (ctypes.c_int*num_delays)()
        
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
        woden_settings.hdf5_beam_path = create_string_buffer(args.hdf5_beam_path.encode('utf-8'))
        
    elif args.primary_beam == 'EDA2':
        woden_settings.beamtype = BeamTypes.ANALY_DIPOLE.value
        
    elif args.primary_beam == 'MWA_FEE_interp':
        woden_settings.beamtype = BeamTypes.FEE_BEAM_INTERP.value
        woden_settings.hdf5_beam_path = create_string_buffer(args.hdf5_beam_path.encode('utf-8'))
        
    elif args.primary_beam == 'MWA_analy':
        woden_settings.beamtype = BeamTypes.MWA_ANALY.value
        
    elif args.primary_beam == 'everybeam_OSKAR':
        woden_settings.beamtype = BeamTypes.EB_OSKAR.value
        
    if args.no_precession:
        woden_settings.do_precession = 0
    else:
        woden_settings.do_precession = 1
        
    woden_settings.do_autos = args.do_autos
    
    woden_settings.chunking_size = int(args.chunking_size)
    
    woden_settings.num_bands = len(args.band_nums)
    woden_settings.band_nums = (ctypes.c_int*woden_settings.num_bands)()
    
    for ind, band in enumerate(args.band_nums):
        woden_settings.band_nums[ind] = int(band)
    
    ##Are we using dipole amplitudes?
    woden_settings.use_dipamps = args.use_MWA_dipamps
    
    ##If so, populate the array
    if args.use_MWA_dipamps:
        ##Assign le-memory
        woden_settings.mwa_dipole_amps = (ctypes.c_double*len(args.dipamps))()
        ##again, this is ctypes, so we need to populate it with a loop
        for ind, amplitude in enumerate(args.dipamps):
                woden_settings.mwa_dipole_amps[ind] = amplitude
    
    return woden_settings
    
def setup_lsts_and_phase_centre(woden_settings : Woden_Settings) -> np.ndarray: # type: ignore
    """
    Calculate the Local Sidereal Time (LST) for each time step of an observation,
    and set the phase centre coordinates. If `woden_settings.do_precession == True`,
    the returned LSTs are precessed back to J2000.

    Parameters
    ----------
    woden_settings : Woden_Settings
        A populated Woden_Settings object containing the observation parameters.

    Returns
    -------
    lsts : np.ndarray:
        An array of LST values, one for each time step of the observation.

    """

    ##Used for calculating l,m,n for components
    woden_settings.sdec0 = np.sin(woden_settings.dec0)
    woden_settings.cdec0 = np.cos(woden_settings.dec0)

    print("Setting phase centre RA,DEC {:.5f}deg {:.5f}deg".format(woden_settings.ra0/DD2R, woden_settings.dec0/DD2R))

    
    ##Calculate all lsts for this observation
    ##Used in some python calcs later, and by the C code, so store a ctypes
    ##array as well as a numpy one
    lsts = np.empty(woden_settings.num_time_steps)

    num_time_array = woden_settings.num_time_steps*c_double


    woden_settings.lsts = num_time_array()
    woden_settings.latitudes = num_time_array()
    woden_settings.mjds = num_time_array()
    
    mjd = woden_settings.jd_date - 2400000.5

    for time_step in range(woden_settings.num_time_steps):
        
        ##Add on the angle accrued by current time step to the base LST
        lst_current = woden_settings.lst_obs_epoch_base + time_step*woden_settings.time_res*SOLAR2SIDEREAL*DS2R

        ##Add half a time_res so we are sampling centre of each time step
        lst_current += 0.5*woden_settings.time_res*SOLAR2SIDEREAL*DS2R

        if (woden_settings.do_precession):
            #
            ##Move the mjd to the time of the current step
            ##Want the LST at centre of time step so 0.5 adds half a time step extra
            mjd_current = mjd + ((time_step + 0.5)*woden_settings.time_res)/(24.0*60.0*60.0)
            woden_settings.mjds[time_step] = mjd_current
            
            lst_J2000, latitude_J2000 = RTS_Precess_LST_Lat_to_J2000(
                                 lst_current,
                                 woden_settings.latitude_obs_epoch_base,
                                 mjd_current)

            lsts[time_step] = lst_J2000
            woden_settings.lsts[time_step] = lst_J2000
            woden_settings.latitudes[time_step] = latitude_J2000

            if (time_step == 0):
                print("Obs epoch initial LST was {:.10f} deg".format(lst_current/D2R) )
                print("Setting initial J2000 LST to {:.10f} deg".format(lst_J2000/D2R) )
                print("Setting initial mjd to {:.10f}".format(woden_settings.mjds[time_step]) )
                print("After precession initial latitude of the array is {:.10f} deg".format(latitude_J2000/D2R) )
                woden_settings.lst_base = lst_J2000
                woden_settings.latitude = latitude_J2000
    
        else:
            
            lsts[time_step] = lst_current
            woden_settings.lsts[time_step] = lst_current
            woden_settings.latitudes[time_step] = woden_settings.latitude_obs_epoch_base
            if time_step == 0:
                print("Obs epoch initial LST was {:.10f} deg\n".format(lst_current/D2R) )
            
    
    return lsts