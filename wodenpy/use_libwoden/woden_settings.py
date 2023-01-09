import ctypes 
# import importlib_resources
from wodenpy.woden_lib import *
import numpy as np
import argparse

from ctypes import c_double, c_float, c_int, POINTER, c_ulong, c_char, create_string_buffer, c_uint, c_long, c_uint64

D2R = np.pi/180.0

class Woden_Settings_Double(ctypes.Structure):
    """A class structured equivalently to a `visi_set` struct, used by 
    the C and CUDA code in libwoden_double.so
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
                ("gauss_beam_FWHM", c_double),
                ("gauss_beam_ref_freq", c_double),
                ("chunking_size", c_uint64),
                ("hdf5_beam_path", POINTER(c_char)),
                ("jd_date", c_double),
                ("array_layout_file", c_int),
                ("array_layout_file_path", POINTER(c_char)),
                ("latitude", c_double),
                ("latitude_obs_epoch_base", c_double),
                ("longitude", c_double),
                ("FEE_ideal_delays", (c_double*16)),
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
                ("do_autos", c_int)]

##TODO gotta be a way to set the float or double fields via some kind of
##variable instead of making two different classes
class Woden_Settings_Float(ctypes.Structure):
    """A class structured equivalently to a `visi_set` struct, used by 
    the C and CUDA code in libwoden_float.so
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
                ("gauss_beam_FWHM", c_float),
                ("gauss_beam_ref_freq", c_double),
                ("chunking_size", c_ulong),
                ("hdf5_beam_path", POINTER(c_char)),
                ("jd_date", c_double),
                ("array_layout_file", c_int),
                ("array_layout_file_path", POINTER(c_char)),
                ("latitude", c_double),
                ("latitude_obs_epoch_base", c_double),
                ("longitude", c_float),
                ("FEE_ideal_delays", (c_float*16)),
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
                ("do_autos", c_int)]
    
def create_woden_settings(args : argparse.Namespace,
                          jd_date : float, lst : float) -> ctypes.Structure:
    """Given the parsed and checked arguments in `args`, populate a
    `woden_settings` ctypes.Structure that can be passed into
    libwoden_float.so or libwoden_double.so, depending on the desired
    `precision`.

    Parameters
    ----------
    args : argparse.Namespace
        The populated arguments `args = parser.parse_args()`` as returned from
        the parser given by :func:`~wodenpy.run_setup.check_args`
    jd_date : float
        The julian date of the first time step of the observation
    lst : float
        The LST of the first time step of the observation (degrees)

    Returns
    -------
    woden_settings : ctypes.Structure
       Populated ctype struct that can be passed into the C/CUDA code
    """
    
    if args.precision == 'float':
        woden_settings = Woden_Settings_Float()
    else:
        woden_settings = Woden_Settings_Double()
        
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
        delays = np.array(args.MWA_FEE_delays.strip('[]').split(','))
                          
        for ind, delay in enumerate(delays):
            woden_settings.FEE_ideal_delays[ind] = float(delay)
        
    if args.primary_beam == 'none':
        woden_settings.beamtype = 0
    
    elif args.primary_beam == 'Gaussian':
        woden_settings.beamtype = 1
        woden_settings.gauss_beam_FWHM = float(args.gauss_beam_FWHM)
        woden_settings.gauss_beam_ref_freq = float(args.gauss_beam_ref_freq)
        woden_settings.gauss_ra_point = float(args.gauss_ra_point)
        woden_settings.gauss_dec_point = float(args.gauss_dec_point)

    elif args.primary_beam == 'MWA_FEE':
        woden_settings.beamtype = 2
        woden_settings.hdf5_beam_path = create_string_buffer(args.hdf5_beam_path.encode('utf-8'))
        
    elif args.primary_beam == 'EDA2':
        woden_settings.beamtype = 3
        
    elif args.primary_beam == 'MWA_FEE_interp':
        woden_settings.beamtype = 4
        woden_settings.hdf5_beam_path = create_string_buffer(args.hdf5_beam_path.encode('utf-8'))
        
    elif args.primary_beam == 'MWA_analy':
        woden_settings.beamtype = 5
        
    if args.no_precession:
        woden_settings.do_precession = 0
    else:
        woden_settings.do_precession = 1
        
    woden_settings.do_autos = args.do_autos
    
    woden_settings.chunking_size = int(args.chunking_size)
    
    ##TODO array layout will be used directly inside python code, no need to
    ##write text file-----------------------------------------------------------
    if args.array_layout == 'from_the_metafits':
        json_band_str =  '-'.join(map(str, args.band_nums))
        band_array_layout = f"WODEN_array_layout_band{json_band_str}.txt"
        command(f"cp WODEN_array_layout.txt {band_array_layout}")
    else:
        band_array_layout = args.array_layout_name

    woden_settings.array_layout_file = 1
    woden_settings.array_layout_file_path = create_string_buffer(band_array_layout.encode('utf-8'))
    ##--------------------------------------------------------------------------

    woden_settings.num_bands = len(args.band_nums)
    woden_settings.band_nums = (ctypes.c_int*woden_settings.num_bands)()
    
    for ind, band in enumerate(args.band_nums):
            woden_settings.band_nums[ind] = int(band)
    
    # woden_settings.sdec0
    # woden_settings.cdec0
    # woden_settings.num_baselines
    # woden_settings.num_ants
    
    # woden_settings.longitude
    # 
    # woden_settings.num_cross
    # woden_settings.num_autos
    # woden_settings.num_visis
    # woden_settings.base_band_freq
    # woden_settings.lsts
    # woden_settings.latitudes
    # woden_settings.mjds
        
    return woden_settings
    
