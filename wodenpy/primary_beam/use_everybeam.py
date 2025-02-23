import numpy as np
from astropy.coordinates import ITRS, SkyCoord, AltAz, EarthLocation
from astropy.time import Time, TimeDelta
import astropy.units as u
import argparse
from wodenpy.use_libwoden.create_woden_struct_classes import Woden_Struct_Classes
import erfa
from typing import Union, Tuple
import concurrent.futures
from line_profiler import profile
import os
import astropy
import wodenpy
import importlib_resources
from wodenpy.use_libwoden.skymodel_structs import c_double_complex
from ctypes import c_char_p, c_int, c_double, POINTER, c_bool
import ctypes
from wodenpy.wodenpy_setup.woden_logger import simple_logger
from logging import Logger


##This call is so we can use it as a type annotation
woden_struct_classes = Woden_Struct_Classes()
Source_Catalogue = woden_struct_classes.Source_Catalogue
Woden_Settings = woden_struct_classes.Woden_Settings

def check_ms_telescope_type_matches_element_response(ms_path : str,
                                                     element_response_model : str = 'default',
                                                     logger : Logger = False) -> str:
    
    
    if not logger:
        logger = simple_logger()
    
    woden_path = importlib_resources.files(wodenpy).joinpath(f"libuse_everybeam.so")
    woden_lib = ctypes.cdll.LoadLibrary(woden_path)
    
    check_ms_telescope_type = woden_lib.check_ms_telescope_type
    check_ms_telescope_type.argtypes = [c_char_p]
    check_ms_telescope_type.restype = c_char_p
    
    ms_path_ctypes = ctypes.c_char_p(ms_path.encode('utf-8'))
    
    telescope_type = check_ms_telescope_type(ms_path_ctypes).decode('utf-8')
    
    use_element_response_model = False
    
    if telescope_type == 'MWA':
        if element_response_model == 'default':
            use_element_response_model = "MWA"
        else:
            if element_response_model != 'MWA':
                logger.warning(f"Measurement set telescope type is MWA, but element_response_model was set to {element_response_model}. Changing to 'MWA'")
                use_element_response_model = "MWA"
            else:
                use_element_response_model = element_response_model
                
    elif telescope_type == 'LOFAR':
        if element_response_model == 'default':
            use_element_response_model = "hamaker"
        else:
            if element_response_model not in ['hamaker', 'hamakerlba', 'lobes']:
                logger.warning(f"Measurement set telescope type is LOFAR, but "
                               f"element_response_model was set to {element_response_model}, "
                               "which is not one of ['hamaker', 'hamakerlba','lobes']. "
                               "Defaulting to 'hamaker'")
                use_element_response_model = "hamaker"
            else:
                use_element_response_model = element_response_model
                
    elif telescope_type == 'OSKAR':
        if element_response_model == 'default':
            use_element_response_model = "skala40_wave"
        else:
            if element_response_model != 'skala40_wave':
                logger.warning(f"Measurement set telescope type is OSKAR, but element_response_model "
                               f"was set to {element_response_model}. Changing to 'skala40_wave'")
                use_element_response_model = "skala40_wave"
            else:
                use_element_response_model = element_response_model
    else:
        use_element_response_model = element_response_model
        exit_message = f"Measurement set telescope type is {telescope_type}. "
        exit_message += "WODEN currently supports LOFAR, MWA, OSKAR EveryBeam."
        exit_message += "Cannot proceed as unknown behaviour will happen in C++ code. "
        exit_message += "Exiting now."
        logger.error(exit_message)
        exit(exit_message)
        
    
    return telescope_type, use_element_response_model


def convert_common_args_to_everybeam_args(ms_path : str, coeff_path : str,
                                          element_response_model : str,
                                          station_idxs : np.ndarray,
                                          freqs : np.ndarray,
                                          mjd_sec_times : np.ndarray):
    
    num_stations = len(station_idxs)
    num_times = len(mjd_sec_times)
    num_freqs = len(freqs)
                                          
    mjd_sec_times_ctypes = (ctypes.c_double * num_times)()
    for i in range(num_times):
        mjd_sec_times_ctypes[i] = mjd_sec_times[i]
        
    freqs_ctypes = (ctypes.c_double * num_freqs)()
    for i in range(num_freqs):
        freqs_ctypes[i] = freqs[i]
    
    station_idxs_ctypes = (ctypes.c_int * num_stations)()
    for i in range(num_stations):
        station_idxs_ctypes[i] = station_idxs[i]
    
    
    ms_path_ctypes = ctypes.c_char_p(ms_path.encode('utf-8'))
    element_response_model_ctypes = ctypes.c_char_p(element_response_model.encode('utf-8'))
    coeff_path_ctypes = ctypes.c_char_p(coeff_path.encode('utf-8'))
    
    return ms_path_ctypes, coeff_path_ctypes, element_response_model_ctypes, station_idxs_ctypes, freqs_ctypes, mjd_sec_times_ctypes

def run_everybeam(ms_path : str, coeff_path : str,
                  ras: np.ndarray, decs: np.ndarray,
                  beam_ra0: float, beam_dec0: float,
                  j2000_latitudes: np.ndarray, j2000_lsts: np.ndarray,
                  times: np.ndarray, freqs: np.ndarray,
                  station_ids: np.ndarray,
                  element_response_model='default',
                  apply_beam_norms: bool = True,
                  iau_order: bool = False,
                  element_only: bool = False,
                  parallactic_rotate: bool = False,
                  logger : Logger = False) -> np.ndarray:
    """
    Calculate the Jones matrices for a given set of coordinates, times,
    frequencies, and station ids using the EveryBeam library.
    `j2000_latitudes` should be the array latitude as precessed back to J2000,
    with `j2000_lsts` being the matching LST in J2000. 
    
    
    Parameters
    ------------
    ras : np.ndarray
        Right ascensions of the coordinates in radians.
    decs : np.ndarray
         Declinations of the coordinates in radians.
    beam_ra0 : float
        Right ascension of the beam center in radians.
    beam_dec0 : float
        Declination of the beam center in radians.
    j2000_latitudes : np.ndarray
        Latitudes in J2000 coordinates.
    j2000_lsts : np.ndarray
        Local sidereal times in J2000 coordinates.
    times : np.ndarray
        Array of observation times.
    freqs : np.ndarray
        Array of frequencies.
    telescope : (eb.Telescope):
        Telescope object from the EveryBeam library.
    station_ids : np.ndarray
        Array of station IDs.
    apply_beam_norms : bool, optional
        Whether to apply beam normalisation. Defaults to True. Achieved by
        calculating the beam response at beam centre, and multiplying all
        Jones by the inverse of this central beam response.
    iau_order : bool, optional
        If True, use IAU polarisation ordering, so set jones[0,0] to the NS dipole and jones[1,1] to EW. If False, jones[0,0] is EW.
    element_only : bool, optional
        Whether to use only the element response. Defaults to False. Use this to
        look at the dipole response only, not the beam formed response.
    parallactic_rotate : bool, optional
        Whether to apply parallactic angle rotation. Defaults to False.
        
    Returns
    --------
    np.ndarray
        The calculated Jones matrices with shape (num_stations, num_times, num_freqs, num_coords, 2, 2).
    """
    
    if not logger:
        logger = simple_logger()
        
    telescope_type, checked_element_response_model = check_ms_telescope_type_matches_element_response(ms_path,
                                                                           element_response_model,
                                                                           logger)
    
    num_stations = len(station_ids)
    num_times = len(times)
    num_freqs = len(freqs)
    num_coords = len(ras)
    
    mjd_sec_times = np.array([time.mjd * 86400.0 for time in times])
    
    if telescope_type == 'MWA':
        jones = run_mwa_beam(ms_path, checked_element_response_model,
                             coeff_path, station_ids,
                             ras, decs, mjd_sec_times,
                             j2000_lsts, j2000_latitudes, freqs, 
                             apply_beam_norms=apply_beam_norms,
                             parallactic_rotate=parallactic_rotate,
                             iau_order=iau_order,
                             element_only=element_only)
        
    elif telescope_type == 'LOFAR':
        jones = run_lofar_beam(ms_path, checked_element_response_model,
                               coeff_path, station_ids,
                               beam_ra0, beam_dec0,
                               ras, decs, mjd_sec_times, freqs,
                               apply_beam_norms=apply_beam_norms,
                               parallactic_rotate=parallactic_rotate,
                               iau_order=iau_order,
                               element_only=element_only)
        
    else:
        jones = False
    
    return jones


def run_everybeam_thread(num_threads : int, thread_id : int,
                         ms_path : str, coeff_path : str,
                         ras : np.ndarray, decs : np.ndarray,
                         ra0 : float, dec0 : float,
                         j2000_latitudes : np.ndarray, j2000_lsts : np.ndarray,
                         times : np.ndarray, freqs : np.ndarray,
                         station_ids : np.ndarray,
                         full_accuracy : bool = True,
                         apply_beam_norms : bool = True,
                         iau_order : bool = False,
                         element_only : bool = False,
                         parallactic_rotate : bool = True,
                         element_response_model='default',
                         use_local_mwa : bool = True) -> Tuple[np.ndarray, int]:
    """
    Thread function called by `run_everybeam_over_threads` to calculate the
    EveryBeam response in parrallel. Calls `run_everybeam` with a subset of
    the coordinates; see `run_everybeam` for more details of the parameters.
    
    Creates a new EveryBeam telescope object from `ms_path` for each thread.
    This has to be done because `concurrent.futures.ProcessPoolExecutor` has
    to pickle the function and all it's arguments, and EveryBeam objects can't
    be pickled. This is somewhat wasteful but I can't work out a better way
    to make things parallel.
    
    Parameters
    ------------
    num_threads : int
        Number of threads being in call by `run_everybeam_over_threads`.
    thread_id : int
        ID of the current thread. Useds to work out what chunk of `ras` and `decs`
        to process.
    ms_path : str
        Path to the measurement set to load the EveryBeam telescope from.
    ras : np.ndarray
        Right ascensions of the coordinates in radians.
    decs : np.ndarray
         Declinations of the coordinates in radians.
    beam_ra0 : float
        Right ascension of the beam center in radians.
    beam_dec0 : float
        Declination of the beam center in radians.
    j2000_latitudes : np.ndarray
        Latitudes in J2000 coordinates.
    j2000_lsts : np.ndarray
        Local sidereal times in J2000 coordinates.
    times : np.ndarray
        Array of observation times.
    freqs : np.ndarray
        Array of frequencies.
    station_ids : np.ndarray
        Array of station IDs.
    full_accuracy : bool, optional
        Whether to use the full accuracy of the EveryBeam library. Defaults to True.
        If False, the array factor and element response are calculated separately
        (the two values that multiply to give the full response). The array factor
        is only calculated at the middle time and frequency, under the assumption
        that the range of frequencies and times is small enough that the array factor
        will not change significantly. Magnitude of the differences vary by
        frequency and direction so use with caution.
    apply_beam_norms : bool, optional
        Whether to apply beam normalisation. Defaults to True. Achieved by
        calculating the beam response at beam centre, and multiplying all
        Jones by the inverse of this central beam response.
    iau_order : bool, optional
        If True, use IAU polarisation ordering, so set jones[0,0] to the NS dipole and jones[1,1] to EW. If False, jones[0,0] is EW.
    element_only : bool, optional
        Whether to use only the element response. Defaults to False. Use this to
        look at the dipole response only, not the beam formed response.
    parallactic_rotate : bool, optional
        Whether to apply parallactic angle rotation. Defaults to False.
    element_response_model : str, optional
        The Everybeam element response model to use. Defaults to 'hamaker'.
        Avaible options are 'hamaker' (LOFAR), 'skala40_wave' (OSKAR), and 'MWA' (MWA).
    use_local_mwa : bool, optional
        Whether to use the local MWA model. Defaults to True. The local MWA model
        takes za/az instead of RA/Dec.
        
    Returns
    --------
    Tuple[np.ndarray, int]
        The calculated Jones matrices with shape
        (num_stations, num_times, num_freqs, num_coords_in_thread, 2, 2), as
        well as the thread ID. Use the thread ID to insert this thread output
        into the correct place in the final Jones matrix.
        
    """
    
    num_coords = len(ras)
    coords_per_thread = int(np.ceil(num_coords / num_threads))
    
    low_coord = thread_id * coords_per_thread
    high_coord = (thread_id + 1) * coords_per_thread
    
    print(f"Thread {thread_id} processing coords {low_coord} to {high_coord}")
    
    jones = run_everybeam(ms_path, coeff_path,
                  ras[low_coord:high_coord],
                  decs[low_coord:high_coord],
                  ra0, dec0, j2000_latitudes, j2000_lsts,
                  times, freqs,
                  station_ids,
                  element_response_model=element_response_model,
                  apply_beam_norms=apply_beam_norms,
                  iau_order=iau_order,
                  element_only=element_only,
                  parallactic_rotate=parallactic_rotate)
    
    print(f"Thread {thread_id} finished")
    
    return jones, thread_id

def run_everybeam_over_threads(num_threads : int,
                               ms_path : str,
                               coeff_path : str,
                               ras : np.ndarray, decs : np.ndarray,
                               ra0 : float, dec0 : float,
                               j2000_latitudes : np.ndarray, j2000_lsts : np.ndarray,
                               times : np.ndarray, freqs : np.ndarray,
                               station_ids : np.ndarray,
                               full_accuracy : bool = True,
                               apply_beam_norms : bool = True,
                               iau_order : bool = False,
                               element_only : bool = False,
                               parallactic_rotate : bool = True,
                               use_local_mwa : bool = True,
                               element_response_model='default'):
    """
    Runs `run_everybeam` in parallel over `num_threads` threads, using
    `concurrent.futures.ProcessPoolExecutor`. See `run_everybeam` for more
    details of what each parameter does.
    
    Creates a new EveryBeam telescope object from `ms_path` for each thread.
    This has to be done because `concurrent.futures.ProcessPoolExecutor` has
    to pickle the function and all it's arguments, and EveryBeam objects can't
    be pickled. This is somewhat wasteful but I can't work out a better way
    to make things parallel.
    
    Parameters
    ------------
    num_threads : int
        Number of threads being in call by `run_everybeam_over_threads`.
    ms_path : str
        Path to the measurement set to load the EveryBeam telescope from.
    ras : np.ndarray
        Right ascensions of the coordinates in radians.
    decs : np.ndarray
         Declinations of the coordinates in radians.
    beam_ra0 : float
        Right ascension of the beam center in radians.
    beam_dec0 : float
        Declination of the beam center in radians.
    j2000_latitudes : np.ndarray
        Latitudes in J2000 coordinates.
    j2000_lsts : np.ndarray
        Local sidereal times in J2000 coordinates.
    times : np.ndarray
        Array of observation times.
    freqs : np.ndarray
        Array of frequencies.
    station_ids : np.ndarray
        Array of station IDs.
    full_accuracy : bool, optional
        Whether to use the full accuracy of the EveryBeam library. Defaults to True.
        If False, the array factor and element response are calculated separately
        (the two values that multiply to give the full response). The array factor
        is only calculated at the middle time and frequency, under the assumption
        that the range of frequencies and times is small enough that the array factor
        will not change significantly. Magnitude of the differences vary by
        frequency and direction so use with caution.
    apply_beam_norms : bool, optional
        Whether to apply beam normalisation. Defaults to True. Achieved by
        calculating the beam response at beam centre, and multiplying all
        Jones by the inverse of this central beam response.
    iau_order : bool, optional
        If True, use IAU polarisation ordering, so set jones[0,0] to the NS dipole and jones[1,1] to EW. If False, jones[0,0] is EW.
    element_only : bool, optional
        Whether to use only the element response. Defaults to False. Use this to
        look at the dipole response only, not the beam formed response.
    parallactic_rotate : bool, optional
        Whether to apply parallactic angle. Defaults to False.
    para_angle_offset : float, optional
        Offset to add to the parallactic angle. Defaults to 0.
    element_response_model : str, optional
        The Everybeam element response model to use. Defaults to 'hamaker'.
        Avaible options are 'hamaker' (LOFAR), 'skala40_wave' (OSKAR), and 'MWA' (MWA).
    use_local_mwa : bool, optional
        Whether to use the local MWA model. Defaults to True. The local MWA model
        takes za/az instead of RA/Dec.
        
    Returns
    --------
    np.ndarray
        The calculated Jones matrices with shape
        (num_stations, num_times, num_freqs, num_coord, 2, 2)
    """
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
        
            future_data = [executor.submit(run_everybeam_thread,
                                           num_threads, thread_id,
                                           ms_path, 
                                           coeff_path,
                                           ras, decs,
                                           ra0, dec0,
                                           j2000_latitudes, j2000_lsts,
                                           times, freqs,
                                           station_ids,
                                           full_accuracy=full_accuracy,
                                           apply_beam_norms=apply_beam_norms,
                                           iau_order=iau_order,
                                           element_only=element_only,
                                           parallactic_rotate=parallactic_rotate,
                                           element_response_model=element_response_model,
                                           use_local_mwa=use_local_mwa)
                                    for thread_id in range(num_threads)]
            
            all_jones_chunks = []
            all_thread_ids = []
            for future in concurrent.futures.as_completed(future_data):
                jones_chunk, thread_id = future.result()
                all_jones_chunks.append(jones_chunk)
                all_thread_ids.append(thread_id)
                
    num_stations = len(station_ids)
    num_times = len(times)
    num_freqs = len(freqs)
    num_coords = len(ras)
    
    all_jones = np.zeros((num_stations, num_times, num_freqs, num_coords, 2, 2), dtype=np.complex128)*np.nan
    
    coords_per_thread = int(np.ceil(num_coords / num_threads))
    
    for jones_chunk, thread_id in zip(all_jones_chunks, all_thread_ids):
        
        low_coord = thread_id * coords_per_thread
        high_coord = (thread_id + 1) * coords_per_thread
        
        all_jones[:, :, :, low_coord:high_coord, :, :] = jones_chunk
    
    return all_jones




def run_lofar_beam(ms_path : str, element_response_model : bool,
                   coeff_path : str,
                   station_idxs : np.ndarray,
                   beam_ra0 : float, beam_dec0 : float,
                   ras : np.ndarray, decs : np.ndarray,
                   mjd_sec_times : np.ndarray,
                   freqs : np.ndarray,
                   apply_beam_norms : bool,
                   parallactic_rotate : bool,
                   iau_order : bool = True,
                   element_only : bool = False):
    
    woden_path = importlib_resources.files(wodenpy).joinpath(f"libuse_everybeam.so")
    woden_lib = ctypes.cdll.LoadLibrary(woden_path)
    
    load_and_run_lofar_beam = woden_lib.load_and_run_lofar_beam
    
    num_stations = len(station_idxs)
    num_dirs = len(ras)
    num_freqs = len(freqs)
    num_times = len(mjd_sec_times)
    
    ras_ctypes = (ctypes.c_double * num_dirs)()
    decs_ctypes = (ctypes.c_double * num_dirs)()
    for i in range(num_dirs):
        ras_ctypes[i] = ras[i]
        decs_ctypes[i] = decs[i]
        
    mjd_sec_times_ctypes = (ctypes.c_double * num_times)()
    for i in range(num_times):
        mjd_sec_times_ctypes[i] = mjd_sec_times[i]
        
    freqs_ctypes = (ctypes.c_double * num_freqs)()
    for i in range(num_freqs):
        freqs_ctypes[i] = freqs[i]
    
    station_idxs_ctypes = (ctypes.c_int * num_stations)()
    for i in range(num_stations):
        station_idxs_ctypes[i] = station_idxs[i]
    
    
    ms_path_ctypes = ctypes.c_char_p(ms_path.encode('utf-8'))
    element_response_model_ctypes = ctypes.c_char_p(element_response_model.encode('utf-8'))
    coeff_path_ctypes = ctypes.c_char_p(coeff_path.encode('utf-8'))
    
    jones = ((num_stations*num_times*num_freqs*num_dirs*4)*c_double_complex)()
    
    load_and_run_lofar_beam.argtypes = [c_char_p, c_char_p, c_char_p,
                                        c_int, POINTER(c_int),
                                        c_int, c_double, c_double,
                                        POINTER(c_double), POINTER(c_double),
                                        c_int, POINTER(c_double),
                                        c_int, POINTER(c_double),
                                        c_bool, c_bool, c_bool, c_bool,
                                        POINTER(c_double_complex)]
    
    load_and_run_lofar_beam(ms_path_ctypes,
                            element_response_model_ctypes,
                            coeff_path_ctypes,
                            num_stations, station_idxs_ctypes,
                            num_dirs,
                            beam_ra0, beam_dec0,
                            ras_ctypes, decs_ctypes,
                            num_times, mjd_sec_times_ctypes,
                            num_freqs, freqs_ctypes,
                            apply_beam_norms, parallactic_rotate,
                            element_only, iau_order,
                            jones)
    
    # print(jones)
    
    jones_py = np.ctypeslib.as_array(jones, shape=(num_stations*num_times*num_freqs*num_dirs*4))
    jones_py = jones_py['real'] + 1j*jones_py['imag']
    
    jones_py = jones_py.reshape(num_stations, num_times, num_freqs, num_dirs, 2, 2)
    
    
    return jones_py



def run_mwa_beam(ms_path : str, element_response_model : bool,
                   coeff_path : str,
                   station_idxs : np.ndarray,
                   ras : np.ndarray, decs : np.ndarray,
                   mjd_sec_times : np.ndarray,
                   j2000_lsts : np.ndarray, j2000_latitudes : np.ndarray,
                   freqs : np.ndarray,
                   apply_beam_norms : bool, parallactic_rotate : bool,
                   iau_order : bool = True,
                   element_only : bool = False):
    
    num_stations = len(station_idxs)
    num_dirs = len(ras)
    num_freqs = len(freqs)
    num_times = len(mjd_sec_times)
    
    azs = np.empty(num_dirs*num_times)
    zas = np.empty(num_dirs*num_times)
    para_angles = np.empty(num_dirs*num_times)
    
    for comp_ind in range(num_dirs):
        comp_has = j2000_lsts - ras[comp_ind]
        these_azs, these_els = erfa.hd2ae(comp_has, decs[comp_ind], j2000_latitudes)
        these_zas = np.pi/2 - these_els
        
        these_para_angles = erfa.hd2pa(comp_has, decs[comp_ind], j2000_latitudes)
    
        azs[comp_ind*num_times:(comp_ind+1)*num_times] = these_azs
        zas[comp_ind*num_times:(comp_ind+1)*num_times] = these_zas
        para_angles[comp_ind*num_times:(comp_ind+1)*num_times] = these_para_angles
        
    woden_path = importlib_resources.files(wodenpy).joinpath(f"libuse_everybeam.so")
    
    woden_path = "/home/jack-line/software/WODEN_dev/build/cmake_testing/GPU_or_C_code/libuse_everybeam.so"
    
    woden_lib = ctypes.cdll.LoadLibrary(woden_path)
    
    load_and_run_mwa_beam = woden_lib.load_and_run_mwa_beam
    
    zas_ctypes = (ctypes.c_double*(num_dirs*num_times))()
    azs_ctypes = (ctypes.c_double*(num_dirs*num_times))()
    para_angles_ctypes = (ctypes.c_double*(num_dirs*num_times))()
    for i in range(num_dirs*num_times):
        zas_ctypes[i] = zas[i]
        azs_ctypes[i] = azs[i]
        para_angles_ctypes[i] = para_angles[i]
    
    ms_path_ctypes, coeff_path_ctypes, \
    element_response_model_ctypes, \
    station_idxs_ctypes, freqs_ctypes, \
    mjd_sec_times_ctypes = convert_common_args_to_everybeam_args(ms_path, 
                                    coeff_path, element_response_model,
                                    station_idxs, freqs, mjd_sec_times)
    
    jones = ((num_stations*num_times*num_freqs*num_dirs*4)*c_double_complex)()
    
    load_and_run_mwa_beam.argtypes = [c_char_p, c_char_p, c_char_p,
                                        c_int, POINTER(c_int),
                                        c_int, 
                                        POINTER(c_double), POINTER(c_double),
                                        POINTER(c_double),
                                        c_int, POINTER(c_double),
                                        c_int, POINTER(c_double),
                                        c_bool, c_bool, c_bool, c_bool,
                                        POINTER(c_double_complex)]
    
    load_and_run_mwa_beam(ms_path_ctypes,
                            element_response_model_ctypes,
                            coeff_path_ctypes,
                            num_stations, station_idxs_ctypes,
                            num_dirs,
                            azs_ctypes, zas_ctypes,
                            para_angles_ctypes,
                            num_times, mjd_sec_times_ctypes,
                            num_freqs, freqs_ctypes,
                            apply_beam_norms, parallactic_rotate, element_only,
                            iau_order,
                            jones)
    
    jones_py = np.ctypeslib.as_array(jones, shape=(num_stations*num_times*num_freqs*num_dirs*4))
    jones_py = jones_py['real'] + 1j*jones_py['imag']
    
    jones_py = jones_py.reshape(num_stations, num_times, num_freqs, num_dirs, 2, 2)
    
    
    return jones_py