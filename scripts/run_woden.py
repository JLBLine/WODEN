#!/usr/bin/env python3
"""Wrapper script to run the WODEN simulator. Author: J.L.B. Line
"""
import numpy as np
import sys
import os
import importlib_resources
import wodenpy
from time import time
import argparse
from threading import Thread
from ctypes import POINTER, c_double, c_float, c_int, create_string_buffer, string_at
import ctypes
from typing import Union, Tuple, List, Callable
from astropy.table import Table
from astropy.io import fits
from wodenpy.use_libwoden.beam_settings import BeamTypes, BeamGroups
from wodenpy.use_libwoden.woden_settings import fill_woden_settings_python, setup_lsts_and_phase_centre, Woden_Settings_Python, convert_woden_settings_to_ctypes
from wodenpy.use_libwoden.visibility_set import setup_visi_set_array, load_visibility_set, Visi_Set_Python
from wodenpy.wodenpy_setup.run_setup import get_parser, check_args, get_code_version
from wodenpy.use_libwoden.use_libwoden import load_in_run_woden
from wodenpy.observational.calc_obs import get_uvfits_date_and_position_constants, calc_jdcal
from wodenpy.skymodel.woden_skymodel import crop_below_horizon
from wodenpy.skymodel.read_skymodel import get_skymodel_tables, read_radec_count_components, create_source_catalogue_from_python_sources
from wodenpy.skymodel.chunk_sky_model import create_skymodel_chunk_map, reshape_chunked_skymodel_map_sets, Skymodel_Chunk_Map
from wodenpy.array_layout.create_array_layout import calc_XYZ_diffs, enh2xyz, Array_Layout_Python, convert_array_layout_to_ctypes, setup_array_layout_python
from wodenpy.use_libwoden.array_layout_struct import Array_Layout_Ctypes
from wodenpy.uvfits.wodenpy_uvfits import make_antenna_table, make_baseline_date_arrays, create_uvfits
from wodenpy.phase_rotate.remove_phase_track import remove_phase_tracking
from wodenpy.use_libwoden.shapelets import create_sbf
from wodenpy.use_libwoden.create_woden_struct_classes import Woden_Struct_Classes
from wodenpy.skymodel.read_fits_skymodel import read_fits_skymodel_chunks
from wodenpy.use_libwoden.skymodel_structs import setup_source_catalogue, Source_Python
import concurrent.futures
from logging import Logger
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed, Future
from line_profiler import profile, LineProfiler
from wodenpy.primary_beam.use_everybeam import run_everybeam
from wodenpy.wodenpy_setup.woden_logger import get_log_callback, set_logger_header, simple_logger, log_chosen_beamtype, summarise_input_args, set_woden_logger
from copy import deepcopy
from multiprocessing import Process, Queue
from datetime import timedelta
import traceback
import shutil
from wodenpy.primary_beam.use_uvbeam import setup_MWA_uvbeams, setup_HERA_uvbeams_from_CST, setup_HERA_uvbeams_from_single_file

##Constants
R2D = 180.0 / np.pi
D2R = np.pi / 180.0
MWA_LAT = -26.703319405555554
MWA_LONG = 116.67081523611111
MWA_HEIGHT = 377.827
VELC = 299792458.0
SOLAR2SIDEREAL = 1.00274

sbf_N = 101
sbf_L = 10001
sbf_c = 5000
sbf_dx = 0.01

##This call is so we can use it as a type annotation
woden_struct_classes = Woden_Struct_Classes()
Woden_Settings = woden_struct_classes.Woden_Settings
Visi_Set = woden_struct_classes.Visi_Set
Source_Catalogue = woden_struct_classes.Source_Catalogue

def woden_worker(thread_ind : int, 
                all_loaded_python_sources : List[List[Source_Python]],
                all_loaded_sources_orders : List[int], round_num : int,
                woden_settings_python : Woden_Settings_Python,
                array_layout_python : Array_Layout_Python,
                visi_sets_python : List[Visi_Set_Python],
                beamtype : int, logging_level : int = logging.DEBUG,
                precision : str = "double", 
                logger : Logger = False,
                profile = False) -> Tuple[List[Visi_Set_Python], int, int]:
    """
    This worker function runs WODEN CPU/GPU code, to process chunked source catalogues
    from a queue until the queue is empty.
    
    Parameters
    ----------
    thread_ind : int
        The index of the thread this worker is running in. Just used for logging here.
    all_loaded_python_sources : List[List[Source_Python]]
        A list of lists of `Source_Python` to be processed, as ouput by `read_skymodel_worker`,
        where each list of `Source_Python` matches a chunk_map in `chunked_skymodel_map_set`.
        Each entry is a list itself as the Shapelet chunking can result in multiple sources
        per thread per round, as the chunking happens by sky direction as well as number of
        shapelet coefficients.
    all_loaded_sources_orders : List[int]
        The order of the sources in `all_loaded_python_sources` as matched to `chunked_skymodel_map_set`;
        orders also output by `read_skymodel_worker` (different threads finish in different times so
        the order of the sources in `all_loaded_python_sources` may not match the order of the chunk_maps)
    round_num : int
        The round number of the processing
    woden_settings_python : Woden_Settings_Python
        The WODEN settings to be used.
    array_layout_python : Array_Layout_Python
        The array layout to be used.
    visi_sets_python : List[Visi_Set_Python]
        The existing visibility sets information. The results of this worker
        will be added to this and returned. Should have length of `num_bands`,
        where `num_bands` is the number of frequency bands.
    beamtype : int
        The value of the type of beam being used, e.g. BeamTypes.FEE_BEAM.value
    logging_level : int
        The logging level to use for logging. Default is logging.DEBUG.
    precision : str
        The precision to use for the WODEN processing. Default is "double".
        Should be either "float" or "double".
    logger : Logger
        The logger to use for logging. If not provided, a default logger
        will be created.
    profile : bool
        Whether to profile the function. This is a VERY experimental
        way of doing things, and only includes some specific functions. Use
        with caution. Default is False.
    """
    
    try:
    
        if woden_settings_python.do_gpu == 1:
            device = 'GPU'
        else:
            device = 'CPU'
        
        woden_struct_classes = Woden_Struct_Classes(precision)
        
        ##Depending on what precision was selected by the user, load in the
        # ##C/CUDA library and return the `run_woden` function
        # run_woden = load_in_woden_library(woden_struct_classes)
        
        woden_lib_path = importlib_resources.files(wodenpy).joinpath(f"libwoden_{woden_struct_classes.precision}.so")
        woden_lib = ctypes.cdll.LoadLibrary(woden_lib_path)
        
        ##This lets the WODEN C/C++/GPU code write to the logger
        c_log_callback = get_log_callback(logger, logging_level)
        woden_lib.set_log_callback(c_log_callback)
        
        run_woden = load_in_run_woden(woden_lib, woden_struct_classes)
        
        python_sources = []
        
        ##grab all the python sources, and reorder/flatten them to match the order of
        ##chunked_skymodel_maps. Not actually sure we need to do this, but if we
        ##change this order it'll mean we have to change the order of a bunch of
        ##tests. I don't think it costs much, and saves a bunch of developing time
        ordering = np.argsort(all_loaded_sources_orders)
        for order in ordering:
            sources = all_loaded_python_sources[order]
            for source in sources:
                python_sources.append(source)
                
        if len(python_sources) == 0:
            logger.warning(f"Visibility processing thread {thread_ind} has no work to do")
            # print(f"Visibility processing thread {thread_ind} has no work to do")
            return visi_sets_python, thread_ind, round_num
        
        if woden_settings_python.beamtype in BeamGroups.python_calc_beams:
            start = time()
            logger.info(f"Coverting sky model and beam values into `ctypes` for sky set {round_num} chunk {thread_ind}")
        ##Create a ctypes Source_Catalogue from the python sources to feed the GPU
        source_catalogue = create_source_catalogue_from_python_sources(python_sources,
                                                                    woden_struct_classes,
                                                                    beamtype, precision)
        if woden_settings_python.beamtype in BeamGroups.python_calc_beams:
            end = time()
            logger.info(f"`ctype` catalogue for round {round_num} chunk {thread_ind} created in {end-start:.1f} seconds")
            # print(f"Created source catalogue for round {round_num} in {end-start:.1f} seconds")
        
        woden_settings = woden_struct_classes.Woden_Settings()
        woden_settings = convert_woden_settings_to_ctypes(woden_settings_python, woden_settings)
        
        visibility_set = setup_visi_set_array(woden_struct_classes.Visi_Set,
                                                woden_settings.num_bands,
                                                woden_settings.num_visis,
                                                precision=precision)
        
        array_layout = Array_Layout_Ctypes()
        array_layout = convert_array_layout_to_ctypes(array_layout_python, array_layout)
        
        sbf = create_sbf(precision=precision)
        
        ##We send all sky models chunks to the GPU at once, and then process them
        if device == 'GPU':
            logger.info(f"Sending Sky set {round_num} to the {device}")
        ##We send each chunk to a separate CPU thread, so report the thread ind
        else:
            logger.info(f"Sending Sky set {round_num} chunk {thread_ind} to the {device}")
        # print(f"Sending Sky set {round_num} thread {thread_ind} to CPU")
        start = time()
                
        run_woden(woden_settings, visibility_set, source_catalogue, array_layout,
                    sbf)
        end = time()
        
        
        if device == 'GPU':
            logger.info(f"Sky set {round_num} has returned from the {device} after {end-start:.1f} seconds")
        else:
            logger.info(f"Sky set {round_num} chunk {thread_ind} has returned from the {device} after {end-start:.1f} seconds")
            
        for band_ind in range(woden_settings_python.num_bands):
            visi_sets_python[band_ind].us_metres = deepcopy(np.ctypeslib.as_array(visibility_set[band_ind].us_metres, shape=(woden_settings_python.num_visis,)))
            visi_sets_python[band_ind].vs_metres = deepcopy(np.ctypeslib.as_array(visibility_set[band_ind].vs_metres, shape=(woden_settings_python.num_visis,)))
            visi_sets_python[band_ind].ws_metres = deepcopy(np.ctypeslib.as_array(visibility_set[band_ind].ws_metres, shape=(woden_settings_python.num_visis,)))
            
            visi_sets_python[band_ind].sum_visi_XX_real += np.ctypeslib.as_array(visibility_set[band_ind].sum_visi_XX_real, shape=(woden_settings_python.num_visis,))
            visi_sets_python[band_ind].sum_visi_XX_imag += np.ctypeslib.as_array(visibility_set[band_ind].sum_visi_XX_imag, shape=(woden_settings_python.num_visis,))
            visi_sets_python[band_ind].sum_visi_XY_real += np.ctypeslib.as_array(visibility_set[band_ind].sum_visi_XY_real, shape=(woden_settings_python.num_visis,))
            visi_sets_python[band_ind].sum_visi_XY_imag += np.ctypeslib.as_array(visibility_set[band_ind].sum_visi_XY_imag, shape=(woden_settings_python.num_visis,))
            visi_sets_python[band_ind].sum_visi_YX_real += np.ctypeslib.as_array(visibility_set[band_ind].sum_visi_YX_real, shape=(woden_settings_python.num_visis,))
            visi_sets_python[band_ind].sum_visi_YX_imag += np.ctypeslib.as_array(visibility_set[band_ind].sum_visi_YX_imag, shape=(woden_settings_python.num_visis,))
            visi_sets_python[band_ind].sum_visi_YY_real += np.ctypeslib.as_array(visibility_set[band_ind].sum_visi_YY_real, shape=(woden_settings_python.num_visis,))
            visi_sets_python[band_ind].sum_visi_YY_imag += np.ctypeslib.as_array(visibility_set[band_ind].sum_visi_YY_imag, shape=(woden_settings_python.num_visis,))
            
        return visi_sets_python, thread_ind, round_num
    except Exception as e:
        return f"woden_worker: {e}\n{traceback.format_exc()}"

def read_skymodel_worker(thread_id : int, num_threads : int, 
                         chunked_skymodel_map_sets: List[Skymodel_Chunk_Map],
                         lsts : np.ndarray, latitudes : np.ndarray,
                         args : argparse.Namespace,
                         beamtype : int,
                         main_table : Table, shape_table : Table,
                         v_table : Table = False, q_table : Table = False,
                         u_table : Table = False, p_table : Table = False,
                         logger : Logger = False,
                         uvbeam_objs : np.ndarray = None) -> Tuple[List[Source_Python], int]:
    """
    Reads a list of `Skymodel_Chunk_Map` types in `chunked_skymodel_map_sets`
    from the given astropy tables, into a list of `Source_Python` types.
    Depending on the type of primary beam specified by `beamtype`
    
    Parameters
    ==========
    thread_id : int
        The ID of the current thread.
    num_threads : int
        The total number of threads.
    chunked_skymodel_map_sets : List[Skymodel_Chunk_Map]
        A list of chunked skymodel map sets.
    lsts : np.ndarray
        Local Sidereal Times (LSTs) for all time steps.
    latitudes : np.ndarray
        Latitudes of the array at all time steps; when doing array precession,
        these will all be slightly different.
    args : argparse.Namespace
        Command-line arguments namespace.
    beamtype : int
        The value of the type of beam being used, e.g. BeamTypes.FEE_BEAM.value
    main_table : Table
        Main table containing skymodel data.
    shape_table : Table
        Table containing Shapelet data.
    v_table : Table, optional
        Table containing V polarization data (default is False).
    q_table : Table, optional
        Table containing Q polarization data (default is False).
    u_table : Table, optional
        Table containing U polarization data (default is False).
    p_table : Table, optional
        Table containing P polarization data (default is False).
    logger : Logger, optional
        The logger to use for logging (default is False). If not provided,
        a default logger will be created.
    uvbeam_objs : np.ndarray, optional
        The UVBeam objects to use for the primary beam calculations. Only needed
        if the beamtype is in `BeamGroups.uvbeam_beams`.
        
    Returns
    =======
    tuple : Tuple[List[Source_Python], int]
        A tuple containing the list of `Source_Python`s and the thread number,
        where thead number is `thread_id % num_threads`; this can be used to
        match the order of recovered data to the order of the chunked maps.
    
    """
    if not logger:
        logger = simple_logger(args.log_level)
    
    try:
        set_ind = thread_id // num_threads
        thread_num = thread_id % num_threads
        
        start = time()
        
        chunk_maps = chunked_skymodel_map_sets[set_ind][thread_num]
        
        n_p, n_g, n_s, n_c = 0, 0, 0, 0
        for ind, chunk_map in enumerate(chunk_maps):
            n_p += chunk_map.n_points
            n_g += chunk_map.n_gauss
            n_s += chunk_map.n_shapes
            n_c += chunk_map.n_shape_coeffs
            
        msg = f"From sky set {set_ind} thread num {thread_num} reading {n_p} points, {n_g} gauss, {n_s} shape, {n_c} shape coeffs"
        # msg = f'cmon man {thread_id}'
        logger.info(msg)
                    
        if len(chunk_maps) == 0:
            logger.warning(f"Sky model reading thread {thread_id} has no work to do")
            # print(f"Sky model reading thread {thread_id} has no work to do")
            return [], thread_num
        
        if args.profile:
            profiler = LineProfiler()
            ##Add whatever functions that are called in `read_fits_skymodel_chunks`
            ##here to profile them
            profiler.add_function(read_fits_skymodel_chunks)
            profiler.add_function(run_everybeam)
            profiler.enable()
        
        python_sources = read_fits_skymodel_chunks(args, main_table, shape_table,
                                chunk_maps,
                                args.num_freq_channels, args.num_time_steps,
                                beamtype, lsts, latitudes,
                                v_table, q_table,
                                u_table, p_table,
                                args.precision, uvbeam_objs,
                                logger=logger)
        
        end = time()
        
        if beamtype in BeamGroups.python_calc_beams:
            logger.info(f"Finshed sky set {set_ind} reading thread num {thread_num} in {end-start:.1f} seconds")
        else:
            logger.debug(f"Finshed sky set {set_ind} reading thread num {thread_num} in {end-start:.1f} seconds")
        
        if args.profile:
            profiler.disable()
            profile_filename = f"wod_read_skymodel_worker_{os.getpid()}_{set_ind:04d}.lprof"
            logger.info(f"Dumping profile for set {set_ind} thread {thread_num} to {profile_filename}")
            profiler.dump_stats(profile_filename)
            
        return python_sources, thread_num
    except Exception as e:
        return f"read_skymodel_worker: {e}\n{traceback.format_exc()}"


def sum_components_in_chunked_skymodel_map_sets(chunked_skymodel_map_sets):
    """
    Sums the number of components in a list of `Skymodel_Chunk_Map` types.
    
    Parameters
    ==========
    chunked_skymodel_map_sets : List[Skymodel_Chunk_Map]
        A list of chunked skymodel map sets.
    
    Returns
    =======
    int : The total number of components in the list of `Skymodel_Chunk_Map` types.
    
    
    """
    
    n_points, n_gauss, n_shapes, n_shape_coeffs = 0, 0, 0, 0
    for chunk_map_set in chunked_skymodel_map_sets:
        for chunk_maps in chunk_map_set:
            if len(chunk_maps) > 0:
                for chunk_map in chunk_maps:
                
                    n_points += chunk_map.n_points
                    n_gauss += chunk_map.n_gauss
                    n_shapes += chunk_map.n_shapes
                    n_shape_coeffs += chunk_map.n_shape_coeffs
        
    return n_points, n_gauss, n_shapes, n_shape_coeffs

def get_future_result(future : Future, logger : Logger,
                      function : Callable):
    """
    This function gets the result of a future, and handles any exceptions
    that may occur. It logs the error and exits the program if an exception
    occurs.
    Parameters
    ----------
    future : Future
        The future to get the result from, as returned by the `ProcessPoolExecutor`.
    logger : Logger
        The logger to use for logging.
    function : Callable
        The function that was run in the future. Used for logging.
    Returns
    -------
    Any
        The result of the future, or None if an exception occurred.
    """
    try:
        result = future.result()
        if isinstance(result, str):
            logger.error(f"Future error running function `{function}`: {result}" )
            sys.exit(result)
        else:
            return result
    except Exception as e:
        
        msg = f"Exception running function `{function}`: {e}\n"
        msg += f"Error traceback is as follows:\n{traceback.format_exc()}"
        logger.error(msg)
        
        sys.exit(msg)
        
def woden_worker_into_queue(q : Queue,
                            all_loaded_python_sources : List[List[Source_Python]],
                            all_loaded_sources_orders : List[int],
                            round_num : int,
                            woden_settings_python : Woden_Settings_Python,
                            array_layout_python : Array_Layout_Python,
                            visi_sets_python : List[Visi_Set_Python],
                            beamtype : int,
                            args : argparse.Namespace,
                            logger : Logger):
    """
    This function launches the :func:`woden_worker` function in a separate process,
    putting outputs into the queue `q`. Specifically, it places
    
     - visi_set_python
     - completed_round
     
    into the queue. 
    
    Do this to keep the casacore version linked to `libwoden_*.so` separate from
    the one linked to `python-casacore`. When both casacore versions are initialised
    in the main process, they clash, and can cause segfaults.
    
    This wrapper is only intended for use when running WODEN in serial mode. As
    such, we pass 0 as the thread index and `visi_sets_python[0,:]`
    into :func:`woden_worker`, as there should only be one thread's worth of
    visibility sets.
    
    Parameters
    ----------
    q : Queue
        The queue to put the outputs into.
    all_loaded_python_sources : List[List[Source_Python]]
        A list of lists of `Source_Python` to be processed, as ouput by `read_skymodel_worker`,
        where each list of `Source_Python` matches a chunk_map in `chunked_skymodel_map_set`.
        Each entry is a list itself as the Shapelet chunking can result in multiple sources
        per thread per round, as the chunking happens by sky direction as well as number of
        shapelet coefficients.
    all_loaded_sources_orders : List[int]
        The order of the sources in `all_loaded_python_sources` as matched to `chunked_skymodel_map_set`;
        orders also output by `read_skymodel_worker` (different threads finish in different times so
        the order of the sources in `all_loaded_python_sources` may not match the order of the chunk_maps)
    round_num : int
        The round number of the processing
    woden_settings_python : Woden_Settings_Python
        The WODEN settings to be used.
    array_layout_python : Array_Layout_Python
        The array layout to be used.
    visi_sets_python : List[Visi_Set_Python]
        The exisiting visibility sets information. The results of this round
        of processing will be added to this and put in the queue; `visi_sets_python`
        will then need to be updated in turn with this `visi_sets_python` in
        the main process.
    beamtype : int
        The value of the type of beam being used, e.g. BeamTypes.FEE_BEAM.value
    args : argparse.Namespace
        Command-line arguments namespace as returned by :func:`get_parser`.
    logger : Logger
        The logger to use for logging.
    """
    
    visi_set_python, _, completed_round = woden_worker(0, all_loaded_python_sources,
                                                       all_loaded_sources_orders,
                                                       round_num,
                                                       woden_settings_python,
                                                       array_layout_python, 
                                                       visi_sets_python[0,:],
                                                       beamtype, args.log_level,
                                                       args.precision,
                                                       profile=args.profile,
                                                       logger=logger)
    
    q.put((visi_set_python, completed_round))
    return
    

def run_woden_processing(num_threads, num_rounds, chunked_skymodel_map_sets,
                 lsts, latitudes, args, beamtype,
                 main_table, shape_table, v_table, q_table, u_table, p_table,
                 woden_settings_python, array_layout_python, visi_sets_python,
                 logger = False, serial_mode = False, uvbeam_objs = None):
    """
    This function runs the WODEN processing, either in serial or parallel mode.
    It sets up the WODEN settings, array layout, and visibility set arrays,
    and then runs the WODEN processing in either serial or parallel mode.
    It also handles the logging and profiling of the processing.
    Parameters
    ----------
    num_threads : int
        The number of threads to use for processing.
    num_rounds : int
        The number of rounds of processing to perform.
    chunked_skymodel_map_sets : List[Skymodel_Chunk_Map]
        A list of chunked skymodel map sets.
    lsts : np.ndarray
        Local Sidereal Times (LSTs) for all time steps.
    latitudes : np.ndarray
        Latitudes of the array at all time steps; when doing array precession,
        these will all be slightly different.
    args : argparse.Namespace
        Command-line arguments namespace.
    beamtype : int
        The value of the type of beam being used, e.g. BeamTypes.FEE_BEAM.value
    main_table : Table
        Main table containing skymodel data.
    shape_table : Table
        Table containing Shapelet data.
    v_table : Table
        Table containing V polarisation data (if not needed, can be boolean False)
    q_table : Table
        Table containing Q polarisation data (if not needed, can be boolean False)
    u_table : Table
        Table containing U polarisation data (if not needed, can be boolean False)
    p_table : Table
        Table containing P polarisation data (if not needed, can be boolean False)
    woden_settings_python : Woden_Settings_Python
        The WODEN settings to be used.
    array_layout_python : Array_Layout_Python
        The array layout to be used.
    visi_sets_python : List[List[Visi_Set_Python]]
        Place to output the visibilities to. Should be of shape (num_visi_threads, num_bands)
        where num_visi_threads is the number of threads to use for processing
        visibilities, and num_bands is the number of frequency bands.
    logger : Logger
        The logger to use for logging. If not provided, a default logger
        will be created.
    serial_mode : bool
        Whether to run in serial mode (True) or parallel mode (False). Default is False.
    uvbeam_objs : np.ndarray
        The UVBeam objects to use for the primary beam calculations. Only needed
        if the beamtype is in `BeamGroups.uvbeam_beams`.
    """
    
    if not logger:
        logger = simple_logger(args.log_level)
        
    total_n_points, total_n_gauss, total_n_shapes, total_n_shape_coeffs = sum_components_in_chunked_skymodel_map_sets(chunked_skymodel_map_sets)
    total_comps = total_n_points + total_n_gauss + total_n_shapes + total_n_shape_coeffs
    
    
    if serial_mode:
        for round_num in range(num_rounds):
            if woden_settings_python.beamtype in BeamGroups.python_calc_beams:
                logger.info(f"Reading (and calculating primary beam values for) set {round_num} sky models")
            else:
                logger.info(f"Reading set {round_num} sky models")
            
            all_loaded_python_sources = []
            all_loaded_sources_orders = []
            
            for i in range(num_threads):
                outputs = read_skymodel_worker(i + round_num * num_threads,
                                                            num_threads, 
                                                            chunked_skymodel_map_sets,
                                                            lsts, latitudes,
                                                            args, beamtype,
                                                            main_table, shape_table,
                                                            v_table, q_table, u_table, p_table,
                                                            logger, uvbeam_objs)
                ##If everything ran fine, should have two outputs
                if len(outputs) == 2:
                    python_sources, order = outputs
                    all_loaded_python_sources.append(python_sources)
                    all_loaded_sources_orders.append(order)
                    
                    q = Queue()
                    p = Process(target=woden_worker_into_queue,
                                args=(q,
                                    all_loaded_python_sources,
                                    all_loaded_sources_orders,
                                    round_num,
                                    woden_settings_python,
                                    array_layout_python,
                                    visi_sets_python,
                                    beamtype,
                                    args,
                                    logger))
                    p.start()
                    
                    visi_set_python, completed_round = q.get()
                    p.join()
                    p.terminate()
                    q.close()
                ##Otherwise, something went wrong, so log the error
                else:
                    logger.error(f"{outputs}")
                    sys.exit(f"{outputs}")
            
            visi_sets_python[0, :] = visi_set_python
            done_n_points, done_n_gauss, done_n_shapes, done_n_shape_coeffs = sum_components_in_chunked_skymodel_map_sets(chunked_skymodel_map_sets[:completed_round + 1])
            done_comps = done_n_points + done_n_gauss + done_n_shapes + done_n_shape_coeffs
            
            logger.info(f"Have completed {done_comps} of {total_comps} components calcs ({(done_comps/total_comps)*100:.1f}%)")
            logger.debug(f"\t{done_n_points} of {total_n_points} points\n"
                        f"\t{done_n_gauss} of {total_n_gauss} gauss\n"
                        f"\t{done_n_shapes} of {total_n_shapes} shapelets\n"
                        f"\t{done_n_shape_coeffs} of {total_n_shape_coeffs} shape coeffs")
    else:
        if args.cpu_mode:
            num_visi_threads = num_threads
            device = 'CPU'
        else:
            num_visi_threads = 1
            device = 'GPU'
            
        ##Currently, multi-threading doesn't work across UVBeam. So read in
        ##sky model with one thread for now
        if woden_settings_python.beamtype in BeamGroups.uvbeam_beams:
            num_sky_model_threads = 1
        else:
            num_sky_model_threads = num_threads
            
        if device == 'CPU':
        
            with ProcessPoolExecutor(max_workers=num_sky_model_threads) as sky_model_executor, ProcessPoolExecutor(max_workers=num_visi_threads) as visi_executor:
            
                # # Loop through multiple rounds of data reading and calculation
                for round_num in range(num_rounds):
                    if woden_settings_python.beamtype in BeamGroups.python_calc_beams:
                        logger.info(f"Reading (and calculating primary beam values for) set {round_num} sky models")
                    else:
                        logger.info(f"Reading set {round_num} sky models")
                        
                        
                    future_data_sky = [sky_model_executor.submit(read_skymodel_worker,
                                                i + round_num * num_threads, num_threads,
                                                chunked_skymodel_map_sets,
                                                lsts, latitudes,
                                                args, beamtype,
                                                main_table, shape_table,
                                                v_table, q_table, u_table, p_table,
                                                logger, uvbeam_objs)
                                        for i in range(num_threads)]
                    
                    all_loaded_python_sources = [0]*num_threads
                    all_loaded_sources_orders = [0]*num_threads
                    for future in future_data_sky:
                        python_sources, order = get_future_result(future, logger, "read_skymodel_worker")
                        
                        all_loaded_python_sources[order] = python_sources
                        all_loaded_sources_orders[order] = order
                        
                    future_data_visi = [visi_executor.submit(woden_worker, i,
                                                [all_loaded_python_sources[i]], #all_loaded_python_sources
                                                [all_loaded_sources_orders[i]], #all_loaded_sources_orders
                                                round_num,
                                                woden_settings_python,
                                                array_layout_python, 
                                                visi_sets_python[i, :],
                                                beamtype, args.log_level,
                                                args.precision, logger,
                                                args.profile)
                                        for i in range(num_visi_threads)]
                    
                    completed = 0
                    
                    for future in concurrent.futures.as_completed(future_data_visi):
                        visi_set_python, thread_ind, completed_round = get_future_result(future, logger, "woden_worker")
                        
                        visi_sets_python[thread_ind, :] = visi_set_python
                        
                        completed += 1
                        
                        if completed == num_visi_threads:
                        
                            done_n_points, done_n_gauss, done_n_shapes, done_n_shape_coeffs = sum_components_in_chunked_skymodel_map_sets(chunked_skymodel_map_sets[:completed_round + 1])
                            done_comps = done_n_points + done_n_gauss + done_n_shapes + done_n_shape_coeffs
                            logger.info(f"Have completed {done_comps} of {total_comps} components calcs ({(done_comps/total_comps)*100:.1f}%)")
                            logger.debug(f"\t{done_n_points} of {total_n_points} points\n"
                                    f"\t{done_n_gauss} of {total_n_gauss} gauss\n"
                                    f"\t{done_n_shapes} of {total_n_shapes} shapelets\n"
                                    f"\t{done_n_shape_coeffs} of {total_n_shape_coeffs} shape coeffs")
                        
        else:
                        
            with concurrent.futures.ProcessPoolExecutor(max_workers=num_sky_model_threads) as sky_model_executor:
                
                visi_executor = concurrent.futures.ProcessPoolExecutor(max_workers=1)
                gpu_calc = None  # To store the future of the calculation thread
                
            # # Loop through multiple rounds of data reading and calculation
                for round_num in range(num_rounds):
                    if woden_settings_python.beamtype in BeamGroups.python_calc_beams:
                        logger.info(f"Reading (and calculating primary beam values for) set {round_num} sky models")
                    else:
                        logger.info(f"Reading set {round_num} sky models")
                    future_data_sky = [sky_model_executor.submit(read_skymodel_worker,
                                                i + round_num * num_threads, num_threads,
                                                chunked_skymodel_map_sets,
                                                lsts, latitudes,
                                                args, beamtype,
                                                main_table, shape_table,
                                                v_table, q_table, u_table, p_table,
                                                logger, uvbeam_objs)
                                        for i in range(num_threads)]
                    
                    all_loaded_python_sources = []
                    all_loaded_sources_orders = []
                    for future in concurrent.futures.as_completed(future_data_sky):
                        python_sources, order = get_future_result(future, logger, "read_skymodel_worker")
                        all_loaded_python_sources.append(python_sources)
                        all_loaded_sources_orders.append(order)
                        
                    # Wait for previous calculation to complete (if there was one)
                    if gpu_calc is not None:
                        # visi_set_python, _, completed_round = gpu_calc.result()
                        visi_set_python, _, completed_round = get_future_result(gpu_calc, logger, "woden_worker")
                        
                        visi_sets_python[0, :] = visi_set_python
                        done_n_points, done_n_gauss, done_n_shapes, done_n_shape_coeffs = sum_components_in_chunked_skymodel_map_sets(chunked_skymodel_map_sets[:completed_round + 1])
                        done_comps = done_n_points + done_n_gauss + done_n_shapes + done_n_shape_coeffs
                        
                        logger.info(f"Have completed {done_comps} of {total_comps} components calcs ({(done_comps/total_comps)*100:.1f}%)")
                        logger.debug(f"\t{done_n_points} of {total_n_points} points\n"
                                    f"\t{done_n_gauss} of {total_n_gauss} gauss\n"
                                    f"\t{done_n_shapes} of {total_n_shapes} shapelets\n"
                                    f"\t{done_n_shape_coeffs} of {total_n_shape_coeffs} shape coeffs")

                    gpu_calc = visi_executor.submit(woden_worker, 0,
                                                        all_loaded_python_sources,
                                                        all_loaded_sources_orders,
                                                        round_num,
                                                        woden_settings_python,
                                                        array_layout_python, 
                                                        visi_sets_python[0,:],
                                                        beamtype, args.log_level,
                                                        args.precision, logger,
                                                        args.profile)
                    
                # # Wait for the final calculation to complete
                if gpu_calc is not None:
                    # visi_set_python, _, completed_round = gpu_calc.result()
                    visi_set_python, _, completed_round = get_future_result(gpu_calc, logger, "woden_worker")
                    visi_sets_python[0, :] = visi_set_python
                    done_n_points, done_n_gauss, done_n_shapes, done_n_shape_coeffs = sum_components_in_chunked_skymodel_map_sets(chunked_skymodel_map_sets[:completed_round + 1])
                    done_comps = done_n_points + done_n_gauss + done_n_shapes + done_n_shape_coeffs
                    
                    logger.info(f"Have completed {done_comps} of {total_comps} components calcs ({(done_comps/total_comps)*100:.1f}%)")
                    logger.debug(f"\t{done_n_points} of {total_n_points} points\n"
                                    f"\t{done_n_gauss} of {total_n_gauss} gauss\n"
                                    f"\t{done_n_shapes} of {total_n_shapes} shapelets\n"
                                    f"\t{done_n_shape_coeffs} of {total_n_shape_coeffs} shape coeffs")
            
        logger.info("Finished all rounds of processing")
        
    return visi_sets_python

@profile
def main(argv=None):
    """Runs the WODEN simulator, given command line inputs. Does fancy things
    like reading in the sky model and running the GPU code in parallel; the
    sky model lazy load allows us to simulate sky models that cannot fit into
    RAM.

    Parameters
    ----------
    argv : _type_, optional
        Will be parsed from the command line if run in main script. Can be
        passed in explicitly for easy testing
    """
    
    woden_start = time()
    
    ##Find out what git/release version we are using, and where the code lives
    gitlabel = get_code_version()
    
    if not argv:
        argv = sys.argv[1:]

    if '--version' in argv:
        print(f"You are using wodenpy version/git hash: {gitlabel}")
        sys.exit()
    
    ##Grab the parser and parse some args
    parser = get_parser()
    args = parser.parse_args(argv)
    
    
    
    ##Check that the input arguments make sense
    args = check_args(args)
    
    if args.num_threads == 1:
        serial_mode = True
    else:
        serial_mode = False

    main_logger = set_woden_logger(args.log_level, args.log_file_name)
    
    set_logger_header(main_logger, gitlabel)
    summarise_input_args(main_logger, args)
    
    lst_deg, gst0_deg, degpdy, ut1utc = get_uvfits_date_and_position_constants(latitude=args.latitude, longitude=args.longitude,
                        height=args.array_height, date=args.date)

    jd_midnight, float_jd = calc_jdcal(args.date)
    jd_date = jd_midnight + float_jd

    ##Generate the woden_struct_classes, which are used to mimic the C
    ##structs. Must be made dynamically, as the precision can change
    woden_struct_classes = Woden_Struct_Classes(args.precision)
    
    num_baselines = int(((args.num_antennas - 1)*args.num_antennas) / 2)

    if args.do_autos:
        num_visis = args.num_time_steps*args.num_freq_channels*(num_baselines + args.num_antennas)
    else:
        num_visis = args.num_time_steps*args.num_freq_channels*num_baselines

    woden_settings_python = fill_woden_settings_python(args, jd_date, lst_deg)
    
    ##fill the lst and mjds fields, precessing if necessary
    lsts, latitudes = setup_lsts_and_phase_centre(woden_settings_python, main_logger)
    
    ##Precessing the array can be expensive if there are 1000s of time steps
    ##so only calculate it if we need to
    if args.dry_run:
        array_layout_python = setup_array_layout_python(woden_settings_python, args)
    else:
    ##calculate the array layout
        array_layout_python = calc_XYZ_diffs(woden_settings_python, args, main_logger)
    
    ##report what beam type we are using
    log_chosen_beamtype(main_logger, woden_settings_python, args)

    # ##read in and chunk the sky model=======================================
    main_logger.info("Doing the initial mapping of sky model")
    t_before = time()
    
    comp_counter = read_radec_count_components(args.cat_filename)
    t_after = time()
    main_logger.info(f"Sky model mapping took {(t_after - t_before)/60.0:.1f} mins")
    
    crop_by_component = True
    if args.sky_crop_sources:
        crop_by_component = False
        
    comp_counter = crop_below_horizon(woden_settings_python.lsts[0],
                                        woden_settings_python.latitude,
                                        comp_counter, main_logger,
                                        crop_by_component=crop_by_component)
    
    main_logger.debug(comp_counter.info_string())
    
    if serial_mode:
        num_threads = 1
    else:
        num_threads = args.num_threads
    
    if args.cpu_mode:
        max_chunks_per_set=num_threads
    else:
        max_chunks_per_set=32
        
    # max_chunks_per_set=args.num_threads*3
        
    chunked_skymodel_map_sets = create_skymodel_chunk_map(comp_counter,
                                                    args.chunking_size,
                                                    woden_settings_python.num_baselines,
                                                    args.num_freq_channels,
                                                    args.num_time_steps,
                                                    num_threads=num_threads,
                                                    max_chunks_per_set=max_chunks_per_set,
                                                    max_dirs=args.max_sky_directions,
                                                    beamtype_value=woden_settings_python.beamtype)
    
    if args.dry_run:
        ##User only wants to check whether the arguments would have worked or not
        main_logger.info("As --dry-run was selected, exiting now before calculating visibilities")
    else:
        
        woden_lib_path = importlib_resources.files(wodenpy).joinpath(f"libwoden_{args.precision}.so")
        main_logger.info(f"Will load libwoden from {woden_lib_path}")
        
        ##Create the shapelet basis functions
        sbf = create_sbf(precision=args.precision)

        ##Create an array of visibility_sets, which get fed into run_woden
        ##and store the output visibilities
        visi_set_array = setup_visi_set_array(woden_struct_classes.Visi_Set,
                                            len(args.band_nums), num_visis,
                                            precision=args.precision)
        
        if woden_settings_python.do_gpu or serial_mode:
            num_visi_threads = 1
        else:
            num_visi_threads = args.num_threads
        
        ##Something to hold all visibility outputs for all threads, for all bands
        visi_sets_python = [[Visi_Set_Python() for _ in range(len(args.band_nums))] for _ in range(num_visi_threads)]
        visi_sets_python = np.array(visi_sets_python)
        
        ##initialise to zeros
        for thread_ind in range(num_visi_threads):
            for band_ind in range(len(args.band_nums)):
                visi_sets_python[thread_ind, band_ind].sum_visi_XX_real = np.zeros(num_visis)
                visi_sets_python[thread_ind, band_ind].sum_visi_XX_imag = np.zeros(num_visis)
                visi_sets_python[thread_ind, band_ind].sum_visi_XY_real = np.zeros(num_visis)
                visi_sets_python[thread_ind, band_ind].sum_visi_XY_imag = np.zeros(num_visis)
                visi_sets_python[thread_ind, band_ind].sum_visi_YX_real = np.zeros(num_visis)
                visi_sets_python[thread_ind, band_ind].sum_visi_YX_imag = np.zeros(num_visis)
                visi_sets_python[thread_ind, band_ind].sum_visi_YY_real = np.zeros(num_visis)
                visi_sets_python[thread_ind, band_ind].sum_visi_YY_imag = np.zeros(num_visis)
                visi_sets_python[thread_ind, band_ind].us_metres = np.zeros(num_visis)
                visi_sets_python[thread_ind, band_ind].vs_metres = np.zeros(num_visis)
                visi_sets_python[thread_ind, band_ind].ws_metres = np.zeros(num_visis)
            

        ##Setup UVBeam objects; only have to do this once, so do it
        ##before looping over all the sky model chunks
        ##TODO pull this out and whack into a function in wodenpy.primary_beams.use_uvbeam.py
        
        if woden_settings_python.beamtype in BeamGroups.uvbeam_beams:
            main_logger.info("Creating UVBeam objects. This may take a while")
            
            band_num = args.band_nums[0]
            base_band_freq = ((band_num - 1)*float(args.coarse_band_width)) + args.lowest_channel_freq
            top_freq = base_band_freq + args.num_freq_channels*args.freq_res
            freqs = np.array([base_band_freq, top_freq])
            
            if woden_settings_python.beamtype == BeamTypes.UVB_MWA.value:
            
                if args.use_MWA_dipamps:
                    dipole_amps = woden_settings_python.mwa_dipole_amps
                else:
                    dipole_amps = np.ones(32)
                
                uvbeam_objs = setup_MWA_uvbeams(args.hdf5_beam_path, freqs,
                                                woden_settings_python.FEE_ideal_delays[:16],
                                                dipole_amps, pixels_per_deg = 5)
                
            elif woden_settings_python.beamtype == BeamTypes.UVB_HERA.value:
                
                ##We only give args the attribute cst_paths if we are
                ##using CST files, so check for that and run accordingly
                if hasattr(args, 'cst_paths'):
                
                    uvbeam_objs = setup_HERA_uvbeams_from_CST(args.cst_paths,
                                                              args.cst_freqs,
                                                              main_logger)
                else:
                    uvbeam_objs = setup_HERA_uvbeams_from_single_file(args.uvbeam_file_path,
                                                                      main_logger)
            
            main_logger.info("UVBeam objects have been initialised")
            
        else:
            uvbeam_objs = False

        ###---------------------------------------------------------------------
        ### heavy lifting area - here we setup running the sky model reading
        ### and running CPU/GPU code at the same time. Means we can limit the
        ### absolute amount of RAM used, and save time doing IO at the same
        ### time as compute
        
        # num_threads = args.num_threads
        
        num_rounds = len(chunked_skymodel_map_sets)
        main_table, shape_table, v_table, q_table, u_table, p_table = get_skymodel_tables(args.cat_filename)
        
        msg = ""
        
        if serial_mode:
            para_mode = "serial"
        else:
            para_mode = "parallel"
        
        if args.cpu_mode:
            if woden_settings_python.beamtype in BeamGroups.uvbeam_beams:
                msg = f"Running in {para_mode} on CPU mode with {num_threads} threads\n"
                msg += "Will read sky model using one thread (UVBeam models seem to do their own threading, trying to thread leads to GIL issues)"
            else:
                msg = f"Running in {para_mode} on CPU mode with {num_threads} threads"
        else:
            if woden_settings_python.beamtype in BeamGroups.uvbeam_beams:
                msg = f"Running in {para_mode} on GPU.\n"
                msg += "Will read sky model using one thread (UVBeam models seem to do their own threading, trying to thread leads to GIL issues)"
            else:
                msg = f"Running in {para_mode} on GPU.\nWill read sky model using {num_threads} threads"
        msg += f"\nThere are {num_rounds} sets of sky models to run"
        
        main_logger.info(msg)
        
        run_woden_processing(num_threads, num_rounds, chunked_skymodel_map_sets,
            lsts, latitudes, args, woden_settings_python.beamtype,
            main_table, shape_table, v_table, q_table, u_table, p_table,
            woden_settings_python, array_layout_python, visi_sets_python,
            main_logger, serial_mode, uvbeam_objs)
            
        ### we've now calculated all the visibilities
        ###---------------------------------------------------------------------
        
        visi_sets_python_combined = visi_sets_python[0, :]
        
        for thread_ind in range(1, num_visi_threads):
            for band_ind in range(len(args.band_nums)):
                
                visi_sets_python_combined[band_ind].sum_visi_XX_real += visi_sets_python[thread_ind, band_ind].sum_visi_XX_real
                visi_sets_python_combined[band_ind].sum_visi_XX_imag += visi_sets_python[thread_ind, band_ind].sum_visi_XX_imag
                visi_sets_python_combined[band_ind].sum_visi_XY_real += visi_sets_python[thread_ind, band_ind].sum_visi_XY_real
                visi_sets_python_combined[band_ind].sum_visi_XY_imag += visi_sets_python[thread_ind, band_ind].sum_visi_XY_imag
                visi_sets_python_combined[band_ind].sum_visi_YX_real += visi_sets_python[thread_ind, band_ind].sum_visi_YX_real
                visi_sets_python_combined[band_ind].sum_visi_YX_imag += visi_sets_python[thread_ind, band_ind].sum_visi_YX_imag
                visi_sets_python_combined[band_ind].sum_visi_YY_real += visi_sets_python[thread_ind, band_ind].sum_visi_YY_real
                visi_sets_python_combined[band_ind].sum_visi_YY_imag += visi_sets_python[thread_ind, band_ind].sum_visi_YY_imag
        
        ##I think we want to X,Y,Z to be in the current frame for writing
        ##out to the uvfits, so calculate again
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

        ##Loop over coarse frequency band and convert visibilities output
        ##from woden.c into uvfits files
        for band_ind, band in enumerate(args.band_nums):

            output_uvfits_name = args.output_uvfits_prepend + '_band%02d.uvfits' %band

            band_low_freq = args.lowest_channel_freq + (band - 1)*args.coarse_band_width
            central_freq_chan_value = band_low_freq + central_freq_chan*args.freq_res

            uus,vvs,wws,v_container = load_visibility_set(visibility_set=visi_sets_python_combined[band_ind],
                                                          num_baselines=num_baselines,
                                                          num_freq_channels=args.num_freq_channels,
                                                          num_time_steps=args.num_time_steps,
                                                          precision=args.precision,
                                                          do_autos=args.do_autos,
                                                          num_ants=args.num_antennas)
            
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
                        array_height=args.array_height,
                        ant_names=args.ant_names)

            baselines_array, date_array = make_baseline_date_arrays(args.num_antennas,
                                        args.date, args.num_time_steps, args.time_res,
                                        do_autos=args.do_autos)

            create_uvfits(v_container=v_container, freq_cent=central_freq_chan_value,
                        ra_point=args.ra0, dec_point=args.dec0,
                        output_uvfits_name=output_uvfits_name,
                        uu=uus, vv=vvs, ww=wws, baselines_array=baselines_array,
                        date_array=date_array,
                        central_freq_chan=central_freq_chan,
                        ch_width=args.freq_res, jd_midnight=jd_midnight,
                        hdu_ant=hdu_ant, gitlabel=gitlabel,
                        longitude=args.longitude, latitude=args.latitude,
                        array_height=args.array_height,
                        telescope_name=args.telescope_name,
                        IAU_order=args.IAU_order,
                        comment=args.command)
            
    if args.pointed_ms_file_name and not args.dry_run:
        main_logger.info(f"Deleting {args.pointed_ms_file_name}...")
        shutil.rmtree(args.pointed_ms_file_name)
        main_logger.info(f"Done")
        
            
    woden_end = time()
    time_passed = timedelta(seconds=woden_end - woden_start)
    main_logger.info(f"Full run took {time_passed}")
    main_logger.info("WODEN is done. Closing the log. S'later")
    
    ##Make sure we clean out the handlers. If we don't, when running
    ##integration tests, old messages hang around and it gets confusing as
    for handler in main_logger.handlers:
        handler.close()
        main_logger.removeHandler(handler)
    logging.shutdown()
    
if __name__ == "__main__":
    main()