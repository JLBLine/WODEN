#!/usr/bin/env python3
"""Wrapper script to run the GPU WODEN code. Author: J.L.B. Line
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
from typing import Union, Tuple, List
from astropy.table import Table
from astropy.io import fits
from wodenpy.use_libwoden.beam_settings import BeamTypes, BeamGroups
from wodenpy.use_libwoden.woden_settings import fill_woden_settings_python, setup_lsts_and_phase_centre, Woden_Settings_Python, convert_woden_settings_to_ctypes
from wodenpy.use_libwoden.visibility_set import setup_visi_set_array, load_visibility_set
from wodenpy.wodenpy_setup.run_setup import get_parser, check_args, get_code_version
from wodenpy.use_libwoden.use_libwoden import load_in_run_woden
from wodenpy.observational.calc_obs import get_uvfits_date_and_position_constants, calc_jdcal
from wodenpy.skymodel.woden_skymodel import crop_below_horizon
from wodenpy.skymodel.read_skymodel import get_skymodel_tables, read_radec_count_components, create_source_catalogue_from_python_sources
from wodenpy.skymodel.chunk_sky_model import create_skymodel_chunk_map, reshape_chunked_skymodel_map_sets, Skymodel_Chunk_Map
from wodenpy.array_layout.create_array_layout import calc_XYZ_diffs, enh2xyz, Array_Layout_Python, convert_array_layout_to_ctypes
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
from concurrent.futures import ProcessPoolExecutor, as_completed
from line_profiler import profile, LineProfiler
from wodenpy.primary_beam.use_everybeam import run_everybeam
from wodenpy.wodenpy_setup.woden_logger import get_log_callback, set_logger_header, simple_logger, log_chosen_beamtype, summarise_input_args, set_woden_logger
from multiprocessing import Manager, set_start_method
from copy import deepcopy
import multiprocessing
from datetime import timedelta
import traceback

# ##Miiiiight give us a speed up on Linux??
# try:
#     print("WE BE DOING A FORKSERVER")
#     set_start_method("forkserver")  # Alternatives: "spawn" or "forkserver"
# except RuntimeError:
#     # Start method is already set in some environments like Jupyter
#     pass

# import faulthandler
# faulthandler.enable()

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


class Visi_Set_Python(object):
    def __init__(self):
        self.us_metres = None
        self.vs_metres = None
        self.ws_metres = None
        self.sum_visi_XX_real = None
        self.sum_visi_XX_imag = None
        self.sum_visi_XY_real = None
        self.sum_visi_XY_imag = None
        self.sum_visi_YX_real = None
        self.sum_visi_YX_imag = None
        self.sum_visi_YY_real = None
        self.sum_visi_YY_imag = None
        
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
    This function runs WODEN CPU code on a separate thread, processing source catalogues from a queue until the queue is empty.
    
    Parameters
    ----------
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
    run_woden : _NamedFuncPointer
        A pointer to the WODEN GPU function to be run.
    woden_settings : Woden_Settings
        The WODEN settings to be used.
    array_layout : Array_Layout_Ctypes
        The array layout to be used.
    woden_struct_classes : Woden_Struct_Classes
        The WODEN struct classes to be used
    sbf : np.ndarray
        The shapelet basis function array
    beamtype : int
        The value of the type of beam being used, e.g. BeamTypes.FEE_BEAM.value
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
        
        ##Create a ctypes Source_Catalogue from the python sources to feed the GPU
        source_catalogue = create_source_catalogue_from_python_sources(python_sources,
                                                                    woden_struct_classes,
                                                                    beamtype, precision)
        
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
            
        # if profile:
        #     profiler.disable()
        #     profile_filename = f"wod_woden_worker_{os.getpid()}_{round_num:04d}.lprof"
        #     # logger.info("Dumping profile to", profile_filename)
        #     profiler.dump_stats(profile_filename)
            
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
                         logger : Logger = False) -> Tuple[List[Source_Python], int]:
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
    profile : bool, optional
        Whether to profile the function (default is False).
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
                                beamtype,
                                lsts, latitudes,
                                v_table, q_table,
                                u_table, p_table,
                                args.precision)
        
        end = time()
        
        if beamtype in BeamGroups.eb_beam_values:
            logger.info(f"Finshed sky set {set_ind} reading thread num {thread_num} in {end-start:.1f} seconds")
        else:
            logger.debug(f"Finshed sky set {set_ind} reading thread num {thread_num} in {end-start:.1f} seconds")
        
        if args.profile:
            profiler.disable()
            profile_filename = f"wod_read_skymodel_worker_{os.getpid()}_{set_ind:04d}.lprof"
            # logger.info("Dumping profile to", profile_filename)
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

def get_future_result(future, logger):
    try:
        result = future.result()
        if isinstance(result, str):
            logger.error(result)
            sys.exit(result)
        else:
            return result
    except Exception as e:
        msg = f"Exception in getting result from future: {e}"
        logger.error(msg)
        sys.exit(msg)
    

def run_woden_processing(num_threads, num_rounds, chunked_skymodel_map_sets,
                 lsts, latitudes, args, beamtype,
                 main_table, shape_table, v_table, q_table, u_table, p_table,
                 woden_settings_python, array_layout_python, visi_sets_python,
                 logger = False, serial_mode = False):
    
    if not logger:
        logger = simple_logger(args.log_level)
        
    total_n_points, total_n_gauss, total_n_shapes, total_n_shape_coeffs = sum_components_in_chunked_skymodel_map_sets(chunked_skymodel_map_sets)
    total_comps = total_n_points + total_n_gauss + total_n_shapes + total_n_shape_coeffs
    
    
    if serial_mode:
        for round_num in range(num_rounds):
            
            logger.info(f"Reading set {round_num} sky models")
            
            all_loaded_python_sources = []
            all_loaded_sources_orders = []
            
            for i in range(num_threads):
                python_sources, order = read_skymodel_worker(i + round_num * num_threads,
                                                            num_threads, 
                                                            chunked_skymodel_map_sets,
                                                            lsts, latitudes,
                                                            args, beamtype,
                                                            main_table, shape_table,
                                                            v_table, q_table, u_table, p_table,
                                                            logger)
            
                all_loaded_python_sources.append(python_sources)
                all_loaded_sources_orders.append(order)
                
            visi_set_python, _, completed_round = woden_worker(0, 
                                                               all_loaded_python_sources,
                                                               all_loaded_sources_orders,
                                                               round_num,
                                                               woden_settings_python,
                                                               array_layout_python, 
                                                               visi_sets_python[0,:],
                                                               beamtype, args.log_level,
                                                               args.precision,
                                                               profile=args.profile,
                                                               logger=logger)
                
            
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
            
        # logger.info(f"Running WODEN with {num_threads} threads for sky model reading and {num_visi_threads} threads for visibility processing")
        
        if device == 'CPU':
        
            with ProcessPoolExecutor(max_workers=num_threads) as sky_model_executor, ProcessPoolExecutor(max_workers=num_visi_threads) as visi_executor:
            
                # # Loop through multiple rounds of data reading and calculation
                for round_num in range(num_rounds):
                    logger.info(f"Reading set {round_num} sky models")
                    future_data_sky = [sky_model_executor.submit(read_skymodel_worker,
                                                i + round_num * num_threads, num_threads,
                                                chunked_skymodel_map_sets,
                                                lsts, latitudes,
                                                args, beamtype,
                                                main_table, shape_table,
                                                v_table, q_table, u_table, p_table,
                                                logger)
                                        for i in range(num_threads)]
                    
                    all_loaded_python_sources = [0]*num_threads
                    all_loaded_sources_orders = [0]*num_threads
                    for future in future_data_sky:
                        python_sources, order = get_future_result(future, logger)
                        
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
                        visi_set_python, thread_ind, completed_round = get_future_result(future, logger)
                        
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
                        
            with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as sky_model_executor:
                
                visi_executor = concurrent.futures.ThreadPoolExecutor(max_workers=1)
                gpu_calc = None  # To store the future of the calculation thread
                
            # # Loop through multiple rounds of data reading and calculation
                for round_num in range(num_rounds):
                    
                    logger.info(f"Reading set {round_num} sky models")
                    future_data_sky = [sky_model_executor.submit(read_skymodel_worker,
                                                i + round_num * num_threads, num_threads,
                                                chunked_skymodel_map_sets,
                                                lsts, latitudes,
                                                args, beamtype,
                                                main_table, shape_table,
                                                v_table, q_table, u_table, p_table, 
                                                logger)
                                        for i in range(num_threads)]
                    
                    all_loaded_python_sources = []
                    all_loaded_sources_orders = []
                    for future in concurrent.futures.as_completed(future_data_sky):
                        python_sources, order = get_future_result(future, logger)
                        all_loaded_python_sources.append(python_sources)
                        all_loaded_sources_orders.append(order)
                        
                    # Wait for previous calculation to complete (if there was one)
                    if gpu_calc is not None:
                        # visi_set_python, _, completed_round = gpu_calc.result()
                        visi_set_python, _, completed_round = get_future_result(gpu_calc, logger)
                        
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
                    visi_set_python, _, completed_round = get_future_result(gpu_calc, logger)
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
    
    ##Grab the parser and parse some args
    parser = get_parser()
    args = parser.parse_args(argv)
    ##Check that the input arguments make sense
    args = check_args(args)
    
    ##Find out what git/release version we are using, and where the code lives
    gitlabel = get_code_version()
    
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
            msg = f"Running in {para_mode} on CPU mode with {num_threads} threads"
        else:
            msg = f"Running in {para_mode} on GPU.\nWill read sky model using {num_threads} threads"
        msg += f"\nThere are {num_rounds} sets of sky models to run"
        
        main_logger.info(msg)
        
        run_woden_processing(num_threads, num_rounds, chunked_skymodel_map_sets,
            lsts, latitudes, args, woden_settings_python.beamtype,
            main_table, shape_table, v_table, q_table, u_table, p_table,
            woden_settings_python, array_layout_python, visi_sets_python,
            main_logger, serial_mode)
            
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
            
    woden_end = time()
    time_passed = timedelta(seconds=woden_end - woden_start)
    main_logger.info(f"Full run took {time_passed}")
    main_logger.info("WODEN is done. Closing the log. S'later")
    
if __name__ == "__main__":
    main()