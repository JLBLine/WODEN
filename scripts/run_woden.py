#!/usr/bin/env python3
"""Wrapper script to run the GPU WODEN code. Author: J.L.B. Line
"""
import numpy as np
import sys
import os
from time import time
import argparse
from queue import Queue
from threading import Thread
from ctypes import POINTER, c_double, c_float, c_int, create_string_buffer, string_at
from typing import Union, Tuple, List
from astropy.table import Table
from astropy.io import fits
from wodenpy.use_libwoden.beam_settings import BeamTypes
from wodenpy.use_libwoden.woden_settings import create_woden_settings, setup_lsts_and_phase_centre
from wodenpy.use_libwoden.visibility_set import setup_visi_set_array, load_visibility_set
from wodenpy.wodenpy_setup.run_setup import get_parser, check_args, get_code_version
from wodenpy.use_libwoden.use_libwoden import load_in_woden_library
from wodenpy.observational.calc_obs import get_uvfits_date_and_position_constants, calc_jdcal
from wodenpy.skymodel.woden_skymodel import crop_below_horizon
from wodenpy.skymodel.read_skymodel import get_skymodel_tables, read_radec_count_components, create_source_catalogue_from_python_sources
from wodenpy.skymodel.chunk_sky_model import create_skymodel_chunk_map, reshape_chunked_skymodel_map_sets, Skymodel_Chunk_Map
from wodenpy.array_layout.create_array_layout import calc_XYZ_diffs, enh2xyz
from wodenpy.use_libwoden.array_layout_struct import Array_Layout
from wodenpy.uvfits.wodenpy_uvfits import make_antenna_table, make_baseline_date_arrays, create_uvfits
from wodenpy.phase_rotate.remove_phase_track import remove_phase_tracking
from wodenpy.use_libwoden.shapelets import create_sbf
from wodenpy.use_libwoden.create_woden_struct_classes import Woden_Struct_Classes
from wodenpy.skymodel.read_fits_skymodel import read_fits_skymodel_chunks, calc_everybeam_for_components
from wodenpy.use_libwoden.skymodel_structs import setup_source_catalogue, Source_Python
import concurrent.futures
# import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed

from line_profiler import profile, LineProfiler
from wodenpy.primary_beam.use_everybeam import run_everybeam

from copy import deepcopy

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

def _run_run_woden(run_woden, woden_settings : Woden_Settings,
                   visibility_set : Visi_Set, source_catalogue : Source_Catalogue,
                   array_layout : Array_Layout, sbf : np.ndarray):
    """This is silly, but for profiling purposes, we need to wrap the `run_woden`
    C/GPU function in another function so it's easy to pick out when profiling.
    We can search for just `_run_run_woden` in the profiling output to see how
    long we're spending on the GPU (just searching for `run_woden` doesn't
    seem to work, at least in `yappi`.)
    """
    
    run_woden(woden_settings, visibility_set, source_catalogue, array_layout,
                  sbf)
    
    
class Woden_Settings_Python(object):
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
    
    
def convert_woden_settings_to_python(woden_settings_ctypes : Woden_Settings):
    
    woden_settings_python = Woden_Settings_Python()
    
    woden_settings_python.lst_base = woden_settings_ctypes.lst_base
    woden_settings_python.lst_obs_epoch_base = woden_settings_ctypes.lst_obs_epoch_base
    woden_settings_python.ra0 = woden_settings_ctypes.ra0
    woden_settings_python.dec0 = woden_settings_ctypes.dec0
    woden_settings_python.sdec0 = woden_settings_ctypes.sdec0
    woden_settings_python.cdec0 = woden_settings_ctypes.cdec0
    woden_settings_python.num_baselines = woden_settings_ctypes.num_baselines
    woden_settings_python.num_ants = woden_settings_ctypes.num_ants
    woden_settings_python.num_freqs = woden_settings_ctypes.num_freqs
    woden_settings_python.frequency_resolution = woden_settings_ctypes.frequency_resolution
    woden_settings_python.base_low_freq = woden_settings_ctypes.base_low_freq
    woden_settings_python.num_time_steps = woden_settings_ctypes.num_time_steps
    woden_settings_python.time_res = woden_settings_ctypes.time_res
    woden_settings_python.num_bands = woden_settings_ctypes.num_bands
    woden_settings_python.beamtype = woden_settings_ctypes.beamtype
    woden_settings_python.gauss_beam_FWHM = woden_settings_ctypes.gauss_beam_FWHM
    woden_settings_python.gauss_beam_ref_freq = woden_settings_ctypes.gauss_beam_ref_freq
    woden_settings_python.chunking_size = woden_settings_ctypes.chunking_size
    woden_settings_python.jd_date = woden_settings_ctypes.jd_date
    woden_settings_python.latitude = woden_settings_ctypes.latitude
    woden_settings_python.latitude_obs_epoch_base = woden_settings_ctypes.latitude_obs_epoch_base
    woden_settings_python.longitude = woden_settings_ctypes.longitude
    woden_settings_python.coarse_band_width = woden_settings_ctypes.coarse_band_width
    woden_settings_python.gauss_ra_point = woden_settings_ctypes.gauss_ra_point
    woden_settings_python.gauss_dec_point = woden_settings_ctypes.gauss_dec_point
    woden_settings_python.num_cross = woden_settings_ctypes.num_cross
    woden_settings_python.num_autos = woden_settings_ctypes.num_autos
    woden_settings_python.num_visis = woden_settings_ctypes.num_visis
    woden_settings_python.base_band_freq = woden_settings_ctypes.base_band_freq
    woden_settings_python.do_precession = woden_settings_ctypes.do_precession
    woden_settings_python.do_autos = woden_settings_ctypes.do_autos
    woden_settings_python.use_dipamps = woden_settings_ctypes.use_dipamps
    woden_settings_python.single_everybeam_station = woden_settings_ctypes.single_everybeam_station
    woden_settings_python.off_cardinal_dipoles = woden_settings_ctypes.off_cardinal_dipoles
    woden_settings_python.do_gpu = woden_settings_ctypes.do_gpu
    
    woden_settings_python.mjds = np.ctypeslib.as_array(woden_settings_ctypes.mjds, shape=(woden_settings_ctypes.num_time_steps,))
    woden_settings_python.latitudes = np.ctypeslib.as_array(woden_settings_ctypes.latitudes, shape=(woden_settings_ctypes.num_time_steps,))
    woden_settings_python.lsts = np.ctypeslib.as_array(woden_settings_ctypes.lsts, shape=(woden_settings_ctypes.num_time_steps,))
    woden_settings_python.band_nums = np.ctypeslib.as_array(woden_settings_ctypes.band_nums, shape=(woden_settings_ctypes.num_bands,))
    
    if woden_settings_ctypes.beamtype == BeamTypes.FEE_BEAM.value or \
        woden_settings_ctypes.beamtype == BeamTypes.FEE_BEAM_INTERP.value or \
        woden_settings_ctypes.beamtype == BeamTypes.MWA_ANALY.value:
    
        woden_settings_python.FEE_ideal_delays = np.ctypeslib.as_array(woden_settings_ctypes.FEE_ideal_delays, shape=(16,))
        woden_settings_python.hdf5_beam_path = string_at(woden_settings_ctypes.hdf5_beam_path).decode('utf-8')
        print('PYTHON TRANSLATED IT TO', woden_settings_python.hdf5_beam_path)
    
    if woden_settings_ctypes.use_dipamps:
        woden_settings_python.mwa_dipole_amps = np.ctypeslib.as_array(woden_settings_ctypes.mwa_dipole_amps, shape=(32*woden_settings_ctypes.num_ants,))
    
    # woden_settings_python.array_layout_file = woden_settings_ctypes.array_layout_file
    # woden_settings_python.array_layout_file_path = woden_settings_ctypes.array_layout_file_path
    # woden_settings_python.band_nums = woden_settings_ctypes.band_nums
    # woden_settings_python.hdf5_beam_path = woden_settings_ctypes.hdf5_beam_path
    # woden_settings_python.sky_crop_type = woden_settings_ctypes.sky_crop_type
    
    return woden_settings_python

def convert_woden_settings_to_ctypes(woden_settings_python : Woden_Settings_Python,
                                     woden_settings_ctypes : Woden_Settings) -> Woden_Settings:
    
    woden_settings_ctypes.lst_base = woden_settings_python.lst_base
    woden_settings_ctypes.lst_obs_epoch_base = woden_settings_python.lst_obs_epoch_base
    
    # if woden_settings_python.lst_base is None:
    #     pass
    # else:
        
    # if woden_settings_python.lst_obs_epoch_base is None:
    #     pass
    # else:
        
        
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
    woden_settings_ctypes.gauss_beam_FWHM = woden_settings_python.gauss_beam_FWHM
    woden_settings_ctypes.gauss_beam_ref_freq = woden_settings_python.gauss_beam_ref_freq
    woden_settings_ctypes.chunking_size = woden_settings_python.chunking_size
    woden_settings_ctypes.jd_date = woden_settings_python.jd_date
    woden_settings_ctypes.latitude = woden_settings_python.latitude
    woden_settings_ctypes.latitude_obs_epoch_base = woden_settings_python.latitude_obs_epoch_base
    woden_settings_ctypes.longitude = woden_settings_python.longitude
    woden_settings_ctypes.coarse_band_width = woden_settings_python.coarse_band_width
    woden_settings_ctypes.gauss_ra_point = woden_settings_python.gauss_ra_point
    woden_settings_ctypes.gauss_dec_point = woden_settings_python.gauss_dec_point
    woden_settings_ctypes.num_cross = woden_settings_python.num_cross
    woden_settings_ctypes.num_autos = woden_settings_python.num_autos
    woden_settings_ctypes.num_visis = woden_settings_python.num_visis
    woden_settings_ctypes.base_band_freq = woden_settings_python.base_band_freq
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
    
    if woden_settings_ctypes.beamtype == BeamTypes.FEE_BEAM.value or \
        woden_settings_ctypes.beamtype == BeamTypes.FEE_BEAM_INTERP.value or \
        woden_settings_ctypes.beamtype == BeamTypes.MWA_ANALY.value:
    
        woden_settings_ctypes.FEE_ideal_delays = woden_settings_python.FEE_ideal_delays.ctypes.data_as(POINTER(c_int))
        woden_settings_ctypes.hdf5_beam_path = create_string_buffer(woden_settings_python.hdf5_beam_path.encode('utf-8'))
    
    if woden_settings_ctypes.use_dipamps:
        woden_settings_ctypes.mwa_dipole_amps = woden_settings_python.mwa_dipole_amps.ctypes.data_as(POINTER(c_double))
    
    # woden_settings_python.array_layout_file = woden_settings_ctypes.array_layout_file
    # woden_settings_python.array_layout_file_path = woden_settings_ctypes.array_layout_file_path
    
    # woden_settings_python.hdf5_beam_path = woden_settings_ctypes.hdf5_beam_path
    # woden_settings_python.sky_crop_type = woden_settings_ctypes.sky_crop_type
    
    return woden_settings_ctypes


class Array_Layout_Python(object):
    def __init__(self):
        self.ant_X = None
        self.ant_Y = None
        self.ant_Z = None
        self.X_diff_metres = None
        self.Y_diff_metres = None
        self.Z_diff_metres = None
        self.ant_east = None
        self.ant_north = None
        self.ant_height = None
        self.latitude = None
        self.num_baselines = None
        self.num_tiles = None
        self.lst_base = None
        
def convert_array_layout_to_python(array_layout_ctypes : Array_Layout, num_times : int):
    array_layout_python = Array_Layout_Python()
    
    # array_layout_python.ant_X = np.ctypeslib.as_array(array_layout_ctypes.ant_X, shape=(array_layout_ctypes.num_tiles,))
    # array_layout_python.ant_Y = np.ctypeslib.as_array(array_layout_ctypes.ant_Y, shape=(array_layout_ctypes.num_tiles,))
    # array_layout_python.ant_Z = np.ctypeslib.as_array(array_layout_ctypes.ant_Z, shape=(array_layout_ctypes.num_tiles,))
    
    
    num_diffs = array_layout_ctypes.num_baselines * num_times
    
    array_layout_python.X_diff_metres = np.ctypeslib.as_array(array_layout_ctypes.X_diff_metres, shape=(num_diffs,))
    array_layout_python.Y_diff_metres = np.ctypeslib.as_array(array_layout_ctypes.Y_diff_metres, shape=(num_diffs,))
    array_layout_python.Z_diff_metres = np.ctypeslib.as_array(array_layout_ctypes.Z_diff_metres, shape=(num_diffs,))
    
    array_layout_python.latitude = array_layout_ctypes.latitude
    array_layout_python.num_baselines = array_layout_ctypes.num_baselines
    array_layout_python.num_tiles = array_layout_ctypes.num_tiles
    array_layout_python.lst_base = array_layout_ctypes.lst_base
    
    return array_layout_python

def convert_array_layout_to_ctypes(array_layout_python : Array_Layout_Python,
                                   array_layout_ctypes : Array_Layout) -> Array_Layout:
    
    # array_layout_ctypes.ant_X = array_layout_python.ant_X.ctypes.data_as(POINTER(c_double))
    # array_layout_ctypes.ant_Y = array_layout_python.ant_Y.ctypes.data_as(POINTER(c_double))
    # array_layout_ctypes.ant_Z = array_layout_python.ant_Z.ctypes.data_as(POINTER(c_double))
    array_layout_ctypes.X_diff_metres = array_layout_python.X_diff_metres.ctypes.data_as(POINTER(c_double))
    array_layout_ctypes.Y_diff_metres = array_layout_python.Y_diff_metres.ctypes.data_as(POINTER(c_double))
    array_layout_ctypes.Z_diff_metres = array_layout_python.Z_diff_metres.ctypes.data_as(POINTER(c_double))
    
    array_layout_ctypes.latitude = array_layout_python.latitude
    array_layout_ctypes.num_baselines = array_layout_python.num_baselines
    array_layout_ctypes.num_tiles = array_layout_python.num_tiles
    array_layout_ctypes.lst_base = array_layout_python.lst_base
    
    return array_layout_ctypes
    
    
    
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

def woden_thread(all_loaded_python_sources : List[List[Source_Python]],
                 all_loaded_sources_orders : List[int], round_num : int,
                 run_woden, woden_settings_python : Woden_Settings_Python,
                 array_layout_python : Array_Layout_Python,
                 visi_sets_python : List[Visi_Set_Python],
                 woden_struct_classes : Woden_Struct_Classes,
                 beamtype : int,
                 precision : str = "double"):
    """
    This function runs WODEN C/GPU code on a separate thread, processing source catalogues from a queue until the queue is empty.
    
    Parameters
    ----------
    all_loaded_python_sources : List[List[Source_Python]]
        A list of lists of `Source_Python` to be processed, as ouput by `read_skymodel_thread`,
        where each list of `Source_Python` matches a chunk_map in `chunked_skymodel_map_set`.
        Each entry is a list itself as the Shapelet chunking can result in multiple sources
        per thread per round, as the chunking happens by sky direction as well as number of
        shapelet coefficients.
    all_loaded_sources_orders : List[int]
        The order of the sources in `all_loaded_python_sources` as matched to `chunked_skymodel_map_set`;
        orders also output by `read_skymodel_thread` (different threads finish in different times so
        the order of the sources in `all_loaded_python_sources` may not match the order of the chunk_maps)
    round_num : int
        The round number of the processing
    run_woden : _NamedFuncPointer
        A pointer to the WODEN GPU function to be run.
    woden_settings : Woden_Settings
        The WODEN settings to be used.
    array_layout : Array_Layout
        The array layout to be used.
    woden_struct_classes : Woden_Struct_Classes
        The WODEN struct classes to be used
    sbf : np.ndarray
        The shapelet basis function array
    beamtype : int
        The value of the type of beam being used, e.g. BeamTypes.FEE_BEAM.value
    """
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
            
    # print("THIS LENGTH", len(python_sources))
    
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
    
    array_layout = Array_Layout()
    array_layout = convert_array_layout_to_ctypes(array_layout_python, array_layout)
    
    sbf = create_sbf(precision=precision)
    
    if woden_settings.do_gpu:
        print(f"Sending set {round_num} to GPU")
    else:
        print(f"Sending set {round_num} to CPU")
    start = time()
            
    #Run the GPU code in a wrapper function so we can easily profile it
    # _run_run_woden(run_woden, woden_settings, visibility_set,
    #                source_catalogue, array_layout, sbf)
    
    run_woden(woden_settings, visibility_set, source_catalogue, array_layout,
                  sbf)
    
    
    end = time()
    
    if woden_settings.do_gpu:
        print(f"Set {round_num} has returned from the GPU after {end-start:.1f} seconds")
    else:
        print(f"Set {round_num} has returned from the CPU after {end-start:.1f} seconds")
    
    for set_ind in range(woden_settings_python.num_bands):
        visi_set = visi_sets_python[set_ind]
        visi_set.us_metres = deepcopy(np.ctypeslib.as_array(visibility_set[set_ind].us_metres, shape=(woden_settings_python.num_visis,)))
        visi_set.vs_metres = deepcopy(np.ctypeslib.as_array(visibility_set[set_ind].vs_metres, shape=(woden_settings_python.num_visis,)))
        visi_set.ws_metres = deepcopy(np.ctypeslib.as_array(visibility_set[set_ind].ws_metres, shape=(woden_settings_python.num_visis,)))
        visi_set.sum_visi_XX_real += np.ctypeslib.as_array(visibility_set[set_ind].sum_visi_XX_real, shape=(woden_settings_python.num_visis,))
        visi_set.sum_visi_XX_imag += np.ctypeslib.as_array(visibility_set[set_ind].sum_visi_XX_imag, shape=(woden_settings_python.num_visis,))
        visi_set.sum_visi_XY_real += np.ctypeslib.as_array(visibility_set[set_ind].sum_visi_XY_real, shape=(woden_settings_python.num_visis,))
        visi_set.sum_visi_XY_imag += np.ctypeslib.as_array(visibility_set[set_ind].sum_visi_XY_imag, shape=(woden_settings_python.num_visis,))
        visi_set.sum_visi_YX_real += np.ctypeslib.as_array(visibility_set[set_ind].sum_visi_YX_real, shape=(woden_settings_python.num_visis,))
        visi_set.sum_visi_YX_imag += np.ctypeslib.as_array(visibility_set[set_ind].sum_visi_YX_imag, shape=(woden_settings_python.num_visis,))
        visi_set.sum_visi_YY_real += np.ctypeslib.as_array(visibility_set[set_ind].sum_visi_YY_real, shape=(woden_settings_python.num_visis,))
        visi_set.sum_visi_YY_imag += np.ctypeslib.as_array(visibility_set[set_ind].sum_visi_YY_imag, shape=(woden_settings_python.num_visis,))
        
    return visi_sets_python


def woden_multithread(thread_ind : int,
                      all_loaded_python_sources : List[List[Source_Python]],
                      all_loaded_sources_orders : List[int], round_num : int,
                      woden_settings_python : Woden_Settings_Python,
                      array_layout_python : Array_Layout_Python,
                      visi_sets_python : List[Visi_Set_Python],
                      beamtype : int,
                      precision : str = "double"):
    """
    This function runs WODEN C/GPU code on a separate thread, processing source catalogues from a queue until the queue is empty.
    
    Parameters
    ----------
    all_loaded_python_sources : List[List[Source_Python]]
        A list of lists of `Source_Python` to be processed, as ouput by `read_skymodel_thread`,
        where each list of `Source_Python` matches a chunk_map in `chunked_skymodel_map_set`.
        Each entry is a list itself as the Shapelet chunking can result in multiple sources
        per thread per round, as the chunking happens by sky direction as well as number of
        shapelet coefficients.
    all_loaded_sources_orders : List[int]
        The order of the sources in `all_loaded_python_sources` as matched to `chunked_skymodel_map_set`;
        orders also output by `read_skymodel_thread` (different threads finish in different times so
        the order of the sources in `all_loaded_python_sources` may not match the order of the chunk_maps)
    round_num : int
        The round number of the processing
    run_woden : _NamedFuncPointer
        A pointer to the WODEN GPU function to be run.
    woden_settings : Woden_Settings
        The WODEN settings to be used.
    array_layout : Array_Layout
        The array layout to be used.
    woden_struct_classes : Woden_Struct_Classes
        The WODEN struct classes to be used
    sbf : np.ndarray
        The shapelet basis function array
    beamtype : int
        The value of the type of beam being used, e.g. BeamTypes.FEE_BEAM.value
    """
    
    woden_struct_classes = Woden_Struct_Classes(precision)
    
    ##Depending on what precision was selected by the user, load in the
    ##C/CUDA library and return the `run_woden` function
    run_woden = load_in_woden_library(woden_struct_classes)
    
    
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
            
    # ##We only use one GPU, so send all the sources to the GPU
    # if woden_settings_python.do_gpu:
    #     sources_to_use = python_sources
    # ##otherwise, only send the sources to the CPU that are meant for this thread
    # else:
    #     sources_to_use = [python_sources[thread_ind]]
    
    sources_to_use = python_sources
    
    ##Create a ctypes Source_Catalogue from the python sources to feed the GPU
    source_catalogue = create_source_catalogue_from_python_sources(sources_to_use,
                                                                   woden_struct_classes,
                                                                   beamtype, precision)
    
    woden_settings = woden_struct_classes.Woden_Settings()
    woden_settings = convert_woden_settings_to_ctypes(woden_settings_python, woden_settings)
    
    visibility_set = setup_visi_set_array(woden_struct_classes.Visi_Set,
                                              woden_settings.num_bands,
                                              woden_settings.num_visis,
                                              precision=precision)
    
    array_layout = Array_Layout()
    array_layout = convert_array_layout_to_ctypes(array_layout_python, array_layout)
    
    sbf = create_sbf(precision=precision)
    
    if woden_settings.do_gpu:
        print(f"Sending Sky {thread_ind} set {round_num} to GPU")
    else:
        print(f"Sending Sky {thread_ind} set {round_num} to CPU")
    start = time()
            
    run_woden(woden_settings, visibility_set, source_catalogue, array_layout,
                  sbf)
    end = time()
    
    if woden_settings.do_gpu:
        print(f"Sky {thread_ind} Set {round_num} has returned from the GPU after {end-start:.1f} seconds")
    else:
        print(f"Sky {thread_ind} Set {round_num} has returned from the CPU after {end-start:.1f} seconds")
        
        
    # maybe_visi_set_python = [Visi_Set_Python() for band in range(woden_settings_python.num_bands)]
    
    for band_ind in range(woden_settings_python.num_bands):
        # visi_set = visi_sets_python[band_ind]
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
        
        # maybe_visi_set_python[band_ind].us_metres = deepcopy(np.ctypeslib.as_array(visibility_set[band_ind].us_metres, shape=(woden_settings_python.num_visis,)))
        # maybe_visi_set_python[band_ind].vs_metres = deepcopy(np.ctypeslib.as_array(visibility_set[band_ind].vs_metres, shape=(woden_settings_python.num_visis,)))
        # maybe_visi_set_python[band_ind].ws_metres = deepcopy(np.ctypeslib.as_array(visibility_set[band_ind].ws_metres, shape=(woden_settings_python.num_visis,)))
        
        # maybe_visi_set_python[band_ind].sum_visi_XX_real = deepcopy(np.ctypeslib.as_array(visibility_set[band_ind].sum_visi_XX_real, shape=(woden_settings_python.num_visis,)))
        # maybe_visi_set_python[band_ind].sum_visi_XX_imag = deepcopy(np.ctypeslib.as_array(visibility_set[band_ind].sum_visi_XX_imag, shape=(woden_settings_python.num_visis,)))
        # maybe_visi_set_python[band_ind].sum_visi_XY_real = deepcopy(np.ctypeslib.as_array(visibility_set[band_ind].sum_visi_XY_real, shape=(woden_settings_python.num_visis,)))
        # maybe_visi_set_python[band_ind].sum_visi_XY_imag = deepcopy(np.ctypeslib.as_array(visibility_set[band_ind].sum_visi_XY_imag, shape=(woden_settings_python.num_visis,)))
        # maybe_visi_set_python[band_ind].sum_visi_YX_real = deepcopy(np.ctypeslib.as_array(visibility_set[band_ind].sum_visi_YX_real, shape=(woden_settings_python.num_visis,)))
        # maybe_visi_set_python[band_ind].sum_visi_YX_imag = deepcopy(np.ctypeslib.as_array(visibility_set[band_ind].sum_visi_YX_imag, shape=(woden_settings_python.num_visis,)))
        # maybe_visi_set_python[band_ind].sum_visi_YY_real = deepcopy(np.ctypeslib.as_array(visibility_set[band_ind].sum_visi_YY_real, shape=(woden_settings_python.num_visis,)))
        # maybe_visi_set_python[band_ind].sum_visi_YY_imag = deepcopy(np.ctypeslib.as_array(visibility_set[band_ind].sum_visi_YY_imag, shape=(woden_settings_python.num_visis,)))
        
        
    # return visi_sets_python, thread_ind
    return visi_sets_python, thread_ind

def read_skymodel_thread(thread_id : int, num_threads : int,
                         chunked_skymodel_map_sets: List[Skymodel_Chunk_Map],
                         lsts : np.ndarray, latitudes : np.ndarray,
                         args : argparse.Namespace,
                         beamtype : int,
                         main_table : Table, shape_table : Table,
                         v_table : Table = False, q_table : Table = False,
                         u_table : Table = False, p_table : Table = False) -> Tuple[List[Source_Python], int]:
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
        
    print(f"From Set {set_ind} thread num {thread_num} reading {n_p} points, {n_g} gauss, {n_s} shape, {n_c} shape coeffs using thread id {thread_id}")
                
    if len(chunk_maps) == 0:
        print(f"Thread {thread_id} has no work to do")
        return [], thread_num
    
    
    if args.profile:
        profiler = LineProfiler()
        ##Add whatever functions that are called in `read_fits_skymodel_chunks`
        ##here to profile them
        profiler.add_function(read_fits_skymodel_chunks)
        profiler.add_function(calc_everybeam_for_components)
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
    
    print(f"Finshed sky reading thread {thread_id} in {end-start:.1f} seconds")
    
    if args.profile:
        profiler.disable()
        profile_filename = f"line_profile_{os.getpid()}.lprof"
        print("Dumping profile to", profile_filename)
        profiler.dump_stats(profile_filename)
    
    return python_sources, thread_num, set_ind

def run_cpu_mode(num_threads, num_rounds, chunked_skymodel_map_sets,
                 lsts, latitudes, args, beamtype,
                 main_table, shape_table, v_table, q_table, u_table, p_table,
                 woden_settings_python, array_layout_python, visi_sets_python):
    
    with ProcessPoolExecutor(max_workers=num_threads) as sky_model_executor, ProcessPoolExecutor(max_workers=num_threads) as visi_executor:
        
        # # Loop through multiple rounds of data reading and calculation
        for round_num in range(num_rounds):
            future_data_sky = [sky_model_executor.submit(read_skymodel_thread, i + round_num * num_threads,
                                        num_threads, chunked_skymodel_map_sets,
                                        lsts, latitudes,
                                        args, beamtype,
                                        main_table, shape_table,
                                        v_table, q_table, u_table, p_table)
                                for i in range(num_threads)]
            
            all_loaded_python_sources = []
            all_loaded_sources_orders = []
            for future in concurrent.futures.as_completed(future_data_sky):
                python_sources, order, round_num = future.result()
                all_loaded_python_sources.append(python_sources)
                all_loaded_sources_orders.append(order)
                
                
            future_data_visi = [visi_executor.submit(woden_multithread, i,
                                        [all_loaded_python_sources[i]], #all_loaded_python_sources
                                        [i], #all_loaded_sources_orders
                                        round_num,
                                        woden_settings_python,
                                        array_layout_python, 
                                        visi_sets_python[i, :],
                                        beamtype,
                                        args.precision)
                                for i in range(num_threads)]
            
            for future in concurrent.futures.as_completed(future_data_visi):
                visi_set_python, thread_ind = future.result()
                visi_sets_python[thread_ind, :] = visi_set_python
                    
    return visi_sets_python

# def run_cpu_mode_alt(num_threads, num_rounds, chunked_skymodel_map_sets,
#                  lsts, latitudes, args, beamtype,
#                  main_table, shape_table, v_table, q_table, u_table, p_table,
#                  woden_settings_python, array_layout_python, visi_sets_python):
    
#     num_pols = 8
#     num_bands = woden_settings_python.num_bands
#     num_visis = woden_settings_python.num_visis
#     visi_sums = multiprocessing.Array('d', [0] * num_threads * num_pols * num_visis * num_bands)
    
#     with ProcessPoolExecutor(max_workers=num_threads) as sky_model_executor, ProcessPoolExecutor(max_workers=num_threads) as visi_executor:
#         futures_sky_model_exec = []  # To hold futures for executor1
#         futures_visi_exec = []  # To hold futures for executor2

        
#         for round_num in range(num_rounds):
#             for thread_num in range(num_threads):  # Split input into chunks
#                 futures_sky_model_exec.append(sky_model_executor.submit(read_skymodel_thread,
#                                                 thread_num + round_num * num_threads,
#                                                 num_threads, chunked_skymodel_map_sets,
#                                                 lsts, latitudes,
#                                                 args, beamtype,
#                                                 main_table, shape_table,
#                                                 v_table, q_table, u_table, p_table))

#         # Process pipeline
#         while futures_sky_model_exec or futures_visi_exec:
#             # Collect completed tasks from executor1 and submit to executor2
#             ready_executor1 = [fut for fut in futures_sky_model_exec if fut.done()]
#             for future in ready_executor1:
#                 futures_sky_model_exec.remove(future)
#                 try:
#                     python_sources, thread_num, round_num = future.result()
#                     # print(f"executor1 produced: {intermediate_result}")
#                     futures_visi_exec.append(visi_executor.submit(woden_multithread, thread_num,
#                                                 [python_sources],
#                                                 [thread_num],
#                                                 round_num,
#                                                 woden_settings_python,
#                                                 array_layout_python, 
#                                                 visi_sets_python[thread_num, :],
#                                                 beamtype,
#                                                 args.precision))
#                 except Exception as e:
#                     print(f"sky_model_executor task failed: {e}")

#             # Collect completed tasks from executor2
#             ready_executor2 = [fut for fut in futures_visi_exec if fut.done()]
#             for future in ready_executor2:
#                 futures_visi_exec.remove(future)
#                 try:
#                     visi_set_python, thread_num = future.result()
#                     visi_sets_python[thread_num, :] = visi_set_python
                    
#                     with visi_sums.get_lock():  # Synchronize updates to the shared array
#                         for band_ind in range(num_bands):
#                             stripe = band_ind*num_bands*num_pols*num_visis + num_pols*num_visis*thread_num
                        
#                             visi_sums[stripe+0*num_visis:stripe+1*num_visis] += visi_set_python[band_ind].sum_visi_XX_real
#                             visi_sums[stripe+1*num_visis:stripe+2*num_visis] += visi_set_python[band_ind].sum_visi_XX_imag
#                             visi_sums[stripe+2*num_visis:stripe+3*num_visis] += visi_set_python[band_ind].sum_visi_XY_real
#                             visi_sums[stripe+3*num_visis:stripe+4*num_visis] += visi_set_python[band_ind].sum_visi_XY_imag
#                             visi_sums[stripe+4*num_visis:stripe+5*num_visis] += visi_set_python[band_ind].sum_visi_YX_real
#                             visi_sums[stripe+5*num_visis:stripe+6*num_visis] += visi_set_python[band_ind].sum_visi_YX_imag
#                             visi_sums[stripe+6*num_visis:stripe+7*num_visis] += visi_set_python[band_ind].sum_visi_YY_real
#                             visi_sums[stripe+7*num_visis:stripe+8*num_visis] += visi_set_python[band_ind].sum_visi_YY_imag
                    
#                 except Exception as e:
#                     print(f"visi_executor task failed: {e}", thread_num)

#     return visi_sets_python


def run_gpu_mode(num_threads, num_rounds, chunked_skymodel_map_sets,
                 lsts, latitudes, args, beamtype,
                 main_table, shape_table, v_table, q_table, u_table, p_table,
                 woden_settings_python, array_layout_python, visi_sets_python, 
                 run_woden, woden_struct_classes):
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as sky_model_executor:
        
        visi_executor = concurrent.futures.ThreadPoolExecutor(max_workers=1)
        gpu_calc = None  # To store the future of the calculation thread
        

    # # Loop through multiple rounds of data reading and calculation
        for round_num in range(num_rounds):
            future_data_sky = [sky_model_executor.submit(read_skymodel_thread, i + round_num * num_threads,
                                        num_threads, chunked_skymodel_map_sets,
                                        lsts, latitudes,
                                        args, beamtype,
                                        main_table, shape_table,
                                        v_table, q_table, u_table, p_table)
                                for i in range(num_threads)]
            
            all_loaded_python_sources = []
            all_loaded_sources_orders = []
            for future in concurrent.futures.as_completed(future_data_sky):
                python_sources, order, _ = future.result()
                all_loaded_python_sources.append(python_sources)
                all_loaded_sources_orders.append(order)
                
                
            
            # Wait for previous calculation to complete (if there was one)
            if gpu_calc is not None:
                visi_sets_python = gpu_calc.result()

            gpu_calc = visi_executor.submit(woden_thread,
                                        all_loaded_python_sources,
                                        all_loaded_sources_orders,
                                        round_num,
                                        run_woden, woden_settings_python,
                                        array_layout_python, 
                                        visi_sets_python,
                                        woden_struct_classes, 
                                        beamtype,
                                        args.precision)
            
        # # Wait for the final calculation to complete
        if gpu_calc is not None:
            visi_sets_python = gpu_calc.result()
            
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
    ##Grab the parser and parse some args
    parser = get_parser()
    args = parser.parse_args(argv)
    
    ##Check that the input arguments make sense
    args = check_args(args)
    
    ##Find out what git/release version we are using, and where the code lives
    gitlabel = get_code_version()

    lst_deg, gst0_deg, degpdy, ut1utc = get_uvfits_date_and_position_constants(latitude=args.latitude, longitude=args.longitude,
                        height=args.array_height, date=args.date)

    int_jd, float_jd = calc_jdcal(args.date)
    jd_date = int_jd + float_jd

    ##Check the uvfits prepend to make sure we end in .uvfits
    output_uvfits_prepend = args.output_uvfits_prepend
    if output_uvfits_prepend[-7:] == '.uvfits': output_uvfits_prepend = output_uvfits_prepend[-7:]

    ##Generate the woden_struct_classes, which are used to mimic the C
    ##structs. Must be made dynamically, as the precision can change
    woden_struct_classes = Woden_Struct_Classes(args.precision)
    
    ##Depending on what precision was selected by the user, load in the
    ##C/CUDA library and return the `run_woden` function
    run_woden = load_in_woden_library(woden_struct_classes)
    
    num_baselines = int(((args.num_antennas - 1)*args.num_antennas) / 2)

    if args.do_autos:
        num_visis = args.num_time_steps*args.num_freq_channels*(num_baselines + args.num_antennas)
    else:
        num_visis = args.num_time_steps*args.num_freq_channels*num_baselines

    ##populates a ctype equivalent of woden_settings struct to pass
    ##to the C library
    woden_settings = create_woden_settings(woden_struct_classes.Woden_Settings(),
                                            args, jd_date, lst_deg)
    
    ##fill the lst and mjds fields, precessing if necessary
    lsts, latitudes = setup_lsts_and_phase_centre(woden_settings)
    

    ##calculate the array layout
    array_layout = calc_XYZ_diffs(woden_settings, args)

    # ##read in and chunk the sky model=======================================
    print("Doing the initial reading/mapping of sky model into chunks")
    t_before = time()
    
    comp_counter = read_radec_count_components(args.cat_filename)
    t_after = time()
    print(f"Mapping took {(t_after - t_before)/60.0:.1f} mins")
    
    crop_by_component = True
    if args.sky_crop_sources:
        crop_by_component = False
        
    comp_counter = crop_below_horizon(woden_settings.lsts[0],
                                        woden_settings.latitude,
                                        comp_counter,
                                        crop_by_component=crop_by_component)
    max_chunks_per_set=args.num_threads
    
    chunked_skymodel_map_sets = create_skymodel_chunk_map(comp_counter,
                                                      args.chunking_size, woden_settings.num_baselines,
                                                      args.num_freq_channels,
                                                      args.num_time_steps,
                                                      num_threads=args.num_threads,
                                                      max_chunks_per_set=max_chunks_per_set,
                                                      max_dirs=args.max_sky_directions)
    
    woden_settings_python = convert_woden_settings_to_python(woden_settings)
    array_layout_python = convert_array_layout_to_python(array_layout, woden_settings.num_time_steps)
    
    if args.dry_run:
        ##User only wants to check whether the arguments would have worked or not
        pass
    else:
        ##Create the shapelet basis functions
        sbf = create_sbf(precision=args.precision)

        ##Create an array of visibility_sets, which get fed into run_woden
        ##and store the output visibilities
        visi_set_array = setup_visi_set_array(woden_struct_classes.Visi_Set,
                                              len(args.band_nums), num_visis,
                                              precision=args.precision)
        
        if woden_settings_python.do_gpu:
            num_visi_threads = 1
        else:
            num_visi_threads = args.num_threads
        
        ##Something to hold all visibility outputs for all threads, for all bands
        visi_sets_python = [[Visi_Set_Python() for _ in range(len(args.band_nums))] for _ in range(num_visi_threads)]
        visi_sets_python = np.array(visi_sets_python)
        
        print("DIS", visi_sets_python.shape)
        
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
            

        ###---------------------------------------------------------------------
        ### heavy lifting area - here we setup running the sky model reading
        ### and running CPU/GPU code at the same time. Means we can limit the
        ### absolute amount of RAM used, and save time doing IO at the same
        ### time as compute
        
        num_threads = args.num_threads
        
        num_rounds = len(chunked_skymodel_map_sets)
        main_table, shape_table, v_table, q_table, u_table, p_table = get_skymodel_tables(args.cat_filename)
        
        if args.cpu_mode:
            run_cpu_mode(num_threads, num_rounds, chunked_skymodel_map_sets,
                 lsts, latitudes, args, woden_settings.beamtype,
                 main_table, shape_table, v_table, q_table, u_table, p_table,
                 woden_settings_python, array_layout_python, visi_sets_python)
            # run_cpu_mode_alt(num_threads, num_rounds, chunked_skymodel_map_sets,
            #      lsts, latitudes, args, woden_settings.beamtype,
            #      main_table, shape_table, v_table, q_table, u_table, p_table,
            #      woden_settings_python, array_layout_python, visi_sets_python)
        else:
            run_gpu_mode(num_threads, num_rounds, chunked_skymodel_map_sets,
                 lsts, latitudes, args, woden_settings.beamtype,
                 main_table, shape_table, v_table, q_table, u_table, p_table,
                 woden_settings_python, array_layout_python, visi_sets_python[0], 
                 run_woden, woden_struct_classes)
            
        
        
        ### we've now calculated all the visibilities
        ###---------------------------------------------------------------------
        
        visi_sets_python_combined = visi_sets_python[0, :]
        
        for thread_ind in range(1, num_visi_threads):
            # visi_set_thread = visi_sets_python[thread_ind, :]
            for band_ind in range(len(args.band_nums)):
                # visi_set = visi_set_thread[thread_ind, band_ind]
                # visi_set_combined = visi_sets_python_combined[0, band_ind]
                
                visi_sets_python_combined[band_ind].sum_visi_XX_real += visi_sets_python[thread_ind, band_ind].sum_visi_XX_real
                visi_sets_python_combined[band_ind].sum_visi_XX_imag += visi_sets_python[thread_ind, band_ind].sum_visi_XX_imag
                visi_sets_python_combined[band_ind].sum_visi_XY_real += visi_sets_python[thread_ind, band_ind].sum_visi_XY_real
                visi_sets_python_combined[band_ind].sum_visi_XY_imag += visi_sets_python[thread_ind, band_ind].sum_visi_XY_imag
                visi_sets_python_combined[band_ind].sum_visi_YX_real += visi_sets_python[thread_ind, band_ind].sum_visi_YX_real
                visi_sets_python_combined[band_ind].sum_visi_YX_imag += visi_sets_python[thread_ind, band_ind].sum_visi_YX_imag
                visi_sets_python_combined[band_ind].sum_visi_YY_real += visi_sets_python[thread_ind, band_ind].sum_visi_YY_real
                visi_sets_python_combined[band_ind].sum_visi_YY_imag += visi_sets_python[thread_ind, band_ind].sum_visi_YY_imag
                
                
        # print("CMON", visi_sets_python_combined[0].us_metres)
        
        ##I think we want to X,Y,Z to be in the current frame for writing
        ##out to the uvfits, so calculate again
        X,Y,Z = enh2xyz(args.east, args.north, args.height, args.latitude*D2R)
        ##X,Y,Z are stored in a 2D array in units of seconds in the uvfits file
        XYZ_array = np.empty((args.num_antennas,3))
        XYZ_array[:,0] = X
        XYZ_array[:,1] = Y
        XYZ_array[:,2] = Z

        ##If we were to grab them out of the code used in the simulation,
        ##need to do this
        # expec_num = woden_settings.num_time_steps*args.num_antennas
        # XYZ_array[:,0] = np.ctypeslib.as_array(array_layout.ant_X, shape=(expec_num, ))[:args.num_antennas]
        # XYZ_array[:,1] = np.ctypeslib.as_array(array_layout.ant_Y, shape=(expec_num, ))[:args.num_antennas]
        # XYZ_array[:,2] = np.ctypeslib.as_array(array_layout.ant_Z, shape=(expec_num, ))[:args.num_antennas]

        ##Get the central frequency channels, used in the uvfits header
        central_freq_chan = int(np.floor(args.num_freq_channels / 2.0))
        ##Useful number
        num_baselines = int(((args.num_antennas - 1)*args.num_antennas) / 2)

        ##Loop over coarse frequency band and convert visibilities output
        ##from woden.c into uvfits files
        for band_ind, band in enumerate(args.band_nums):

            output_uvfits_name = output_uvfits_prepend + '_band%02d.uvfits' %band

            band_low_freq = args.lowest_channel_freq + (band - 1)*args.coarse_band_width
            central_freq_chan_value = band_low_freq + central_freq_chan*args.freq_res

            uus,vvs,wws,v_container = load_visibility_set(visibility_set=visi_sets_python_combined[band_ind],
                                                num_baselines=num_baselines,
                                                num_freq_channels=args.num_freq_channels,
                                                num_time_steps=args.num_time_steps,
                                                precision=args.precision,
                                                do_autos=args.do_autos,
                                                num_ants=args.num_antennas)
            
            # print("THE FINAL FORM", uus)

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
                          ch_width=args.freq_res, int_jd=int_jd,
                          hdu_ant=hdu_ant, gitlabel=gitlabel,
                          longitude=args.longitude, latitude=args.latitude,
                          array_height=args.array_height,
                          telescope_name=args.telescope_name,
                          IAU_order=args.IAU_order,
                          comment=args.command)

if __name__ == "__main__":
    main()