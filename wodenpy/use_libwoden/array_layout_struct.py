import numpy as np
from ctypes import POINTER, c_double, c_int, Structure
import sys
import os
from typing import Union
import argparse
from wodenpy.use_libwoden.create_woden_struct_classes import Woden_Struct_Classes

VELC = 299792458.0

class Array_Layout_Ctypes(Structure):
    """
    A ctypes structure representing the layout of an array of antennas.

    Attributes:
        ant_X (POINTER(c_double)): Pointer to an array of antenna X positions.
        ant_Y (POINTER(c_double)): Pointer to an array of antenna Y positions.
        ant_Z (POINTER(c_double)): Pointer to an array of antenna Z positions.
        X_diff_metres (POINTER(c_double)): Pointer to an array of X position differences in metres.
        Y_diff_metres (POINTER(c_double)): Pointer to an array of Y position differences in metres.
        Z_diff_metres (POINTER(c_double)): Pointer to an array of Z position differences in metres.
        ant_east (POINTER(c_double)): Pointer to an array of antenna east positions.
        ant_north (POINTER(c_double)): Pointer to an array of antenna north positions.
        ant_height (POINTER(c_double)): Pointer to an array of antenna height positions.
        latitude (c_double): The latitude of the array.
        num_baselines (c_int): The number of baselines in the array.
        num_tiles (c_int): The number of tiles in the array.
        lst_base (c_double): The local sidereal time of the array.
    """
    _fields_ = [("ant_X", POINTER(c_double)),
                ("ant_Y", POINTER(c_double)),
                ("ant_Z", POINTER(c_double)),
                ("X_diff_metres", POINTER(c_double)),
                ("Y_diff_metres", POINTER(c_double)),
                ("Z_diff_metres", POINTER(c_double)),
                ("ant_east", POINTER(c_double)),
                ("ant_north", POINTER(c_double)),
                ("ant_height", POINTER(c_double)),
                ("latitude", c_double),
                ("num_baselines", c_int),
                ("num_tiles", c_int),
                ("lst_base", c_double)]
    
##This call is so we can use it as a type annotation
woden_struct_classes = Woden_Struct_Classes()
Woden_Settings = woden_struct_classes.Woden_Settings
    
def setup_array_layout_ctypes(woden_settings : Woden_Settings, args : argparse.Namespace) -> Array_Layout_Ctypes: #type: ignore
    """Given the populated Woden_Settings struct and the command line arguments, set up the `array_layout` struct, and fill it with the correct values, as well as doing the equivalent of a "malloc" for the arrays:
     - array_layout.ant_X
     - array_layout.ant_Y
     - array_layout.ant_Z
     - array_layout.X_diff_metres
     - array_layout.Y_diff_metres
     - array_layout.Z_diff_metres

    Parameters
    ----------
    woden_settings : Woden_Settings
        Populated Woden_Settings struct.
    args : argparse.Namespace
        Input arguments from the command line that have been checked using
        `wodenpy.use_libwoden.check_args.check_args`.

    Returns
    -------
    array_layout : Array_Layout_Ctypes
        Initialised Array_Layout_Ctypes struct.
    """

    array_layout = Array_Layout_Ctypes()

    array_layout.num_tiles = args.num_antennas
    woden_settings.num_ants = args.num_antennas

    array_layout.latitude = woden_settings.latitude
    array_layout.num_baselines = int((array_layout.num_tiles*(array_layout.num_tiles-1)) / 2)
    woden_settings.num_baselines = array_layout.num_baselines

    num_cross = woden_settings.num_baselines * woden_settings.num_time_steps * woden_settings.num_freqs

    woden_settings.num_cross = num_cross

    ##If user asked for auto-correlations, set this number to total autos to work out
    ##Otherwise stick it to zero
    if woden_settings.do_autos:
        num_autos = woden_settings.num_ants * woden_settings.num_time_steps * woden_settings.num_freqs
    else:
        num_autos = 0
    
    woden_settings.num_autos = num_autos
    woden_settings.num_visis = num_cross + num_autos
    
    num_timesants_array = (woden_settings.num_time_steps*args.num_antennas)*c_double
    array_layout.ant_X = num_timesants_array()
    array_layout.ant_Y = num_timesants_array()
    array_layout.ant_Z = num_timesants_array()

    diffs_array = (woden_settings.num_time_steps*array_layout.num_baselines)*c_double
    array_layout.X_diff_metres = diffs_array()
    array_layout.Y_diff_metres = diffs_array()
    array_layout.Z_diff_metres = diffs_array()
    
    return array_layout