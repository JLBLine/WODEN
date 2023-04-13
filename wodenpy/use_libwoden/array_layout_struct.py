import numpy as np
from ctypes import POINTER, c_double, c_int, Structure
import sys
import os
from typing import Union
import argparse

##If we are performing a ctest, this check means we use the code we are
##testing and NOT what has been pip or conda installed
try:
    testdir = os.environ['CMAKE_CURRENT_SOURCE_DIR']
    sys.path.append('{:s}/../../../wodenpy/use_libwoden'.format(testdir))
    from woden_settings import Woden_Settings_Float, Woden_Settings_Double
    
except KeyError:
    from wodenpy.use_libwoden.woden_settings import Woden_Settings_Float, Woden_Settings_Double

VELC = 299792458.0

class Array_Layout(Structure):
    """A class structured equivalently to a `array_layout_t` struct, used by 
    the C and CUDA code in libwoden_double.so or libwoden_float.so
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
    
def setup_array_layout(woden_settings : Union[Woden_Settings_Float, Woden_Settings_Double], args : argparse.Namespace) -> Structure:

    array_layout = Array_Layout()

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
    
    # num_ants_array = args.num_antennas*c_double
    # array_layout.ant_east = num_ants_array()
    # array_layout.ant_north = num_ants_array()
    # array_layout.ant_heigh = num_ants_array()

    # ##TODO copy across the contents of args.east, args.north, args.height

    num_timesants_array = (woden_settings.num_time_steps*args.num_antennas)*c_double
    array_layout.ant_X = num_timesants_array()
    array_layout.ant_Y = num_timesants_array()
    array_layout.ant_Z = num_timesants_array()


    diffs_array = (woden_settings.num_time_steps*array_layout.num_baselines)*c_double
    array_layout.X_diff_metres = diffs_array()
    array_layout.Y_diff_metres = diffs_array()
    array_layout.Z_diff_metres = diffs_array()
    
    return array_layout
