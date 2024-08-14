import ctypes 
import importlib_resources
import wodenpy
import numpy as np
import sys
import os

from wodenpy.use_libwoden.create_woden_struct_classes import Woden_Struct_Classes
from wodenpy.use_libwoden.array_layout_struct import Array_Layout

VELC = 299792458.0

def load_in_woden_library(woden_struct_classes : Woden_Struct_Classes, precision='double'):
    """Load in the WODEN C and CUDA code via a dynamic library, with the
    required `precision` (either load `libwoden_float.so` or `libwoden_double.so`)

    Parameters
    ----------
    woden_struct_classes : Woden_Struct_Classes
        This holds all the various ctype structure classes that are equivalent
        to the C/CUDA structs. Should have been initialised with the correct
        precision ("float" or "double").
    precision : str, optional
        Choose to use "float" precision or "double". By default 'double'

    Returns
    -------
    run_woden : _NamedFuncPointer
        The C wrapper function `run_woden`, which runs the C/CUDA code.
        This function takes the following args when `precision="double"`
         - ctypes.POINTER(Woden_Settings_Double)
         - ctypes.POINTER(Visi_Set_Double)
         - ctypes.POINTER(Source_Catalogue_Double)
         - ctypes.POINTER(Array_Layout)
         - ctypes.POINTER(ctypes.c_double)
        This function takes the following args when `precision="float"`
         - ctypes.POINTER(Woden_Settings_Float)
         - ctypes.POINTER(Visi_Set_Float)
         - ctypes.POINTER(Source_Catalogue_Float)
         - ctypes.POINTER(Array_Layout)
         - ctypes.POINTER(ctypes.c_float)

        
    """
    
    woden_lib = importlib_resources.files(wodenpy).joinpath(f"libwoden_{woden_struct_classes.precision}.so")
    
    print("LOADING IN", woden_lib)

    ## Read in the C library
    libwoden = ctypes.cdll.LoadLibrary(woden_lib)

    # ##Define the input and return types for the `test_RTS_calculate_MWA_analytic_beam` function
    run_woden = libwoden.run_woden

    run_woden.restype = ctypes.c_int
    
    Source_Catalogue = woden_struct_classes.Source_Catalogue
    
    if precision == 'float':
        run_woden.argtypes = [ctypes.POINTER(woden_struct_classes.Woden_Settings),
                              ctypes.POINTER(woden_struct_classes.Visi_Set),
                              ctypes.POINTER(woden_struct_classes.Source_Catalogue),
                              ctypes.POINTER(Array_Layout),
                              ctypes.POINTER(ctypes.c_float)]
    else:
        run_woden.argtypes = [ctypes.POINTER(woden_struct_classes.Woden_Settings),
                              ctypes.POINTER(woden_struct_classes.Visi_Set),
                              ctypes.POINTER(woden_struct_classes.Source_Catalogue),
                              ctypes.POINTER(Array_Layout),
                              ctypes.POINTER(ctypes.c_double)]
        
    return run_woden