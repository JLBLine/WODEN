import ctypes 
import importlib_resources
import wodenpy
import numpy as np
import sys
import os

from wodenpy.use_libwoden.create_woden_struct_classes import Woden_Struct_Classes
from wodenpy.use_libwoden.array_layout_struct import Array_Layout

VELC = 299792458.0

def load_in_woden_library(woden_struct_classes : Woden_Struct_Classes):
    """Load in the WODEN C and CUDA code via a dynamic library, with the
    required `precision` (either load `libwoden_float.so` or `libwoden_double.so`)

    Parameters
    ----------
    woden_struct_classes : Woden_Struct_Classes
        This holds all the various ctype structure classes that are equivalent
        to the C/CUDA structs. Should have been initialised with the correct
        precision ("float" or "double").

    Returns
    -------
    run_woden : _NamedFuncPointer
        The C wrapper function `run_woden`, which runs the C/CUDA code.
        This function takes the following args, where `sbf_pointer` is either
        ctypes.POINTER(ctypes.c_float) or ctypes.POINTER(ctypes.c_double),
        depending on the `woden_struct_classes.precision`:
         - ctypes.POINTER(woden_struct_classes.Woden_Settings)
         - ctypes.POINTER(woden_struct_classes.Visi_Set)
         - ctypes.POINTER(woden_struct_classes.Source_Catalogue)
         - ctypes.POINTER(Array_Layout)
         - sbf_pointer
    """
    
    woden_lib = importlib_resources.files(wodenpy).joinpath(f"libwoden_{woden_struct_classes.precision}.so")
    
    print("LOADING IN", woden_lib)

    ## Read in the C library
    libwoden = ctypes.cdll.LoadLibrary(woden_lib)

    #Select the run_woden function and define the return type
    run_woden = libwoden.run_woden
    run_woden.restype = ctypes.c_int
    
    ##now define the argument types; we have defined the classes needed
    ##in woden_struct_classes. Final argument is the `sbf` array, which depends
    ##on the precision required
    if woden_struct_classes.precision == 'float':
        sbf_pointer = ctypes.POINTER(ctypes.c_float)
    else:
        sbf_pointer = ctypes.POINTER(ctypes.c_double)
        
        
    run_woden.argtypes = [ctypes.POINTER(woden_struct_classes.Woden_Settings),
                            ctypes.POINTER(woden_struct_classes.Visi_Set),
                            ctypes.POINTER(woden_struct_classes.Source_Catalogue),
                            ctypes.POINTER(Array_Layout),
                            sbf_pointer]
    
    return run_woden