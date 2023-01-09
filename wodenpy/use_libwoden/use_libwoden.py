import ctypes 
import importlib_resources
import wodenpy
import numpy as np
from wodenpy.use_libwoden.visibility_set import *
from wodenpy.use_libwoden.woden_settings import *

VELC = 299792458.0

def load_in_woden_library(precision='double'):
    """Load in the WODEN C and CUDA code"""

    woden_lib = importlib_resources.files(wodenpy).joinpath(f"libwoden_{precision}.so")

    ## Read in the C library
    libwoden = ctypes.cdll.LoadLibrary(woden_lib)

    # ##Define the input and return types for the `test_RTS_calculate_MWA_analytic_beam` function
    run_woden = libwoden.run_woden

    run_woden.restype = ctypes.c_int
    
    if precision == 'float':
        run_woden.argtypes = [ctypes.POINTER(Woden_Settings_Float),
                              ctypes.POINTER(Visi_Set_Float)]
    else:
        run_woden.argtypes = [ctypes.POINTER(Woden_Settings_Double),
                              ctypes.POINTER(Visi_Set_Double)]

    return run_woden