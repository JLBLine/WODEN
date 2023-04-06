import ctypes 
import importlib_resources
import wodenpy
import numpy as np
import sys
import os

##If we are performing a ctest, this check means we use the code we are
##testing and NOT what has been pip or conda installed
try:
    testdir = os.environ['CMAKE_CURRENT_SOURCE_DIR']
    sys.path.append('{:s}/../../../wodenpy/use_libwoden'.format(testdir))
    
    # from use_libwoden.visibility_set import *
    # from use_libwoden.woden_settings import *
    from visibility_set import *
    from woden_settings import *
    from skymodel_structs import * 
    
except KeyError:
    from wodenpy.use_libwoden.visibility_set import *
    from wodenpy.use_libwoden.woden_settings import *
    from wodenpy.use_libwoden.skymodel_structs import *

VELC = 299792458.0

def load_in_woden_library(precision='double'):
    """Load in the WODEN C and CUDA code"""

    woden_lib = importlib_resources.files(wodenpy).joinpath(f"libwoden_{precision}.so")
    
    print("LOADING IN", woden_lib)

    ## Read in the C library
    libwoden = ctypes.cdll.LoadLibrary(woden_lib)

    # ##Define the input and return types for the `test_RTS_calculate_MWA_analytic_beam` function
    run_woden = libwoden.run_woden

    run_woden.restype = ctypes.c_int
    
    if precision == 'float':
        run_woden.argtypes = [ctypes.POINTER(Woden_Settings_Float),
                              ctypes.POINTER(Visi_Set_Float),
                              ctypes.POINTER(Source_Catalogue_Float)]
    else:
        run_woden.argtypes = [ctypes.POINTER(Woden_Settings_Double),
                              ctypes.POINTER(Visi_Set_Double),
                              ctypes.POINTER(Source_Catalogue_Double)]
        
    return run_woden