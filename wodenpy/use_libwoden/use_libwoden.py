"""Functions to load in the WODEN C/C++/GPU code via a dynamic library,
with the required `precision` (either load `libwoden_float.so` or `libwoden_double.so`)."""

import ctypes 
import importlib_resources
import wodenpy
import numpy as np
import sys
import os
from ctypes import CFUNCTYPE
from multiprocessing import Queue, Process
from wodenpy.use_libwoden.create_woden_struct_classes import Woden_Struct_Classes
from wodenpy.use_libwoden.array_layout_struct import Array_Layout_Ctypes


VELC = 299792458.0

def load_in_run_woden(woden_lib : ctypes.CDLL,
                      woden_struct_classes : Woden_Struct_Classes):
    """Load in and define the C wrapper function `run_woden`, which runs the C/C++/GPU code.
    `woden_lib` is the ctypes object which has loaded in the WODEN library
    via `ctypes.cdll.LoadLibrary(woden_lib_path)`.
    
    Here `woden_lib_path` is the path to the WODEN library, which is either
    `libwoden_float.so` or `libwoden_double.so`, depending on the precision. 

    Parameters
    ----------
    woden_lib : ctypes.CDLL
        The ctypes object which has loaded in the WODEN library
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
    
    #Select the run_woden function and define the return type
    run_woden = woden_lib.run_woden
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
                            ctypes.POINTER(Array_Layout_Ctypes),
                            sbf_pointer]
    
    return run_woden

def worker_check_for_everybeam(woden_lib_path: str, q : Queue):
    """
    Checks if libwoden*.so has been compiled against EveryBeam (via a flag
    -DHAVE_EVERYBEAM which was set via CMake during compilation). Puts True
    into the queue `q`, if it has, False otherwise.
    
    This function is run in a separate process to avoid clashing with 
    `python-casacore`. `python-casacore` and `libwoden*.so`
    can be linked to different versions of `casacore`. When
    loaded, both create dependency paths and global variables in a persisting
    `c++` state. This causes much explosions and segfaults. Wrapping in a 
    separate process isolates the state of the `c++` code. This is a hot-fix.
    
    .. todo::
        Replace all calls to `python-casacore` with calls to `casacore`
        directly, via libuse_everybeam.so. It's also possible that using 
        pybind11 to call c++ instead of `ctypes` would be better?
        
    Parameters
    ----------
    woden_lib_path : str
        The file path to the WODEN library (either `/some/path/libwoden_float.so`
        or `/some/path/libwoden_double.so`).
    q : Queue
        A multiprocessing queue to store the result of the check.
    
    """
    
    woden_lib = ctypes.cdll.LoadLibrary(woden_lib_path)
    
    check_for_everybeam_compilation = woden_lib.check_for_everybeam_compilation
    
    check_for_everybeam_compilation.restype = ctypes.c_bool
    
    q.put(check_for_everybeam_compilation())


def check_for_everybeam(woden_lib_path: str) -> bool:
    """
    Checks if libwoden*.so has been compiled against EveryBeam (via a flag
    -DHAVE_EVERYBEAM which was set via CMake during compilation).
    
    Returns True if it has, False otherwise.
    
    This function calls :func:~`worker_check_for_everybeam` in a separate
    process to avoid clashing with `python-casacore`. See 
    :func:~`worker_check_for_everybeam` docs for more details.
    
    Parameters
    ----------
    woden_lib_path : str
        The file path to the WODEN library (either `/some/path/libwoden_float.so`
        or `/some/path/libwoden_double.so`).
    
    Returns
    -------
    bool
        True if the 'EveryBeam' feature is compiled in the WODEN library, False otherwise.
    """
    
    q = Queue()
    p = Process(target=worker_check_for_everybeam, args=(woden_lib_path, q))
    p.start()
    p.join()
    
    return q.get()