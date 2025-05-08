from sys import path
import os
import unittest
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
import concurrent
import sys

# from astropy.io import fits
# from unittest import mock
# import argparse
# from astropy.time import Time, TimeDelta
# import numpy as np
# from astropy.coordinates import EarthLocation
# from astropy import units as u
# from astropy.constants import c as speed_of_light

# from wodenpy.array_layout.precession import RTS_Precess_LST_Lat_to_J2000
# from astropy.time import Time, TimeDelta
# import numpy as np
# from astropy.coordinates import EarthLocation
# from astropy import units as u
# import scipy.optimize as opt
# from glob import glob

code_dir = os.path.realpath(__file__)
code_dir = ('/').join(code_dir.split('/')[:-1])

##Append the location of add_woden_uvfits to the sys.path to import from it
path.append('{:s}/../../../scripts'.format(code_dir))


##annoying path hack to find where the C library is. The C library is
##build during ctest runs
test_dir = code_dir + "/../../../build/cmake_testing/scripts/run_woden/"
testing_lib = f"{test_dir}/libdo_simple_functions.so"

##Code we are testing
# from wodenpy.uvfits import wodenpy_uvfits
import run_woden as rw
import numpy.testing as npt
from wodenpy.wodenpy_setup.woden_logger import simple_logger
from ctypes import c_int, c_float, c_double, c_char_p, c_void_p
import ctypes

# from multiprocessing import Manager, set_start_method

# from subprocess import call

def proc_one_works_python(thread_ind):
    return thread_ind, thread_ind*2

def proc_two_works_python(thread_ind, input):
    return thread_ind, input**2

def proc_one_fails_python(thread_ind):
    return thread_ind, 1 / 0

def proc_two_fails_python(thread_ind, input):
    return thread_ind, 1 / 0

def proc_one_works_c(thread_ind):
    
    c_lib = ctypes.cdll.LoadLibrary(testing_lib)
    proc_one_works_c = c_lib.proc_one_works_c
    proc_one_works_c.argtypes = [ctypes.c_int]
    proc_one_works_c.restype = ctypes.c_int
    
    output = proc_one_works_c(thread_ind)
    
    return thread_ind, output

def proc_two_works_c(thread_ind, input):
    
    c_lib = ctypes.cdll.LoadLibrary(testing_lib)
    proc_two_works_c = c_lib.proc_two_works_c
    proc_two_works_c.argtypes = [ctypes.c_int]
    proc_two_works_c.restype = ctypes.c_int
    
    output = proc_two_works_c(input)
    
    return thread_ind, output

def proc_one_fails_c(thread_ind):
    
    c_lib = ctypes.cdll.LoadLibrary(testing_lib)
    proc_one_fails_c = c_lib.proc_one_fails_c
    proc_one_fails_c.argtypes = [ctypes.c_int]
    proc_one_fails_c.restype = ctypes.c_int
    
    # print("WE HERE")
    
    # with open("worker_stderr.log", "w") as f:
    #     os.dup2(f.fileno(), sys.stderr.fileno())
        
    output = proc_one_fails_c(thread_ind)
    
    # try:
    #     output = proc_one_fails_c(thread_ind)
    # except Exception as e:
    #     msg = f"Error in proc_one_fails_c: {e}"
    #     sys.exit(msg)
    
    return thread_ind, output

def proc_two_fails_c(thread_ind, input):
    
    c_lib = ctypes.cdll.LoadLibrary(testing_lib)
    proc_two_fails_c = c_lib.proc_two_fails_c
    proc_two_fails_c.argtypes = [ctypes.c_int]
    proc_two_fails_c.restype = ctypes.c_int
    
    output = proc_two_fails_c(input)
    
    return thread_ind, output


##Vehicle for running tests
class Test(unittest.TestCase):
    
    def run_multi_process(self, num_proc_one, num_proc_two, num_threads,
                                proc_one, proc_two):
        
        logger = simple_logger()
        logger.info("Starting multi-process test")

        with ProcessPoolExecutor(max_workers=num_proc_one) as executor_one, ProcessPoolExecutor(max_workers=num_proc_two) as executor_two:
            
                future_proc_one = [executor_one.submit(proc_one, i) for i in range(num_threads)]
                    
                proc_one_outputs = [0]*num_threads
                all_loaded_sources_orders = [0]*num_threads
                for future in future_proc_one:
                    one_index, proc_one_output = rw.get_future_result(future, logger, "proc_one")
                    
                    proc_one_outputs[one_index] = proc_one_output
                    
                logger.info("Finished proc_one")
                    
                future_proc_two = [executor_one.submit(proc_two, i, proc_one_outputs[i]) for i in range(num_threads)]
                
                proc_two_outputs = [0]*num_threads
                
                for future in concurrent.futures.as_completed(future_proc_two):
                    two_index, proc_two_output = rw.get_future_result(future, logger, "proc_two")
                    
                    proc_two_outputs[two_index] = proc_two_output
                    
                logger.info("Finished proc_two")
            
        logger.info("Finished multi-process test")
            
        for handler in logger.handlers:
            handler.close()
            logger.removeHandler(handler)
            
        return proc_two_outputs
    
    
    
    def test_recovers_outputs_correctly_python(self):
        
        num_proc_one = 4
        num_proc_two = 2
        num_threads = 4
        
        proc_two_outputs = self.run_multi_process(num_proc_one, num_proc_two, num_threads,
                                                  proc_one_works_python, proc_two_works_python)
            
        expec_output = (np.arange(num_threads)*2)**2
        npt.assert_array_equal(proc_two_outputs, expec_output)
        
    def test_recovers_outputs_correctly_c(self):
        
        num_proc_one = 4
        num_proc_two = 2
        num_threads = 4
        
        proc_two_outputs = self.run_multi_process(num_proc_one, num_proc_two, num_threads,
                                                  proc_one_works_c, proc_two_works_c)
            
        expec_output = (np.arange(num_threads)*2)**2
        npt.assert_array_equal(proc_two_outputs, expec_output)
        
    def test_bad_proc_one_python(self):
        
        num_proc_one = 4
        num_proc_two = 2
        num_threads = 4
        
        with self.assertRaises(SystemExit) as cm:
            proc_two_outputs = self.run_multi_process(num_proc_one, num_proc_two, num_threads,
                                                  proc_one_fails_python, proc_two_works_python)
            
    def test_bad_proc_two_python(self):
        
        num_proc_one = 4
        num_proc_two = 2
        num_threads = 4
        
        with self.assertRaises(SystemExit) as cm:
            proc_two_outputs = self.run_multi_process(num_proc_one, num_proc_two, num_threads,
                                                  proc_one_works_python, proc_two_fails_python)
            
    def test_bad_proc_one_c(self):
        
        num_proc_one = 4
        num_proc_two = 2
        num_threads = 4
        
        # proc_one_fails_c(1)
        
        with self.assertRaises(SystemExit) as cm:
            proc_two_outputs = self.run_multi_process(num_proc_one, num_proc_two, num_threads,
                                                  proc_one_fails_c, proc_two_works_c)
            
    def test_bad_proc_two_c(self):
        
        num_proc_one = 4
        num_proc_two = 2
        num_threads = 4
        
        with self.assertRaises(SystemExit) as cm:
            proc_two_outputs = self.run_multi_process(num_proc_one, num_proc_two, num_threads,
                                                      proc_one_works_c, proc_two_fails_c)
        
        
##Run the test
if __name__ == '__main__':
    unittest.main()
