from sys import path
import os
import unittest
import numpy as np
import ctypes
from wodenpy.use_libwoden.use_libwoden import check_for_everybeam


test_dir = os.environ['CMAKE_CURRENT_SOURCE_DIR'] + "/../../../build/cmake_testing/wodenpy/use_libwoden/"

##Vehicle for running tests
class Test(unittest.TestCase):
    """Check the shapelet basis is created correctly by sampling a few values"""
    
    def test_check_for_everybeam(self):
        
        with_flag_lib = test_dir + "libcheck_comp_has_flag_with_flag.so"
        have_everybeam = check_for_everybeam(with_flag_lib)
        self.assertTrue(have_everybeam)
        
        with_flag_lib = test_dir + "libcheck_comp_has_flag_without_flag.so"
        have_everybeam = check_for_everybeam(with_flag_lib)
        self.assertFalse(have_everybeam)
        
        
        
        
##Run the test
if __name__ == '__main__':
   unittest.main()