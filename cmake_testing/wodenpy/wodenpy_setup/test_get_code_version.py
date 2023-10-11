from sys import path
import os
import unittest
import numpy as np

code_dir = os.path.realpath(__file__)
code_dir = ('/').join(code_dir.split('/')[:-1])

# ##Code we are testing
from wodenpy.wodenpy_setup import run_setup



##Vehicle for running tests
class Test(unittest.TestCase):
    """Test whether the args collected by the argument parser are read in
    correctly, and that sanity checks on certain combinations work such that
    we don't feed WODEN arguments that won't work"""

    def test_get_code_version(self):
        """Call `run_setup.code_version` and just check it runs without error.
        Hard to automate a test for this one"""
        version = run_setup.get_code_version()
        print(version)

    

##Run the test
if __name__ == '__main__':
    unittest.main()
    