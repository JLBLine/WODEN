from sys import path
import os
import unittest
import numpy as np
import ctypes
import numpy.testing as npt

##Code we are testing
from wodenpy.use_libwoden.shapelets import create_sbf

##Indexes to test in sbf
test_indexes = np.array([5000, 15011, 25022, 35033, 45044, 55055,
                        65066, 75077, 85088, 95099, 105110, 115121,
                        125132, 135143, 145154, 155165, 165176,
                        175187, 185198, 195209], dtype=int)

##Matching expected values
expected = np.array([1.00000000, 0.14071601, -0.63765672, -0.46694581,
                       0.22280003, 0.58389053, 0.31120116, -0.2327714,
                       -0.52281663, -0.35279668, 0.081828458, 0.42677188,
                       0.45464097, 0.18577076, -0.18055274, -0.42282655,
                       -0.42453389, -0.21182068, 0.09026158, 0.33683096])

##Vehicle for running tests
class Test(unittest.TestCase):
    """Check the shapelet basis is created correctly by sampling a few values"""
    
    def test_create_sbf(self):
        
        sbf = create_sbf()
        
        ##The sbf is a ctypes object, so iterate over values as can't
        ##directly read whole array without fiddling
        
        for index, value in zip(test_indexes, expected):
        
            npt.assert_allclose(sbf[index], value, rtol=1e-7)
        
        
        
##Run the test
if __name__ == '__main__':
   unittest.main()