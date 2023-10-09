from sys import path
import os
import unittest
import numpy as np

##Do some disgusting path finding exercise, there must be a better way
##to do this
fileloc = os.path.realpath(__file__)
path.append('{:s}/../../../wodenpy/array_layout'.format(('/').join(fileloc.split('/')[:-1])))

import create_array_layout

D2R = np.pi / 180.0

##Vehicle for running tests
class Test(unittest.TestCase):
    def test_enh2xyz(self):
        """Tests the `create_array_layout.enh2xyz` function, which should calculate the
        local X,Y,Z coord using the local east, north, height. Test using
        the cases where latitude is 0 and -30 deg, which have specific
        outcomes"""

        ##Dummy array of boringness
        east = np.arange(0,100,10)
        north = np.arange(0,100,10)
        height = np.arange(0,100,10)

        ##Test a latitude of 0.0 - X,Y,Z and e,n,h just map to one another
        latitude = 0.0
        X,Y,Z = create_array_layout.enh2xyz(east, north, height, latitude=latitude)

        self.assertTrue(np.array_equal(height, X))
        self.assertTrue(np.array_equal(east, Y))
        self.assertTrue(np.array_equal(north, Z))

        ##Test at latitude of -30, mapping to X,Y,Z is somewhat simple still
        latitude = -30*D2R
        X,Y,Z = create_array_layout.enh2xyz(east, north, height, latitude=latitude)

        self.assertTrue(np.allclose(0.5*north + (np.sqrt(3)/2)*height, X, atol=1e-15))
        self.assertTrue(np.array_equal(east, Y))
        self.assertTrue(np.allclose((np.sqrt(3)/2)*north + -0.5*height, Z, atol=1e-15))

##Run the test
if __name__ == '__main__':
   unittest.main()
