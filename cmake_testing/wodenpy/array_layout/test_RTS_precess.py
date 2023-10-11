from sys import path
import os
import unittest
import numpy as np

##Do some disgusting path finding exercise, there must be a better way
##to do this

from wodenpy.array_layout import precession

D2R = np.pi / 180.0

##Vehicle for running tests
class Test(unittest.TestCase):
    def test_RTS_Precess_LST_Lat_to_J2000(self):
        
        lst_current = 0.0077340
        latitude_current = -0.4660608377603597

        jd_date = 2457278.201145833
        sec_to_d = 1 / (24.0*60.0*60.0)
        
        mjd = 57277.7010995000600815

        expec_lst = 0.004245887738
        expec_lat = -0.467587211832
        
        lst_J2000, latitude_J2000 = precession.RTS_Precess_LST_Lat_to_J2000(lst_current,
                                                                            latitude_current,
                                                                            mjd)
        
        self.assertAlmostEqual(lst_J2000, expec_lst, delta=1e-12)
        self.assertAlmostEqual(latitude_J2000, expec_lat, delta=1e-12)

##Run the test
if __name__ == '__main__':
   unittest.main()
