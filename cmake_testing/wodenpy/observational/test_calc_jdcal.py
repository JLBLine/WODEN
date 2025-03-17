from sys import path
import os
import unittest
from wodenpy.observational import calc_obs
import numpy.testing as npt

##Vehicle for running tests
class Test(unittest.TestCase):
    def test_cal_jdcal(self):
        """Tests the `rw.calc_jdcal` function, which should calculate the
        Julian Date and split into a day and fractional day. Just test by
        running a couple of dates"""

        ##Run the code with an example date
        date1 = "2019-06-12T13:04:12"
        jd_day1, jd_fraction1 = calc_obs.calc_jdcal(date1)

        ##Test it gives correct outputs
        npt.assert_equal(2458646.5, jd_day1)
        npt.assert_allclose(0.5445833331905305, jd_fraction1)
        
        ##old values from just doing a floor on the jd date
        sum1 = 2458647.0 + 0.04458333319
        npt.assert_allclose(jd_day1 + jd_fraction1, sum1, atol=1e-10)

        ##Do the same again for another date
        date2 = "2013-01-12T09:45:09"
        jd_day2, jd_fraction2 = calc_obs.calc_jdcal(date2)

        npt.assert_equal(2456304.5, jd_day2)
        npt.assert_allclose(0.4063541665673256, jd_fraction2)
        
        ##old values from just doing a floor on the jd date
        sum2 = 2456304.0 + 0.90635416656
        npt.assert_allclose(jd_day2 + jd_fraction2, sum2, atol=1e-10)

##Run the test
if __name__ == '__main__':
   unittest.main()
