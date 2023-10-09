from sys import path
import os
import unittest

##Do some disgusting path finding exercise, there must be a better way
##to do this
fileloc = os.path.realpath(__file__)
path.append('{:s}/../../../wodenpy/observational'.format(('/').join(fileloc.split('/')[:-1])))

import calc_obs

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
        self.assertEqual(2458647.0, jd_day1)
        self.assertAlmostEqual(0.04458333319, jd_fraction1)

        ##Do the same again for another date
        date2 = "2013-01-12T09:45:09"
        jd_day2, jd_fraction2 = calc_obs.calc_jdcal(date2)

        self.assertEqual(2456304.0, jd_day2)
        self.assertAlmostEqual(0.90635416656, jd_fraction2)

##Run the test
if __name__ == '__main__':
   unittest.main()
