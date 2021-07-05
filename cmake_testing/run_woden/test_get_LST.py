from sys import path
import os
import unittest

##Do some disgusting path finding exercise, there must be a better way
##to do this
fileloc = os.path.realpath(__file__)
path.append('{:s}/../../src'.format(('/').join(fileloc.split('/')[:-1])))

##Code we are testing
import run_woden as rw

##Vehicle for running tests
class Test(unittest.TestCase):
    def test_get_LST(self):
        """Tests the `rw.get_LST` function, which should calculate the
        Local Sidereal Time for a given Long/Lat and UTC date"""

        ##Some input test params
        date1 = "2019-06-12T13:04:12"
        long1 = 95.2
        lat1 = 46.0
        height1 = 0.0

        date2 = "2013-01-12T09:45:09"
        long2 = 232.9
        lat2 = -26.7
        height2 = 200.0

        ##Code to be tested
        lst_loc1date1 = rw.get_LST(lat1, long1, date1, height1)
        lst_loc2date1 = rw.get_LST(lat2, long2, date1, height1)
        lst_loc1date2 = rw.get_LST(lat1, long1, date2, height2)
        lst_loc2date2 = rw.get_LST(lat2, long2, date2, height2)

        ##Test answers are as expected
        self.assertAlmostEqual(191.817149246, lst_loc1date1, delta=1e-9)
        self.assertAlmostEqual(329.517149246, lst_loc2date1, delta=1e-9)
        self.assertAlmostEqual(353.542251851, lst_loc1date2, delta=1e-9)
        self.assertAlmostEqual(131.242251851, lst_loc2date2, delta=1e-9)

##Run the test
if __name__ == '__main__':
   unittest.main()
