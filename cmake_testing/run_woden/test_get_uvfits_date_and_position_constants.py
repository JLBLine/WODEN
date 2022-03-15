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
    def test_uvfits_date_and_position_constants(self):
        """Tests the `rw.get_uvfits_date_and_position_constants` function,
        which should calculate the LST, GST0 (greenwich sidereal time at 0 hours
        of the given date), DEGPDY (rotational speed of the Earth) and UT1UTC (
        difference between UT1 and UTC) for a given Long/Lat and UTC date.
        """

        ##Some input test params
        date1 = "2019-06-12T13:04:12"
        long1 = 95.2
        lat1 = 46.0
        height1 = 0.0

        date2 = "2013-01-12T09:45:09"
        long2 = 232.9
        lat2 = -26.7
        height2 = 200.0

        ##Code to be tested=====================================================
        lst_loc1date1, gst0_loc1date1, degpdy_loc1date1, ut1utc_loc1date1 = rw.get_uvfits_date_and_position_constants(lat1, long1, date1, height1)
        lst_loc2date1, gst0_loc2date1, degpdy_loc2date1, ut1utc_loc2date1 = rw.get_uvfits_date_and_position_constants(lat2, long2, date1, height1)
        lst_loc1date2, gst0_loc1date2, degpdy_loc1date2, ut1utc_loc1date2 = rw.get_uvfits_date_and_position_constants(lat1, long1, date2, height2)
        lst_loc2date2, gst0_loc2date2, degpdy_loc2date2, ut1utc_loc2date2 = rw.get_uvfits_date_and_position_constants(lat2, long2, date2, height2)

        ##Test answers are as expected==========================================

        gst0_1 = 260.0303917560829
        degpdy1 = 360.98563614850775
        ut1utc1 = -0.17421478816666666
        gst0_2 = 111.65422226828267
        degpdy2 = 360.98565484299667
        ut1utc2 = 0.26520905895833335

        delta = 3e-7

        ##Test LST
        self.assertAlmostEqual(191.817149246, lst_loc1date1, delta=delta)
        self.assertAlmostEqual(329.517149246, lst_loc2date1, delta=delta)
        self.assertAlmostEqual(353.542251851, lst_loc1date2, delta=delta)
        self.assertAlmostEqual(131.242251851, lst_loc2date2, delta=delta)

        ##Test GST0
        self.assertAlmostEqual(gst0_1, gst0_loc1date1, delta=delta)
        self.assertAlmostEqual(gst0_1, gst0_loc2date1, delta=delta)
        self.assertAlmostEqual(gst0_2, gst0_loc1date2, delta=delta)
        self.assertAlmostEqual(gst0_2, gst0_loc2date2, delta=delta)

        ##Test DEGPDY
        self.assertAlmostEqual(degpdy1, degpdy_loc1date1, delta=delta)
        self.assertAlmostEqual(degpdy1, degpdy_loc2date1, delta=delta)
        self.assertAlmostEqual(degpdy2, degpdy_loc1date2, delta=delta)
        self.assertAlmostEqual(degpdy2, degpdy_loc2date2, delta=delta)

        ##Astropy varies wildly from version to version it seems
        ##Up the delta here
        delta = 4e-5
        ##Test UT1UTC
        self.assertAlmostEqual(ut1utc1, ut1utc_loc1date1, delta=delta)
        self.assertAlmostEqual(ut1utc1, ut1utc_loc2date1, delta=delta)
        self.assertAlmostEqual(ut1utc2, ut1utc_loc1date2, delta=delta)
        self.assertAlmostEqual(ut1utc2, ut1utc_loc2date2, delta=delta)

##Run the test
if __name__ == '__main__':
   unittest.main()
