from sys import path
import os
import unittest
import numpy as np

##Do some disgusting path finding exercise, there must be a better way
##to do this
fileloc = os.path.realpath(__file__)
path.append('{:s}/../../src'.format(('/').join(fileloc.split('/')[:-1])))

##Code we are testing
import run_woden as rw

##Vehicle for running tests
class Test(unittest.TestCase):
    def test_make_antenna_table(self):
        """Tests the `rw.make_antenna_table` function, which should create
        the antenna table that goes into a uvfits file. Test by giving it a
        known set of input parameters and checking those values end up in
        the correct location and format in the output antenna table"""

        ##Some input test params================================================
        date = "2019-06-12T13:04:12"
        num_antennas = 10
        freq_cent = 160e+6
        telescope_name = "MWA"
        longitude = 116.670813889
        latitude = -26.7033194444
        array_height = 377.0
        gst0_deg = 260.0303917560829
        degpdy = 360.98563614850775
        ut1utc = -0.17421478816666666

        XYZ_array = np.arange(num_antennas*3)
        XYZ_array.shape = (num_antennas, 3)

        ##Code being tested=====================================================
        ant_table = rw.make_antenna_table(XYZ_array=XYZ_array,
                      telescope_name=telescope_name,
                      num_antennas=num_antennas, freq_cent=freq_cent,
                      date=date, gst0_deg=gst0_deg, degpdy=degpdy,
                      ut1utc=ut1utc,
                      longitude=longitude, latitude=latitude,
                      array_height=array_height)

        ##Test the type of the output is correct
        from astropy.io.fits.hdu.table import BinTableHDU
        self.assertEqual(BinTableHDU, type(ant_table))

        ##Some expected output values
        annnames = np.array(['00001', '00002', '00003', '00004', '00005',
                             '00006', '00007', '00008', '00009', '00010'])

        ##This is a weirdo empty array type that seems to hang around in a
        ##uvfits file
        nothing_array = np.array([], dtype='|V8')
        nothing_array.shape = (0,)

        ##Test that a bunch of named columns in the output array are equal======
        self.assertTrue(np.array_equal(annnames, ant_table.data['ANNAME']))
        self.assertTrue(np.array_equal(XYZ_array, ant_table.data['STABXYZ']))
        self.assertTrue(np.array_equal(nothing_array, ant_table.data['ORBPARM']))
        self.assertTrue(np.array_equal(np.arange(1, num_antennas+1), ant_table.data['NOSTA']))
        self.assertTrue(np.array_equal(np.zeros(num_antennas), ant_table.data['MNTSTA']))
        self.assertTrue(np.array_equal(np.zeros(num_antennas), ant_table.data['STAXOF']))
        self.assertTrue(np.array_equal(np.array(['X']*num_antennas), ant_table.data['POLTYA']))
        self.assertTrue(np.array_equal(np.zeros(num_antennas), ant_table.data['POLAA']))
        self.assertTrue(np.array_equal(np.zeros(num_antennas), ant_table.data['POLCALA']))
        self.assertTrue(np.array_equal(np.array(['Y']*num_antennas), ant_table.data['POLTYB']))
        self.assertTrue(np.array_equal(np.zeros(num_antennas), ant_table.data['POLAB']))
        self.assertTrue(np.array_equal(np.zeros(num_antennas), ant_table.data['POLCALB']))

        self.assertAlmostEqual(ant_table.header['ARRAYX'], -2559453.6265158253)
        self.assertAlmostEqual(ant_table.header['ARRAYY'], 5095371.541942583)
        self.assertAlmostEqual(ant_table.header['ARRAYZ'], -2849056.817561836)
        self.assertEqual(ant_table.header['FREQ'], freq_cent)

        self.assertEqual(ant_table.header['GSTIA0'], gst0_deg)
        self.assertEqual(ant_table.header['DEGPDY'], degpdy)
        self.assertEqual(ant_table.header['UT1UTC'], ut1utc)
        self.assertEqual(ant_table.header['XYZHAND'], 'RIGHT')
        self.assertEqual(ant_table.header['FRAME'], '????')

        self.assertEqual(ant_table.header['RDATE'], date)
        self.assertEqual(ant_table.header['TIMSYS'], 'UTC')
        self.assertEqual(ant_table.header['ARRNAM'], telescope_name )
        self.assertEqual(ant_table.header['NUMORB'], 0)
        self.assertEqual(ant_table.header['NOPCAL'], 0)
        self.assertEqual(ant_table.header['POLTYPE'], '')
        self.assertEqual(ant_table.header['CREATOR'], 'WODEN_uvfits_writer')

##Run the test
if __name__ == '__main__':
   unittest.main()
