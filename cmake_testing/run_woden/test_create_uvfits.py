from sys import path
import os
import unittest
import numpy as np
from astropy.io import fits

np.random.seed(179)

##Do some disgusting path finding exercise, there must be a better way
##to do this
fileloc = os.path.realpath(__file__)
path.append('{:s}/../../src'.format(('/').join(fileloc.split('/')[:-1])))

##Code we are testing
import run_woden as rw

##Vehicle for running tests
class Test(unittest.TestCase):
    def make_dummy_intputs(self):
        """Makes some inputs used to pass to rw.create_uvfits"""
        self.date = "2019-06-12T13:04:12"
        self.num_antennas = 10
        self.freq_cent = 160e+6
        self.telescope_name = "TESTTELE"
        self.longitude = 116.670813889
        self.latitude = -26.7033194444
        self.array_height = 377.0
        self.central_freq_chan = 1
        self.ch_width = 40e3
        self.ra_point = 0.0
        self.dec_point = -26.7
        self.num_baselines = int(((self.num_antennas - 1)*self.num_antennas) / 2)
        self.num_time_steps = 2
        self.num_freq_channels = 3
        self.int_jd = 2458647.0
        self.gitlabel = 'as987has'
        self.XYZ_array = np.arange(self.num_antennas*3)
        self.XYZ_array.shape = (self.num_antennas, 3)

        self.output_uvfits_name = "unittest_example.uvfits"


    def test_create_uvfits_outputs(self):
        """Tests the `rw.create_uvfits` function, which should take a whole
        heap of inputs and write out a uvits. Test by running funciton,
        reading in the created file and checking contents"""

        ##Some input test params
        self.make_dummy_intputs()

        ##Create an antenna table
        ant_table = rw.make_antenna_table(XYZ_array=self.XYZ_array,
                      telescope_name=self.telescope_name,
                      num_antennas=self.num_antennas, freq_cent=self.freq_cent,
                      date=self.date,
                      longitude=self.longitude, latitude=self.latitude,
                      array_height=self.array_height)

        ##Create some dummy data inputs
        v_container = np.empty((self.num_time_steps*self.num_baselines, 1, 1, self.num_freq_channels, 4, 3))
        random_data = np.random.uniform(-10000, 10000, self.num_time_steps*self.num_baselines*self.num_freq_channels*4*3)
        random_data.shape = (self.num_time_steps*self.num_baselines, self.num_freq_channels, 4, 3)
        v_container[:,0,0,:,:,:] = random_data

        uu = np.arange(self.num_baselines*self.num_time_steps)
        vv = np.arange(self.num_baselines*self.num_time_steps) + self.num_baselines
        ww = np.arange(self.num_baselines*self.num_time_steps) + 2*self.num_baselines
        baselines_array = np.arange(self.num_baselines*self.num_time_steps)
        date_array = np.ones(self.num_time_steps*self.num_baselines)*0.04458333319
        ##Add a second of time onto the second half of the date array
        date_array[self.num_baselines:] += (1 / (24.0*60.0*60.0))

        rw.create_uvfits(freq_cent=self.freq_cent, ch_width=self.ch_width,
                  central_freq_chan=self.central_freq_chan,
                  ra_point=self.ra_point, dec_point=self.dec_point,
                  output_uvfits_name=self.output_uvfits_name,
                  int_jd=self.int_jd, gitlabel=self.gitlabel,
                  v_container=v_container,
                  uu=uu, vv=vv, ww=ww,
                  baselines_array=baselines_array, date_array=date_array,
                  hdu_ant=ant_table)

        self.assertTrue(os.path.isfile(self.output_uvfits_name))

        with fits.open(self.output_uvfits_name) as hdu:
            ant_table = hdu[1]
            data_table = hdu[0]

            ##Some expected output values
            annnames = np.array(['00001', '00002', '00003', '00004', '00005',
                                 '00006', '00007', '00008', '00009', '00010'])

            ##This is a weirdo empty array type that seems to hang around in a
            ##uvfits file
            nothing_array = np.array([], dtype='|V8')
            nothing_array.shape = (0,)

            ##Test that a bunch of named columns in the antenna table are correct
            self.assertTrue(np.array_equal(annnames, ant_table.data['ANNAME']))
            self.assertTrue(np.array_equal(self.XYZ_array, ant_table.data['STABXYZ']))
            self.assertTrue(np.array_equal(nothing_array, ant_table.data['ORBPARM']))
            self.assertTrue(np.array_equal(np.arange(1, self.num_antennas+1), ant_table.data['NOSTA']))
            self.assertTrue(np.array_equal(np.zeros(self.num_antennas), ant_table.data['MNTSTA']))
            self.assertTrue(np.array_equal(np.zeros(self.num_antennas), ant_table.data['STAXOF']))
            self.assertTrue(np.array_equal(np.array(['X']*self.num_antennas), ant_table.data['POLTYA']))
            self.assertTrue(np.array_equal(np.zeros(self.num_antennas), ant_table.data['POLAA']))
            self.assertTrue(np.array_equal(np.zeros(self.num_antennas), ant_table.data['POLCALA']))
            self.assertTrue(np.array_equal(np.array(['Y']*self.num_antennas), ant_table.data['POLTYB']))
            self.assertTrue(np.array_equal(np.zeros(self.num_antennas), ant_table.data['POLAB']))
            self.assertTrue(np.array_equal(np.zeros(self.num_antennas), ant_table.data['POLCALB']))
            ##Check the header of the antenna table
            self.assertAlmostEqual(ant_table.header['ARRAYX'], -2559453.6265158253)
            self.assertAlmostEqual(ant_table.header['ARRAYY'], 5095371.541942583)
            self.assertAlmostEqual(ant_table.header['ARRAYZ'], -2849056.817561836)
            self.assertEqual(ant_table.header['FREQ'], self.freq_cent)
            self.assertEqual(ant_table.header['RDATE'], self.date)
            self.assertEqual(ant_table.header['TIMSYS'], 'UTC')
            self.assertEqual(ant_table.header['ARRNAM'], self.telescope_name )
            self.assertEqual(ant_table.header['NUMORB'], 0)
            self.assertEqual(ant_table.header['NOPCAL'], 0)
            self.assertEqual(ant_table.header['POLTYPE'], '')
            self.assertEqual(ant_table.header['CREATOR'], 'WODEN_uvfits_writer')

            ##Test that a bunch of named columns in the data table are correct
            self.assertTrue(np.array_equal(uu, data_table.data['UU']))
            self.assertTrue(np.array_equal(vv, data_table.data['VV']))
            self.assertTrue(np.array_equal(ww, data_table.data['WW']))
            self.assertTrue(np.array_equal(baselines_array, data_table.data['BASELINE']))
            ##Astropy automatically adds the header value to the DATE array,
            ##so need to subtract before comparison
            self.assertTrue(np.array_equal(date_array, data_table.data['DATE']) - data_table.header['PZERO5'])

            ##Check the actual visisbility values are correct
            self.assertTrue(np.allclose(v_container, data_table.data.data))

            ##Check the header of the data table
            self.assertEqual('COMPLEX', data_table.header['CTYPE2'])
            self.assertEqual(1.0, data_table.header['CRVAL2'])
            self.assertEqual(1.0, data_table.header['CRPIX2'])
            self.assertEqual(1.0, data_table.header['CDELT2'])
            self.assertEqual('STOKES', data_table.header['CTYPE3'])
            self.assertEqual(-5.0, data_table.header['CRVAL3'])
            self.assertEqual( 1.0, data_table.header['CRPIX3'])
            self.assertEqual(-1.0, data_table.header['CDELT3'])
            self.assertEqual('FREQ', data_table.header['CTYPE4'])
            self.assertEqual(self.freq_cent, data_table.header['CRVAL4'])
            self.assertEqual(int(self.central_freq_chan) + 1, data_table.header['CRPIX4'])
            self.assertEqual(self.ch_width, data_table.header['CDELT4'])
            self.assertEqual('RA', data_table.header['CTYPE5'])
            self.assertEqual(self.ra_point, data_table.header['CRVAL5'])
            self.assertEqual(1.0, data_table.header['CRPIX5'])
            self.assertEqual(1.0, data_table.header['CDELT5'])
            self.assertEqual('DEC', data_table.header['CTYPE6'])
            self.assertEqual(self.dec_point, data_table.header['CRVAL6'])
            self.assertEqual(1.0, data_table.header['CRPIX6'])
            self.assertEqual(1.0, data_table.header['CDELT6'])
            self.assertEqual(1.0, data_table.header['PSCAL1'])
            self.assertEqual(0.0, data_table.header['PZERO1'])
            self.assertEqual(1.0, data_table.header['PSCAL2'])
            self.assertEqual(0.0, data_table.header['PZERO2'])
            self.assertEqual(1.0, data_table.header['PSCAL3'])
            self.assertEqual(0.0, data_table.header['PZERO3'])
            self.assertEqual(1.0, data_table.header['PSCAL4'])
            self.assertEqual(0.0, data_table.header['PZERO4'])
            self.assertEqual(1.0, data_table.header['PSCAL5'])
            self.assertEqual(float(self.int_jd), data_table.header['PZERO5'])
            self.assertEqual('Undefined', data_table.header['OBJECT'])
            self.assertEqual(self.ra_point, data_table.header['OBSRA'])
            self.assertEqual(self.dec_point, data_table.header['OBSDEC'])
            self.assertEqual(self.gitlabel, data_table.header['GITLABEL'])

    def test_error_for_unequal_arrays(self):
        """Tests the `rw.create_uvfits` function. Function should fail with
        an error message if some arrays do not have equal dimensions. Check for
        failure and check the error message"""

        ##Some input test params
        self.make_dummy_intputs()

        ##Create an antenna table
        ant_table = rw.make_antenna_table(XYZ_array=self.XYZ_array,
                      telescope_name=self.telescope_name,
                      num_antennas=self.num_antennas, freq_cent=self.freq_cent,
                      date=self.date,
                      longitude=self.longitude, latitude=self.latitude,
                      array_height=self.array_height)

        ##Make unequal inputs to create failure
        v_container = np.zeros((10, 1, 1, 16, 4, 3))
        uu = np.zeros(1)
        vv = np.zeros(2)
        ww = np.zeros(3)
        baselines_array = np.zeros(4)
        date_array = np.zeros(5)

        ##Call the code
        with self.assertRaises(SystemExit) as cm:
            rw.create_uvfits(freq_cent=self.freq_cent, ch_width=self.ch_width,
                      central_freq_chan=self.central_freq_chan,
                      ra_point=self.ra_point, dec_point=self.dec_point,
                      output_uvfits_name=self.output_uvfits_name,
                      int_jd=self.int_jd, gitlabel=self.gitlabel,
                      v_container=v_container,
                      uu=uu, vv=vv, ww=ww,
                      baselines_array=baselines_array, date_array=date_array,
                      hdu_ant=ant_table)
        ##Check the error string is as expected
        expec_fail_string = ("run_woden.create_uvfits: The first dimension of the arrays:\n"
                            "v_container, uu, vv, ww, baselines_array, date_array\n"
                            "must be equal to make a uvfits file. Exiting now.")
        self.assertEqual(cm.exception.code, expec_fail_string)

##Run the test
if __name__ == '__main__':
   unittest.main()
