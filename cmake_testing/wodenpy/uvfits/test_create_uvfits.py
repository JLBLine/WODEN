from sys import path
import os
import unittest
import numpy as np
from astropy.io import fits
from copy import deepcopy

np.random.seed(179)

from wodenpy.uvfits import wodenpy_uvfits

##Vehicle for running tests
class Test(unittest.TestCase):

    def make_dummy_intputs(self):
        """Makes some inputs used to pass to wodenpy_uvfits.create_uvfits"""
        self.date = "2019-06-12T13:04:12"
        self.gst0_deg = 260.0303917560829
        self.degpdy = 360.98563614850775
        self.ut1utc = -0.17421478816666666
        self.num_antennas = 10
        self.freq_cent = 160e+6
        self.telescope_name = "MWA"
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
        self.jd_midnight = 2458646.5
        self.gitlabel = 'as987has'
        self.XYZ_array = np.arange(self.num_antennas*3)
        self.XYZ_array.shape = (self.num_antennas, 3)

        self.output_uvfits_name = "unittest_example.uvfits"
        

    def setup_inputs_and_make_ant_table(self):
        """Setup some simple inputs and make the antenna table"""


         ##Some input test params
        self.make_dummy_intputs()

        ##Create an antenna table
        self.ant_table = wodenpy_uvfits.make_antenna_table(XYZ_array=self.XYZ_array,
                      telescope_name=self.telescope_name,
                      num_antennas=self.num_antennas, freq_cent=self.freq_cent,
                      date=self.date, gst0_deg=self.gst0_deg,
                      degpdy=self.degpdy, ut1utc=self.ut1utc,
                      longitude=self.longitude, latitude=self.latitude,
                      array_height=self.array_height)

    def check_output_uvfits(self, IAU_order=False):
        """One the uvfits file has been written out, test it's been written
        correctly"""
        self.assertTrue(os.path.isfile(self.output_uvfits_name))

        with fits.open(self.output_uvfits_name) as hdu:
            ant_table = hdu[1]
            data_table = hdu[0]

            ##Some expected output values
            annnames = np.array(['00001', '00002', '00003', '00004', '00005',
                                 '00006', '00007', '00008', '00009', '00010'])

            ##ant_table.data['ORBPARM'] is a weirdo empty array type that seems
            ##to hang around in a uvfits file
            ##At some point, the behaviour of astropy changes drastically so
            ##this file gets read in as a totally different data type and shape
            ##So based on what mood astropy is in, create an array to test against
            ##SIIIIIIIIIIIIIIIIIGH

            if ant_table.data['ORBPARM'].dtype == "float64":
                nothing_array = np.array([], dtype='float64')
                nothing_array.shape = (0,0)

            else:
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
            self.assertAlmostEqual(ant_table.header['GSTIA0'], self.gst0_deg, delta=1e-9)
            self.assertAlmostEqual(ant_table.header['DEGPDY'], self.degpdy, delta=1e-9)
            self.assertAlmostEqual(ant_table.header['UT1UTC'], self.ut1utc, delta=1e-9)
            self.assertEqual(ant_table.header['XYZHAND'], 'RIGHT')
            self.assertEqual(ant_table.header['FRAME'], '????')
            self.assertEqual(ant_table.header['RDATE'], self.date)
            self.assertEqual(ant_table.header['TIMSYS'], 'UTC')
            self.assertEqual(ant_table.header['ARRNAM'], self.telescope_name )
            self.assertEqual(ant_table.header['NUMORB'], 0)
            self.assertEqual(ant_table.header['NOPCAL'], 0)
            self.assertEqual(ant_table.header['POLTYPE'], '')
            self.assertEqual(ant_table.header['CREATOR'], 'WODEN_uvfits_writer')

            ##Test that a bunch of named columns in the data table are correct
            self.assertTrue(np.array_equal(self.uu, data_table.data['UU']))
            self.assertTrue(np.array_equal(self.vv, data_table.data['VV']))
            self.assertTrue(np.array_equal(self.ww, data_table.data['WW']))
            self.assertTrue(np.array_equal(self.baselines_array, data_table.data['BASELINE']))
            ##Astropy automatically adds the header value to the DATE array,
            ##so need to subtract before comparison
            self.assertTrue(np.array_equal(self.date_array, data_table.data['DATE']) - data_table.header['PZERO4'])

            if IAU_order:
                
                self.assertTrue(np.allclose(self.v_container, data_table.data.data))
                self.assertEqual(data_table.header['IAUORDER'], True)
            else:
                ##Check the actual visisbility values are correct and
                ##reordered compared to what was input
                self.assertTrue(np.allclose(self.orig_v_container[:,:,:,:,0,:],
                               data_table.data.data[:,:,:,:,1,:]))
                self.assertTrue(np.allclose(self.orig_v_container[:,:,:,:,1,:],
                               data_table.data.data[:,:,:,:,0,:]))
                self.assertTrue(np.allclose(self.orig_v_container[:,:,:,:,2,:],
                               data_table.data.data[:,:,:,:,3,:]))
                self.assertTrue(np.allclose(self.orig_v_container[:,:,:,:,3,:],
                               data_table.data.data[:,:,:,:,2,:]))
                self.assertEqual(data_table.header['IAUORDER'], False)

            

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
            ##This is the date array so has an offset for JD date
            self.assertEqual(float(self.jd_midnight), data_table.header['PZERO4'])
            self.assertEqual(1.0, data_table.header['PSCAL5'])
            self.assertEqual(0.0, data_table.header['PZERO5'])

            self.assertEqual('Undefined', data_table.header['OBJECT'])
            self.assertEqual(self.ra_point, data_table.header['OBSRA'])
            self.assertEqual(self.dec_point, data_table.header['OBSDEC'])
            self.assertEqual(self.gitlabel, data_table.header['GITLABEL'])

            self.assertEqual(data_table.header['TELESCOP'], self.telescope_name)
            self.assertEqual(data_table.header['LAT'], self.latitude)
            self.assertEqual(data_table.header['LON'], self.longitude)
            self.assertEqual(data_table.header['ALT'], self.array_height)
            self.assertEqual(data_table.header['INSTRUME'], self.telescope_name)

    def make_dummy_data(self):
        """"Fill out enough data fields to write out a uvfits. Make a copy
        of the input data to test the IAU reordering"""
        ##Create some dummy data inputs
        self.v_container = np.empty((self.num_time_steps*self.num_baselines, 1, 1, self.num_freq_channels, 4, 3))
        random_data = np.random.uniform(-10000, 10000, self.num_time_steps*self.num_baselines*self.num_freq_channels*4*3)
        random_data.shape = (self.num_time_steps*self.num_baselines, self.num_freq_channels, 4, 3)
        self.v_container[:,0,0,:,:,:] = random_data

        self.orig_v_container = deepcopy(self.v_container)

        self.uu = np.arange(self.num_baselines*self.num_time_steps)
        self.vv = np.arange(self.num_baselines*self.num_time_steps) + self.num_baselines
        self.ww = np.arange(self.num_baselines*self.num_time_steps) + 2*self.num_baselines
        self.baselines_array = np.arange(self.num_baselines*self.num_time_steps)

        self.date_array = np.ones(self.num_time_steps*self.num_baselines)*0.04458333319
        ##Add a second of time onto the second half of the date array
        self.date_array[self.num_baselines:] += (1 / (24.0*60.0*60.0))

    def test_create_uvfits_outputs(self):
        """Tests the `wodenpy_uvfits.create_uvfits` function, which should take a whole
        heap of inputs and write out a uvits. Test by running funciton,
        reading in the created file and checking contents"""

        self.setup_inputs_and_make_ant_table()

        self.make_dummy_data()

        wodenpy_uvfits.create_uvfits(freq_cent=self.freq_cent, ch_width=self.ch_width,
                  central_freq_chan=self.central_freq_chan,
                  ra_point=self.ra_point, dec_point=self.dec_point,
                  output_uvfits_name=self.output_uvfits_name,
                  jd_midnight=self.jd_midnight, gitlabel=self.gitlabel,
                  v_container=self.v_container,
                  uu=self.uu, vv=self.vv, ww=self.ww,
                  baselines_array=self.baselines_array, date_array=self.date_array,
                  hdu_ant=self.ant_table, telescope_name=self.telescope_name,
                  longitude=self.longitude, latitude=self.latitude,
                  array_height=self.array_height)

        self.check_output_uvfits()

    def test_create_uvfits_outputs_IAU_order(self):
        """Tests the `wodenpy_uvfits.create_uvfits` function, which should take a whole
        heap of inputs and write out a uvits. Test by running funciton,
        reading in the created file and checking contents"""

        self.setup_inputs_and_make_ant_table()

        self.make_dummy_data()

        wodenpy_uvfits.create_uvfits(freq_cent=self.freq_cent, ch_width=self.ch_width,
                  central_freq_chan=self.central_freq_chan,
                  ra_point=self.ra_point, dec_point=self.dec_point,
                  output_uvfits_name=self.output_uvfits_name,
                  jd_midnight=self.jd_midnight, gitlabel=self.gitlabel,
                  v_container=self.v_container,
                  uu=self.uu, vv=self.vv, ww=self.ww,
                  baselines_array=self.baselines_array, date_array=self.date_array,
                  hdu_ant=self.ant_table, telescope_name=self.telescope_name,
                  longitude=self.longitude, latitude=self.latitude,
                  array_height=self.array_height,
                  IAU_order=True)

        self.check_output_uvfits(IAU_order=True)

    def test_error_for_unequal_arrays(self):
        """Tests the `wodenpy_uvfits.create_uvfits` function. Function should fail with
        an error message if some arrays do not have equal dimensions. Check for
        failure and check the error message"""

        ##Some input test params
        self.make_dummy_intputs()

        ##Create an antenna table
        self.ant_table = wodenpy_uvfits.make_antenna_table(XYZ_array=self.XYZ_array,
                      telescope_name=self.telescope_name,
                      num_antennas=self.num_antennas, freq_cent=self.freq_cent,
                      date=self.date,
                      longitude=self.longitude, latitude=self.latitude,
                      array_height=self.array_height,
                      gst0_deg=self.gst0_deg, ut1utc=self.ut1utc,
                      degpdy=self.degpdy)

        ##Make unequal inputs to create failure
        self.v_container = np.zeros((10, 1, 1, 16, 4, 3))
        self.uu = np.zeros(1)
        self.vv = np.zeros(2)
        self.ww = np.zeros(3)
        self.baselines_array = np.zeros(4)
        self.date_array = np.zeros(5)

        ##Call the code
        with self.assertRaises(SystemExit) as cm:
            wodenpy_uvfits.create_uvfits(freq_cent=self.freq_cent, ch_width=self.ch_width,
                      central_freq_chan=self.central_freq_chan,
                      ra_point=self.ra_point, dec_point=self.dec_point,
                      output_uvfits_name=self.output_uvfits_name,
                      jd_midnight=self.jd_midnight, gitlabel=self.gitlabel,
                      v_container=self.v_container,
                      uu=self.uu, vv=self.vv, ww=self.ww,
                      baselines_array=self.baselines_array, date_array=self.date_array,
                      hdu_ant=self.ant_table, longitude=self.longitude,
                      latitude=self.latitude, array_height=self.array_height,
                      telescope_name=self.telescope_name)
        ##Check the error string is as expected
        expec_fail_string = ("run_woden.create_uvfits: The first dimension of the arrays:\n"
                            "v_container, uu, vv, ww, baselines_array, date_array\n"
                            "must be equal to make a uvfits file. Exiting now.")
        self.assertEqual(cm.exception.code, expec_fail_string)

##Run the test
if __name__ == '__main__':
   unittest.main()
