from sys import path
import os
import unittest
import numpy as np
from astropy.io import fits
from unittest import mock
import argparse

code_dir = os.path.realpath(__file__)
code_dir = ('/').join(code_dir.split('/')[:-1])

##Append the location of add_woden_uvfits to the sys.path to import from it
path.append('{:s}/../../../scripts'.format(code_dir))

##Code we are testing
from wodenpy.uvfits import wodenpy_uvfits
import concat_woden_uvfits as cwu

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
        self.central_freq_chan = 2
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
        self.output_uvfits_name = "unittest_example1_band01.uvfits"

    def create_uvfits_outputs(self, visi_values=1, num_freqs=3):
        """Tests the `wodenpy_uvfits.create_uvfits` function, which should take a whole
        heap of inputs and write out a uvits. Test by running funciton,
        reading in the created file and checking contents"""

        self.num_freq_channels = num_freqs

        ##Create an antenna table
        ant_table = wodenpy_uvfits.make_antenna_table(XYZ_array=self.XYZ_array,
                      telescope_name=self.telescope_name,
                      num_antennas=self.num_antennas, freq_cent=self.freq_cent,
                      date=self.date, gst0_deg=self.gst0_deg,
                      degpdy=self.degpdy, ut1utc=self.ut1utc,
                      longitude=self.longitude, latitude=self.latitude,
                      array_height=self.array_height)

        ##Create some dummy data inputs
        v_container = np.ones((self.num_time_steps*self.num_baselines, 1, 1, self.num_freq_channels, 4, 3))*visi_values

        uu = np.arange(self.num_baselines*self.num_time_steps)
        vv = np.arange(self.num_baselines*self.num_time_steps) + self.num_baselines
        ww = np.arange(self.num_baselines*self.num_time_steps) + 2*self.num_baselines
        baselines_array = np.arange(self.num_baselines*self.num_time_steps)

        date_array = np.ones(self.num_time_steps*self.num_baselines)*0.04458333319
        ##Add a second of time onto the second half of the date array
        date_array[self.num_baselines:] += (1 / (24.0*60.0*60.0))

        wodenpy_uvfits.create_uvfits(freq_cent=self.freq_cent, ch_width=self.ch_width,
                  central_freq_chan=self.central_freq_chan,
                  ra_point=self.ra_point, dec_point=self.dec_point,
                  output_uvfits_name=self.output_uvfits_name,
                  int_jd=self.int_jd, gitlabel=self.gitlabel,
                  v_container=v_container,
                  uu=uu, vv=vv, ww=ww,
                  baselines_array=baselines_array, date_array=date_array,
                  hdu_ant=ant_table, telescope_name=self.telescope_name,
                  longitude=self.longitude, latitude=self.latitude,
                  array_height=self.array_height)

        self.uu = uu
        self.vv = vv
        self.ww = ww

        self.baselines_array = baselines_array
        self.date_array = date_array

    def check_uvfits_contents(self, uvfits_name, expec_data, expec_num_freqs):

        with fits.open(uvfits_name) as hdu:
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

            ##Check the actual visisbility values are correct
            self.assertTrue(np.allclose(expec_data, data_table.data.data))

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
            self.assertEqual(float(self.int_jd), data_table.header['PZERO4'])
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

    ##Setup some fake command line things for the following test
    @mock.patch('argparse.ArgumentParser.parse_args',
                return_value=argparse.Namespace(num_bands=4,
                            uvfits_prepend="unittest_example1_band",
                            output_name="concatenated_uvfuts.uvfits",
                            swap_pols=False))

    def test_concat_four_uvfits(self, mock_args):
        """Checks that concatenating four uvfits files that have `band{band}.uvfits`
        at the end works, by creating four uvfits with known properties"""

        ##Some input test params
        self.make_dummy_intputs()
        self.create_uvfits_outputs()

        expec_data = np.ones((self.num_time_steps*self.num_baselines, 1, 1, self.num_freq_channels*4, 4, 3))

        for band in range(2,5):
            cent_freq = self.freq_cent

            self.freq_cent += self.num_freq_channels*self.ch_width

            self.output_uvfits_name = f"unittest_example1_band{band:02d}.uvfits"
            self.create_uvfits_outputs(visi_values=band)

            low = self.num_freq_channels*(band - 1)
            high = self.num_freq_channels*band

            expec_data[:, 0, 0, low:high, :, :] = band

        cwu.main()

        ##Stick the frequencies settings to the first uvfits
        self.freq_cent = 160e+6
        self.check_uvfits_contents("concatenated_uvfuts.uvfits",
                              expec_data,
                              self.num_freq_channels*4)

    ##TODO haven't actually setup the code to check the dimensions but
    ##this is the kind of test you would run to check that behaviour if we do
    ##if in the le future
    # ##Setup some fake command line things for the following test
    # @mock.patch('argparse.ArgumentParser.parse_args',
    #             return_value=argparse.Namespace(num_bands=2,
    #                         uvfits_prepend="unittest_example1_band",
    #                         output_name="concatenated_uvfuts.uvfits",
    #                         swap_pols=False))
    # def test_mismatch_data_fails(self, mock_args):
    #     """Check everything fails when the dimensions of data in uvfits1
    #     and uvfits2 do not match"""
    
    # #     self.make_dummy_intputs()
    # #     self.output_uvfits_name = "unittest_diff_dimension.uvfits"
    # #     self.create_uvfits_outputs(num_freqs=2)

    #     ##Some input test params
    #     self.make_dummy_intputs()
    #     self.create_uvfits_outputs()

    #     band = 2
    #     self.output_uvfits_name = f"unittest_example1_band{band:02d}.uvfits"
    #     self.create_uvfits_outputs(visi_values=band, num_freqs=2)
    
    #     ##Contain things when they go wrong
    #     with self.assertRaises(SystemExit) as cm:
    #         ##This call to main should die
    #         cwu.main()
    
    #     expec_fail_string = "Shape of the data in the uvfits are not the same. Data "\
    #          "might have different number of times steps, frequencies, "\
    #          "or baselines."
    #     self.assertEqual(cm.exception.code, expec_fail_string)

##Run the test
if __name__ == '__main__':
   unittest.main()
