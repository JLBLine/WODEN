from sys import path
import os
import unittest
import numpy as np
import ctypes

# ##Code we are testing
import wodenpy.use_libwoden.woden_settings as ws
from wodenpy.wodenpy_setup import run_setup
from wodenpy.use_libwoden.create_woden_struct_classes import Woden_Struct_Classes
import numpy.testing as npt

code_dir = os.path.realpath(__file__)
code_dir = ('/').join(code_dir.split('/')[:-1])

D2R = np.pi/180.0

class PretendArgs():
    """Function `ws.fill_woden_settings_python` takes an arg class out of argparse as an
    input. Write a dummy one here to feed in and test"""
    def make_basic_args(self):
        self.ra0 = 0.0
        self.dec0 = -26.7
        self.chunking_size = 1e+10
        self.coarse_band_width = 1.28e+6
        self.sky_crop_components = True
        self.sky_crop_sources = False
        self.primary_beam = "none"

        self.num_freq_channels = 16
        self.num_time_steps = 14
        self.time_res = 2.0
        self.freq_res = 40e+3

        self.array_layout = 'WODEN_array_layout.txt'
        self.array_layout_name = 'bespoke_array_layout.txt'
        self.lowest_channel_freq = 160e+6
        self.latitude = -26.7
        self.band_nums = [1,8,10]

        ##Set the optional arguments to False
        self.gauss_beam_FWHM = False
        self.gauss_beam_ref_freq = False
        self.gauss_ra_point = False
        self.gauss_dec_point = False
        self.MWA_FEE_delays = False
        self.no_precession = False
        self.do_autos = False
        self.use_MWA_dipamps = False
        self.use_dipflags = False
        self.dipamps = False
        
        self.hdf5_beam_path = ""
        self.station_id = np.nan
        
        self.off_cardinal_dipoles = False
        self.cpu_mode = False
        self.verbose = False
        self.no_beam_normalisation = False
        self.beam_ms_path = ""
        
    def __init__(self):
        self.make_basic_args()

##Vehicle for running tests
class Test(unittest.TestCase):
    """Test the `rw.write_json` function, which writes an input file to feed
    into the WODEN executable"""
    
    def make_basic_inputs(self, precision):
        """Make some basis input arguments for `rw.write_json`"""
        self.jd_date = 2458646.544583333
        # self.json_name = json_name
        self.lst = 10.0
        self.args = PretendArgs()
        self.args.precision = precision
        self.precision = precision
        
    def check_basic_inputs(self, woden_settings):
        """Check the basic outputs are set correctly"""

        if self.args.precision == 'float':
            delta = 1e-8
        else:
            delta = 1e-12
            
        self.assertAlmostEqual(self.args.ra0*D2R, woden_settings.ra0, delta=delta)
        self.assertAlmostEqual(self.args.dec0*D2R, woden_settings.dec0, delta=delta)
        self.assertEqual(self.args.num_freq_channels, woden_settings.num_freqs)
        self.assertEqual(self.args.num_time_steps, woden_settings.num_time_steps)
        self.assertEqual(self.args.time_res, woden_settings.time_res)
        self.assertEqual(self.args.freq_res, woden_settings.frequency_resolution)
        self.assertEqual(self.args.chunking_size, woden_settings.chunking_size)
        self.assertEqual(self.jd_date, woden_settings.jd_date)
        self.assertAlmostEqual(self.lst*D2R, woden_settings.lst_base, delta=delta)
        self.assertAlmostEqual(self.lst*D2R, woden_settings.lst_obs_epoch_base, delta=delta)
        self.assertEqual(self.args.lowest_channel_freq, woden_settings.base_low_freq)
        self.assertAlmostEqual(self.args.latitude*D2R, woden_settings.latitude, delta=delta)
        self.assertAlmostEqual(self.args.latitude*D2R, woden_settings.latitude_obs_epoch_base, delta=delta)
        self.assertEqual(self.args.coarse_band_width, woden_settings.coarse_band_width)
        self.assertEqual(len(self.args.band_nums), woden_settings.num_bands)
        
        npt.assert_allclose(self.args.band_nums, woden_settings.band_nums, atol=1e-8)
        
        # self.data = data

    def call_fill_woden_settings_python(self):
        """Calls the function under test"""
        woden_settings = ws.fill_woden_settings_python(args=self.args,
                                   lst=self.lst, jd_date=self.jd_date)
        return woden_settings
    

    def test_make_minimum_settings(self):
        """Run the bare mimimum set of choices """
        
        for precision in ['float', 'double']:
            self.make_basic_inputs(precision)
            woden_settings = self.call_fill_woden_settings_python()
            self.check_basic_inputs(woden_settings)
        
    def test_write_gaussian_beam(self):
        """Test that the Gaussian primary beam options work correctly"""
        
        def add_extra_gauss_args():
            ##Add in primary_beam as Gaussian to inputs
            self.args.primary_beam = 'Gaussian'
            self.args.gauss_beam_FWHM = 20.0
            self.args.gauss_beam_ref_freq = 180e+6
            self.args.gauss_ra_point = 60.0
            self.args.gauss_dec_point = -10.0
            
        def check_gauss_outputs(woden_settings):
            ##Check the extra arguments have resulting in correct outputs
            ##Optional arguments should now exists
            
            if self.args.precision == 'float':
                delta = 1e-8
            else:
                delta = 1e-10
                
            self.assertEqual(1, woden_settings.beamtype)
            self.assertAlmostEqual(self.args.gauss_ra_point,
                             woden_settings.gauss_ra_point/D2R,
                             delta=delta)
            self.assertAlmostEqual(self.args.gauss_dec_point,
                             woden_settings.gauss_dec_point/D2R,
                             delta=delta)
            self.assertAlmostEqual(self.args.gauss_beam_FWHM,
                             woden_settings.gauss_beam_FWHM, delta=delta)
            self.assertEqual(self.args.gauss_beam_ref_freq,
                             woden_settings.gauss_beam_ref_freq)


        for precision in ['float', 'double']:
            ##This makes fake args with double precision
            self.make_basic_inputs(precision)
            add_extra_gauss_args()
            
            ##This runs `fill_woden_settings_python`
            woden_settings = self.call_fill_woden_settings_python()

            ##This passes woden_settings into C code which write contents to a text
            ## file, reads in that text file, and checks the outputs make sense
            self.check_basic_inputs(woden_settings)
            check_gauss_outputs(woden_settings)

    def check_mwa_beam_delays(self, woden_settings):
        
        ##turn string of list into an actual list
        FEE_delays = [int(delay) for delay in self.args.MWA_FEE_delays.strip("[]").split(',')]
        
        # ##Check the extra arguments have resulting in correct outputs
        npt.assert_array_equal(FEE_delays, woden_settings.FEE_ideal_delays)

    def test_write_MWA_FEE_beam(self):
        """Test that the MWA FEE primary beam options work correctly"""

        def add_extra_MWAFEE_args():
            ##Add in primary_beam as Gaussian to inputs
            self.args.primary_beam = 'MWA_FEE'
            self.args.hdf5_beam_path = "This_is_the_way"
            self.args.MWA_FEE_delays = "[0,2,4,6,0,2,4,6,0,2,4,6,0,2,4,6]"
            
        def check_MWAFEE_outputs():
            self.assertEqual(2, woden_settings.beamtype)
            self.assertEqual(self.args.hdf5_beam_path, woden_settings.hdf5_beam_path)
            self.check_mwa_beam_delays(woden_settings)
            
        for precision in ['float', 'double']:
            
        
            ##This makes fake args with double precision
            self.make_basic_inputs(precision)
            add_extra_MWAFEE_args()
            
            ##This runs `fill_woden_settings_python`
            woden_settings = self.call_fill_woden_settings_python()

            ##This passes woden_settings into C code which write contents to a text
            ## file, reads in that text file, and checks the outputs make sense
            self.check_basic_inputs(woden_settings)
            check_MWAFEE_outputs()
        
    def test_write_MWA_FEE_beam_interp(self):
        """Test that the interpolated MWA FEE primary beam options work correctly"""

        def add_extra_MWAFEE_args():
            ##Add in primary_beam as Gaussian to inputs
            self.args.primary_beam = 'MWA_FEE_interp'
            self.args.hdf5_beam_path = "This_is_the_way_interped"
            self.args.MWA_FEE_delays = "[0,2,4,6,0,2,4,6,0,2,4,6,0,2,4,6]"
            
        def check_MWAFEE_outputs():
            self.assertEqual(4, woden_settings.beamtype)
            self.assertEqual(self.args.hdf5_beam_path, woden_settings.hdf5_beam_path)
            self.check_mwa_beam_delays(woden_settings)
            
        
        for precision in ['float', 'double']:
            ##This makes fake args with double precision
            self.make_basic_inputs(precision)
            add_extra_MWAFEE_args()
            
            ##This runs `fill_woden_settings_python`
            woden_settings = self.call_fill_woden_settings_python()

            ##This passes woden_settings into C code which write contents to a text
            ## file, reads in that text file, and checks the outputs make sense
            self.check_basic_inputs(woden_settings)
            check_MWAFEE_outputs()
        
    def test_write_MWA_analy_beam(self):
        """Test that the interpolated analytic MWA primary beam options work correctly"""

        def add_extra_MWAanaly_args():
            ##Add in primary_beam as Gaussian to inputs
            self.args.primary_beam = 'MWA_analy'
            self.args.MWA_FEE_delays = "[0,2,4,6,0,2,4,6,0,2,4,6,0,2,4,6]"
            
        def check_MWAanaly_outputs():
            self.assertEqual(5, woden_settings.beamtype)
            self.check_mwa_beam_delays(woden_settings)
            
        for precision in ['float', 'double']:
            ##This makes fake args with double precision
            self.make_basic_inputs(precision)
            add_extra_MWAanaly_args()
            
            ##This runs `fill_woden_settings_python`
            woden_settings = self.call_fill_woden_settings_python()

            ##This passes woden_settings into C code which write contents to a text
            ## file, reads in that text file, and checks the outputs make sense
            self.check_basic_inputs(woden_settings)
            check_MWAanaly_outputs()
            
    def test_write_EDA2_beam(self):
        """Test that the EDA2 primary beam options work correctly"""
        
        for precision in ['float', 'double']:
            ##This makes fake args with double precision
            self.make_basic_inputs(precision)
            self.args.primary_beam = 'EDA2'
            
            ##This runs `fill_woden_settings_python`
            woden_settings = self.call_fill_woden_settings_python()

            ##This passes woden_settings into C code which write contents to a text
            ## file, reads in that text file, and checks the outputs make sense
            self.check_basic_inputs(woden_settings)
            self.assertEqual(3, woden_settings.beamtype)

    def test_write_do_autos(self):
        """Test that doing autos gets swithced on. No need to check for
        float or double, this is just changing an int"""

        ##This makes fake args with double precision
        self.make_basic_inputs('double')
        self.args.do_autos = True
        
        ##This runs `fill_woden_settings_python`
        woden_settings = self.call_fill_woden_settings_python()
        
        ##This passes woden_settings into C code which write contents to a text
        ## file, reads in that text file, and checks the outputs make sense
        self.check_basic_inputs(woden_settings)
        self.assertTrue(woden_settings.do_autos)
        
        self.args.do_autos = False
        woden_settings = self.call_fill_woden_settings_python()
        self.check_basic_inputs(woden_settings)
        self.assertFalse(woden_settings.do_autos)


        
    def test_write_do_precssion(self):
        """Test that doing precession gets swithced on/off. No need to check for
        float or double, this is just changing an int"""
        
        ##This makes fake args with double precision
        self.make_basic_inputs('double')
        self.args.no_precession = False
        
        ##This runs `fill_woden_settings_python`
        woden_settings = self.call_fill_woden_settings_python()
        
        ##This passes woden_settings into C code which write contents to a text
        ## file, reads in that text file, and checks the outputs make sense
        self.check_basic_inputs(woden_settings)
        self.assertTrue(woden_settings.do_precession)
        
        self.args.no_precession = True
        woden_settings = self.call_fill_woden_settings_python()
        self.check_basic_inputs(woden_settings)
        self.assertFalse(woden_settings.do_precession)
        
    def test_use_dipamps(self):
        """Test that the use_dipamp flag gets propagated correctly, and if
        true, the dipole amplitude array is correctly passed through"""

        ##This makes fake args with double precision
        self.make_basic_inputs('double')
        # self.args.no_precession = False
        
        ##This runs `fill_woden_settings_python`
        woden_settings = self.call_fill_woden_settings_python()
        
        ##This passes woden_settings into C code which write contents to a text
        ## file, reads in that text file, and checks the outputs make sense
        self.check_basic_inputs(woden_settings)
        
        ##check it defaults to False
        self.assertFalse(woden_settings.use_dipamps)
        
        
        ##Now ask for some dipole amps, provide some (as would happen using
        # run_setup.check_args) and check things are propagated correctly
        self.args.use_MWA_dipamps = True
        self.args.dipamps = np.arange(16)
        ##This runs `fill_woden_settings_python`
        woden_settings = self.call_fill_woden_settings_python()
        
        self.check_basic_inputs(woden_settings)
        ##check it we have a True flag
        self.assertTrue(woden_settings.use_dipamps)
        
        # mwa_dipole_amps = np.array([float(amp) for amp in self.data['mwa_dipole_amps'].split(',')])
        # npt.assert_allclose(mwa_dipole_amps, np.arange(16), atol=1e-8)
        npt.assert_allclose(woden_settings.mwa_dipole_amps, np.arange(16), atol=1e-8)
        
        
    def test_off_cardinal_dipoles(self):
        """Test that the use_dipamp flag gets propagated correctly, and if
        true, the dipole amplitude array is correctly passed through"""

        ##First up, check that the off_cardinal_dipoles flag off by default
        self.make_basic_inputs('double')
        woden_settings = self.call_fill_woden_settings_python()
        self.check_basic_inputs(woden_settings)
        self.assertFalse(woden_settings.off_cardinal_dipoles)
        
        ##Next, turn it on via the command line
        self.make_basic_inputs('double')
        self.args.off_cardinal_dipoles = True
        woden_settings = self.call_fill_woden_settings_python()
        self.check_basic_inputs(woden_settings)
        self.assertTrue(woden_settings.off_cardinal_dipoles)
        
        ##Next, check it gets switched on given the right primary beam
        self.make_basic_inputs('double')
        self.args.primary_beam = 'everybeam_LOFAR'
        woden_settings = self.call_fill_woden_settings_python()
        self.check_basic_inputs(woden_settings)
        self.assertTrue(woden_settings.off_cardinal_dipoles)
        
    def test_cpu_mode(self):
        """Test that the `cpu_mode` flag gets propagated correctly and sets
        the `do_gpu` flag correctly"""

        ##First up, check that the off_cardinal_dipoles flag off by default
        self.make_basic_inputs('double')
        woden_settings = self.call_fill_woden_settings_python()
        self.check_basic_inputs(woden_settings)
        self.assertTrue(woden_settings.do_gpu)
        
        ##Next, turn it on via the command line
        self.make_basic_inputs('double')
        self.args.cpu_mode = True
        woden_settings = self.call_fill_woden_settings_python()
        self.check_basic_inputs(woden_settings)
        self.assertFalse(woden_settings.do_gpu)
        
    def test_no_beam_norm(self):
        """Test that the `no_beam_normalistion` option gets propagated
        and sets the `normalise_primary_beam` flag correctly"""

        ##First up, check that the off_cardinal_dipoles flag off by default
        self.make_basic_inputs('double')
        woden_settings = self.call_fill_woden_settings_python()
        self.check_basic_inputs(woden_settings)
        self.assertTrue(woden_settings.normalise_primary_beam)
        
        ##Next, turn it on via the command line
        self.make_basic_inputs('double')
        self.args.no_beam_normalisation = True
        woden_settings = self.call_fill_woden_settings_python()
        self.check_basic_inputs(woden_settings)
        self.assertFalse(woden_settings.normalise_primary_beam)
        
##Run the test
if __name__ == '__main__':
   unittest.main()