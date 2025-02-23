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

##annoying path hack to find where the C library is
test_dir = os.environ['CMAKE_CURRENT_SOURCE_DIR'] + "/../../../build/cmake_testing/wodenpy/use_libwoden/"

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
        self.num_antennas = 16
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
        
    def __init__(self):
        self.make_basic_args()

##Vehicle for running tests
class Test(unittest.TestCase):
    """Test the `rw.write_json` function, which writes an input file to feed
    into the WODEN executable"""
    
    
    def read_in_C_functions(self):
        """Read in the C library function that reads the ctypes woden_settings
        structure and writes the content to a text file"""
        
        woden_struct_classes = Woden_Struct_Classes(self.precision)
        self.woden_settings = woden_struct_classes.Woden_Settings()
        
        ## Read in the C library for float version
        libwoden = ctypes.cdll.LoadLibrary(f"{test_dir}/libread_woden_settings_{self.precision}.so")
        self.read_woden_settings = libwoden.read_woden_settings
        self.read_woden_settings.argtypes = [ctypes.POINTER(woden_struct_classes.Woden_Settings)]
        
    def make_basic_inputs(self, precision):
        """Make some basis input arguments for `rw.write_json`"""
        self.jd_date = 2458646.544583333
        # self.json_name = json_name
        self.lst = 10.0
        self.args = PretendArgs()
        self.args.precision = precision
        self.precision = precision
        
        ##load in the C library
        self.read_in_C_functions()
        
    def read_C_output_textfile(self):
        
        output_dict = {}
        
        with open('woden_settings.txt', 'r') as infile:
            for line in infile.read().split('\n'):
                if line == '':
                    pass
                else:
                    output_dict[line.split()[0]] = line.split()[1]
            
        return output_dict

    def check_basic_inputs(self, woden_settings):
        """Run C code that reads in the woden_settings struct and prints
        out the contents to a text file"""

        if self.args.precision == 'float':
            # print("Checking the FLOAT")
            # self.read_woden_settings_float(woden_settings)
            delta = 1e-8
        else:
            # print("Checking the DOUBLE")
            # self.read_woden_settings_double(woden_settings)
            delta = 1e-12
            
        self.read_woden_settings(woden_settings)
        data = self.read_C_output_textfile()

        self.assertAlmostEqual(self.args.ra0*D2R, float(data['ra0']), delta=delta)
        self.assertAlmostEqual(self.args.dec0*D2R, float(data['dec0']), delta=delta)
        self.assertEqual(self.args.num_freq_channels, int(data['num_freqs']))
        self.assertEqual(self.args.num_time_steps, int(data['num_time_steps']))
        self.assertEqual(self.args.time_res, float(data['time_res']))
        self.assertEqual(self.args.freq_res, float(data['frequency_resolution']))
        self.assertEqual(self.args.chunking_size, float(data['chunking_size']))
        self.assertEqual(self.jd_date, float(data['jd_date']))
        self.assertAlmostEqual(self.lst*D2R, float(data['lst_base']), delta=delta)
        self.assertAlmostEqual(self.lst*D2R, float(data['lst_obs_epoch_base']), delta=delta)
        self.assertEqual(self.args.lowest_channel_freq, float(data['base_low_freq']))
        self.assertAlmostEqual(self.args.latitude*D2R, float(data['latitude']), delta=delta)
        self.assertAlmostEqual(self.args.latitude*D2R, float(data['latitude_obs_epoch_base']), delta=delta)
        self.assertEqual(self.args.coarse_band_width, float(data['coarse_band_width']))

        self.assertEqual(len(self.args.band_nums), int(data['num_bands']))
        
        ##Loop through the output band numbers and stick into a list
        written_bands = []
        for key in data.keys():
            if key[:7] == 'bandnum':
                written_bands.append(int(data[key]))
                
        ##Check it matches
        self.assertEqual(self.args.band_nums, written_bands)
        
        self.data = data

    def call_fill_woden_settings_python(self):
        """Calls the function under test"""
        woden_settings_python = ws.fill_woden_settings_python(args=self.args,
                                 lst=self.lst, jd_date=self.jd_date)
        
        ##These things normally gets set by `setup_array_layout`, so set them
        ##here
        woden_settings_python.num_baselines = int((self.args.num_antennas*(self.args.num_antennas-1)) / 2)
        woden_settings_python.num_ants = self.args.num_antennas
        num_cross = woden_settings_python.num_baselines * woden_settings_python.num_time_steps * woden_settings_python.num_freqs
        woden_settings_python.num_cross = num_cross
        if woden_settings_python.do_autos:
            num_autos = woden_settings_python.num_ants * woden_settings_python.num_time_steps * woden_settings_python.num_freqs
        else:
            num_autos = 0
        
        woden_settings_python.num_autos = num_autos
        woden_settings_python.num_visis = num_cross + num_autos
        
        _, _ = ws.setup_lsts_and_phase_centre(woden_settings_python)
        ##This can end up doing precession, which is tested elsewhere. For simnplicity
        ##here convert LST back to the input
        
        woden_settings_python.lst_base = self.lst*D2R
        woden_settings_python.latitude = self.args.latitude*D2R
        
        woden_settings = ws.convert_woden_settings_to_ctypes(woden_settings_python,
                                                             self.woden_settings)
        
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
            
        def check_gauss_outputs():
            ##Check the extra arguments have resulting in correct outputs
            ##Optional arguments should now exists
            
            if self.args.precision == 'float':
                delta = 1e-8
            else:
                delta = 1e-10
                
            self.assertEqual(1, int(self.data['beamtype']))
            self.assertAlmostEqual(self.args.gauss_ra_point,
                             float(self.data['gauss_ra_point'])/D2R,
                             delta=delta)
            self.assertAlmostEqual(self.args.gauss_dec_point,
                             float(self.data['gauss_dec_point'])/D2R,
                             delta=delta)
            self.assertAlmostEqual(self.args.gauss_beam_FWHM,
                             float(self.data['gauss_beam_FWHM']), delta=delta)
            self.assertEqual(self.args.gauss_beam_ref_freq,
                             float(self.data['gauss_beam_ref_freq']))


        for precision in ['float', 'double']:
            ##This makes fake args with double precision
            self.make_basic_inputs(precision)
            add_extra_gauss_args()
            
            ##This runs `fill_woden_settings_python`
            woden_settings = self.call_fill_woden_settings_python()

            ##This passes woden_settings into C code which write contents to a text
            ## file, reads in that text file, and checks the outputs make sense
            self.check_basic_inputs(woden_settings)
            check_gauss_outputs()

    def check_mwa_beam_delays(self):
        
        ##turn string of list into an actual list
        FEE_delays = [int(delay) for delay in self.args.MWA_FEE_delays.strip("[]").split(',')]
        
        # ##Check the extra arguments have resulting in correct outputs
        c_delays =  [int(float(delay)) for delay in self.data['FEE_ideal_delays'].split(',')]
        
        self.assertEqual(FEE_delays, c_delays)

    def test_write_MWA_FEE_beam(self):
        """Test that the MWA FEE primary beam options work correctly"""

        def add_extra_MWAFEE_args():
            ##Add in primary_beam as Gaussian to inputs
            self.args.primary_beam = 'MWA_FEE'
            self.args.hdf5_beam_path = "This_is_the_way"
            self.args.MWA_FEE_delays = "[0,2,4,6,0,2,4,6,0,2,4,6,0,2,4,6]"
            
        def check_MWAFEE_outputs():
            self.assertEqual(2, int(self.data['beamtype']))
            self.assertEqual(self.args.hdf5_beam_path, self.data['hdf5_beam_path'])
            self.check_mwa_beam_delays()
            
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
            self.assertEqual(4, int(self.data['beamtype']))
            self.assertEqual(self.args.hdf5_beam_path, self.data['hdf5_beam_path'])
            self.check_mwa_beam_delays()
            
        
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
            self.assertEqual(5, int(self.data['beamtype']))
            self.check_mwa_beam_delays()
            
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
            self.assertEqual(3, int(self.data['beamtype']))

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
        self.assertTrue(int(self.data['do_autos']))
        
        self.args.do_autos = False
        woden_settings = self.call_fill_woden_settings_python()
        self.check_basic_inputs(woden_settings)
        self.assertFalse(int(self.data['do_autos']))


        
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
        self.assertTrue(int(self.data['do_precession']))
        
        self.args.no_precession = True
        woden_settings = self.call_fill_woden_settings_python()
        self.check_basic_inputs(woden_settings)
        self.assertFalse(int(self.data['do_precession']))
        
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
        self.assertFalse(int(self.data['use_dipamps']))
        
        
        ##Now ask for some dipole amps, provide some (as would happen using
        # run_setup.check_args) and check things are propagated correctly
        self.args.use_MWA_dipamps = True
        self.args.dipamps = np.arange(16)
        ##This runs `fill_woden_settings_python`
        woden_settings = self.call_fill_woden_settings_python()
        
        self.check_basic_inputs(woden_settings)
        ##check it we have a True flag
        self.assertTrue(int(self.data['use_dipamps']))
        
        mwa_dipole_amps = np.array([float(amp) for amp in self.data['mwa_dipole_amps'].split(',')])
        npt.assert_allclose(mwa_dipole_amps, np.arange(16), atol=1e-8)
        
        
    def test_off_cardinal_dipoles(self):
        """Test that the use_dipamp flag gets propagated correctly, and if
        true, the dipole amplitude array is correctly passed through"""

        ##First up, check that the off_cardinal_dipoles flag off by default
        self.make_basic_inputs('double')
        woden_settings = self.call_fill_woden_settings_python()
        self.check_basic_inputs(woden_settings)
        self.assertFalse(int(self.data['off_cardinal_dipoles']))
        
        ##Next, turn it on via the command line
        self.make_basic_inputs('double')
        self.args.off_cardinal_dipoles = True
        woden_settings = self.call_fill_woden_settings_python()
        self.check_basic_inputs(woden_settings)
        self.assertTrue(int(self.data['off_cardinal_dipoles']))
        
        ##Next, check it gets switched on given the right primary beam
        self.make_basic_inputs('double')
        self.args.primary_beam = 'everybeam_LOFAR'
        woden_settings = self.call_fill_woden_settings_python()
        self.check_basic_inputs(woden_settings)
        self.assertTrue(int(self.data['off_cardinal_dipoles']))
        
    def test_cpu_mode(self):
        """Test that the `cpu_mode` flag gets propagated correctly and sets
        the `do_gpu` flag correctly"""

        ##First up, check that the off_cardinal_dipoles flag off by default
        self.make_basic_inputs('double')
        woden_settings = self.call_fill_woden_settings_python()
        self.check_basic_inputs(woden_settings)
        self.assertTrue(int(self.data['do_gpu']))
        
        ##Next, turn it on via the command line
        self.make_basic_inputs('double')
        self.args.cpu_mode = True
        woden_settings = self.call_fill_woden_settings_python()
        self.check_basic_inputs(woden_settings)
        self.assertFalse(int(self.data['do_gpu']))
        
    def test_no_beam_norm(self):
        """Test that the `no_beam_normalistion` option gets propagated
        and sets the `normalise_primary_beam` flag correctly"""

        ##First up, check that the off_cardinal_dipoles flag off by default
        self.make_basic_inputs('double')
        woden_settings = self.call_fill_woden_settings_python()
        self.check_basic_inputs(woden_settings)
        self.assertTrue(int(self.data['normalise_primary_beam']))
        
        ##Next, turn it on via the command line
        self.make_basic_inputs('double')
        self.args.no_beam_normalisation = True
        woden_settings = self.call_fill_woden_settings_python()
        self.check_basic_inputs(woden_settings)
        self.assertFalse(int(self.data['normalise_primary_beam']))
        
##Run the test
if __name__ == '__main__':
   unittest.main()