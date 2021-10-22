from sys import path
import os
import unittest
import json

##Do some disgusting path finding exercise, there must be a better way
##to do this
fileloc = os.path.realpath(__file__)
path.append('{:s}/../../src'.format(('/').join(fileloc.split('/')[:-1])))

##Code we are testing
import run_woden as rw

class PretendArgs():
    """Function `rw.write_json` takes an arg class out of argparse as an
    input. Write a dummy one here to feed in and test"""
    def make_basic_args(self):
        self.ra0 = 0.0
        self.dec0 = -26.7
        self.cat_filename = 'srclist_of_dreams.txt'
        self.chunking_size = 1e+10
        self.coarse_band_width = 1.28e+6
        self.sky_crop_components = False
        self.primary_beam = "none"

        self.num_freq_channels = 16
        self.num_time_steps = 14
        self.time_res = 2.0
        self.freq_res = 40e+3

        self.array_layout_name = 'WODEN_array_layout.txt'
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

    def __init__(self):
        self.make_basic_args()

##Vehicle for running tests
class Test(unittest.TestCase):
    """Test the `rw.write_json` function, which writes an input file to feed
    into the WODEN executable"""
    def make_basic_inputs(self, json_name):
        """Make some basis input arguments for `rw.write_json`"""
        self.jd_date = 2458647.044583333
        self.json_name = json_name
        self.lst = 10.0
        self.args = PretendArgs()

    def check_basic_inputs(self):
        """Attempts to read in a .json and checks the information that should
        have been written with the minimum set of inputs"""

        self.assertTrue(os.path.isfile(self.json_name))

        with open(self.json_name) as json_file:
            data = json.load(json_file)

        self.data = data

        self.assertEqual(self.args.ra0 ,data['ra0'])
        self.assertEqual(self.args.dec0 ,data['dec0'])
        self.assertEqual(self.args.num_freq_channels ,data['num_freqs'])
        self.assertEqual(self.args.num_time_steps ,data['num_time_steps'])
        self.assertEqual(self.args.cat_filename ,data['cat_filename'])
        self.assertEqual(self.args.time_res ,data['time_res'])
        self.assertEqual(self.args.freq_res ,data['frequency_resolution'])
        self.assertEqual(self.args.chunking_size ,data['chunking_size'])
        self.assertEqual(self.jd_date ,data['jd_date'])
        self.assertEqual(self.lst ,data['LST'])
        self.assertEqual(self.args.array_layout_name ,data['array_layout'])
        self.assertEqual(self.args.lowest_channel_freq ,data['lowest_channel_freq'])
        self.assertEqual(self.args.latitude ,data['latitude'])
        self.assertEqual(self.args.coarse_band_width ,data['coarse_band_width'])
        self.assertEqual(self.args.band_nums ,data['band_nums'])

    def call_write_json(self):
        """Calls the function under test"""
        rw.write_json(json_name=self.json_name, lst=self.lst,
                   jd_date=self.jd_date, args=self.args)

    def test_write_minimum_json(self):
        """Test that a .json file is written correctly with the minimum set
        of input arguments. Test by calling function, reading in resultant
        .json file and checking the outputs match"""

        ##Some input test params
        self.make_basic_inputs("test_write_minimum_json.json")

        ##Code to be tested
        self.call_write_json()

        ##Check the basic outputs
        self.check_basic_inputs()

    def test_write_gaussian_beam(self):
        """Test that the Gaussian primary beam options work correctly"""

        ##Some input test params
        self.make_basic_inputs("test_write_gaussian_beam.json")
        ##Add in primary_beam as Gaussian to inputs
        self.args.primary_beam = 'Gaussian'
        self.args.gauss_beam_FWHM = 20.0
        self.args.gauss_beam_ref_freq = 180e+6
        self.args.gauss_ra_point = 60.0
        self.args.gauss_dec_point = -10.0

        ##Code to be tested
        self.call_write_json()

        ##Check the basic outputs - need to call this to read in the json file
        ##and update self.data
        self.check_basic_inputs()

        ##Check the extra arguments have resulting in correct outputs
        ##Optional arguments should now exists
        self.assertTrue(self.data['use_gaussian_beam'])
        self.assertEqual(self.args.gauss_ra_point, self.data['gauss_ra_point'])
        self.assertEqual(self.args.gauss_dec_point, self.data['gauss_dec_point'])
        self.assertEqual(self.args.gauss_beam_FWHM, self.data['gauss_beam_FWHM'])
        self.assertEqual(self.args.gauss_beam_ref_freq, self.data['gauss_beam_ref_freq'])

    def test_write_MWA_FEE_beam(self):
        """Test that the MWA FEE primary beam options work correctly"""

        ##Some input test params
        self.make_basic_inputs("test_write_mwafee_beam.json")
        ##Add in primary_beam as Gaussian to inputs
        self.args.primary_beam = 'MWA_FEE'
        self.args.hdf5_beam_path = "This_is_the_way"
        self.args.MWA_FEE_delays = "[0,2,4,6,0,2,4,6,0,2,4,6,0,2,4,6]"
        ##turn string of list into an actual list
        FEE_delays = [int(delay) for delay in self.args.MWA_FEE_delays.strip("[]").split(',')]

        ##Code to be tested
        self.call_write_json()

        ##Check the basic outputs
        self.check_basic_inputs()

        ##Check the extra arguments have resulting in correct outputs
        self.assertTrue(self.data['use_FEE_beam'])
        self.assertEqual(self.args.hdf5_beam_path, self.data['hdf5_beam_path'])
        self.assertEqual(FEE_delays, self.data['FEE_delays'])

    def test_write_EDA2_beam(self):
        """Test that the EDA2 primary beam options work correctly"""

        ##Some input test params
        self.make_basic_inputs("test_write_mwafee_beam.json")
        ##Add in primary_beam as Gaussian to inputs
        self.args.primary_beam = 'EDA2'

        ##Code to be tested
        self.call_write_json()

        ##Check the basic outputs
        self.check_basic_inputs()

        ##Check the extra arguments have resulting in correct outputs
        self.assertTrue(self.data['use_EDA2_beam'])

    def test_write_no_precession(self):
        """Test that no_precession is written when selected in args"""

        ##Some input test params
        self.make_basic_inputs("test_write_np_precess.json")
        ##Add in primary_beam as Gaussian to inputs
        self.args.no_precession = True

        ##Code to be tested
        self.call_write_json()

        ##Check the basic outputs
        self.check_basic_inputs()

        ##Check the extra arguments have resulting in correct outputs
        self.assertTrue(self.data['no_precession'])

##Run the test
if __name__ == '__main__':
   unittest.main()
