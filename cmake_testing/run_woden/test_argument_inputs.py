from sys import path
import os
import unittest
import numpy as np

##This is where our code lives
code_dir = os.environ['CMAKE_CURRENT_SOURCE_DIR']

##Append the location of run_woden.py to the sys.path to import it
path.append('{:s}/../../src'.format(code_dir))

##Code we are testing
import run_woden as rw

##some expected values
east = [-5.85675e+02, -5.75557e+02, -4.89943e+02, -5.18046e+02, -1.02607e+03,
        -1.07960e+03, -6.85027e+02, -6.53467e+02,  3.70650e+01, -5.88010e+01,
         4.42990e+02,  6.53774e+02,  3.26101e+02,  3.25688e+02,  2.57984e+02,
         2.32842e+02, -5.75533e+02, -4.80270e+02, -1.41898e+02, -4.12457e+02,
        -5.68540e+02, -6.33237e+02, -6.31403e+02, -6.74021e+02, -3.88152e+02,
        -2.84457e+02, -1.98194e+02, -1.95980e+02, -3.69608e+02, -3.64524e+02,
        -4.68817e+02, -4.52686e+02,  2.23817e+03,  2.32742e+03,  1.99447e+03,
         1.65437e+03,  3.02946e+03,  2.40567e+03,  2.67033e+03,  2.69858e+03,
         1.36097e+03,  8.79780e+02,  8.75570e+02,  3.63950e+02,  1.22470e+03,
         1.63127e+03,  1.37137e+03,  1.68351e+03,  2.44155e+03,  2.40004e+03,
         2.05493e+03,  1.90759e+03,  2.87228e+03,  3.10564e+03,  2.77082e+03,
         2.34629e+03, -1.75090e+02, -4.55360e+02, -3.30610e+02, -8.29330e+02,
         3.27910e+02,  5.20820e+02,  7.76550e+02, -5.25500e+01,  2.27878e+02,
         2.38310e+02,  2.88066e+02,  7.16393e+02,  1.01500e+00,  3.52240e+01,
        -5.88740e+01,  6.14430e+01, -2.56952e+02, -1.73717e+02, -2.27151e+02,
        -3.33758e+02, -8.89471e+02, -3.00696e+02, -2.49647e+02, -1.95018e+02,
         3.01193e+02,  4.98629e+02,  1.20113e+03,  1.10094e+03,  5.23727e+02,
         5.06793e+02,  5.44324e+02,  3.70648e+02,  5.56122e+02,  1.31361e+03,
         1.17729e+03,  8.39943e+02,  3.81069e+02,  5.28598e+02,  4.88857e+02,
         5.71541e+02,  1.67727e+03,  1.17792e+03,  1.47178e+03,  8.03010e+02,
         1.99131e+03,  2.28548e+03,  1.88127e+03,  1.55938e+03,  8.05810e+02,
         3.38490e+02,  5.11310e+02,  1.65760e+02,  1.66429e+03,  1.33634e+03,
         9.31990e+02,  1.10572e+03, -1.47317e+03, -1.60021e+03, -1.77014e+03,
        -1.99981e+03, -8.89290e+02, -1.10097e+03, -1.17292e+03, -1.36727e+03,
         9.18090e+01, -2.17100e+01,  1.77467e+02, -5.55600e+01, -8.75160e+01,
        -1.18037e+02, -6.75820e+01,  1.53700e+02]
north = [ -101.53,    415.008,  604.568,  632.503,  574.527, 87.018,   212.168,
           142.255,  -405.05,  -505.332, -835.306, -1037.75, -306.307, -318.508,
          -366.453,  -357.891, -863.782, -817.324, -420.261, -409.06,  -257.66,
          -169.362,  -245.107, -638.871, -44.022,  -196.181, 95.15,    189.795,
           266.484,   183.098,  185.187, -15.933,  408.81,   1026.38,  629.74,
           845.05,    1095.29,    23.61,  572.92,  931.59,   2458.85,  1687.73,
           2181.3,    1536.78,  1387.35,  1497.8,  1863.8,   2153.35,  1785.7,
           2214.24,   1775.43,  1247.36,  1455.11, 1640.33,  1905.8,   1367.78,
           2644.3,    2039.45,  2829.62,  2563.46, 1962.28,  2352.17,  2535.7,
           2142.93,   819.724,  768.413,  618.265, 1405.26,  1412.36,  950.242,
           810.53,    815.255,  892.005,  781.638, 898.129,  1660.73,  1088.14,
           704.326,   731.709,  762.893, -48.952,  -264.594, -443.25,  151.76,
           75.494,   -3.216,   -91.152,  -38.415,  244.485,  553.372,  858.632,
           990.227,   652.864,  539.194,  424.98,  359.21,   -478.33,  -691.74,
           179.53,   -812.85,  -940.61,  -251.75,  -106.66,  -932.27,  -1408.62,
           -2065.83, -1486.98, -991.57,  -1355.78, -1867.29, -1975.87, -1139.39,
           -194.32,  -64.26,    53.59,    206.85,  -779.73,  -698.66,  -569.74,
           -453.65,  -80.565,   77.005,  -43.377,  124.801,  125.197,  -39.873,
           -22.361,  -259.929]
height = [375.212, 373.374, 372.907, 372.604, 374.039, 377.009, 374.236, 374.462, 374.034,
          372.913, 372.251, 371.071, 375.128, 375.029, 375.221, 375.262, 375.522, 375.72,
          373.503, 373.969, 375.246, 375.616, 375.544, 375.34,  375.372, 374.901, 377.024,
          377.161, 375.24,  375.739, 374.752, 375.111, 371.5,   374.52,  373.7,   372.22,
          373.69,  369.35,  371.79,  373.23,  380.61,  376.53,  377.05,  374.63,  376.46,
          376.74,  378.97,  380.99,  377.35,  380.84,  378.58,  375.33,  376.45,  376.8,
          378.43,  375.66,  373.96,  370.95,  373.37,  369.52,  374.38,  376.35,  377.7,
          373.19,  375.748, 375.805, 375.134, 374.602, 373.682, 376.275, 375.17,  376.001,
          373.492, 374.217, 373.694, 371.068, 370.381, 373.514, 373.797, 374.093, 374.618,
          373.22,  368.328, 369.065, 372.613, 372.803, 372.812, 374.34,  372.284, 370.354,
          372.122, 372.696, 374.251, 372.789, 372.967, 372.509, 366.65,  368.31,  369.24,
          371.04,  365.15,  367.95,  368.26,  367.41,  367.88,  365.99,  368.75,  370.3,
          365.88,  366.17,  365.39,  370.54,  376.22,  375.54,  374.44,  372.92,  374.41,
          374.12,  374.77,  374.91,  375.76,  376.351, 375.005, 376.803, 377.031, 376.054,
          376.175, 376.476]


##Vehicle for running tests
class Test(unittest.TestCase):
    """Test whether the args collected by the argument parser are read in
    correctly, and that sanity checks on certain combinations work such that
    we don't feed WODEN arguments that won't work"""

    def run_parser_on_inputs(self):
        """Call `rw.get_parser` and run the returned parser using the inputs
        in `self.inputs`. Return the recovered arguments"""
        parser = rw.get_parser()
        args =  parser.parse_args(self.inputs)
        return args

    def assert_parser_errors(self):
        """Assert that the parser returned by `rw.get_parser` errors when
        run with the current set of self.inputs"""
        with self.assertRaises(SystemExit) as cm:
            ##call the argparser with the gives=n inputs
            args = self.run_parser_on_inputs()

    def assert_check_args_errors(self):
        """Assert that `rw.check_args` errors when run on the parser returned
        by `rw.get_parser` is run with the current arguments in `self.inputs`"""
        ##call the argparser with the given inputs
        args = self.run_parser_on_inputs()
        ##Assert the code raisers a sys exit
        with self.assertRaises(SystemExit) as cm:
            rw.check_args(args)

        return args

    def run_parser_and_check_args(self):
        """Runs the parser on the current set of args in self.inputs, then runs
        the output through `rw.check_args`"""
        args = self.run_parser_on_inputs()
        args = rw.check_args(args)
        return args

    def make_required_args(self):
        """These are arguments that if missing, the parser itself should fail"""
        self.inputs = ['--ra0=0.0', '--dec0=-26.7', '--cat_filename=srclist.txt']

    def make_minimum_required_args_without_metafits(self):
        """When not passing a metafits file, these are the minimum set of
        arguments needed to pass onto WODEN. Other arguments either have
        defaults or are optional"""
        self.make_required_args()

        self.inputs.append('--num_time_steps=16')
        self.inputs.append('--time_res=8.0')
        self.inputs.append('--freq_res=40e+3')
        self.inputs.append('--array_layout={:s}/example_array_layout.txt'.format(code_dir))
        self.inputs.append('--lowest_channel_freq=160e+6')
        self.inputs.append('--date="2019-06-12T13:04:12"')

    def test_parser_fails(self):
        """Tests the `rw.get_parser` function, which creates an
        `argparse.ArgumentParser` object to collect the arguments from the
        command. There are three arguments which are always required -
        check things fail if they are missing"""

        self.inputs = []
        self.assert_parser_errors()

        ##Add one of the required args
        self.inputs.append('--ra0=0.0')
        self.assert_parser_errors()

        ##Add second required args
        self.inputs.append('--dec0=-26.7')
        self.assert_parser_errors()

        ##Add final required args
        self.inputs.append('--cat_filename=srclist.txt')
        ##This should run now
        args = self.run_parser_on_inputs()

    def test_missing_args_without_metafits_fails(self):
        """When not providing a metafits file, a number of arguments must
        be given to complete the observational settings. Can't make them
        required by the argparser as only needed if metafits not given,
        so we use `rw.check_args` to make sure everything that is needed
        is present.

        Here, we make a list of every neccessary input argument in a loop,
        deleting a different one each time to make sure if one is missing
        we get an error every time"""

        ##First 3 arguments in self.inputs generated by
        ##self.make_minimum_required_args_without_metafits are required by
        ##the parser, so start loop at 3

        self.make_minimum_required_args_without_metafits()
        num_req_args = len(self.inputs)

        ##Loop through the needed args, delete them from the list of inputs
        ##and make sure it errors for each one
        for arg_ind in range(4, num_req_args):
            ##regenerate input args to replace anything that was deleted
            self.make_minimum_required_args_without_metafits()
            del[self.inputs[arg_ind]]

            ##Check that it errors
            self.assert_check_args_errors()

    def test_metafits_read_fails(self):
        """If the path to the metafits file is not real, `rw.check_args`
        should fail"""
        self.make_required_args()
        ##Add a bad path to the metafits option
        self.inputs.append("--metafits_filename=this_is_not_a_path")
        ##Check it fails with the bad path
        self.assert_check_args_errors()

    def test_read_metafits_succeeds(self):
        self.make_required_args()
        ##Add a path to a metafits to the options
        self.inputs.append("--metafits_filename={:s}/1202815152_metafits_ppds.fits".format(code_dir))

        ##Run the parser and get the args
        args = self.run_parser_on_inputs()
        ##Check the args make sense and read in options from the metafits
        args = rw.check_args(args)

        ##Check the options set by the metafits are correct
        self.assertEqual(169595000.0, args.lowest_channel_freq)
        self.assertEqual(240, args.num_time_steps)
        self.assertEqual(10000.0, args.freq_res)
        self.assertEqual(0.5, args.time_res)
        self.assertEqual("2018-02-16T11:18:54", args.date)
        self.assertEqual("from_the_metafits", args.array_layout)
        self.assertEqual(128, args.num_freq_channels)
        self.assertEqual(range(1, 25), args.band_nums)

        self.assertEqual(-26.7033194444 ,args.latitude)
        self.assertEqual(1280000.0 ,args.coarse_band_width)

        self.assertEqual(54.75681779724241, args.gauss_ra_point)
        self.assertEqual(-39.48422590285089, args.gauss_dec_point)

        self.assertTrue(np.allclose(np.array(east), args.east, atol=1e-10))
        self.assertTrue(np.allclose(np.array(north), args.north, atol=1e-10))
        self.assertTrue(np.allclose(np.array(height), args.height, atol=1e-10))

    def test_EDA2_args_work(self):
        """Check that with a minimal set of arguments, the EDA2 primary
        beam is selected by `rw.check_args` correctly"""
        self.make_minimum_required_args_without_metafits()
        self.inputs.append('--primary_beam=EDA2')

        args = self.run_parser_and_check_args()

        self.assertEqual('EDA2', args.primary_beam)

    def test_GaussBeam_args_work(self):
        """Check that the Gaussian primary beam is handled `ra.check_args`
        correctly, and that optional arguments are added correctly"""
        self.make_minimum_required_args_without_metafits()
        self.inputs.append('--primary_beam=Gaussian')

        args = self.run_parser_and_check_args()

        self.assertEqual('Gaussian', args.primary_beam)

        ##Should be pointed at phase centre with False FWHM and ref freq
        ##as none of these have been specified
        self.assertEqual(args.ra0, args.gauss_ra_point)
        self.assertEqual(args.dec0, args.gauss_dec_point)
        self.assertFalse(args.gauss_beam_FWHM)
        self.assertFalse(args.gauss_beam_ref_freq)

        ##Keep adding in optional arguments and checking they end up manifesting
        self.inputs.append('--gauss_ra_point=234.0')
        args = self.run_parser_and_check_args()
        self.assertEqual(234.0, args.gauss_ra_point)
        self.assertEqual(args.dec0, args.gauss_dec_point)
        self.assertFalse(args.gauss_beam_FWHM)
        self.assertFalse(args.gauss_beam_ref_freq)

        self.inputs.append('--gauss_dec_point=19.0')
        args = self.run_parser_and_check_args()
        self.assertEqual(234.0, args.gauss_ra_point)
        self.assertEqual(19.0, args.gauss_dec_point)
        self.assertFalse(args.gauss_beam_FWHM)
        self.assertFalse(args.gauss_beam_ref_freq)

        self.inputs.append('--gauss_beam_FWHM=64.0')
        args = self.run_parser_and_check_args()
        self.assertEqual(234.0, args.gauss_ra_point)
        self.assertEqual(19.0, args.gauss_dec_point)
        self.assertEqual(64.0, args.gauss_beam_FWHM)
        self.assertFalse(args.gauss_beam_ref_freq)

        self.inputs.append('--gauss_beam_ref_freq=291e+6')
        args = self.run_parser_and_check_args()
        self.assertEqual(234.0, args.gauss_ra_point)
        self.assertEqual(19.0, args.gauss_dec_point)
        self.assertEqual(64.0, args.gauss_beam_FWHM)
        self.assertEqual(291e+6, args.gauss_beam_ref_freq)

    def _check_MWA_FEE_generic(self, beam_name, beam_env_var):
        """Run the tests common to `test_MWAFEEBeam_args_work` and
        `test_MWAFEEBeamInterp_args_work`. The function should error out
        if certain paths to the hdf5 file that holds the spherical harmonic
        information is missing, and if the delays have been specified
        incorrectly"""

        ##We want to test that the argument fails if the environment key
        ##MWA_FEE_HDF5 is not set. Here we check if it is set and delete it
        ##if so
        try:
            del os.environ[beam_env_var]
        except KeyError:
            pass

        ##Trying to run an MWA_FEE simulation with no metafits, no 'MWA_FEE_HDF5'
        ##or --hdf5_beam_path should fail
        self.make_minimum_required_args_without_metafits()
        self.inputs.append(f'--primary_beam={beam_name}')

        ##Check the primary_beam has been selected and that `rw.check_args` fails
        args = self.assert_check_args_errors()
        self.assertEqual(beam_name, args.primary_beam)

        ##Set the path to the hdf5 file. Just point it at a text file.
        ##TODO - have `check_args` actually test reading the hdf5 file. At the
        ##mo just tests that the file exists
        ##This should still fail because we have no delays set
        self.inputs.append('--hdf5_beam_path={:s}/example_array_layout.txt'.format(code_dir))
        self.assert_check_args_errors()

        ##now set the delays to something that should produce a failure
        self.inputs.append('--MWA_FEE_delays=oijasdoiasjd')
        self.assert_check_args_errors()

        ##Set the delays to something that should work
        self.inputs.append('--MWA_FEE_delays=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]')
        self.run_parser_and_check_args()

        ##Reset the arguments try to run using a metafits file, but still
        ##have no path to the hdf5 file. Should fail
        self.make_minimum_required_args_without_metafits()
        self.inputs.append(f'--primary_beam={beam_name}')
        self.inputs.append("--metafits_filename={:s}/1202815152_metafits_ppds.fits".format(code_dir))
        args = self.assert_check_args_errors()

        ##This time, set the environment variable to the hdf5 file. Should
        ##pass (just use a file we know exists somewhere)
        os.environ[beam_env_var] = '{:s}/example_array_layout.txt'.format(code_dir)
        args = self.run_parser_and_check_args()
        ##Assert the delays were read in correctly
        self.assertEqual("[6, 4, 2, 0, 8, 6, 4, 2, 10, 8, 6, 4, 12, 10, 8, 6]", args.MWA_FEE_delays)

        ##Setting the delays manually should override the delays in the
        ##metafits
        self.inputs.append('--MWA_FEE_delays=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]')
        args = self.run_parser_and_check_args()
        ##Check the command line delays are there instead of metafits
        self.assertEqual("[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]", args.MWA_FEE_delays)

    def test_MWAFEEBeam_args_work(self):
        """Check that the MWA FEE primary beam is handled `ra.check_args`
        correctly. The function should error out if certain paths to the
        hdf5 file that holds the spherical harmonic information is missing,
        and if the delays have been specified incorrectly"""

        self._check_MWA_FEE_generic('MWA_FEE', 'MWA_FEE_HDF5')


    def test_MWAFEEBeamInterp_args_work(self):
        """Check that the interpolated MWA FEE primary beam is handled `ra.check_args`
        correctly. The function should error out if certain paths to the
        hdf5 file that holds the spherical harmonic information is missing,
        and if the delays have been specified incorrectly"""

        self._check_MWA_FEE_generic('MWA_FEE_interp', 'MWA_FEE_HDF5_INTERP')


    def test_MWAAnalyBeam_args_work(self):
        """Check `ra.check_args` works correctly for the analytic MWA beam.
        Should fail if the delays have been specified incorrectly"""

        ##Trying to run an MWA_FEE simulation with no metafits, no 'MWA_FEE_HDF5'
        ##or --hdf5_beam_path should fail
        self.make_minimum_required_args_without_metafits()
        self.inputs.append(f'--primary_beam=MWA_analy')

        ##Check the primary_beam has been selected and that `rw.check_args` fails
        args = self.assert_check_args_errors()
        self.assertEqual('MWA_analy', args.primary_beam)

        ##now set the delays to something that should produce a failure
        self.inputs.append('--MWA_FEE_delays=oijasdoiasjd')
        self.assert_check_args_errors()

        ##Set the delays to something that should work
        self.inputs.append('--MWA_FEE_delays=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]')
        self.run_parser_and_check_args()

        ##Reset the arguments try to run using delays from a metafits file.
        ##Should pass
        self.make_minimum_required_args_without_metafits()
        self.inputs.append(f'--primary_beam=MWA_analy')
        self.inputs.append("--metafits_filename={:s}/1202815152_metafits_ppds.fits".format(code_dir))
        self.run_parser_and_check_args()

        ##Setting the delays manually should override the delays in the
        ##metafits
        self.inputs.append('--MWA_FEE_delays=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]')
        args = self.run_parser_and_check_args()
        ##Check the command line delays are there instead of metafits
        self.assertEqual("[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]", args.MWA_FEE_delays)


##Run the test
if __name__ == '__main__':
    unittest.main()
