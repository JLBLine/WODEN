from sys import path
import os
import unittest
import numpy as np
import numpy.testing as npt

code_dir = os.path.realpath(__file__)
code_dir = ('/').join(code_dir.split('/')[:-1])

# ##Code we are testing
from wodenpy.wodenpy_setup import run_setup

##some expected values
east = [-5.55600e+01, 1.77467e+02, -2.17100e+01, 9.18090e+01, 1.53700e+02,
        -6.75820e+01, -1.18037e+02, -8.75160e+01, -1.95980e+02, -1.98194e+02,
        -2.84457e+02, -3.88152e+02, -4.52686e+02, -4.68817e+02, -3.64524e+02,
        -3.69608e+02, -5.18046e+02, -4.89943e+02, -5.75557e+02, -5.85675e+02,
        -6.53467e+02, -6.85027e+02, -1.07960e+03, -1.02607e+03, -3.33758e+02,
        -2.27151e+02, -1.73717e+02, -2.56952e+02, -1.95018e+02, -2.49647e+02,
        -3.00696e+02, -8.89471e+02, 7.16393e+02, 2.88066e+02, 2.38310e+02,
        2.27878e+02, 6.14430e+01, -5.88740e+01, 3.52240e+01, 1.01500e+00,
        8.39943e+02, 1.17729e+03, 1.31361e+03, 5.56122e+02, 5.71541e+02,
        4.88857e+02, 5.28598e+02, 3.81069e+02, 1.10094e+03, 1.20113e+03,
        4.98629e+02, 3.01193e+02, 3.70648e+02, 5.44324e+02, 5.06793e+02,
        5.23727e+02, 6.53774e+02, 4.42990e+02, -5.88010e+01, 3.70650e+01,
        2.32842e+02, 2.57984e+02, 3.25688e+02, 3.26101e+02, -4.12457e+02,
        -1.41898e+02, -4.80270e+02, -5.75533e+02, -6.74021e+02, -6.31403e+02,
        -6.33237e+02, -5.68540e+02, -1.99981e+03, -1.77014e+03, -1.60021e+03,
        -1.47317e+03, -1.36727e+03, -1.17292e+03, -1.10097e+03, -8.89290e+02,
        -8.29330e+02, -3.30610e+02, -4.55360e+02, -1.75090e+02, -5.25500e+01,
        7.76550e+02, 5.20820e+02, 3.27910e+02, 3.63950e+02, 8.75570e+02,
        8.79780e+02, 1.36097e+03, 1.68351e+03, 1.37137e+03, 1.63127e+03,
        1.22470e+03, 1.90759e+03, 2.05493e+03, 2.40004e+03, 2.44155e+03,
        2.34629e+03, 2.77082e+03, 3.10564e+03, 2.87228e+03, 1.65437e+03,
        1.99447e+03, 2.32742e+03, 2.23817e+03, 2.69858e+03, 2.67033e+03,
        2.40567e+03, 3.02946e+03, 8.03010e+02, 1.47178e+03, 1.17792e+03,
        1.67727e+03, 1.55938e+03, 1.88127e+03, 2.28548e+03, 1.99131e+03,
        1.65760e+02, 5.11310e+02, 3.38490e+02, 8.05810e+02, 1.10572e+03,
        9.31990e+02, 1.33634e+03, 1.66429e+03]

north = [ 124.801, -43.377, 77.005, -80.565, -259.929, -22.361, -39.873,
        125.197, 189.795, 95.15, -196.181, -44.022, -15.933, 185.187,
        183.098, 266.484, 632.503, 604.568, 415.008, -101.53, 142.255,
        212.168, 87.018, 574.527, 1660.73, 898.129, 781.638, 892.005,
        762.893, 731.709, 704.326, 1088.14, 1405.26, 618.265, 768.413,
        819.724, 815.255, 810.53, 950.242, 1412.36, 990.227, 858.632,
        553.372, 244.485, 359.21, 424.98, 539.194, 652.864, 151.76,
        -443.25, -264.594, -48.952, -38.415, -91.152, -3.216, 75.494,
        -1037.75, -835.306, -505.332, -405.05, -357.891, -366.453, -318.508,
        -306.307, -409.06, -420.261, -817.324, -863.782, -638.871, -245.107,
        -169.362, -257.66, 206.85, 53.59, -64.26, -194.32, -453.65,
        -569.74, -698.66, -779.73, 2563.46, 2829.62, 2039.45, 2644.3,
        2142.93, 2535.7, 2352.17, 1962.28, 1536.78, 2181.3, 1687.73,
        2458.85, 2153.35, 1863.8, 1497.8, 1387.35, 1247.36, 1775.43,
        2214.24, 1785.7, 1367.78, 1905.8, 1640.33, 1455.11, 845.05,
        629.74, 1026.38, 408.81, 931.59, 572.92, 23.61, 1095.29,
        -812.85, 179.53, -691.74, -478.33, -932.27, -106.66, -251.75,
        -940.61, -991.57, -1486.98, -2065.83, -1408.62, -1139.39, -1975.87,
        -1867.29, -1355.78]

height = [376.803, 375.005, 376.351, 375.76, 376.476,
          376.175, 376.054, 377.031, 377.161,
          377.024, 374.901, 375.372, 375.111, 374.752, 375.739, 375.24, 372.604, 372.907,
          373.374, 375.212, 374.462, 374.236, 377.009, 374.039, 371.068, 373.694, 374.217,
          373.492, 374.093, 373.797, 373.514, 370.381, 374.602, 375.134, 375.805, 375.748,
          376.001, 375.17, 376.275, 373.682, 372.696, 372.122, 370.354, 372.284, 372.509,
          372.967, 372.789, 374.251, 369.065, 368.328, 373.22, 374.618, 374.34, 372.812,
          372.803, 372.613, 371.071, 372.251, 372.913, 374.034, 375.262, 375.221, 375.029,
          375.128, 373.969, 373.503, 375.72, 375.522, 375.34, 375.544, 375.616, 375.246,
          372.92, 374.44, 375.54, 376.22, 374.91, 374.77, 374.12, 374.41, 369.52,
          373.37, 370.95, 373.96, 373.19, 377.7, 376.35, 374.38, 374.63, 377.05,
          376.53, 380.61, 380.99, 378.97, 376.74, 376.46, 375.33, 378.58, 380.84,
          377.35, 375.66, 378.43, 376.8, 376.45, 372.22, 373.7, 374.52, 371.5,
          373.23, 371.79, 369.35, 373.69, 371.04, 369.24, 368.31, 366.65, 367.41,
          368.26, 367.95, 365.15, 370.3, 368.75, 365.99, 367.88, 370.54, 365.39,
          366.17, 365.88]

names = ['Tile051', 'Tile052', 'Tile053', 'Tile054',
         'Tile055', 'Tile056', 'Tile057',
'Tile058', 'Tile071', 'Tile072', 'Tile073', 'Tile074', 'Tile075', 'Tile076', 
'Tile077', 'Tile078', 'Tile101', 'Tile102', 'Tile103', 'Tile104', 'Tile105', 
'Tile106', 'Tile107', 'Tile108', 'Tile111', 'Tile112', 'Tile113', 'Tile114', 
'Tile115', 'Tile116', 'Tile117', 'Tile118', 'Tile121', 'Tile122', 'Tile123', 
'Tile124', 'Tile125', 'Tile126', 'Tile127', 'Tile128', 'Tile131', 'Tile132', 
'Tile133', 'Tile134', 'Tile135', 'Tile136', 'Tile137', 'Tile138', 'Tile141', 
'Tile142', 'Tile143', 'Tile144', 'Tile145', 'Tile146', 'Tile147', 'Tile148', 
'Tile151', 'Tile152', 'Tile153', 'Tile154', 'Tile155', 'Tile156', 'Tile157', 
'Tile158', 'Tile161', 'Tile162', 'Tile163', 'Tile164', 'Tile165', 'Tile166', 
'Tile167', 'Tile168', 'LBA1', 'LBA2', 'LBA3', 'LBA4', 'LBA5', 'LBA6', 'LBA7', 
'LBA8', 'LBB1', 'LBB2', 'LBB3', 'LBB4', 'LBB5', 'LBB6', 'LBB7', 'LBB8', 'LBC1', 
'LBC2', 'LBC3', 'LBC4', 'LBC5', 'LBC6', 'LBC7', 'LBC8', 'LBD1', 'LBD2', 'LBD3', 
'LBD4', 'LBD5', 'LBD6', 'LBD7', 'LBD8', 'LBE1', 'LBE2', 'LBE3', 'LBE4', 'LBE5', 
'LBE6', 'LBE7', 'LBE8', 'LBF1', 'LBF2', 'LBF3', 'LBF4', 'LBF5', 'LBF6', 'LBF7', 
'LBF8', 'LBG1', 'LBG2', 'LBG3', 'LBG4', 'LBG5', 'LBG6', 'LBG7', 'LBG8']


dipflags = np.array([  26,   59,   97,  104,  107,  118,  126,  197,  198,  205,  207,
        214,  218,  220,  221,  447,  534,  616, 1032, 1044, 1059, 1074,
       1092, 1099, 1101, 1122, 1148, 1163, 1179, 1215, 1217, 1252, 1350,
       1386, 1413, 1583, 1588, 1629, 1716, 1740, 1744, 1811, 2172, 2307,
       2327, 2398, 2418, 2496, 2569, 2629, 2652, 2657, 2667, 2703, 2719,
       2745, 2816, 2842, 2864, 2981, 3102, 3247, 3260, 3281, 3308, 3359,
       3368, 3404, 3452, 3480, 3561, 3617, 3634, 3679, 3700, 3756, 3791,
       3823, 3881, 3896, 3985, 4019, 4073])

dipamp_indexes = np.array([ 26, 56, 79, 100, 345, 768, 2780, 3678, 4000])
dipamps = np.array([0.90376163, 1., 0.85784072, 0.70395702, 1., 1.,
                    0.95163655, 0.9428432, 0.89301419])

dipamps_flagged = np.array([0.0, 1., 0.85784072, 0.70395702, 1., 1.,
                    0.95163655, 0.9428432, 0.89301419])

##Vehicle for running tests
class Test(unittest.TestCase):
    """Test whether the args collected by the argument parser are read in
    correctly, and that sanity checks on certain combinations work such that
    we don't feed WODEN arguments that won't work"""

    def run_parser_on_inputs(self):
        """Call `run_setup.get_parser` and run the returned parser using the inputs
        in `self.inputs`. Return the recovered arguments"""
        parser = run_setup.get_parser()
        args =  parser.parse_args(self.inputs)
        return args

    def assert_parser_errors(self):
        """Assert that the parser returned by `run_setup.get_parser` errors when
        run with the current set of self.inputs"""
        with self.assertRaises(SystemExit) as cm:
            ##call the argparser with the gives=n inputs
            args = self.run_parser_on_inputs()

    def assert_check_args_errors(self):
        """Assert that `run_setup.check_args` errors when run on the parser returned
        by `run_setup.get_parser` is run with the current arguments in `self.inputs`"""
        ##call the argparser with the given inputs
        args = self.run_parser_on_inputs()
        ##Assert the code raisers a sys exit
        with self.assertRaises(SystemExit) as cm:
            run_setup.check_args(args)

        return args

    def run_parser_and_check_args(self):
        """Runs the parser on the current set of args in self.inputs, then runs
        the output through `run_setup.check_args`"""
        args = self.run_parser_on_inputs()
        args = run_setup.check_args(args)
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
        """Tests the `run_setup.get_parser` function, which creates an
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
        so we use `run_setup.check_args` to make sure everything that is needed
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
        """If the path to the metafits file is not real, `run_setup.check_args`
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
        args = run_setup.check_args(args)

        ##Check the options set by the metafits are correct
        self.assertEqual(169595000.0, args.lowest_channel_freq)
        self.assertEqual(240, args.num_time_steps)
        self.assertEqual(10000.0, args.freq_res)
        self.assertEqual(0.5, args.time_res)
        self.assertEqual("2018-02-16T11:18:54", args.date)
        self.assertEqual("from_the_metafits", args.array_layout)
        self.assertEqual(128, args.num_freq_channels)
        self.assertEqual(range(1, 25), args.band_nums)

        self.assertEqual(-26.703319405555554 ,args.latitude)
        self.assertEqual(1280000.0 ,args.coarse_band_width)

        self.assertEqual(54.75681779724241, args.gauss_ra_point)
        self.assertEqual(-39.48422590285089, args.gauss_dec_point)

        ##Check east, north, height coords read correctly
        self.assertTrue(np.allclose(np.array(east), args.east, atol=1e-5))
        self.assertTrue(np.allclose(np.array(north), args.north, atol=1e-5))
        self.assertTrue(np.allclose(np.array(height), args.height, atol=1e-5))
        
        ##Check the tile (antenna) names are read correctly
        self.assertTrue(np.char.equal(np.array(names), args.ant_names).all())

    def test_EDA2_args_work(self):
        """Check that with a minimal set of arguments, the EDA2 primary
        beam is selected by `run_setup.check_args` correctly"""
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

        ##Should be pointed at phase centre with default 20.0 FWHM and ref freq 150e+6
        ##as none of these have been specified
        self.assertEqual(args.ra0, args.gauss_ra_point)
        self.assertEqual(args.dec0, args.gauss_dec_point)
        self.assertEqual(20.0, args.gauss_beam_FWHM)
        self.assertEqual(150e+6, args.gauss_beam_ref_freq)

        ##Keep adding in optional arguments and checking they end up manifesting
        self.inputs.append('--gauss_ra_point=234.0')
        args = self.run_parser_and_check_args()
        self.assertEqual(234.0, args.gauss_ra_point)
        self.assertEqual(args.dec0, args.gauss_dec_point)
        self.assertEqual(20.0, args.gauss_beam_FWHM)
        self.assertEqual(150e+6, args.gauss_beam_ref_freq)

        self.inputs.append('--gauss_dec_point=19.0')
        args = self.run_parser_and_check_args()
        self.assertEqual(234.0, args.gauss_ra_point)
        self.assertEqual(19.0, args.gauss_dec_point)
        self.assertEqual(20.0, args.gauss_beam_FWHM)
        self.assertEqual(150e+6, args.gauss_beam_ref_freq)

        self.inputs.append('--gauss_beam_FWHM=64.0')
        args = self.run_parser_and_check_args()
        self.assertEqual(234.0, args.gauss_ra_point)
        self.assertEqual(19.0, args.gauss_dec_point)
        self.assertEqual(64.0, args.gauss_beam_FWHM)
        self.assertEqual(150e+6, args.gauss_beam_ref_freq)

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

        ##Check the primary_beam has been selected and that `run_setup.check_args` fails
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

        ##Check the primary_beam has been selected and that `run_setup.check_args` fails
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


    def test_do_autos_work(self):
        """Check that with a minimal set of arguments, the parser defaults to
        no asking for autocorrelations, and check that this changes when
        requested"""

        ##Make minimum arguements, run parser, check do_autos == False
        self.make_minimum_required_args_without_metafits()
        args = self.run_parser_and_check_args()
        self.assertEqual(False, args.do_autos)

        ##Append the --do_autos flag, and check do_autos == True
        self.inputs.append('--do_autos')
        args = self.run_parser_and_check_args()
        self.assertEqual(True, args.do_autos)

    def test_sky_cropping_works(self):
        """Check that with a minimal set of arguments, the parser defaults to
        crop by sky component, and check that this changes when
        requested via --sky_crop_sources"""

        ##Make minimum arguements, run parser, check do_autos == False
        self.make_minimum_required_args_without_metafits()
        args = self.run_parser_and_check_args()
        self.assertEqual(True, args.sky_crop_components)

        ##Append the --do_autos flag, and check do_autos == True
        self.inputs.append('--sky_crop_sources')
        args = self.run_parser_and_check_args()
        self.assertEqual(False, args.sky_crop_components)
        
    def test_use_MWA_dipamps_works(self):
        """Check `ra.check_args` works correctly for the `--use_MWA_dipamps` flag.
        Should fail for a number of cases"""

        ##Make minimum arguements
        self.make_minimum_required_args_without_metafits()
        
        ##Set flag we want to test
        self.inputs.append('--use_MWA_dipamps')
        
        # ##Try running without selecting a FEE beam - should fail
        self.assert_check_args_errors()
        
        ##Put in a FEE beam, but not a metafits file - should fail
        self.inputs.append('--primary_beam=MWA_FEE')
        self.assert_check_args_errors()
        
        ##now try adding a metafits file that does not have the `Dipamps` key
        ##this should also fail
        self.inputs.append("--metafits_filename={:s}/1202815152_metafits_ppds.fits".format(code_dir))
        self.assert_check_args_errors()
        
        ##Add a metafits file that has the `Dipamps` key. However input
        ##args have a bespoke array layout which has 8 antennas, whereas this
        ##metafits has 128 antennas. Should fail
        
        self.inputs.append("--metafits_filename={:s}/1088285600_DipAmps.metafits".format(code_dir))
        self.assert_check_args_errors()
        
        ##reset to remove incorrect --array_layout flags
        self.make_required_args()
        self.inputs.append('--use_MWA_dipamps')
        self.inputs.append('--primary_beam=MWA_FEE')
        self.inputs.append("--metafits_filename={:s}/1088285600_DipAmps.metafits".format(code_dir))
        
        args = self.run_parser_and_check_args()
        
        ##check outputs are good
        self.assertEqual(args.use_MWA_dipamps, True)
        ##Using atol=1e-8 as only stored to 8 decimal places
        npt.assert_allclose(args.dipamps[dipamp_indexes], dipamps, atol=1e-8)
        
    def test_use_MWA_dipflags_works(self):
        """Check `ra.check_args` works correctly for the `--use_MWA_dipamps` flag.
        Should fail for a number of cases"""

        ##Make minimum arguements
        self.make_minimum_required_args_without_metafits()
        
        ##Set flag we want to test
        self.inputs.append('--use_MWA_dipflags')
        
        # ##Try running without selecting a FEE beam - should fail
        self.assert_check_args_errors()
        
        ##Put in a FEE beam, but not a metafits file - should fail
        self.inputs.append('--primary_beam=MWA_FEE')
        self.assert_check_args_errors()
        
        ##Do give it a metafits, but combined with the `--array_layout` flag
        ##this should crash as we have wrong number of tiles
        self.inputs.append("--metafits_filename={:s}/1088285600_DipAmps.metafits".format(code_dir))
        self.assert_check_args_errors()
        
        #reset to remove incorrect --array_layout flags
        self.make_required_args()
        self.inputs.append('--use_MWA_dipflags')
        self.inputs.append('--primary_beam=MWA_FEE')
        
        ##Try using a metafits file that has no dipole flagging at all
        ##As we haven't asked for dipole amplitudes here, we should stick
        ##in a warning saying we won't flag any dipoles and will run with
        ##perfect FEE beams (which will be way faster)
        self.inputs.append("--metafits_filename={:s}/1088285600_DipAmps.metafits".format(code_dir))
        args = self.run_parser_and_check_args()
        self.assertEqual(args.use_MWA_dipflags, False)
        
        ##Finally run with a metafites that does have flags and read them in
        self.inputs.append("--metafits_filename={:s}/1088285600_DipAmps_withflags.metafits".format(code_dir))
        args = self.run_parser_and_check_args()
        
        ##Righto, so our fina result should live in args.dipamps and we
        ##should have switched --use_MWA_dipamps to True
        self.assertEqual(args.use_MWA_dipflags,True)
        self.assertEqual(args.use_MWA_dipamps, True)
        
        ##Check we have a zero where we expect
        npt.assert_array_equal(dipflags, np.where(args.dipamps == 0)[0])
        
    def test_use_both_MWA_dipflags_dipamps_works(self):
        """Check `ra.check_args` works correctly for the `--use_MWA_dipflags` and
        `--use_MWA_dipamps` flags. Make sure it combines the arrays"""
        
        self.make_required_args()
        self.inputs.append('--use_MWA_dipflags')
        self.inputs.append('--use_MWA_dipamps')
        self.inputs.append('--primary_beam=MWA_FEE')
        self.inputs.append("--metafits_filename={:s}/1088285600_DipAmps_withflags.metafits".format(code_dir))
        args = self.run_parser_and_check_args()
        
        ##Righto, so our fina result should live in args.dipamps and we
        ##should have switched --use_MWA_dipamps to True
        self.assertEqual(args.use_MWA_dipflags,True)
        self.assertEqual(args.use_MWA_dipamps, True)
        
        ##Check we have a zero where we expect
        npt.assert_array_equal(dipflags, np.where(args.dipamps == 0)[0])
        npt.assert_allclose(args.dipamps[dipamp_indexes], dipamps_flagged, atol=1e-8)
        

##Run the test
if __name__ == '__main__':
    unittest.main()
    