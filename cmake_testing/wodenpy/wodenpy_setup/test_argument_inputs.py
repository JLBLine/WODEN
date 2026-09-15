from sys import path
import os
import unittest
import numpy as np
import numpy.testing as npt
import importlib_resources
import shutil
from casacore.tables import table
from pyuvdata.telescopes import known_telescope_location

code_dir = os.path.realpath(__file__)
code_dir = ('/').join(code_dir.split('/')[:-1])

# ##Code we are testing
from wodenpy.wodenpy_setup import run_setup
from wodenpy.use_libwoden.use_libwoden import check_for_everybeam


woden_lib_path = importlib_resources.files('wodenpy').joinpath('libwoden_double.so')
HAVE_EVERYBEAM = check_for_everybeam(woden_lib_path)

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


east_ms = [-1.50306139e+02, -9.58861392e+01, -8.90331391e+01, -7.92191387e+01,
 -8.74611383e+01, -1.02846138e+02, -9.38671381e+01, -1.48851137e+02,
 -1.33460140e+02, -1.40388139e+02, -1.24188140e+02, -9.34101397e+01,
 -7.83561397e+01, -8.39971395e+01, -7.75091393e+01, -1.76527140e+02,
 -2.41411364e+01, -2.43231372e+01, -4.10101378e+01, -5.50481379e+01,
 -6.53721384e+01, -7.05961388e+01, -6.05781389e+01, -5.11971390e+01,
 -6.85111393e+01, -6.05791393e+01, -7.14061396e+01, -4.60611393e+01,
 -5.36771395e+01, -2.31801389e+01, -5.56181398e+01, -1.54211393e+01,
 -5.60811332e+01,  1.76945874e+02, -2.22311312e+01,  9.12878753e+01,
    1.53178883e+02, -6.81031271e+01, -1.18558126e+02, -8.80371332e+01,
 -1.05802137e+02, -9.90511376e+01, -8.23661375e+01, -8.00621379e+01,
 -7.56411383e+01, -7.15381378e+01, -6.30881375e+01, -5.12121372e+01,
 -1.96501136e+02, -1.98715132e+02, -2.84978120e+02, -3.88673126e+02,
 -4.53207127e+02, -4.69338136e+02, -3.65045136e+02, -3.70129139e+02,
 -1.60986151e+02, -1.28567142e+02, -7.02921402e+01, -7.92521403e+01,
 -1.04462140e+02, -1.00693140e+02, -3.95973143e+02, -2.64068144e+02,
    3.25619864e+02,  1.72903864e+02, -3.03614012e+00, -4.70751399e+01,
 -5.66621403e+01, -3.42261425e+01,  4.96058499e+01,  8.43588506e+01,
 -5.18567154e+02, -4.90464153e+02, -5.76078145e+02, -5.86196124e+02,
 -6.53988134e+02, -6.85548137e+02, -1.08012113e+03, -1.02659115e+03,
 -3.34279196e+02, -2.27672165e+02, -1.74238160e+02, -2.57473164e+02,
 -1.95539159e+02, -2.50168158e+02, -3.01217157e+02, -8.89992172e+02,
    7.15871814e+02,  2.87544847e+02,  2.37788840e+02,  2.27356838e+02,
    6.09218385e+01, -5.93951612e+01,  3.47028329e+01,  4.93814166e-01,
    8.39421832e+02,  1.17676884e+03,  1.31308885e+03,  5.55600862e+02,
    5.71019858e+02,  4.88335855e+02,  5.28076850e+02,  3.80547845e+02,
    1.10041887e+03,  1.20060889e+03,  4.98107883e+02,  3.00671874e+02,
    3.70126874e+02,  5.43802876e+02,  5.06271872e+02,  5.23205869e+02,
    6.53252915e+02,  4.42468907e+02, -5.93221070e+01,  3.65438888e+01,
    2.32320887e+02,  2.57462887e+02,  3.25166885e+02,  3.25579885e+02,
 -4.12978111e+02, -1.42419111e+02, -4.80791094e+02, -5.76054093e+02,
 -6.74542102e+02, -6.31924118e+02, -6.33758121e+02, -5.69061117e+02]

north_ms = [  265.85893427,   270.23293651,   266.06593679,   258.47593719,
    248.49693685,   244.13293622,   242.32393659,   220.31593433,
    282.65593496,   275.54293468,   284.89893534,   281.40293661,
    282.88393723,   277.607937,     272.68193726,   289.94793319,
    201.59093945,   222.39793944,   235.21493876,   237.56293818,
    250.85793776,   259.44693755,   263.72093796,   264.59093834,
    271.14793763,   271.94993796,   278.86293751,   273.28993855,
    276.31293824,   263.73893949,   284.32993816,   273.91493981,
    124.84593814,   -43.3320523,     77.04993952,   -80.52005582,
  -259.88405327,   -22.31606236,   -39.82806443,   125.24193683,
    217.2249361,    230.20693638,   229.38193706,   238.17293716,
    247.04693734,   235.97393751,   229.03193785,   221.92993834,
    189.83993237,    95.19493228,  -196.13607128,   -43.97707553,
    -15.88807818,   185.23192115,   183.14292544,   266.52892523,
    573.25193382,   350.53893516,   294.07593756,   297.30393719,
    300.95893615,   288.94093631,   371.19692416,   389.31892959,
    203.65295379,   193.63994753,   292.88594032,   287.51393851,
    295.96693812,   351.09293904,   536.20394247,   519.1859439,
    632.54791912,   604.61292027,   415.05291676,  -101.48508364,
    142.29991357,   212.21291227,    87.06289609,   574.57189826,
  1660.77492667,   898.17393107,   781.68293327,   892.04992984,
    762.93793239,   731.75393015,   704.37092805,  1088.18490385,
  1405.30496982,   618.30995224,   768.4579502,    819.76894977,
    815.29994294,   810.57493799,   950.28694186,  1412.40494044,
    990.27197488,   858.67698873,   553.41699432,   244.52996323,
    359.25496386,   425.02496047,   539.2389621,    652.90895605,
    151.80498558,  -443.20501031,  -264.54903913,   -48.90704723,
    -38.37004438,   -91.10703726,    -3.1710388,     75.5389619,
 -1037.70503277,  -835.26104142,  -505.28706202,  -405.00505808,
  -357.84605003,  -366.408049,    -318.46304622,  -306.2620462,
  -409.01507654,  -420.21606543,  -817.27907931,  -863.73708322,
  -638.82608727,  -245.06208552,  -169.31708559,  -257.61508294]


height_ms = [754.01098587, 754.28499028, 754.30899087, 754.29099173, 754.33899112,
 754.2579899, 754.26999064, 754.21398631, 754.00198713, 754.01798661,
 754.04998787, 754.25399041, 754.25299162, 754.2999912, 754.29799177,
 753.62298356, 753.80499663, 753.87199646, 753.96099501, 754.12299385,
 754.23899291, 754.26599242, 754.17999321, 754.09999397, 754.23599251,
 754.20399315, 754.25199222, 754.06599433, 754.08999369, 753.94799626,
 754.18099347, 753.90999683, 753.80299456, 752.00501478, 753.35099766,
 752.76000805, 753.47601437, 753.17499462, 753.05399063, 754.03099195,
 754.25498985, 754.25999031, 754.26999168, 754.2699918, 754.2549921,
 754.23899251, 754.20999325, 754.07799427, 754.16098264, 754.02398313,
 751.90097815, 752.3719686, 752.11096314, 751.7519604, 752.73896892,
 752.23996792, 751.77298282, 753.56998705, 754.1919922, 754.19599145,
 754.04998936, 754.15698976, 751.51496507, 752.32197571, 750.94502517,
 752.03201277, 753.8919977, 754.11499415, 754.1559933, 753.97599474,
 753.18100028, 753.31400324, 749.60395321, 749.9069557, 750.37395005,
 752.21195288, 751.46194562, 751.23594255, 754.00891123, 751.03891215,
 748.06796098, 750.69397508, 751.21698027, 750.49197269, 751.09297866,
 750.79697442, 750.51397045, 747.38091966, 751.60204852, 752.13401912,
 752.805014, 752.74801278, 753.00099923, 752.16998944, 753.27499613,
 750.68199007, 749.69606154, 749.12209002, 747.35410331, 749.28404365,
 749.5090441, 749.96703688, 749.78903932, 751.25102647, 746.06508879,
 745.32810118, 750.22004256, 751.61802492, 751.34003051, 749.81204506,
 749.80304138, 749.6130422, 748.0710607, 749.25104206, 749.91299876,
 751.03400587, 752.26202152, 752.22102364, 752.02902882, 752.12802877,
 750.9689692, 750.50299137, 752.71996656, 752.52195911, 752.33994948,
 752.54395017, 752.61594948, 752.24595539]

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


dipflags = np.array([10, 43, 102, 110, 113, 120, 123, 198, 202, 204, 205, 213, 214, 221,
                    223,  431,  534,  616, 1028, 1048, 1059, 1074, 1092, 1099, 1101, 1122, 1148, 1163,
                    1179, 1215, 1217, 1268, 1350, 1386, 1429, 1572, 1599, 1613, 1700, 1740, 1744, 1811,
                    2156, 2311, 2323, 2398, 2418, 2496, 2569, 2629, 2652, 2657, 2667, 2703, 2719, 2745,
                    2816, 2842, 2864, 2997, 3102, 3247, 3260, 3265, 3308, 3343, 3384, 3420, 3436, 3480,
                    3561, 3618, 3633, 3663, 3700, 3772, 3791, 3823, 3880, 3897, 3969, 4003, 4089])

dipamp_indexes = np.array([ 10, 43, 102, 100, 345, 768, 2780, 3678, 4000])
dipamps = np.array([0.90376163, 0.95701027, 0.44699562, 0.95701027,
                    1., 0.97801661, 0.95163655, 1., 1.])

dipamps_flagged = np.array([0.0, 0.0, 0.0, 0.95701027,
                    1., 0.97801661, 0.95163655, 1., 1.])


MWA_LAT = -26.703319405555554
MWA_LONG = 116.67081523611111
MWA_HEIGHT = 377.827

LOFAR_LAT = 52.905329712
LOFAR_LONG = 6.867996528
LOFAR_HEIGHT = 0.0

##Vehicle for running tests
class Test(unittest.TestCase):
    """Test whether the args collected by the argument parser are read in
    correctly, and that sanity checks on certain combinations work such that
    we don't feed WODEN arguments that won't work"""
    
    def setUp(self):
        """Grab two environment variables. We're going to break these during
        tests to check things fail, so we need to restore them afterwards
        We do that via the `addCleanup` method, which will run the function
        passed to it after the test has finished. This is a unittest feature"""
        self.hdf5_path = os.environ['MWA_FEE_HDF5']
        self.hdf5_interp_path = os.environ['MWA_FEE_HDF5_INTERP']
            
        self.addCleanup(self.reset_env_var)  # Ensure cleanup runs

    def reset_env_var(self):
        """Restore the environment variables to their original values"""
        os.environ["MWA_FEE_HDF5"] = self.hdf5_path  # Restore original value
        os.environ["MWA_FEE_HDF5_INTERP"] = self.hdf5_interp_path  # Restore original value

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
        
    def make_required_args(self):
        """These are arguments that if missing, the parser itself should fail"""
        self.inputs = ['--ra0=0.0', '--dec0=-26.7', '--cat_filename=srclist.txt']

    def make_minimum_required_args_minus_radec(self):
        """When not passing a metafits file, these are the minimum set of
        arguments needed to pass onto WODEN. Other arguments either have
        defaults or are optional"""
        self.inputs = ['--cat_filename=srclist.txt']

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

        actual_hdf5_path = os.environ[beam_env_var]
        # self.actual_hdf5_path = actual_hdf5_path

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
        
        ##UVBeam MWA can only work with one band number, so add that arg here.
        ##Doesn't change the test outputs for the other beam models
        self.inputs.append('--band_nums=1')
        
        ##Check the primary_beam has been selected and that `run_setup.check_args` fails
        args = self.assert_check_args_errors()
        self.assertEqual(beam_name, args.primary_beam)
        
        
        self.inputs.append('--hdf5_beam_path=nuthin')
        self.assert_check_args_errors()

        ##Set the path to the hdf5 file. Just point it at a text file.
        ##TODO - have `check_args` actually test reading the hdf5 file. At the
        ##mo just tests that the file exists
        ##This should still fail because we have no delays set
        self.inputs.pop(-1)
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
        ##UVBeam MWA can only work with one band number, so add that arg here.
        ##Doesn't change the test outputs for the other beam models
        self.inputs.append('--band_nums=1')
        self.inputs.append(f'--primary_beam={beam_name}')
        self.inputs.append("--metafits_filename={:s}/1202815152_metafits_ppds.fits".format(code_dir))
        args = self.assert_check_args_errors()
        
        
        ##This time, set the environment variable to a non-existant file.
        ##Should faiul
        os.environ[beam_env_var] = 'this_aint_real'
        self.assert_check_args_errors()

        ##This time, set the environment variable to the hdf5 file. Should
        ##pass (just use a file we know exists somewhere)
        os.environ[beam_env_var] = actual_hdf5_path
        
        args = self.run_parser_and_check_args()
        ##Assert the delays were read in correctly
        self.assertEqual("6,4,2,0,8,6,4,2,10,8,6,4,12,10,8,6", args.MWA_FEE_delays)

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
        
    def test_UVBeamMWA_args_work(self):
        """Check that the UVBeam MWA primary beam is handled `ra.check_args`
        correctly. The function should error out if certain paths to the
        hdf5 file that holds the spherical harmonic information is missing,
        and if the delays have been specified incorrectly"""

        self._check_MWA_FEE_generic('uvbeam_MWA', 'MWA_FEE_HDF5')


    def test_MWAAnalyBeam_args_work(self):
        """Check `ra.check_args` works correctly for the analytic MWA beam.
        Should fail if the delays have been specified incorrectly"""

        ##Trying to run an MWA_FEE simulation with no metafits or
        ##delays should fail
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
        
        # self.make_required_args()
        self.make_minimum_required_args_without_metafits()
        self.inputs.append('--use_MWA_dipflags')
        self.inputs.append('--use_MWA_dipamps')
        self.inputs.append('--primary_beam=MWA_FEE')
        self.inputs.append('--MWA_FEE_delays=0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0')
        
        ##should error if we don't give it a metafits file
        self.assert_check_args_errors()
        
        ##reset args. give it a metafits file without Delay column, should error
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
        
    def test_use_off_cardinal_dipoles(self):
        """Check `ra.check_args` works correctly for the `--use_MWA_dipflags` and
        `--use_MWA_dipamps` flags. Make sure it combines the arrays"""
        
        ##First, run as normal, and assert that off cardinal dipoles are not used
        self.make_minimum_required_args_without_metafits()
        self.inputs.append('--primary_beam=None')
        args = self.run_parser_and_check_args()
        self.assertEqual(args.off_cardinal_dipoles,False)
        
        ##Next, add --off_cardinal_dipoles and check it is set to True
        self.make_minimum_required_args_without_metafits()
        self.inputs.append('--primary_beam=None')
        self.inputs.append('--off_cardinal_dipoles')
        args = self.run_parser_and_check_args()
        self.assertEqual(args.off_cardinal_dipoles,True)
        
    def test_use_everybeam_MWA(self):
        """Check `ra.check_args` works correctly for the `--primary_beam=everybeam_MWA`
        
        However, can only run these tests if we have EveryBeam. If not, we
        expect and error, so adjust the test accordingly
        """
        
        if HAVE_EVERYBEAM:
        
            self.make_minimum_required_args_without_metafits()
            self.inputs.append('--primary_beam=everybeam_MWA')
            ##check it errors as we haven't provided a hdf5 path or metafits
            self.assert_check_args_errors()
            
            ##check it errors with bad path to hdf5
            self.inputs.append('--hdf5_beam_path=hdf5_beam_path')
            self.assert_check_args_errors()
            
            
            ##reset input args
            self.make_minimum_required_args_without_metafits()
            self.inputs.append('--primary_beam=everybeam_MWA')
            self.inputs.append("--metafits_filename={:s}/1202815152_metafits_ppds.fits".format(code_dir))
            
            args = self.run_parser_and_check_args()
            ##Make sure delays were read in from the metafits, and the 
            ##amplitudes were set to 1.0
            self.assertEqual("6,4,2,0,8,6,4,2,10,8,6,4,12,10,8,6", args.MWA_FEE_delays)
            
            npt.assert_array_equal(np.ones(16), args.dipamps)

            
        ##If we don't have EveryBeam, we expect an error when someone asks
        ##for it
        else:
            self.make_minimum_required_args_without_metafits()
            self.inputs.append('--primary_beam=everybeam_MWA')
            self.assert_check_args_errors()
        
    def test_bad_primary_beam(self):
        """Check `ra.check_args` works correctly for the `--primary_beam=nonsense`"""
        
        self.make_minimum_required_args_without_metafits()
        self.inputs.append('--primary_beam=nonsense')
        self.assert_check_args_errors()
        
    def test_bad_no_beam_ms_path(self):
        """Check things fail if using no primary beam, but trying to get the
        array layout from a measurement set"""
        
        if HAVE_EVERYBEAM:
        
            self.make_minimum_required_args_without_metafits()
            self.inputs.append('--primary_beam=none')
            self.inputs.append('--beam_ms_path=lkjsdfsjhfkjshdf')
            self.assert_check_args_errors()
            
            ##now put in a proper path to the ms and check all good
            self.inputs.pop(-1)
            ms_path = f"{code_dir}/../../../test_installation/everybeam/MWA-single-timeslot.ms"
            self.inputs.append(f'--beam_ms_path={ms_path}')
            self.run_parser_and_check_args()
        
    def test_num_freq_chans(self):
        """Check that the number of frequency channels is read in correctly"""
        
        self.make_minimum_required_args_without_metafits()
        self.inputs.append('--num_freq_channels=128')
        args = self.run_parser_and_check_args()
        self.assertEqual(128, args.num_freq_channels)
        
    def test_bad_band_nums(self):
        """Check that the band numbers fail if non-sensical"""
        
        self.make_minimum_required_args_without_metafits()
        self.inputs.append('--band_nums=kjhasdkjahsd')
        self.assert_check_args_errors()
        
    def test_wrong_length_delays(self):
        """Check that the delays fail if the wrong length"""
        
        self.make_minimum_required_args_without_metafits()
        self.inputs.append('--MWA_FEE_delays=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15')
        self.assert_check_args_errors()
        
    def test_bad_precision(self):
        """Check that the precision defaults to double if nonsense is given"""
        
        self.make_minimum_required_args_without_metafits()
        self.inputs.append('--precision=yo_momma')
        args = self.run_parser_and_check_args()
        
        self.assertEqual(args.precision, 'double')
        
    def test_bad_beam_index(self):
        """Check that the beam index fails if too large"""
        
        if HAVE_EVERYBEAM:
        
            self.make_minimum_required_args_without_metafits()
            self.inputs.append('--primary_beam=everybeam_MWA')
            self.inputs.append('--band_nums=1')
            ms_path = f"{code_dir}/../../../test_installation/everybeam/MWA-single-timeslot.ms"
            self.inputs.append(f'--beam_ms_path={ms_path}')
            self.inputs.append('--station_id=9393')
            ##we need to get rid of the --array_layout flag, as it will
            ##override reading the measurement set array layout
            for ind, arg in enumerate(self.inputs):
                if arg.startswith('--array_layout'):
                    self.inputs.pop(ind)
                    break
            self.assert_check_args_errors()
        
    def test_bad_enh_text_file(self):
        """Check that the enhancement text file fails if it doesn't exist"""
        
        self.make_minimum_required_args_without_metafits()
        for ind, arg in enumerate(self.inputs):
            if arg.startswith('--array_layout'):
                self.inputs[ind] = f'--array_layout=non_existent_file.txt'
        self.assert_check_args_errors()
        
    def check_ms_pointing(self, ms_path, ra0, dec0):
        with table(ms_path+'::FIELD') as field_table:
        
            delay_dir = np.squeeze(field_table.getcol('DELAY_DIR'))
            ref_dir = np.squeeze(field_table.getcol('REFERENCE_DIR'))
            
            self.assertAlmostEqual(delay_dir[0], ra0, 6)
            self.assertAlmostEqual(delay_dir[1], dec0, 6)
            self.assertAlmostEqual(ref_dir[0], ra0, 6)
            self.assertAlmostEqual(ref_dir[1], dec0, 6)
        
    def test_read_radec_from_MS(self):
        """Check we can read the ra,dec pointing from a measurement set,
        and also check we can override with a new phase centre and/or beam
        pointing
        
        However, can only run these tests if we have EveryBeam. If not, we
        expect and error, so adjust the test accordingly
        """
        
        if HAVE_EVERYBEAM:
        
            self.make_minimum_required_args_minus_radec()
            self.inputs.append('--primary_beam=everybeam_LOFAR')
            ms_path = f"{code_dir}/../../../test_installation/everybeam/LOFAR_LBA_MOCK.ms"
            self.inputs.append(f'--beam_ms_path={ms_path}')
            self.inputs.append(f'--band_nums=1')
            
            ##we need to get rid of the --array_layout flag, as it will
            ##override reading the measurement set array layout
            for ind, arg in enumerate(self.inputs):
                if arg.startswith('--array_layout'):
                    self.inputs.pop(ind)
                    break
            
            args = self.run_parser_and_check_args()
            # print(args.ra0, args.dec0)
            
            self.assertAlmostEqual(args.ra0, -82.61757960000001, 6)
            self.assertAlmostEqual(args.dec0, 48.7461556, 6)
            
            self.assertAlmostEqual(args.eb_ra_point, -82.61757960000001, 6)
            self.assertAlmostEqual(args.eb_dec_point, 48.7461556, 6)
            
            ##Now check we can override the defaults
            self.inputs.append('--ra0=0.0')
            args = self.run_parser_and_check_args()
            
            self.assertAlmostEqual(args.ra0, 0.0, 6)
            self.assertAlmostEqual(args.dec0, 48.7461556, 6)
            
            self.assertAlmostEqual(args.eb_ra_point, -82.61757960000001, 6)
            self.assertAlmostEqual(args.eb_dec_point, 48.7461556, 6)
            
            ##Now check we can override the defaults
            self.inputs.append('--dec0=-26.7')
            args = self.run_parser_and_check_args()
            
            self.assertAlmostEqual(args.ra0, 0.0, 6)
            self.assertAlmostEqual(args.dec0, -26.7, 6)
            
            self.assertAlmostEqual(args.eb_ra_point, -82.61757960000001, 6)
            self.assertAlmostEqual(args.eb_dec_point, 48.7461556, 6)
            
            #The above has changed the phase centre; now check we can
            #set the pointing to be the same as the phase centre
            self.inputs.append('--eb_point_to_phase')
            args = self.run_parser_and_check_args()
            
            self.assertAlmostEqual(args.ra0, 0.0, 6)
            self.assertAlmostEqual(args.dec0, -26.7, 6)
            
            self.assertAlmostEqual(args.eb_ra_point, 0.0, 6)
            self.assertAlmostEqual(args.eb_dec_point, -26.7, 6)
            
            self.check_ms_pointing("pointed_output_band01.ms", 0.0, np.radians(-26.7))
            
            self.inputs.append('--eb_ra_point=20.5')
            args = self.run_parser_and_check_args()
            
            self.assertAlmostEqual(args.ra0, 0.0, 6)
            self.assertAlmostEqual(args.dec0, -26.7, 6)
            
            self.assertAlmostEqual(args.eb_ra_point, 20.5, 6)
            self.assertAlmostEqual(args.eb_dec_point, -26.7, 6)
            
            self.check_ms_pointing("pointed_output_band01.ms", np.radians(20.5), np.radians(-26.7))
            
            self.inputs.append('--eb_dec_point=-53.7')
            args = self.run_parser_and_check_args()
            
            self.assertAlmostEqual(args.ra0, 0.0, 6)
            self.assertAlmostEqual(args.dec0, -26.7, 6)
            
            self.assertAlmostEqual(args.eb_ra_point, 20.5, 6)
            self.assertAlmostEqual(args.eb_dec_point, -53.7, 6)
            
            self.check_ms_pointing("pointed_output_band01.ms", np.radians(20.5), np.radians(-53.7))
            
            shutil.rmtree("pointed_output_band01.ms")
            
    def test_read_inputs_from_metafits(self):
        
        # ##First double check it fails without ra,dec
        
        self.inputs = ['--cat_filename=srclist.txt']
        self.inputs.append('--array_layout={:s}/example_array_layout.txt'.format(code_dir))
        self.assert_check_args_errors()
        
        ##Now check it works with the metafits file
        self.inputs.append("--metafits_filename={:s}/1202815152_metafits_ppds.fits".format(code_dir))
        
        args = self.run_parser_and_check_args()
        
        self.assertAlmostEqual(args.ra0, 54.75681779724241, 6)
        self.assertAlmostEqual(args.dec0, -39.48422590285089, 6)
        
        self.assertAlmostEqual(args.lowest_channel_freq, 169595000.0, 2)
        self.assertAlmostEqual(args.num_time_steps, 240, 2)
        self.assertAlmostEqual(args.freq_res, 10000.0, 2)
        self.assertAlmostEqual(args.time_res, 0.5, 2)
        self.assertEqual("2018-02-16T11:18:54", args.date)
        
    def test_output_uvfits_prepend(self):
        """Check that the output_uvfits_prepend is set correctly, a.k.a
        if the user actually includes '.uvfits' that we strip it out"""
        
        self.make_minimum_required_args_without_metafits()
        
        self.inputs.append('--output_uvfits_prepend=some_cool_name.uvfits')
        
        args = self.run_parser_and_check_args()
        self.assertEqual("some_cool_name", args.output_uvfits_prepend)
        
    def test_setting_lon_lat_height(self):
        """Check latitude / longitude / height are set correctly based
        on requested beam type.
        
        However, can only run these tests if we have EveryBeam. If not, we
        expect and error, so adjust the test accordingly
        """
        
        if HAVE_EVERYBEAM:
            ##Doing LOFAR should stick the beam at the LOFAR lat/long/height
            self.make_minimum_required_args_minus_radec()
            self.inputs.append('--primary_beam=everybeam_LOFAR')
            ms_path = f"{code_dir}/../../../test_installation/everybeam/LOFAR_HBA_MOCK.ms"
            self.inputs.append(f'--beam_ms_path={ms_path}')
            self.inputs.append(f'--band_nums=1')
            
            args = self.run_parser_and_check_args()
            
            self.assertEqual(args.orig_lat, LOFAR_LAT)
            self.assertEqual(args.orig_long, LOFAR_LONG)
            self.assertEqual(args.orig_height, LOFAR_HEIGHT)
            
            self.assertEqual(args.latitude, LOFAR_LAT)
            self.assertEqual(args.longitude, LOFAR_LONG)
            self.assertEqual(args.array_height, LOFAR_HEIGHT)
            
            ##If we do it again, but set the lat/long/height manually, it should
            ##use those instead. However original lat/long/height should
            ##still be the LOFAR ones
            
            self.inputs.append(f'--latitude={MWA_LAT}')
            self.inputs.append(f'--longitude={MWA_LONG}')
            self.inputs.append(f'--array_height={MWA_HEIGHT}')
            
            args = self.run_parser_and_check_args()
            
            self.assertEqual(args.orig_lat, LOFAR_LAT)
            self.assertEqual(args.orig_long, LOFAR_LONG)
            self.assertEqual(args.orig_height, LOFAR_HEIGHT)
            
            self.assertEqual(args.latitude, MWA_LAT)
            self.assertEqual(args.longitude, MWA_LONG)
            self.assertEqual(args.array_height, MWA_HEIGHT)
            
        ##All other beams currently default to the MWA lat/long/height
        self.make_minimum_required_args_without_metafits()
        self.inputs.append('--primary_beam=Gaussian')
        
        args = self.run_parser_and_check_args()
        self.assertEqual(args.orig_lat, MWA_LAT)
        self.assertEqual(args.orig_long, MWA_LONG)
        self.assertEqual(args.orig_height, MWA_HEIGHT)
        
        self.assertEqual(args.latitude, MWA_LAT)
        self.assertEqual(args.longitude, MWA_LONG)
        self.assertEqual(args.array_height, MWA_HEIGHT)
        
    def test_check_asking_for_log_does_something(self):
        """Check that asking for a log file actually does something"""
        
        self.make_minimum_required_args_without_metafits()
        self.inputs.append('--save_log')
        
        args = self.run_parser_and_check_args()
        
        ##Check that the log file name is set. It will be False if not set.
        self.assertTrue(args.log_file_name)
        
    def test_use_uvdata_HERA_beam(self):
        """Check that asking for a HERA beam demands you give a
        text file listing CST files. Should also error if bad file is linked,
        or if it can't read the file correctly"""
        
        self.make_minimum_required_args_without_metafits()
        self.inputs.append('--primary_beam=uvbeam_HERA')
        ##should die if --cst_file_list not set
        self.assert_check_args_errors()
        
        
        self.inputs.append('--cst_file_list=note_a_file.txt')
        ##should die if --cst_file_list isn't a file
        self.assert_check_args_errors()
        
        ##Get rid of old --cst_file_list arg, replace with a bad file
        ##check it fails
        self.inputs.pop(-1)
        self.inputs.append(f'--cst_file_list={code_dir}/bad_eg_cst_file_list.csv')
        self.assert_check_args_errors()
        
        ##should run now. Check certain things are set successfully
        self.inputs.pop(-1)
        self.inputs.append(f'--cst_file_list={code_dir}/eg_cst_file_list.csv')
        
        args = self.run_parser_and_check_args()
        
        self.assertEqual(args.primary_beam, 'uvbeam_HERA')
        self.assertEqual(args.cst_paths, ['/hey/look/here/file 100 [1].txt', '/hey/look/here/file 150 [1].txt', '/hey/look/here/file 200 [1].txt'])
        self.assertEqual(args.cst_freqs, [100e+6, 150e+6, 200e+6])
        
        location = known_telescope_location("HERA")
        
        self.assertEqual(args.latitude, location.lat.value)
        self.assertEqual(args.longitude, location.lon.value)
        self.assertEqual(args.array_height, location.height.value)
        
        
        
        ##Now check the single file settings, rather than the CST list
        ##remove the cst option, and add the single file option
        ##add a borked path to check that this fails
        self.inputs.pop(-1)
        self.inputs.append(f'--uvbeam_file_path=note_a_file.txt')
        self.assert_check_args_errors()
        
        ##parser only checks if file exists, not if it is a valid (UVBeam can do that)
        ##so we can check things run here via an existing path
        self.inputs.pop(-1)
        self.inputs.append(f'--uvbeam_file_path={code_dir}/eg_cst_file_list.csv')
        
        args = self.run_parser_and_check_args()
        
        self.assertEqual(args.primary_beam, 'uvbeam_HERA')
        self.assertEqual(args.uvbeam_file_path, f'{code_dir}/eg_cst_file_list.csv')
        self.assertEqual(args.latitude, location.lat.value)
        self.assertEqual(args.longitude, location.lon.value)
        self.assertEqual(args.array_height, location.height.value)
        
        ##have made behaviour to default to single file if both
        ##--cst_file_list and --uvbeam_file_path are set, so check that
        
        self.inputs.append(f'--cst_file_list={code_dir}/eg_cst_file_list.csv')
        
        self.assertEqual(args.primary_beam, 'uvbeam_HERA')
        self.assertEqual(args.uvbeam_file_path, f'{code_dir}/eg_cst_file_list.csv')
        self.assertEqual(args.latitude, location.lat.value)
        self.assertEqual(args.longitude, location.lon.value)
        self.assertEqual(args.array_height, location.height.value)
        
        ##assert that args.cst_paths doesn't exist
        self.assertFalse(hasattr(args, 'cst_paths'))
        self.assertFalse(hasattr(args, 'cst_freqs'))
        

        
##Run the test
if __name__ == '__main__':
    unittest.main()
    