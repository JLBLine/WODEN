from sys import path
import os
import unittest
import numpy as np
from ctypes import c_double
import argparse

##This is where our code lives
from wodenpy.use_libwoden.woden_settings import Woden_Settings_Double
from wodenpy.use_libwoden.array_layout_struct import setup_array_layout
from wodenpy.array_layout.create_array_layout import RTS_PrecessXYZtoJ2000

D2R = np.pi / 180.0


LST_BEFORE = 0.0077340
LST_AFTER = 0.0042458877378725
JD_DATE = 2457278.2010995

lsts = np.array([0.0045377974422681, 0.0051216168443364])
mjds = np.array([57277.7011457963599241, 57277.7012383889523335])


expec_X_noprec = np.array([53.009559921431, 93.452927296631, 173.879499021789,
                           106.473886729789, 258.350407484468, 270.515794526912,
                           31.455952402933, 242.617035144687,
                           53.009559921431, 93.452927296631, 173.879499021789,
                           106.473886729789, 258.350407484468, 270.515794526912,
                           31.455952402933, 242.617035144687])

expec_Y_noprec = np.array([84.00000, 12.00000, 86.00000, 780.00000,
                           813.00000, 899.00000, 460.00000, 810.00000,
                           84.00000, 12.00000, 86.00000, 780.00000,
                           813.00000, 899.00000, 460.00000, 810.00000])

expec_Z_noprec = np.array([98.706567952372, 179.107650254506, 334.544346566986,
                           200.542542730095,  498.021151109679,  535.557844598951,
                           62.534175123882, 464.518002081299,
                           98.706567952372, 179.107650254506, 334.544346566986,
                           200.542542730095,  498.021151109679,  535.557844598951,
                           62.534175123882, 464.518002081299])


expec_X_prec = np.array([53.162354565825, 93.726514308806, 174.392177908297,
                          106.800243325144, 259.131508165121, 271.356422590930,
                          31.563384834262, 243.346937978264, 53.162389323963,
                          93.726514245527, 174.392206378199, 106.800587900124,
                          259.131858510815, 271.356810457851, 31.563589743626,
                          243.347287996261])

expec_Y_prec = np.array([83.993507108853, 11.988289841444, 85.978144795225,
                          779.986839276406, 812.967476912920, 898.965215900110,
                          459.995941258304, 809.969621844301, 83.993395385724,
                          11.988088348935, 85.977768738665, 779.986612802000,
                          812.966917276866, 898.964617357155, 459.995871407581,
                          809.969099115337])

expec_Z_prec = np.array([98.629888069350, 178.965419626805, 334.283004207471,
                          200.420155127401, 497.668306161755, 535.190837477925,
                          62.509893326505, 464.189055736995, 98.629964478068,
                          178.965433157135, 334.283086077152, 200.420852891734,
                          497.669037933927, 535.191646196057, 62.510303874448,
                          464.189784358892])

##Vehicle for running tests
class Test(unittest.TestCase):
    def test_RTS_PrecessXYZtoJ2000(self):
        
        woden_settings = Woden_Settings_Double()
        ##Set up where woden_settings correctly
        woden_settings.latitude = -0.46606083776035967
        woden_settings.jd_date = 2457278.2010995
        woden_settings.lst_obs_epoch_base = LST_BEFORE
        woden_settings.lst_base = LST_AFTER
        woden_settings.do_precession = 1
        woden_settings.num_time_steps = 2
        woden_settings.time_res = 8

        ##These are ctype double arrays, and as ever we have to iterate
        ##over them to populate
        num_times_array = woden_settings.num_time_steps*c_double
        woden_settings.lsts = num_times_array()
        woden_settings.mjds = num_times_array()

        for time_ind in range(woden_settings.num_time_steps):
            woden_settings.lsts[time_ind] = lsts[time_ind]
            woden_settings.mjds[time_ind] = mjds[time_ind]

        ##Some things get setup after checking the input arguments, so
        ##quickly knock up a reduced version for this test
        parser = argparse.ArgumentParser()
        parser.add_argument('--num_antennas', type=int)
        # parser.add_argument('--east',)
        # parser.add_argument('--north',)
        # parser.add_argument('--height',)
        args = parser.parse_args(['--num_antennas=8'])

        # args.east = np.zeros(8)
        # args.north = np.zeros(8)
        # args.height = np.zeros(8)

        array_layout = setup_array_layout(woden_settings, args)
        # print(array_layout.num_baselines)

        for xyz_ind in range(len(expec_X_noprec)):
            array_layout.ant_X[xyz_ind] = expec_X_noprec[xyz_ind]
            array_layout.ant_Y[xyz_ind] = expec_Y_noprec[xyz_ind]
            array_layout.ant_Z[xyz_ind] = expec_Z_noprec[xyz_ind]


        ##code we are actually testing
        array_layout = RTS_PrecessXYZtoJ2000(array_layout, woden_settings)

        expec_num = woden_settings.num_time_steps*args.num_antennas

        found_X_prec = np.ctypeslib.as_array(array_layout.ant_X, shape=(expec_num, ))
        found_Y_prec = np.ctypeslib.as_array(array_layout.ant_Y, shape=(expec_num, ))
        found_Z_prec = np.ctypeslib.as_array(array_layout.ant_Z, shape=(expec_num, ))

        self.assertTrue(np.allclose(found_X_prec, expec_X_prec, atol=1e-10))
        self.assertTrue(np.allclose(found_Y_prec, expec_Y_prec, atol=1e-10))
        self.assertTrue(np.allclose(found_Z_prec, expec_Z_prec, atol=1e-10))

        



##Run the test
if __name__ == '__main__':
   unittest.main()
