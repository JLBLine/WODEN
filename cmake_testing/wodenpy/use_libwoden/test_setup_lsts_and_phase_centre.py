from sys import path
import os
import unittest
import numpy as np
import ctypes

# ##Code we are testing
import wodenpy.use_libwoden.woden_settings as ws
import wodenpy.wodenpy_setup.run_setup
# import woden_lib
from wodenpy.use_libwoden.create_woden_struct_classes import Woden_Struct_Classes

woden_struct_classes = Woden_Struct_Classes()
Woden_Settings = woden_struct_classes.Woden_Settings

D2R = np.pi/180.0
SOLAR2SIDEREAL = 1.00274
DS2R  = 7.2722052166430399038487115353692196393452995355905e-5

##Vehicle for running tests
class Test(unittest.TestCase):
    """Test the `rw.write_json` function, which writes an input file to feed
    into the WODEN executable"""
    
    def run_test(self, woden_settings, expec_prec_lsts):
        
        ##do it with precession
        woden_settings.do_precession = 1
        
        lsts, latitudes = ws.setup_lsts_and_phase_centre(woden_settings)
        
        self.assertTrue(np.array_equal(expec_prec_lsts, lsts))
        
        ##do it again with no precession
        
        woden_settings.do_precession = 0
        expected_lsts = np.empty(woden_settings.num_time_steps)
        
        for time_step in range(woden_settings.num_time_steps):
        
            expected_lsts[time_step] = woden_settings.lst_obs_epoch_base + (time_step + 0.5)*woden_settings.time_res*SOLAR2SIDEREAL*DS2R
        
        lsts, latitudes = ws.setup_lsts_and_phase_centre(woden_settings)
        self.assertTrue(np.allclose(expected_lsts, lsts, atol=1e-8))
        
    
    def test_setup_lsts_and_phase_centre_1(self):
        
        woden_settings = Woden_Settings()
        
        woden_settings.time_res = 1.0
        woden_settings.latitude = -0.46606083776035967
        woden_settings.latitude_obs_epoch_base = -0.46606083776035967
        
        woden_settings.lst_base = 0.0
        woden_settings.lst_obs_epoch_base = 0.0
        
        woden_settings.jd_date = 2457278.201145833
        woden_settings.num_time_steps = 3
        
        expec_prec_lsts = np.array([6.2797277312015503,
                                    6.2798007086408179,
                                    6.2798736860799780])
        
        self.run_test(woden_settings, expec_prec_lsts)
        
    def test_setup_lsts_and_phase_centre_2(self):
        
        woden_settings = Woden_Settings()
        
        woden_settings.time_res = 1.0
        woden_settings.latitude = -0.46606083776035967
        woden_settings.latitude_obs_epoch_base = -0.46606083776035967
        
        woden_settings.lst_base = 1.2345
        woden_settings.lst_obs_epoch_base = 1.2345
        
        woden_settings.jd_date = 2457278.201145833
        woden_settings.num_time_steps = 8
        
        expec_prec_lsts = np.array([1.2317545128009186,
                                    1.2318274511154514,
                                    1.2319003894260832,
                                    1.2319733277328140,
                                    1.2320462660356444,
                                    1.2321192043345734,
                                    1.2321921426296016,
                                    1.2322650809207287])
        
        self.run_test(woden_settings, expec_prec_lsts)
        
##Run the test
if __name__ == '__main__':
   unittest.main()