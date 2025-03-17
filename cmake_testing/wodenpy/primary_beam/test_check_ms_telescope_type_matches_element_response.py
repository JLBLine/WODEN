"""
Test `wodenpy.primary_beam.use_everybeam.run_everybeam_over_threads`, which
should run `wodenpy.primary_beam.use_everybeam.run_everybeam` over multiple
threads. This is a test of the parallelisation of the code.
"""

from sys import path
import os
import unittest
import numpy as np
import numpy.testing as npt

from wodenpy.use_libwoden.skymodel_structs import Components_Python
from astropy.coordinates import EarthLocation
from astropy import units as u
from astropy.time import Time
from wodenpy.primary_beam.use_everybeam import check_ms_telescope_type_matches_element_response

from subprocess import call

##Location of this file; use it to find test measurement sets
code_dir = os.path.realpath(__file__)
code_dir = ('/').join(code_dir.split('/')[:-1])


# path.append('{:s}/../../../scripts'.format(code_dir))

# import woden_uv2ms

##Vehicle for running tests
class Test(unittest.TestCase):
    """Vehicle for running tests"""
    
    def run_the_test(self, ms_path, expected_telescope_type, default_element,
                     accepted_elements):
        
                ##Default element response
        telescope_type, checked_response = check_ms_telescope_type_matches_element_response(ms_path)
        self.assertEqual(telescope_type, expected_telescope_type)
        self.assertEqual(checked_response, default_element)
        
        ##Specify certain allowed beam models; check they come back as expected
        
        for model in accepted_elements:
        
            telescope_type, checked_response = check_ms_telescope_type_matches_element_response(ms_path,
                                                                                                model)
            self.assertEqual(telescope_type, expected_telescope_type)
            self.assertEqual(checked_response, model)
        
        
        ##Ask for something wrong, should change it to hamaker with a warning
        telescope_type, checked_response = check_ms_telescope_type_matches_element_response(ms_path,
                                                                                            "wibbly")
        self.assertEqual(telescope_type, expected_telescope_type)
        self.assertEqual(checked_response, default_element)
    
    def test_check_MWA(self):
        """Check we get expected behaviour with an MWA MS"""
        
        
        ms_path = f'{code_dir}/../../../test_installation/everybeam/MWA-single-timeslot.ms'
        
        self.run_the_test(ms_path, "MWA", "MWA", ["MWA"])

        
    def test_check_LOFAR(self):
        """Check we get expected behaviour with a LOFAR MS"""
        
        
        ms_path = f'{code_dir}/../../../test_installation/everybeam/LOFAR_HBA_MOCK.ms'
        
        self.run_the_test(ms_path, "LOFAR", "hamaker",
                          ['hamaker', 'hamakerlba', 'lobes'])
        
    def test_check_OSKAR(self):
        """Check we get expected behaviour with an OSKAR MS"""
        
        ms_path = f'{code_dir}/../../../test_installation/everybeam/create_OSKAR-SKA_ms/OSKAR-SKA-layout.ms'
        self.run_the_test(ms_path, "OSKAR", "skala40_wave", ["skala40_wave"])
        
    def test_it_errors_with_bad_MS(self):
        """Check we get an error if we give it a bad MS"""
        
        cmd = "run_woden.py"
        
        cmd += f" --band_nums=1 "
        cmd += "--lowest_channel_freq=180e6 "
        cmd += "--coarse_band_width=1.28e+6 "
        cmd += "--freq_res=80e+3 "
        cmd += "--num_freq_channels='obs' "
        cmd += "--num_time_steps=2 "
        cmd += "--time_res=2 "
        cmd += "--ra0=0.0 --dec0=0.0 "
        cmd += "--date 2000-01-01T12:00:00 "
        cmd += "--latitude=0.0 --longitude=79.53789457709424 --array_height=0.0 "
        cmd += f"--cat_filename={code_dir}/../../scripts/woden_uv2ms/simple_sky.txt "
        cmd += "--output_uvfits_prepend=unittest "
        cmd += "--primary_beam='none' "
        cmd += "--no_precession "
        cmd += f"--array_layout={code_dir}/../../scripts/woden_uv2ms/simple_array.txt "
        cmd += "--precision='double' "
        cmd += "--num_threads=1 "
        cmd += "--telescope_name=what_is_this "
        
        call(cmd, shell=True)
        
        cmd = "woden_uv2ms.py --single_uvfits unittest_band01.uvfits"
        call(cmd, shell=True)
        
        with self.assertRaises(SystemExit) as cm:
            telescope_type, checked_response = check_ms_telescope_type_matches_element_response("unittest_band01.ms")
            
        call("rm unittest_band01.uvfits", shell=True)
        call("rm -r unittest_band01.ms", shell=True)
        
##Run the test
if __name__ == '__main__':
    unittest.main()