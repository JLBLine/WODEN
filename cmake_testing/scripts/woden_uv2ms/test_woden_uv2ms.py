from sys import path
import os
import unittest
import numpy as np
from astropy.io import fits
from unittest import mock
import argparse
from astropy.constants import c as speed_of_light
from subprocess import call

code_dir = os.path.realpath(__file__)
code_dir = ('/').join(code_dir.split('/')[:-1])

##Append the location of add_woden_uvfits to the sys.path to import from it
path.append('{:s}/../../../scripts'.format(code_dir))

##Code we are testing
from wodenpy.uvfits import wodenpy_uvfits
import woden_uv2ms
import numpy.testing as npt

##Vehicle for running tests
class Test(unittest.TestCase):
    def call_woden(self, band):
        """Call WODEN on the command line to make a uvfits file"""

        cmd = "run_woden.py"
        
        cmd += f" --band_nums={band} "
        cmd += "--lowest_channel_freq=180e6 "
        cmd += "--coarse_band_width=1.28e+6 "
        cmd += "--freq_res=80e+3 "
        cmd += "--num_freq_channels='obs' "
        cmd += "--num_time_steps=2 "
        cmd += "--time_res=2 "
        cmd += "--ra0=0.0 --dec0=0.0 "
        cmd += "--date 2000-01-01T12:00:00 "
        cmd += "--latitude=0.0 --longitude=79.53789457709424 --array_height=0.0 "
        cmd += f"--cat_filename={code_dir}/simple_sky.txt "
        cmd += "--output_uvfits_prepend=unittest "
        cmd += "--primary_beam='none' "
        cmd += "--no_precession "
        cmd += f"--array_layout={code_dir}/simple_array.txt "
        cmd += "--precision='double' "
        
        call(cmd, shell=True)

        
    ##Setup some fake command line things for the following test
    @mock.patch('argparse.ArgumentParser.parse_args',
                return_value=argparse.Namespace(single_uvfits=False,
                            uvfits_prepend="unittest_band",
                            band_nums="1,2,3",
                            no_delete=False))
        
    def test_convert_three_uvfits(self, mock_args):
        """Checks we can convert three uvfits files to measurement sets,
        by first generating the uvfits with run_woden.py. Check it works by
        asserting a measurement set has been created"""

        for band in range(1,4):
            self.call_woden(band)
        
        woden_uv2ms.main()
        
        for band in range(1,4):
            self.assertTrue(os.path.exists(f"unittest_band{band:02d}.ms"))
            
        call("rm -r unittest_band*.ms", shell=True)
        
    ##Setup some fake command line things for the following test
    @mock.patch('argparse.ArgumentParser.parse_args',
                return_value=argparse.Namespace(single_uvfits="unittest_band01.uvfits",
                            uvfits_prepend=False,
                            no_delete=False))
        
    def test_convert_single_uvfits(self, mock_args):
        """Checks we can convert three uvfits files to measurement sets,
        by first generating the uvfits with run_woden.py. Check it works by
        asserting a measurement set has been created"""

        self.call_woden(1)
        
        woden_uv2ms.main()
        
        self.assertTrue(os.path.exists(f"unittest_band01.ms"))
            
        call("rm -r unittest_band*.ms", shell=True)
        
    ##Setup some fake command line things for the following test
    @mock.patch('argparse.ArgumentParser.parse_args',
                return_value=argparse.Namespace(single_uvfits="unittest_band01.uvfits",
                            uvfits_prepend=False,
                            no_delete=True))
    
    def test_no_delete(self, mock_args):
        """Checks we can convert three uvfits files to measurement sets,
        by first generating the uvfits with run_woden.py. Check it works by
        asserting a measurement set has been created"""

        self.call_woden(1)
        
        woden_uv2ms.main()
        
        self.assertTrue(os.path.exists(f"unittest_band01.ms"))
            
        call("rm -r unittest_band*.ms", shell=True)
        
    ##Setup some fake command line things for the following test
    @mock.patch('argparse.ArgumentParser.parse_args',
                return_value=argparse.Namespace(single_uvfits="asdasdasd",
                            uvfits_prepend=False,
                            no_delete=True))
    
    def test_bad_path(self, mock_args):
        """Checks we can convert three uvfits files to measurement sets,
        by first generating the uvfits with run_woden.py. Check it works by
        asserting a measurement set has been created"""

        woden_uv2ms.main()
        
    ##Setup some fake command line things for the following test
    @mock.patch('argparse.ArgumentParser.parse_args',
                return_value=argparse.Namespace(single_uvfits=False,
                            uvfits_prepend="unittest_band",
                            band_nums="asdasdasd",
                            no_delete=False))
    
    def test_bad_band_nums(self, mock_args):
        """WHen band numbers are poopy, it should exit"""

        with self.assertRaises(SystemExit) as cm:

            woden_uv2ms.main()
        

##Run the test
if __name__ == '__main__':
   unittest.main()
