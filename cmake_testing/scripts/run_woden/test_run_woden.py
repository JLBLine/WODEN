from sys import path
import os
import unittest
import numpy as np
from astropy.io import fits
from unittest import mock
import argparse
from astropy.time import Time, TimeDelta
import numpy as np
from astropy.coordinates import EarthLocation
from astropy import units as u
from astropy.constants import c as speed_of_light

from wodenpy.array_layout.precession import RTS_Precess_LST_Lat_to_J2000
from astropy.time import Time, TimeDelta
import numpy as np
from astropy.coordinates import EarthLocation
from astropy import units as u
import scipy.optimize as opt
from glob import glob

code_dir = os.path.realpath(__file__)
code_dir = ('/').join(code_dir.split('/')[:-1])

##Append the location of add_woden_uvfits to the sys.path to import from it
path.append('{:s}/../../../scripts'.format(code_dir))

##Code we are testing
from wodenpy.uvfits import wodenpy_uvfits
import run_woden as rw
import numpy.testing as npt

# from multiprocessing import Manager, set_start_method

from subprocess import call

def get_lon(inputs, lon):
    """Every time they update astropy, the exact values of LST seems to change
    This is worrying and annoying as hell. So here I come up with a way to 
    optimise the lat, lon to get the LST I want. This is a bit of a hack, but
    it works. """
    
    # lat, lon = inputs
    
    lat, height = inputs
    date = "2000-01-01T12:00:00.0"
    
    ##Setup location
    observing_location = EarthLocation(lat=lat*u.deg, lon=lon*u.deg, height=height)
    ##Setup time at that location

    lst_type = 'mean'

    observing_time = Time(date, scale='utc', location=observing_location)
    ##Grab the LST
    LST = observing_time.sidereal_time(lst_type)
    LST_deg = LST.value*15.0
    
    lst_current = LST_deg
    
    return lst_current
    
def estimate_zero_lat_long():
    """Every time they update astropy, the exact values of LST seems to change
    This is worrying and annoying as hell. So here I come up with a way to 
    optimise the lat, lon to get the LST I want. This is a bit of a hack, but
    it works. """

    intitial_guess = [79.6423588359480874]
    
    desried_lst = [0.0]
    
    height=0.0
    lat=0.0
    
    inputs = [lat, height]
    
    popt, pcov = opt.curve_fit(get_lon, inputs, desried_lst, p0=intitial_guess)
    
    return popt[0]

##Vehicle for running tests
class Test(unittest.TestCase):
    
    @classmethod
    def tearDown(self):
        """Remove the uvfits files and log files that are created after all
        tests are run"""
        
        uvfits_files = glob("*.uvfits")
        log_files = glob("*.log")
        for uvfits_file in uvfits_files:
            os.remove(uvfits_file)
        for log_file in log_files:
            os.remove(log_file)        
    
    def check_uvfits_contents(self, uvfits_name):

        with fits.open(uvfits_name) as hdu:
            ant_table = hdu[1]
            data_table = hdu[0]
            
            # ##Test that a bunch of named columns in the data table are correct
            npt.assert_allclose(self.uu, data_table.data['UU']*speed_of_light.value, atol=1e-12)
            npt.assert_allclose(np.zeros(self.num_baselines), data_table.data['VV']*speed_of_light.value, atol=1e-12)
            npt.assert_allclose(np.zeros(self.num_baselines), data_table.data['WW']*speed_of_light.value, atol=1e-12)
            
            npt.assert_allclose(self.expec_data, data_table.data.data, atol=1e-12)
            
    def expec_values(self, autos=False, expec_gain=5.0):
        """Makes some expected outcomes """
        
        longitude = estimate_zero_lat_long()
        latitude = 0.0
        
        self.date = "2000-01-01T12:00:00.0"
        self.lowest = 180e+6
        self.telescope_name = "test"
        self.longitude = longitude
        self.latitude = latitude
        self.height = 0.0
        
        observing_location = EarthLocation(lat=self.latitude*u.deg,
                                           lon=self.longitude*u.deg,
                                           height=self.height)
        ##Setup time at that location

        lst_type = 'mean'

        observing_time = Time(self.date, scale='utc', location=observing_location)
        ##Grab the LST
        LST = observing_time.sidereal_time(lst_type)
        LST_deg = LST.value*15.0
        
        self.central_freq_chan = 8
        self.ch_width = 80e3
        self.ra_point = 0.0
        self.dec_point = 0.0
        self.num_antennas = 3
        if autos:
            self.num_baselines = int(((self.num_antennas - 1)*self.num_antennas) / 2) + self.num_antennas
        else:
            self.num_baselines = int(((self.num_antennas - 1)*self.num_antennas) / 2)
            
        self.num_time_steps = 1
        self.num_freq_channels = 16
        # self.jd_midnight = 2458646.5
        # self.gitlabel = 'as987has'
        self.XYZ_array = np.zeros(self.num_antennas*3)
        self.XYZ_array.shape = (self.num_antennas, 3)
        
        self.XYZ_array[:, 1] = np.array([0, 10, 20])
        # self.output_uvfits_name = "test_run_woden_band01.uvfits"
        
        if autos:
            self.uu = np.array([0, -10, -20, 0, -10, 0])
        else:
            self.uu = np.array([-10, -20, -10])
        
        
        
        expec_data = np.zeros((self.num_time_steps*self.num_baselines, 1, 1, self.num_freq_channels, 4, 3))
        
        expec_data[:,0,0,:,0,0] = expec_gain
        expec_data[:,0,0,:,1,0] = expec_gain
        expec_data[:,0,0,:,:,2] = 1.0
        
        self.expec_data = expec_data

    def test_nobeam_phase_centre_double(self):
        """Checks that concatenating four uvfits files that have `band{band}.uvfits`
        at the end works, by creating four uvfits with known properties"""
        
        self.expec_values()
        
        args = []
        args.append("--band_nums=1")
        args.append("--lowest_channel_freq=180e6")
        args.append("--coarse_band_width=1.28e+6")
        args.append("--freq_res=80e+3")
        args.append("--num_time_steps=1")
        args.append("--time_res=1e-16")
        args.append("--ra0=0.0")
        args.append("--dec0=0.0")
        args.append("--date=2000-01-01T12:00:00")
        args.append("--array_height=0.0")
        args.append(f"--longitude={self.longitude}")
        args.append(f"--latitude=0.0")
        args.append(f"--cat_filename={code_dir}/simple_sky.txt")
        args.append("--output_uvfits_prepend=test_run_woden")
        args.append("--primary_beam=none")
        args.append(f"--array_layout={code_dir}/simple_array.txt")
        args.append("--telescope_name=test")
        args.append("--no_precession")
        args.append("--num_threads=1")

        rw.main(args)
        
        self.check_uvfits_contents("test_run_woden_band01.uvfits")
        
        # ##do it again in CPU mode
        args.append("--cpu_mode")
        rw.main(args)
        self.check_uvfits_contents("test_run_woden_band01.uvfits")
        

    def test_nobeam_phase_centre_float(self):
        
        self.expec_values()
        
        args = []
        args.append("--band_nums=1")
        args.append("--lowest_channel_freq=180e6")
        args.append("--coarse_band_width=1.28e+6")
        args.append("--freq_res=80e+3")
        args.append("--num_time_steps=1")
        args.append("--time_res=1e-16")
        args.append("--ra0=0.0")
        args.append("--dec0=0.0")
        args.append("--date=2000-01-01T12:00:00")
        args.append("--array_height=0.0")
        args.append(f"--longitude={self.longitude}")
        args.append(f"--latitude=0.0")
        args.append(f"--cat_filename={code_dir}/simple_sky.txt")
        args.append("--output_uvfits_prepend=test_run_woden")
        args.append("--primary_beam=none")
        args.append(f"--array_layout={code_dir}/simple_array.txt")
        args.append("--telescope_name=test")
        args.append("--no_precession")
        args.append("--precision=float")
        args.append("--num_threads=1")

        rw.main(args)
        
        self.check_uvfits_contents("test_run_woden_band01.uvfits")
        
        ##do it again in CPU mode
        args.append("--cpu_mode")
        rw.main(args)
        self.check_uvfits_contents("test_run_woden_band01.uvfits")
    
    def test_nobeam_phase_centre_float_autos(self):

        self.expec_values(autos=True)
        
        args = []
        args.append("--band_nums=1")
        args.append("--lowest_channel_freq=180e6")
        args.append("--coarse_band_width=1.28e+6")
        args.append("--freq_res=80e+3")
        args.append("--num_time_steps=1")
        args.append("--time_res=1e-16")
        args.append("--ra0=0.0")
        args.append("--dec0=0.0")
        args.append("--date=2000-01-01T12:00:00")
        args.append("--array_height=0.0")
        args.append(f"--longitude={self.longitude}")
        args.append(f"--latitude=0.0")
        args.append(f"--cat_filename={code_dir}/simple_sky.txt")
        args.append("--output_uvfits_prepend=test_run_woden")
        args.append("--primary_beam=none")
        args.append(f"--array_layout={code_dir}/simple_array.txt")
        args.append("--telescope_name=test")
        args.append("--no_precession")
        args.append("--do_autos")
        args.append("--num_threads=1")

        rw.main(args)
        
        self.check_uvfits_contents("test_run_woden_band01.uvfits")
        
        ##do it again in CPU mode
        args.append("--cpu_mode")
        rw.main(args)
        self.check_uvfits_contents("test_run_woden_band01.uvfits")
        
        
    def test_nobeam_phase_centre_float_autos_multithread_multiband(self):

        self.expec_values(autos=True)
        
        args = []
        args.append("--band_nums=1,2,3")
        args.append("--lowest_channel_freq=180e6")
        args.append("--coarse_band_width=1.28e+6")
        args.append("--freq_res=80e+3")
        args.append("--num_time_steps=1")
        args.append("--time_res=1e-16")
        args.append("--ra0=0.0")
        args.append("--dec0=0.0")
        args.append("--date=2000-01-01T12:00:00")
        args.append("--array_height=0.0")
        args.append(f"--longitude={self.longitude}")
        args.append(f"--latitude=0.0")
        args.append(f"--cat_filename={code_dir}/simple_sky.txt")
        args.append("--output_uvfits_prepend=test_run_woden")
        args.append("--primary_beam=none")
        args.append(f"--array_layout={code_dir}/simple_array.txt")
        args.append("--telescope_name=test")
        args.append("--no_precession")
        args.append("--do_autos")
        args.append("--num_threads=3")

        rw.main(args)
        
        self.check_uvfits_contents("test_run_woden_band01.uvfits")
        self.check_uvfits_contents("test_run_woden_band02.uvfits")
        self.check_uvfits_contents("test_run_woden_band03.uvfits")
        
        ##do it again in CPU mode
        args.append("--cpu_mode")
        rw.main(args)
        self.check_uvfits_contents("test_run_woden_band01.uvfits")
        self.check_uvfits_contents("test_run_woden_band02.uvfits")
        self.check_uvfits_contents("test_run_woden_band03.uvfits")
        
        
    def test_runs_with_profiler_on(self):

        self.expec_values()
        
        args = []
        args.append("--band_nums=1")
        args.append("--lowest_channel_freq=180e6")
        args.append("--coarse_band_width=1.28e+6")
        args.append("--freq_res=80e+3")
        args.append("--num_time_steps=1")
        args.append("--time_res=1e-16")
        args.append("--ra0=0.0")
        args.append("--dec0=0.0")
        args.append("--date=2000-01-01T12:00:00")
        args.append("--array_height=0.0")
        args.append(f"--longitude={self.longitude}")
        args.append(f"--latitude=0.0")
        args.append(f"--cat_filename={code_dir}/simple_sky.txt")
        args.append("--output_uvfits_prepend=test_run_woden")
        args.append("--primary_beam=none")
        args.append(f"--array_layout={code_dir}/simple_array.txt")
        args.append("--telescope_name=test")
        args.append("--no_precession")
        args.append("--num_threads=1")
        args.append("--profile")

        rw.main(args)
        
        call("rm *.lprof", shell=True)
        self.check_uvfits_contents("test_run_woden_band01.uvfits")
        
    def test_verbose_logger(self):
        
        self.expec_values(autos=True)
        
        args = []
        args.append("--band_nums=1,2,3")
        args.append("--lowest_channel_freq=180e6")
        args.append("--coarse_band_width=1.28e+6")
        args.append("--freq_res=80e+3")
        args.append("--num_time_steps=1")
        args.append("--time_res=1e-16")
        args.append("--ra0=0.0")
        args.append("--dec0=0.0")
        args.append("--date=2000-01-01T12:00:00")
        args.append("--array_height=0.0")
        args.append(f"--longitude={self.longitude}")
        args.append(f"--latitude=0.0")
        args.append(f"--cat_filename={code_dir}/simple_sky.txt")
        args.append("--output_uvfits_prepend=test_run_woden")
        args.append("--primary_beam=none")
        args.append(f"--array_layout={code_dir}/simple_array.txt")
        args.append("--telescope_name=test")
        args.append("--no_precession")
        args.append("--do_autos")
        args.append("--num_threads=1")
        args.append("--save_log")
        args.append("--verbose")

        rw.main(args)
        self.check_uvfits_contents("test_run_woden_band01.uvfits")
        
        ##look for a log file to check it was saved
        log_files = glob("*test_run_woden*.log")
        
        self.assertTrue(len(log_files) > 0)
        
    def test_empty_model_thread(self):
        """Check if there are more threads than sources in sky things still work"""
        
        self.expec_values(autos=True)
        
        args = []
        args.append("--band_nums=1,2,3")
        args.append("--lowest_channel_freq=180e6")
        args.append("--coarse_band_width=1.28e+6")
        args.append("--freq_res=80e+3")
        args.append("--num_time_steps=1")
        args.append("--time_res=1e-16")
        args.append("--ra0=0.0")
        args.append("--dec0=0.0")
        args.append("--date=2000-01-01T12:00:00")
        args.append("--array_height=0.0")
        args.append(f"--longitude={self.longitude}")
        args.append(f"--latitude=0.0")
        args.append(f"--cat_filename={code_dir}/simple_sky.txt")
        args.append("--output_uvfits_prepend=test_run_woden")
        args.append("--primary_beam=none")
        args.append(f"--array_layout={code_dir}/simple_array.txt")
        args.append("--telescope_name=test")
        args.append("--no_precession")
        args.append("--do_autos")
        args.append("--num_threads=6")
        args.append("--cpu_mode")

        rw.main(args)
        self.check_uvfits_contents("test_run_woden_band01.uvfits")
        
        
    def test_runs_with_uvbeam_MWA(self):
        
        args = []
        args.append("--band_nums=1")
        args.append("--num_time_steps=1")
        args.append("--num_freq_channels=1")
        args.append(f"--metafits_filename={code_dir}/../../../examples/metafits/1125949472_metafits.fits")
        args.append(f"--cat_filename={code_dir}/simple_sky.txt")
        args.append("--output_uvfits_prepend=test_run_woden")
        args.append("--primary_beam=uvbeam_MWA")
        args.append(f"--array_layout={code_dir}/simple_array.txt")
        args.append("--num_threads=1")
        args.append("--cpu_mode")
        
        rw.main(args)
        
        ##we don't actually care what the outputs are here; this is just
        ##to check that it runs without error. All the myriad other tests
        ##get into the nitty gritty of the output values.
            
    def test_runs_multiple_GPU_rounds(self):
        """Want to check that multiple rounds of processing work correctly.
        This normally happens when there are more chunks than threads. Size
        of chunks is normally set my maximum number of visibilities that can
        be calculated. Instead, we set --max_sky_directions=2. As there are
        5 components in the sky model, and we set --num_threads=2, we get
        4 chunks in first round, and 1 in the second. """
        
        self.expec_values(expec_gain=65.0)
        
        args = []
        args.append("--band_nums=1")
        args.append("--lowest_channel_freq=180e6")
        args.append("--coarse_band_width=1.28e+6")
        args.append("--freq_res=80e+3")
        args.append("--num_time_steps=1")
        args.append("--time_res=1e-16")
        args.append("--ra0=0.0")
        args.append("--dec0=0.0")
        args.append("--date=2000-01-01T12:00:00")
        args.append("--array_height=0.0")
        args.append(f"--longitude={self.longitude}")
        args.append(f"--latitude=0.0")
        args.append(f"--cat_filename={code_dir}/simple_sky_two.txt")
        args.append("--output_uvfits_prepend=test_run_woden")
        args.append("--primary_beam=none")
        args.append(f"--array_layout={code_dir}/simple_array.txt")
        args.append("--telescope_name=test")
        args.append("--no_precession")
        args.append("--chunking_size=1")
        args.append("--num_threads=2")
        
        rw.main(args)
        self.check_uvfits_contents("test_run_woden_band01.uvfits")
        
        
    def test_runs_with_everybeam_moved_location(self):
        """Check things work if we run using the LOFAR everybeam model,
        but move the location to the MRO"""
        
        args = []
        
        args.append("--ra0=60.0")
        args.append("--dec0=-27.0")
        args.append("--latitude=-26.703319405555554")
        args.append("--longitude=116.67081523611111")
        args.append("--height=377.0")
        args.append("--num_freq_channels=1")
        args.append("--num_time_steps=1")
        args.append("--freq_res=80e+3")
        args.append("--time_res=8.0")
        args.append(f"--cat_filename={code_dir}/simple_sky.txt")
        args.append("--date=2016-01-09T13:11:19")
        args.append("--lowest_channel_freq=167.035e+6")
        args.append("--coarse_band_width=1.28e+6")
        args.append("--band_nums=2")
        args.append("--output_uvfits_prepend=woden_LOFAR-MRO_EoR1")
        args.append(f"--beam_ms_path={code_dir}/../../../test_installation/everybeam/LOFAR_HBA_MOCK.ms")
        args.append("--primary_beam=everybeam_LOFAR")
        args.append("--off_cardinal_dipoles")
        args.append("--move_array_to_latlon")
        args.append("--eb_point_to_phase")
        args.append("--num_threads=5")
        args.append("--cpu_mode")
        
        rw.main(args)
        
        # with fits.open("woden_LOFAR-MRO_EoR1_band02.uvfits") as hdu:
        #     ant_table = hdu[1]
        #     data_table = hdu[0]
            
            # print(data_table.data.data[0,0,0,0,:,:])
            
            # ##Test that a bunch of named columns in the data table are correct
            # npt.assert_allclose(self.uu, data_table.data['UU']*speed_of_light.value, atol=1e-12)
            # npt.assert_allclose(np.zeros(self.num_baselines), data_table.data['VV']*speed_of_light.value, atol=1e-12)
            # npt.assert_allclose(np.zeros(self.num_baselines), data_table.data['WW']*speed_of_light.value, atol=1e-12)
            
            # npt.assert_allclose(self.expec_data, data_table.data.data, atol=1e-12)
            
    def test_runs_with_uvbeam_HERA_FITS(self):
        """Check it runs with the HERA uvbeam FITS file. Run test at MRO so
        we can use same obsevational settings and sky model as other tests"""
        
        args = []
        args.append("--band_nums=1")
        args.append("--num_time_steps=1")
        args.append("--num_freq_channels=1")
        args.append("--latitude=-26.703319405555554")
        args.append("--longitude=116.67081523611111")
        args.append("--height=377.0")
        args.append(f"--metafits_filename={code_dir}/../../../examples/metafits/1125949472_metafits.fits")
        args.append(f"--cat_filename={code_dir}/simple_sky.txt")
        args.append("--output_uvfits_prepend=test_run_woden")
        args.append("--primary_beam=uvbeam_HERA")
        args.append(f"--uvbeam_file_path={code_dir}/../../wodenpy/primary_beam/NF_HERA_Dipole_small.fits")
        args.append(f"--array_layout={code_dir}/simple_array.txt")
        args.append("--num_threads=1")
        args.append("--cpu_mode")
        
        rw.main(args)
        
        
##Run the test
if __name__ == '__main__':
    unittest.main()
