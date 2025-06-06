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
import wodenpy
from wodenpy.use_libwoden.skymodel_structs import Components_Python
from astropy.coordinates import EarthLocation
from astropy import units as u
from astropy.time import Time
from wodenpy.primary_beam.use_everybeam import run_everybeam_over_threads, run_everybeam
from wodenpy.use_libwoden.use_libwoden import check_for_everybeam
import importlib_resources

##Location of this file; use it to find test measurement sets
code_dir = os.path.realpath(__file__)
code_dir = ('/').join(code_dir.split('/')[:-1])

##Vehicle for running tests
class Test(unittest.TestCase):
    """Vehicle for running tests"""
    
    def do_test(self):
        """First of all, run `run_everybeam` to get a set of expected values
        for given inputs. Then run `run_everybeam_over_threads` for a number
        of threads, and check the outputs match."""
        
        arr_latitude = np.radians(-26.703319405555554)
        arr_long = np.radians(116.67081523611111)
        date = "2024-07-21T20:13:00"
        ms_path = f'{code_dir}/../../../test_installation/everybeam/MWA-single-timeslot.ms'
        
        coeff_path=os.environ['MWA_FEE_HDF5']
        
        location = EarthLocation(lat=arr_latitude*u.rad, 
                                 lon=arr_long*u.rad)

        observing_time = Time(date, scale='utc', location=location)

        ##Grab the LST
        LST = observing_time.sidereal_time('mean')
        LST_deg = LST.value*15
        # print(f"LST: {LST_deg} deg")
        
        components = Components_Python()
        
        beam_ra0 = np.radians(LST_deg)
        beam_dec0 = arr_latitude
        
        ##Setup components ra/decs in a way where they are all above the horizon
        num_comps = 25
        ras = np.linspace(beam_ra0-np.pi/4, beam_ra0+np.pi/4, num_comps)
        ras[ras < 0] += np.pi*2
        decs = np.linspace(arr_latitude - np.pi/4, arr_latitude + np.pi/4, len(ras))
        decs[decs > np.pi/2] = np.pi - decs[decs > np.pi/2]
        
        ##No need to precess for this simple test
        
        j2000_latitudes = np.array([arr_latitude])
        j2000_lsts = np.array([np.radians(LST_deg)])
        
        all_times = np.array([observing_time])
        all_freqs = np.array([160e+6])
        station_ids = np.array([0])
        
        mwa_jones =  run_everybeam(ms_path, coeff_path,
                                ras, decs, beam_ra0, beam_dec0,
                                j2000_latitudes, j2000_lsts,
                                all_times, all_freqs,
                                station_ids,
                                apply_beam_norms=False,
                                iau_order=True,
                                element_only=False,
                                parallactic_rotate=True)
        
        
        ms_path = f'{code_dir}/../../../test_installation/everybeam/create_OSKAR-MWA_ms/OSKAR-MWA-layout.ms'
        oskar_jones =  run_everybeam(ms_path, coeff_path,
                                ras, decs, beam_ra0, beam_dec0,
                                j2000_latitudes, j2000_lsts,
                                all_times, all_freqs,
                                station_ids,
                                apply_beam_norms=True,
                                iau_order=True,
                                element_only=False,
                                parallactic_rotate=True)
        
        ##These are same array layout, but different dipoles (bow-tie vs.
        # analytic log dipole) so test abs value matches to within 15% is good enough
        npt.assert_allclose(np.abs(mwa_jones), np.abs(oskar_jones), atol=0.15)
        
            
    def test_run_OSKAR_MWA_vs_MWA(self):
        woden_lib_path = importlib_resources.files(wodenpy).joinpath(f"libwoden_double.so")
        have_everybeam = check_for_everybeam(woden_lib_path)
        if have_everybeam:
            self.do_test()
        else:
            print("Skipping test_run_everybeam_over_threads as everybeam not installed")

##Run the test
if __name__ == '__main__':
    unittest.main()