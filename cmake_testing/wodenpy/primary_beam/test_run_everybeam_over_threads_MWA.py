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
from wodenpy.primary_beam.use_everybeam import run_everybeam_over_threads, run_everybeam
from wodenpy.wodenpy_setup.run_setup import check_for_library

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
        all_freqs = np.array([100e+6])
        station_ids = np.array([0,10])
        
        # ##check the code runs when we manually add normalisation, even
        # ##though EveryBeam already seems to normalise the MWA beam
        # serial_jones =  run_everybeam(ms_path, coeff_path,
        #                         ras, decs, beam_ra0, beam_dec0,
        #                         j2000_latitudes, j2000_lsts,
        #                         all_times, all_freqs,
        #                         station_ids,
        #                         apply_beam_norms=True,
        #                         iau_order=True,
        #                         element_only=False,
        #                         parallactic_rotate=True)
        
        serial_jones =  run_everybeam(ms_path, coeff_path,
                                ras, decs, beam_ra0, beam_dec0,
                                j2000_latitudes, j2000_lsts,
                                all_times, all_freqs,
                                station_ids,
                                apply_beam_norms=False,
                                iau_order=True,
                                element_only=False,
                                parallactic_rotate=True)
        
        for num_threads in [1, 3, 4, 8]:
        
            parallel_jones =  run_everybeam_over_threads(num_threads,
                                ms_path, coeff_path,
                                ras, decs, beam_ra0, beam_dec0,
                                j2000_latitudes, j2000_lsts,
                                all_times, all_freqs,
                                station_ids,
                                apply_beam_norms=False,
                                iau_order=True,
                                element_only=False,
                                parallactic_rotate=True)
            
            npt.assert_allclose(serial_jones, parallel_jones, atol=1e-8)
            
    def test_run_everybeam_over_threads(self):
        have_everybeam = check_for_library('everybeam')
        if have_everybeam:
            self.do_test()
        else:
            print("Skipping test_run_everybeam_over_threads as everybeam not installed")

##Run the test
if __name__ == '__main__':
    unittest.main()