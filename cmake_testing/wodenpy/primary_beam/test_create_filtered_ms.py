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
from astropy.coordinates import EarthLocation
from astropy import units as u
from astropy.time import Time
from wodenpy.primary_beam.use_everybeam import run_everybeam_over_threads, create_filtered_ms
from wodenpy.use_libwoden.use_libwoden import check_for_everybeam
import importlib_resources
from erfa import ae2hd

##Location of this file; use it to find test measurement sets
code_dir = os.path.realpath(__file__)
code_dir = ('/').join(code_dir.split('/')[:-1])

MWA_LAT = -26.703319405555554
MWA_LONG = 116.67081523611111
MWA_HEIGHT = 377.827

LOFAR_LAT = 52.905329712
LOFAR_LONG = 6.867996528
LOFAR_HEIGHT = 0.0

##Vehicle for running tests
class Test(unittest.TestCase):
    """Vehicle for running tests"""
    
    def run_beam_given_inputs(self, arr_latitude, LST_deg, ms_path, observing_time):
        
        ##set the beam centre pointing to zenith
        beam_ra0 = np.radians(LST_deg)
        beam_dec0 = arr_latitude
        
        ##No need to precess for this simple test
        
        j2000_latitudes = np.array([arr_latitude])
        j2000_lsts = np.array([np.radians(LST_deg)])
        
        all_times = np.array([observing_time])
        all_freqs = np.array([100e+6])
        station_ids = np.array([0,10])
        
        
        ##go for az,za close to zenith
        
        az_range = np.arange(0, 2*np.pi, np.pi/4)
        za_range = np.arange(np.radians(1), np.radians(4), np.radians(1))
        azs, zas = np.meshgrid(az_range, za_range)
        azs = azs.flatten()
        zas = zas.flatten()
        
        ##convert to ra,dec as that's what the beam function uses
        has, decs = ae2hd(azs, np.pi/2 - zas, arr_latitude)
        ras = np.radians(LST_deg) - has
        
        ras[ras < 0] += 2*np.pi
        
        # print(np.degrees(ras))
        # print(np.degrees(decs))
        
        coeff_path = ""
        num_threads = 1
        jones =  run_everybeam_over_threads(num_threads, ms_path, coeff_path,
                                ras, decs, beam_ra0, beam_dec0,
                                j2000_latitudes, j2000_lsts,
                                all_times, all_freqs,
                                station_ids,
                                apply_beam_norms=False,
                                iau_order=True,
                                element_only=False,
                                parallactic_rotate=False)
        
        return jones
    
    def test_move_LOFAR_to_MWA(self):
        """Try moving the array centre to a different location. Test by running
        with a set of az,za at original location. Then try and move the array
        and compare the same az,za to the new location. """
        
        # date = "2024-07-21T03:35:00"
        date = "2014-07-01T21:33:04"
        
        arr_latitude = np.radians(LOFAR_LAT)
        arr_long = np.radians(LOFAR_LONG)
        orig_ms_path = f'{code_dir}/../../../test_installation/everybeam/LOFAR_HBA_MOCK.ms'
        
        location = EarthLocation(lat=arr_latitude*u.rad, 
                                 lon=arr_long*u.rad)

        observing_time = Time(date, scale='utc', location=location)

        ##Grab the LST
        LST = observing_time.sidereal_time('mean')
        LST_deg = LST.value*15
        
        
        ms_path = "test.ms"
        
        ##First, use the function just to change the pointing of the measurement set
        create_filtered_ms(orig_ms_path, ms_path,
                           np.radians(LST_deg), arr_latitude)        
        
        in_situ_jones = self.run_beam_given_inputs(arr_latitude, LST_deg, ms_path, observing_time)
        # print(np.abs(in_situ_jones[0,0,0,:,0,0]))
        # ##Next, use the function to move the array location as well
        
        arr_latitude = np.radians(MWA_LAT)
        arr_long = np.radians(MWA_LONG)
        
        orig_ms_path = f'{code_dir}/../../../test_installation/everybeam/LOFAR_HBA_MOCK.ms'
        
        location = EarthLocation(lat=arr_latitude*u.rad, 
                                 lon=arr_long*u.rad)

        observing_time = Time(date, scale='utc', location=location)

        ##Grab the LST
        LST = observing_time.sidereal_time('mean')
        LST_deg = LST.value*15
        
        recentre_array=True
        create_filtered_ms(orig_ms_path, ms_path,
                           np.radians(LST_deg), arr_latitude,
                           recentre_array,
                           np.radians(LOFAR_LAT), np.radians(LOFAR_LONG),
                           np.radians(MWA_LAT), np.radians(MWA_LONG))  
        
        recentred_jones = self.run_beam_given_inputs(arr_latitude, LST_deg, ms_path, observing_time)
        # print(np.abs(recentred_jones[0,0,0,:,0,0]))
        
        npt.assert_allclose(np.abs(in_situ_jones),
                            np.abs(recentred_jones),
                            atol=0.16)
        
        npt.assert_allclose(np.angle(in_situ_jones),
                            np.angle(recentred_jones),
                            atol=np.radians(10))
        
        
    def test_move_LOFAR_to_LOFAR(self):
        """Try moving a LOFAR array to what we think as the LOFAR array centre.
        Do this to check how well our method is working; if we can do the
        correct rotations, we should be able to get the same beam pattern."""
        
        date = "2014-07-01T21:33:04"
        
        arr_latitude = np.radians(LOFAR_LAT)
        arr_long = np.radians(LOFAR_LONG)
        orig_ms_path = f'{code_dir}/../../../test_installation/everybeam/LOFAR_HBA_MOCK.ms'
        
        location = EarthLocation(lat=arr_latitude*u.rad, 
                                 lon=arr_long*u.rad)

        observing_time = Time(date, scale='utc', location=location)

        ##Grab the LST
        LST = observing_time.sidereal_time('mean')
        LST_deg = LST.value*15
        
        ms_path = "test.ms"
        
        ##First, use the function just to change the pointing of the measurement set
        create_filtered_ms(orig_ms_path, ms_path,
                           np.radians(LST_deg), arr_latitude)        
        
        in_situ_jones = self.run_beam_given_inputs(arr_latitude, LST_deg, ms_path, observing_time)
        # print(np.abs(in_situ_jones[0,0,0,:,0,0]))
        
        # # ##Next, use the function to move the array location as well
        
        
        
        recentre_array=True
        create_filtered_ms(orig_ms_path, ms_path,
                           np.radians(LST_deg), arr_latitude,
                           recentre_array,
                           np.radians(LOFAR_LAT), np.radians(LOFAR_LONG),
                           np.radians(LOFAR_LAT), np.radians(LOFAR_LONG))  
        
        recentred_jones = self.run_beam_given_inputs(arr_latitude, LST_deg, ms_path, observing_time)
        # print(np.abs(in_situ_jones[0,0,0,:,0,0]))
        
        npt.assert_allclose(np.abs(in_situ_jones),
                            np.abs(recentred_jones),
                            atol=0.01)
        
        npt.assert_allclose(np.angle(in_situ_jones),
                            np.angle(recentred_jones),
                            atol=np.radians(1))
            

##Run the test
if __name__ == '__main__':
    unittest.main()