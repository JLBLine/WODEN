"""
Test `wodenpy.skymodel.read_fits_skymodel.calc_uvbeam_for_components`, which should
calculate uvbeam jones matrices for a set of components.
"""

from sys import path
import os
import unittest
import numpy as np
import numpy.testing as npt

from wodenpy.primary_beam.use_uvbeam import calc_uvbeam_for_components, setup_MWA_uvbeams
from wodenpy.use_libwoden.skymodel_structs import Components_Python
from astropy.coordinates import EarthLocation
from astropy import units as u
from astropy.time import Time
from wodenpy.array_layout.precession import RTS_Precess_LST_Lat_to_J2000
# from pyuvdata import ShortDipoleBeam, BeamInterface, UVBeam

##Location of this file; use it to find test measurement sets
code_dir = os.path.realpath(__file__)
code_dir = ('/').join(code_dir.split('/')[:-1])

FREQ = 99840000

##Vehicle for running tests
class Test(unittest.TestCase):
    """Vehicle for running tests"""
    
    def run_calc_uvbeam_for_components(self, uvbeam_objs, arr_latitude, arr_long,
                                          date, expec_gxs, expec_Dxs,
                                          expec_Dys, expec_gys):
        """Run calc_uvbeam_for_components for a minimal set of inputs,
        given the telescope, array location, date. Test outputs
        against expected values"""
        
        
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
        num_comps = 5
        ras = np.linspace(beam_ra0-np.pi/4, beam_ra0+np.pi/4, num_comps)
        ras[ras < 0] += np.pi*2
        decs = np.linspace(arr_latitude - np.pi/4, arr_latitude + np.pi/4, len(ras))
        decs[decs > np.pi/2] = np.pi - decs[decs > np.pi/2]
        
        ##No need to precess for this simple test
        
        j2000_latitudes = np.array([arr_latitude])
        j2000_lsts = np.array([np.radians(LST_deg)])
        
        all_times = np.array([observing_time])
        all_freqs = np.array([FREQ])
        station_id = 0
        
        components.ras = ras
        components.decs = decs
        
        calc_uvbeam_for_components(components, uvbeam_objs,
                                   all_freqs,
                                   j2000_latitudes, j2000_lsts)
        
        ##I've intentionally put the central coord at beam centre, so we should
        ##get gains of 1 and leakages of 0
        atol=4e-4
        npt.assert_allclose(np.abs(components.gxs[2]), 1, atol=atol)
        npt.assert_allclose(np.abs(components.Dxs[2]), 0, atol=atol)
        npt.assert_allclose(np.abs(components.Dys[2]), 0, atol=atol)
        npt.assert_allclose(np.abs(components.gys[2]), 1, atol=atol)
        
        ##Now just check against the expected values
        npt.assert_allclose(components.gxs, expec_gxs, rtol=1e-4, atol=atol)
        npt.assert_allclose(components.Dxs, expec_Dxs, rtol=1e-4, atol=atol)
        npt.assert_allclose(components.Dys, expec_Dys, rtol=1e-4, atol=atol)
        npt.assert_allclose(components.gys, expec_gys, rtol=1e-4, atol=atol)
        
    def test_MWA_zenith(self):
        """Run test using an MWA beam. Skip if we can't find the hdf5 file"""
        
        try:
            mwa_coeff = os.environ['MWA_FEE_HDF5']
        
            mwa_lat = np.radians(-26.703319405555554)
            mwa_long = np.radians(116.67081523611111)
            date = "2024-07-21T20:13:00"
            
            
            expec_gxs = [-2.11355046e-02+0.01884995j, 2.10453938e-01-0.28730457j,
                         6.32404804e-01-0.77463806j, 1.83363035e-01-0.24991899j,
                         3.43799517e-04-0.00361672j]
            expec_Dxs = [-0.01458456+8.32428025e-03j, 0.04058375-5.83309369e-02j,
                         -0.00030464-7.05280627e-06j, -0.03390884+4.28678427e-02j,
                         -0.0021708 -2.71366936e-04j]
            expec_Dys = [0.01672783-2.89442996e-02j, -0.0684786 +8.83901277e-02j,
                        -0.00030525+2.06055097e-06j,  0.00571018-1.14795095e-02j,
                        -0.00270409-9.71247947e-04j]
            expec_gys = [-0.01732507+0.02962904j, 0.2206394 -0.29093557j,
                         0.63274342-0.77436155j, 0.18302489-0.25105522j,
                         -0.00100178-0.00475868j]
            
            freqs = np.array([FREQ])
            delays = np.zeros(16, dtype=int)
            gains = np.ones(32, dtype=float)
            
            uvbeam_objs = setup_MWA_uvbeams(mwa_coeff, freqs,
                        delays, gains,
                        pixels_per_deg = 5)
            
            self.run_calc_uvbeam_for_components(uvbeam_objs, mwa_lat, mwa_long,
                                                date, expec_gxs, expec_Dxs,
                                                expec_Dys, expec_gys)
        except KeyError:
            print("MWA_FEE_HDF5 not set. Skipping test on MWA beam as cannot access FEE hdf5 file")
            
        
##Run the test
if __name__ == '__main__':
    
    unittest.main()