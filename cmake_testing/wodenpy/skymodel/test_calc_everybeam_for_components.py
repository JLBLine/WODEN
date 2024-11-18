"""
Test `wodenpy.skymodel.read_fits_skymodel.calc_everybeam_for_components`, which should
calculate EveryBeam jones matrices for a set of components.
"""

from sys import path
import os
import unittest
import numpy as np
import numpy.testing as npt

from wodenpy.skymodel.read_fits_skymodel import calc_everybeam_for_components
from wodenpy.use_libwoden.skymodel_structs import Components_Python
from wodenpy.primary_beam.use_everybeam import load_LOFAR_telescope, load_MWA_telescope, load_OSKAR_telescope
from astropy.coordinates import EarthLocation
from astropy import units as u
from astropy.time import Time
from wodenpy.array_layout.precession import RTS_Precess_LST_Lat_to_J2000

from wodenpy.wodenpy_setup.run_setup import check_for_library
have_everybeam = check_for_library('everybeam')

##Location of this file; use it to find test measurement sets
code_dir = os.path.realpath(__file__)
code_dir = ('/').join(code_dir.split('/')[:-1])

##Vehicle for running tests
class Test(unittest.TestCase):
    """Vehicle for running tests"""
    
    def run_calc_everybeam_for_components(self, telescope, arr_latitude, arr_long,
                                          date, expec_gxs, expec_Dxs,
                                          expec_Dys, expec_gys):
        """Run calc_everybeam_for_components for a minimal set of inputs,
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
        all_freqs = np.array([100e+6])
        station_id = 0
        
        components.ras = ras
        components.decs = decs
        
        calc_everybeam_for_components(beam_ra0, beam_dec0, num_comps,
                                components, telescope,
                                all_times, all_freqs,
                                j2000_latitudes, j2000_lsts,
                                arr_latitude, arr_long,
                                station_id=station_id)
        
        ##I've intention put the central coord at beam centre, so we should
        ##get gains of 1 and leakages of 0
        atol=4e-4
        npt.assert_allclose(np.abs(components.gxs[2]), 1, atol=atol)
        npt.assert_allclose(np.abs(components.Dxs[2]), 0, atol=atol)
        npt.assert_allclose(np.abs(components.Dys[2]), 0, atol=atol)
        npt.assert_allclose(np.abs(components.gys[2]), 1, atol=atol)
        
        ##Now just check against the expected values
        npt.assert_allclose(components.gxs, expec_gxs, rtol=1e-5, atol=1e-5)
        npt.assert_allclose(components.Dxs, expec_Dxs, rtol=1e-5, atol=1e-5)
        npt.assert_allclose(components.Dys, expec_Dys, rtol=1e-5, atol=1e-5)
        npt.assert_allclose(components.gys, expec_gys, rtol=1e-5, atol=1e-5)
        
    
    def do_LOFAR_beam(self):
        """Run test using a LOFAR beam"""
        
        lofar_lat = np.radians(52.905329712)
        lofar_long = np.radians(6.867996528)
        date = "2024-07-21T03:35:00"
        ms_path = f'{code_dir}/../../../test_installation/everybeam/lba.MS'
        telescope = load_LOFAR_telescope(ms_path)
        
        expec_gxs = [0.05270926+0.16701497j, -0.09573493-0.13976068j,  1.+0.j,
                     0.04727046+0.01778293j, -0.03659679-0.01156235j]
        expec_Dxs = [-1.01143394e-01-0.13900946j, 5.24021239e-02+0.06176052j,
                     -1.11022302e-16+0.j, 1.40412645e-02+0.00620097j,
                     -2.80333957e-02-0.01121601j]
        expec_Dys = [-0.04592168+3.29134486e-02j, -0.00383622-1.83244478e-02j,
                     0.0+2.71050543e-20j, -0.02328202-5.89294827e-03j,
                    0.04240413+5.82550985e-03j]
        expec_gys = [0.07605183+0.20740728j, -0.10956174-0.1518494j, 1.+0.j,
                     0.05155368+0.0144725j, -0.03881918-0.00673202j]
        
        self.run_calc_everybeam_for_components(telescope, lofar_lat, lofar_long,
                                               date, expec_gxs, expec_Dxs,
                                               expec_Dys, expec_gys)

    
    def do_MWA_beam(self):
        """Run test using a LOFAR beam. Skip if we can't find the hdf5 file"""
        
        try:
            mwa_coeff = os.environ['MWA_FEE_HDF5']
        
            mwa_lat = np.radians(-26.703319405555554)
            mwa_long = np.radians(116.67081523611111)
            date = "2024-07-21T20:13:00"
            
            ms_path = f'{code_dir}/../../../test_installation/everybeam/MWA-single-timeslot.ms'
            telescope = load_MWA_telescope(ms_path, coeff_path=mwa_coeff)
            
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
            
            self.run_calc_everybeam_for_components(telescope, mwa_lat, mwa_long,
                                                date, expec_gxs, expec_Dxs,
                                                expec_Dys, expec_gys)
        except KeyError:
            print("MWA_FEE_HDF5 not set. Skipping test on MWA beam as cannot access FEE hdf5 file")
            
    def do_OSKAR_beam(self):
        """Run test using a OSKAR beam"""
        
        mwa_lat = np.radians(-26.703319405555554)
        mwa_long = np.radians(116.67081523611111)
        date = "2024-07-21T20:13:00"
        ms_path = f'{code_dir}/OSKAR-MWA-layout.ms'
        telescope = load_OSKAR_telescope(ms_path)
        
        expec_gxs = [-2.43335832e-02+2.31609975e-02j,  3.34658510e-01-9.41421237e-02j,
                     1.00000000e+00-3.46096504e-19j,  2.91801131e-01-9.99791100e-02j,
                     8.64411217e-04-1.25713046e-03j]
        expec_Dxs = [-0.01155318+0.01050554j,  0.06288374-0.01461055j,
                     0.+0.j, -0.05672406+0.02249763j, -0.00048496+0.00078277j]
        expec_Dys = [0.03609715-3.35324770e-02j, -0.11445022+3.51491847e-02j,
                     0.-1.35525272e-20j,  0.00361486+4.41763924e-03j,
                     -0.00033303+7.26006209e-04j]
        expec_gys = [-3.56786286e-02+3.06814048e-02j, 3.55409668e-01-9.18447232e-02j,
                     1.00000000e+00-3.46097331e-19j, 2.90453908e-01-9.59879126e-02j,
                     7.07749933e-04-1.09958296e-03j]
        
        self.run_calc_everybeam_for_components(telescope, mwa_lat, mwa_long,
                                                date, expec_gxs, expec_Dxs,
                                                expec_Dys, expec_gys)
        
    def test_all_beams(self):
        
        if have_everybeam:
            self.do_LOFAR_beam()
            self.do_MWA_beam()
            self.do_OSKAR_beam()
            
        else:
            print("Skipping test_calc_everybeam_for_components as everybeam not installed")

##Run the test
if __name__ == '__main__':
    
    unittest.main()