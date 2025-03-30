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
from wodenpy.primary_beam.use_uvbeam import run_uvbeam
import importlib_resources
import mwa_hyperbeam
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import TimeDelta
import erfa
from pyuvdata import ShortDipoleBeam, BeamInterface, UVBeam
import matplotlib.pyplot as plt

##Location of this file; use it to find test measurement sets
code_dir = os.path.realpath(__file__)
code_dir = ('/').join(code_dir.split('/')[:-1])


D2R = np.pi / 180.0
VELC = 299792458.0
MWA_LAT = -26.703319405555554
MWA_LAT_RAD = MWA_LAT * D2R
MWA_LONG = 116.67081523611111
MWA_LONG_RAD = MWA_LONG * D2R
MWA_HEIGHT = 377.83

NSIDE=41

##Vehicle for running tests
class Test(unittest.TestCase):
    """Vehicle for running tests"""
    
    def do_MWA_beam_test(self, both_delays):
        """First of all, run `run_everybeam` to get a set of expected values
        for given inputs. Then run `run_everybeam_over_threads` for a number
        of threads, and check the outputs match."""
        
        
        # freqs = np.array([1.6768e+08, 1.8944e+08])
        freqs = np.array([1.87e+08, 1.8944e+08])
        num_freqs = len(freqs)
        num_times = 2
        num_beams = 1
        
        date = "2013-05-16T10:48:48"
        
        location = EarthLocation(lat=MWA_LAT_RAD*u.rad, 
                                 lon=MWA_LONG_RAD*u.rad)

        observing_time = Time(date, scale='utc', location=location)

        time_inc = 3600

        times = [Time(date, scale='utc', location=location) + TimeDelta(time_inc*i, format='sec') for i in range(num_times)]
        lsts = [time.sidereal_time('mean').rad for time in times]
        latitudes = [MWA_LAT_RAD]*num_times ##no need to precess for this simple test
        
        ##Setup a dummy FITS header with appropriate settings
        header = fits.Header()

        ##Give it 301 pixel for each axis
        nside = NSIDE

        ##This resolution seems to cover the full sky nicely
        cpix = int(nside // 2) + 1
        cdelt = 0.35
        cdelt = 80 / nside

        header['NAXIS']   = 2
        header['NAXIS1']  = nside
        header['NAXIS2']  = nside
        header['CTYPE1']  = 'RA---SIN'
        header['CRPIX1']  = cpix
        header['CRVAL1']  = np.degrees(lsts[0])
        header['CDELT1']  = cdelt
        header['CUNIT1']  = 'deg     '
        header['CTYPE2']  = 'DEC--SIN'
        header['CRPIX2']  = cpix
        header['CRVAL2']  = MWA_LAT
        header['CDELT2']  = cdelt
        header['CUNIT2']  = 'deg     '

        ##Make a world coord system
        wcs = WCS(header)

        ##Set up x/y pixels that cover the whole image
        x_mesh, y_mesh = np.meshgrid(np.arange(nside), np.arange(nside))

        x_pixels = x_mesh.flatten()
        y_pixels = y_mesh.flatten()

        ##convert to ra, dec
        ras, decs = wcs.all_pix2world(x_pixels, y_pixels, 0.0)
        ras = np.radians(ras)
        decs = np.radians(decs)
        
        coeff_path=os.environ['MWA_FEE_HDF5']
        hyp_beam = mwa_hyperbeam.FEEBeam(hdf5_file=coeff_path)
        
        norm_to_zenith = True
        iau_order = True
        
        gains = [1]*16
        
        full_hyp_jones = np.zeros((num_beams, num_times, num_freqs, nside*nside, 2, 2), dtype=complex)
        
        beam_ind = 0
        for time_ind, time in enumerate(times):
            for freq_ind, freq in enumerate(freqs):
                
                lst = lsts[time_ind]
                
                ##Then use erfa to convert these values into azs, els
                has = lst - ras

                ##use this erfa function to convert to azimuth and elevation
                ##were using degrees, but erfa uses rads, so convert here
                azs, els = erfa.hd2ae(has, decs, MWA_LAT_RAD)
                
                zas = np.pi/2 - els
                
                hyp_jones = hyp_beam.calc_jones_array(azs, zas, freq,
                                                  both_delays[0], gains, norm_to_zenith,
                                                  MWA_LAT_RAD, iau_order)
                
                
                # hyp_jones = hyp_beam.calc_jones_array(azs, zas, freq,
                #                                   both_delays[0], gains, norm_to_zenith)
                
                full_hyp_jones[beam_ind, time_ind, freq_ind, :, 0, 0] = hyp_jones[:, 0]
                full_hyp_jones[beam_ind, time_ind, freq_ind, :, 0, 1] = hyp_jones[:, 1]
                full_hyp_jones[beam_ind, time_ind, freq_ind, :, 1, 0] = hyp_jones[:, 2]
                full_hyp_jones[beam_ind, time_ind, freq_ind, :, 1, 1] = hyp_jones[:, 3]
          
        uvbeam_mwa = UVBeam.from_file(coeff_path, pixels_per_deg=5,
                                      delays=both_delays,
                                      freq_range=[freqs.min()-2*1.28e6, freqs.max()+2*1.28e+6])
        
        # print(uvbeam_mwa.bandpass_array)
        
        uvbeam_mwa.peak_normalize()
                
        uvbeam_obs = [uvbeam_mwa]
        
        uvbeam_jones = run_uvbeam(uvbeam_obs, ras, decs,
                                  latitudes, lsts,
                                  MWA_LAT_RAD, MWA_LONG_RAD,
                                  times, freqs, iau_order=True,
                                  parallactic_rotate=True)
        
        fig, axs = plt.subplots(4, 3, figsize=(8, 10), layout='constrained')
        
        this_hyp = np.abs(full_hyp_jones)
        this_uvbeam = np.abs(uvbeam_jones)
        
        for col, jones in enumerate([this_hyp, this_uvbeam, this_hyp - this_uvbeam]):
        # for col, jones in zip(range(1), [full_hyp_jones]):
            for row in range(4):
                
                beam = jones[0, 0, 0, :, row//2, row%2]
                beam.shape = (nside, nside)
                
                im = axs[row, col].imshow(beam, origin='lower')
                
                axs[row, col].set_xticks([])
                axs[row, col].set_yticks([])
                
                plt.colorbar(im, ax=axs[row, col])
                
        axs[0, 0].set_title('hyperbeam')
        axs[0, 1].set_title('pyuvdata')
        axs[0, 2].set_title('hyperbeam - pyuvdata')
                
        fig.savefig('test_uvbeam.png',bbox_inches='tight')
        plt.close(fig)
        
        npt.assert_allclose(uvbeam_jones, full_hyp_jones,
                            rtol=1e-1, atol=1e-4)
            
            
    
    def test_zenith(self):
        """Test the zenith beam"""
        
        # Set up the test parameters
        both_delays = np.zeros((2, 16), dtype='int')
        
        self.do_MWA_beam_test(both_delays)
        
        
    # def test_offzenith(self):
    #     """Test the zenith beam"""
        
    #     # Set up the test parameters
    #     both_delays = np.zeros((2, 16), dtype='int')
        
    #     both_delays[0, :] = np.tile(np.arange(0,8,2), 4)
    #     both_delays[1, :] = np.tile(np.arange(0,8,2), 4)
        
    #     self.do_MWA_beam_test(both_delays)

##Run the test
if __name__ == '__main__':
    unittest.main()