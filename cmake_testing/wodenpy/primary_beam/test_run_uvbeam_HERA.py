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
from wodenpy.primary_beam.use_uvbeam import run_uvbeam, setup_HERA_uvbeams
import importlib_resources
import mwa_hyperbeam
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import TimeDelta
import erfa
from pyuvdata import ShortDipoleBeam, BeamInterface, UVBeam
from pyuvdata.telescopes import known_telescope_location
import matplotlib.pyplot as plt

location = known_telescope_location("HERA")
# import matplotlib.pyplot as plt

##Location of this file; use it to find the beam files
code_dir = os.path.realpath(__file__)
code_dir = ('/').join(code_dir.split('/')[:-1])


D2R = np.pi / 180.0
VELC = 299792458.0
HERA_LAT = location.lat.value
HERA_LAT_RAD = HERA_LAT * D2R
HERA_LONG = location.lon.value
HERA_LONG_RAD = HERA_LONG * D2R
HERA_HEIGHT = location.height.value

NSIDE=51

def plot_beam_response(all_jones, station_ind, time_ind, freq_ind, nside):
    
    fig, axs = plt.subplots(2, 2, figsize=(10, 10), layout='constrained')
    
    gx = all_jones[station_ind, time_ind, freq_ind, :, 0, 0]
    Dx = all_jones[station_ind, time_ind, freq_ind, :, 0, 1]
    Dy = all_jones[station_ind, time_ind, freq_ind, :, 1, 0]
    gy = all_jones[station_ind, time_ind, freq_ind, :, 1, 1]
    
    gx.shape = (nside, nside)
    Dx.shape = (nside, nside)
    Dy.shape = (nside, nside)
    gy.shape = (nside, nside)
    
    im = axs[0, 0].imshow(np.abs(gx), origin='lower')
    axs[0, 0].set_title('gx')
    plt.colorbar(im, ax=axs[0, 0])
    im = axs[0, 1].imshow(np.abs(Dx), origin='lower')
    axs[0, 1].set_title('Dx')
    plt.colorbar(im, ax=axs[0, 1])
    im = axs[1, 0].imshow(np.abs(Dy), origin='lower')
    axs[1, 0].set_title('Dy')
    plt.colorbar(im, ax=axs[1, 0])
    im = axs[1, 1].imshow(np.abs(gy), origin='lower')
    axs[1, 1].set_title('gy')
    plt.colorbar(im, ax=axs[1, 1])
    
    
    for ax in axs.flatten():
        ax.set_xticks([])
        ax.set_yticks([])
    
    fig.suptitle(f'HERA beam {station_ind} Time {time_ind} Frequency {freq_ind}')
    
    fig.savefig(f'hera_beam_station{station_ind:03d}_time{time_ind:03d}_freq{freq_ind:03d}.png', bbox_inches='tight')
    # fig.savefig(f'hera_beam_station{station_ind:03d}_time{time_ind:03d}_freq{freq_ind:03d}_no-rotation.png', bbox_inches='tight')
    plt.close()

##Vehicle for running tests
class Test(unittest.TestCase):
    """Vehicle for running tests"""
    
    # @classmethod
    # def setUpClass(cls):
    #     """Setup the """
    
    def do_HERA_beam_test(self, filenames, freqs):
        """First of all, run `run_everybeam` to get a set of expected values
        for given inputs. Then run `run_everybeam_over_threads` for a number
        of threads, and check the outputs match."""
        
        num_freqs = len(freqs)
        num_times = 2
        
        num_beams = 1
        
        date = "2013-05-16T10:48:48"
        
        # location = EarthLocation(lat=HERA_LAT_RAD*u.rad, 
        #                          lon=HERA_LONG_RAD*u.rad)

        observing_time = Time(date, scale='utc', location=location)

        time_inc = 60*60*3.0 ##3 hours

        times = [Time(date, scale='utc', location=location) + TimeDelta(time_inc*i, format='sec') for i in range(num_times)]
        lsts = [time.sidereal_time('mean').rad for time in times]
        latitudes = [HERA_LAT_RAD]*num_times ##no need to precess for this simple test
        
        ##Setup a dummy FITS header with appropriate settings
        header = fits.Header()

        ##Give it 301 pixel for each axis
        nside = NSIDE

        ##This resolution seems to cover the full sky nicely
        cpix = int(nside // 2) + 1
        cdelt = 0.35
        cdelt = 120 / nside
        # cdelt = 80 / nside

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
        header['CRVAL2']  = HERA_LAT
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
        
        uvbeam_objs = setup_HERA_uvbeams(filenames, freqs)
        
        uvbeam_jones = run_uvbeam(uvbeam_objs, ras, decs,
                                  latitudes, lsts,
                                  freqs, iau_order=True,
                                  parallactic_rotate=True)
        
        station_ind = 0
        ##super minimal test. Just check that the maximum values is close to 1,
        ##and sits at zenith
        for time_ind in range(num_times):
            for freq_ind in range(num_freqs):
                # plot_beam_response(uvbeam_jones, station_ind, time_ind, freq_ind, nside)
                
                lst = lsts[time_ind]
                
                x, y = wcs.all_world2pix(np.degrees(lst), HERA_LAT, 0.0)
                ##central pixel after flattening a 2D image
                loc_location = np.round(y)*nside + np.round(x)
                
                gx = np.abs(uvbeam_jones[station_ind, time_ind, freq_ind, :, 0, 0])
                gy = np.abs(uvbeam_jones[station_ind, time_ind, freq_ind, :, 1, 1])
                
                for g in [gx, gy]:
                    max_loc = np.nanargmax(g)
                    max_val = np.nanmax(g)
                    
                    self.assertAlmostEqual(max_val, 1.0, delta=0.05)
                    self.assertAlmostEqual(max_loc, loc_location, delta=1)
                
                
    def test_bad_filepaths(self):
        """Test things error when you give bad file paths"""
        
        freqs = np.array([100e+6, 125e+5])
        filenames = ['hoogyboogy.txt', 'hoogyboogy2.txt']
        
        with self.assertRaises(FileNotFoundError) as cm:
            self.do_HERA_beam_test(filenames, freqs)
            
    def test_unequal_files_and_freqs(self):
        """Test things error when you uneven number of files and freqs"""
        
        freqs = np.array([100e+6])
        filenames = ['hoogyboogy.txt', 'hoogyboogy2.txt']
        
        with self.assertRaises(ValueError) as cm:
            self.do_HERA_beam_test(filenames, freqs)
        
    def test_make_beams(self):
        """Test it works if everything is correct"""
        
        freqs = np.array([100e+6, 125e+6, 150e+6])
        filenames = [f'{code_dir}/HERA_4.9m_E-pattern_{int(freq/1e+6)}MHz.txt' for freq in freqs]
        interp_freqs = np.array([1050e+6, 120e+6, 137e+6])
        
        self.do_HERA_beam_test(filenames, interp_freqs)
        
    def test_make_beams_less_than_three(self):
        """Test that things run OK if we give less than three frequencies,
        but we get a warning"""
        
        freqs = np.array([100e+6, 125e+6])
        filenames = [f'{code_dir}/HERA_4.9m_E-pattern_{int(freq/1e+6)}MHz.txt' for freq in freqs]
        interp_freqs = np.array([105e+6, 120e+6])
        
        ##Use the assertLogs context manager to capture log messages
        with self.assertLogs() as cm:
        
            ##This is the message we expect to see
            msg = "Number of CST files submitted to HERA pyvudata UVBeam is less than 3.\n"
            msg += "WODEN will proceed, but will switch off frequency interpolation.\n"
            msg += "UVBeam needs at least three frequencies to perform default interpolation.\n"
            msg += "Further warnings of this type will be suppressed.\n"
            
            self.do_HERA_beam_test(filenames, interp_freqs)
            
            ##should only have one logging message
            self.assertEqual(len(cm.output), 1)
            
            ##The output message will have 'WARNING:wodenpy.wodenpy_setup.woden_logger:' at the start
            ##so we need to remove that
            output_msg = cm.output[0].split('WARNING:wodenpy.wodenpy_setup.woden_logger:')[1]
            self.assertEqual(output_msg, msg)
        
            ##If we run the same command again, we should not see the message
            ##repeated
            self.do_HERA_beam_test(filenames, interp_freqs)
            self.assertEqual(len(cm.output), 1)
            output_msg = cm.output[0].split('WARNING:wodenpy.wodenpy_setup.woden_logger:')[1]
            self.assertEqual(output_msg, msg)
            

##Run the test
if __name__ == '__main__':
    unittest.main()