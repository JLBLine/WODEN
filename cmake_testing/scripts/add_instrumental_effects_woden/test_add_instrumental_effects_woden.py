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

import importlib_resources
import wodenpy

code_dir = os.path.realpath(__file__)
code_dir = ('/').join(code_dir.split('/')[:-1])

##Append the location of add_woden_uvfits to the sys.path to import from it
path.append('{:s}/../../../scripts'.format(code_dir))

##Code we are testing
from wodenpy.uvfits import wodenpy_uvfits
import add_instrumental_effects_woden as aiew
import numpy.testing as npt

##Vehicle for running tests
class Test(unittest.TestCase):

    def create_uvfits_outputs(self, num_freqs=3, do_autos=False,
                              num_ants=10, ch_width=40e3):
        """Tests the `wodenpy_uvfits.create_uvfits` function, which should take a whole
        heap of inputs and write out a uvits. Test by running funciton,
        reading in the created file and checking contents"""
        
        self.date = "2019-06-12T13:04:12"
        self.gst0_deg = 260.0303917560829
        self.degpdy = 360.98563614850775
        self.ut1utc = -0.17421478816666666
        self.num_antennas = num_ants
        self.freq_cent = 160e+6
        self.telescope_name = "MWA"
        self.longitude = 116.670813889
        self.latitude = -26.7033194444
        self.array_height = 377.0
        self.central_freq_chan = 2
        self.ch_width = ch_width
        self.ra_point = 0.0
        self.dec_point = -26.7
        
        if do_autos:
            self.num_baselines = int(((self.num_antennas - 1)*self.num_antennas) / 2) + self.num_antennas
        else:
            self.num_baselines = int(((self.num_antennas - 1)*self.num_antennas) / 2)
            
        self.num_time_steps = 2
        self.int_jd = 2458647.0
        self.gitlabel = 'as987has'
        self.XYZ_array = np.arange(self.num_antennas*3)
        self.XYZ_array.shape = (self.num_antennas, 3)
        self.output_uvfits_name = "unittest_example1_band01.uvfits"

        self.num_freq_channels = num_freqs
        self.num_freqs = num_freqs

        ##Some input test params
        # self.make_dummy_intputs()

        ##Create an antenna table
        ant_table = wodenpy_uvfits.make_antenna_table(XYZ_array=self.XYZ_array,
                      telescope_name=self.telescope_name,
                      num_antennas=self.num_antennas, freq_cent=self.freq_cent,
                      date=self.date, gst0_deg=self.gst0_deg,
                      degpdy=self.degpdy, ut1utc=self.ut1utc,
                      longitude=self.longitude, latitude=self.latitude,
                      array_height=self.array_height)

        ##Create some dummy data inputs
        v_container = np.ones((self.num_time_steps*self.num_baselines, 1, 1, self.num_freq_channels, 4, 3))
        
        v_container[:, :, :, :, :, 1] = 0.0

        uu = np.arange(self.num_baselines*self.num_time_steps)
        vv = np.arange(self.num_baselines*self.num_time_steps) + self.num_baselines
        ww = np.arange(self.num_baselines*self.num_time_steps) + 2*self.num_baselines
        
        self.time_res = 1.0
        baselines_array, date_array = wodenpy_uvfits.make_baseline_date_arrays(self.num_antennas,
                                self.date,
                                self.num_time_steps, self.time_res,
                                do_autos=do_autos)

        wodenpy_uvfits.create_uvfits(freq_cent=self.freq_cent, ch_width=self.ch_width,
                  central_freq_chan=self.central_freq_chan,
                  ra_point=self.ra_point, dec_point=self.dec_point,
                  output_uvfits_name=self.output_uvfits_name,
                  int_jd=self.int_jd, gitlabel=self.gitlabel,
                  v_container=v_container,
                  uu=uu, vv=vv, ww=ww,
                  baselines_array=baselines_array, date_array=date_array,
                  hdu_ant=ant_table, telescope_name=self.telescope_name,
                  longitude=self.longitude, latitude=self.latitude,
                  array_height=self.array_height)
        
        self.do_autos = do_autos
        
    def get_expec_visi(self, gxs, Dxs, Dys, gys, ant1, ant2):
        """Apply the gains to the visibilities for a specific baseline"""
        gx1 = gxs[ant1]
        Dx1 = Dxs[ant1]
        Dy1 = Dys[ant1]
        gy1 = gys[ant1]
        gx2 = np.conjugate(gxs[ant2])
        Dx2 = np.conjugate(Dxs[ant2])
        Dy2 = np.conjugate(Dys[ant2])
        gy2 = np.conjugate(gys[ant2])
        
        xx = gx1*(gx2 + Dy2) + Dx1*(gx2 + Dy2)
        xy = gx1*(Dx2 + gy2) + Dx1*(Dx2 + gy2)
        yx = Dy1*(gx2 + Dy2) + gy1*(gx2 + Dy2)
        yy = Dy1*(Dx2 + gy2) + gy1*(Dx2 + gy2)
        
        return xx, xy, yx, yy
        
    def add_visi(self, visi_ind, xx, xy, yx, yy):
        """put the visibilities into the expected visibilities array,
           assuming that the gains don't change with time"""
        
        for time in range(self.num_time_steps):
        
            self.expected_visis[visi_ind + time*self.num_baselines, :, 0] = xx
            self.expected_visis[visi_ind + time*self.num_baselines, :, 1] = yy
            self.expected_visis[visi_ind + time*self.num_baselines, :, 2] = xy
            self.expected_visis[visi_ind + time*self.num_baselines, :, 3] = yx
            
    def expected_visis(self, gxs, Dxs, Dys, gys):
        """Given the input params and expected gains, calculate the expected
        visibilities, assuming the sky model is complex(1,0) everywhere.
        Do it explicitly for each baseline and time step to check that the
        matrix multiplcation in the script is correct
        
        Gains/leakages should have shape (num_antennas, num_freqs)"""
        
        self.expected_visis = np.ones((self.num_baselines*self.num_time_steps, self.num_freq_channels, 4), dtype=complex)
        
        
        visi_ind = 0
        for ant1 in range(self.num_antennas):
            for ant2 in range(ant1, self.num_antennas):
                if ant1 == ant2:
                    if self.do_autos:
                        xx, xy, yx, yy = self.get_expec_visi(gxs, Dxs, Dys, gys, ant1, ant2)
                        self.add_visi(visi_ind, xx, xy, yx, yy)
                        visi_ind += 1
                
                else:
                    xx, xy, yx, yy = self.get_expec_visi(gxs, Dxs, Dys, gys, ant1, ant2)
                    self.add_visi(visi_ind, xx, xy, yx, yy)
                    visi_ind += 1
                    
    def run_code_test_ant_gains(self, args, expec_gx, expec_Dx, expec_Dy, expec_gy):
        """Runs the main script with the given args, and checks the output
        aingst the expected values"""
        
        ##Actually run the code!!
        aiew.main(args)
        ##read in the resuts
        uvfits = aiew.UVFITS('instrumental.uvfits')
        
        self.expected_visis(expec_gx, expec_Dx, expec_Dy, expec_gy)
        
        # diffs = np.where(np.abs(uvfits.visibilities - self.expected_visis) > 1e-6)
        ##check everything matches what we expect
        npt.assert_allclose(uvfits.visibilities, self.expected_visis, atol=1e-10)
        
        data = np.load("gains_applied_woden.npz")
        gx = data['gx']
        Dx = data['Dx']
        Dy = data['Dy']
        gy = data['gy']
        
        npt.assert_allclose(gx, expec_gx, atol=1e-6)
        npt.assert_allclose(Dx, expec_Dx, atol=1e-6)
        npt.assert_allclose(Dy, expec_Dy, atol=1e-6)
        npt.assert_allclose(gy, expec_gy, atol=1e-6)
        
    def test_add_gain_amps(self):
        """"""
        
        self.create_uvfits_outputs(do_autos=False)
        
        amp_err = 0.05
        seed = 983745
        
        args = []
        args.append(f"--numpy_seed={seed}")
        args.append("--uvfits=unittest_example1_band01.uvfits")
        args.append(f"--ant_gain_amp_error={amp_err}")
        
        ##Reset the random seed, and make what should have been created
        np.random.seed(seed)
        
        expec_gx = 1 + np.random.uniform(-amp_err, amp_err, self.num_antennas) + 1j*np.zeros(self.num_antennas)
        expec_gx = np.repeat(expec_gx, self.num_freqs)
        expec_gx.shape = (self.num_antennas, self.num_freqs)
        
        ##the main script calls random to generate phases, so we need to do that here too
        for ant in range(self.num_antennas):
            phase = np.random.uniform(-0.0*(np.pi/180.0), 0.0*(np.pi/180.0))
        
        expec_gy = 1 + np.random.uniform(-amp_err, amp_err, self.num_antennas) + 1j*np.zeros(self.num_antennas)
        expec_gy = np.repeat(expec_gy, self.num_freqs)
        expec_gy.shape = (self.num_antennas, self.num_freqs)
        
        expec_gx[0] = complex(1,0)
        expec_gy[0] = complex(1,0)
        
        expec_Dx = np.zeros((self.num_antennas, self.num_freqs), dtype=complex)
        expec_Dy = np.zeros((self.num_antennas, self.num_freqs), dtype=complex)
        
        
        self.run_code_test_ant_gains(args, expec_gx, expec_Dx, expec_Dy, expec_gy)
        
    def test_add_phase_error(self):
        """"""
        
        self.create_uvfits_outputs(do_autos=False)
        
        phase_err = 10
        seed = 34987
        
        args = []
        args.append(f"--numpy_seed={seed}")
        args.append("--uvfits=unittest_example1_band01.uvfits")
        args.append(f"--ant_gain_phase_error={phase_err}")
        
        ##Set the random seed, and make what should have been created
        np.random.seed(seed)
        
        ##the main script calls random to generate amps, so we need to do that here too
        amp_err = 0.0
        expec_gx = 1 + np.random.uniform(-amp_err, amp_err, self.num_antennas) + 1j*np.zeros(self.num_antennas)
        expec_gx = np.repeat(expec_gx, self.num_freqs)
        expec_gx.shape = (self.num_antennas, self.num_freqs)
        
        
        for ant in range(self.num_antennas):
            phase = np.random.uniform(-phase_err*(np.pi/180.0), phase_err*(np.pi/180.0))
            phase_grad = np.linspace(-phase, phase, self.num_freqs)
            expec_gx[ant] = expec_gx[ant]*np.exp(1j*phase_grad)
        
        expec_gy = 1 + np.random.uniform(-amp_err, amp_err, self.num_antennas) + 1j*np.zeros(self.num_antennas)
        expec_gy = np.repeat(expec_gy, self.num_freqs)
        expec_gy.shape = (self.num_antennas, self.num_freqs)
        
        for ant in range(self.num_antennas):
            phase = np.random.uniform(-phase_err*(np.pi/180.0), phase_err*(np.pi/180.0))
            phase_grad = np.linspace(-phase, phase, self.num_freqs)
            expec_gy[ant] = expec_gy[ant]*np.exp(1j*phase_grad)
        
        expec_gx[0] = complex(1,0)
        expec_gy[0] = complex(1,0)
        
        expec_Dx = np.zeros((self.num_antennas, self.num_freqs), dtype=complex)
        expec_Dy = np.zeros((self.num_antennas, self.num_freqs), dtype=complex)
        
        
        self.run_code_test_ant_gains(args, expec_gx, expec_Dx, expec_Dy, expec_gy)
        
    def test_add_amp_and_phase_error(self):
        """"""
        
        self.create_uvfits_outputs(do_autos=False)
        
        phase_err = 40
        amp_err = 0.5
        seed = 3984
        
        args = []
        args.append(f"--numpy_seed={seed}")
        args.append("--uvfits=unittest_example1_band01.uvfits")
        args.append(f"--ant_gain_phase_error={phase_err}")
        args.append(f"--ant_gain_amp_error={amp_err}")
        
        ##Set the random seed, and make what should have been created
        np.random.seed(seed)
        
        ##the main script calls random to generate amps, so we need to do that here too
        
        expec_gx = 1 + np.random.uniform(-amp_err, amp_err, self.num_antennas) + 1j*np.zeros(self.num_antennas)
        expec_gx = np.repeat(expec_gx, self.num_freqs)
        expec_gx.shape = (self.num_antennas, self.num_freqs)
        
        
        for ant in range(self.num_antennas):
            phase = np.random.uniform(-phase_err*(np.pi/180.0), phase_err*(np.pi/180.0))
            phase_grad = np.linspace(-phase, phase, self.num_freqs)
            expec_gx[ant] = expec_gx[ant]*np.exp(1j*phase_grad)
        
        expec_gy = 1 + np.random.uniform(-amp_err, amp_err, self.num_antennas) + 1j*np.zeros(self.num_antennas)
        expec_gy = np.repeat(expec_gy, self.num_freqs)
        expec_gy.shape = (self.num_antennas, self.num_freqs)
        
        for ant in range(self.num_antennas):
            phase = np.random.uniform(-phase_err*(np.pi/180.0), phase_err*(np.pi/180.0))
            phase_grad = np.linspace(-phase, phase, self.num_freqs)
            expec_gy[ant] = expec_gy[ant]*np.exp(1j*phase_grad)
        
        expec_gx[0] = complex(1,0)
        expec_gy[0] = complex(1,0)
        
        expec_Dx = np.zeros((self.num_antennas, self.num_freqs), dtype=complex)
        expec_Dy = np.zeros((self.num_antennas, self.num_freqs), dtype=complex)
        
        self.run_code_test_ant_gains(args, expec_gx, expec_Dx, expec_Dy, expec_gy)
        
    def test_add_amp_and_phase_autos_error(self):
        """"""
        
        self.create_uvfits_outputs(do_autos=True)
        
        phase_err = 40
        amp_err = 0.5
        seed = 3984
        
        args = []
        args.append(f"--numpy_seed={seed}")
        args.append("--uvfits=unittest_example1_band01.uvfits")
        args.append(f"--ant_gain_phase_error={phase_err}")
        args.append(f"--ant_gain_amp_error={amp_err}")
        
        ##Set the random seed, and make what should have been created
        np.random.seed(seed)
        
        ##the main script calls random to generate amps, so we need to do that here too
        
        expec_gx = 1 + np.random.uniform(-amp_err, amp_err, self.num_antennas) + 1j*np.zeros(self.num_antennas)
        expec_gx = np.repeat(expec_gx, self.num_freqs)
        expec_gx.shape = (self.num_antennas, self.num_freqs)
        
        
        for ant in range(self.num_antennas):
            phase = np.random.uniform(-phase_err*(np.pi/180.0), phase_err*(np.pi/180.0))
            phase_grad = np.linspace(-phase, phase, self.num_freqs)
            expec_gx[ant] = expec_gx[ant]*np.exp(1j*phase_grad)
        
        expec_gy = 1 + np.random.uniform(-amp_err, amp_err, self.num_antennas) + 1j*np.zeros(self.num_antennas)
        expec_gy = np.repeat(expec_gy, self.num_freqs)
        expec_gy.shape = (self.num_antennas, self.num_freqs)
        
        for ant in range(self.num_antennas):
            phase = np.random.uniform(-phase_err*(np.pi/180.0), phase_err*(np.pi/180.0))
            phase_grad = np.linspace(-phase, phase, self.num_freqs)
            expec_gy[ant] = expec_gy[ant]*np.exp(1j*phase_grad)
        
        expec_gx[0] = complex(1,0)
        expec_gy[0] = complex(1,0)
        
        expec_Dx = np.zeros((self.num_antennas, self.num_freqs), dtype=complex)
        expec_Dy = np.zeros((self.num_antennas, self.num_freqs), dtype=complex)
        
        self.run_code_test_ant_gains(args, expec_gx, expec_Dx, expec_Dy, expec_gy)
        

    def test_add_amp_and_phase_and_leakage_error(self):
        """"""
        
        self.create_uvfits_outputs(do_autos=True)
        
        phase_err = 40
        amp_err = 0.5
        seed = 3984
        psi_err = 0.02
        chi_err = 0.05
        
        args = []
        args.append(f"--numpy_seed={seed}")
        args.append("--uvfits=unittest_example1_band01.uvfits")
        args.append(f"--ant_gain_phase_error={phase_err}")
        args.append(f"--ant_gain_amp_error={amp_err}")
        args.append(f"--ant_leak_errs")
        args.append(f"{psi_err}")
        args.append(f"{chi_err}")
        
        ##Set the random seed, and make what should have been created
        np.random.seed(seed)
        
        ##the main script calls random to generate amps, so we need to do that here too
        
        expec_gx = 1 + np.random.uniform(-amp_err, amp_err, self.num_antennas) + 1j*np.zeros(self.num_antennas)
        expec_gx = np.repeat(expec_gx, self.num_freqs)
        expec_gx.shape = (self.num_antennas, self.num_freqs)
        
        
        for ant in range(self.num_antennas):
            phase = np.random.uniform(-phase_err*(np.pi/180.0), phase_err*(np.pi/180.0))
            phase_grad = np.linspace(-phase, phase, self.num_freqs)
            expec_gx[ant] = expec_gx[ant]*np.exp(1j*phase_grad)
        
        expec_gy = 1 + np.random.uniform(-amp_err, amp_err, self.num_antennas) + 1j*np.zeros(self.num_antennas)
        expec_gy = np.repeat(expec_gy, self.num_freqs)
        expec_gy.shape = (self.num_antennas, self.num_freqs)
        
        for ant in range(self.num_antennas):
            phase = np.random.uniform(-phase_err*(np.pi/180.0), phase_err*(np.pi/180.0))
            phase_grad = np.linspace(-phase, phase, self.num_freqs)
            expec_gy[ant] = expec_gy[ant]*np.exp(1j*phase_grad)
        
        expec_gx[0] = complex(1,0)
        expec_gy[0] = complex(1,0)
        
        psi_err *= np.pi/180.0
        chi_err *= np.pi/180.0
        
        expec_Dx = np.repeat(np.random.uniform(0, psi_err, self.num_antennas), self.num_freqs) - 1j*np.repeat(np.random.uniform(0, chi_err, self.num_antennas), self.num_freqs)
        expec_Dx.shape = (self.num_antennas, self.num_freqs)
        
        expec_Dy = -np.repeat(np.random.uniform(0, psi_err, self.num_antennas), self.num_freqs) + 1j*np.repeat(np.random.uniform(0, chi_err, self.num_antennas), self.num_freqs)
        expec_Dy.shape = (self.num_antennas, self.num_freqs)
        
        self.run_code_test_ant_gains(args, expec_gx, expec_Dx, expec_Dy, expec_gy)
        
    def run_code_test_noise(self, args, int_time=False, freq_reso=False):
        
        
        ##Actually run the code!!
        aiew.main(args)
        ##read in the resuts
        uvfits = aiew.UVFITS('instrumental.uvfits')
        
        if not int_time:
            int_time = uvfits.time_res
        if not freq_reso:
            freq_reso = uvfits.freq_res
        
        ##don't check auto correlations when testing cross-correlations
        baseline_cross = np.where(uvfits.b1s != uvfits.b2s)[0]
        baseline_auto = np.where(uvfits.b1s == uvfits.b2s)[0]
        
        for find, freq in enumerate(uvfits.all_freqs):
            
            expec_std = aiew.visibility_noise_stddev(freq, int_time, freq_reso)
            
            found_std_xx_re = np.std(np.real(uvfits.visibilities[baseline_cross, find, 0]))
            found_std_xx_im = np.std(np.imag(uvfits.visibilities[baseline_cross, find, 0]))
            found_std_yy_re = np.std(np.real(uvfits.visibilities[baseline_cross, find, 1]))
            found_std_yy_im = np.std(np.imag(uvfits.visibilities[baseline_cross, find, 1]))
            
            npt.assert_allclose(found_std_xx_re, expec_std, rtol=0.05)
            npt.assert_allclose(found_std_xx_im, expec_std, rtol=0.05)
            npt.assert_allclose(found_std_yy_re, expec_std, rtol=0.05)
            npt.assert_allclose(found_std_yy_im, expec_std, rtol=0.05)
            
            if len(baseline_auto) > 0:
            
                found_std_xx_re = np.std(np.real(uvfits.visibilities[baseline_auto, find, 0]))
                found_std_xx_im = np.std(np.imag(uvfits.visibilities[baseline_auto, find, 0]))
                found_std_yy_re = np.std(np.real(uvfits.visibilities[baseline_auto, find, 1]))
                found_std_yy_im = np.std(np.imag(uvfits.visibilities[baseline_auto, find, 1]))
                
                npt.assert_allclose(found_std_xx_re, np.sqrt(2)*expec_std, rtol=0.10)
                npt.assert_allclose(found_std_xx_im, np.sqrt(2)*expec_std, rtol=0.10)
                npt.assert_allclose(found_std_yy_re, np.sqrt(2)*expec_std, rtol=0.10)
                npt.assert_allclose(found_std_yy_im, np.sqrt(2)*expec_std, rtol=0.10)
            
    def test_add_visi_noise(self):
        """"""
        
        self.create_uvfits_outputs(do_autos=False, num_ants=100)
        
        args = []
        args.append("--uvfits=unittest_example1_band01.uvfits")
        args.append(f"--add_visi_noise")
        
        self.run_code_test_noise(args)
        
    def test_add_visi_noise_set_reso(self):
        """"""
        
        self.create_uvfits_outputs(do_autos=True, num_ants=100)
        
        int_time = 1000*60*60.0
        freq_reso = 1e+6
        
        args = []
        args.append("--uvfits=unittest_example1_band01.uvfits")
        args.append(f"--visi_noise_int_time={int_time}")
        args.append(f"--visi_noise_freq_reso={freq_reso}")
        
        self.run_code_test_noise(args, int_time, freq_reso)
        
    def run_code_test_reflections(self, args, metafits):
        
        
        ##Actually run the code!!
        aiew.main(args)
        ##read in the resuts
        uvfits = aiew.UVFITS('instrumental.uvfits')
        
        with fits.open(metafits) as f:
            antenna_order = np.argsort(f[1].data['Tile'])
            
            selection = np.arange(0,len(antenna_order),2)
            ##this has both XX and YY pol, assume they have the same cable lengths
            antenna_order = antenna_order[selection]
            
            ##the length causing reflections is in the flavour column
            flavours = f[1].data['Flavors'][antenna_order]
            cable_lengths = np.array([float(flavour.split('_')[1]) for flavour in flavours])
        
        meta_delays = aiew.get_cable_delay(cable_lengths)
        
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(4, 2, figsize=(10, 6))
        
        autos = np.where(uvfits.b1s == uvfits.b2s)[0]
        crosses = np.where(uvfits.b1s != uvfits.b2s)[0]
        
        delays = np.fft.fftshift(np.fft.fftfreq(len(uvfits.all_freqs), uvfits.freq_res))
        
        # cut_delays = np.where(delays > 0 )[0]
        use = np.where(delays != 0)
        delays = delays[use]
        
        # for ind, ant in enumerate([0, 2, 50, 100]):
        
        base_plot = 0
        # for baseline in crosses[:8128]:
        for baseline in range(8256):
            
            b1 = uvfits.b1s[baseline]
            b2 = uvfits.b2s[baseline]
            
            delay1, delay2 = meta_delays[b1], meta_delays[b2]
            
            # print(b1, b2, delay1, delay2, delay1 + delay2, np.abs(delay2 - delay1))
            
            
            
            
            fft_gain = np.fft.fftshift(np.fft.fft(uvfits.visibilities[baseline, :,  0]))
            fft_gain = fft_gain[use]
            
            idx = np.argsort(np.abs(fft_gain))[::-1]
            
            
            ##OK, so times the peak delay sits across two samples, so we don't
            ##
            from scipy.signal import find_peaks
            idx, props = find_peaks(np.abs(fft_gain), distance=3)
            # print(idx)
            
            idx1 = np.argsort(np.abs(fft_gain[idx]))[::-1][:2]
            
            found_delays = np.sort(np.abs(delays[idx[idx1]]))
            
            # found_delays = np.sort([np.abs(delays[idx[0]]), np.abs(delays[idx[1]])])
            expec_delays = np.sort([delay1, delay2])
            
            
            
            try:
                npt.assert_allclose(found_delays, expec_delays, rtol=0.05)
            except:
                # print(expec_delays, found_delays, delay2 - delay1)
                print(expec_delays, found_delays)
                print(delays[idx[idx1]], np.abs(fft_gain[idx[idx1]]))
                
                plt.plot(delays/1e-6, np.abs(fft_gain))
                plt.savefig('ting.png')
                plt.close()
                exit()
            
            # print(delays[idx[1]]/1e-6)
            # print()
            
            # auto = autos[ant]
            # axs[ind, 0].plot(uvfits.all_freqs/1e+6, np.abs(uvfits.visibilities[auto, :,  0]))
            # fft_gain = np.fft.fft(uvfits.visibilities[auto, :,  0])
            
            base_plots = [0, 2, 50, 100]
            if baseline in base_plots:
            
            # if base_plot <= 3:
                
                base_plot = base_plots.index(baseline)
            
                axs[base_plot, 0].plot(uvfits.all_freqs/1e+6, np.real(uvfits.visibilities[baseline, :,  0]))
                
                # axs[ind, 1].plot(delays[cut_delays]/1e-6, np.real(fft_gain[cut_delays]))
                axs[base_plot, 1].plot(delays/1e-6, np.abs(fft_gain), 'k')
                
                for ind, delay in enumerate(expec_delays):
                    axs[base_plot, 1].axvline(delay/1e-6, color=f'C{ind:d}', label=f'Delay {delay/1e-6:.1f}', linestyle='--',alpha=0.5)
                    axs[base_plot, 1].axvline(-delay/1e-6, color=f'C{ind:d}', linestyle='--',alpha=0.5)
                    
                axs[base_plot, 1].axvline(-delay/1e-6, color=f'C{ind:d}', linestyle='--',alpha=0.5)
                    
                axs[base_plot, 1].set_xlim(-5, 5)
                axs[base_plot, 1].legend()
                
        axs[0,0].set_title('Real part visibility')
        axs[0,1].set_title('FT visi along frequency')
          
        axs[3, 0].set_xlabel("Frequency MHz")
        axs[3, 1].set_xlabel("Delay $\mu$s")
            
        plt.tight_layout()
        fig.savefig("reflections_test.png", bbox_inches='tight')
        plt.close()
        
    def test_add_cable_reflect(self):
        """"""
        
        np.random.seed(93208)
        self.create_uvfits_outputs(do_autos=True, num_ants=128, 
                                   num_freqs=384, ch_width=80e+3)
        
        metafits = f"{code_dir}/1088285720_metafits.fits"
        
        args = []
        args.append("--uvfits=unittest_example1_band01.uvfits")
        args.append(f"--cable_reflection_from_metafits={metafits}")
        args.append(f"--cable_reflection_coeff_amp=0.05")
        
        self.run_code_test_reflections(args, metafits)
        
    def test_add_bandpass(self):
        """"""
        
        # np.random.seed(93208)
        # self.create_uvfits_outputs(do_autos=True, num_ants=128, 
        #                            num_freqs=384, ch_width=80e+3)
        
        # metafits = f"{code_dir}/1088285720_metafits.fits"
        
        # args = []
        # args.append("--uvfits=unittest_example1_band01.uvfits")
        # args.append(f"--cable_reflection_from_metafits={metafits}")
        # args.append(f"--cable_reflection_coeff_amp=0.05")
        
        # self.run_code_test_reflections(args, metafits)
        
        import matplotlib.pyplot as plt
        
        
        freq_res = 80e+3
        
        
        fig, axs = plt.subplots(1, 1, figsize=(6, 6))
        
        bandpass = importlib_resources.files(wodenpy).joinpath("bandpass_1kHz.txt")
        bandpass = np.loadtxt(bandpass)
        
        axs.plot(np.arange(5e+2, 1.28e+6, 1e+3)/1e+6, bandpass, 'k-', label='Jake Jones bandpass')
        
        for freq_res in [10e+3, 20e+3, 40e+3, 80e+3]:
            freqs = np.arange(freq_res/2, 1.28e+6, freq_res) / 1e+6
            single_bandpass, full_bandpass = aiew.calculate_bandpass(freq_res)
            
            axs.plot(freqs, single_bandpass, label=f'{freq_res/1e+3:.1f} kHz',
                     linestyle='none', mfc='none', marker='o')
            
            
        axs.set_xlabel('Frequency MHz')
        axs.set_ylabel('Bandpass')
            
        axs.legend()
        plt.tight_layout()
        fig.savefig("bandpass_test.png", bbox_inches='tight')
        plt.close()
        

##Run the test
if __name__ == '__main__':
   unittest.main()
