from sys import path
import os
import unittest
import numpy as np
from copy import deepcopy
np.random.seed(7687)

##Do some disgusting path finding exercise, there must be a better way
##to do this
##Do some disgusting path finding exercise, there must be a better way
##to do this
fileloc = os.path.realpath(__file__)
path.append('{:s}/../../../wodenpy/phase_rotate/'.format(('/').join(fileloc.split('/')[:-1])))

from wodenpy.phase_rotate import remove_phase_track

##Constants
D2R = np.pi / 180.0
SOLAR2SIDEREAL = 1.00274
MWA_LAT_RAD = -0.46606083776035967
VELC = 299792458.0

##Vehicle for running tests
class Test(unittest.TestCase):

    def enh2xyz(self, east, north, height, latitude):
        """
        Takes local east, north, height coords for a given latitude (radians)
        and returns local X,Y,Z coords. This is a copy of what is in run_woden.py,
        but trying to isolate `remove_phase_tracking` so copy functionality here
        """

        sl = np.sin(latitude)
        cl = np.cos(latitude)
        X = -north*sl + height*cl
        Y = east
        Z = north*cl + height*sl

        return X,Y,Z

    def get_uvw(self, xdiff, ydiff, zdiff, dec, ha):
        """Get uvw coords from baseline lengths in X,Y,Z (xdiff, ydiff, zdiff)
        towards hour angle (ha) declination (dec)"""
        u = np.sin(ha)*xdiff + np.cos(ha)*ydiff
        v = -np.sin(dec)*np.cos(ha)*xdiff + np.sin(dec)*np.sin(ha)*ydiff + np.cos(dec)*zdiff
        w = np.cos(dec)*np.cos(ha)*xdiff - np.cos(dec)*np.sin(ha)*ydiff + np.sin(dec)*zdiff

        return u,v,w

    def get_lmn(self, ra, ra0, dec, dec0):
        """Get l,m,n at a given ra,dec with a given phase centre ra0,dec0"""

        ##RTS way of doing it
        cdec0 = np.cos(dec0)
        sdec0 = np.sin(dec0)
        cdec = np.cos(dec)
        sdec = np.sin(dec)
        cdra = np.cos(ra-ra0)
        sdra = np.sin(ra-ra0)
        l = cdec*sdra
        m = sdec*cdec0 - cdec*sdec0*cdra
        n = sdec*sdec0 + cdec*cdec0*cdra

        return l,m,n

    def get_xyz_diffs(self, X, Y, Z):
        num_ants = len(X)

        xdiffs = []
        ydiffs = []
        zdiffs = []

        for ant1 in np.arange(num_ants - 1):
            for ant2 in np.arange(ant1 + 1, num_ants):
                xdiffs.append(X[ant1] - X[ant2])
                ydiffs.append(Y[ant1] - Y[ant2])
                zdiffs.append(Z[ant1] - Z[ant2])

        return np.array(xdiffs), np.array(ydiffs), np.array(zdiffs)

    def create_uvw_and_vcontainer(self, xdiffs, ydiffs, zdiffs, frequencies,
                                  ra_source, dec_source, initial_lst,
                                  num_time_steps, time_res,
                                  ra_phase=None, dec_phase=None, phase_tracked=True):
        """Data is stored in the uvfits in ``v_container` in shape:

        (num_time_steps*num_baselines,1,1,num_freq_channels,4,3)

        Calculate the u,v,w and l,m,n to calculate visibilities for all time
        steps. Just simulate a single source at location ra_source,dec_source
        with phase centre ra_phase,dec_phase. For no phase tracking, just
        point the u,v,w at zenith (which means sourc l,m,n changes with time)

        u,v,w change with time/freq so need to know number of times, time
        res to move lst along, and frequencies"""

        num_baselines = len(xdiffs)
        num_freqs = len(frequencies)

        v_container = np.empty((num_baselines*num_time_steps, 1, 1, num_freqs, 4, 3))
        wws_seconds = np.empty(num_baselines*num_time_steps)

        ##Set the weights to one for everything
        v_container[:,0,0,:,:,2] = 1.0

        ##Fill 'em up
        for time_ind in np.arange(num_time_steps):
            ##Adjust the lst for sky rotation. In main code usually sample at
            ##the centre of the time sample (so +0.5*time_res) but I want
            ##source to be at phase centre at first time step here so sample at
            ##beginning instead. That means first time step is purely real
            ##in the non-phase track case
            lst = initial_lst + time_ind*time_res*SOLAR2SIDEREAL*(15.0/3600.0)
            ##u,v,w centred on phase centre when phase tracking, zenith if not
            if phase_tracked:
                ha_phase = lst - ra_phase
                uu_metres, vv_metres, ww_metres = self.get_uvw(xdiffs, ydiffs, zdiffs, dec=dec_phase, ha=ha_phase)
            else:
                uu_metres, vv_metres, ww_metres = self.get_uvw(xdiffs, ydiffs, zdiffs, dec=MWA_LAT_RAD, ha=0.0)

            ##u,v,w are stored in a uvfits in units of seconds
            wws_seconds[time_ind*num_baselines:(time_ind + 1)*num_baselines] = ww_metres / VELC

            for freq_ind, freq in enumerate(frequencies):

                wavelength = VELC / freq
                uu_waves = uu_metres / wavelength
                vv_waves = vv_metres / wavelength
                ww_waves = ww_metres / wavelength

                if phase_tracked:
                    ls, ms, ns = self.get_lmn(ra_source, ra_phase, dec_source, dec_phase)
                    ##Phase tracking case - has w*(n-1) in the exponential
                    visis = np.exp(2j*np.pi*(uu_waves*ls + vv_waves*ms + ww_waves*(ns - 1)))

                else:
                    ##If not phase tracking, centre of the coord system is zenith
                    ls, ms, ns = self.get_lmn(ra_source, lst, dec_source, MWA_LAT_RAD)
                    ##No phase tracking case - has w*n in the exponential
                    visis = np.exp(2j*np.pi*(uu_waves*ls + vv_waves*ms + ww_waves*ns))

                ##Stick a different gain on each polarisation to check correct
                ##mapping from remove_phase_track.remove_phase_tracking`
                v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,0,0] = np.real(visis)
                v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,0,1] = np.imag(visis)

                v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,1,0] = 2*np.real(visis)
                v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,1,1] = 2*np.imag(visis)

                v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,2,0] = 3*np.real(visis)
                v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,2,1] = 3*np.imag(visis)

                v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,3,0] = 4*np.real(visis)
                v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,3,1] = 4*np.imag(visis)

        return wws_seconds, v_container

    def test_remove_phase_tracking(self):
        """Tests the `remove_phase_track.remove_phase_tracking` function, which should remove
        the phase tracking applied to visibilities. The original MWA correlator
        did not phase track, so the RTS expects no phase tracking in it's inputs.
        First, set up a random array layout and calculate some phase tracked and
        non-phase tracked visibilities. The unwrap the phase tracking from the
        phase tracked case, and check it matches the non-phase tracked case"""

        num_antennas = 50

        ##Make a random array layout
        east = np.random.uniform(-1000, 1000, num_antennas)
        north = np.random.uniform(-1000, 1000, num_antennas)
        height = np.random.uniform(0, 10, num_antennas)

        ##Stick at the MRO, convert to X,Y,Z and find the baseline lengths
        X, Y, Z = self.enh2xyz(east, north, height, MWA_LAT_RAD)
        xdiffs, ydiffs, zdiffs = self.get_xyz_diffs(X, Y, Z)

        print(f'Number of baselines: {len(xdiffs)}')

        ##Setup some overserving variables for zenith

        initial_lst = 0.0
        ra_source = 10.0*D2R
        dec_source = -15.0*D2R

        ##Random phase centre
        ra_phase = 40*D2R
        dec_phase = -50.0*D2R

        ##Other inputs params
        num_time_steps = 10
        time_res = 2.0
        num_freqs = 10

        frequencies = np.arange(num_freqs)*10e+6 + 100e+6

        ##Make data with phase tracking
        w_seconds_phased, v_container_phased = self.create_uvw_and_vcontainer(xdiffs,
                                      ydiffs, zdiffs, frequencies,
                                      ra_source, dec_source, initial_lst,
                                      num_time_steps, time_res,
                                      ra_phase=ra_phase, dec_phase=dec_phase,
                                      phase_tracked=True)

        ##Make data with no phase tracking
        w_seconds_nophased, v_container_nophased = self.create_uvw_and_vcontainer(xdiffs,
                                      ydiffs, zdiffs, frequencies,
                                      ra_source, dec_source, initial_lst,
                                      num_time_steps, time_res,
                                      phase_tracked=False)

        ##Check we haven't accidentally made the same set of visis twice
        self.assertFalse(np.allclose(v_container_nophased, v_container_phased, atol=1e-3))

        ## Take out the phase tracking from the phase tracked case
        ## using the code under test
        v_container_unwrapped = remove_phase_track.remove_phase_tracking(frequencies=frequencies,
                                  wws_seconds=w_seconds_phased, num_time_steps=num_time_steps,
                                  v_container=v_container_phased,
                                  num_baselines=len(xdiffs))

        ## Require the no phase track case to be withing 1e-10 precision of
        ## the unwrapped phase track case
        self.assertTrue(np.allclose(v_container_nophased, v_container_unwrapped, atol=1e-10))

##Run the test
if __name__ == '__main__':
   unittest.main()
