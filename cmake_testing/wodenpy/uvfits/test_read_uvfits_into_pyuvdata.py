from sys import path
import os
import unittest
import numpy as np
from astropy.io import fits

np.random.seed(179)

from wodenpy.uvfits import wodenpy_uvfits

##Vehicle for running tests
class Test(unittest.TestCase):
    def make_dummy_intputs(self):
        """Makes some inputs used to pass to wodenpy_uvfits.create_uvfits"""
        self.date = "2019-06-12T13:04:12"
        self.gst0_deg = 260.0303917560829
        self.degpdy = 360.98563614850775
        self.ut1utc = -0.17421478816666666
        self.num_antennas = 10
        self.freq_cent = 160e+6
        self.telescope_name = "MWA"
        self.longitude = 116.670813889
        self.latitude = -26.7033194444
        self.array_height = 377.0
        self.central_freq_chan = 1
        self.ch_width = 40e3
        self.ra_point = 0.0
        self.dec_point = -26.7
        self.num_baselines = int(((self.num_antennas - 1)*self.num_antennas) / 2)
        self.num_time_steps = 2
        self.num_freq_channels = 3
        self.jd_midnight = 2458646.5
        self.gitlabel = 'as987has'
        self.XYZ_array = np.arange(self.num_antennas*3)
        self.XYZ_array.shape = (self.num_antennas, 3)

        self.output_uvfits_name = "unittest_example.uvfits"


    def create_uvfits_file(self):
        """Runs the `wodenpy_uvfits.create_uvfits` function, which should take a whole
        heap of inputs and write out a uvits."""

        ##Some input test params
        self.make_dummy_intputs()

        ##Create an antenna table
        ant_table = wodenpy_uvfits.make_antenna_table(XYZ_array=self.XYZ_array,
                      telescope_name=self.telescope_name,
                      num_antennas=self.num_antennas, freq_cent=self.freq_cent,
                      date=self.date, gst0_deg=self.gst0_deg,
                      degpdy=self.degpdy,
                      longitude=self.longitude, latitude=self.latitude,
                      array_height=self.array_height)

        ##Create some dummy data inputs
        v_container = np.empty((self.num_time_steps*self.num_baselines, 1, 1, self.num_freq_channels, 4, 3))
        random_data = np.random.uniform(-100, 100, self.num_time_steps*self.num_baselines*self.num_freq_channels*4*3)
        random_data.shape = (self.num_time_steps*self.num_baselines, self.num_freq_channels, 4, 3)
        v_container[:,0,0,:,:,:] = random_data

        uu = (np.arange(self.num_baselines*self.num_time_steps) + 1) / 3e+8
        vv = (np.arange(self.num_baselines*self.num_time_steps) + 1) / 3e+8
        ww = (np.arange(self.num_baselines*self.num_time_steps) + 1) / 3e+8

        ##Create all baseline pairs and encode using the ol' miriad ways
        baselines = np.empty(self.num_baselines)

        base_ind = 0
        for ant1 in range(1,self.num_antennas):
            for ant2 in range(ant1 + 1,self.num_antennas + 1):
                baselines[base_ind] = wodenpy_uvfits.RTS_encode_baseline(ant1, ant2)
                base_ind += 1

        ##Repeat the baseline pairs for as many timesteps as necessary
        baselines_array = np.tile(baselines, self.num_time_steps)

        date_array = np.ones(self.num_time_steps*self.num_baselines)*0.54458333319
        ##Add a second of time onto the second half of the date array
        date_array[self.num_baselines:] += (1 / (24.0*60.0*60.0))

        wodenpy_uvfits.create_uvfits(freq_cent=self.freq_cent, ch_width=self.ch_width,
                  central_freq_chan=self.central_freq_chan,
                  ra_point=self.ra_point, dec_point=self.dec_point,
                  output_uvfits_name=self.output_uvfits_name,
                  jd_midnight=self.jd_midnight, gitlabel=self.gitlabel,
                  v_container=v_container,
                  uu=uu, vv=vv, ww=ww,
                  baselines_array=baselines_array, date_array=date_array,
                  hdu_ant=ant_table, telescope_name=self.telescope_name,
                  longitude=self.longitude, latitude=self.latitude,
                  array_height=self.array_height)


    def test_pyuvdata_reads_uvfits_file(self):
        """Runs the `self.create_uvfits_file` function, which should take a whole
        heap of inputs and write out a uvits. It then tries to read the
        created uvfits file into a pyuvdata `UVData` object"""

        try:

            from pyuvdata import UVData
            ##Make a uvfits file with dummy data
            self.create_uvfits_file()
            ##Setup a UVData object and read in the uvfits file. This will print
            ##out a warning about u,v,w coords as I'm not setting up real coords
            ##here, but shouldn't error if all is well
            UV = UVData()
            UV.read('unittest_example.uvfits')

        except ModuleNotFoundError:
            print('There is no pyuvdata module, so cannot test uvfits work with pyuvdata')

##Run the test
if __name__ == '__main__':
   unittest.main()
