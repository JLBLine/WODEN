from sys import path
import os
import unittest
import numpy as np

np.random.seed(20398)

VELC = 299792458.0

##Do some disgusting path finding exercise, there must be a better way
##to do this
fileloc = os.path.realpath(__file__)
path.append('{:s}/../../src'.format(('/').join(fileloc.split('/')[:-1])))

##Code we are testing
import run_woden as rw

##Vehicle for running tests
class Test(unittest.TestCase):
    """Tests the `rw.load_data` function, which should read in a binary
    file as output by WODEN, into various arrays. Test by writing out
    a binary file with known input params, reading it in with the funtion
    under test, and compare the inputs to outputs"""

    def reshape_v_container(self, v_slice, num_baselines, num_freq_channels,
                            num_time_steps):
        """Reshape a real/imag single polarisation slice of v_container to match
        the order of data in the WODEN binary file"""
        out_slice = np.empty(num_baselines*num_freq_channels*num_time_steps)

        for time_ind in np.arange(num_time_steps):
            for freq_ind in np.arange(num_freq_channels):
                index_in = time_ind*num_baselines
                index_out = time_ind*num_baselines*num_freq_channels + freq_ind*num_baselines
                out_slice[index_out:index_out+num_baselines] = v_slice[index_in:index_in+num_baselines, freq_ind]

        return out_slice

    def test_load_data(self):
        """Does the binary file writing/reading and comparisons"""

        num_baselines = 20
        num_freq_channels = 4
        num_time_steps = 2

        num_visis = num_baselines*num_freq_channels*num_time_steps

        ##In the real WODEN output, the u values are repeated for each frequency
        ##as they are in metres. First make two different random sets of u,v,w
        ##one for each time step

        input_u_t1 = np.random.uniform(-1000, 1000, num_baselines)
        input_v_t1 = np.random.uniform(-1000, 1000, num_baselines)
        input_w_t1 = np.random.uniform(-100, 100, num_baselines)

        input_u_t2 = np.random.uniform(-1000, 1000, num_baselines)
        input_v_t2 = np.random.uniform(-1000, 1000, num_baselines)
        input_w_t2 = np.random.uniform(-100, 100, num_baselines)

        ##Tile them up in the way the would come out of WODEN
        input_u = np.empty(num_visis)
        input_v = np.empty(num_visis)
        input_w = np.empty(num_visis)

        input_u[:num_baselines*num_freq_channels] = np.tile(input_u_t1, num_freq_channels)
        input_v[:num_baselines*num_freq_channels] = np.tile(input_v_t1, num_freq_channels)
        input_w[:num_baselines*num_freq_channels] = np.tile(input_w_t1, num_freq_channels)

        input_u[num_baselines*num_freq_channels:] = np.tile(input_u_t2, num_freq_channels)
        input_v[num_baselines*num_freq_channels:] = np.tile(input_v_t2, num_freq_channels)
        input_w[num_baselines*num_freq_channels:] = np.tile(input_w_t2, num_freq_channels)

        ##make some fake visibilities
        input_xx_re = np.random.uniform(-1000, 1000, num_visis)
        input_xx_im = np.random.uniform(-1000, 1000, num_visis)
        input_xy_re = np.random.uniform(-1000, 1000, num_visis)
        input_xy_im = np.random.uniform(-1000, 1000, num_visis)
        input_yx_re = np.random.uniform(-1000, 1000, num_visis)
        input_yx_im = np.random.uniform(-1000, 1000, num_visis)
        input_yy_re = np.random.uniform(-1000, 1000, num_visis)
        input_yy_im = np.random.uniform(-1000, 1000, num_visis)

        ##Essentially when writing the binary all the data goes into a massive
        ##1D array, so assemble in correct order and write out

        binary_array = np.empty(11*num_visis)
        binary_array[0*num_visis:1*num_visis] = input_u
        binary_array[1*num_visis:2*num_visis] = input_v
        binary_array[2*num_visis:3*num_visis] = input_w
        binary_array[3*num_visis:4*num_visis] = input_xx_re
        binary_array[4*num_visis:5*num_visis] = input_xx_im
        binary_array[5*num_visis:6*num_visis] = input_xy_re
        binary_array[6*num_visis:7*num_visis] = input_xy_im
        binary_array[7*num_visis:8*num_visis] = input_yx_re
        binary_array[8*num_visis:9*num_visis] = input_yx_im
        binary_array[9*num_visis:10*num_visis] = input_yy_re
        binary_array[10*num_visis:11*num_visis] = input_yy_im

        ##Write them out to a binary file
        filename = "test_load_data.dat"
        binary_array.astype('float32').tofile(filename)

        ##Read in binary file using code under test
        output_u, output_v, output_w, v_container = rw.load_data(filename=filename,
                                                 num_baselines=num_baselines,
                                                 num_freq_channels=num_freq_channels,
                                                 num_time_steps=num_time_steps)

        ##Check the u,v,w positions are to millimetre precision
        mm_accuracy = 1e-3 / VELC

        ##First half of the output u,v,w should equal first time step u,v,w
        self.assertTrue(np.allclose(input_u_t1 / VELC, output_u[:num_baselines], atol=mm_accuracy))
        self.assertTrue(np.allclose(input_v_t1 / VELC, output_v[:num_baselines], atol=mm_accuracy))
        self.assertTrue(np.allclose(input_w_t1 / VELC, output_w[:num_baselines], atol=mm_accuracy))

        ##Same for the second half u,v,w
        self.assertTrue(np.allclose(input_u_t2 / VELC, output_u[num_baselines:], atol=mm_accuracy))
        self.assertTrue(np.allclose(input_v_t2 / VELC, output_v[num_baselines:], atol=mm_accuracy))
        self.assertTrue(np.allclose(input_w_t2 / VELC, output_w[num_baselines:], atol=mm_accuracy))

        ##Reshape the data read in from binary file to what was input
        output_xx_re = self.reshape_v_container(v_container[:,0,0,:,0,0],
                            num_baselines, num_freq_channels, num_time_steps)
        output_xx_im = self.reshape_v_container(v_container[:,0,0,:,0,1],
                            num_baselines, num_freq_channels, num_time_steps)
        output_xy_re = self.reshape_v_container(v_container[:,0,0,:,2,0],
                            num_baselines, num_freq_channels, num_time_steps)
        output_xy_im = self.reshape_v_container(v_container[:,0,0,:,2,1],
                            num_baselines, num_freq_channels, num_time_steps)
        output_yx_re = self.reshape_v_container(v_container[:,0,0,:,3,0],
                            num_baselines, num_freq_channels, num_time_steps)
        output_yx_im = self.reshape_v_container(v_container[:,0,0,:,3,1],
                            num_baselines, num_freq_channels, num_time_steps)
        output_yy_re = self.reshape_v_container(v_container[:,0,0,:,1,0],
                            num_baselines, num_freq_channels, num_time_steps)
        output_yy_im = self.reshape_v_container(v_container[:,0,0,:,1,1],
                            num_baselines, num_freq_channels, num_time_steps)

        ##Test things are what they're supposed to be
        self.assertTrue(np.allclose(input_xx_re, output_xx_re, atol=1e-10))
        self.assertTrue(np.allclose(input_xx_im, output_xx_im, atol=1e-10))
        self.assertTrue(np.allclose(input_xy_re, output_xy_re, atol=1e-10))
        self.assertTrue(np.allclose(input_xy_im, output_xy_im, atol=1e-10))
        self.assertTrue(np.allclose(input_yx_re, output_yx_re, atol=1e-10))
        self.assertTrue(np.allclose(input_yx_im, output_yx_im, atol=1e-10))
        self.assertTrue(np.allclose(input_yy_re, output_yy_re, atol=1e-10))
        self.assertTrue(np.allclose(input_yy_im, output_yy_im, atol=1e-10))

        ##The weights should all be one
        self.assertTrue((np.all(v_container[:,0,0,:,0,2] == 1.0)))
        self.assertTrue((np.all(v_container[:,0,0,:,1,2] == 1.0)))
        self.assertTrue((np.all(v_container[:,0,0,:,2,2] == 1.0)))
        self.assertTrue((np.all(v_container[:,0,0,:,3,2] == 1.0)))

##Run the test
if __name__ == '__main__':
   unittest.main()
