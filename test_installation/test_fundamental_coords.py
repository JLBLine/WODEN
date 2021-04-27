import pycuda.driver as cuda
import pycuda.autoinit
from pycuda.compiler import SourceModule
import numpy as np
import pycuda.gpuarray as gpuarray
from os import getcwd


D2R = np.pi / 180.0
VELC = 299792458.0

# Save the path to this directory
pwd = getcwd()

# include_dir =

##When running np.isclose, have a relative tolerance of 5e-7 and an absolute
##tolerance of 0.0, which means if abs(array1 - array2) <= 5e-7 * abs(array1)
##is True, results pass
rtol = 5e-7
atol = 0.0


def make_cuda_float_array(a):
    """Take a numpy array `a`, convert to float32 and make a copy on the GPU"""
    a = a.astype(np.float32)
    a_gpu = gpuarray.to_gpu(a)

    return a, a_gpu

class TestFundamentalCoords():
    def _get_kernel(self, kernel_name):
        """Read in fundamental_coords.cu CUDA code and use pycuda to create
        a SourceModule with the relevant include path. Compile and return the
        kernel specified by `kernel_name`"""
        ##TODO sort out hard-coded paths
        loaded_from_source = open("../src/fundamental_coords.cu",'r').read()
        mod = SourceModule(loaded_from_source, include_dirs=['/home/jline/software/WODEN/include/'])
        kernel = mod.get_function(kernel_name)

        return kernel


    def _get_uvw_python(self, X_diff, Y_diff, Z_diff, sdec0, cdec0, cha0s, sha0s,
                       wavelengths):
        """python equivalent to `fundamental_coords::kern_calc_uvw`. Only
        difference is that `kern_calc_uvw` performs a modulus on each thread
        to access correct `X_diff, Y_diff, Z_diff`. Here we'll just tile
        those arrays"""

        num_tiles = int(len(wavelengths) / len(X_diff))

        X_diff_tile = np.tile(X_diff, num_tiles).astype(np.float32)
        Y_diff_tile = np.tile(Y_diff, num_tiles).astype(np.float32)
        Z_diff_tile = np.tile(Z_diff, num_tiles).astype(np.float32)

        u_metres = sha0s*X_diff_tile + cha0s*Y_diff_tile
        v_metres = -sdec0*cha0s*X_diff_tile + sdec0*sha0s*Y_diff_tile + cdec0*Z_diff_tile
        w_metres = cdec0*cha0s*X_diff_tile - cdec0*sha0s*Y_diff_tile + sdec0*Z_diff_tile

        u = u_metres / wavelengths
        v = v_metres / wavelengths
        w = w_metres / wavelengths

        return u, v, w, u_metres, v_metres, w_metres

    def test_kern_calc_uvw(self):
        """Call a `kern_calc_uvw` with 100 baselines, two frequencies, and two
        time steps. Check against an equivalent python function."""

        ##Grab the kernel from
        kernel_calc_uvw = self._get_kernel("kern_calc_uvw")

        num_baselines = 100
        num_times = 2
        num_freqs = 2

        num_visis = num_baselines * num_times * num_freqs

        ##`make_cuda_float_array` returns a np.float32 version of the array on
        ##the host, and also a float device version

        ##Input X,Y,Z baseline lengths
        X_diff, d_X_diff = make_cuda_float_array(np.arange(num_baselines)+1)
        Y_diff, d_Y_diff = make_cuda_float_array(np.arange(num_baselines)+1)
        Z_diff, d_Z_diff = make_cuda_float_array(np.arange(num_baselines)+1)

        ##kern_calc_uvw expects d_sha0s, d_cha0s, d_wavelengths to include
        ##values for all time and frequency steps, so do some repeating/tiling
        ha0s = np.empty(num_visis, dtype=np.float32)
        ##Set first time step to zero, for all frequencies and baselines
        ha0s[:num_baselines*num_freqs] = 0.0
        ##Let's say the sky has moved 10 arcminute for the second time step
        ha0s[num_baselines*num_freqs:] = (10/60.0)*D2R

        ##Setup enough frequencies to cover the first time step
        wavelengths = np.empty(num_baselines*num_freqs, dtype=np.float32)
        ##Set first frequency to 150 MHz
        wavelengths[:num_baselines] = VELC / 150e+6
        ##Set second frequency to 200 MHz
        wavelengths[num_baselines:] = VELC / 200e+6

        ##Need to repeat this for both time steps, so tile it
        wavelengths = np.tile(wavelengths, num_times)

        ##Phase centre stays constant in declination so only need one value
        dec0 = -30.0*D2R
        sdec0, cdec0 = np.float32(np.sin(dec0)), np.float32(np.cos(dec0))

        ##Stick relevant arrays onto the GPU
        sha0s, d_sha0s = make_cuda_float_array(np.sin(ha0s))
        cha0s, d_cha0s = make_cuda_float_array(np.cos(ha0s))
        wavelengths, d_wavelengths = make_cuda_float_array(wavelengths)

        ##Arrays to hold outputs
        u, d_u = make_cuda_float_array(np.zeros(num_visis))
        v, d_v = make_cuda_float_array(np.zeros(num_visis))
        w, d_w = make_cuda_float_array(np.zeros(num_visis))

        u_metres, d_u_metres = make_cuda_float_array(np.zeros(num_visis))
        v_metres, d_v_metres = make_cuda_float_array(np.zeros(num_visis))
        w_metres, d_w_metres = make_cuda_float_array(np.zeros(num_visis))

        ##Setup the grid/block dimensions needed for kernel call
        threads_x = 128
        grid_x = int(np.ceil(num_visis / threads_x))

        ##Call device code
        kernel_calc_uvw(d_X_diff, d_Y_diff, d_Z_diff,
                d_u_metres, d_v_metres, d_w_metres,
                d_u, d_v, d_w, d_wavelengths,
                sdec0, cdec0, d_cha0s, d_sha0s,
                np.int32(num_visis), np.int32(num_baselines),
                block=(threads_x,1,1), grid=(grid_x, 1, 1))

        ##Call host code
        u, v, w, u_metres, v_metres, w_metres =  self._get_uvw_python(X_diff, Y_diff, Z_diff,
                                                                sdec0, cdec0, cha0s, sha0s,
                                                                wavelengths)


        ##Check the GPU function matches the python version and append to errors
        ##if not
        ##Check they match to a relative tolerance defined at top of tile, rtol
        ##Floating point error / whatever numpy is doing under the hood means
        ##they never match perfectly
        errors = []

        u_close = np.isclose(u, d_u.get(),rtol=rtol, atol=atol)
        v_close = np.isclose(v, d_v.get(),rtol=rtol, atol=atol)
        w_close = np.isclose(w, d_w.get(),rtol=rtol, atol=atol)

        u_m_close = np.isclose(u_metres, d_u_metres.get(),rtol=rtol, atol=atol)
        v_m_close = np.isclose(v_metres, d_v_metres.get(),rtol=rtol, atol=atol)
        w_m_close = np.isclose(w_metres, d_w_metres.get(),rtol=rtol, atol=atol)

        ##If you want to see the maximum rtol found, run
        ## $ pytest -s
        ##The -s option ensures the print statements are included
        u_ratio = abs((u / d_u.get()) - 1)
        v_ratio = abs((v / d_v.get()) - 1)
        w_ratio = abs((w / d_w.get()) - 1)
        print('\nu ratios rtol max', u_ratio.max())
        print('v ratios rtol max', v_ratio.max())
        print('w ratios rtol max', w_ratio.max())

        u_m_ratio = abs((u_metres / d_u_metres.get()) - 1)
        v_m_ratio = abs((v_metres / d_v_metres.get()) - 1)
        w_m_ratio = abs((w_metres / d_w_metres.get()) - 1)
        print('\nu_metres ratios rtol max', u_m_ratio.max())
        print('v_metres ratios rtol max', v_m_ratio.max())
        print('w_metres ratios rtol max', w_m_ratio.max())

        ##Make absolutely sure the device memory is being freed
        d_u.gpudata.free()
        d_v.gpudata.free()
        d_w.gpudata.free()
        d_X_diff.gpudata.free()
        d_Y_diff.gpudata.free()
        d_Z_diff.gpudata.free()
        d_u_metres.gpudata.free()
        d_v_metres.gpudata.free()
        d_w_metres.gpudata.free()

        # Build up a list of conditions we want to have met
        if False in u_close:
            errors.append("d_u returned by kern_calc_uvw is incorrect")
        if False in v_close:
            errors.append("d_v returned by kern_calc_uvw is incorrect")
        if False in w_close:
            errors.append("d_w returned by kern_calc_uvw is incorrect")
        if False in u_m_close:
            errors.append("d_u_metres returned by kern_calc_uvw is incorrect")
        if False in v_m_close:
            errors.append("d_v_metres returned by kern_calc_uvw is incorrect")
        if False in w_m_close:
            errors.append("d_w_metres returned by kern_calc_uvw is incorrect")

        # assert no error message has been occurred, else print messages
        assert not errors, "errors occured:\n{}".format("\n".join(errors))

if __name__ == '__main__':
    thing = TestFundamentalCoords()
    thing.test_kern_calc_uvw()
