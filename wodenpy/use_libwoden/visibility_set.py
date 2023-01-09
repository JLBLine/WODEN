import ctypes 
import importlib_resources
import wodenpy
import numpy as np

VELC = 299792458.0

class Visi_Set_Double(ctypes.Structure):
    """A class structured equivalently to a `visi_set` struct, used by 
    the C and CUDA code in libwoden_double.so
    """
    
    _fields_ = [("us_metres", ctypes.POINTER(ctypes.c_double)),
                ("vs_metres", ctypes.POINTER(ctypes.c_double)),
                ("ws_metres", ctypes.POINTER(ctypes.c_double)),
                ("allsteps_sha0s", ctypes.POINTER(ctypes.c_double)),
                ("allsteps_cha0s", ctypes.POINTER(ctypes.c_double)),
                ("allsteps_lsts", ctypes.POINTER(ctypes.c_double)),
                ("allsteps_wavelengths", ctypes.POINTER(ctypes.c_double)),
                ("channel_frequencies", ctypes.POINTER(ctypes.c_double)),
                ("sum_visi_XX_real", ctypes.POINTER(ctypes.c_double)),
                ("sum_visi_XX_imag", ctypes.POINTER(ctypes.c_double)),
                ("sum_visi_XY_real", ctypes.POINTER(ctypes.c_double)),
                ("sum_visi_XY_imag", ctypes.POINTER(ctypes.c_double)),
                ("sum_visi_YX_real", ctypes.POINTER(ctypes.c_double)),
                ("sum_visi_YX_imag", ctypes.POINTER(ctypes.c_double)),
                ("sum_visi_YY_real", ctypes.POINTER(ctypes.c_double)),
                ("sum_visi_YY_imag", ctypes.POINTER(ctypes.c_double))]
    
class Visi_Set_Float(ctypes.Structure):
    """A class structured equivalent to a `visi_set` struct, used by 
    the C and CUDA code in libwoden_float.so
    """
    
    _fields_ = [("us_metres", ctypes.POINTER(ctypes.c_float)),
                ("vs_metres", ctypes.POINTER(ctypes.c_float)),
                ("ws_metres", ctypes.POINTER(ctypes.c_float)),
                ("allsteps_sha0s", ctypes.POINTER(ctypes.c_float)),
                ("allsteps_cha0s", ctypes.POINTER(ctypes.c_float)),
                ("allsteps_lsts", ctypes.POINTER(ctypes.c_float)),
                ("allsteps_wavelengths", ctypes.POINTER(ctypes.c_float)),
                ("channel_frequencies", ctypes.POINTER(ctypes.c_float)),
                ("sum_visi_XX_real", ctypes.POINTER(ctypes.c_float)),
                ("sum_visi_XX_imag", ctypes.POINTER(ctypes.c_float)),
                ("sum_visi_XY_real", ctypes.POINTER(ctypes.c_float)),
                ("sum_visi_XY_imag", ctypes.POINTER(ctypes.c_float)),
                ("sum_visi_YX_real", ctypes.POINTER(ctypes.c_float)),
                ("sum_visi_YX_imag", ctypes.POINTER(ctypes.c_float)),
                ("sum_visi_YY_real", ctypes.POINTER(ctypes.c_float)),
                ("sum_visi_YY_imag", ctypes.POINTER(ctypes.c_float))]
    
def setup_visi_set(num_visis : int, precision='double') -> ctypes.Structure:
    """Sets up a ctypes structure class to contain the visibility outputs.
    This class is compatible with the C/CUDA code, and will allocate the
    correct amount of memory, based on whether the precision is either
    'double' or 'float'.
    
    Parameters
    ----------
    num_visis : int
        Number of visibilities to assign memory for
    precision : str, optional
        Precision to be used, either 'float' or 'double. Defaults to 'double'

    Returns
    -------
    visibility_set : ctypes.Structure
        An initialised `wodenpy.use_libwoden.use_ctypes.Visi_Set_Float` or
        `wodenpy.use_libwoden.use_ctypes.Visi_Set_Double` class,
        compatible with libwoden_float.so or libwoden_double.so.
    """
    
    if precision == 'float':
    
        visibility_set = Visi_Set_Float()
        num_visi_array = ctypes.c_float*num_visis
        
    else:
        visibility_set = Visi_Set_Double()
        num_visi_array = ctypes.c_double*num_visis
        
    visibility_set.us_metres = num_visi_array()
    visibility_set.vs_metres = num_visi_array()
    visibility_set.ws_metres = num_visi_array()
    visibility_set.sum_visi_XX_real = num_visi_array()
    visibility_set.sum_visi_XX_imag = num_visi_array()
    visibility_set.sum_visi_XY_real = num_visi_array()
    visibility_set.sum_visi_XY_imag = num_visi_array()
    visibility_set.sum_visi_YX_real = num_visi_array()
    visibility_set.sum_visi_YX_imag = num_visi_array()
    visibility_set.sum_visi_YY_real = num_visi_array()
    visibility_set.sum_visi_YY_imag = num_visi_array()
    
    return visibility_set

def load_visibility_set(visibility_set=None,num_baselines=None,num_freq_channels=None,
              num_time_steps=None, precision=None, do_autos=False, num_ants=0):
    """
    Read the WODEN ctype output and shove into a numpy arrays, ready to be put
    into a uvfits file. By default, WODEN only outputs cross-correlations.
    In this case, the output binary is ordered by baseline (fastest
    changing), frequency, and time (slowest changing). Visibility coords and
    data are read in, with the visi data output into an array of
    `shape=(num_time_steps*num_baselines,1,1,num_freq_channels,4,3))`, which is
    appropriate for a uvfits file. Needs to know whether WODEN was run with
    'float' (32 bit) or 'double' (64 bit) precision to read in the data
    correctly.

    If WODEN was run with `do_autos=True`, then auto correlations are
    also in the output binary. These are stored AFTER the cross-correlations,
    and orders by antenna (fastest changing), frequency, and time (slowest changing).
    In this case the data are output into an array of
    `shape=(num_time_steps*(num_baselines+num_ants),1,1,num_freq_channels,4,3))`,
    where we say a baseline is only defined between to different antennas.
    The visibilities are output to match the BASELINE array, which orders
    the autos and crosses via antenna pairs as (1,1), (1,2), (1,3) .. (2,2),
    (2,3) etc etc meaning the autos and crosses are mixed.

    Parameters
    ----------
    filename : string
        Name of WODEN binary file to read from
    num_baselines : int
        Number of baselines in the binary file
    num_freq_channels : int
        Number of frequencies in the binary file
    num_time_steps : int
        Number of time steps in the binary file
    precision : string
        Precision WODEN was run with - either 'float' or 'double'
    do_autos : Boolean
        if True, data has auto-correlations in
    do_autos : int
        Number of antennas in the array

    Returns
    -------
    uus : float array
        The :math:`u` coordinates (seconds). These are zero for auto-correlations.
    vvs : float array
        The :math:`v` coordinates (seconds). These are zero for auto-correlations.
    wws : float array
        The :math:`w` coordinates (seconds). These are zero for auto-correlations.
    v_container : float array
        Visibility data with
        `shape=(num_time_steps*num_baselines,1,1,num_freq_channels,4,3))`
    """

    ##If not doing autos, ensure this number is zero
    if do_autos == False:
        num_ants = 0

    n_data = num_time_steps * (num_baselines + num_ants)
    uus = np.zeros(n_data)
    vvs = np.zeros(n_data)
    wws = np.zeros(n_data)
    v_container = np.zeros((n_data,1,1,num_freq_channels,4,3))

    num_visi = num_time_steps * num_freq_channels * (num_baselines + num_ants)
    num_cross = num_time_steps * num_freq_channels * num_baselines

    ##Righto, this converts from the ctype POINTER into a numpy array
    ##This is grabbing all the lovely things calculated by the GPU
    us_metres = np.ctypeslib.as_array(visibility_set.us_metres, shape=(num_visi,))
    vs_metres = np.ctypeslib.as_array(visibility_set.vs_metres, shape=(num_visi,))
    ws_metres = np.ctypeslib.as_array(visibility_set.ws_metres, shape=(num_visi,))
    visi_XX_real = np.ctypeslib.as_array(visibility_set.sum_visi_XX_real, shape=(num_visi,))
    visi_XX_imag = np.ctypeslib.as_array(visibility_set.sum_visi_XX_imag, shape=(num_visi,))
    visi_XY_real = np.ctypeslib.as_array(visibility_set.sum_visi_XY_real, shape=(num_visi,))
    visi_XY_imag = np.ctypeslib.as_array(visibility_set.sum_visi_XY_imag, shape=(num_visi,))
    visi_YX_real = np.ctypeslib.as_array(visibility_set.sum_visi_YX_real, shape=(num_visi,))
    visi_YX_imag = np.ctypeslib.as_array(visibility_set.sum_visi_YX_imag, shape=(num_visi,))
    visi_YY_real = np.ctypeslib.as_array(visibility_set.sum_visi_YY_real, shape=(num_visi,))
    visi_YY_imag = np.ctypeslib.as_array(visibility_set.sum_visi_YY_imag, shape=(num_visi,))

    ##If doing auto-correlations, need some mapping arrays so we can
    ##shove the correct data into the correct spots
    if do_autos:
        cross_map = []
        auto_map = []

        visi_map = 0
        for b1 in np.arange(num_ants):
            for b2 in np.arange(b1,num_ants):
                if b1 == b2:
                    auto_map.append(visi_map)
                else:
                    cross_map.append(visi_map)
                visi_map += 1

        cross_map = np.array(cross_map, dtype=int)
        auto_map = np.array(auto_map, dtype=int)

    ##We only care about the u,v,w for the cross-correlations, so fill
    ##them only
    for time_ind in np.arange(num_time_steps):

        time_step = num_baselines * time_ind * num_freq_channels

        ##TODO now that we are reading the u,v,w directly from memory out of
        ##GPU, we don't have to store copies of u,v,w for every frequency like
        ##we did when writing to file. Will be a moderate save on GPU memory
        if do_autos:
            time_base = int(time_ind*(num_baselines + num_ants))
            this_cross_map = cross_map + time_base

            ##The baseline length
            uus[this_cross_map] = us_metres[time_step:time_step+num_baselines] / VELC
            vvs[this_cross_map] = vs_metres[time_step:time_step+num_baselines] / VELC
            wws[this_cross_map] = ws_metres[time_step:time_step+num_baselines] / VELC

        else:

            uus[time_ind*num_baselines:(time_ind + 1)*num_baselines] = us_metres[time_step:time_step+num_baselines] / VELC
            vvs[time_ind*num_baselines:(time_ind + 1)*num_baselines] = vs_metres[time_step:time_step+num_baselines] / VELC
            wws[time_ind*num_baselines:(time_ind + 1)*num_baselines] = ws_metres[time_step:time_step+num_baselines] / VELC

    for time_ind in np.arange(num_time_steps):
        for freq_ind in np.arange(num_freq_channels):

            freq_step = num_baselines * (time_ind * num_freq_channels + freq_ind)

            cross_XX_re = visi_XX_real[freq_step:freq_step+num_baselines]
            cross_XX_im = visi_XX_imag[freq_step:freq_step+num_baselines]
            cross_YY_re = visi_YY_real[freq_step:freq_step+num_baselines]
            cross_YY_im = visi_YY_imag[freq_step:freq_step+num_baselines]
            cross_XY_re = visi_XY_real[freq_step:freq_step+num_baselines]
            cross_XY_im = visi_XY_imag[freq_step:freq_step+num_baselines]
            cross_YX_re = visi_YX_real[freq_step:freq_step+num_baselines]
            cross_YX_im = visi_YX_imag[freq_step:freq_step+num_baselines]

            ##If doing auto-correlations, load up the autos and do some fancy
            ##mapping
            if do_autos:

                time_base = int(time_ind*(num_baselines + num_ants))
                this_cross_map = cross_map + time_base

                v_container[this_cross_map,0,0,freq_ind,0,0] = cross_XX_re
                v_container[this_cross_map,0,0,freq_ind,0,1] = cross_XX_im
                v_container[this_cross_map,0,0,freq_ind,1,0] = cross_YY_re
                v_container[this_cross_map,0,0,freq_ind,1,1] = cross_YY_im
                v_container[this_cross_map,0,0,freq_ind,2,0] = cross_XY_re
                v_container[this_cross_map,0,0,freq_ind,2,1] = cross_XY_im
                v_container[this_cross_map,0,0,freq_ind,3,0] = cross_YX_re
                v_container[this_cross_map,0,0,freq_ind,3,1] = cross_YX_im

                freq_step = num_ants * (time_ind * num_freq_channels + freq_ind)

                real_XX_ind = num_cross + freq_step
                imag_XX_ind = num_cross + freq_step
                real_YY_ind = num_cross + freq_step
                imag_YY_ind = num_cross + freq_step
                real_XY_ind = num_cross + freq_step
                imag_XY_ind = num_cross + freq_step
                real_YX_ind = num_cross + freq_step
                imag_YX_ind = num_cross + freq_step

                auto_XX_re = visi_XX_real[real_XX_ind:real_XX_ind+num_ants]
                auto_XX_im = visi_XX_imag[imag_XX_ind:imag_XX_ind+num_ants]
                auto_YY_re = visi_YY_real[real_YY_ind:real_YY_ind+num_ants]
                auto_YY_im = visi_YY_imag[imag_YY_ind:imag_YY_ind+num_ants]
                auto_XY_re = visi_XY_real[real_XY_ind:real_XY_ind+num_ants]
                auto_XY_im = visi_XY_imag[imag_XY_ind:imag_XY_ind+num_ants]
                auto_YX_re = visi_YX_real[real_YX_ind:real_YX_ind+num_ants]
                auto_YX_im = visi_YX_imag[imag_YX_ind:imag_YX_ind+num_ants]

                this_auto_map = auto_map + time_base

                v_container[this_auto_map,0,0,freq_ind,0,0] = auto_XX_re
                v_container[this_auto_map,0,0,freq_ind,0,1] = auto_XX_im
                v_container[this_auto_map,0,0,freq_ind,1,0] = auto_YY_re
                v_container[this_auto_map,0,0,freq_ind,1,1] = auto_YY_im
                v_container[this_auto_map,0,0,freq_ind,2,0] = auto_XY_re
                v_container[this_auto_map,0,0,freq_ind,2,1] = auto_XY_im
                v_container[this_auto_map,0,0,freq_ind,3,0] = auto_YX_re
                v_container[this_auto_map,0,0,freq_ind,3,1] = auto_YX_im

            ##Otherwise, everything is a cross-correlation so just bung em in
            else:
                v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,0,0] = cross_XX_re
                v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,0,1] = cross_XX_im
                v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,1,0] = cross_YY_re
                v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,1,1] = cross_YY_im
                v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,2,0] = cross_XY_re
                v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,2,1] = cross_XY_im
                v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,3,0] = cross_YX_re
                v_container[time_ind*num_baselines:(time_ind + 1)*num_baselines,0,0,freq_ind,3,1] = cross_YX_im

    ##Set the weights for everything to one
    v_container[:,0,0,:,:,2] = 1.0

    return uus, vvs, wws, v_container