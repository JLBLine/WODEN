import ctypes 
import importlib_resources
import wodenpy
import numpy as np

VELC = 299792458.0

# typedef struct _components_t {

#   //Instrinsic to COMPONENT values
#   double *ras; /*!< COMPONENT right ascensions (radians) */
#   double *decs; /*!< COMPONENT declinations (radians) */

#   //power law params
#   double *power_ref_freqs; /*!< COMPONENT Flux density reference frequencies (Hz) */
#   user_precision_t *power_ref_stokesI; /*!< COMPONENT Stokes I reference flux density (Jy) */
#   user_precision_t *power_ref_stokesQ; /*!< COMPONENT Stokes Q reference flux density (Jy) */
#   user_precision_t *power_ref_stokesU; /*!< COMPONENT Stokes U reference flux density (Jy) */
#   user_precision_t *power_ref_stokesV; /*!< COMPONENT Stokes V reference flux density (Jy) */
#   user_precision_t *power_SIs; /*!<  COMPONENT spectral indexes */

#   // curved power law params
#   double *curve_ref_freqs; /*!< COMPONENT Flux density reference frequencies (Hz) */
#   user_precision_t *curve_ref_stokesI; /*!< COMPONENT Stokes I reference flux density (Jy) */
#   user_precision_t *curve_ref_stokesQ; /*!< COMPONENT Stokes Q reference flux density (Jy) */
#   user_precision_t *curve_ref_stokesU; /*!< COMPONENT Stokes U reference flux density (Jy) */
#   user_precision_t *curve_ref_stokesV; /*!< COMPONENT Stokes V reference flux density (Jy) */
#   user_precision_t *curve_SIs; /*!<  COMPONENT spectral indexes */
#   user_precision_t *curve_qs; /*!<  COMPONENT curvature */

#   int *power_comp_inds; /*!< The indexes of all power-law models w.r.t ra,dec */
#   int *curve_comp_inds; /*!< The indexes of all curved power-law models w.r.t ra,dec */
#   int *list_comp_inds; /*!< The indexes of all list models w.r.t ra,dec */

#   //list flux params
#   double *list_freqs; /*!< COMPONENT Flux density references frequencies (Hz) */
#   user_precision_t *list_stokesI; /*!< COMPONENT Stokes I list flux density (Jy) */
#   user_precision_t *list_stokesQ; /*!< COMPONENT Stokes Q list flux density (Jy) */
#   user_precision_t *list_stokesU; /*!< COMPONENT Stokes U list flux density (Jy) */
#   user_precision_t *list_stokesV; /*!< COMPONENT Stokes V list flux density (Jy) */
#   int *num_list_values; /*!< How many freq/flux values are in each COMPONENT list*/
#   int *list_start_indexes; /*!< How many freq/flux values are in each COMPONENT list*/

#   int total_num_flux_entires; /*!< The total number of freq/flux values are in all lists combined*/

#   //something to store extrapolated output fluxes in
#   user_precision_t *extrap_stokesI; /*!< extrapolated COMPONENT Stokes I flux densities (Jy) */
#   user_precision_t *extrap_stokesQ; /*!< extrapolated COMPONENT Stokes Q flux densities (Jy) */
#   user_precision_t *extrap_stokesU; /*!< extrapolated COMPONENT Stokes U flux densities (Jy) */
#   user_precision_t *extrap_stokesV; /*!< extrapolated COMPONENT Stokes V flux densities (Jy) */

#   //SHAPELET / GAUSSIAN params
#   user_precision_t *shape_coeffs; /*!< Scaling coefficients for SHAPELET basis functions */
#   user_precision_t *n1s; /*!< 1st basis function order for SHAPELET basis functions */
#   user_precision_t *n2s; /*!< 2nd basis function order for SHAPELET basis functions */
#   user_precision_t *majors; /*!< GAUSSIAN/SHAPELET major axis (beta1, radians) */
#   user_precision_t *minors; /*!< GAUSSIAN/SHAPELET minor axis (beta2, radians) */
#   user_precision_t *pas; /*!< GAUSSIAN/SHAPELET position angles (radians) */
#   user_precision_t *param_indexes; /*!< An index value to match each coeff, n1, and n2
#   to the correct ra, dec, major, minor, pa for a SHAPELET */

#   //Specific to observation settings for these COMPONENTs
#   user_precision_t *azs; /*!< SHAPELET source azimuth angles for all time steps */
#   user_precision_t *zas; /*!< SHAPELET source zenith angles for all time steps */
#   double *beam_has; /*!< Hour angle of COMPONENTs for all time steps, used for
#    beam calculations */
#   double *beam_decs; /*!< Declinations of COMPONENTs for all time steps, used for
#    beam calculations */
#   int num_primarybeam_values; /*!< Number of beam calculations needed for
#   COMPONENTs */
#   //float *angular_seps; /*!< Angular separation of each COMPONENTS from
#   //zenith, for all time steps. Used in MWA analytic beam  */

#   /*
#   These beam values are two dimensional, in anticipation of one day having
#   different primary beams for different tiles. At the moment everything is
#   the same for each tile, so only have length one for first dimension
#   */
#   user_precision_complex_t *gxs; /*!< North-South Beam gain values for all directions,
#   frequencies, and times for these COMPONENTS*/
#   user_precision_complex_t *Dxs; /*!< North-South Beam leakage values for all directions,
#   frequencies, and times for these COMPONENTS*/
#   user_precision_complex_t *Dys; /*!< East-West Beam leakage values for all directions,
#   frequencies, and times for these COMPONENTS*/
#   user_precision_complex_t *gys; /*!< East-West Beam gain values for all directions,
#   frequencies, and times for these COMPONENTS*/

#   //Leave off the d_ from these device values, as the components_t struct
#   //itself will have the d_ label if doing things on the GPU
#   double *ls; /*!< Device memory l cosine direction coords for these COMPONENTs*/
#   double *ms; /*!< Device memory m cosine direction coords for these COMPONENTs*/
#   double *ns; /*!< Device memory n cosine direction coords for these COMPONENTs*/

# } components_t;



# /*!
# A struct to contain sky model values for a single SOURCE
# */
# typedef struct _source_t {
#   //General source info
#   char name[32]; /*!< Source name */
#   int n_comps; /*!< Total number of COMPONENTs in source  */

#   int n_points; /*!< Number of POINT source COMPONENTs  */
#   int n_point_lists; /*!< Number of POINTs with LIST type flux */
#   int n_point_powers; /*!< Number of POINTs with POWER_LAW type flux */
#   int n_point_curves; /*!< Number of POINTs with CURVED_POWER_LAW type flux */

#   int n_gauss; /*!< Number of GAUSSIAN source COMPONENTs */
#   int n_gauss_lists; /*!< Number of GAUSSIANs with LIST type flux */
#   int n_gauss_powers; /*!< Number of GAUSSIANs with POWER_LAW type flux */
#   int n_gauss_curves; /*!< Number of GAUSSIANs with CURVED_POWER_LAW type flux */

#   int n_shapes; /*!< Number of SHAPELET source COMPONENTs */
#   int n_shape_lists; /*!< Number of SHAPELETs with LIST type flux */
#   int n_shape_powers; /*!< Number of SHAPELETs with POWER_LAW type flux */
#   int n_shape_curves; /*!< Number of SHAPELETs with CURVED_POWER_LAW type flux */
#   int n_shape_coeffs; /*!< Total number of SHAPELET coefficients */

#   components_t point_components; /*!< `components_t` holding component
#   information for all POINT COMPONENTs in this SOURCE.*/
#   components_t gauss_components; /*!< `components_t` holding component
#   information for all GAUSSIAN COMPONENTs in this SOURCE.*/
#   components_t shape_components; /*!< `components_t` holding component
#   information for all SHAPELET COMPONENTs in this SOURCE.*/

#   //Device versions
#   components_t d_point_components; /*!< `components_t` holding component
#   information for all POINT COMPONENTs in this SOURCE.*/
#   components_t d_gauss_components; /*!< `components_t` holding component
#   information for all GAUSSIAN COMPONENTs in this SOURCE.*/
#   components_t d_shape_components; /*!< `components_t` holding component
#   information for all SHAPELET COMPONENTs in this SOURCE.*/

# } source_t;

class Components_Double(ctypes.Structure):
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
    
class Components_Float(ctypes.Structure):
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
    
def setup_source(chunk_map : Skymodel_Chunk_Map, precision='double') -> ctypes.Structure:
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