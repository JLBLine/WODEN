/*! \file
  Device methods to calculate interferometric coorindates \f$u,v,w\f$ and
  \f$l,m,n\f$.
*/

#include <stdio.h>
#include <stdlib.h>
#include "gpu_macros.h"
#include <complex.h>
#include <math.h>
#include "constants.h"
#include "woden_precision_defs.h"
#include "woden_struct_defs.h"

/**
@brief Use the given baseline lengths in local `X,Y,Z` coords, and a given
phase centre, calculate `u,v,w` coordinates

@details Performs Equation 4.1 from TMS
<https://link.springer.com/book/10.1007/978-3-319-44431-4>. This `__device__`
function is called by `kern_calc_uvw`, which calculates `u,v,w` for multiple
time steps. Using the index `iBaseline`, `num_baselines`, `num_times`, `num_freqs`,
this function does a modulus to index `d_X_diff`, `d_Y_diff`, `d_Z_diff` correctly,
as the diff arrays do not change with frequency (but change with time).

@param[in] d_X_diff Baseline length in the `X` direction (metres)
@param[in] d_Y_diff Baseline length in the `Y` direction (metres)
@param[in] d_Z_diff Baseline length in the `Z` direction (metres)
@param[in] sdec0 Sine of the declination of the phase centre
@param[in] cdec0 Cosine of the declination of the phase centre
@param[in] sha0 Sine of the hour angle of the phase centre
@param[in] cha0 Cosine of the hour angle of the phase centre
@param[in] iBaseline Index passed from `kern_calc_uvw`
@param[in] num_baselines Number of baselines for a single time step
@param[in] num_times Number of time steps
@param[in] num_freqs Number of frequency steps
@param[in,out] u Output `u` coord (metres)
@param[in,out] v Output `v` coord (metres)
@param[in,out] w Output `w` coord (metres)
*/
__device__ void calc_uvw(double *d_X_diff, double *d_Y_diff,
                         double *d_Z_diff,
                         double sdec0, double cdec0,
                         double sha0, double cha0,
                         int iBaseline, int num_baselines, int num_times,
                         int num_freqs,
                         user_precision_t * u, user_precision_t * v,
                         user_precision_t * w);

/**
@brief Calculate the \f$u,v,w\f$ in metres and in wavelengths for multiple
time steps and frequencies

@details `d_X_diff`, `d_Y_diff`, `d_Z_diff` contain baseline lengths in local
`X,Y,Z` coords, which can change with time (when precessed back to J2000),
so should be of length `num_baselines`*`num_times`. To keep indexing simple, the wavelengths in `d_wavelengths` and hour angles in `d_cha0s, d_sha0s` should contain
values for all baselines, all frequencies, and all time steps. .

When called with `dim3 grid, threads`, kernel should be called with `grid.x`
 set, where:
 - grid.x * threads.x >= `num_visis`

@param[in] d_X_diff Baseline lengths in the `X` direction (metres)
@param[in] d_Y_diff Baseline lengths in the `Y` direction (metres)
@param[in] d_Z_diff Baseline lengths in the `Z` direction (metres)
@param[in,out] d_u_metres Output `u` coords (metres)
@param[in,out] d_v_metres Output `v` coords (metres)
@param[in,out] d_w_metres Output `w` coords (metres)
@param[in,out] d_u Output `u` coord (wavelengths)
@param[in,out] d_v Output `v` coord (wavelengths)
@param[in,out] d_w Output `w` coord (wavelengths)
@param[in] d_wavelengths Wavelengths for all baselines, frequencies, and time
steps (metres)
@param[in] sdec0 Sine of the declination of the phase centre
@param[in] cdec0 Cosine of the declination of the phase centre
@param[in] d_cha0s Cosine of the hour angle of the phase centre for all
baselines, frequencies, and time steps
@param[in] d_sha0s Sine of the hour angle of the phase centre for all
baselines, frequencies, and time steps
@param[in] num_cross Total number of `u,v,w` coords to be calculated (number of cross correlations)
@param[in] num_baselines Number of baselines for a single time step
@param[in] num_times Number of time steps
@param[in] num_freqs Number of frequency steps
*/
__global__ void kern_calc_uvw(double *d_X_diff, double *d_Y_diff,
           double *d_Z_diff, user_precision_t *d_u_metres,
           user_precision_t *d_v_metres, user_precision_t *d_w_metres,
           user_precision_t *d_u, user_precision_t *d_v, user_precision_t *d_w, user_precision_t *d_wavelengths,
           double sdec0, double cdec0,
           double *d_cha0s, double *d_sha0s,
           int num_cross, int num_baselines, int num_times, int num_freqs);

/**
@brief This is just a wrapper to launch `kern_calc_uvw` to calculate `u,v,w` coords
on the GPU. Wrapper means we can call this function from the CPU. All arrays
should be initialised and filled on the GPU before calling this function.

@param[in] d_X_diff Baseline lengths in the `X` direction (metres)
@param[in] d_Y_diff Baseline lengths in the `Y` direction (metres)
@param[in] d_Z_diff Baseline lengths in the `Z` direction (metres)
@param[in,out] d_u_metres Output `u` coords (metres)
@param[in,out] d_v_metres Output `v` coords (metres)
@param[in,out] d_w_metres Output `w` coords (metres)
@param[in,out] d_us Output `u` coord (wavelengths)
@param[in,out] d_vs Output `v` coord (wavelengths)
@param[in,out] d_ws Output `w` coord (wavelengths)
@param[in] d_allsteps_wavelengths Wavelengths for all baselines, frequencies, and time
steps (metres)
@param[in] d_allsteps_cha0s Cosine of the hour angle of the phase centre for all
baselines, frequencies, and time steps
@param[in] d_allsteps_sha0s Sine of the hour angle of the phase centre for all
baselines, frequencies, and time steps
@param[in] woden_settings Struct containing simulation settings
*/
extern "C" void calc_uvw_gpu(double *d_X_diff, double *d_Y_diff, double *d_Z_diff,
                             user_precision_t *d_u_metres,
                             user_precision_t *d_v_metres, user_precision_t *d_w_metres,
                             user_precision_t *d_us, user_precision_t *d_vs,
                             user_precision_t *d_ws, user_precision_t *d_allsteps_wavelengths,
                             double *d_allsteps_cha0s, double *d_allsteps_sha0s,
                             woden_settings_t *woden_settings);


/**
@brief The SHAPELET visibility envelope calculation uses a `u,v` coordinate
system where the phase centre is set to the RA/Dec of that particular
SHAPELET COMPONENT. This kernel calculates multiple \f$u,v\f$ coodinates
with varying phase centres to be used later by
`source_components.kern_calc_visi_shapelets`.

@details The output arrays `d_u_shapes`, `d_v_shapes`,
contain coordinates for all baselines (fastest changing), all times, and
all SHAPELET components (slowest changing), in units of metres.

When called with `dim3 grid, threads`, kernel should be called with both
`grid.x` and `grid.y` set, where:
 - grid.x * threads.x >= `num_baselines*num_times`
 - grid.y * threads.y >= `num_shapes`

@param[in] d_X_diff Baseline lengths in the `X` direction (metres)
@param[in] d_Y_diff Baseline lengths in the `Y` direction (metres)
@param[in] d_Z_diff Baseline lengths in the `Z` direction (metres)
@param[in,out] d_u_shapes Output `u` coords with various phase centres for
SHAPELET components (metres)
@param[in,out] d_v_shapes Output `v` coords with various phase centres for
SHAPELET components (metres)
@param[in] d_lsts The local sidereal times of all time steps in simulation
(radians)
@param[in] d_ras Array of SHAPELET Right Ascensions (radians)
@param[in] d_decs Array of SHAPELET Declinations (radians)
@param[in] num_baselines Number of baselines for a single time step
@param[in] num_times Number of time steps
@param[in] num_shapes Number of SHAPELET COMPONENTs

*/
__global__ void kern_calc_uv_shapelet(double *d_X_diff,
      double *d_Y_diff, double *d_Z_diff,
      user_precision_t *d_u_shapes, user_precision_t *d_v_shapes,
      double *d_lsts, double *d_ras, double *d_decs,
      const int num_baselines, const int num_times, const int num_shapes);


/**
@brief This is just a wrapper to launch `kern_calc_uv_shapelet` to calculate
`u,v` coords on the GPU. Wrapper means we can call this function from the CPU.

@details All `d_*` arrays should be initialised and filled on the GPU before calling this.
The function will copy `woden_settings->lsts` to the GPU, use it, then free it.
Must also have `woden_settings->num_baselines`, `woden_settings->num_time_steps`
intialised.

@param[in,out] d_u_shapes Output `u` coords with various phase centres for SHAPELET components (wavelengths)
@param[in,out] d_v_shapes Output `v` coords with various phase centres for SHAPELET components (wavelengths)
@param[in] num_shapes Number of SHAPELET COMPONENTs
@param[in] d_X_diff Baseline lengths in the `X` direction (metres)
@param[in] d_Y_diff Baseline lengths in the `Y` direction (metres)
@param[in] d_Z_diff Baseline lengths in the `Z` direction (metres)
@param[in] d_ras Array of SHAPELET Right Ascensions (radians)
@param[in] d_decs Array of SHAPELET Declinations (radians)
@param[in] woden_settings Struct containing simulation settings
 */
extern "C" void calc_uv_shapelet_gpu(user_precision_t *d_u_shapes, user_precision_t *d_v_shapes,
                                    int num_shapes,
                                    double *d_X_diff, double *d_Y_diff, double *d_Z_diff,
                                    double *d_ras, double *d_decs,
                                    woden_settings_t *woden_settings);


/**
@brief Calculate interferometric \f$l,m,n\f$ image coords for a given RA,Dec
and phase centre

@details Here it is assumed the phase centre is not changing with time

@param[in] ra0 Right Ascension of the phase centre (radians)
@param[in] sdec0 Sine of Declination of the phase centre
@param[in] cdec0 Cosine of Declination of the phase centre
@param[in] ra Right Ascension (radians)
@param[in] dec Declination (radians)
@param[in] l Output \f$l\f$ coodinate
@param[in] m Output \f$m\f$ coodinate
@param[in] n Output \f$n\f$ coodinate

*/
__device__ void calc_lmn(double ra0, double sdec0,
                         double cdec0,
                         double ra, double dec,
                         double * l, double * m, double * n);

/**
@brief Calculate interferometric \f$l,m,n\f$ image coords for a set of
RA,Dec coodinates, with a single RA/Dec phase centre.

@details When called with `dim3 grid, threads`, kernel should be called with
`grid.x` set, where:
 - grid.x * threads.x >= `num_components`

 @param[in] ra0 Right Ascension of the phase centre (radians)
 @param[in] sdec0 Sine of Declination of the phase centre
 @param[in] cdec0 Cosine of Declination of the phase centre
 @param[in] d_ras Array of Right Ascensions (radians)
 @param[in] d_decs Array of Declinations (radians)
 @param[in,out] d_l Output \f$l\f$ coodinates
 @param[in,out] d_m Output \f$m\f$ coodinates
 @param[in,out] d_n Output \f$n\f$ coodinates
 @param[in] num_components Number of RA,Dec coords
*/
__global__ void kern_calc_lmn(double ra0, double sdec0,
                              double cdec0,
                              double *d_ras, double *d_decs,
                              double *d_l, double *d_m, double *d_n,
                              int num_components);


/**
@brief Calculates the \f$l,m,n\f$ image coords for an intialised `d_components`
struct. Assumes d_components->ras, d_components->decs have been allocated
and set on the device.

@details Runs `kern_calc_lmn` using d_components->ras, d_components->decs and
woden_settings->ra0, woden_settings->sdec0, woden_settings->cdec0 as inputs.
This function allocates device memory for `d_components->ls`, `d_components->ms`,
and `d_components->ns`, and puts the results in these arrays.

 @param[in,out] d_components Initialised `components_t` struct containing the
 ra and dec arrays
 @param[in] num_components Number of RA,Dec coords in `d_components`
 @param[in] woden_settings Woden settings struct containing the phase centre
*/
extern "C" void calc_lmn_for_components_gpu(components_t *d_components,
                                           int num_components,
                                           woden_settings_t *woden_settings);


/**
 @brief Internally to WODEN C/GPU, all cross-correlations are stored first,
 then the autos after. All autos have zero length, so use this function to
 set all autos uvw to zero.

 @details `num_cross` is the number of cross-correlations, so fill any uvw after
 this by adding `iAuto = threadIdx.x + (blockDim.x*blockIdx.x)`, meaning you
 should run this kernel with `grid.x` set, where:
  - grid.x * threads.x >= `num_autos`

 @param[in] num_cross Number of cross-correlations
 @param[in] num_autos Number of auto-correlations
 @param[in,out] d_u u coords
 @param[in,out] d_v v coords
 @param[in,out] d_w w coords
*/
__global__ void kern_set_auto_uvw_to_zero(int num_cross, int num_autos,
                                     user_precision_t *d_u,
                                     user_precision_t *d_v,
                                     user_precision_t *d_w);


/**
 @brief This is just a wrapper to launch `kern_set_auto_uvw_to_zero` to set
 all auto-correlations to zero. Wrapper means we can call this function from
 the CPU.

 @details All `d_*` arrays should be initialised and filled on the GPU before
 calling this.

 @param[in] num_cross Number of cross-correlations
 @param[in] num_autos Number of auto-correlations
 @param[in,out] d_u u coords
 @param[in,out] d_v v coords
 @param[in,out] d_w w coords
*/
extern "C" void set_auto_uvw_to_zero_gpu(int num_cross, int num_autos,
                                     user_precision_t *d_u,
                                     user_precision_t *d_v,
                                     user_precision_t *d_w);