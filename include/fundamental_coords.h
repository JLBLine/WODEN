/*! \file
  Device methods to calculate interferometric coorindates \f$u,v,w\f$ and
  \f$l,m,n\f$.
*/

#include <math.h>
#include <stdint.h>
#include "woden_precision_defs.h"

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
