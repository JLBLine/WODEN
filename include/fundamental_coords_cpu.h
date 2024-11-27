/*! \file
  Device methods to calculate interferometric coorindates \f$u,v,w\f$ and
  \f$l,m,n\f$ on the CPU.
*/

#include <math.h>
#include <stdint.h>
#include "woden_precision_defs.h"

/**
@brief Calculate the \f$u,v,w\f$ in metres and in wavelengths for multiple
time steps and frequencies

@details `X_diff`, `Y_diff`, `Z_diff` contain baseline lengths in local
`X,Y,Z` coords, which can change with time (when precessed back to J2000),
so should be of length `num_baselines`*`num_times`. To keep indexing simple, the wavelengths in `wavelengths` and hour angles in `cha0s, sha0s` should contain
values for all baselines, all frequencies, and all time steps. .

@param[in] X_diff Baseline lengths in the `X` direction (metres)
@param[in] Y_diff Baseline lengths in the `Y` direction (metres)
@param[in] Z_diff Baseline lengths in the `Z` direction (metres)
@param[in,out] u_metres Output `u` coords (metres)
@param[in,out] v_metres Output `v` coords (metres)
@param[in,out] w_metres Output `w` coords (metres)
@param[in,out] us Output `u` coord (wavelengths)
@param[in,out] vs Output `v` coord (wavelengths)
@param[in,out] ws Output `w` coord (wavelengths)
@param[in] wavelengths Wavelengths for all baselines, frequencies, and time
steps (metres)
@param[in] sdec0 Sine of the declination of the phase centre
@param[in] cdec0 Cosine of the declination of the phase centre
@param[in] cha0s Cosine of the hour angle of the phase centre for all
baselines, frequencies, and time steps
@param[in] sha0s Sine of the hour angle of the phase centre for all
baselines, frequencies, and time steps
@param[in] num_cross Total number of `u,v,w` coords to be calculated (number of cross correlations)
@param[in] num_baselines Number of baselines for a single time step
@param[in] num_times Number of time steps
@param[in] num_freqs Number of frequency steps
*/
void calc_uvw_cpu(double *X_diff, double *Y_diff, double *Z_diff,
                  user_precision_t *u_metres, user_precision_t *v_metres, user_precision_t *w_metres,
                  user_precision_t *us, user_precision_t *vs, user_precision_t *ws,
                  user_precision_t *wavelengths,
                  double sdec0, double cdec0,
                  double *cha0s, double *sha0s,
                  int num_cross, int num_baselines, int num_times, int num_freqs);

/**
@brief The SHAPELET visibility envelope calculation uses a `u,v` coordinate
system where the phase centre is set to the RA/Dec of that particular
SHAPELET COMPONENT. This kernel calculates multiple \f$u,v\f$ coodinates
with varying phase centres to be used later by
`source_components.kern_calc_visi_shapelets`.

@details The output arrays `u_shapes`, `v_shapes`,
contain coordinates for all baselines (fastest changing), all times, and
all SHAPELET components (slowest changing), in units of metres.

@param[in] X_diff Baseline lengths in the `X` direction (metres)
@param[in] Y_diff Baseline lengths in the `Y` direction (metres)
@param[in] Z_diff Baseline lengths in the `Z` direction (metres)
@param[in,out] u_shapes Output `u` coords with various phase centres for
SHAPELET components (metres)
@param[in,out] v_shapes Output `v` coords with various phase centres for
SHAPELET components (metres)
@param[in] lsts The local sidereal times of all time steps in simulation
(radians)
@param[in] ras Array of SHAPELET Right Ascensions (radians)
@param[in] decs Array of SHAPELET Declinations (radians)
@param[in] num_baselines Number of baselines for a single time step
@param[in] num_times Number of time steps
@param[in] num_shapes Number of SHAPELET COMPONENTs

*/
void calc_uv_shapelet_cpu(double *X_diff, double *Y_diff, double *Z_diff,
      user_precision_t *u_shapes, user_precision_t *v_shapes,
      double *lsts, double *ras, double *decs,
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
@param[in] ls Output \f$l\f$ coodinate
@param[in] ms Output \f$m\f$ coodinate
@param[in] ns Output \f$n\f$ coodinate

*/
void calc_lmn_cpu(double ra0, double sdec0, double cdec0, 
                  double *ras, double *decs,
                  double *ls, double *ms, double *ns, int num_components);

/**
@brief Internally to WODEN , all cross-correlations are stored first,
then the autos after. All autos have zero length, so use this function to
set all autos uvw to zero.

@details `num_cross` is the number of cross-correlations, so fill any uvw after
this by adding `num_autos` worth of zeros

 @param[in] num_cross Number of cross-correlations
 @param[in] num_autos Number of auto-correlations
 @param[in,out] us u coords
 @param[in,out] vs v coords
 @param[in,out] ws w coords
*/
void set_auto_uvw_to_zero(int num_cross, int num_autos,
                          user_precision_t *us, user_precision_t *vs,
                          user_precision_t *ws);