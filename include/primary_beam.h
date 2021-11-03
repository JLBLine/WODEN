/*! \file
Functions to calculate the parallactic angle, and fill in primary beam settings
requested for the simulation, based on settings in `woden_settings` and the sky
model `cropped_src`
*/


#include "constants.h"
#include "woden_struct_defs.h"

/**
@brief Given the component positions in the sky model `cropped_src`,
calculate the parallactic angle.

@details The MWA FEE beam is stored in instrumental polarisations (aligned
with azimuth and elevation), which must be rotated to align with celestial
based XX/YY. This is achieved by rotating by parallactic angle. Parallactic
angle is calculated using `<erfa> eraHd2pa`. The `sine` and `cosine` of the
parallactic angles are stored separately for each COMPONENT type in

      cropped_src->sin_point_para_angs
      cropped_src->cos_point_para_angs
      cropped_src->sin_gauss_para_angs
      cropped_src->cos_gauss_para_angs
      cropped_src->sin_shape_para_angs
      cropped_src->cos_shape_para_angs

as they are used later by `kern_rotate_FEE_beam` to perform the rotation on the
device.

@param[in,out] *cropped_src A populated `catsource_t` sky model
@param[in] *lsts All local sidereal times in the simulation
@param[in] latitude Latitude of the array (radians)
@param[in] num_time_steps Number of time steps in the simulation
*/
void calc_para_angle(catsource_t *cropped_src, user_precision_t *lsts,
                     double latitude, int num_time_steps);


/**
@brief Given the settings specified in `woden_settings`, populate a
`beam_settings_t` and `catsource_t` struct with appropriate attributes to be
used in primary beam modelling.

@details If using a `GAUSS_BEAM` primary beam, calculate the hour angle and
declination of the pointing centre for all times steps, and for all components in
sky model `cropped_src`. If using `FEE_BEAM`, calculate the parallactic angle
for all times steps, for all components.

@param[in,out] *cropped_src A populated `catsource_t` sky model
@param[in] *woden_settings A populated `woden_settings_t` struct
@param[in] *lsts All local sidereal times in the simulation

@return `beam_settings` - a populated `beam_settings_t` struct containing
attributes necessary to simulate the requested beam response

@todo For all sky simulations, there can be millions of components, so
calculating parallactic angle can be expensive. Consider multi-threading this
somehow

*/
beam_settings_t * fill_primary_beam_settings(woden_settings_t *woden_settings,
                            catsource_t *cropped_src, user_precision_t *lsts);
