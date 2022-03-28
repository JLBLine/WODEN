/*! \file
Functions to calculate the parallactic angle, and fill in primary beam settings
requested for the simulation, based on settings in `woden_settings` and the sky
model `cropped_src`
*/


#include "constants.h"
#include "woden_struct_defs.h"

/**
@brief Given the settings specified in `woden_settings`, populate a
`beam_settings_t` and `source_t` struct with appropriate attributes to be
used in primary beam modelling.

@details If using a `GAUSS_BEAM` or `MWA_ANALY` primary beam, calculate the hour angle and
declination of the pointing centre for all times steps, and for all components in
sky model `cropped_src`.

@param[in,out] *cropped_src A populated `source_t` sky model
@param[in] *woden_settings A populated `woden_settings_t` struct
@param[in] *lsts All local sidereal times in the simulation

@return `beam_settings` - a populated `beam_settings_t` struct containing
attributes necessary to simulate the requested beam response

@todo For all sky simulations, there can be millions of components, so
calculating parallactic angle can be expensive. Consider multi-threading this
somehow

*/
beam_settings_t * fill_primary_beam_settings(woden_settings_t *woden_settings,
                            source_t *cropped_src, double *lsts);
