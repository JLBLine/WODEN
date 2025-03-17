/*! \file
Functions to calculate the parallactic angle, and fill in primary beam settings
requested for the simulation, based on settings in `woden_settings` and the sky
model `cropped_src`
*/


#include "constants.h"
#include "woden_struct_defs.h"
#include "stdlib.h"
#include "logger.h"

/**
@brief Given the settings specified in `woden_settings`, populate a
`beam_settings_t` struct with appropriate attributes to be
used in primary beam modelling.

@details Not much doing here, just sets up the `beam_settings_t` struct with some useful values.

@param[in] *woden_settings A populated `woden_settings_t` struct
@param[in] *lsts All local sidereal times in the simulation

@return `beam_settings` - a populated `beam_settings_t` struct containing
attributes necessary to simulate the requested beam response

*/
beam_settings_t * fill_primary_beam_settings(woden_settings_t *woden_settings,
                                             double *lsts);
