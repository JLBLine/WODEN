/*! \file
  Methods to read in simulation parameters from a .json file and prepare
  settings into a `woden_settings_t` struct

  @author J.L.B. Line
*/
#pragma once
#include <math.h>
#include <stdint.h>
#include "woden_struct_defs.h"

/**
@brief Takes a path to .json WODEN parameter file, and populates a
`woden_settings_t` struct with the contents of `filename`.

@details For what can be included in the .json file, see the documentation
for
@todo Work out how to link the print_help function here

@param[in] *filename Path to a WODEN *.json settings file
@param[in,out] *woden_settings A pointer to a populated `woden_settings_t` struct
 @return Integer where 0 if read was successful, 1 if failed
 */
int read_json_settings(const char *filename, woden_settings_t *woden_settings);

/**
@brief Creates and returns an array of LST values for all time steps in
the simulation based on the settings in `woden_settings`

@details Also sets `woden_settings->num_visis` to the overall number of
visiblities (baselines*freqs*time_steps) and takes the sine and cosine of the
declination of the phase centre and stores in `woden_settings->sdec0`,
`woden_settings->cdec0`

@param[in] *filename Path to a WODEN *.json settings file
@param[in,out] *woden_settings A pointer to a populated `woden_settings_t` struct
@return `lsts` Array of all LSTs in the simulation (radians)
*/
float * setup_lsts_and_phase_centre(woden_settings_t *woden_settings);
