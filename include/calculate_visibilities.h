/*! \file
  Primary function to take the simulation settings and sky models generated
  on the host, transfer everything over to the device, run the simulation
  through the CUDA code, and grab everything back off the device and onto the
  host.
  @author J.L.B. Line
*/

#pragma once
#include <math.h>
#include "woden_struct_defs.h"

/**
@brief Given the telescope, settings, and sky models detailed in `array_layout`,
 `cropped_sky_models`, and `woden_settings`, run the GPU simulation, and store
 the output visibilities in `visibility_set`

@details Uses the telescope model described in `array_layout`, an
array of sky models and beam settings described by `cropped_sky_models`,
and the simulation settings described in `woden_settings` to run the GPU
simulation. `sbf` is the shapelet basis function, which will be copied in device
memory if SHAPELET components are present in the sky models. `num_chunks`
corresponds to the length of the `cropped_sky_models` array. The function will
run each sky model in `cropped_sky_models` as a separate simulation, and sum
the resultant visibilities into `visibility_set`.
@param[in] *array_layout Pointer to an `array_layout_t` struct
@param[in] *cropped_sky_models Pointer to a populated `source_catalogue_t`
struct
@param[in] *woden_settings Pointer to a populated `woden_settings_t` struct
@param[in,out] *visibility_set Pointer to an initialised visibility_set_t struct
@param[in] *sbf An array of gridded shapelet basis function values as created by
`create_sbf`

      NOTE

      The following attributes *must* be set for `calculate_visibilities` to work:

      array_layout->X_diff_metres
      array_layout->Y_diff_metres
      array_layout->Z_diff_metres
      cropped_sky_models->num_sources
      cropped_sky_models->catsources
      cropped_sky_models->beam_settings
      woden_settings->num_baselines
      woden_settings->num_time_steps
      woden_settings->num_visis
      woden_settings->num_freqs
      woden_settings->beamtype
      woden_settings->ra0
      woden_settings->sdec0
      woden_settings->cdec0
      visibility_set->allsteps_cha0s
      visibility_set->allsteps_sha0s
      visibility_set->allsteps_lsts
      visibility_set->allsteps_wavelengths
      visibility_set->channel_frequencies

      The following must be malloc-ed for `calculate_visibilities` to work:

      visibility_set->sum_visi_real
      visibility_set->sum_visi_imag
      visibility_set->us_metres
      visibility_set->vs_metres
      visibility_set->ws_metres
      visibility_set->sum_visi_XX_real
      visibility_set->sum_visi_XX_imag
      visibility_set->sum_visi_XY_real
      visibility_set->sum_visi_XY_imag
      visibility_set->sum_visi_YX_real
      visibility_set->sum_visi_YX_imag
      visibility_set->sum_visi_YY_real
      visibility_set->sum_visi_YY_imag

      Beyond that, depending on your primary beam settings,
      other attributes may be required.
*/
extern "C" void calculate_visibilities(array_layout_t * array_layout,
  source_catalogue_t *cropped_sky_models, woden_settings_t *woden_settings,
  visibility_set_t *visibility_set, float *sbf);
