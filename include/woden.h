#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "woden_precision_defs.h"
#include "constants.h"
#include "woden_struct_defs.h"
#include "beam_settings.h"
#include "visibility_set.h"
#include "hyperbeam_error.h"
#include "calculate_visibilities_common.h"
#include "logger.h"


/**
 * @brief Runs the WODEN simulation with the provided settings and data.
 * 
 * @details `woden_settings` contains settings indicating what kind of simulation
 * is to be run. `cropped_sky_models` contains the sky model, chunked into
 * pieces that should fit on the CPU/GPU. `array_layout` hold the array layout,
 *  and `sbf` holds the shapelet basis functions. All should be fully
 * initialised before calling this function. The intention is
 * to call this `run_woden` function from `run_woden.py`, so see that for
 * how to set up the input parameters.
 * 
 * `visibility_sets` should be an array of `visibility_set_t` structs, one for
 * each coarse band. Within each `visibility_set`, the following arrays should be allocated:
 *  - visibility_set->us_metres
 *  - visibility_set->vs_metres
 *  - visibility_set->ws_metres
 *  - visibility_set->sum_visi_XX_real
 *  - visibility_set->sum_visi_XX_imag
 *  - visibility_set->sum_visi_XY_real
 *  - visibility_set->sum_visi_XY_imag
 *  - visibility_set->sum_visi_YX_real
 *  - visibility_set->sum_visi_YX_imag
 *  - visibility_set->sum_visi_YY_real
 *  - visibility_set->sum_visi_YY_imag
 * 
 * as these hold the outputs of the simulation.
 * 
 * All other members of `visibility_set` will be populated/freed by this function.
 *
 * @param woden_settings Pointer to the settings structure for WODEN.
 * @param visibility_sets Pointer to the visibility sets structure.
 * @param cropped_sky_models Pointer to the source catalogue structure containing cropped sky models.
 * @param array_layout Pointer to the array layout structure.
 * @param sbf Pointer to the shapelet basis functions array
 *
 * @return An integer indicating the success or failure of the simulation.
 */
int run_woden(woden_settings_t *woden_settings, visibility_set_t *visibility_sets,
             source_catalogue_t *cropped_sky_models, array_layout_t * array_layout,
             user_precision_t *sbf);