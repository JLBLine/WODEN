#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

#include "woden_precision_defs.h"
#include "woden_struct_defs.h"
#include "constants.h"

#include "source_components_common.h"
#include "hyperbeam_error.h"
#include "visibility_set.h"



/**
 * @brief Creates and initialises the calc_visi_inouts_t structure for CPU calculations.
 * This function gathers all inputs into one struct, and allocates memory
 * for intermediate outputs.
 * 
 * @details Specifically, we make the following pointers:
 * 
 * - calc_visi_inouts->allsteps_sha0s = visibility_set->allsteps_sha0s;
   - calc_visi_inouts->allsteps_cha0s = visibility_set->allsteps_cha0s;
   - calc_visi_inouts->allsteps_wavelengths = visibility_set->allsteps_wavelengths;
   - calc_visi_inouts->freqs = visibility_set->channel_frequencies;

 * and depending on the number of shapelets and whether we using a unique
 * beam for each station, we allocate memory for:
 * 
 * - calc_visi_inouts->u_metres
 * - calc_visi_inouts->v_metres
 * - calc_visi_inouts->w_metres
 * - calc_visi_inouts->us
 * - calc_visi_inouts->vs
 * - calc_visi_inouts->ws
 * - calc_visi_inouts->u_shapes
 * - calc_visi_inouts->v_shapes
 * - calc_visi_inouts->ant1_to_baseline_map
 * - calc_visi_inouts->ant2_to_baseline_map
 *
 * @param visibility_set Pointer to the visibility set structure.
 * @param sbf Pointer to the user precision structure.
 * @param woden_settings Pointer to the WODEN settings structure.
 * @param num_shapelets Number of shapelets to be used in the calculation.
 * @param use_twobeams Flag indicating whether to use two beams per baseline (1) or not (0) (assumes all primary beams are unique)
 * 
 * @return Pointer to the initialised calc_visi_inouts_t structure.
 */
calc_visi_inouts_t * create_calc_visi_inouts_cpu(array_layout_t *array_layout,
                visibility_set_t *visibility_set, user_precision_t *sbf, 
                woden_settings_t *woden_settings,
                int num_shapelets, int use_twobeams);

/**
 * @brief Sets the visibility set to zero.
 * 
 * @details This function sets the real and imaginary parts of all visibilities
 * in the visibility set to zero. Call before looping over sources to accumulate
 * visibilities.
 * 
 * @param visibility_set Pointer to the visibility set structure.
 * @param num_visis Number of visibilities to be set to zero.
 */
void set_visi_set_to_zero_cpu(visibility_set_t *visibility_set, int num_visis);

/**
 * @brief Frees the memory allocated for the calc_visi_inouts_t structure.
 * 
 * @details Frees everything that was allocated in `create_calc_visi_inouts_cpu`.
 * 
 * @param calc_visi_inouts Pointer to the calc_visi_inouts_t structure.
 * @param num_shapelets Number of shapelets in the original source_catalogue.
 * @param use_twobeams Flag indicating whether to use two beams per baseline (1) or not (0) (assumes all primary beams are unique)
 */
void free_calc_visi_inouts_cpu(calc_visi_inouts_t *calc_visi_inouts,
                               int num_shapelets, int use_twobeams);

/**
 * @brief Frees any arrays that were allocated when calculating `lmn` coords
 * and extrapolating Stokes fluxes; doesn't free the whole `components_t` struct
 * as that should have been created and freed elsewhere (Python via ctypes
 * in `run_woden.py`)
 * 
 * @param components Pointer to the components structure.
 * 
*/
void free_components_cpu(components_t *components);

/**
 * @brief Frees the memory allocated for the beam gains structure.
 * 
 * @details Depending on the beamtype, either frees beams gains or beamg gains
 * and leakages.
 * 
 * @param beam_gains Pointer to the beam gains structure.
 * @param beamtype The type of beam to be freed (see e_beamtype in woden_struct_defs.h).
 */
void free_beam_gains_cpu(beam_gains_t *beam_gains, e_beamtype beamtype);