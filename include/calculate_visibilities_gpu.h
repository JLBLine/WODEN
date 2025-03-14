/*! \file
  Primary function to take the simulation settings and sky models generated
  on the host, transfer everything over to the device, run the simulation
  through the GPU code, and grab everything back off the device and onto the
  host.
  @author J.L.B. Line
*/

#pragma once
#include <math.h>
#include "woden_struct_defs.h"

/**
 * @brief Creates and initialises the `calc_visi_inouts_t` structure for GPU calculations.
 * This function gathers all inputs into one struct, and allocates memory
 * for intermediate outputs. Makes it easy to pass all the inputs and outputs
 * to the GPU functions.
 * 
 * @details Specifically, we copy the following into GPU mempory:
 * 
 * - `d_calc_visi_inouts->X_diff` from `array_layout->X_diff_metres`
 * - `d_calc_visi_inouts->Y_diff` from `array_layout->Y_diff_metres`
 * - `d_calc_visi_inouts->Z_diff` from `array_layout->Z_diff_metres`
 * - `d_calc_visi_inouts->allsteps_sha0s` from `visibility_set->allsteps_sha0s`
 * - `d_calc_visi_inouts->allsteps_cha0s` from `visibility_set->allsteps_cha0s`
 * - `d_calc_visi_inouts->allsteps_wavelengths` from `visibility_set->allsteps_wavelengths`
 * - `d_calc_visi_inouts->freqs` from `visibility_set->channel_frequencies`
 * - `d_calc_visi_inouts->sbf` from `sbf` (if we have shapelets in our sky model)
 *
 * and depending on the number of shapelets and whether we using a unique
 * beam for each station, we allocate GPU memory for:
 * 
 * - `d_calc_visi_inouts->u_metres`
 * - `d_calc_visi_inouts->v_metres`
 * - `d_calc_visi_inouts->w_metres`
 * - `d_calc_visi_inouts->us`
 * - `d_calc_visi_inouts->vs`
 * - `d_calc_visi_inouts->ws`
 * - `d_calc_visi_inouts->d_u_shapes`
 * - `d_calc_visi_inouts->d_v_shapes`
 * - `d_calc_visi_inouts->d_ant1_to_baseline_map`
 * - `d_calc_visi_inouts->d_ant2_to_baseline_map`
 * - `d_visibility_set->sum_visi_XX_real`
 * - `d_visibility_set->sum_visi_XX_imag`
 * - `d_visibility_set->sum_visi_XY_real`
 * - `d_visibility_set->sum_visi_XY_imag`
 * - `d_visibility_set->sum_visi_YX_real`
 * - `d_visibility_set->sum_visi_YX_imag`
 * - `d_visibility_set->sum_visi_YY_real`
 * - `d_visibility_set->sum_visi_YY_imag`
 *
 * @param array_layout Pointer to the array layout structure.
 * @param visibility_set Pointer to the visibility set structure.
 * @param d_visibility_set Pointer to the GPU visibility set structure (shouldn't have any attribute initialised)
 * @param sbf Pointer to shapelet basis function array
 * @param woden_settings Pointer to the WODEN settings structure.
 * @param num_shapelets Number of shapelets in the source_catalogue
 * @param use_twobeams Flag indicating whether to use two beams per baseline (1) or not (0) (assumes all primary beams are unique)
 * 
 * @return Pointer to the initialised calc_visi_inouts_t structure.
 */
extern "C" calc_visi_inouts_t * create_calc_visi_inouts_gpu(array_layout_t *array_layout,
  visibility_set_t *visibility_set, visibility_set_t *d_visibility_set,
  user_precision_t *sbf, woden_settings_t *woden_settings,
  int num_shapelets, int use_twobeams);

/**
 * @brief Sets the visibility set to zero.
 * 
 * @details This function sets the real and imaginary parts of all visibilities
 * in the GPU visibility set to zero, by copying the zeroed arrays in the
 * `chunk_visibility_set` to the `d_visibility_set`.
 * 
 * @param d_visibility_set Pointer to the GPU visibility set structure.
 * @param chunk_visibility_set Pointer to the CPU visibility set structure. This 
 * should have it's visibilities set to zero
 * @param num_visis Number of visibilities to be set to zero.
 */
extern "C" void set_visibilities_to_zero_gpu(visibility_set_t *d_visibility_set,
    visibility_set_t *chunk_visibility_set, int num_visis);

/**
 * @brief Copies the GPU visibilities and `uvw` coords set back to the CPU
 * visibility set.
 * 
 * @details This function copies the GPU visibility arrays out of `d_visibility_set`,
 * and the `uvw` coords out of `d_calc_visi_inouts` into `chunk_visibility_set`.
 * 
 * @param d_visibility_set Pointer to the GPU visibility set structure.
 * @param d_calc_visi_inouts Pointer to the GPU calc_visi_inouts structure.
 * @param chunk_visibility_set Pointer to the CPU visibility set structure.
 * @param num_visis Number of visibilities to be copied.
 */
extern "C" void copy_gpu_visi_set_to_host(visibility_set_t *d_visibility_set,
      calc_visi_inouts_t *d_calc_visi_inouts,
      visibility_set_t *chunk_visibility_set,
      int num_visis);

/**
 * @brief Frees the GPU memory allocated for the `calc_visi_inouts_t` structure.
 * 
 * @details Frees everything that was allocated in `create_calc_visi_inouts_gpu`.
 * 
 * @param d_calc_visi_inouts Pointer to the GPU `calc_visi_inouts_t` structure.
 * @param d_visibility_set Pointer to the GPU `visibility_set_t` structure.
 * @param num_shapelets Number of shapelets in the original source_catalogue.
 * @param use_twobeams Flag indicating whether to use two beams per baseline (1) or not (0) (assumes all primary beams are unique)
 */
extern "C"  void free_calc_visi_inouts_gpu(calc_visi_inouts_t *d_calc_visi_inouts,
        visibility_set_t *d_visibility_set,
        int num_shapelets, int use_twobeams);