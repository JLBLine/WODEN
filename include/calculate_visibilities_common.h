#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

#include "woden_precision_defs.h"
#include "woden_struct_defs.h"
#include "constants.h"

#include "calculate_visibilities_cpu.h"
#include "source_components_common.h"
#include "source_components_cpu.h"
#include "hyperbeam_error.h"
#include "visibility_set.h"
#include "logger.h"


//Function defined in calculate_visibilities_common.c---------------------------

/**
 * @brief Calculate the visibilities for a given component type, on either the
 * CPU or GPU, for a single chunk of the sky model
 *
 * This function computes primary beam reponses and calculates the visibilities
 * for a specified component type (POINT, GAUSSIAN, or SHAPELET) using the
 * provided input and output structures, frequency channels, and settings.
 * Any input with `mem_` prefix needs to be allocated on the GPU if `do_gpu`
 * is set to 1. Otherwise it should be allocated on the CPU.
 * 
 * Because EveryBeam only works on the CPU, even when running in GPU mode,
 * need access to the CPU version of the `source`, hence we have to input
 * `source` and `mem_chunked_source` separately. Stops multiple GPU copying
 * doing the copy in `mem_chunked_source` outside of this function.
 *
 * @param comptype The type of component for which to calculate visibilities.
 * @param mem_calc_visi_inouts Pointer to the structure containing input and output data for visibility calculation.
 * @param cpu_channel_freqs Pointer to the array of frequency channels.
 * @param woden_settings Pointer to the WODEN settings structure.
 * @param beam_settings Pointer to the beam settings structure.
 * @param source Pointer to the source structure.
 * @param mem_chunked_source Pointer to the chunked source structure.
 * @param mem_visibility_set Pointer to the visibility set structure.
 * @param num_beams The number of beams to be used in the calculation.
 * @param use_twobeams Flag indicating whether to use two beams (1) or not (0). This means use two different beams per baseline, so assumes a unique beam for each station.
 * @param do_gpu Flag indicating whether to perform the calculation on the GPU (1) or CPU (0).
 */
void calculate_component_visis(e_component_type comptype,
  calc_visi_inouts_t *mem_calc_visi_inouts,
  double *cpu_channel_freqs,
  woden_settings_t *woden_settings,
  beam_settings_t *beam_settings,
  source_t *source, source_t *mem_chunked_source,
  visibility_set_t *mem_visibility_set,
  int num_beams, int use_twobeams,
  int do_gpu);

/**
 * @brief Calculates the visibilities for a given array layout and source catalogue.
 *
 * This function computes the visibilities based on the provided array layout,
 * cropped sky models, beam settings, and WODEN settings. The results are stored
 * in the visibility set.
 * 
 * `woden_settings` contains settings indicating what kind of simulation
 * is to be run. `cropped_sky_models` contains the sky model, chunked into
 * pieces that should fit on the CPU/GPU. `array_layout` hold the array layout,
 *  and `sbf` holds the shapelet basis functions. `beam_settings` holds the primary
 * some primary beam settings, and initialised `hyperbeam` objects. `visibility_set`
 * should be an array of `visibility_set_t` structs, one for each coarse band. 
 * All should be fully initialised before calling this function. 
 * 
  NOTE The following attributes *must* be set for `calculate_visibilities` to work:

   - array_layout->X_diff_metres
   - array_layout->Y_diff_metres
   - array_layout->Z_diff_metres
   - cropped_sky_models->num_sources
   - cropped_sky_models->sources
   - beam_settings->beamtype
   - woden_settings->num_baselines
   - woden_settings->num_time_steps
   - woden_settings->num_visis
   - woden_settings->num_freqs
   - woden_settings->beamtype
   - woden_settings->ra0
   - woden_settings->sdec0
   - woden_settings->cdec0
   - visibility_set->allsteps_cha0s
   - visibility_set->allsteps_sha0s
   - visibility_set->allsteps_lsts
   - visibility_set->allsteps_wavelengths
   - visibility_set->channel_frequencies

  The following must be malloc-ed for `calculate_visibilities` to work:

   - visibility_set->us_metres
   - visibility_set->vs_metres
   - visibility_set->ws_metres
   - visibility_set->sum_visi_XX_real
   - visibility_set->sum_visi_XX_imag
   - visibility_set->sum_visi_XY_real
   - visibility_set->sum_visi_XY_imag
   - visibility_set->sum_visi_YX_real
   - visibility_set->sum_visi_YX_imag
   - visibility_set->sum_visi_YY_real
   - visibility_set->sum_visi_YY_imag

  for all `visibility_set` in `visibility_sets`. Beyond that, depending on your
  settings, other attributes may be required.
 *
 * @param array_layout Pointer to the array layout structure.
 * @param cropped_sky_models Pointer to the cropped sky models structure.
 * @param beam_settings Pointer to the beam settings structure.
 * @param woden_settings Pointer to the WODEN settings structure.
 * @param visibility_set Pointer to the visibility set structure where results will be stored.
 * @param sbf Pointer to the shapelet basis functions array
 */
void calculate_visibilities(array_layout_t *array_layout,
  source_catalogue_t *cropped_sky_models, beam_settings_t *beam_settings,
  woden_settings_t *woden_settings, visibility_set_t *visibility_set,
  user_precision_t *sbf);



//External GPU code we're linking in--------------------------------------------

/**
 * External GPU code linked in from calculate_visibilities_gpu.c, see calculate_visibilities_gpu.h for full details.
 */
extern calc_visi_inouts_t * create_calc_visi_inouts_gpu(array_layout_t *array_layout,
        visibility_set_t *visibility_set, visibility_set_t *d_visibility_set,
        user_precision_t *sbf, woden_settings_t *woden_settings,
        int num_shapelets, int use_twobeams);

/**
 * External GPU code linked in from calculate_visibilities_gpu.c, see calculate_visibilities_gpu.h for full details.
 */
extern void set_visibilities_to_zero_gpu(visibility_set_t *d_visibility_set,
                              visibility_set_t *chunk_visibility_set, int num_visis);

/**
 * External GPU code linked in from calculate_visibilities_gpu.c, see calculate_visibilities_gpu.h for full details.
 */
extern void copy_gpu_visi_set_to_host(visibility_set_t *d_visibility_set,
  calc_visi_inouts_t *d_calc_visi_inouts,
  visibility_set_t *chunk_visibility_set,
  int num_visis);

/**
* External GPU code linked in from calculate_visibilities_gpu.c, see calculate_visibilities_gpu.h for full details.
*/
extern void free_calc_visi_inouts_gpu(calc_visi_inouts_t *d_calc_visi_inouts,
      visibility_set_t *d_visibility_set,
      int num_shapelets, int use_twobeams);

// /**
// * External GPU code linked in from source_components_gpu.cpp, see source_components_gpu.h for full details.
// */
// extern void fill_ant_to_baseline_mapping_gpu(int num_ants, int *d_ant1_to_baseline_map,
//         int *d_ant2_to_baseline_map);

/**
* External GPU code linked in from source_components_gpu.cpp, see source_components_gpu.h for full details.
*/
extern source_t * copy_chunked_source_to_GPU(source_t *chunked_source, woden_settings_t *woden_settings);

/**
* External GPU code linked in from source_components_gpu.cpp, see source_components_gpu.h for full details.
*/
extern void free_extrapolated_flux_arrays_gpu(components_t *d_components);

/**
* External GPU code linked in from source_components_gpu.cpp, see source_components_gpu.h for full details.
*/
extern void free_components_gpu(source_t *d_chunked_source, e_component_type comptype);

/**
* External GPU code linked in from source_components_gpu.cpp, see source_components_gpu.h for full details.
*/
extern void free_beam_gains_gpu(beam_gains_t *d_beam_gains, e_beamtype beamtype);

/**
* External GPU code linked in from source_components_gpu.cpp, see source_components_gpu.h for full details.
*/
extern void calc_visi_point_or_gauss_gpu(components_t d_components,
                                        beam_gains_t d_component_beam_gains,
                                        calc_visi_inouts_t *d_calc_visi_inouts,
                                        visibility_set_t *d_visibility_set,
                                        int num_components, e_beamtype beamtype,
                                        e_component_type comptype,
                                        woden_settings_t *woden_settings);
/**
* External GPU code linked in from source_components_gpu.cpp, see source_components_gpu.h for full details.
*/
extern void calc_visi_shapelets_gpu(components_t d_components,
                                        beam_gains_t d_component_beam_gains,
                                        calc_visi_inouts_t *d_calc_visi_inouts,
                                        visibility_set_t *d_visibility_set,
                                        int num_shapes, int num_shape_coeffs,
                                        e_beamtype beamtype,
                                        woden_settings_t *woden_settings);


/**
* External GPU code linked in from fundamental_coords_gpu.cpp, see fundamental_coords_gpu.h for full details.
*/
extern void calc_uvw_gpu(double *d_X_diff, double *d_Y_diff, double *d_Z_diff,
                             user_precision_t *d_u_metres,
                             user_precision_t *d_v_metres, user_precision_t *d_w_metres,
                             user_precision_t *d_us, user_precision_t *d_vs,
                             user_precision_t *d_ws, user_precision_t *d_allsteps_wavelengths,
                             double *d_allsteps_cha0s, double *d_allsteps_sha0s,
                             woden_settings_t *woden_settings);

/**
* External GPU code linked in from fundamental_coords_gpu.cpp, see fundamental_coords_gpu.h for full details.
*/
extern void calc_uv_shapelet_gpu(user_precision_t *d_u_shapes, user_precision_t *d_v_shapes,
                                    int num_shapes,
                                    double *d_X_diff, double *d_Y_diff, double *d_Z_diff,
                                    double *d_ras, double *d_decs,
                                    woden_settings_t *woden_settings);
                                    
/**
* External GPU code linked in from fundamental_coords_gpu.cpp, see fundamental_coords_gpu.h for full details.
*/
extern void set_auto_uvw_to_zero_gpu(int num_cross, int num_autos,
                              user_precision_t *d_u, user_precision_t *d_v,
                              user_precision_t *d_w);