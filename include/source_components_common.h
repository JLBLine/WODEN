#pragma once
#include "constants.h"
#include "woden_struct_defs.h"
#include "woden_precision_defs.h"
#include "primary_beam_cpu.h"
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

//Code defined in source_components_common.c------------------------------------

/**
 * @brief Performs necessary calculations that are common to all POINT,
 * GAUSSIAN, and SHAPELET component types: calculating the 
 * `lmn` coordinates; extrapolating flux densities to given frequencies;
 * calculating primary beam values.

 * @details This function either calls the CPU or GPU version of the
 * functions, depending on the value of `woden_settings->do_gpu`. The `mem_` prefix
 * indicates that the variable should be allocated on the GPU if `do_gpu` is set to 1.
 * Otherwise, it should be allocated on the CPU.
 *
 * If woden_settings->do_autos is 1, the function will calculate the auto-correlations
 *
 * @param woden_settings Pointer to the WODEN settings structure.
 * @param beam_settings Pointer to the beam settings structure.
 * @param cpu_freqs Pointer to the array of CPU frequencies.
 * @param mem_freqs Pointer to the array of memory frequencies.
 * @param chunked_source Pointer to the chunked source structure.
 * @param mem_chunked_source Pointer to the memory chunked source structure.
 * @param mem_component_beam_gains Pointer to the memory component beam gains structure.
 * @param comptype The type of component being processed.
 * @param mem_visibility_set Pointer to the memory visibility set structure.
 */
void source_component_common(woden_settings_t *woden_settings,
           beam_settings_t *beam_settings,
           double *cpu_freqs, double *mem_freqs,
           source_t *chunked_source, source_t *mem_chunked_source,
           beam_gains_t *mem_component_beam_gains,
           e_component_type comptype,
           visibility_set_t *mem_visibility_set);


/**
 * @brief Extrapolates Stokes IQUV parameters for a given `comptype` to a
 * requested set of frequencies `mem_extrap_freqs` for a given populated
 * `mem_chunked_source`. This function will extrapolate all available SED model
 * for the source.
 *
 * This function performs extrapolation of Stokes parameters for a given source
 * over a specified number of frequencies. The extrapolation can be performed
 * either on the CPU or GPU depending on the value of the `do_gpu` parameter.
 *
 * @param mem_chunked_source Pointer to the source data structure.
 * @param mem_extrap_freqs Pointer to an array of frequencies for extrapolation.
 * @param num_extrap_freqs Number of frequencies for extrapolation.
 * @param comptype The type of component for which the extrapolation is performed (POINT, GAUSSIAN, or SHAPELET)
 * @param do_gpu Flag indicating whether to perform the extrapolation on the GPU (1) or CPU (0).
 */
void extrapolate_Stokes(source_t *mem_chunked_source, double *mem_extrap_freqs,
                        int num_extrap_freqs, e_component_type comptype,
                        int do_gpu);

/**
@brief Fill the `ant1_to_baseline_map` and `ant2_to_baseline_map` arrays
with indexes corresponding to ant1 and ant2 for all unique baselines in an
array of `num_ants` antennas.

@details The `ant1_to_baseline_map` and `ant2_to_baseline_map` should
already have their memory allocated

@param[in] num_ants Number of antennas in the array
@param[in,out] *ant1_to_baseline_map Memory-allocated array of size `((num_ants - 1)*num_ants) / 2`
@param[in,out] *ant2_to_baseline_map Memory-allocated array of size `((num_ants - 1)*num_ants) / 2`

*/
void fill_ant_to_baseline_mapping_cpu(int num_ants, int *ant1_to_baseline_map,
                                      int *ant2_to_baseline_map);

//Put all the GPU code we are linking in here
//Only include this file in the final executable if we are compiling for GPU?

//External GPU code to link in--------------------------------------------------
#if !defined (__NVCC__) && !defined (__HIPCC__)

/**
* External GPU code linked in from source_components_gpu.cpp, see source_components_gpu.h for full details.
*/
extern void malloc_extrapolated_flux_arrays_gpu(components_t *d_components, int num_comps,
                                                int num_freqs);
/**
* External GPU code linked in from source_components_gpu.cpp, see source_components_gpu.h for full details.
*/
extern void extrap_power_laws_stokesI_gpu(components_t d_components,
                                   int n_powers, double *d_extrap_freqs,
                                   int num_extrap_freqs);
/**
* External GPU code linked in from source_components_gpu.cpp, see source_components_gpu.h for full details.
*/
extern void extrap_curved_power_laws_stokesI_gpu(components_t d_components,
                                   int n_curves, double *d_extrap_freqs,
                                   int num_extrap_freqs);
/**
* External GPU code linked in from source_components_gpu.cpp, see source_components_gpu.h for full details.
*/
extern void extrap_list_fluxes_stokesI_gpu(components_t d_components,
                                   int n_lists, double *d_extrap_freqs,
                                   int num_extrap_freqs);
/**
* External GPU code linked in from source_components_gpu.cpp, see source_components_gpu.h for full details.
*/
extern void extrap_power_laws_stokesV_gpu(components_t d_components,
                                              double *d_extrap_freqs,
                                              int num_extrap_freqs);
/**
* External GPU code linked in from source_components_gpu.cpp, see source_components_gpu.h for full details.
*/
extern void extrap_curved_power_laws_stokesV_gpu(components_t d_components,
                                                     double *d_extrap_freqs,
                                                     int num_extrap_freqs);
/**
* External GPU code linked in from source_components_gpu.cpp, see source_components_gpu.h for full details.
*/
extern void polarisation_fraction_stokesV_gpu(components_t d_components,
                                                     double *d_extrap_freqs,
                                                     int num_extrap_freqs);
/**
* External GPU code linked in from source_components_gpu.cpp, see source_components_gpu.h for full details.
*/
extern void extrap_list_fluxes_stokesV_gpu(components_t d_components,
                                               double *d_extrap_freqs,
                                               int num_extrap_freqs);
/**
* External GPU code linked in from source_components_gpu.cpp, see source_components_gpu.h for full details.
*/
extern void extrap_power_laws_linpol_gpu(components_t d_components,
                                               double *d_extrap_freqs,
                                               int num_extrap_freqs);
/**
* External GPU code linked in from source_components_gpu.cpp, see source_components_gpu.h for full details.
*/
extern void extrap_curved_power_laws_linpol_gpu(components_t d_components,
                                               double *d_extrap_freqs,
                                               int num_extrap_freqs);
/**
* External GPU code linked in from source_components_gpu.cpp, see source_components_gpu.h for full details.
*/
extern void polarisation_fraction_linpol_gpu(components_t d_components,
                                               double *d_extrap_freqs,
                                               int num_extrap_freqs);
/**
* External GPU code linked in from source_components_gpu.cpp, see source_components_gpu.h for full details.
*/
extern void extrap_list_fluxes_linpol_gpu(components_t d_components,
                                               double *d_extrap_freqs,
                                               int num_extrap_freqs);
/**
* External GPU code linked in from source_components_gpu.cpp, see source_components_gpu.h for full details.
*/
extern void extrap_p_list_fluxes_linpol_gpu(components_t d_components,
                                                double *d_extrap_freqs,
                                                int num_extrap_freqs);
/**
* External GPU code linked in from source_components_gpu.cpp, see source_components_gpu.h for full details.
*/
extern void apply_rotation_measure_gpu(components_t d_components,
                                                double *d_extrap_freqs,
                                                int num_extrap_freqs);
/**
* External GPU code linked in from source_components_gpu.cpp, see source_components_gpu.h for full details.
*/
extern void malloc_beam_gains_gpu(beam_gains_t *d_component_beam_gains,
                                     int beamtype, int num_gains);
/**
* External GPU code linked in from source_components_gpu.cpp, see source_components_gpu.h for full details.
*/
extern void calc_autos_gpu(components_t *d_components,
                               beam_settings_t *beam_settings,
                               beam_gains_t *d_component_beam_gains,
                               visibility_set_t *d_visibility_set,
                               woden_settings_t *woden_settings,
                               int num_components, int use_twobeams);

/**
 * External GPU code linked in from source_components_gpu.cpp, see source_components_gpu.h for full details.
 */
extern void copy_CPU_component_gains_to_GPU_beam_gains(components_t *components,
                                       beam_gains_t *d_beam_gains, int num_gains);

/**
 * External GPU code linked in from source_components_gpu.cpp, see source_components_gpu.h for full details.
 */
extern void copy_CPU_beam_gains_to_GPU_beam_gains(beam_gains_t *beam_gains,
                                     beam_gains_t *d_beam_gains, int num_gains);

/**
 * External GPU code linked in from fundamental_coods_gpu.cpp, see fundamental_coods_gpu.h for full details.
 */
extern void calc_lmn_for_components_gpu(components_t *d_components,
                                        int num_components,
                                        woden_settings_t *woden_settings);

/**
 * External GPU code linked in from primary_beam_gpu.cpp, see primary_beam_gpu.h for full details.
 */
extern void wrapper_calculate_gaussian_beam_gpu(int num_components,
                          user_precision_t cos_theta,
                          user_precision_t sin_theta, user_precision_t sin_2theta,
                          user_precision_t fwhm_lm,
                          woden_settings_t *woden_settings,
                          beam_settings_t *beam_settings,
                          components_t *components,
                          beam_gains_t *d_component_beam_gains,
                          double *d_freqs);

/**
 * External GPU code linked in from primary_beam_gpu.cpp, see primary_beam_gpu.h for full details.
 */
extern void wrapper_calculate_analytic_dipole_beam_gpu(int num_components,
                          components_t *components,
                          beam_gains_t *d_component_beam_gains,
                          double *d_freqs, woden_settings_t *woden_settings);

/**
 * External GPU code linked in from primary_beam_gpu.cpp, see primary_beam_gpu.h for full details.
 */
extern void wrapper_run_hyperbeam_gpu(int num_components,
                          beam_settings_t *beam_settings,
                          int num_beams, int parallactic,
                          double *reordered_azs, double *reordered_zas,
                          beam_gains_t *d_component_beam_gains,
                          double *d_freqs, woden_settings_t *woden_settings);

/**
 * External GPU code linked in from primary_beam_gpu.cpp, see primary_beam_gpu.h for full details.
 */
extern void wrapper_calculate_RTS_MWA_analytic_beam_gpu(int num_components,
                          components_t *components, int norm,
                          beam_gains_t *d_component_beam_gains,
                          double *d_freqs, woden_settings_t *woden_settings);



#endif