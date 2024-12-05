#pragma once
#include "constants.h"
#include "woden_struct_defs.h"
#include "woden_precision_defs.h"

#if defined (__NVCC__) || defined (__HIPCC__)
#include "gpu_macros.h"
#endif

/*!
A struct to contain primary beam values for a give COMPONENT. `d_gxs,d_Dxs,d_Dys,d_gys`
should be used when all antennas have the same primary beam, and `d_gxs,d_Dxs,d_Dys,d_gys` used when all primary beams are different.
*/
typedef struct _d_beam_gains_t {

  int *d_ant1_to_baseline_map; /*!< The index of antenna 1 in all unique pairs of
antennas. Used to map iBaseline to the correct antenna 1 */
  int *d_ant2_to_baseline_map; /*!< The index of antenna 2 in all unique pairs of
antennas. Used to map iBaseline to the correct antenna 2 */
  int use_twobeams; /*!< The beam gains were made with unique primary beams so
  should use two antenna patterns per visibility */

  user_precision_complex_t *gxs; /*!< Host copy of North-South Beam gain values
  for all beams, directions, frequencies, and times for these COMPONENTS*/
  user_precision_complex_t *Dxs; /*!< Host copy of North-South Beam leakage values
  for all beams, directions, frequencies, and times for these COMPONENTS*/
  user_precision_complex_t *Dys; /*!< Host copy of East-West Beam leakage values
  for all beams, directions, frequencies, and times for these COMPONENTS*/
  user_precision_complex_t *gys; /*!< Host copy of East-West Beam gain values
  for all beams, directions, frequencies, and times for these COMPONENTS*/

  user_precision_complex_t *d_gxs; /*!< Device copy of North-South Beam gain values
  for all beams, directions, frequencies, and times for these COMPONENTS*/
  user_precision_complex_t *d_Dxs; /*!< Device copy of North-South Beam leakage values
  for all beams, directions, frequencies, and times for these COMPONENTS*/
  user_precision_complex_t *d_Dys; /*!< Device copy of East-West Beam leakage values
  for all beams, directions, frequencies, and times for these COMPONENTS*/
  user_precision_complex_t *d_gys; /*!< Device copy of East-West Beam gain values
  for all beams, directions, frequencies, and times for these COMPONENTS*/

  // //OK, so when this is compiled for CPU, the compiler doesn't have access
  // //to GPU libraries. So wrap the GPU specific stuff in an ifdef
  // //that only triggers when we're compiling for GPU
  // #if defined(__NVCC__) || defined(__HIPCC__)
  // gpuUserComplex *d_gxs; /*!< Device copy of North-South Beam gain values
  // for all beams, directions, frequencies, and times for these COMPONENTS*/
  // gpuUserComplex *d_Dxs; /*!< Device copy of North-South Beam leakage values
  // for all beams, directions, frequencies, and times for these COMPONENTS*/
  // gpuUserComplex *d_Dys; /*!< Device copy of East-West Beam leakage values
  // for all beams, directions, frequencies, and times for these COMPONENTS*/
  // gpuUserComplex *d_gys; /*!< Device copy of East-West Beam gain values
  // for all beams, directions, frequencies, and times for these COMPONENTS*/
  // #endif


} d_beam_gains_t;

void source_component_common(woden_settings_t *woden_settings,
           beam_settings_t *beam_settings, double *d_freqs,
           source_t *chunked_source, source_t *d_chunked_source,
           d_beam_gains_t *d_component_beam_gains,
           e_component_type comptype,
           visibility_set_t *d_visibility_set);

void extrapolate_Stokes(source_t *d_chunked_source, double *d_extrap_freqs,
                        int num_extrap_freqs, e_component_type comptype,
                        int do_gpu);

//Put all the GPU code we are linking in here
//Only include this file in the final executable if we are compiling for GPU?

#if !defined (__NVCC__) && !defined (__HIPCC__)

extern void malloc_extrapolated_flux_arrays_gpu(components_t *d_components, int num_comps,
                                        int num_freqs);

extern void extrap_power_laws_stokesI_gpu(components_t d_components,
                                   int n_powers, double *d_extrap_freqs,
                                   int num_extrap_freqs);

extern void extrap_curved_power_laws_stokesI_gpu(components_t d_components,
                                   int n_curves, double *d_extrap_freqs,
                                   int num_extrap_freqs);

extern void extrap_list_fluxes_stokesI_gpu(components_t d_components,
                                   int n_lists, double *d_extrap_freqs,
                                   int num_extrap_freqs);

extern void extrap_power_laws_stokesV_gpu(components_t d_components,
                                              double *d_extrap_freqs,
                                              int num_extrap_freqs);

extern void extrap_curved_power_laws_stokesV_gpu(components_t d_components,
                                                     double *d_extrap_freqs,
                                                     int num_extrap_freqs);

extern void polarisation_fraction_stokesV_gpu(components_t d_components,
                                                     double *d_extrap_freqs,
                                                     int num_extrap_freqs);

extern void extrap_list_fluxes_stokesV_gpu(components_t d_components,
                                               double *d_extrap_freqs,
                                               int num_extrap_freqs);

extern void extrap_power_laws_linpol_gpu(components_t d_components,
                                               double *d_extrap_freqs,
                                               int num_extrap_freqs);

extern void extrap_curved_power_laws_linpol_gpu(components_t d_components,
                                               double *d_extrap_freqs,
                                               int num_extrap_freqs);

extern void polarisation_fraction_linpol_gpu(components_t d_components,
                                               double *d_extrap_freqs,
                                               int num_extrap_freqs);

extern void extrap_list_fluxes_linpol_gpu(components_t d_components,
                                               double *d_extrap_freqs,
                                               int num_extrap_freqs);

extern void extrap_p_list_fluxes_linpol_gpu(components_t d_components,
                                                double *d_extrap_freqs,
                                                int num_extrap_freqs);

extern void apply_rotation_measure_gpu(components_t d_components,
                                                double *d_extrap_freqs,
                                                int num_extrap_freqs);

extern void malloc_beam_gains_gpu(d_beam_gains_t *d_component_beam_gains,
                                     int beamtype, int num_gains);

extern void calc_lmn_for_components_gpu(components_t *d_components,
                                        int num_components,
                                        woden_settings_t *woden_settings);

extern void wrapper_calculate_gaussian_beam_gpu(int num_components,
                          user_precision_t cos_theta,
                          user_precision_t sin_theta, user_precision_t sin_2theta,
                          user_precision_t fwhm_lm,
                          woden_settings_t *woden_settings,
                          beam_settings_t *beam_settings,
                          components_t *components,
                          d_beam_gains_t *d_component_beam_gains,
                          double *d_freqs);

extern void wrapper_calculate_analytic_dipole_beam_gpu(int num_components,
                          components_t *components,
                          d_beam_gains_t *d_component_beam_gains,
                          double *d_freqs, woden_settings_t *woden_settings);

extern void wrapper_run_hyperbeam_gpu(int num_components,
                          components_t *components, beam_settings_t *beam_settings,
                          int num_beams, int parallactic,
                          double *reordered_azs, double *reordered_zas,
                          d_beam_gains_t *d_component_beam_gains,
                          double *d_freqs, woden_settings_t *woden_settings);

extern void wrapper_calculate_RTS_MWA_analytic_beam_gpu(int num_components,
                          components_t *components, int norm,
                          d_beam_gains_t *d_component_beam_gains,
                          double *d_freqs, woden_settings_t *woden_settings);

extern void calc_autos_gpu(components_t *d_components,
                               beam_settings_t *beam_settings,
                               d_beam_gains_t *d_component_beam_gains,
                               visibility_set_t *d_visibility_set,
                               woden_settings_t *woden_settings,
                               int num_components, int use_twobeams);

#endif