#pragma once
#include "constants.h"
#include "woden_struct_defs.h"
#include "woden_precision_defs.h"
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

void source_component_common(woden_settings_t *woden_settings,
           beam_settings_t *beam_settings, double *d_freqs,
           source_t *chunked_source, source_t *d_chunked_source,
           beam_gains_t *d_component_beam_gains,
           e_component_type comptype,
           visibility_set_t *d_visibility_set);

void extrapolate_Stokes(source_t *mem_chunked_source, double *mem_extrap_freqs,
                        int num_extrap_freqs, e_component_type comptype,
                        int do_gpu);

void fill_ant_to_baseline_mapping_cpu(int num_ants, int *ant1_to_baseline_map,
                                      int *ant2_to_baseline_map);

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

extern void malloc_beam_gains_gpu(beam_gains_t *d_component_beam_gains,
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
                          beam_gains_t *d_component_beam_gains,
                          double *d_freqs);

extern void wrapper_calculate_analytic_dipole_beam_gpu(int num_components,
                          components_t *components,
                          beam_gains_t *d_component_beam_gains,
                          double *d_freqs, woden_settings_t *woden_settings);

extern void wrapper_run_hyperbeam_gpu(int num_components,
                          components_t *components, beam_settings_t *beam_settings,
                          int num_beams, int parallactic,
                          double *reordered_azs, double *reordered_zas,
                          beam_gains_t *d_component_beam_gains,
                          double *d_freqs, woden_settings_t *woden_settings);

extern void wrapper_calculate_RTS_MWA_analytic_beam_gpu(int num_components,
                          components_t *components, int norm,
                          beam_gains_t *d_component_beam_gains,
                          double *d_freqs, woden_settings_t *woden_settings);

extern void calc_autos_gpu(components_t *d_components,
                               beam_settings_t *beam_settings,
                               beam_gains_t *d_component_beam_gains,
                               visibility_set_t *d_visibility_set,
                               woden_settings_t *woden_settings,
                               int num_components, int use_twobeams);

#endif