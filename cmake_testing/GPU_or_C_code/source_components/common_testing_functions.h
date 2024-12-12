#pragma once
#include "constants.h"
#include "woden_struct_defs.h"
#include "woden_precision_defs.h"
#include <math.h>
#include <unity.h>
#include <string.h>

//Linear interpolation between list flux values - go through the list
//and find out which points the desired frequency lies between, and then
//interpolate between the fluxes for that point
void extrap_stokes_list_flux(int *arr_num_list_values, int *list_start_indexes,
            user_precision_t *list_fluxes, double *list_freqs,
           double *extrap_freqs, int iFluxComp, int iFreq,
           double * extrap_flux);


//Extrapolate the power law flux of the iPowerComp component in components
//and the iFreq frequency in extrap_freqs
void extrap_stokes_power_law(user_precision_t *power_refs, user_precision_t *power_SIs,
           double *extrap_freqs, int iPowerComp, int iFreq,
           double * flux);

//Extrapolate the curved power law flux of the iCurveComp component in components
//and the iFreq frequency in extrap_freqs
void extrap_stokes_curved_power_law(user_precision_t *curve_refs,
           user_precision_t *curve_SIs, user_precision_t *curve_qs,
           double *extrap_freqs, int iCurveComp, int iFreq,
           double * flux);

//Given the components in *comps, extrapolate the expected fluxes
//using the CPU code to compare to the GPU
void CPU_extrapolate_fluxes_in_components(components_t *comps, int num_powers,
  int num_curves, int num_lists,
  double *extrap_freqs, int num_extrap_freqs,
  double *expec_flux_I, double *expec_flux_Q,
  double *expec_flux_U, double *expec_flux_V);


//Stick the components into a source_t struct
source_t * put_components_into_source(components_t components,
                                      e_component_type comptype,
                                      int n_powers, int n_curves, int n_lists,
                                      int num_shape_coeffs);

void copy_outputs_source_component_common_cpu(int num_of_each_flux_type,
           source_t *mem_chunked_source, beam_gains_t *mem_beam_gains,
           woden_settings_t *woden_settings,
           beam_settings_t *beam_settings,
           user_precision_complex_t *gxs, user_precision_complex_t *Dxs,
           user_precision_complex_t *Dys, user_precision_complex_t *gys,
           user_precision_t *extrap_flux_I, user_precision_t *extrap_flux_Q,
           user_precision_t *extrap_flux_U, user_precision_t *extrap_flux_V,
           double *ls, double *ms, double *ns,
           e_component_type comptype);