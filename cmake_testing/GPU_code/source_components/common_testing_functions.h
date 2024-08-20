#pragma once
#include "constants.h"
#include "woden_struct_defs.h"
#include "woden_precision_defs.h"
#include <math.h>
#include <unity.h>

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
