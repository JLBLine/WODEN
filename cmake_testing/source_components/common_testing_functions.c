#include "common_testing_functions.h"

//Takes selected indexes and does the linear interpolation between them
double calc_gradient_extrap_list(user_precision_t *list_fluxes,
          double *list_freqs, double desired_freq, int low_ind_1, int low_ind_2) {

  double gradient = (list_fluxes[low_ind_2] - list_fluxes[low_ind_1]) / (list_freqs[low_ind_2] - list_freqs[low_ind_1]);
  double extrap_flux = list_fluxes[low_ind_1] + gradient*(desired_freq - list_freqs[low_ind_1]);

  return extrap_flux;
}


//Linear interpolation between list flux values - go through the list
//and find out which points the desired frequency lies between, and then
//interpolate between the fluxes for that point
void extrap_stokes_list_flux(components_t *components,
           double *extrap_freqs, int iFluxComp, int iFreq,
           double * flux_I, double * flux_Q,
           double * flux_U, double * flux_V){

  int num_list_values = components->num_list_values[iFluxComp];
  int list_start_ind = components->list_start_indexes[iFluxComp];

  double extrap_freq = extrap_freqs[iFreq];

  int low_ind_1 = 0;
  int low_ind_2 = 0;

  double low_val_1 = 1e16;
  double low_val_2 = 1e16;

  double ref_freq;
  double abs_diff_freq;

  for (int i = 0; i < num_list_values; i++) {
    ref_freq = components->list_freqs[list_start_ind + i];
    abs_diff_freq = abs(ref_freq - extrap_freq);

    if (abs_diff_freq < low_val_1) {
      low_val_1 = abs_diff_freq;
      low_ind_1 = i;
    }
  }

  for (int i = 0; i < num_list_values; i++) {
    ref_freq = components->list_freqs[list_start_ind + i];
    abs_diff_freq = abs(ref_freq - extrap_freq);

    if (abs_diff_freq < low_val_2 && i != low_ind_1) {
      low_val_2 = abs_diff_freq;
      low_ind_2 = i;
    }
  }

  if (low_ind_1 == low_ind_2){
      low_ind_2 = num_list_values - 1;
      low_ind_1 = num_list_values - 2;
  }

  * flux_I = calc_gradient_extrap_list(components->list_stokesI,
            components->list_freqs, extrap_freq,
            list_start_ind + low_ind_1, list_start_ind + low_ind_2);
  * flux_Q = calc_gradient_extrap_list(components->list_stokesQ,
            components->list_freqs, extrap_freq,
            list_start_ind + low_ind_1, list_start_ind + low_ind_2);
  * flux_U = calc_gradient_extrap_list(components->list_stokesU,
            components->list_freqs, extrap_freq,
            list_start_ind + low_ind_1, list_start_ind + low_ind_2);
  * flux_V = calc_gradient_extrap_list(components->list_stokesV,
            components->list_freqs, extrap_freq,
            list_start_ind + low_ind_1, list_start_ind + low_ind_2);

}

//Extrapolate the power law flux of the iPowerComp component in components
//and the iFreq frequency in extrap_freqs
void extrap_stokes_power_law(components_t *components,
           double *extrap_freqs, int iPowerComp, int iFreq,
           double * flux_I, double * flux_Q,
           double * flux_U, double * flux_V) {

  double extrap_freq = extrap_freqs[iFreq];
  double flux_ratio = pow(extrap_freq / components->power_ref_freqs[iPowerComp], components->power_SIs[iPowerComp]);

  * flux_I = components->power_ref_stokesI[iPowerComp] * flux_ratio;
  * flux_Q = components->power_ref_stokesQ[iPowerComp] * flux_ratio;
  * flux_U = components->power_ref_stokesU[iPowerComp] * flux_ratio;
  * flux_V = components->power_ref_stokesV[iPowerComp] * flux_ratio;

}

//Extrapolate the curved power law flux of the iCurveComp component in components
//and the iFreq frequency in extrap_freqs
void extrap_stokes_curved_power_law(components_t *components,
           double *extrap_freqs, int iCurveComp, int iFreq,
           double * flux_I, double * flux_Q,
           double * flux_U, double * flux_V) {

  // printf("%.1e %.1e %.1e %.1e\n", components->curve_ref_freqs[iCurveComp],
  //                                 components->curve_SIs[iCurveComp],
  //                                 components->curve_qs[iCurveComp],
  //                                 components->curve_ref_stokesI[iCurveComp]);

  double extrap_freq = extrap_freqs[iFreq];
  double si_ratio = pow(extrap_freq / components->curve_ref_freqs[iCurveComp],
                                             components->curve_SIs[iCurveComp]);

  double exp_numer = exp((double)components->curve_qs[iCurveComp]*log(extrap_freq)*log(extrap_freq));
  double exp_denom = exp((double)components->curve_qs[iCurveComp]*log(components->curve_ref_freqs[iCurveComp])*log(components->curve_ref_freqs[iCurveComp]));

  double exp_ratio = exp_numer / exp_denom;

  // printf("%d %.4f %.4f %.4f\n",iCurveComp, si_ratio, exp_ratio );

  * flux_I = components->curve_ref_stokesI[iCurveComp] * exp_ratio * si_ratio;
  * flux_Q = components->curve_ref_stokesQ[iCurveComp] * exp_ratio * si_ratio;
  * flux_U = components->curve_ref_stokesU[iCurveComp] * exp_ratio * si_ratio;
  * flux_V = components->curve_ref_stokesV[iCurveComp] * exp_ratio * si_ratio;

  // if (iCurveComp == 9) {
  //   printf("In C %d %.1f %.1f %.1f %.1f %.5f\n", iFreq,
  //                    components->curve_ref_freqs[iCurveComp],
  //                    components->curve_SIs[iCurveComp],
  //                    components->curve_qs[iCurveComp],
  //                    components->curve_ref_stokesI[iCurveComp],
  //                    * flux_I);
  // }

}



void CPU_extrapolate_fluxes_in_components(components_t *comps, int num_powers,
  int num_curves, int num_lists,
  double *extrap_freqs, int num_extrap_freqs,
  double *expec_flux_I, double *expec_flux_Q,
  double *expec_flux_U, double *expec_flux_V){

  double flux_I;
  double flux_Q;
  double flux_U;
  double flux_V;

  //Fill values with what should have been found for POWER_LAW
  int ind = 0;
  for (int comp = 0; comp < num_powers; comp++) {
    for (int extrap = 0; extrap < num_extrap_freqs; extrap++) {

      extrap_stokes_power_law(comps, extrap_freqs, comp, extrap,
                                            &flux_I, &flux_Q, &flux_U, &flux_V);

      expec_flux_I[ind] = flux_I;
      expec_flux_Q[ind] = flux_Q;
      expec_flux_U[ind] = flux_U;
      expec_flux_V[ind] = flux_V;
      // printf("ind, flux %d %.3f\n",ind, ref_stokesI[comp] * flux_ratio );
      ind ++;

      // if (extrap == 0) {
      //     printf("In C %d %.5f\n",comp, flux_I);
      // }

    }
  }

  //Fill values with what should have been found for CURVED_POWER_LAW
  ind = 0;
  for (int comp = 0; comp < num_curves; comp++) {
    for (int extrap = 0; extrap < num_extrap_freqs; extrap++) {

      extrap_stokes_curved_power_law(comps, extrap_freqs, comp, extrap,
                                            &flux_I, &flux_Q, &flux_U, &flux_V);

      expec_flux_I[ind + num_extrap_freqs*num_powers] = flux_I;
      expec_flux_Q[ind + num_extrap_freqs*num_powers] = flux_Q;
      expec_flux_U[ind + num_extrap_freqs*num_powers] = flux_U;
      expec_flux_V[ind + num_extrap_freqs*num_powers] = flux_V;
      // printf("ind, flux %d %.3f\n",ind + num_extrap_freqs*num_powers,
      //                              expec_flux_I[ind + num_extrap_freqs*num_powers]);
      ind ++;

      // if (extrap == 0) {
      //     printf("In C %d %.5f\n",15 + comp, flux_I);
      // }
    }
  }

  //Fill values with what should have been found for LIST fluxes
  ind = 0;
  for (int comp = 0; comp < num_lists; comp++) {
    for (int extrap = 0; extrap < num_extrap_freqs; extrap++) {


      extrap_stokes_list_flux(comps, extrap_freqs,
                     comp, extrap,
                     &flux_I, &flux_Q, &flux_U, &flux_V);

      expec_flux_I[ind + num_extrap_freqs*(num_powers + num_curves)] = flux_I;
      expec_flux_Q[ind + num_extrap_freqs*(num_powers + num_curves)] = flux_Q;
      expec_flux_U[ind + num_extrap_freqs*(num_powers + num_curves)] = flux_U;
      expec_flux_V[ind + num_extrap_freqs*(num_powers + num_curves)] = flux_V;
      // printf("ind, flux %d %.3f\n",ind + num_extrap_freqs*num_powers,
      //                              expec_flux_I[ind + num_extrap_freqs*num_powers]);

      ind ++;
    }
  }
}