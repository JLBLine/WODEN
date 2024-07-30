#include "common_testing_functions.h"

//Takes selected indexes and does the linear interpolation between them
double calc_gradient_extrap_list(user_precision_t *list_fluxes,
          double *list_freqs, double desired_freq, int low_ind_1, int low_ind_2) {

  double gradient;
  double extrap_flux;

  //If one is negative, do interpolation in linear space
  if (list_fluxes[low_ind_1] <= 0 || list_fluxes[low_ind_2] <= 0) {
    gradient = (list_fluxes[low_ind_2] - list_fluxes[low_ind_1]) / (list_freqs[low_ind_2] - list_freqs[low_ind_1]);
    extrap_flux = list_fluxes[low_ind_1] + gradient*(desired_freq - list_freqs[low_ind_1]);
  }

  else {

    double logflux1, logflux2, logfreq1, logfreq2, log_des_freq;

    // printf("what a do %d %d %.3e %.3e %.3e %.3e %.3e\n",low_ind_1, low_ind_2,
    // list_fluxes[low_ind_1],
    // list_fluxes[low_ind_2], list_freqs[low_ind_1], list_freqs[low_ind_2], desired_freq );

    logflux1 = log10((double)list_fluxes[low_ind_1]);
    logflux2 = log10((double)list_fluxes[low_ind_2]);
    logfreq1 = log10(list_freqs[low_ind_1]);
    logfreq2 = log10(list_freqs[low_ind_2]);
    log_des_freq = log10(desired_freq);

    gradient = (logflux2 - logflux1) / (logfreq2 - logfreq1);
    extrap_flux = logflux1 + gradient*(log_des_freq - logfreq1);

    extrap_flux = pow(10, extrap_flux);

  }

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

  int low_ind_1 = -1;
  int low_ind_2 = -1;

  double low_val_1 = 1e16;
  // double low_val_2 = 1e16;

  double ref_freq;
  double abs_diff_freq;

  if (num_list_values == 1) {
    * flux_I = components->list_stokesI[list_start_ind];
    * flux_Q = components->list_stokesQ[list_start_ind];
    * flux_U = components->list_stokesU[list_start_ind];
    * flux_V = components->list_stokesV[list_start_ind];
    return;
  }

  //First loop finds the absolute closest frequency
  for (int i = 0; i < num_list_values; i++) {
    ref_freq = components->list_freqs[list_start_ind + i];
    abs_diff_freq = abs(ref_freq - extrap_freq);

    if (abs_diff_freq < low_val_1) {
      low_val_1 = abs_diff_freq;
      low_ind_1 = i;
    }
  }

  //Depending on the closest frequency, we either want to search above or
  //below the target frequency to find points either side of the target freq

  //We happen to need the reference frequency; just return the refs
  if (components->list_freqs[list_start_ind + low_ind_1] == extrap_freq) {
    // if (iFluxComp == 5 && iFreq == 13){
      // printf("We are heeeeere iFreq %d\n", iFreq);
    // }
    * flux_I = components->list_stokesI[list_start_ind + low_ind_1];
    * flux_Q = components->list_stokesQ[list_start_ind + low_ind_1];
    * flux_U = components->list_stokesU[list_start_ind + low_ind_1];
    * flux_V = components->list_stokesV[list_start_ind + low_ind_1];
    return;
  }
  else {
    //The closest freq is the first index, so set the second index to the second
    if (low_ind_1 == 0) {
      low_ind_2 = 1;
    }
    //closest freq the highest list entry - set second index to one below
    //(order of indexes doesn't matter, as the calculated gradient is pos/neg
    //as needed)
    else if (low_ind_1 == num_list_values - 1){
      low_ind_2 = low_ind_1 - 1;
    }
    else {
      //closest freq is higher than desired - set second index to one below
      //(order of indexes doesn't matter, as the calculated gradient is pos/neg
      //as needed)
      if (components->list_freqs[list_start_ind + low_ind_1] > extrap_freq){
        low_ind_2 = low_ind_1 - 1;
      }
      else {
        low_ind_2 = low_ind_1 + 1;
      }
        //We are extrapolating to a frequency that is lower than all list entries
        //so just stick low_ind_2 to one above low_ind_1
    }
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
void extrap_stokes_power_law(user_precision_t *power_refs, user_precision_t *power_SIs,
           double *extrap_freqs, int iPowerComp, int iFreq,
           double * flux) {

  double extrap_freq = extrap_freqs[iFreq];
  double flux_ratio = pow(extrap_freq / REF_FREQ, power_SIs[iPowerComp]);

  * flux = power_refs[iPowerComp] * flux_ratio;
}

//Extrapolate the curved power law flux of the iCurveComp component in components
//and the iFreq frequency in extrap_freqs
void extrap_stokes_curved_power_law(user_precision_t *curve_refs,
           user_precision_t *curve_SIs, user_precision_t *curve_qs,
           double *extrap_freqs, int iCurveComp, int iFreq,
           double * flux) {

  // printf("%.1e %.1e %.1e %.1e\n", components->curve_ref_freqs[iCurveComp],
  //                                 components->curve_SIs[iCurveComp],
  //                                 components->curve_qs[iCurveComp],
  //                                 components->curve_ref_stokesI[iCurveComp]);


  double extrap_freq = extrap_freqs[iFreq];
  // double ref_freq = components->curve_ref_freqs[iCurveComp];

  double si_ratio = pow(extrap_freq / REF_FREQ, curve_SIs[iCurveComp]);
  double log_ratio = log(extrap_freq / REF_FREQ);

  double exp_ratio = exp((double)curve_qs[iCurveComp]*log_ratio*log_ratio);

  * flux = curve_refs[iCurveComp] * exp_ratio * si_ratio;
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

  //Fill values with what should have been found for Stokes I POWER_LAW
  //Use this opportunity the set the polarisation values to 0
  //We'll fill in what should have values later
  int ind = 0;
  for (int comp = 0; comp < num_powers; comp++) {
    for (int extrap = 0; extrap < num_extrap_freqs; extrap++) {

      extrap_stokes_power_law(comps->power_ref_stokesI, comps->power_SIs,
                              extrap_freqs, comp, extrap, &flux_I);

      expec_flux_I[ind] = flux_I;
      expec_flux_Q[ind] = 0.0;
      expec_flux_U[ind] = 0.0;
      expec_flux_V[ind] = 0.0;
      ind ++;

    }
  }

  //Fill values with what should have been found for CURVED_POWER_LAW
  ind = 0;
  for (int comp = 0; comp < num_curves; comp++) {
    for (int extrap = 0; extrap < num_extrap_freqs; extrap++) {

      extrap_stokes_curved_power_law(comps->curve_ref_stokesI, comps->curve_SIs,
                                     comps->curve_qs ,
                                      extrap_freqs, comp, extrap, &flux_I);

      expec_flux_I[ind + num_extrap_freqs*num_powers] = flux_I;
      expec_flux_Q[ind + num_extrap_freqs*num_powers] = 0.0;
      expec_flux_U[ind + num_extrap_freqs*num_powers] = 0.0;
      expec_flux_V[ind + num_extrap_freqs*num_powers] = 0.0;
      ind ++;
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
      // expec_flux_Q[ind + num_extrap_freqs*(num_powers + num_curves)] = flux_Q;
      // expec_flux_U[ind + num_extrap_freqs*(num_powers + num_curves)] = flux_U;
      // expec_flux_V[ind + num_extrap_freqs*(num_powers + num_curves)] = flux_V;
      expec_flux_Q[ind + num_extrap_freqs*(num_powers + num_curves)] = 0.0;
      expec_flux_U[ind + num_extrap_freqs*(num_powers + num_curves)] = 0.0;
      expec_flux_V[ind + num_extrap_freqs*(num_powers + num_curves)] = 0.0;
      ind ++;
    }
  }
  
  //Polarisation times
  int vind;
  int lind;
  int vcomp;
  int lcomp;

  for (int comp = 0; comp < num_powers; comp++) {
    for (int extrap = 0; extrap < num_extrap_freqs; extrap++) {

      extrap_stokes_power_law(comps->stokesV_power_ref_flux,
                              comps->stokesV_power_SIs,
                              extrap_freqs, comp, extrap, &flux_V);

      vcomp = comps->stokesV_power_comp_inds[comp];
      vind = vcomp*num_extrap_freqs + extrap;
      expec_flux_V[vind] = flux_V;
      // expec_flux_V[vind] = 0.0;

      extrap_stokes_power_law(comps->linpol_power_ref_flux,
                              comps->linpol_power_SIs,
                              extrap_freqs, comp, extrap, &flux_Q);

      lcomp = comps->linpol_power_comp_inds[comp];
      lind = lcomp*num_extrap_freqs + extrap;
      expec_flux_Q[lind] = flux_Q;

      ind ++;

    }
  }

  for (int comp = 0; comp < num_curves; comp++) {
    for (int extrap = 0; extrap < num_extrap_freqs; extrap++) {

      extrap_stokes_curved_power_law(comps->stokesV_curve_ref_flux,
                              comps->stokesV_curve_SIs,
                              comps->stokesV_curve_qs,
                              extrap_freqs, comp, extrap, &flux_V);

      vcomp = comps->stokesV_curve_comp_inds[comp];
      vind = vcomp*num_extrap_freqs + extrap;
      expec_flux_V[vind] = flux_V;

      extrap_stokes_curved_power_law(comps->linpol_curve_ref_flux,
                              comps->linpol_curve_SIs,
                              comps->linpol_curve_qs,
                              extrap_freqs, comp, extrap, &flux_Q);

      lcomp = comps->linpol_curve_comp_inds[comp];
      lind = lcomp*num_extrap_freqs + extrap;
      expec_flux_Q[lind] = flux_Q;


      ind ++;

    }
  }

  user_precision_t pol_frac;
  for (int comp = 0; comp < comps->n_stokesV_pol_frac; comp++) {
    for (int extrap = 0; extrap < num_extrap_freqs; extrap++) {

      vcomp = comps->stokesV_pol_frac_comp_inds[comp];
      vind = vcomp*num_extrap_freqs + extrap;

      pol_frac = comps->stokesV_pol_fracs[comp];
      expec_flux_V[vind] = expec_flux_I[vind]*pol_frac;

      lcomp = comps->linpol_pol_frac_comp_inds[comp];
      lind = lcomp*num_extrap_freqs + extrap;

      pol_frac = comps->linpol_pol_fracs[comp];

      // printf("lind, pol_frac, expec_flux_I[lind] %d %f %f\n", lind, pol_frac, expec_flux_I[lind]);

      expec_flux_Q[lind] = expec_flux_I[lind]*pol_frac;

      ind ++;

    }
  }
}
