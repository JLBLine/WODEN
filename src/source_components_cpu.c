#define _GNU_SOURCE  // Enable GNU extensions
#include <math.h>
#include "source_components_cpu.h"
// #include <omp.h>

//_GNU_SOURCE above means we enable GNU extensions, which includes sincos
//function. sincos is in theory faster than sin and cos separately

user_precision_complex_t calc_measurement_equation_cpu(user_precision_t u,
           user_precision_t v, user_precision_t w,
           double l, double m, double n, double offset){

  user_precision_complex_t visi;

  double temp = 2*M_PI*( u*l + v*m + w*(n-1) + offset);
  double temp_im, temp_re;

  sincos(temp, &temp_im, &temp_re);

  visi = temp_re + I*temp_im;

  return visi;
}

void calc_measurement_equation_arrays_cpu(int num_cross, int num_components,
           user_precision_t *us, user_precision_t *vs, user_precision_t *ws,
           double *ls, double *ms, double *ns, user_precision_complex_t *visis){

  user_precision_complex_t visi;

  double u, v, w;
  double l, m, n;

  for (int iBaseline = 0; iBaseline < num_cross; iBaseline++) {

    u = (double)us[iBaseline];
    v = (double)vs[iBaseline];
    w = (double)ws[iBaseline];

    for (int iComponent = 0; iComponent < num_components; iComponent++) {
      l = ls[iComponent];
      m = ms[iComponent];
      n = ns[iComponent];

      visi = calc_measurement_equation_cpu(u, v, w, l, m, n, 0);

      visis[num_components*iBaseline + iComponent] = visi;
    }
  }
}


void apply_beam_gains_stokesI_on_cardinal_cpu(user_precision_complex_t g1x, user_precision_complex_t D1x,
          user_precision_complex_t D1y, user_precision_complex_t g1y,
          user_precision_complex_t g2x, user_precision_complex_t D2x,
          user_precision_complex_t D2y, user_precision_complex_t g2y,
          user_precision_t flux_I, user_precision_complex_t visi_component,
          user_precision_complex_t * visi_XX, user_precision_complex_t * visi_XY,
          user_precision_complex_t * visi_YX, user_precision_complex_t * visi_YY) {

  //Conjugate the second beam gains
  user_precision_complex_t g2x_conj = conj(g2x);
  user_precision_complex_t D2x_conj = conj(D2x);
  user_precision_complex_t D2y_conj = conj(D2y);
  user_precision_complex_t g2y_conj = conj(g2y);

  user_precision_complex_t visi_I = (flux_I + I*0.0)*visi_component;
  user_precision_complex_t this_XX;
  user_precision_complex_t this_XY;
  user_precision_complex_t this_YX;
  user_precision_complex_t this_YY;

  this_XX = (g1x*g2x_conj + D1x*D2x_conj)*visi_I;
  this_XY = (g1x*D2y_conj + D1x*g2y_conj)*visi_I;
  this_YX = (D1y*g2x_conj + g1y*D2x_conj)*visi_I;
  this_YY = (D1y*D2y_conj + g1y*g2y_conj)*visi_I;

  * visi_XX = this_XX;
  * visi_XY = this_XY;
  * visi_YX = this_YX;
  * visi_YY = this_YY;

}

void apply_beam_gains_stokesIQUV_on_cardinal_cpu(user_precision_complex_t g1x, user_precision_complex_t D1x,
          user_precision_complex_t D1y, user_precision_complex_t g1y,
          user_precision_complex_t g2x, user_precision_complex_t D2x,
          user_precision_complex_t D2y, user_precision_complex_t g2y,
          user_precision_t flux_I, user_precision_t flux_Q,
          user_precision_t flux_U, user_precision_t flux_V,
          user_precision_complex_t visi_component,
          user_precision_complex_t * visi_XX, user_precision_complex_t * visi_XY,
          user_precision_complex_t * visi_YX, user_precision_complex_t * visi_YY) {

  //Conjugate the second beam gains
  user_precision_complex_t g2x_conj = conj(g2x);
  user_precision_complex_t D2x_conj = conj(D2x);
  user_precision_complex_t D2y_conj = conj(D2y);
  user_precision_complex_t g2y_conj = conj(g2y);

  //Create the Stokes visibilities
  user_precision_complex_t visi_I = (flux_I + I*0.0)*visi_component;
  user_precision_complex_t visi_Q = (flux_Q + I*0.0)*visi_component;
  user_precision_complex_t visi_U = (flux_U + I*0.0)*visi_component;
  user_precision_complex_t visi_V = (flux_V + I*0.0)*visi_component;

  user_precision_complex_t this_XX;
  user_precision_complex_t this_XY;
  user_precision_complex_t this_YX;
  user_precision_complex_t this_YY;

  user_precision_complex_t temp1 = g1x*g2x_conj;
  user_precision_complex_t temp2 = D1x*D2x_conj;
  user_precision_complex_t temp3 = g1x*D2x_conj;
  user_precision_complex_t temp4 = D1x*g2x_conj;

  this_XX = (temp1 + temp2)*visi_I;
  this_XX += (temp1 - temp2)*visi_Q;
  this_XX += (temp3 + temp4)*visi_U;
  this_XX += ((0.0 + I*1.0)*visi_V)*(temp3 - temp4);

  temp1 = g1x*D2y_conj;
  temp2 = D1x*g2y_conj;
  temp3 = g1x*g2y_conj;
  temp4 = D1x*D2y_conj;

  this_XY = (temp1 + temp2)*visi_I;
  this_XY += (temp1 - temp2)*visi_Q;
  this_XY += (temp3 + temp4)*visi_U;
  this_XY += ((0.0 + I*1.0)*visi_V)*(temp3 - temp4);

  temp1 = D1y*g2x_conj;
  temp2 = g1y*D2x_conj;
  temp3 = D1y*D2x_conj;
  temp4 = g1y*g2x_conj;

  this_YX = (temp1 + temp2)*visi_I;
  this_YX += (temp1 - temp2)*visi_Q;
  this_YX += (temp3 + temp4)*visi_U;
  this_YX += ((0.0 + I*1.0)*visi_V)*(temp3 - temp4);

  temp1 = D1y*D2y_conj;
  temp2 = g1y*g2y_conj;
  temp3 = D1y*g2y_conj;
  temp4 = g1y*D2y_conj;

  this_YY = (temp1 + temp2)*visi_I;
  this_YY += (temp1 - temp2)*visi_Q;
  this_YY += (temp3 + temp4)*visi_U;
  this_YY += ((0.0 + I*1.0)*visi_V)*(temp3 - temp4);

  // this_XX = (g1x*g2x_conj + D1x*D2x_conj)*visi_I;
  // this_XX += (g1x*g2x_conj - D1x*D2x_conj)*visi_Q;
  // this_XX += (g1x*D2x_conj + D1x*g2x_conj)*visi_U;
  // this_XX += ((0.0 + I*1.0)*visi_V)*(g1x*D2x_conj - D1x*g2x_conj);

  // this_XY = (g1x*D2y_conj + D1x*g2y_conj)*visi_I;
  // this_XY += (g1x*D2y_conj - D1x*g2y_conj)*visi_Q;
  // this_XY += (g1x*g2y_conj + D1x*D2y_conj)*visi_U;
  // this_XY += ((0.0 + I*1.0)*visi_V)*(g1x*g2y_conj - D1x*D2y_conj);

  // this_YX = (D1y*g2x_conj + g1y*D2x_conj)*visi_I;
  // this_YX += (D1y*g2x_conj - g1y*D2x_conj)*visi_Q;
  // this_YX += (D1y*D2x_conj + g1y*g2x_conj)*visi_U;
  // this_YX += ((0.0 + I*1.0)*visi_V)*(D1y*D2x_conj - g1y*g2x_conj);

  // this_YY = (D1y*D2y_conj + g1y*g2y_conj)*visi_I;
  // this_YY += (D1y*D2y_conj - g1y*g2y_conj)*visi_Q;
  // this_YY += (D1y*g2y_conj + g1y*D2y_conj)*visi_U;
  // this_YY += ((0.0 + I*1.0)*visi_V)*(D1y*g2y_conj - g1y*D2y_conj);

  * visi_XX = this_XX;
  * visi_XY = this_XY;
  * visi_YX = this_YX;
  * visi_YY = this_YY;

}


void apply_beam_gains_arrays_cpu(int num_gains,
          user_precision_complex_t *g1xs, user_precision_complex_t *D1xs,
          user_precision_complex_t *D1ys, user_precision_complex_t *g1ys,
          user_precision_complex_t *g2xs, user_precision_complex_t *D2xs,
          user_precision_complex_t *D2ys, user_precision_complex_t *g2ys,
          user_precision_t *flux_Is, user_precision_t *flux_Qs,
          user_precision_t *flux_Us, user_precision_t *flux_Vs,
          user_precision_complex_t *visi_components,
          user_precision_complex_t *visi_XXs, user_precision_complex_t *visi_XYs,
          user_precision_complex_t *visi_YXs, user_precision_complex_t *visi_YYs,
          int off_cardinal_dipoles, int do_QUV){

  user_precision_complex_t this_XX, this_XY, this_YX, this_YY;

  for (int gain = 0; gain < num_gains; gain++)
  {

    if (do_QUV == 1) {
      if (off_cardinal_dipoles == 1) {
        apply_beam_gains_stokesIQUV_off_cardinal_cpu(g1xs[gain], D1xs[gain],
          D1ys[gain], g1ys[gain], g2xs[gain], D2xs[gain], D2ys[gain], g2ys[gain],
          flux_Is[gain], flux_Qs[gain], flux_Us[gain], flux_Vs[gain],
          visi_components[gain], &this_XX, &this_XY, &this_YX, &this_YY);
      } else {
        apply_beam_gains_stokesIQUV_on_cardinal_cpu(g1xs[gain], D1xs[gain],
          D1ys[gain], g1ys[gain], g2xs[gain], D2xs[gain], D2ys[gain], g2ys[gain],
          flux_Is[gain], flux_Qs[gain], flux_Us[gain], flux_Vs[gain],
          visi_components[gain], &this_XX, &this_XY, &this_YX, &this_YY);
      }
    } else {
      if (off_cardinal_dipoles == 1) {
        apply_beam_gains_stokesI_off_cardinal_cpu(g1xs[gain], D1xs[gain],
          D1ys[gain], g1ys[gain], g2xs[gain], D2xs[gain], D2ys[gain], g2ys[gain],
          flux_Is[gain],
          visi_components[gain], &this_XX, &this_XY, &this_YX, &this_YY);
      } else {
        apply_beam_gains_stokesI_on_cardinal_cpu(g1xs[gain], D1xs[gain],
          D1ys[gain], g1ys[gain], g2xs[gain], D2xs[gain], D2ys[gain], g2ys[gain],
          flux_Is[gain],
          visi_components[gain], &this_XX, &this_XY, &this_YX, &this_YY);
          // printf("%d XX: %f + %f i\n", gain, creal(this_XX), cimag(this_XX));
      }
    }
    
    visi_XXs[gain] = this_XX;
    visi_XYs[gain] = this_XY;
    visi_YXs[gain] = this_YX;
    visi_YYs[gain] = this_YY;
  }
}

void apply_beam_gains_stokesI_off_cardinal_cpu(user_precision_complex_t g1x, user_precision_complex_t D1x,
          user_precision_complex_t D1y, user_precision_complex_t g1y,
          user_precision_complex_t g2x, user_precision_complex_t D2x,
          user_precision_complex_t D2y, user_precision_complex_t g2y,
          user_precision_t flux_I, user_precision_complex_t visi_component,
          user_precision_complex_t * visi_XX, user_precision_complex_t * visi_XY,
          user_precision_complex_t * visi_YX, user_precision_complex_t * visi_YY) {

  //Conjugate the second beam gains
  user_precision_complex_t g2x_conj = conj(g2x);
  user_precision_complex_t D2x_conj = conj(D2x);
  user_precision_complex_t D2y_conj = conj(D2y);
  user_precision_complex_t g2y_conj = conj(g2y);

  //Create the Stokes visibilities
  user_precision_complex_t visi_I = (flux_I + I*0.0)*visi_component;

  user_precision_complex_t this_XX;
  user_precision_complex_t this_XY;
  user_precision_complex_t this_YX;
  user_precision_complex_t this_YY;

  this_XX = (g1x*g2x_conj + D1x*D2x_conj)*visi_I;
  this_XY = (g1x*D2y_conj + D1x*g2y_conj)*visi_I;
  this_YX = (D1y*g2x_conj + g1y*D2x_conj)*visi_I;
  this_YY = (D1y*D2y_conj + g1y*g2y_conj)*visi_I;

  * visi_XX = this_XX;
  * visi_XY = this_XY;
  * visi_YX = this_YX;
  * visi_YY = this_YY;

}

void apply_beam_gains_stokesIQUV_off_cardinal_cpu(user_precision_complex_t g1x, user_precision_complex_t D1x,
          user_precision_complex_t D1y, user_precision_complex_t g1y,
          user_precision_complex_t g2x, user_precision_complex_t D2x,
          user_precision_complex_t D2y, user_precision_complex_t g2y,
          user_precision_t flux_I, user_precision_t flux_Q,
          user_precision_t flux_U, user_precision_t flux_V,
          user_precision_complex_t visi_component,
          user_precision_complex_t * visi_XX, user_precision_complex_t * visi_XY,
          user_precision_complex_t * visi_YX, user_precision_complex_t * visi_YY) {

  //Conjugate the second beam gains
  user_precision_complex_t g2x_conj = conj(g2x);
  user_precision_complex_t D2x_conj = conj(D2x);
  user_precision_complex_t D2y_conj = conj(D2y);
  user_precision_complex_t g2y_conj = conj(g2y);

  //Create the Stokes visibilities
  user_precision_complex_t visi_I = (flux_I + I*0.0)*visi_component;
  user_precision_complex_t visi_Q = (flux_Q + I*0.0)*visi_component;
  user_precision_complex_t visi_U = (flux_U + I*0.0)*visi_component;
  user_precision_complex_t visi_V = (flux_V + I*0.0)*visi_component;

  user_precision_complex_t this_XX;
  user_precision_complex_t this_XY;
  user_precision_complex_t this_YX;
  user_precision_complex_t this_YY;

  this_XX = (g1x*g2x_conj + D1x*D2x_conj)*visi_I;
  this_XX -= (g1x*D2x_conj + D1x*g2x_conj)*visi_Q;
  this_XX += (g1x*g2x_conj + D1x*D2x_conj)*visi_U;
  this_XX += ((0.0 + I*1.0)*visi_V)*(g1x*D2x_conj - D1x*g2x_conj);

  this_XY = (g1x*D2y_conj + D1x*g2y_conj)*visi_I;
  this_XY -= (g1x*g2y_conj + D1x*D2y_conj)*visi_Q;
  this_XY += (g1x*D2y_conj - D1x*g2y_conj)*visi_U;
  this_XY += ((0.0 + I*1.0)*visi_V)*(g1x*g2y_conj - D1x*D2y_conj);

  this_YX = (D1y*g2x_conj + g1y*D2x_conj)*visi_I;
  this_YX -= (D1y*D2x_conj + g1y*g2x_conj)*visi_Q;
  this_YX += (D1y*g2x_conj - g1y*D2x_conj)*visi_U;
  this_YX += ((0.0 + I*1.0)*visi_V)*(D1y*D2x_conj - g1y*g2x_conj);

  this_YY = (D1y*D2y_conj + g1y*g2y_conj)*visi_I;
  this_YY -= (D1y*g2y_conj + g1y*D2y_conj)*visi_Q;
  this_YY += (D1y*D2y_conj - g1y*g2y_conj)*visi_U;
  this_YY += ((0.0 + I*1.0)*visi_V)*(D1y*g2y_conj - g1y*D2y_conj);

  * visi_XX = this_XX;
  * visi_XY = this_XY;
  * visi_YX = this_YX;
  * visi_YY = this_YY;

}

void malloc_extrapolated_flux_arrays_cpu(components_t *components, int num_comps,
                                     int num_freqs){
  components->extrap_stokesI = malloc(num_comps*num_freqs*sizeof(user_precision_t));

  if (components->do_QUV == 1)
  {
      // printf("Doing full polarisation\n");
      components->extrap_stokesQ = malloc(num_comps*num_freqs*sizeof(user_precision_t));
      components->extrap_stokesU = malloc(num_comps*num_freqs*sizeof(user_precision_t));
      components->extrap_stokesV = malloc(num_comps*num_freqs*sizeof(user_precision_t));

      //set everything to zero at not every component has to have full polarisation

      for (int extra_ind = 0; extra_ind < num_comps*num_freqs; extra_ind++)
      {
          components->extrap_stokesQ[extra_ind] = 0.0;
          components->extrap_stokesU[extra_ind] = 0.0;
          components->extrap_stokesV[extra_ind] = 0.0;
      }
  }
}

void free_extrapolated_flux_arrays_cpu(components_t *components){
  free(components->extrap_stokesI);

  if (components->do_QUV) {
    free(components->extrap_stokesQ);
    free(components->extrap_stokesU);
    free(components->extrap_stokesV);
  }
}


void extrap_stokes_power_law_cpu(user_precision_t ref_flux, user_precision_t SI,
                             double extrap_freq,  user_precision_t * extrap_flux){

  user_precision_t flux_ratio = pow(extrap_freq / REF_FREQ, SI);

  *extrap_flux = ref_flux*flux_ratio;
}


void extrap_power_laws_stokesI_cpu(components_t components, int num_components,
                                   double *extrap_freqs, int num_extrap_freqs){

  user_precision_t flux_I;

  for (int iFluxComp = 0; iFluxComp < num_components; iFluxComp++) {
    for (int iFreq = 0; iFreq < num_extrap_freqs; iFreq++) {

      extrap_stokes_power_law_cpu(components.power_ref_stokesI[iFluxComp],
                              components.power_SIs[iFluxComp],
                              extrap_freqs[iFreq], &flux_I);

      int iComponent = components.power_comp_inds[iFluxComp];
      int extrap_ind = num_extrap_freqs*iComponent + iFreq;

      components.extrap_stokesI[extrap_ind] = flux_I;
    }
  }
}

void extrap_power_laws_stokesV_cpu(components_t components, double *extrap_freqs,
                                   int num_extrap_freqs){

  user_precision_t flux_V;

  for (int iFluxComp = 0; iFluxComp < components.n_stokesV_power; iFluxComp++) {
    for (int iFreq = 0; iFreq < num_extrap_freqs; iFreq++) {

      extrap_stokes_power_law_cpu(components.stokesV_power_ref_flux[iFluxComp],
                              components.stokesV_power_SIs[iFluxComp],
                              extrap_freqs[iFreq], &flux_V);

      int iComponent = components.stokesV_power_comp_inds[iFluxComp];
      int extrap_ind = num_extrap_freqs*iComponent + iFreq;

      components.extrap_stokesV[extrap_ind] = flux_V;
    }
  }
}

void extrap_power_laws_linpol_cpu(components_t components, double *extrap_freqs,
                                   int num_extrap_freqs){

  user_precision_t flux_linpol;

  for (int iFluxComp = 0; iFluxComp < components.n_linpol_power; iFluxComp++) {
    for (int iFreq = 0; iFreq < num_extrap_freqs; iFreq++) {

      extrap_stokes_power_law_cpu(components.linpol_power_ref_flux[iFluxComp],
                              components.linpol_power_SIs[iFluxComp],
                              extrap_freqs[iFreq], &flux_linpol);

      int iComponent = components.linpol_power_comp_inds[iFluxComp];
      int extrap_ind = num_extrap_freqs*iComponent + iFreq;

      components.extrap_stokesQ[extrap_ind] = flux_linpol;
    }
  }
}



void extrap_stokes_curved_power_law_cpu(user_precision_t ref_flux,
           user_precision_t SI, user_precision_t q,
           double extrap_freq, user_precision_t * extrap_flux){

  user_precision_t si_ratio = pow(extrap_freq / REF_FREQ, SI);
  double log_freq_ratio = log(extrap_freq / REF_FREQ);
  double exp_bit = exp(q*log_freq_ratio*log_freq_ratio);

  user_precision_t flux_ratio = si_ratio * exp_bit;

  * extrap_flux = ref_flux * flux_ratio;
}

void extrap_curved_power_laws_stokesI_cpu(components_t components, int num_components,
                                   double *extrap_freqs, int num_extrap_freqs){

  user_precision_t flux_I;

  for (int iFluxComp = 0; iFluxComp < num_components; iFluxComp++) {
    for (int iFreq = 0; iFreq < num_extrap_freqs; iFreq++) {

      extrap_stokes_curved_power_law_cpu(components.curve_ref_stokesI[iFluxComp],
                              components.curve_SIs[iFluxComp],
                              components.curve_qs[iFluxComp],
                              extrap_freqs[iFreq], &flux_I);

      int iComponent = components.curve_comp_inds[iFluxComp];
      int extrap_ind = num_extrap_freqs*iComponent + iFreq;

      components.extrap_stokesI[extrap_ind] = flux_I;
    }
  }
}


void extrap_curved_power_laws_stokesV_cpu(components_t components,
                                    double *extrap_freqs, int num_extrap_freqs){

  user_precision_t flux_V;

  for (int iFluxComp = 0; iFluxComp < components.n_stokesV_curve; iFluxComp++) {
    for (int iFreq = 0; iFreq < num_extrap_freqs; iFreq++) {

      extrap_stokes_curved_power_law_cpu(components.stokesV_curve_ref_flux[iFluxComp],
                              components.stokesV_curve_SIs[iFluxComp],
                              components.stokesV_curve_qs[iFluxComp],
                              extrap_freqs[iFreq], &flux_V);

      int iComponent = components.stokesV_curve_comp_inds[iFluxComp];
      int extrap_ind = num_extrap_freqs*iComponent + iFreq;

      components.extrap_stokesV[extrap_ind] = flux_V;
    }
  }
}

void extrap_curved_power_laws_linpol_cpu(components_t components,
                                    double *extrap_freqs, int num_extrap_freqs){

  user_precision_t flux_linpol;

  for (int iFluxComp = 0; iFluxComp < components.n_linpol_curve; iFluxComp++) {
    for (int iFreq = 0; iFreq < num_extrap_freqs; iFreq++) {

      extrap_stokes_curved_power_law_cpu(components.linpol_curve_ref_flux[iFluxComp],
                              components.linpol_curve_SIs[iFluxComp],
                              components.linpol_curve_qs[iFluxComp],
                              extrap_freqs[iFreq], &flux_linpol);

      int iComponent = components.linpol_curve_comp_inds[iFluxComp];
      int extrap_ind = num_extrap_freqs*iComponent + iFreq;

      components.extrap_stokesQ[extrap_ind] = flux_linpol;
    }
  }
}

double calc_gradient_extrap_list_cpu(user_precision_t *list_fluxes,
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
void extrap_stokes_list_flux_cpu(int *arr_num_list_values, int *list_start_indexes,
            user_precision_t *list_fluxes, double *list_freqs,
           double *extrap_freqs, int iFluxComp, int iFreq,
           user_precision_t * extrap_flux){

  int num_list_values = arr_num_list_values[iFluxComp];
  int list_start_ind = list_start_indexes[iFluxComp];

  double extrap_freq = extrap_freqs[iFreq];

  int low_ind_1 = -1;
  int low_ind_2 = -1;

  double low_val_1 = 1e16;
  // double low_val_2 = 1e16;

  double ref_freq;
  double abs_diff_freq;

  if (num_list_values == 1) {
    * extrap_flux = list_fluxes[list_start_ind];
    return;
  }

  //First loop finds the absolute closest frequency
  for (int i = 0; i < num_list_values; i++) {
    ref_freq = list_freqs[list_start_ind + i];
    abs_diff_freq = abs(ref_freq - extrap_freq);

    if (abs_diff_freq < low_val_1) {
      low_val_1 = abs_diff_freq;
      low_ind_1 = i;
    }
  }

  //Depending on the closest frequency, we either want to search above or
  //below the target frequency to find points either side of the target freq

  //We happen to need the reference frequency; just return the refs
  if (list_freqs[list_start_ind + low_ind_1] == extrap_freq) {
    // if (iFluxComp == 5 && iFreq == 13){
      // printf("We are heeeeere iFreq %d\n", iFreq);
    // }
    * extrap_flux = list_fluxes[list_start_ind + low_ind_1];
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
      if (list_freqs[list_start_ind + low_ind_1] > extrap_freq){
        low_ind_2 = low_ind_1 - 1;
      }
      else {
        low_ind_2 = low_ind_1 + 1;
      }
        //We are extrapolating to a frequency that is lower than all list entries
        //so just stick low_ind_2 to one above low_ind_1
    }
  }

  * extrap_flux = calc_gradient_extrap_list_cpu(list_fluxes,
            list_freqs, extrap_freq,
            list_start_ind + low_ind_1, list_start_ind + low_ind_2);

}



void extrap_stokes_list_flux_arrays_cpu(user_precision_t *list_stokes, double *list_freqs,
                                int *num_list_values, int *list_start_indexes,
                                int *list_comp_inds,
                                int num_extrap_freqs, double *extrap_freqs,
                                int num_comps, user_precision_t *extrap_stokes){

  for (int iFluxComp = 0; iFluxComp < num_comps; iFluxComp++) {
    for (int iFreq = 0; iFreq < num_extrap_freqs; iFreq++) {

      int iComponent = list_comp_inds[iFluxComp];
      int extrap_ind = num_extrap_freqs*iComponent + iFreq;

      extrap_stokes_list_flux_cpu(num_list_values, list_start_indexes,
                                  list_stokes, list_freqs,
                                  extrap_freqs, iFluxComp, iFreq,
                                  &extrap_stokes[extrap_ind]);
    }
  }
}

void extrap_list_fluxes_stokesI_cpu(components_t components,
                                   int n_lists, double *extrap_freqs,
                                   int num_extrap_freqs) {
  
  extrap_stokes_list_flux_arrays_cpu(components.list_stokesI, components.list_freqs,
                       components.num_list_values, components.list_start_indexes,
                       components.list_comp_inds,
                       num_extrap_freqs, extrap_freqs,
                       n_lists, components.extrap_stokesI);

}

void extrap_list_fluxes_stokesV_cpu(components_t components, double *extrap_freqs,
                                   int num_extrap_freqs) {
  
  extrap_stokes_list_flux_arrays_cpu(components.stokesV_list_ref_flux,
                      components.stokesV_list_ref_freqs,
                      components.stokesV_num_list_values,
                      components.stokesV_list_start_indexes,
                      components.stokesV_list_comp_inds,
                      num_extrap_freqs, extrap_freqs,
                      components.n_stokesV_list, components.extrap_stokesV);

}

void extrap_list_fluxes_linpol_cpu(components_t components, double *extrap_freqs,
                                   int num_extrap_freqs) {
  
  extrap_stokes_list_flux_arrays_cpu(components.stokesQ_list_ref_flux,
                      components.stokesQ_list_ref_freqs,
                      components.stokesQ_num_list_values,
                      components.stokesQ_list_start_indexes,
                      components.stokesQ_list_comp_inds,
                      num_extrap_freqs, extrap_freqs,
                      components.n_linpol_list, components.extrap_stokesQ);

  extrap_stokes_list_flux_arrays_cpu(components.stokesU_list_ref_flux,
                      components.stokesU_list_ref_freqs,
                      components.stokesU_num_list_values,
                      components.stokesU_list_start_indexes,
                      components.stokesU_list_comp_inds,
                      num_extrap_freqs, extrap_freqs,
                      components.n_linpol_list, components.extrap_stokesU);

}

void extrap_p_list_fluxes_linpol_cpu(components_t components, double *extrap_freqs,
                                   int num_extrap_freqs) {
  
  extrap_stokes_list_flux_arrays_cpu(components.linpol_p_list_ref_flux,
                      components.linpol_p_list_ref_freqs,
                      components.linpol_p_num_list_values,
                      components.linpol_p_list_start_indexes,
                      components.linpol_p_list_comp_inds,
                      num_extrap_freqs, extrap_freqs,
                      components.n_linpol_p_list, components.extrap_stokesQ);

}


void polarisation_fraction_stokesV_cpu(components_t components, double *extrap_freqs,
                                   int num_extrap_freqs){

  for (int iFluxComp = 0; iFluxComp < components.n_stokesV_pol_frac; iFluxComp++) {
    for (int iFreq = 0; iFreq < num_extrap_freqs; iFreq++) {

      user_precision_t pol_frac;
      pol_frac = components.stokesV_pol_fracs[iFluxComp];

      int iComponent = components.stokesV_pol_frac_comp_inds[iFluxComp];
      int extrap_ind = num_extrap_freqs*iComponent + iFreq;

      components.extrap_stokesV[extrap_ind] = pol_frac*components.extrap_stokesI[extrap_ind];
    }
  }
}

void polarisation_fraction_linpol_cpu(components_t components, double *extrap_freqs,
                                   int num_extrap_freqs){

  for (int iFluxComp = 0; iFluxComp < components.n_linpol_pol_frac; iFluxComp++) {
    for (int iFreq = 0; iFreq < num_extrap_freqs; iFreq++) {

      user_precision_t pol_frac;
      pol_frac = components.linpol_pol_fracs[iFluxComp];

      int iComponent = components.linpol_pol_frac_comp_inds[iFluxComp];
      int extrap_ind = num_extrap_freqs*iComponent + iFreq;

      components.extrap_stokesQ[extrap_ind] = pol_frac*components.extrap_stokesI[extrap_ind];
    }
  }
                                  
}


void apply_rotation_measure_cpu(components_t components,
                                double *extrap_freqs, int num_extrap_freqs) {

  // Start by computing which baseline we're going to do
  
  for (int iFluxComp = 0; iFluxComp < components.n_linpol_angles; iFluxComp++) {
    for (int iFreq = 0; iFreq < num_extrap_freqs; iFreq++) {

      int iComponent = components.linpol_angle_inds[iFluxComp];
      int extrap_ind = num_extrap_freqs*iComponent + iFreq;

      user_precision_t rm = components.rm_values[iFluxComp];
      user_precision_t intr_pol_angle = components.intr_pol_angle[iFluxComp];

      //We should have calculated the linear polarisation flux before, and
      //shoved it in the Stokes Q extrap array. No point in wasting precious
      //memory as we'll immediate overwrite it
      user_precision_t linpol_flux = components.extrap_stokesQ[extrap_ind];

      double wavelength = VELC / extrap_freqs[iFreq];
      double angle = 2*(intr_pol_angle + rm*wavelength*wavelength);

      components.extrap_stokesQ[extrap_ind] = linpol_flux*cos(angle);
      components.extrap_stokesU[extrap_ind] = linpol_flux*sin(angle);

    }
  }
}


void get_beam_gains_cpu(int iBaseline, int iComponent, int num_freqs,
           int num_baselines, int num_components, int num_times, int beamtype,
           user_precision_complex_t *gxs, user_precision_complex_t *Dxs,
           user_precision_complex_t *Dys, user_precision_complex_t *gys,
           user_precision_complex_t * g1x, user_precision_complex_t * D1x,
           user_precision_complex_t * D1y, user_precision_complex_t * g1y,
           user_precision_complex_t * g2x, user_precision_complex_t * D2x,
           user_precision_complex_t * D2y, user_precision_complex_t * g2y){

  int beam_ind = 0;
  int time_ind = 0;
  int freq_ind = 0;

  time_ind = (int)floorf( (float)iBaseline / ((float)num_baselines * (float)num_freqs));
  freq_ind = (int)floorf( ((float)iBaseline - ((float)time_ind*(float)num_baselines * (float)num_freqs)) / (float)num_baselines);
  beam_ind = num_freqs*time_ind*num_components + (num_components*freq_ind) + iComponent;

    //Set gains to one if no beam
  if (beamtype == NO_BEAM) {
    * g1x = 1.0 + I*0.0;
    * g2x = 1.0 + I*0.0;
    * g1y = 1.0 + I*0.0;
    * g2y = 1.0 + I*0.0;
  }

  //Get gains if using a beam
  else {
    * g1x = gxs[beam_ind];
    * g2x = gxs[beam_ind];
    * g1y = gys[beam_ind];
    * g2y = gys[beam_ind];

  }

  // Set leakage to zero if no leakage
  if (beamtype == NO_BEAM || beamtype == GAUSS_BEAM || beamtype == ANALY_DIPOLE) {
    * D1x = 0.0 + I*0.0;
    * D2x = 0.0 + I*0.0;
    * D1y = 0.0 + I*0.0;
    * D2y = 0.0 + I*0.0;
  }
  
  else {
    * D1x = Dxs[beam_ind];
    * D2x = Dxs[beam_ind];
    * D1y = Dys[beam_ind];
    * D2y = Dys[beam_ind];
  }

  // printf("beam_ind g1x %d %.3e %.3e\n", beam_ind, creal(*g1x), cimag(*g1x));
  // printf("beam_ind g2x %d %.3e %.3e\n", beam_ind, creal(*g2x), cimag(*g2x));
  // printf("beam_ind g1y %d %.3e %.3e\n", beam_ind, creal(*g1y), cimag(*g1y));
  // printf("beam_ind g2y %d %.3e %.3e\n", beam_ind, creal(*g2y), cimag(*g2y));

} //end get_beam_gains_cpu


void get_beam_gains_multibeams_cpu(int iBaseline, int iComponent, int num_freqs,
           int num_baselines, int num_components, int num_times, int beamtype,
           user_precision_complex_t *gxs, user_precision_complex_t *Dxs,
           user_precision_complex_t *Dys, user_precision_complex_t *gys,
           int *ant1_to_baseline_map, int *ant2_to_baseline_map,
           user_precision_complex_t * g1x, user_precision_complex_t * D1x,
           user_precision_complex_t * D1y, user_precision_complex_t * g1y,
           user_precision_complex_t * g2x, user_precision_complex_t * D2x,
           user_precision_complex_t * D2y, user_precision_complex_t * g2y){

  int time_ind = 0;
  int freq_ind = 0;

  time_ind = (int)floorf( (float)iBaseline / ((float)num_baselines * (float)num_freqs));
  freq_ind = (int)floorf( ((float)iBaseline - ((float)time_ind*(float)num_baselines * (float)num_freqs)) / (float)num_baselines);
  // beam_ind = num_freqs*time_ind*num_components + (num_components*freq_ind) + iComponent;

  int baseline_ind = iBaseline % num_baselines;

  int ant1_ind = ant1_to_baseline_map[baseline_ind];
  int ant2_ind = ant2_to_baseline_map[baseline_ind];

  int beam1 = ant1_ind*num_freqs*num_components*num_times + num_freqs*time_ind*num_components + (num_components*freq_ind) + iComponent;
  int beam2 = ant2_ind*num_freqs*num_components*num_times + num_freqs*time_ind*num_components + (num_components*freq_ind) + iComponent;

    //Set gains to one if no beam
  if (beamtype == NO_BEAM) {
    * g1x = 1.0 + I*0.0;
    * g2x = 1.0 + I*0.0;
    * g1y = 1.0 + I*0.0;
    * g2y = 1.0 + I*0.0;
  }

  //Get gains if using a beam
  else {
    * g1x = gxs[beam1];
    * g2x = gxs[beam2];
    * g1y = gys[beam1];
    * g2y = gys[beam2];

  }

  //These models have no leakage terms
  if (beamtype == NO_BEAM || beamtype == GAUSS_BEAM || beamtype == ANALY_DIPOLE) {
    * D1x = 0.0 + I*0.0;
    * D2x = 0.0 + I*0.0;
    * D1y = 0.0 + I*0.0;
    * D2y = 0.0 + I*0.0;
  }
  // Set leakage to zero if no leakage
  else {
    * D1x = Dxs[beam1];
    * D2x = Dxs[beam2];
    * D1y = Dys[beam1];
    * D2y = Dys[beam2];
  }

} //end __device__ get_beam_gains_multibeams_cpu

void update_sum_visis_stokesI_cpu(int iBaseline, int iComponent,
    int num_freqs, int num_baselines, int num_components, int num_times,
    int beamtype, int off_cardinal_dipoles,
    user_precision_complex_t *gxs, user_precision_complex_t *Dxs,
    user_precision_complex_t *Dys, user_precision_complex_t *gys,
    int *ant1_to_baseline_map, int *ant2_to_baseline_map, int use_twobeams,
    user_precision_complex_t visi_component,  user_precision_t flux_I,
    user_precision_t *sum_visi_XX_real, user_precision_t *sum_visi_XX_imag,
    user_precision_t *sum_visi_XY_real, user_precision_t *sum_visi_XY_imag,
    user_precision_t *sum_visi_YX_real, user_precision_t *sum_visi_YX_imag,
    user_precision_t *sum_visi_YY_real, user_precision_t *sum_visi_YY_imag){

    user_precision_complex_t g1x;
    user_precision_complex_t D1x;
    user_precision_complex_t D1y;
    user_precision_complex_t g1y;
    user_precision_complex_t g2x;
    user_precision_complex_t D2x;
    user_precision_complex_t D2y;
    user_precision_complex_t g2y;

    if (use_twobeams == 1){
      get_beam_gains_multibeams_cpu(iBaseline, iComponent, num_freqs,
               num_baselines, num_components, num_times, beamtype,
               gxs, Dxs,
               Dys, gys,
               ant1_to_baseline_map, ant2_to_baseline_map,
               &g1x, &D1x, &D1y, &g1y, &g2x, &D2x, &D2y, &g2y);
    }
    else {
      get_beam_gains_cpu(iBaseline, iComponent, num_freqs,
               num_baselines, num_components, num_times, beamtype,
               gxs, Dxs,
               Dys, gys,
               &g1x, &D1x, &D1y, &g1y, &g2x, &D2x, &D2y, &g2y);
    }

    // if (iBaseline == 0){
    //   printf("Beam gains: %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f\n", g1x.x, g1x.y, D1x.x, D1x.y, D1y.x, D1y.y, g1y.x, g1y.y);
    //   printf("Beam gains: %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f\n", g2x.x, g2x.y, D2x.x, D2x.y, D2y.x, D2y.y, g2y.x, g2y.y);
    // }

    user_precision_complex_t visi_XX;
    user_precision_complex_t visi_XY;
    user_precision_complex_t visi_YX;
    user_precision_complex_t visi_YY;

    // printf("iComponent IQUV: %d %f %f %f %f\n", iComponent, flux_I, flux_Q, flux_U, flux_V);

    if (off_cardinal_dipoles == 1) {
      apply_beam_gains_stokesI_off_cardinal_cpu(g1x, D1x, D1y, g1y, g2x, D2x, D2y, g2y,
                    flux_I, visi_component, &visi_XX, &visi_XY, &visi_YX, &visi_YY);
    } else {
      apply_beam_gains_stokesI_on_cardinal_cpu(g1x, D1x, D1y, g1y, g2x, D2x, D2y, g2y,
                    flux_I, visi_component, &visi_XX, &visi_XY, &visi_YX, &visi_YY);
    }

    // if (iBaseline == 0){
    //   printf("Visibilities: %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f\n", visi_XX.x, visi_XX.y, visi_XY.x, visi_XY.y,
    //     visi_YX.x, visi_YX.y, visi_YY.x, visi_YY.y);
    // }

    sum_visi_XX_real[iBaseline] += creal(visi_XX);
    sum_visi_XX_imag[iBaseline] += cimag(visi_XX);

    sum_visi_XY_real[iBaseline] += creal(visi_XY);
    sum_visi_XY_imag[iBaseline] += cimag(visi_XY);

    sum_visi_YX_real[iBaseline] += creal(visi_YX);
    sum_visi_YX_imag[iBaseline] += cimag(visi_YX);

    sum_visi_YY_real[iBaseline] += creal(visi_YY);
    sum_visi_YY_imag[iBaseline] += cimag(visi_YY);

    // if (iBaseline == 0){
    //   printf("Visibilities: %.3e %.3e %.3e %.3e\n", visi_XX.x, visi_XX.y,
    //                         sum_visi_XX_real[iBaseline], sum_visi_XX_imag[iBaseline]);
    // }
}

void update_sum_visis_stokesIQUV_cpu(int iBaseline, int iComponent,
    int num_freqs, int num_baselines, int num_components, int num_times,
    int beamtype, int off_cardinal_dipoles,
    user_precision_complex_t *gxs, user_precision_complex_t *Dxs,
    user_precision_complex_t *Dys, user_precision_complex_t *gys,
    int *ant1_to_baseline_map, int *ant2_to_baseline_map, int use_twobeams,
    user_precision_complex_t visi_component,
    user_precision_t flux_I, user_precision_t flux_Q,
    user_precision_t flux_U, user_precision_t flux_V,
    user_precision_t *sum_visi_XX_real, user_precision_t *sum_visi_XX_imag,
    user_precision_t *sum_visi_XY_real, user_precision_t *sum_visi_XY_imag,
    user_precision_t *sum_visi_YX_real, user_precision_t *sum_visi_YX_imag,
    user_precision_t *sum_visi_YY_real, user_precision_t *sum_visi_YY_imag){

    user_precision_complex_t g1x;
    user_precision_complex_t D1x;
    user_precision_complex_t D1y;
    user_precision_complex_t g1y;
    user_precision_complex_t g2x;
    user_precision_complex_t D2x;
    user_precision_complex_t D2y;
    user_precision_complex_t g2y;

    if (use_twobeams == 1){
      get_beam_gains_multibeams_cpu(iBaseline, iComponent, num_freqs,
               num_baselines, num_components, num_times, beamtype,
               gxs, Dxs,
               Dys, gys,
               ant1_to_baseline_map, ant2_to_baseline_map,
               &g1x, &D1x, &D1y, &g1y, &g2x, &D2x, &D2y, &g2y);
    }
    else {
      get_beam_gains_cpu(iBaseline, iComponent, num_freqs,
               num_baselines, num_components, num_times, beamtype,
               gxs, Dxs,
               Dys, gys,
               &g1x, &D1x, &D1y, &g1y, &g2x, &D2x, &D2y, &g2y);
    }

    // if (iBaseline == 0){
    //   printf("Beam gains: %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f\n", g1x.x, g1x.y, D1x.x, D1x.y, D1y.x, D1y.y, g1y.x, g1y.y);
    //   printf("Beam gains: %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f\n", g2x.x, g2x.y, D2x.x, D2x.y, D2y.x, D2y.y, g2y.x, g2y.y);
    // }

    user_precision_complex_t visi_XX;
    user_precision_complex_t visi_XY;
    user_precision_complex_t visi_YX;
    user_precision_complex_t visi_YY;

    // printf("iComponent IQUV: %d %f %f %f %f\n", iComponent, flux_I, flux_Q, flux_U, flux_V);

    if (off_cardinal_dipoles == 1) {
      apply_beam_gains_stokesIQUV_off_cardinal_cpu(g1x, D1x, D1y, g1y, g2x, D2x, D2y, g2y,
                    flux_I, flux_Q, flux_U, flux_V,
                    visi_component, &visi_XX, &visi_XY, &visi_YX, &visi_YY);
    } else {
      apply_beam_gains_stokesIQUV_on_cardinal_cpu(g1x, D1x, D1y, g1y, g2x, D2x, D2y, g2y,
                    flux_I, flux_Q, flux_U, flux_V,
                    visi_component, &visi_XX, &visi_XY, &visi_YX, &visi_YY);
    }

    // if (iBaseline == 0){
    //   printf("Visibilities: %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f\n", visi_XX.x, visi_XX.y, visi_XY.x, visi_XY.y,
    //     visi_YX.x, visi_YX.y, visi_YY.x, visi_YY.y);
    // }

    sum_visi_XX_real[iBaseline] += creal(visi_XX);
    sum_visi_XX_imag[iBaseline] += cimag(visi_XX);

    sum_visi_XY_real[iBaseline] += creal(visi_XY);
    sum_visi_XY_imag[iBaseline] += cimag(visi_XY);

    sum_visi_YX_real[iBaseline] += creal(visi_YX);
    sum_visi_YX_imag[iBaseline] += cimag(visi_YX);

    sum_visi_YY_real[iBaseline] += creal(visi_YY);
    sum_visi_YY_imag[iBaseline] += cimag(visi_YY);

    // if (iBaseline == 0){
    //   printf("Visibilities: %.3e %.3e %.3e %.3e\n", visi_XX.x, visi_XX.y,
    //                         sum_visi_XX_real[iBaseline], sum_visi_XX_imag[iBaseline]);
    // }
}

void calc_autos_cpu(components_t components, beam_gains_t component_beam_gains,
                    int beamtype, int num_components, int num_baselines,
                    int num_freqs, int num_times, int num_ants,
                    user_precision_t *sum_visi_XX_real, user_precision_t *sum_visi_XX_imag,
                    user_precision_t *sum_visi_XY_real, user_precision_t *sum_visi_XY_imag,
                    user_precision_t *sum_visi_YX_real, user_precision_t *sum_visi_YX_imag,
                    user_precision_t *sum_visi_YY_real, user_precision_t *sum_visi_YY_imag,
                    int use_twobeams, int *ant1_to_auto_map, int *ant2_to_auto_map,
                    int off_cardinal_dipoles){

  for (int iAnt = 0; iAnt < num_ants; iAnt++) {
    for (int iTimeFreq = 0; iTimeFreq < num_times*num_freqs; iTimeFreq++) {

    int time_ind = (int)floorf( (float)iTimeFreq / (float)num_freqs);
    int freq_ind = iTimeFreq - time_ind*num_freqs;

    //Set up iBaseline to be a cross-pol of the correct time
    //and frequency step, that also correpsonds to the correct antenna
    //get_beam_gains_gpu and get_beam_gains_multibeams_gpu will use this to access the
    //correct beam gains.
    int iBaseline = num_baselines*num_freqs*time_ind + num_baselines*freq_ind + iAnt;

    int num_visis = num_baselines*num_freqs*num_times;
    int iAuto = num_visis + num_ants*num_freqs*time_ind + num_ants*freq_ind + iAnt;

    user_precision_complex_t auto_XX, auto_XY, auto_YX, auto_YY;
    user_precision_complex_t g1x, D1x, D1y, g1y, g2x, D2x, D2y, g2y;

      for (int iComponent = 0; iComponent < num_components; iComponent++) {

        if (use_twobeams == 1){
          get_beam_gains_multibeams_cpu(iBaseline, iComponent, num_freqs,
                  num_baselines, num_components, num_times, beamtype,
                  component_beam_gains.gxs, component_beam_gains.Dxs,
                  component_beam_gains.Dys, component_beam_gains.gys,
                  ant1_to_auto_map, ant2_to_auto_map,
                  &g1x, &D1x, &D1y, &g1y, &g2x, &D2x, &D2y, &g2y);
        }
        else {
          // printf("WE IS ONLY DOING THIS TING\n");
          get_beam_gains_cpu(iBaseline, iComponent, num_freqs,
                  num_baselines, num_components, num_times, beamtype,
                  component_beam_gains.gxs, component_beam_gains.Dxs,
                  component_beam_gains.Dys, component_beam_gains.gys,
                  &g1x, &D1x, &D1y, &g1y, &g2x, &D2x, &D2y, &g2y);
        }

        user_precision_complex_t visi_component = 1 + I*0.0;

        int extrap_ind = num_freqs*iComponent + freq_ind;

        user_precision_t flux_I = components.extrap_stokesI[extrap_ind];

        if (components.do_QUV == 1) {
          user_precision_t flux_Q = components.extrap_stokesQ[extrap_ind];
          user_precision_t flux_U = components.extrap_stokesU[extrap_ind];
          user_precision_t flux_V = components.extrap_stokesV[extrap_ind];

          if (off_cardinal_dipoles == 1) {
            apply_beam_gains_stokesIQUV_off_cardinal_cpu(g1x, D1x, D1y, g1y,
                                      g2x, D2x, D2y, g2y,
                                      flux_I, flux_Q, flux_U, flux_V,
                                      visi_component,
                                      &auto_XX, &auto_XY, &auto_YX, &auto_YY);
          } else {
            apply_beam_gains_stokesIQUV_on_cardinal_cpu(g1x, D1x, D1y, g1y,
                                      g2x, D2x, D2y, g2y,
                                      flux_I, flux_Q, flux_U, flux_V,
                                      visi_component,
                                      &auto_XX, &auto_XY, &auto_YX, &auto_YY);
          }
        } else {
          if (off_cardinal_dipoles == 1) {
            apply_beam_gains_stokesI_off_cardinal_cpu(g1x, D1x, D1y, g1y,
                                  g2x, D2x, D2y, g2y,
                                  flux_I, visi_component,
                                  &auto_XX, &auto_XY, &auto_YX, &auto_YY);
          } else {
            apply_beam_gains_stokesI_on_cardinal_cpu(g1x, D1x, D1y, g1y,
                                  g2x, D2x, D2y, g2y,
                                  flux_I, visi_component,
                                  &auto_XX, &auto_XY, &auto_YX, &auto_YY);
          }
        }
        
        sum_visi_XX_real[iAuto] += creal(auto_XX);
        sum_visi_XX_imag[iAuto] += cimag(auto_XX);

        sum_visi_XY_real[iAuto] += creal(auto_XY);
        sum_visi_XY_imag[iAuto] += cimag(auto_XY);

        sum_visi_YX_real[iAuto] += creal(auto_YX);
        sum_visi_YX_imag[iAuto] += cimag(auto_YX);

        sum_visi_YY_real[iAuto] += creal(auto_YY);
        sum_visi_YY_imag[iAuto] += cimag(auto_YY);
      }
    }
  }
}


void malloc_beam_gains_cpu(beam_gains_t *component_beam_gains,
                           int beamtype, int num_gains){

  //Only some models have leakage terms
  if (beamtype == FEE_BEAM || beamtype == MWA_ANALY || beamtype == FEE_BEAM_INTERP ||
      beamtype == EB_LOFAR || beamtype == EB_OSKAR || beamtype == EB_MWA || 
      beamtype == UVB_MWA || beamtype == UVB_HERA) {
    component_beam_gains->Dxs = malloc(num_gains*sizeof(user_precision_complex_t));
    component_beam_gains->Dys = malloc(num_gains*sizeof(user_precision_complex_t));
  }
  component_beam_gains->gxs = malloc(num_gains*sizeof(user_precision_complex_t));
  component_beam_gains->gys = malloc(num_gains*sizeof(user_precision_complex_t));
}

void calc_visi_point_or_gauss_cpu(components_t components,
                                  beam_gains_t component_beam_gains,
                                  calc_visi_inouts_t *calc_visi_inouts,
                                  visibility_set_t *visibility_set, 
                                  int num_components, e_beamtype beamtype,
                                  e_component_type comptype,
                                  woden_settings_t *woden_settings) {

  int num_freqs = woden_settings->num_freqs;
  int num_cross = woden_settings->num_cross;
  int num_baselines = woden_settings->num_baselines;
  int num_ants = woden_settings->num_ants;
  int num_times = woden_settings->num_time_steps;
  int off_cardinal_dipoles = woden_settings->off_cardinal_dipoles;

  int use_twobeams = component_beam_gains.use_twobeams;

  user_precision_t flux_I;
  user_precision_t flux_Q;
  user_precision_t flux_U;
  user_precision_t flux_V;

  user_precision_complex_t visi_comp;
  user_precision_complex_t V_envelop;

  user_precision_t pa, sinpa, cospa, u, v, x, y, invsig_x, invsig_y;

  int time_ind;
  int freq_ind;
  int extrap_ind;

  // int num_threads = 8;
  // omp_set_num_threads(num_threads);
  // #pragma omp parallel for collapse(2)
  // #pragma omp parallel for
  for (int iBaseline = 0; iBaseline < num_cross; iBaseline++) {
    time_ind = (int)floorf( (float)iBaseline / ((float)num_baselines * (float)num_freqs));
    freq_ind = (int)floorf( ((float)iBaseline - ((float)time_ind*(float)num_baselines * (float)num_freqs)) / (float)num_baselines);
    int baseline_ind = iBaseline % num_baselines;
    int ant1 = time_ind*num_ants + calc_visi_inouts->ant1_to_baseline_map[baseline_ind];
    int ant2 = time_ind*num_ants + calc_visi_inouts->ant2_to_baseline_map[baseline_ind];

    for (int iComponent = 0; iComponent < num_components; iComponent++) {
    
      extrap_ind = num_freqs*iComponent + freq_ind;
      // printf("INSIDE KERN components.extrap_stokesI %p\n", components.extrap_stokesI);

      flux_I = components.extrap_stokesI[extrap_ind];
      if (components.do_QUV == 1) {
        flux_Q = components.extrap_stokesQ[extrap_ind];
        flux_U = components.extrap_stokesU[extrap_ind];
        flux_V = components.extrap_stokesV[extrap_ind];
      }

      double offset = calc_ionospheric_phase_offset_cpu(calc_visi_inouts->ant_X[ant1],
                                                  calc_visi_inouts->ant_Y[ant1],
                                                  calc_visi_inouts->ant_Z[ant1],
                                                  calc_visi_inouts->ant_X[ant2],
                                                  calc_visi_inouts->ant_Y[ant2],
                                                  calc_visi_inouts->ant_Z[ant2],
                                                  components.azs[time_ind*num_components + iComponent],
                                                  components.zas[time_ind*num_components + iComponent]);
      
      visi_comp = calc_measurement_equation_cpu(calc_visi_inouts->us[iBaseline],
                                                calc_visi_inouts->vs[iBaseline],
                                                calc_visi_inouts->ws[iBaseline],
                                                components.ls[iComponent],
                                                components.ms[iComponent],
                                                components.ns[iComponent],
                                                offset);

      if (comptype == GAUSSIAN) {

        V_envelop = 1.0 + I*0.0;

        pa = components.pas[iComponent];
        sinpa = sin(pa);
        cospa = cos(pa);
        u = calc_visi_inouts->us[iBaseline];
        v = calc_visi_inouts->vs[iBaseline];

        x =  cospa*v + sinpa*u; // major axis
        y = -sinpa*v + cospa*u; // minor axis
        invsig_x = components.majors[iComponent];
        invsig_y = components.minors[iComponent];

        V_envelop = exp( -0.5 * ( x*x*invsig_x*invsig_x*M_PI_2_2_LN_2 + y*y*invsig_y*invsig_y*M_PI_2_2_LN_2 )) + I*0.0;
        visi_comp = visi_comp*V_envelop;

      }

      if (components.do_QUV == 1)
      {
        update_sum_visis_stokesIQUV_cpu(iBaseline, iComponent, num_freqs,
             num_baselines, num_components, num_times, beamtype, off_cardinal_dipoles,
             component_beam_gains.gxs, component_beam_gains.Dxs,
             component_beam_gains.Dys, component_beam_gains.gys,
             component_beam_gains.ant1_to_baseline_map,
             component_beam_gains.ant2_to_baseline_map, use_twobeams,
             visi_comp, flux_I, flux_Q, flux_U, flux_V,
             visibility_set->sum_visi_XX_real, visibility_set->sum_visi_XX_imag,
             visibility_set->sum_visi_XY_real, visibility_set->sum_visi_XY_imag,
             visibility_set->sum_visi_YX_real, visibility_set->sum_visi_YX_imag,
             visibility_set->sum_visi_YY_real, visibility_set->sum_visi_YY_imag);
      } else {
        update_sum_visis_stokesI_cpu(iBaseline, iComponent, num_freqs,
             num_baselines, num_components, num_times, beamtype, off_cardinal_dipoles,
             component_beam_gains.gxs, component_beam_gains.Dxs,
             component_beam_gains.Dys, component_beam_gains.gys,
             component_beam_gains.ant1_to_baseline_map,
             component_beam_gains.ant2_to_baseline_map, use_twobeams,
             visi_comp, flux_I,
             visibility_set->sum_visi_XX_real, visibility_set->sum_visi_XX_imag,
             visibility_set->sum_visi_XY_real, visibility_set->sum_visi_XY_imag,
             visibility_set->sum_visi_YX_real, visibility_set->sum_visi_YX_imag,
             visibility_set->sum_visi_YY_real, visibility_set->sum_visi_YY_imag);
      }
    }
  }
}


void calc_visi_shapelets_cpu(components_t components,
                             beam_gains_t component_beam_gains,
                             calc_visi_inouts_t *calc_visi_inouts,
                             visibility_set_t *visibility_set,
                             int num_shapes, int num_shape_coeffs,
                             e_beamtype beamtype,
                             woden_settings_t *woden_settings) {
  int num_freqs = woden_settings->num_freqs;
  int num_cross = woden_settings->num_cross;
  int num_baselines = woden_settings->num_baselines;
  int num_ants = woden_settings->num_ants;
  int num_times = woden_settings->num_time_steps;
  int off_cardinal_dipoles = woden_settings->off_cardinal_dipoles;

  for (int iBaseline = 0; iBaseline < num_cross; iBaseline++) {
    int use_twobeams = component_beam_gains.use_twobeams;

    user_precision_t shape_flux_I;
    user_precision_t shape_flux_Q;
    user_precision_t shape_flux_U;
    user_precision_t shape_flux_V;
    user_precision_complex_t visi_shape;

    int mobaseline = iBaseline - num_baselines*floorf((float)iBaseline / (float)num_baselines);

    //Find out what time and freq index this baseline corresponds to
    int time_ind = (int)floorf( (float)iBaseline / ((float)num_baselines * (float)num_freqs));
    int freq_ind = (int)floorf( ((float)iBaseline - ((float)time_ind*(float)num_baselines * (float)num_freqs)) / (float)num_baselines);
    int baseline_ind = iBaseline % num_baselines;
    int ant1 = time_ind*num_ants + calc_visi_inouts->ant1_to_baseline_map[baseline_ind];
    int ant2 = time_ind*num_ants + calc_visi_inouts->ant2_to_baseline_map[baseline_ind];

    for (int iCoeff = 0; iCoeff < num_shape_coeffs; iCoeff++) {

      //We have multiple coefficients per SHAPELET component - reference
      //them via this array. We chunk over coeffs so might have any
      //number of components here
      int iComponent = components.param_indexes[iCoeff];
      int extrap_ind = num_freqs*iComponent + freq_ind;

      // if (iBaseline == 0) {
      //   printf("iComponent %d iCoeff %d extrap_ind %d\n", iComponent, iCoeff, extrap_ind);
      // }
      // printf("components.extrap_stokesI %p\n", components.extrap_stokesI);

      shape_flux_I = components.extrap_stokesI[extrap_ind];

      if (components.do_QUV == 1) {
        shape_flux_Q = components.extrap_stokesQ[extrap_ind];
        shape_flux_U = components.extrap_stokesU[extrap_ind];
        shape_flux_V = components.extrap_stokesV[extrap_ind];
      }

      double offset = calc_ionospheric_phase_offset_cpu(calc_visi_inouts->ant_X[ant1],
                                                  calc_visi_inouts->ant_Y[ant1],
                                                  calc_visi_inouts->ant_Z[ant1],
                                                  calc_visi_inouts->ant_X[ant2],
                                                  calc_visi_inouts->ant_Y[ant2],
                                                  calc_visi_inouts->ant_Z[ant2],
                                                  components.azs[time_ind*num_shapes + iComponent],
                                                  components.zas[time_ind*num_shapes + iComponent]);

      visi_shape = calc_measurement_equation_cpu(calc_visi_inouts->us[iBaseline],
                                                calc_visi_inouts->vs[iBaseline],
                                                calc_visi_inouts->ws[iBaseline],
                                                components.ls[iComponent],
                                                components.ms[iComponent],
                                                components.ns[iComponent],
                                                offset);

      user_precision_t pa = components.pas[iComponent];
      user_precision_t sinpa = sin(pa);
      user_precision_t cospa = cos(pa);

      int uv_stripe = num_baselines*num_times*iComponent + time_ind*num_baselines + mobaseline;

      user_precision_t u_shape = calc_visi_inouts->u_shapes[uv_stripe] / calc_visi_inouts->allsteps_wavelengths[iBaseline];
      user_precision_t v_shape = calc_visi_inouts->v_shapes[uv_stripe] / calc_visi_inouts->allsteps_wavelengths[iBaseline];

      user_precision_t x = (cospa*v_shape + sinpa*u_shape); // major axis
      user_precision_t y = (-sinpa*v_shape + cospa*u_shape); // minor axis

      //Scales the FWHM to std to match basis functions, and account for the
      //basis functions being stored with beta = 1.0
      //Basis functions have been stored in such a way that x is in the same
      //direction as on sky, but y is opposite, so include negative here
      user_precision_t const_x = (components.majors[iComponent]*SQRT_M_PI_2_2_LN_2)/sbf_dx;
      user_precision_t const_y = -(components.minors[iComponent]*SQRT_M_PI_2_2_LN_2)/sbf_dx;

      // I^(n1+n2) = Ipow_lookup[(n1+n2) % 4]
      user_precision_complex_t Ipow_lookup[] = { 1.0 + I*0.0,
                                                 0.0 + I*1.0,
                                                -1.0 + I*0.0,
                                                 0.0 + I*-1.0 };

      user_precision_t xlow, xhigh, ylow, yhigh, u_value, v_value, f_hat, *sbf_n;

      // find the indices in the basis functions for u*beta_u and v*beta_v

      user_precision_t xpos = x*const_x + sbf_c;
      user_precision_t ypos = y*const_y + sbf_c;

      int xindex = (int)floor(xpos);
      int yindex = (int)floor(ypos);
      //
      int n1 = (int)components.n1s[iCoeff];
      int n2 = (int)components.n2s[iCoeff];

      f_hat = components.shape_coeffs[iCoeff];

      sbf_n = &calc_visi_inouts->sbf[n1*sbf_L];
      xlow  = sbf_n[xindex];
      xhigh = sbf_n[xindex+1];
      u_value = xlow + (xhigh-xlow)*(xpos-xindex);

      sbf_n = &calc_visi_inouts->sbf[n2*sbf_L];
      ylow  = sbf_n[yindex];
      yhigh = sbf_n[yindex+1];
      v_value = ylow + (yhigh-ylow)*(ypos-yindex);

      // accumulate the intensity model for baseline pair (u,v)
      user_precision_complex_t V_envelop = 0.0 + I*0.0;
      V_envelop = V_envelop + Ipow_lookup[(n1+n2) % 4] * f_hat * u_value*v_value;

      visi_shape = visi_shape*V_envelop;

      if (components.do_QUV == 1) {
        update_sum_visis_stokesIQUV_cpu(iBaseline, iComponent, num_freqs,
             num_baselines, num_shapes, num_times, beamtype, off_cardinal_dipoles,
             component_beam_gains.gxs, component_beam_gains.Dxs,
             component_beam_gains.Dys, component_beam_gains.gys,
             component_beam_gains.ant1_to_baseline_map,
             component_beam_gains.ant2_to_baseline_map, use_twobeams,
             visi_shape, shape_flux_I, shape_flux_Q, shape_flux_U, shape_flux_V,
             visibility_set->sum_visi_XX_real, visibility_set->sum_visi_XX_imag,
             visibility_set->sum_visi_XY_real, visibility_set->sum_visi_XY_imag,
             visibility_set->sum_visi_YX_real, visibility_set->sum_visi_YX_imag,
             visibility_set->sum_visi_YY_real, visibility_set->sum_visi_YY_imag);
      } else {
        update_sum_visis_stokesI_cpu(iBaseline, iComponent, num_freqs,
             num_baselines, num_shapes, num_times, beamtype, off_cardinal_dipoles,
             component_beam_gains.gxs, component_beam_gains.Dxs,
             component_beam_gains.Dys, component_beam_gains.gys,
             component_beam_gains.ant1_to_baseline_map,
             component_beam_gains.ant2_to_baseline_map, use_twobeams,
             visi_shape, shape_flux_I,
             visibility_set->sum_visi_XX_real, visibility_set->sum_visi_XX_imag,
             visibility_set->sum_visi_XY_real, visibility_set->sum_visi_XY_imag,
             visibility_set->sum_visi_YX_real, visibility_set->sum_visi_YX_imag,
             visibility_set->sum_visi_YY_real, visibility_set->sum_visi_YY_imag);
      }
    }
  }
}