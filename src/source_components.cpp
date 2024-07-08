#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include "gpucomplex.h"
#include "fundamental_coords.h"
#include "constants.h"
#include "source_components.h"
#include "woden_struct_defs.h"
#include "primary_beam_gpu.h"
#include "woden_precision_defs.h"
#include "gpu_macros.h"

__device__  gpuUserComplex calc_measurement_equation(user_precision_t *d_us,
           user_precision_t *d_vs, user_precision_t *d_ws,
           double *d_ls, double *d_ms, double *d_ns,
           const int iBaseline, const int iComponent){

  gpuUserComplex visi;

  double u, v, w;
  double l, m, n;

  u = (double)d_us[iBaseline];
  v = (double)d_vs[iBaseline];
  w = (double)d_ws[iBaseline];

  l = d_ls[iComponent];
  m = d_ms[iComponent];
  n = d_ns[iComponent];

  //Not sure why, but get match with OSKAR/RTS sims, and correct location
  //on sky through WSClean, without negative infront on 2pi
  double temp = 2*M_PI*( u*l + v*m + w*(n-1) );

  visi.y = (user_precision_t)sin(temp);
  visi.x = (user_precision_t)cos(temp);

  return visi;
}

__device__ void apply_beam_gains_stokesIQUV(gpuUserComplex g1x, gpuUserComplex D1x,
          gpuUserComplex D1y, gpuUserComplex g1y,
          gpuUserComplex g2x, gpuUserComplex D2x,
          gpuUserComplex D2y, gpuUserComplex g2y,
          user_precision_t flux_I, user_precision_t flux_Q,
          user_precision_t flux_U, user_precision_t flux_V,
          gpuUserComplex visi_component,
          gpuUserComplex * visi_XX, gpuUserComplex * visi_XY,
          gpuUserComplex * visi_YX, gpuUserComplex * visi_YY) {

  //Conjugate the second beam gains
  gpuUserComplex g2x_conj = make_gpuUserComplex(g2x.x,-g2x.y);
  gpuUserComplex D2x_conj = make_gpuUserComplex(D2x.x,-D2x.y);
  gpuUserComplex D2y_conj = make_gpuUserComplex(D2y.x,-D2y.y);
  gpuUserComplex g2y_conj = make_gpuUserComplex(g2y.x,-g2y.y);

  //Create the Stokes visibilities
  gpuUserComplex visi_I = make_gpuUserComplex(flux_I, 0.0)*visi_component;
  gpuUserComplex visi_Q = make_gpuUserComplex(flux_Q, 0.0)*visi_component;
  gpuUserComplex visi_U = make_gpuUserComplex(flux_U, 0.0)*visi_component;
  gpuUserComplex visi_V = make_gpuUserComplex(flux_V, 0.0)*visi_component;

  gpuUserComplex this_XX;
  gpuUserComplex this_XY;
  gpuUserComplex this_YX;
  gpuUserComplex this_YY;

  this_XX = (g1x*g2x_conj + D1x*D2x_conj)*visi_I;
  this_XX += (g1x*g2x_conj - D1x*D2x_conj)*visi_Q;
  this_XX += (g1x*D2x_conj + D1x*g2x_conj)*visi_U;
  this_XX += (make_gpuUserComplex(0.0,1.0)*visi_V)*(g1x*D2x_conj - D1x*g2x_conj);

  this_XY = (g1x*D2y_conj + D1x*g2y_conj)*visi_I;
  this_XY += (g1x*D2y_conj - D1x*g2y_conj)*visi_Q;
  this_XY += (g1x*g2y_conj + D1x*D2y_conj)*visi_U;
  this_XY += (make_gpuUserComplex(0.0,1.0)*visi_V)* (g1x*g2y_conj - D1x*D2y_conj);

  this_YX = (D1y*g2x_conj + g1y*D2x_conj)*visi_I;
  this_YX += (D1y*g2x_conj - g1y*D2x_conj)*visi_Q;
  this_YX += (D1y*D2x_conj + g1y*g2x_conj)*visi_U;
  this_YX += (make_gpuUserComplex(0.0,1.0)*visi_V)* (D1y*D2x_conj - g1y*g2x_conj);

  this_YY = (D1y*D2y_conj + g1y*g2y_conj)*visi_I;
  this_YY += (D1y*D2y_conj - g1y*g2y_conj)*visi_Q;
  this_YY += (D1y*g2y_conj + g1y*D2y_conj)*visi_U;
  this_YY += (make_gpuUserComplex(0.0,1.0)*visi_V)* (D1y*g2y_conj - g1y*D2y_conj);

  * visi_XX = this_XX;
  * visi_XY = this_XY;
  * visi_YX = this_YX;
  * visi_YY = this_YY;

}

__device__ void get_beam_gains(int iBaseline, int iComponent, int num_freqs,
           int num_baselines, int num_components, int num_times, int beamtype,
           gpuUserComplex *d_gxs, gpuUserComplex *d_Dxs,
           gpuUserComplex *d_Dys, gpuUserComplex *d_gys,
           gpuUserComplex * g1x, gpuUserComplex * D1x,
           gpuUserComplex * D1y, gpuUserComplex * g1y,
           gpuUserComplex * g2x, gpuUserComplex * D2x,
           gpuUserComplex * D2y, gpuUserComplex * g2y){

  int beam_ind = 0;
  int time_ind = 0;
  int freq_ind = 0;

  time_ind = (int)floorf( (float)iBaseline / ((float)num_baselines * (float)num_freqs));
  freq_ind = (int)floorf( ((float)iBaseline - ((float)time_ind*(float)num_baselines * (float)num_freqs)) / (float)num_baselines);
  beam_ind = num_freqs*time_ind*num_components + (num_components*freq_ind) + iComponent;

    //Set gains to one if no beam
  if (beamtype == NO_BEAM) {
    * g1x = make_gpuUserComplex(1.0, 0.0);
    * g2x = make_gpuUserComplex(1.0, 0.0);
    * g1y = make_gpuUserComplex(1.0, 0.0);
    * g2y = make_gpuUserComplex(1.0, 0.0);
  }

  //Get gains if using a beam
  else {
    * g1x = d_gxs[beam_ind];
    * g2x = d_gxs[beam_ind];
    * g1y = d_gys[beam_ind];
    * g2y = d_gys[beam_ind];

  }

  //Only MWA models have leakge terms at the moment
  if (beamtype == FEE_BEAM || beamtype == FEE_BEAM_INTERP || beamtype == MWA_ANALY) {
    * D1x = d_Dxs[beam_ind];
    * D2x = d_Dxs[beam_ind];
    * D1y = d_Dys[beam_ind];
    * D2y = d_Dys[beam_ind];
  }
  // Set leakage to zero if no leakage
  else {
    * D1x = make_gpuUserComplex(0.0, 0.0);
    * D2x = make_gpuUserComplex(0.0, 0.0);
    * D1y = make_gpuUserComplex(0.0, 0.0);
    * D2y = make_gpuUserComplex(0.0, 0.0);
  }
} //end __device__ get_beam_gains


__device__ void get_beam_gains_multibeams(int iBaseline, int iComponent, int num_freqs,
           int num_baselines, int num_components, int num_times, int beamtype,
           gpuUserComplex *d_gxs, gpuUserComplex *d_Dxs,
           gpuUserComplex *d_Dys, gpuUserComplex *d_gys,
           int *d_ant1_to_baseline_map, int *d_ant2_to_baseline_map,
           gpuUserComplex * g1x, gpuUserComplex * D1x,
           gpuUserComplex * D1y, gpuUserComplex * g1y,
           gpuUserComplex * g2x, gpuUserComplex * D2x,
           gpuUserComplex * D2y, gpuUserComplex * g2y){

  int time_ind = 0;
  int freq_ind = 0;

  time_ind = (int)floorf( (float)iBaseline / ((float)num_baselines * (float)num_freqs));
  freq_ind = (int)floorf( ((float)iBaseline - ((float)time_ind*(float)num_baselines * (float)num_freqs)) / (float)num_baselines);
  // beam_ind = num_freqs*time_ind*num_components + (num_components*freq_ind) + iComponent;

  int baseline_ind = iBaseline % num_baselines;

  int ant1_ind = d_ant1_to_baseline_map[baseline_ind];
  int ant2_ind = d_ant2_to_baseline_map[baseline_ind];

  int beam1 = ant1_ind*num_freqs*num_components*num_times + num_freqs*time_ind*num_components + (num_components*freq_ind) + iComponent;

  int beam2 = ant2_ind*num_freqs*num_components*num_times + num_freqs*time_ind*num_components + (num_components*freq_ind) + iComponent;

    //Set gains to one if no beam
  if (beamtype == NO_BEAM) {
    * g1x = make_gpuUserComplex(1.0, 0.0);
    * g2x = make_gpuUserComplex(1.0, 0.0);
    * g1y = make_gpuUserComplex(1.0, 0.0);
    * g2y = make_gpuUserComplex(1.0, 0.0);
  }

  //Get gains if using a beam
  else {
    * g1x = d_gxs[beam1];
    * g2x = d_gxs[beam2];
    * g1y = d_gys[beam1];
    * g2y = d_gys[beam2];

  }

  //Only MWA models have leakge terms at the moment
  if (beamtype == FEE_BEAM || beamtype == FEE_BEAM_INTERP || beamtype == MWA_ANALY) {
    * D1x = d_Dxs[beam1];
    * D2x = d_Dxs[beam2];
    * D1y = d_Dys[beam1];
    * D2y = d_Dys[beam2];
  }
  // Set leakage to zero if no leakage
  else {
    * D1x = make_gpuUserComplex(0.0, 0.0);
    * D2x = make_gpuUserComplex(0.0, 0.0);
    * D1y = make_gpuUserComplex(0.0, 0.0);
    * D2y = make_gpuUserComplex(0.0, 0.0);
  }
} //end __device__ get_beam_gains_multibeams

__device__ void apply_beam_gains_stokesI(gpuUserComplex g1x, gpuUserComplex D1x,
          gpuUserComplex D1y, gpuUserComplex g1y,
          gpuUserComplex g2x, gpuUserComplex D2x,
          gpuUserComplex D2y, gpuUserComplex g2y,
          user_precision_t flux_I,
          gpuUserComplex visi_component,
          gpuUserComplex * visi_XX, gpuUserComplex * visi_XY,
          gpuUserComplex * visi_YX, gpuUserComplex * visi_YY) {

  //Conjugate the second beam gains
  gpuUserComplex g2x_conj = make_gpuUserComplex(g2x.x,-g2x.y);
  gpuUserComplex D2x_conj = make_gpuUserComplex(D2x.x,-D2x.y);
  gpuUserComplex D2y_conj = make_gpuUserComplex(D2y.x,-D2y.y);
  gpuUserComplex g2y_conj = make_gpuUserComplex(g2y.x,-g2y.y);

  //Create the Stokes visibilities
  gpuUserComplex visi_I = make_gpuUserComplex(flux_I, 0.0)*visi_component;

  gpuUserComplex this_XX;
  gpuUserComplex this_XY;
  gpuUserComplex this_YX;
  gpuUserComplex this_YY;

  this_XX = (g1x*g2x_conj + D1x*D2x_conj)*visi_I;
  this_XY = (g1x*D2y_conj + D1x*g2y_conj)*visi_I;
  this_YX = (D1y*g2x_conj + g1y*D2x_conj)*visi_I;
  this_YY = (D1y*D2y_conj + g1y*g2y_conj)*visi_I;

  * visi_XX = this_XX;
  * visi_XY = this_XY;
  * visi_YX = this_YX;
  * visi_YY = this_YY;

}

__device__ void update_sum_visis_stokesIQUV(int iBaseline, int iComponent, int num_freqs,
    int num_baselines, int num_components, int num_times, int beamtype,
    gpuUserComplex *d_gxs, gpuUserComplex *d_Dxs,
    gpuUserComplex *d_Dys, gpuUserComplex *d_gys,
    int *d_ant1_to_baseline_map, int *d_ant2_to_baseline_map, int use_twobeams,
    gpuUserComplex visi_component,
    user_precision_t flux_I, user_precision_t flux_Q,
    user_precision_t flux_U, user_precision_t flux_V,
    user_precision_t *d_sum_visi_XX_real, user_precision_t *d_sum_visi_XX_imag,
    user_precision_t *d_sum_visi_XY_real, user_precision_t *d_sum_visi_XY_imag,
    user_precision_t *d_sum_visi_YX_real, user_precision_t *d_sum_visi_YX_imag,
    user_precision_t *d_sum_visi_YY_real, user_precision_t *d_sum_visi_YY_imag){

    gpuUserComplex g1x;
    gpuUserComplex D1x;
    gpuUserComplex D1y;
    gpuUserComplex g1y;
    gpuUserComplex g2x;
    gpuUserComplex D2x;
    gpuUserComplex D2y;
    gpuUserComplex g2y;

    if (use_twobeams == 1){
      get_beam_gains_multibeams(iBaseline, iComponent, num_freqs,
               num_baselines, num_components, num_times, beamtype,
               d_gxs, d_Dxs,
               d_Dys, d_gys,
               d_ant1_to_baseline_map, d_ant2_to_baseline_map,
               &g1x, &D1x, &D1y, &g1y, &g2x, &D2x, &D2y, &g2y);
    }
    else {
      get_beam_gains(iBaseline, iComponent, num_freqs,
               num_baselines, num_components, num_times, beamtype,
               d_gxs, d_Dxs,
               d_Dys, d_gys,
               &g1x, &D1x, &D1y, &g1y, &g2x, &D2x, &D2y, &g2y);
    }

    gpuUserComplex visi_XX;
    gpuUserComplex visi_XY;
    gpuUserComplex visi_YX;
    gpuUserComplex visi_YY;

    apply_beam_gains_stokesIQUV(g1x, D1x, D1y, g1y, g2x, D2x, D2y, g2y,
                    flux_I, flux_Q, flux_U, flux_V,
                    visi_component, &visi_XX, &visi_XY, &visi_YX, &visi_YY);

    d_sum_visi_XX_real[iBaseline] += visi_XX.x;
    d_sum_visi_XX_imag[iBaseline] += visi_XX.y;

    d_sum_visi_XY_real[iBaseline] += visi_XY.x;
    d_sum_visi_XY_imag[iBaseline] += visi_XY.y;

    d_sum_visi_YX_real[iBaseline] += visi_YX.x;
    d_sum_visi_YX_imag[iBaseline] += visi_YX.y;

    d_sum_visi_YY_real[iBaseline] += visi_YY.x;
    d_sum_visi_YY_imag[iBaseline] += visi_YY.y;
}

__device__ void update_sum_visis_stokesI(int iBaseline, int iComponent, int num_freqs,
    int num_baselines, int num_components, int num_times, int beamtype,
    gpuUserComplex *d_gxs, gpuUserComplex *d_Dxs,
    gpuUserComplex *d_Dys, gpuUserComplex *d_gys,
    int *d_ant1_to_baseline_map, int *d_ant2_to_baseline_map, int use_twobeams,
    gpuUserComplex visi_component,
    user_precision_t flux_I,
    user_precision_t *d_sum_visi_XX_real, user_precision_t *d_sum_visi_XX_imag,
    user_precision_t *d_sum_visi_XY_real, user_precision_t *d_sum_visi_XY_imag,
    user_precision_t *d_sum_visi_YX_real, user_precision_t *d_sum_visi_YX_imag,
    user_precision_t *d_sum_visi_YY_real, user_precision_t *d_sum_visi_YY_imag){

    gpuUserComplex g1x;
    gpuUserComplex D1x;
    gpuUserComplex D1y;
    gpuUserComplex g1y;
    gpuUserComplex g2x;
    gpuUserComplex D2x;
    gpuUserComplex D2y;
    gpuUserComplex g2y;

    if (use_twobeams == 1){
      get_beam_gains_multibeams(iBaseline, iComponent, num_freqs,
               num_baselines, num_components, num_times, beamtype,
               d_gxs, d_Dxs,
               d_Dys, d_gys,
               d_ant1_to_baseline_map, d_ant2_to_baseline_map,
               &g1x, &D1x, &D1y, &g1y, &g2x, &D2x, &D2y, &g2y);
    }
    else {
      get_beam_gains(iBaseline, iComponent, num_freqs,
               num_baselines, num_components, num_times, beamtype,
               d_gxs, d_Dxs,
               d_Dys, d_gys,
               &g1x, &D1x, &D1y, &g1y, &g2x, &D2x, &D2y, &g2y);
    }

    gpuUserComplex visi_XX;
    gpuUserComplex visi_XY;
    gpuUserComplex visi_YX;
    gpuUserComplex visi_YY;

    apply_beam_gains_stokesI(g1x, D1x, D1y, g1y, g2x, D2x, D2y, g2y,
                    flux_I,
                    visi_component, &visi_XX, &visi_XY, &visi_YX, &visi_YY);

    d_sum_visi_XX_real[iBaseline] += visi_XX.x;
    d_sum_visi_XX_imag[iBaseline] += visi_XX.y;

    d_sum_visi_XY_real[iBaseline] += visi_XY.x;
    d_sum_visi_XY_imag[iBaseline] += visi_XY.y;

    d_sum_visi_YX_real[iBaseline] += visi_YX.x;
    d_sum_visi_YX_imag[iBaseline] += visi_YX.y;

    d_sum_visi_YY_real[iBaseline] += visi_YY.x;
    d_sum_visi_YY_imag[iBaseline] += visi_YY.y;

}


//Allocate space for the extrapolated Stokes parameters
void malloc_extrapolated_flux_arrays(components_t *d_components, int num_comps,
                                     int num_freqs, int do_QUV){
  d_components->extrap_stokesI = NULL;
  ( gpuMalloc( (void**)&d_components->extrap_stokesI,
                                   num_comps*num_freqs*sizeof(double) ));

  if (do_QUV == 1)
  {
      d_components->extrap_stokesQ = NULL;
      ( gpuMalloc( (void**)&d_components->extrap_stokesQ,
                                      num_comps*num_freqs*sizeof(double) ));
      d_components->extrap_stokesU = NULL;
      ( gpuMalloc( (void**)&d_components->extrap_stokesU,
                                      num_comps*num_freqs*sizeof(double) ));
      d_components->extrap_stokesV = NULL;
      ( gpuMalloc( (void**)&d_components->extrap_stokesV,
                                      num_comps*num_freqs*sizeof(double) ));

  }
}

__device__ void extrap_stokes_power_law_stokesI(components_t d_components,
           double *d_extrap_freqs, int iFluxComp, int iFreq,
           user_precision_t * flux_I){

  double d_freq = d_extrap_freqs[iFreq];
  double d_ref_freq = d_components.power_ref_freqs[iFluxComp];

  user_precision_t flux_ratio = pow(d_freq / d_ref_freq, d_components.power_SIs[iFluxComp]);

  * flux_I = d_components.power_ref_stokesI[iFluxComp] * flux_ratio;
}

__global__ void kern_extrap_power_laws_stokesI(int num_extrap_freqs, double *d_extrap_freqs,
                                       int num_comps, components_t d_components) {

  // Start by computing which baseline we're going to do
  const int iFluxComp = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iFreq = threadIdx.y + (blockDim.y*blockIdx.y);
  if(iFluxComp < num_comps && iFreq < num_extrap_freqs) {

    user_precision_t flux_I;

    extrap_stokes_power_law_stokesI(d_components, d_extrap_freqs,
                 iFluxComp, iFreq, &flux_I);

    int iComponent = d_components.power_comp_inds[iFluxComp];
    int extrap_ind = num_extrap_freqs*iComponent + iFreq;

    d_components.extrap_stokesI[extrap_ind] = flux_I;
  }
}

__device__ void extrap_stokes_curved_power_law_stokesI(components_t d_components,
           double *d_extrap_freqs, int iFluxComp, int iFreq,
           user_precision_t * flux_I){

  double d_freq = d_extrap_freqs[iFreq];
  double d_ref_freq = d_components.curve_ref_freqs[iFluxComp];

  user_precision_t si_ratio = pow(d_freq / d_ref_freq, d_components.curve_SIs[iFluxComp]);

  double log_freq_ratio = log(d_freq / d_ref_freq);

  double q = (double)d_components.curve_qs[iFluxComp];
  double exp_bit = exp(q*log_freq_ratio*log_freq_ratio);

  user_precision_t flux_ratio = si_ratio * exp_bit;

  * flux_I = d_components.curve_ref_stokesI[iFluxComp] * flux_ratio;
}

__global__ void kern_extrap_curved_power_laws_stokesI(int num_extrap_freqs, double *d_extrap_freqs,
                                              int num_comps, components_t d_components) {

  // Start by computing which baseline we're going to do
  const int iFluxComp = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iFreq = threadIdx.y + (blockDim.y*blockIdx.y);

  if(iFluxComp < num_comps && iFreq < num_extrap_freqs) {

    user_precision_t flux_I;

    extrap_stokes_curved_power_law_stokesI(d_components, d_extrap_freqs,
                 iFluxComp, iFreq,
                 &flux_I);

    int iComponent = d_components.curve_comp_inds[iFluxComp];
    int extrap_ind = num_extrap_freqs*iComponent + iFreq;

    d_components.extrap_stokesI[extrap_ind] = flux_I;

  }
}

__device__ user_precision_t calc_gradient_extrap_list(user_precision_t *list_fluxes,
          double *list_freqs, double desired_freq, int low_ind_1, int low_ind_2) {

  user_precision_t gradient;
  user_precision_t extrap_flux;

  //If both zero, just stick to zero
  if (list_fluxes[low_ind_1] == 0 && list_fluxes[low_ind_2] == 0) {
   extrap_flux = 0.0;
  }

  //If one is negative, do interpolation in linear space
  else if (list_fluxes[low_ind_1] <= 0 || list_fluxes[low_ind_2] <= 0) {
    gradient = (list_fluxes[low_ind_2] - list_fluxes[low_ind_1]) / (list_freqs[low_ind_2] - list_freqs[low_ind_1]);
    extrap_flux = list_fluxes[low_ind_1] + gradient*(desired_freq - list_freqs[low_ind_1]);
  }

  else {

    user_precision_t logflux1, logflux2, logfreq1, logfreq2, log_des_freq;

    logflux1 = log10(list_fluxes[low_ind_1]);
    logflux2 = log10(list_fluxes[low_ind_2]);
    logfreq1 = log10(list_freqs[low_ind_1]);
    logfreq2 = log10(list_freqs[low_ind_2]);
    log_des_freq = log10(desired_freq);

    gradient = (logflux2 - logflux1) / (logfreq2 - logfreq1);
    extrap_flux = logflux1 + gradient*(log_des_freq - logfreq1);

    extrap_flux = pow(10, extrap_flux);

  }
  return extrap_flux;
}


__device__ void extrap_stokes_list_flux_stokesIQUV(components_t d_components,
           double *d_extrap_freqs, int iFluxComp, int iFreq,
           user_precision_t * flux_I, user_precision_t * flux_Q,
           user_precision_t * flux_U, user_precision_t * flux_V){

  int num_list_values = d_components.num_list_values[iFluxComp];
  int list_start_ind = d_components.list_start_indexes[iFluxComp];

  double d_extrap_freq = d_extrap_freqs[iFreq];

  int low_ind_1 = -1;
  int low_ind_2 = -1;

  double low_val_1 = 1e16;
  // double low_val_2 = 1e16;

  double ref_freq;
  double abs_diff_freq;

  if (num_list_values == 1) {
    * flux_I = d_components.list_stokesI[list_start_ind];
    * flux_Q = d_components.list_stokesQ[list_start_ind];
    * flux_U = d_components.list_stokesU[list_start_ind];
    * flux_V = d_components.list_stokesV[list_start_ind];
    return;
  }

  //First loop finds the absolute closest frequency
  for (int i = 0; i < num_list_values; i++) {
    ref_freq = d_components.list_freqs[list_start_ind + i];
    abs_diff_freq = abs(ref_freq - d_extrap_freq);

    if (abs_diff_freq < low_val_1) {
      low_val_1 = abs_diff_freq;
      low_ind_1 = i;
    }
  }

  //Depending on the closest frequency, we either want to search above or
  //below the target frequency to find points either side of the target freq

  //We happen to need the reference frequency; just return the refs
  if (d_components.list_freqs[list_start_ind + low_ind_1] == d_extrap_freq) {
    * flux_I = d_components.list_stokesI[list_start_ind + low_ind_1];
    * flux_Q = d_components.list_stokesQ[list_start_ind + low_ind_1];
    * flux_U = d_components.list_stokesU[list_start_ind + low_ind_1];
    * flux_V = d_components.list_stokesV[list_start_ind + low_ind_1];
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
      if (d_components.list_freqs[list_start_ind + low_ind_1] > d_extrap_freq){
        low_ind_2 = low_ind_1 - 1;
      }
      else {
        low_ind_2 = low_ind_1 + 1;
      }
        //We are extrapolating to a frequency that is lower than all list entries
        //so just stick low_ind_2 to one above low_ind_1
    }
  }

  * flux_I = calc_gradient_extrap_list(d_components.list_stokesI,
            d_components.list_freqs, d_extrap_freq,
            list_start_ind + low_ind_1, list_start_ind + low_ind_2);
  * flux_Q = calc_gradient_extrap_list(d_components.list_stokesQ,
            d_components.list_freqs, d_extrap_freq,
            list_start_ind + low_ind_1, list_start_ind + low_ind_2);
  * flux_U = calc_gradient_extrap_list(d_components.list_stokesU,
            d_components.list_freqs, d_extrap_freq,
            list_start_ind + low_ind_1, list_start_ind + low_ind_2);
  * flux_V = calc_gradient_extrap_list(d_components.list_stokesV,
            d_components.list_freqs, d_extrap_freq,
            list_start_ind + low_ind_1, list_start_ind + low_ind_2);

  if (low_ind_2 == -1){

    printf("wrong range %.3e %.3e iFreq %d %.3e low %d %.3e\n", d_components.list_freqs[list_start_ind],
    d_components.list_freqs[list_start_ind + num_list_values-1],
    iFreq, d_extrap_freq,
    low_ind_1, d_components.list_freqs[list_start_ind + low_ind_1]);
    printf("The flooxes %.3e %.3e %.3e %.3e\n",* flux_I, * flux_Q, * flux_U, * flux_V );
  }
}

__device__ void extrap_stokes_list_flux_stokesI(components_t d_components,
           double *d_extrap_freqs, int iFluxComp, int iFreq,
           user_precision_t * flux_I){

  int num_list_values = d_components.num_list_values[iFluxComp];
  int list_start_ind = d_components.list_start_indexes[iFluxComp];

  double d_extrap_freq = d_extrap_freqs[iFreq];

  int low_ind_1 = -1;
  int low_ind_2 = -1;

  double low_val_1 = 1e16;
  // double low_val_2 = 1e16;

  double ref_freq;
  double abs_diff_freq;

  if (num_list_values == 1) {
    * flux_I = d_components.list_stokesI[list_start_ind];
    return;
  }

  //First loop finds the absolute closest frequency
  for (int i = 0; i < num_list_values; i++) {
    ref_freq = d_components.list_freqs[list_start_ind + i];
    abs_diff_freq = abs(ref_freq - d_extrap_freq);

    if (abs_diff_freq < low_val_1) {
      low_val_1 = abs_diff_freq;
      low_ind_1 = i;
    }
  }

  //Depending on the closest frequency, we either want to search above or
  //below the target frequency to find points either side of the target freq

  //We happen to need the reference frequency; just return the refs
  if (d_components.list_freqs[list_start_ind + low_ind_1] == d_extrap_freq) {
    * flux_I = d_components.list_stokesI[list_start_ind + low_ind_1];
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
      if (d_components.list_freqs[list_start_ind + low_ind_1] > d_extrap_freq){
        low_ind_2 = low_ind_1 - 1;
      }
      else {
        low_ind_2 = low_ind_1 + 1;
      }
        //We are extrapolating to a frequency that is lower than all list entries
        //so just stick low_ind_2 to one above low_ind_1
    }
  }

  * flux_I = calc_gradient_extrap_list(d_components.list_stokesI,
            d_components.list_freqs, d_extrap_freq,
            list_start_ind + low_ind_1, list_start_ind + low_ind_2);

  if (low_ind_2 == -1){

    printf("wrong range %.3e %.3e iFreq %d %.3e low %d %.3e\n", d_components.list_freqs[list_start_ind],
    d_components.list_freqs[list_start_ind + num_list_values-1],
    iFreq, d_extrap_freq,
    low_ind_1, d_components.list_freqs[list_start_ind + low_ind_1]);
    printf("The flooxes %.3e \n",* flux_I);
  }
}

__global__ void kern_extrap_list_fluxes_stokesIQUV(int num_extrap_freqs, double *d_extrap_freqs,
                                        int num_comps, components_t d_components) {

  // Start by computing which baseline we're going to do
  const int iFluxComp = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iFreq = threadIdx.y + (blockDim.y*blockIdx.y);
  if(iFluxComp < num_comps && iFreq < num_extrap_freqs) {

    user_precision_t flux_I;
    user_precision_t flux_Q;
    user_precision_t flux_U;
    user_precision_t flux_V;

    extrap_stokes_list_flux_stokesIQUV(d_components, d_extrap_freqs,
                 iFluxComp, iFreq,
                 &flux_I, &flux_Q, &flux_U, &flux_V);

    int iComponent = d_components.list_comp_inds[iFluxComp];
    int extrap_ind = num_extrap_freqs*iComponent + iFreq;

    d_components.extrap_stokesI[extrap_ind] = flux_I;
    d_components.extrap_stokesQ[extrap_ind] = flux_Q;
    d_components.extrap_stokesU[extrap_ind] = flux_U;
    d_components.extrap_stokesV[extrap_ind] = flux_V;

  }
}

__global__ void kern_extrap_list_fluxes_stokesI(int num_extrap_freqs, double *d_extrap_freqs,
                                        int num_comps, components_t d_components) {

  // Start by computing which baseline we're going to do
  const int iFluxComp = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iFreq = threadIdx.y + (blockDim.y*blockIdx.y);
  if(iFluxComp < num_comps && iFreq < num_extrap_freqs) {

    user_precision_t flux_I;

    extrap_stokes_list_flux_stokesI(d_components, d_extrap_freqs,
                 iFluxComp, iFreq,
                 &flux_I);

    int iComponent = d_components.list_comp_inds[iFluxComp];
    int extrap_ind = num_extrap_freqs*iComponent + iFreq;

    d_components.extrap_stokesI[extrap_ind] = flux_I;

  }
}


extern "C" void extrapolate_Stokes(source_t *d_chunked_source,
                                   double *d_extrap_freqs, int num_extrap_freqs,
                                   e_component_type comptype,
                                   int do_QUV){

  components_t d_components;
  // int n_comps = 0;
  int n_powers = 0;
  int n_curves = 0;
  int n_lists = 0;

  //Choose the right components to extrapolate for
  if (comptype == POINT) {
    d_components = d_chunked_source->point_components;
    // n_comps = d_chunked_source->n_points;
    n_powers = d_chunked_source->n_point_powers;
    n_curves = d_chunked_source->n_point_curves;
    n_lists = d_chunked_source->n_point_lists;
  }
  else if (comptype == GAUSSIAN) {
    d_components = d_chunked_source->gauss_components;
    // n_comps = d_chunked_source->n_gauss;
    n_powers = d_chunked_source->n_gauss_powers;
    n_curves = d_chunked_source->n_gauss_curves;
    n_lists = d_chunked_source->n_gauss_lists;
  // } else if (comptype == SHAPELET) {
  } else {
    d_components = d_chunked_source->shape_components;
    // n_comps = d_chunked_source->n_shapes;
    n_powers = d_chunked_source->n_shape_powers;
    n_curves = d_chunked_source->n_shape_curves;
    n_lists = d_chunked_source->n_shape_lists;
  }

  dim3 grid, threads;

  threads.x = 16;
  threads.y = 16;

  //First up, do the POWER_LAW types
  grid.y = (int)ceilf( (float)num_extrap_freqs / (float)threads.y );

  if (n_powers > 0) {
    grid.x = (int)ceilf( (float)n_powers / (float)threads.x );
    gpuErrorCheckKernel("kern_extrap_power_laws_stokesI",
                          kern_extrap_power_laws_stokesI, grid, threads,
                          num_extrap_freqs, d_extrap_freqs,
                          n_powers, d_components);
  }
  //Next up, do the CURVED_POWER_LAW types
  if (n_curves > 0) {
    grid.x = (int)ceilf( (float)n_curves / (float)threads.x );
    gpuErrorCheckKernel("kern_extrap_curved_power_laws_stokesI",
                      kern_extrap_curved_power_laws_stokesI, grid, threads,
                      num_extrap_freqs, d_extrap_freqs,
                      n_curves, d_components);
  }

  //Finally, do any list flux peeps
  if (n_lists > 0) {
    grid.x = (int)ceilf( (float)n_lists / (float)threads.x );
    if (do_QUV == 1) {
      gpuErrorCheckKernel("kern_extrap_list_fluxes_stokesIQUV",
                          kern_extrap_list_fluxes_stokesIQUV, grid, threads,
                          num_extrap_freqs, d_extrap_freqs,
                          n_lists, d_components);
    }
    else {
      gpuErrorCheckKernel("kern_extrap_list_fluxes_stokesI",
                          kern_extrap_list_fluxes_stokesI, grid, threads,
                          num_extrap_freqs, d_extrap_freqs,
                          n_lists, d_components);
    }
    
  }
}

extern "C" void source_component_common(woden_settings_t *woden_settings,
           beam_settings_t *beam_settings, double *d_freqs,
           source_t *chunked_source, source_t *d_chunked_source,
           d_beam_gains_t *d_component_beam_gains,
           e_component_type comptype,
           visibility_set_t *d_visibility_set){

  int do_QUV = woden_settings->do_QUV;

  //Here we see if we a single primary beam for all (num_beams = 1) or
  //a primary beam per antenna (num_beams = num_ants)
  //This can be expanded in the future to have a primary beam per tile
  //for different options
  int num_beams;
  int use_twobeams = 0;
    if (woden_settings->use_dipamps == 1) {
      num_beams = woden_settings->num_ants;
      use_twobeams = 1;
    } else {
      num_beams = 1;
  }

  int num_components = 0;
  components_t *components = NULL;
  components_t *d_components = NULL;

  if (comptype == POINT) {
    num_components = d_chunked_source->n_points;
    components = &chunked_source->point_components;
    d_components = &d_chunked_source->point_components;
  } else if (comptype == GAUSSIAN) {
    num_components = d_chunked_source->n_gauss;
    components = &chunked_source->gauss_components;
    d_components = &d_chunked_source->gauss_components;
  } else if (comptype == SHAPELET) {
    num_components = d_chunked_source->n_shapes;
    components = &chunked_source->shape_components;
    d_components = &d_chunked_source->shape_components;
  }

  //Will need this later
  malloc_extrapolated_flux_arrays(d_components, num_components,
                                  woden_settings->num_freqs, do_QUV);

  extrapolate_Stokes(d_chunked_source, d_freqs,
                     woden_settings->num_freqs, comptype, do_QUV);

  int num_gains = d_components->num_primarybeam_values*num_beams;
  
  //Only the MWA beams currently yields cross pol values, so only malloc what
  //we need here
  //TODO in the future, this might need to be a loop over all primary beams,
  //if we have different beams for different tiles
  if (beam_settings->beamtype == FEE_BEAM || beam_settings->beamtype == MWA_ANALY || beam_settings->beamtype == FEE_BEAM_INTERP) {
    ( gpuMalloc( (void**)&d_component_beam_gains->d_Dxs,
                    num_gains*sizeof(gpuUserComplex) ));
    ( gpuMalloc( (void**)&d_component_beam_gains->d_Dys,
                    num_gains*sizeof(gpuUserComplex) ));
  }
  ( gpuMalloc( (void**)&d_component_beam_gains->d_gxs,
                    num_gains*sizeof(gpuUserComplex) ));
  ( gpuMalloc( (void**)&d_component_beam_gains->d_gys,
                    num_gains*sizeof(gpuUserComplex) ));
  //
  ( gpuMalloc( (void**)&d_components->ls,
                                               num_components*sizeof(double) ) );
  ( gpuMalloc( (void**)&d_components->ms,
                                               num_components*sizeof(double) ) );
  ( gpuMalloc( (void**)&d_components->ns,
                                               num_components*sizeof(double) ) );


  dim3 grid, threads;

  threads.x = 128;
  threads.y = 1;
  threads.z = 1;
  grid.x = (int)ceil( (float)num_components / (float)threads.x );
  grid.y = 1;
  grid.z = 1;

  gpuErrorCheckKernel("kern_calc_lmn",
                        kern_calc_lmn, grid, threads,
                        woden_settings->ra0,
                        woden_settings->sdec0, woden_settings->cdec0,
                        d_components->ras, d_components->decs,
                        d_components->ls, d_components->ms, d_components->ns, num_components);

  //If using a gaussian primary beam, calculate beam values for all freqs,
  //lsts and point component locations
  if (beam_settings->beamtype == GAUSS_BEAM) {

    //TODO currently hardcoded to have beam position angle = 0.
    //Should this change with az/za?
    user_precision_t cos_theta = 1.0;
    user_precision_t sin_theta = 0.0;
    user_precision_t sin_2theta = 0.0;
    user_precision_t fwhm_lm = sin(beam_settings->beam_FWHM_rad);

    printf("\tDoing Gaussian Beam\n");

    calculate_gaussian_beam(num_components,
         woden_settings->num_time_steps, woden_settings->num_freqs,
         beam_settings->gauss_ha, beam_settings->gauss_sdec,
         beam_settings->gauss_cdec,
         fwhm_lm, cos_theta, sin_theta, sin_2theta,
         beam_settings->beam_ref_freq, d_freqs,
         components->beam_has,
         components->beam_decs,
         d_component_beam_gains->d_gxs, d_component_beam_gains->d_gys);

  }// end if beam == GAUSS

  else if (beam_settings->beamtype == FEE_BEAM || beam_settings->beamtype == FEE_BEAM_INTERP) {

    if (beam_settings->beamtype == FEE_BEAM_INTERP) {
      printf("\tDoing the hyperbeam (interpolated)\n");
    } else {
      printf("\tDoing the hyperbeam\n");
    }

    //Have to reorder the az/za from comp ind, time ind to time ind, comp ind
    //before feeding into mwa_hyperbeam
    int num_azza = woden_settings->num_time_steps*num_components;
    double *reordered_azs = (double *)malloc(num_azza*sizeof(double));
    double *reordered_zas = (double *)malloc(num_azza*sizeof(double));

    int stripe_new, stripe_old;

    for (int time_ind = 0; time_ind < woden_settings->num_time_steps; time_ind++) {
      for (int comp_ind = 0; comp_ind < num_components; comp_ind++) {
        stripe_new = time_ind*num_components + comp_ind;
        stripe_old = comp_ind*woden_settings->num_time_steps + time_ind;
        reordered_azs[stripe_new] = (double)components->azs[stripe_old];
        reordered_zas[stripe_new] = (double)components->zas[stripe_old];
      }
    }

    //Always be doing parallatic angle rotation
    uint8_t parallactic = 1;

    run_hyperbeam_cuda(num_components,
           woden_settings->num_time_steps, woden_settings->num_freqs,
           num_beams, parallactic,
           beam_settings->cuda_fee_beam,
           reordered_azs, reordered_zas,
           woden_settings->latitudes,
           d_component_beam_gains->d_gxs, d_component_beam_gains->d_Dxs,
           d_component_beam_gains->d_Dys, d_component_beam_gains->d_gys);

    free(reordered_azs);
    free(reordered_zas);
  }

  else if (beam_settings->beamtype == ANALY_DIPOLE) {
    printf("\tDoing analytic_dipole (EDA2 beam)\n");

    calculate_analytic_dipole_beam(num_components,
         woden_settings->num_time_steps, woden_settings->num_freqs,
         components->azs, components->zas, d_freqs,
         d_component_beam_gains->d_gxs, d_component_beam_gains->d_gys);
  }

  else if (beam_settings->beamtype == MWA_ANALY) {

    //Always normalise to zenith
    int norm = 1;

    printf("\tDoing analytic MWA beam\n");

    calculate_RTS_MWA_analytic_beam(num_components,
         woden_settings->num_time_steps, woden_settings->num_freqs,
         components->azs, components->zas,
         woden_settings->FEE_ideal_delays, woden_settings->latitude,
         norm, components->beam_has, components->beam_decs,
         d_freqs, d_component_beam_gains->d_gxs, d_component_beam_gains->d_Dxs,
         d_component_beam_gains->d_Dys, d_component_beam_gains->d_gys);
  }

  //Now we've calculated the beams, we can calculate the auto-correlations,
  //if so required

  if (woden_settings->do_autos){

    int num_freqs = woden_settings->num_freqs;
    int num_times = woden_settings->num_time_steps;
    int num_ants = woden_settings->num_ants;
    int num_baselines = woden_settings->num_baselines;

    threads.x = 64;
    threads.y = 2;
    threads.z = 1;
    grid.x = (int)ceil( (float)(num_freqs*num_times) / (float)threads.x );
    grid.y = (int)ceil( (float)(num_ants) / (float)threads.y );
    grid.z = 1;

    int *d_ant_to_auto_map = NULL;

    if (use_twobeams == 1) {
      int *ant_to_auto_map = NULL;
      ant_to_auto_map = (int *)malloc(num_ants*sizeof(int));
      for (int ant = 0; ant < num_ants; ant++){
          ant_to_auto_map[ant] = ant;
      }
      ( gpuMalloc( (void**)&d_ant_to_auto_map,
                                    num_ants*sizeof(int) ));
      ( gpuMemcpy(d_ant_to_auto_map, ant_to_auto_map,
                                      num_ants*sizeof(int), gpuMemcpyHostToDevice ));
      free(ant_to_auto_map);
    }

    gpuErrorCheckKernel("kern_calc_autos",
                  kern_calc_autos, grid, threads,
                  *d_components, *d_component_beam_gains,
                  beam_settings->beamtype,
                  num_components, num_baselines,
                  num_freqs, num_times, num_ants,
                  d_visibility_set->sum_visi_XX_real,
                  d_visibility_set->sum_visi_XX_imag,
                  d_visibility_set->sum_visi_XY_real,
                  d_visibility_set->sum_visi_XY_imag,
                  d_visibility_set->sum_visi_YX_real,
                  d_visibility_set->sum_visi_YX_imag,
                  d_visibility_set->sum_visi_YY_real,
                  d_visibility_set->sum_visi_YY_imag,
                  do_QUV, use_twobeams, d_ant_to_auto_map,
                  d_ant_to_auto_map);

    if (use_twobeams == 1) {
    (  gpuFree( d_ant_to_auto_map ) );
    }
  }

} //END source_component_common


__global__ void kern_calc_visi_point_or_gauss(components_t d_components,
           d_beam_gains_t d_component_beam_gains,
           user_precision_t *d_us, user_precision_t *d_vs, user_precision_t *d_ws,
           user_precision_t *d_sum_visi_XX_real, user_precision_t *d_sum_visi_XX_imag,
           user_precision_t *d_sum_visi_XY_real, user_precision_t *d_sum_visi_XY_imag,
           user_precision_t *d_sum_visi_YX_real, user_precision_t *d_sum_visi_YX_imag,
           user_precision_t *d_sum_visi_YY_real, user_precision_t *d_sum_visi_YY_imag,
           int num_components, int num_baselines, int num_freqs, int num_cross,
           int num_times, e_beamtype beamtype, e_component_type comptype,
           int do_QUV) {

  // Start by computing which baseline we're going to do
  const int iBaseline = threadIdx.x + (blockDim.x*blockIdx.x);
  if(iBaseline < num_cross) {

    int use_twobeams = d_component_beam_gains.use_twobeams;

    user_precision_t flux_I;
    user_precision_t flux_Q;
    user_precision_t flux_U;
    user_precision_t flux_V;

    gpuUserComplex visi_comp;
    gpuUserComplex V_envelop;

    user_precision_t pa, sinpa, cospa, u, v, x, y, invsig_x, invsig_y;

    //Find out what time and freq index this baseline corresponds to
    int time_ind = (int)floorf( (float)iBaseline / ((float)num_baselines * (float)num_freqs));
    int freq_ind = (int)floorf( ((float)iBaseline - ((float)time_ind*(float)num_baselines * (float)num_freqs)) / (float)num_baselines);

    for (int iComponent = 0; iComponent < num_components; iComponent++) {


      int extrap_ind = num_freqs*iComponent + freq_ind;

      flux_I = d_components.extrap_stokesI[extrap_ind];
      if (do_QUV == 1) {
        flux_Q = d_components.extrap_stokesQ[extrap_ind];
        flux_U = d_components.extrap_stokesU[extrap_ind];
        flux_V = d_components.extrap_stokesV[extrap_ind];
      }
      


      visi_comp = calc_measurement_equation(d_us, d_vs, d_ws,
                             d_components.ls, d_components.ms, d_components.ns,
                             iBaseline, iComponent);

      if (comptype == GAUSSIAN) {

        V_envelop = make_gpuUserComplex( 1.0, 0.0 );

        pa = d_components.pas[iComponent];
        sinpa = sin(pa);
        cospa = cos(pa);
        u = d_us[iBaseline];
        v = d_vs[iBaseline];

        x =  cospa*v + sinpa*u; // major axis
        y = -sinpa*v + cospa*u; // minor axis
        invsig_x = d_components.majors[iComponent];
        invsig_y = d_components.minors[iComponent];

        V_envelop = make_gpuUserComplex( exp( -0.5 * ( x*x*invsig_x*invsig_x*M_PI_2_2_LN_2 + y*y*invsig_y*invsig_y*M_PI_2_2_LN_2 ) ), 0.0 );

        visi_comp = visi_comp*V_envelop;
      }

      if (do_QUV == 1)
      {
        update_sum_visis_stokesIQUV(iBaseline, iComponent, num_freqs,
             num_baselines, num_components, num_times, beamtype,
             d_component_beam_gains.d_gxs, d_component_beam_gains.d_Dxs,
             d_component_beam_gains.d_Dys, d_component_beam_gains.d_gys,
             d_component_beam_gains.d_ant1_to_baseline_map,
             d_component_beam_gains.d_ant2_to_baseline_map, use_twobeams,
             visi_comp, flux_I, flux_Q, flux_U, flux_V,
             d_sum_visi_XX_real, d_sum_visi_XX_imag,
             d_sum_visi_XY_real, d_sum_visi_XY_imag,
             d_sum_visi_YX_real, d_sum_visi_YX_imag,
             d_sum_visi_YY_real, d_sum_visi_YY_imag);
      } else {
        update_sum_visis_stokesI(iBaseline, iComponent, num_freqs,
             num_baselines, num_components, num_times, beamtype,
             d_component_beam_gains.d_gxs, d_component_beam_gains.d_Dxs,
             d_component_beam_gains.d_Dys, d_component_beam_gains.d_gys,
             d_component_beam_gains.d_ant1_to_baseline_map,
             d_component_beam_gains.d_ant2_to_baseline_map, use_twobeams,
             visi_comp, flux_I,
             d_sum_visi_XX_real, d_sum_visi_XX_imag,
             d_sum_visi_XY_real, d_sum_visi_XY_imag,
             d_sum_visi_YX_real, d_sum_visi_YX_imag,
             d_sum_visi_YY_real, d_sum_visi_YY_imag);
      }
    }
  }
}

__global__ void kern_calc_visi_shapelets(components_t d_components,
      d_beam_gains_t d_component_beam_gains,
      user_precision_t *d_us, user_precision_t *d_vs, user_precision_t *d_ws,
      user_precision_t *d_allsteps_wavelengths,
      user_precision_t *d_u_shapes, user_precision_t *d_v_shapes,
      user_precision_t *d_sum_visi_XX_real, user_precision_t *d_sum_visi_XX_imag,
      user_precision_t *d_sum_visi_XY_real, user_precision_t *d_sum_visi_XY_imag,
      user_precision_t *d_sum_visi_YX_real, user_precision_t *d_sum_visi_YX_imag,
      user_precision_t *d_sum_visi_YY_real, user_precision_t *d_sum_visi_YY_imag,
      user_precision_t *d_sbf,
      int num_shapes, int num_baselines, int num_freqs, int num_cross,
      const int num_coeffs, int num_times, e_beamtype beamtype,
      int do_QUV) {

  // Start by computing which baseline we're going to do
  const int iBaseline = threadIdx.x + (blockDim.x*blockIdx.x);

  if (iBaseline < num_cross) {

    int use_twobeams = d_component_beam_gains.use_twobeams;

    user_precision_t shape_flux_I;
    user_precision_t shape_flux_Q;
    user_precision_t shape_flux_U;
    user_precision_t shape_flux_V;
    gpuUserComplex visi_shape;

    int mod_baseline = iBaseline - num_baselines*floorf((float)iBaseline / (float)num_baselines);

    //Find out what time and freq index this baseline corresponds to
    int time_ind = (int)floorf( (float)iBaseline / ((float)num_baselines * (float)num_freqs));
    int freq_ind = (int)floorf( ((float)iBaseline - ((float)time_ind*(float)num_baselines * (float)num_freqs)) / (float)num_baselines);

    for (int iCoeff = 0; iCoeff < num_coeffs; iCoeff++) {

      //We have multiple coefficients per SHAPELET component - reference
      //them via this array. We chunk over coeffs so might have any
      //number of components here
      int iComponent = d_components.param_indexes[iCoeff];
      int extrap_ind = num_freqs*iComponent + freq_ind;

      shape_flux_I = d_components.extrap_stokesI[extrap_ind];

      if (do_QUV == 1) {
        shape_flux_Q = d_components.extrap_stokesQ[extrap_ind];
        shape_flux_U = d_components.extrap_stokesU[extrap_ind];
        shape_flux_V = d_components.extrap_stokesV[extrap_ind];
      }

      visi_shape = calc_measurement_equation(d_us, d_vs, d_ws,
                            d_components.ls, d_components.ms, d_components.ns,
                            iBaseline, iComponent);

      user_precision_t pa = d_components.pas[iComponent];
      user_precision_t sinpa = sin(pa);
      user_precision_t cospa = cos(pa);

      int uv_stripe = num_baselines*num_times*iComponent + time_ind*num_baselines + mod_baseline;

      user_precision_t u_shape = d_u_shapes[uv_stripe] / d_allsteps_wavelengths[iBaseline];
      user_precision_t v_shape = d_v_shapes[uv_stripe] / d_allsteps_wavelengths[iBaseline];

      user_precision_t x = (cospa*v_shape + sinpa*u_shape); // major axis
      user_precision_t y = (-sinpa*v_shape + cospa*u_shape); // minor axis

      //Scales the FWHM to std to match basis functions, and account for the
      //basis functions being stored with beta = 1.0
      //Basis functions have been stored in such a way that x is in the same
      //direction as on sky, but y is opposite, so include negative here
      user_precision_t const_x = (d_components.majors[iComponent]*SQRT_M_PI_2_2_LN_2)/sbf_dx;
      user_precision_t const_y = -(d_components.minors[iComponent]*SQRT_M_PI_2_2_LN_2)/sbf_dx;

      // I^(n1+n2) = Ipow_lookup[(n1+n2) % 4]
      gpuUserComplex Ipow_lookup[] = { make_gpuUserComplex(  1.0,  0.0 ),
                                       make_gpuUserComplex(  0.0,  1.0 ),
                                       make_gpuUserComplex( -1.0,  0.0 ),
                                       make_gpuUserComplex(  0.0, -1.0 ) };

      user_precision_t xlow, xhigh, ylow, yhigh, u_value, v_value, f_hat, *sbf_n;

      // find the indices in the basis functions for u*beta_u and v*beta_v

      user_precision_t xpos = x*const_x + sbf_c;
      user_precision_t ypos = y*const_y + sbf_c;

      int xindex = (int)floor(xpos);
      int yindex = (int)floor(ypos);
      //
      int n1 = (int)d_components.n1s[iCoeff];
      int n2 = (int)d_components.n2s[iCoeff];

      f_hat = d_components.shape_coeffs[iCoeff];

      sbf_n = &d_sbf[n1*sbf_L];
      xlow  = sbf_n[xindex];
      xhigh = sbf_n[xindex+1];
      u_value = xlow + (xhigh-xlow)*(xpos-xindex);

      sbf_n = &d_sbf[n2*sbf_L];
      ylow  = sbf_n[yindex];
      yhigh = sbf_n[yindex+1];
      v_value = ylow + (yhigh-ylow)*(ypos-yindex);

      // accumulate the intensity model for baseline pair (u,v)
      gpuUserComplex V_envelop = make_gpuUserComplex( 0.0, 0.0 );
      V_envelop = V_envelop + Ipow_lookup[(n1+n2) % 4] * f_hat * u_value*v_value;

      visi_shape = visi_shape*V_envelop;

      if (do_QUV == 1) {
        update_sum_visis_stokesIQUV(iBaseline, iComponent, num_freqs,
             num_baselines, num_shapes, num_times, beamtype,
             d_component_beam_gains.d_gxs, d_component_beam_gains.d_Dxs,
             d_component_beam_gains.d_Dys, d_component_beam_gains.d_gys,
             d_component_beam_gains.d_ant1_to_baseline_map,
             d_component_beam_gains.d_ant2_to_baseline_map, use_twobeams,
             visi_shape,
             shape_flux_I, shape_flux_Q, shape_flux_U, shape_flux_V,
             d_sum_visi_XX_real, d_sum_visi_XX_imag,
             d_sum_visi_XY_real, d_sum_visi_XY_imag,
             d_sum_visi_YX_real, d_sum_visi_YX_imag,
             d_sum_visi_YY_real, d_sum_visi_YY_imag);
      } else {
        update_sum_visis_stokesI(iBaseline, iComponent, num_freqs,
             num_baselines, num_shapes, num_times, beamtype,
             d_component_beam_gains.d_gxs, d_component_beam_gains.d_Dxs,
             d_component_beam_gains.d_Dys, d_component_beam_gains.d_gys,
             d_component_beam_gains.d_ant1_to_baseline_map,
             d_component_beam_gains.d_ant2_to_baseline_map, use_twobeams,
             visi_shape, shape_flux_I,
             d_sum_visi_XX_real, d_sum_visi_XX_imag,
             d_sum_visi_XY_real, d_sum_visi_XY_imag,
             d_sum_visi_YX_real, d_sum_visi_YX_imag,
             d_sum_visi_YY_real, d_sum_visi_YY_imag);
      }
    }
  }
}

//Copy the sky model info from a set of components from the CPU to the GPU
void copy_components_to_GPU(source_t *chunked_source, source_t *d_chunked_source,
                            e_component_type comptype) {

  components_t *components=NULL;
  components_t *d_components=NULL;
  int num_comps = 0, num_shape_coeffs = 0;
  int num_powers = 0, num_curves = 0, num_lists = 0;

  if (comptype == POINT) {
    components = &chunked_source->point_components;
    d_components = &d_chunked_source->point_components;

    num_comps = chunked_source->n_points;
    num_shape_coeffs = 0;
    num_powers = chunked_source->n_point_powers;
    num_curves = chunked_source->n_point_curves;
    num_lists = chunked_source->n_point_lists;

  }
  else if (comptype == GAUSSIAN) {
    components = &chunked_source->gauss_components;
    d_components = &d_chunked_source->gauss_components;

    num_comps = chunked_source->n_gauss;
    num_shape_coeffs = 0;
    num_powers = chunked_source->n_gauss_powers;
    num_curves = chunked_source->n_gauss_curves;
    num_lists = chunked_source->n_gauss_lists;

  }
  // else if (comptype == SHAPELET) {
  else {
    components = &chunked_source->shape_components;
    d_components = &d_chunked_source->shape_components;

    num_comps = chunked_source->n_shapes;
    num_shape_coeffs = chunked_source->n_shape_coeffs;
    num_powers = chunked_source->n_shape_powers;
    num_curves = chunked_source->n_shape_curves;
    num_lists = chunked_source->n_shape_lists;

  }

  //Common attributes between all flux types and components types
  ( gpuMalloc( (void**)&d_components->ras,
                      num_comps*sizeof(double) ) );
  ( gpuMemcpy( d_components->ras, components->ras,
                      num_comps*sizeof(double), gpuMemcpyHostToDevice ) );

  ( gpuMalloc( (void**)&d_components->decs,
                      num_comps*sizeof(double) ) );
  ( gpuMemcpy( d_components->decs, components->decs,
                      num_comps*sizeof(double), gpuMemcpyHostToDevice ) );

  d_components->num_primarybeam_values = components->num_primarybeam_values;

  //GAUSSIAN and SHAPELET only attributes
  if (comptype == GAUSSIAN || comptype == SHAPELET ) {
    ( gpuMalloc( (void**)&d_components->pas,
                        num_comps*sizeof(user_precision_t) ) );
    ( gpuMemcpy( d_components->pas, components->pas,
                        num_comps*sizeof(user_precision_t), gpuMemcpyHostToDevice ) );

    ( gpuMalloc( (void**)&d_components->majors,
                        num_comps*sizeof(user_precision_t) ) );
    ( gpuMemcpy( d_components->majors, components->majors,
                        num_comps*sizeof(user_precision_t), gpuMemcpyHostToDevice ) );

    ( gpuMalloc( (void**)&d_components->minors,
                        num_comps*sizeof(user_precision_t) ) );
    ( gpuMemcpy( d_components->minors, components->minors,
                        num_comps*sizeof(user_precision_t), gpuMemcpyHostToDevice ) );
  }

  //SHAPELET only attributes
  if (comptype == SHAPELET) {
    ( gpuMalloc( (void**)&d_components->shape_coeffs,
                        num_shape_coeffs*sizeof(user_precision_t) ) );
    ( gpuMemcpy( d_components->shape_coeffs, components->shape_coeffs,
                        num_shape_coeffs*sizeof(user_precision_t),
                        gpuMemcpyHostToDevice ) );

    ( gpuMalloc( (void**)&d_components->n1s,
                        num_shape_coeffs*sizeof(user_precision_t) ) );
    ( gpuMemcpy( d_components->n1s, components->n1s,
                        num_shape_coeffs*sizeof(user_precision_t),
                        gpuMemcpyHostToDevice ) );

    ( gpuMalloc( (void**)&d_components->n2s,
                        num_shape_coeffs*sizeof(user_precision_t) ) );
    ( gpuMemcpy( d_components->n2s, components->n2s,
                        num_shape_coeffs*sizeof(user_precision_t),
                        gpuMemcpyHostToDevice ) );

    ( gpuMalloc( (void**)&d_components->param_indexes,
                        num_shape_coeffs*sizeof(user_precision_t) ) );
    ( gpuMemcpy( d_components->param_indexes, components->param_indexes,
                        num_shape_coeffs*sizeof(user_precision_t),
                        gpuMemcpyHostToDevice ) );
  }

  //POWER_LAW flux things
  if (num_powers > 0) {
    ( gpuMalloc( (void**)&d_components->power_comp_inds,
                        num_powers*sizeof(int) ) );
    ( gpuMemcpy( d_components->power_comp_inds, components->power_comp_inds,
                        num_powers*sizeof(int), gpuMemcpyHostToDevice ) );

    ( gpuMalloc( (void**)&d_components->power_ref_freqs,
                        num_powers*sizeof(double) ) );
    ( gpuMemcpy( d_components->power_ref_freqs, components->power_ref_freqs,
                        num_powers*sizeof(double), gpuMemcpyHostToDevice ) );

    ( gpuMalloc( (void**)&d_components->power_ref_stokesI,
                        num_powers*sizeof(user_precision_t) ) );
    ( gpuMemcpy( d_components->power_ref_stokesI, components->power_ref_stokesI,
                        num_powers*sizeof(user_precision_t), gpuMemcpyHostToDevice ) );

    ( gpuMalloc( (void**)&d_components->power_ref_stokesQ,
                        num_powers*sizeof(user_precision_t) ) );
    ( gpuMemcpy( d_components->power_ref_stokesQ, components->power_ref_stokesQ,
                        num_powers*sizeof(user_precision_t), gpuMemcpyHostToDevice ) );

    ( gpuMalloc( (void**)&d_components->power_ref_stokesU,
                        num_powers*sizeof(user_precision_t) ) );
    ( gpuMemcpy( d_components->power_ref_stokesU, components->power_ref_stokesU,
                        num_powers*sizeof(user_precision_t), gpuMemcpyHostToDevice ) );

    ( gpuMalloc( (void**)&d_components->power_ref_stokesV,
                        num_powers*sizeof(user_precision_t) ) );
    ( gpuMemcpy( d_components->power_ref_stokesV, components->power_ref_stokesV,
                        num_powers*sizeof(user_precision_t), gpuMemcpyHostToDevice ) );

    ( gpuMalloc( (void**)&d_components->power_SIs,
                        num_powers*sizeof(user_precision_t) ) );
    ( gpuMemcpy( d_components->power_SIs, components->power_SIs,
                        num_powers*sizeof(user_precision_t), gpuMemcpyHostToDevice ) );
  }

  //CURVED_POWER_LAW things
  if (num_curves > 0) {
    ( gpuMalloc( (void**)&d_components->curve_comp_inds,
                        num_curves*sizeof(int) ) );
    ( gpuMemcpy( d_components->curve_comp_inds, components->curve_comp_inds,
                        num_curves*sizeof(int), gpuMemcpyHostToDevice ) );

    ( gpuMalloc( (void**)&d_components->curve_ref_freqs,
                        num_curves*sizeof(double) ) );
    ( gpuMemcpy( d_components->curve_ref_freqs, components->curve_ref_freqs,
                        num_curves*sizeof(double), gpuMemcpyHostToDevice ) );

    ( gpuMalloc( (void**)&d_components->curve_ref_stokesI,
                        num_curves*sizeof(user_precision_t) ) );
    ( gpuMemcpy( d_components->curve_ref_stokesI, components->curve_ref_stokesI,
                        num_curves*sizeof(user_precision_t), gpuMemcpyHostToDevice ) );

    ( gpuMalloc( (void**)&d_components->curve_ref_stokesQ,
                        num_curves*sizeof(user_precision_t) ) );
    ( gpuMemcpy( d_components->curve_ref_stokesQ, components->curve_ref_stokesQ,
                        num_curves*sizeof(user_precision_t), gpuMemcpyHostToDevice ) );

    ( gpuMalloc( (void**)&d_components->curve_ref_stokesU,
                        num_curves*sizeof(user_precision_t) ) );
    ( gpuMemcpy( d_components->curve_ref_stokesU, components->curve_ref_stokesU,
                        num_curves*sizeof(user_precision_t), gpuMemcpyHostToDevice ) );

    ( gpuMalloc( (void**)&d_components->curve_ref_stokesV,
                        num_curves*sizeof(user_precision_t) ) );
    ( gpuMemcpy( d_components->curve_ref_stokesV, components->curve_ref_stokesV,
                        num_curves*sizeof(user_precision_t), gpuMemcpyHostToDevice ) );

    ( gpuMalloc( (void**)&d_components->curve_SIs,
                        num_curves*sizeof(user_precision_t) ) );
    ( gpuMemcpy( d_components->curve_SIs, components->curve_SIs,
                        num_curves*sizeof(user_precision_t), gpuMemcpyHostToDevice ) );

    ( gpuMalloc( (void**)&d_components->curve_qs,
                        num_curves*sizeof(user_precision_t) ) );
    ( gpuMemcpy( d_components->curve_qs, components->curve_qs,
                        num_curves*sizeof(user_precision_t), gpuMemcpyHostToDevice ) );
  }

  //LIST things
  if (num_lists > 0) {
    int num_list_values = components->total_num_flux_entires;

    ( gpuMalloc( (void**)&d_components->list_comp_inds,
                        num_lists*sizeof(int) ) );
    ( gpuMemcpy( d_components->list_comp_inds,
                        components->list_comp_inds,
                        num_lists*sizeof(int), gpuMemcpyHostToDevice ) );

    ( gpuMalloc( (void**)&d_components->num_list_values,
                        num_lists*sizeof(int) ) );
    ( gpuMemcpy( d_components->num_list_values,
                        components->num_list_values,
                        num_lists*sizeof(int), gpuMemcpyHostToDevice ) );

    ( gpuMalloc( (void**)&d_components->list_start_indexes,
                        num_lists*sizeof(int) ) );
    ( gpuMemcpy( d_components->list_start_indexes,
                        components->list_start_indexes,
                        num_lists*sizeof(int), gpuMemcpyHostToDevice ) );

    ( gpuMalloc( (void**)&d_components->list_freqs,
                        num_list_values*sizeof(double) ) );
    ( gpuMemcpy( d_components->list_freqs, components->list_freqs,
                        num_list_values*sizeof(double), gpuMemcpyHostToDevice ) );

    ( gpuMalloc( (void**)&d_components->list_stokesI,
                        num_list_values*sizeof(user_precision_t) ) );
    ( gpuMemcpy( d_components->list_stokesI, components->list_stokesI,
                        num_list_values*sizeof(user_precision_t), gpuMemcpyHostToDevice ) );

    ( gpuMalloc( (void**)&d_components->list_stokesQ,
                        num_list_values*sizeof(user_precision_t) ) );
    ( gpuMemcpy( d_components->list_stokesQ, components->list_stokesQ,
                        num_list_values*sizeof(user_precision_t), gpuMemcpyHostToDevice ) );

    ( gpuMalloc( (void**)&d_components->list_stokesU,
                        num_list_values*sizeof(user_precision_t) ) );
    ( gpuMemcpy( d_components->list_stokesU, components->list_stokesU,
                        num_list_values*sizeof(user_precision_t), gpuMemcpyHostToDevice ) );

    ( gpuMalloc( (void**)&d_components->list_stokesV,
                        num_list_values*sizeof(user_precision_t) ) );
    ( gpuMemcpy( d_components->list_stokesV, components->list_stokesV,
                        num_list_values*sizeof(user_precision_t), gpuMemcpyHostToDevice ) );

  }
}

source_t * copy_chunked_source_to_GPU(source_t *chunked_source){

  source_t *d_chunked_source = (source_t*)malloc(sizeof(source_t));

  if (chunked_source->n_points > 0) {
    copy_components_to_GPU(chunked_source, d_chunked_source, POINT);
  }
  if (chunked_source->n_gauss > 0) {
    copy_components_to_GPU(chunked_source, d_chunked_source, GAUSSIAN);
  }
  if (chunked_source->n_shapes > 0) {
    copy_components_to_GPU(chunked_source, d_chunked_source, SHAPELET);
  }

  //copy across the component counters

  d_chunked_source->n_points = chunked_source->n_points;
  d_chunked_source->n_point_lists = chunked_source->n_point_lists;
  d_chunked_source->n_point_powers = chunked_source->n_point_powers;
  d_chunked_source->n_point_curves = chunked_source->n_point_curves;

  d_chunked_source->n_gauss = chunked_source->n_gauss;
  d_chunked_source->n_gauss_lists = chunked_source->n_gauss_lists;
  d_chunked_source->n_gauss_powers = chunked_source->n_gauss_powers;
  d_chunked_source->n_gauss_curves = chunked_source->n_gauss_curves;

  d_chunked_source->n_shapes = chunked_source->n_shapes;
  d_chunked_source->n_shape_lists = chunked_source->n_shape_lists;
  d_chunked_source->n_shape_powers = chunked_source->n_shape_powers;
  d_chunked_source->n_shape_curves = chunked_source->n_shape_curves;
  d_chunked_source->n_shape_coeffs = chunked_source->n_shape_coeffs;

  return d_chunked_source;
}

void free_extrapolated_flux_arrays(components_t *d_components, int do_QUV){
  ( gpuFree( d_components->extrap_stokesI ) );

  if (do_QUV) {
    ( gpuFree( d_components->extrap_stokesQ ) );
    ( gpuFree( d_components->extrap_stokesU ) );
    ( gpuFree( d_components->extrap_stokesV ) );
  }
}



extern "C" void free_d_components(source_t *d_chunked_source,
                                  e_component_type comptype){
  components_t d_components;
  int n_powers = 0;
  int n_curves = 0;
  int n_lists = 0;

  if (comptype == POINT) {
    d_components = d_chunked_source->point_components;
    n_powers = d_chunked_source->n_point_powers;
    n_curves = d_chunked_source->n_point_curves;
    n_lists = d_chunked_source->n_point_lists;
  }
  else if (comptype == GAUSSIAN) {
    d_components = d_chunked_source->gauss_components;
    n_powers = d_chunked_source->n_gauss_powers;
    n_curves = d_chunked_source->n_gauss_curves;
    n_lists = d_chunked_source->n_gauss_lists;
  }
  else {
    d_components = d_chunked_source->shape_components;
    n_powers = d_chunked_source->n_shape_powers;
    n_curves = d_chunked_source->n_shape_curves;
    n_lists = d_chunked_source->n_shape_lists;
  }

  ( gpuFree( d_components.decs) );
  ( gpuFree( d_components.ras) );

  ( gpuFree( d_components.ls) );
  ( gpuFree( d_components.ms) );
  ( gpuFree( d_components.ns) );

  //The az,za,beam_has,beam_decs are handled by other functions

  if (n_powers > 0) {
    ( gpuFree( d_components.power_ref_freqs ) );
    ( gpuFree( d_components.power_ref_stokesI ) );
    ( gpuFree( d_components.power_ref_stokesQ ) );
    ( gpuFree( d_components.power_ref_stokesU ) );
    ( gpuFree( d_components.power_ref_stokesV ) );
    ( gpuFree( d_components.power_SIs ) );
    ( gpuFree( d_components.power_comp_inds ) );
  }

  if (n_curves > 0) {
    ( gpuFree( d_components.curve_ref_freqs ) );
    ( gpuFree( d_components.curve_ref_stokesI ) );
    ( gpuFree( d_components.curve_ref_stokesQ ) );
    ( gpuFree( d_components.curve_ref_stokesU ) );
    ( gpuFree( d_components.curve_ref_stokesV ) );
    ( gpuFree( d_components.curve_SIs ) );
    ( gpuFree( d_components.curve_qs ) );
    ( gpuFree( d_components.curve_comp_inds ) );
  }
  if (n_lists > 0) {
    ( gpuFree( d_components.list_comp_inds ) );
    ( gpuFree( d_components.list_freqs ) );
    ( gpuFree( d_components.list_stokesI ) );
    ( gpuFree( d_components.list_stokesQ ) );
    ( gpuFree( d_components.list_stokesU ) );
    ( gpuFree( d_components.list_stokesV ) );
    ( gpuFree( d_components.num_list_values ) );
    ( gpuFree( d_components.list_start_indexes ) );
  }

  if (comptype == GAUSSIAN || comptype == SHAPELET) {
    ( gpuFree( d_components.pas ) );
    ( gpuFree( d_components.majors ) );
    ( gpuFree( d_components.minors ) );
  }

  if (comptype == SHAPELET) {
    ( gpuFree( d_components.shape_coeffs ) );
    ( gpuFree( d_components.n1s ) );
    ( gpuFree( d_components.n2s ) );
    ( gpuFree( d_components.param_indexes ) );
  }
}

extern "C" void free_beam_gains(d_beam_gains_t d_beam_gains, e_beamtype beamtype){

  ( gpuFree( d_beam_gains.d_gxs) );
  ( gpuFree( d_beam_gains.d_gys) );

  if (beamtype == FEE_BEAM || beamtype == FEE_BEAM_INTERP || beamtype == MWA_ANALY){
    ( gpuFree( d_beam_gains.d_Dxs ) );
    ( gpuFree( d_beam_gains.d_Dys ) );
  }

}

//Calculate auto-correlations
__global__ void kern_calc_autos(components_t d_components,
                                d_beam_gains_t d_component_beam_gains,
                                int beamtype,
                                int num_components, int num_baselines,
                                int num_freqs, int num_times, int num_ants,
                                user_precision_t *d_sum_visi_XX_real,
                                user_precision_t *d_sum_visi_XX_imag,
                                user_precision_t *d_sum_visi_XY_real,
                                user_precision_t *d_sum_visi_XY_imag,
                                user_precision_t *d_sum_visi_YX_real,
                                user_precision_t *d_sum_visi_YX_imag,
                                user_precision_t *d_sum_visi_YY_real,
                                user_precision_t *d_sum_visi_YY_imag,
                                int do_QUV, int use_twobeams,
                                int *d_ant1_to_auto_map,
                                int *d_ant2_to_auto_map) {

  // Start by computing which baseline we're going to do
  const int iTimeFreq = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iAnt = threadIdx.y + (blockDim.y*blockIdx.y);

  //TODO One day we might have different primary beams for each tile,
  //then we'll have to use iAuto to reference different primary
  //beams - at the mo we get grab the same primary beam for all antennas
  // const int iAuto = threadIdx.z + (blockDim.z*blockIdx.z);
  if(iAnt < num_ants && iTimeFreq < num_times*num_freqs) {

    int time_ind = (int)floorf( (float)iTimeFreq / (float)num_freqs);
    int freq_ind = iTimeFreq - time_ind*num_freqs;

    //Set up iBaseline to be a cross-pol of the correct time
    //and frequency step, that also correpsonds to the correct antenna
    //get_beam_gains and get_beam_gains_multibeams will use this to access the
    //correct beam gains.
    int iBaseline = num_baselines*num_freqs*time_ind + num_baselines*freq_ind + iAnt;

    int num_visis = num_baselines*num_freqs*num_times;
    int iAuto = num_visis + num_ants*num_freqs*time_ind + num_ants*freq_ind + iAnt;

    gpuUserComplex auto_XX, auto_XY, auto_YX, auto_YY;
    gpuUserComplex g1x, D1x, D1y, g1y, g2x, D2x, D2y, g2y;

    for (int iComponent = 0; iComponent < num_components; iComponent++) {

      if (use_twobeams == 1){
        get_beam_gains_multibeams(iBaseline, iComponent, num_freqs,
                num_baselines, num_components, num_times, beamtype,
                d_component_beam_gains.d_gxs, d_component_beam_gains.d_Dxs,
                d_component_beam_gains.d_Dys, d_component_beam_gains.d_gys,
                d_ant1_to_auto_map, d_ant2_to_auto_map,
                &g1x, &D1x, &D1y, &g1y, &g2x, &D2x, &D2y, &g2y);
      }
      else {
        // printf("WE IS ONLY DOING THIS TING\n");
        get_beam_gains(iBaseline, iComponent, num_freqs,
                num_baselines, num_components, num_times, beamtype,
                d_component_beam_gains.d_gxs, d_component_beam_gains.d_Dxs,
                d_component_beam_gains.d_Dys, d_component_beam_gains.d_gys,
                &g1x, &D1x, &D1y, &g1y, &g2x, &D2x, &D2y, &g2y);
      }

      gpuUserComplex visi_component;
      visi_component = make_gpuUserComplex(1.0, 0.0);

      int extrap_ind = num_freqs*iComponent + freq_ind;

      user_precision_t flux_I = d_components.extrap_stokesI[extrap_ind];

      if (do_QUV == 1) {
        user_precision_t flux_Q = d_components.extrap_stokesQ[extrap_ind];
        user_precision_t flux_U = d_components.extrap_stokesU[extrap_ind];
        user_precision_t flux_V = d_components.extrap_stokesV[extrap_ind];

        apply_beam_gains_stokesIQUV(g1x, D1x, D1y, g1y, g2x, D2x, D2y, g2y,
                                    flux_I, flux_Q, flux_U, flux_V,
                                    visi_component,
                                    &auto_XX, &auto_XY, &auto_YX, &auto_YY);

      } else {
        apply_beam_gains_stokesI(g1x, D1x, D1y, g1y, g2x, D2x, D2y, g2y,
                                 flux_I, visi_component,
                                 &auto_XX, &auto_XY, &auto_YX, &auto_YY);
      }
      
      d_sum_visi_XX_real[iAuto] += auto_XX.x;
      d_sum_visi_XX_imag[iAuto] += auto_XX.y;

      d_sum_visi_XY_real[iAuto] += auto_XY.x;
      d_sum_visi_XY_imag[iAuto] += auto_XY.y;

      d_sum_visi_YX_real[iAuto] += auto_YX.x;
      d_sum_visi_YX_imag[iAuto] += auto_YX.y;

      d_sum_visi_YY_real[iAuto] += auto_YY.x;
      d_sum_visi_YY_imag[iAuto] += auto_YY.y;

    }
  }
}

extern "C" void fill_ant_to_baseline_mapping(int num_ants, int *d_ant1_to_baseline_map,
                                               int *d_ant2_to_baseline_map){

  int num_baselines = ((num_ants - 1)*num_ants) / 2;

  int *ant1_to_baseline_map = NULL;
  int *ant2_to_baseline_map = NULL;

  ant1_to_baseline_map = (int *)malloc(num_baselines*sizeof(int));
  ant2_to_baseline_map = (int *)malloc(num_baselines*sizeof(int));

  //These functions only do cross correlations, so create all combos of antennas
  //that make up all the crosses
  int cross_index = 0;
  for (int ant1 = 0; ant1 < num_ants-1; ant1++)
  {
    for (int ant2 = ant1 + 1; ant2 < num_ants; ant2++)
    {
      ant1_to_baseline_map[cross_index] = ant1;
      ant2_to_baseline_map[cross_index] = ant2;

      cross_index += 1;
    }
  }

  ( gpuMemcpy(d_ant1_to_baseline_map, ant1_to_baseline_map,
                                  num_baselines*sizeof(int), gpuMemcpyHostToDevice ));
  ( gpuMemcpy(d_ant2_to_baseline_map, ant2_to_baseline_map,
                                  num_baselines*sizeof(int), gpuMemcpyHostToDevice ));

  free(ant1_to_baseline_map);
  free(ant2_to_baseline_map);

}

/*******************************************************************************
                 Functions below to be used in unit tests
*******************************************************************************/

extern "C" void test_extrap_stokes_all_models(source_t *chunked_source,
           int num_extrap_freqs, double *extrap_freqs,
           user_precision_t *extrap_flux_I, user_precision_t *extrap_flux_Q,
           user_precision_t *extrap_flux_U, user_precision_t *extrap_flux_V){

  int do_QUV = 1;

  source_t *d_chunked_source = copy_chunked_source_to_GPU(chunked_source);

  double *d_extrap_freqs = NULL;
  ( gpuMalloc( (void**)&d_extrap_freqs,
                                   num_extrap_freqs*sizeof(double) ));
  ( gpuMemcpy(d_extrap_freqs, extrap_freqs,
             num_extrap_freqs*sizeof(double), gpuMemcpyHostToDevice ));

  malloc_extrapolated_flux_arrays(&d_chunked_source->point_components,
                                  d_chunked_source->n_points,
                                  num_extrap_freqs, do_QUV);

  extrapolate_Stokes(d_chunked_source, d_extrap_freqs, num_extrap_freqs, POINT, do_QUV);


  components_t d_components = d_chunked_source->point_components;

  ( gpuMemcpy(extrap_flux_I, d_components.extrap_stokesI,
            d_chunked_source->n_points*num_extrap_freqs*sizeof(user_precision_t),
                                                      gpuMemcpyDeviceToHost ));
  ( gpuMemcpy(extrap_flux_Q, d_components.extrap_stokesQ,
            d_chunked_source->n_points*num_extrap_freqs*sizeof(user_precision_t),
                                                      gpuMemcpyDeviceToHost ));
  ( gpuMemcpy(extrap_flux_U, d_components.extrap_stokesU,
            d_chunked_source->n_points*num_extrap_freqs*sizeof(user_precision_t),
                                                      gpuMemcpyDeviceToHost ));
  ( gpuMemcpy(extrap_flux_V, d_components.extrap_stokesV,
            d_chunked_source->n_points*num_extrap_freqs*sizeof(user_precision_t),
                                                      gpuMemcpyDeviceToHost ));
  //
  ( gpuFree( d_extrap_freqs ) );
  free_extrapolated_flux_arrays(&d_chunked_source->point_components, do_QUV);
}


__global__ void kern_calc_measurement_equation(int num_components, int num_baselines,
          user_precision_t *d_us, user_precision_t *d_vs, user_precision_t *d_ws,
          double *d_ls, double *d_ms, double *d_ns, gpuUserComplex *d_visis) {

  // Start by computing which baseline we're going to do
  const int iBaseline = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iComponent = threadIdx.y + (blockDim.y*blockIdx.y);

  if(iComponent < num_components && iBaseline < num_baselines) {

    gpuUserComplex visi;
    visi = calc_measurement_equation(d_us, d_vs, d_ws, d_ls, d_ms, d_ns,
                                     iBaseline, iComponent);

    int visi_ind = num_components*iBaseline + iComponent;
    d_visis[visi_ind] = visi;

  }
}

extern "C" void test_kern_calc_measurement_equation(int num_components,
          int num_baselines,
          user_precision_t *us, user_precision_t *vs, user_precision_t *ws,
          double *ls, double *ms, double *ns, user_precision_complex_t *visis){

  user_precision_t *d_us = NULL;
  user_precision_t *d_vs = NULL;
  user_precision_t *d_ws = NULL;
  ( gpuMalloc( (void**)&d_us, num_baselines*sizeof(user_precision_t) ));
  ( gpuMalloc( (void**)&d_vs, num_baselines*sizeof(user_precision_t) ));
  ( gpuMalloc( (void**)&d_ws, num_baselines*sizeof(user_precision_t) ));
  ( gpuMemcpy(d_us, us, num_baselines*sizeof(user_precision_t),
                                                        gpuMemcpyHostToDevice ));
  ( gpuMemcpy(d_vs, vs, num_baselines*sizeof(user_precision_t),
                                                        gpuMemcpyHostToDevice ));
  ( gpuMemcpy(d_ws, ws, num_baselines*sizeof(user_precision_t),
                                                        gpuMemcpyHostToDevice ));

  double *d_ls = NULL;
  double *d_ms = NULL;
  double *d_ns = NULL;
  ( gpuMalloc( (void**)&d_ls, num_components*sizeof(double) ));
  ( gpuMalloc( (void**)&d_ms, num_components*sizeof(double) ));
  ( gpuMalloc( (void**)&d_ns, num_components*sizeof(double) ));
  ( gpuMemcpy(d_ls, ls, num_components*sizeof(double),
                                                      gpuMemcpyHostToDevice ));
  ( gpuMemcpy(d_ms, ms, num_components*sizeof(double),
                                                      gpuMemcpyHostToDevice ));
  ( gpuMemcpy(d_ns, ns, num_components*sizeof(double),
                                                      gpuMemcpyHostToDevice ));

  user_precision_complex_t *d_visis = NULL;
  ( gpuMalloc( (void**)&d_visis, num_baselines*num_components*sizeof(user_precision_complex_t) ));

  dim3 grid, threads;

  threads.x = 16;
  threads.y = 16;
  grid.x = (int)ceilf( (float)num_baselines / (float)threads.x );
  grid.y = (int)ceilf( (float)num_components / (float)threads.y );

  gpuErrorCheckKernel("kern_calc_measurement_equation",
                      kern_calc_measurement_equation, grid, threads,
                      num_components, num_baselines,
                      d_us, d_vs, d_ws,
                      d_ls, d_ms, d_ns,
                      (gpuUserComplex*)d_visis );

  ( gpuMemcpy(visis, (user_precision_complex_t*)d_visis, num_components*num_baselines*sizeof(user_precision_complex_t),gpuMemcpyDeviceToHost ));

  ( gpuFree( d_us ) );
  ( gpuFree( d_vs ) );
  ( gpuFree( d_ws ) );
  ( gpuFree( d_ls ) );
  ( gpuFree( d_ms ) );
  ( gpuFree( d_ns ) );
  ( gpuFree(d_visis ) );

}

__global__ void kern_apply_beam_gains_stokesIQUV(int num_gains, gpuUserComplex *d_g1xs,
          gpuUserComplex *d_D1xs,
          gpuUserComplex *d_D1ys, gpuUserComplex *d_g1ys,
          gpuUserComplex *d_g2xs, gpuUserComplex *d_D2xs,
          gpuUserComplex *d_D2ys, gpuUserComplex *d_g2ys,
          user_precision_t *d_flux_Is, user_precision_t *d_flux_Qs,
          user_precision_t *d_flux_Us, user_precision_t *d_flux_Vs,
          gpuUserComplex *d_visi_components,
          gpuUserComplex *d_visi_XXs, gpuUserComplex *d_visi_XYs,
          gpuUserComplex *d_visi_YXs, gpuUserComplex *d_visi_YYs) {

  const int iGain = threadIdx.x + (blockDim.x*blockIdx.x);
  if (iGain < num_gains) {

    gpuUserComplex visi_XX;
    gpuUserComplex visi_XY;
    gpuUserComplex visi_YX;
    gpuUserComplex visi_YY;

    apply_beam_gains_stokesIQUV(d_g1xs[iGain], d_D1xs[iGain],
             d_D1ys[iGain], d_g1ys[iGain],
             d_g2xs[iGain], d_D2xs[iGain],
             d_D2ys[iGain], d_g2ys[iGain],
             d_flux_Is[iGain], d_flux_Qs[iGain],
             d_flux_Us[iGain], d_flux_Vs[iGain],
             d_visi_components[iGain],
             &visi_XX, &visi_XY,
             &visi_YX, &visi_YY);

    d_visi_XXs[iGain] = visi_XX;
    d_visi_XYs[iGain] = visi_XY;
    d_visi_YXs[iGain] = visi_YX;
    d_visi_YYs[iGain] = visi_YY;

  }
}

extern "C" void test_kern_apply_beam_gains(int num_gains, user_precision_complex_t *g1xs,
          user_precision_complex_t *D1xs,
          user_precision_complex_t *D1ys, user_precision_complex_t *g1ys,
          user_precision_complex_t *g2xs, user_precision_complex_t *D2xs,
          user_precision_complex_t *D2ys, user_precision_complex_t *g2ys,
          user_precision_t *flux_Is, user_precision_t *flux_Qs,
          user_precision_t *flux_Us, user_precision_t *flux_Vs,
          user_precision_complex_t *visi_components,
          user_precision_complex_t *visi_XXs, user_precision_complex_t *visi_XYs,
          user_precision_complex_t *visi_YXs, user_precision_complex_t *visi_YYs){

  user_precision_complex_t *d_g1xs = NULL;
  user_precision_complex_t *d_D1xs = NULL;
  user_precision_complex_t *d_D1ys = NULL;
  user_precision_complex_t *d_g1ys = NULL;
  user_precision_complex_t *d_g2xs = NULL;
  user_precision_complex_t *d_D2xs = NULL;
  user_precision_complex_t *d_D2ys = NULL;
  user_precision_complex_t *d_g2ys = NULL;
  user_precision_t *d_flux_Is = NULL;
  user_precision_t *d_flux_Qs = NULL;
  user_precision_t *d_flux_Us = NULL;
  user_precision_t *d_flux_Vs = NULL;
  user_precision_complex_t *d_visi_components = NULL;
  user_precision_complex_t *d_visi_XXs = NULL;
  user_precision_complex_t *d_visi_XYs = NULL;
  user_precision_complex_t *d_visi_YXs = NULL;
  user_precision_complex_t *d_visi_YYs = NULL;

  ( gpuMalloc( (void**)&d_g1xs,
                                  num_gains*sizeof(user_precision_complex_t) ));
  ( gpuMalloc( (void**)&d_D1xs,
                                  num_gains*sizeof(user_precision_complex_t) ));
  ( gpuMalloc( (void**)&d_D1ys,
                                  num_gains*sizeof(user_precision_complex_t) ));
  ( gpuMalloc( (void**)&d_g1ys,
                                  num_gains*sizeof(user_precision_complex_t) ));
  ( gpuMalloc( (void**)&d_g2xs,
                                  num_gains*sizeof(user_precision_complex_t) ));
  ( gpuMalloc( (void**)&d_D2xs,
                                  num_gains*sizeof(user_precision_complex_t) ));
  ( gpuMalloc( (void**)&d_D2ys,
                                  num_gains*sizeof(user_precision_complex_t) ));
  ( gpuMalloc( (void**)&d_g2ys,
                                  num_gains*sizeof(user_precision_complex_t) ));
  ( gpuMalloc( (void**)&d_flux_Is,
                                          num_gains*sizeof(user_precision_t) ));
  ( gpuMalloc( (void**)&d_flux_Qs,
                                          num_gains*sizeof(user_precision_t) ));
  ( gpuMalloc( (void**)&d_flux_Us,
                                          num_gains*sizeof(user_precision_t) ));
  ( gpuMalloc( (void**)&d_flux_Vs,
                                          num_gains*sizeof(user_precision_t) ));
  ( gpuMalloc( (void**)&d_visi_components,
                                  num_gains*sizeof(user_precision_complex_t) ));
  ( gpuMalloc( (void**)&d_visi_XXs,
                                  num_gains*sizeof(user_precision_complex_t) ));
  ( gpuMalloc( (void**)&d_visi_XYs,
                                  num_gains*sizeof(user_precision_complex_t) ));
  ( gpuMalloc( (void**)&d_visi_YXs,
                                  num_gains*sizeof(user_precision_complex_t) ));
  ( gpuMalloc( (void**)&d_visi_YYs,
                                  num_gains*sizeof(user_precision_complex_t) ));

  ( gpuMemcpy(d_g1xs, g1xs,
          num_gains*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice ));
  ( gpuMemcpy(d_D1xs, D1xs,
          num_gains*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice ));
  ( gpuMemcpy(d_D1ys, D1ys,
          num_gains*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice ));
  ( gpuMemcpy(d_g1ys, g1ys,
          num_gains*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice ));
  ( gpuMemcpy(d_g2xs, g2xs,
          num_gains*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice ));
  ( gpuMemcpy(d_D2xs, D2xs,
          num_gains*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice ));
  ( gpuMemcpy(d_D2ys, D2ys,
          num_gains*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice ));
  ( gpuMemcpy(d_g2ys, g2ys,
          num_gains*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice ));
  ( gpuMemcpy(d_visi_components, visi_components,
          num_gains*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice ));
  ( gpuMemcpy(d_visi_XXs, visi_XXs,
          num_gains*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice ));
  ( gpuMemcpy(d_visi_XYs, visi_XYs,
          num_gains*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice ));
  ( gpuMemcpy(d_visi_YXs, visi_YXs,
          num_gains*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice ));
  ( gpuMemcpy(d_visi_YYs, visi_YYs,
          num_gains*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice ));

  ( gpuMemcpy(d_flux_Is, flux_Is,
                             num_gains*sizeof(user_precision_t), gpuMemcpyHostToDevice ));
  ( gpuMemcpy(d_flux_Qs, flux_Qs,
                             num_gains*sizeof(user_precision_t), gpuMemcpyHostToDevice ));
  ( gpuMemcpy(d_flux_Us, flux_Us,
                             num_gains*sizeof(user_precision_t), gpuMemcpyHostToDevice ));
  ( gpuMemcpy(d_flux_Vs, flux_Vs,
                             num_gains*sizeof(user_precision_t), gpuMemcpyHostToDevice ));

  dim3 grid, threads;

  threads.x = 128;
  grid.x = (int)ceil( (user_precision_t)num_gains / (user_precision_t)threads.x );

  gpuErrorCheckKernel("kern_apply_beam_gains_stokesIQUV",
                      kern_apply_beam_gains_stokesIQUV, grid, threads,
                      num_gains,
                      (gpuUserComplex *)d_g1xs, (gpuUserComplex *)d_D1xs,
                      (gpuUserComplex *)d_D1ys, (gpuUserComplex *)d_g1ys,
                      (gpuUserComplex *)d_g2xs, (gpuUserComplex *)d_D2xs,
                      (gpuUserComplex *)d_D2ys, (gpuUserComplex *)d_g2ys,
                      d_flux_Is, d_flux_Qs,
                      d_flux_Us, d_flux_Vs,
                      (gpuUserComplex *)d_visi_components,
                      (gpuUserComplex *)d_visi_XXs, (gpuUserComplex *)d_visi_XYs,
                      (gpuUserComplex *)d_visi_YXs, (gpuUserComplex *)d_visi_YYs );

  ( gpuMemcpy(visi_XXs, d_visi_XXs,
           num_gains*sizeof(user_precision_complex_t),gpuMemcpyDeviceToHost ));
  ( gpuMemcpy(visi_XYs, d_visi_XYs,
           num_gains*sizeof(user_precision_complex_t),gpuMemcpyDeviceToHost ));
  ( gpuMemcpy(visi_YXs, d_visi_YXs,
           num_gains*sizeof(user_precision_complex_t),gpuMemcpyDeviceToHost ));
  ( gpuMemcpy(visi_YYs, d_visi_YYs,
           num_gains*sizeof(user_precision_complex_t),gpuMemcpyDeviceToHost ));

  ( gpuFree( d_g1xs ) );
  ( gpuFree( d_D1xs ) );
  ( gpuFree( d_D1ys ) );
  ( gpuFree( d_g1ys ) );
  ( gpuFree( d_g2xs ) );
  ( gpuFree( d_D2xs ) );
  ( gpuFree( d_D2ys ) );
  ( gpuFree( d_g2ys ) );
  ( gpuFree( d_flux_Is ) );
  ( gpuFree( d_flux_Qs ) );
  ( gpuFree( d_flux_Us ) );
  ( gpuFree( d_flux_Vs ) );
  ( gpuFree( d_visi_components ) );
  ( gpuFree( d_visi_XXs ) );
  ( gpuFree( d_visi_XYs ) );
  ( gpuFree( d_visi_YXs ) );
  ( gpuFree( d_visi_YYs ) );

}

__global__ void kern_get_beam_gains(int num_components, int num_baselines,
           int num_freqs, int num_cross, int num_times, int beamtype,
           gpuUserComplex *d_g1xs, gpuUserComplex *d_D1xs,
           gpuUserComplex *d_D1ys, gpuUserComplex *d_g1ys,
           gpuUserComplex *d_recov_g1x, gpuUserComplex *d_recov_D1x,
           gpuUserComplex *d_recov_D1y, gpuUserComplex *d_recov_g1y,
           gpuUserComplex *d_recov_g2x, gpuUserComplex *d_recov_D2x,
           gpuUserComplex *d_recov_D2y, gpuUserComplex *d_recov_g2y,
           int use_twobeams, int num_ants,
           int *d_ant1_to_baseline_map, int *d_ant2_to_baseline_map) {

  // Start by computing which baseline we're going to do
  // This iBaseline means all baselines for all times and freqs
  const int iBaseline = threadIdx.x + (blockDim.x*blockIdx.x);
  if(iBaseline < num_cross) {

    for (int iComponent = 0; iComponent < num_components; iComponent++) {

      gpuUserComplex g1x;
      gpuUserComplex D1x;
      gpuUserComplex D1y;
      gpuUserComplex g1y;
      gpuUserComplex g2x;
      gpuUserComplex D2x;
      gpuUserComplex D2y;
      gpuUserComplex g2y;

      if (use_twobeams == 1) {
        get_beam_gains_multibeams(iBaseline, iComponent, num_freqs,
                 num_baselines, num_components, num_times, beamtype,
                 d_g1xs, d_D1xs,
                 d_D1ys, d_g1ys,
                 d_ant1_to_baseline_map, d_ant2_to_baseline_map,
                 &g1x, &D1x, &D1y, &g1y, &g2x, &D2x, &D2y, &g2y);
      } else {
        get_beam_gains(iBaseline, iComponent, num_freqs,
                 num_baselines, num_components, num_times, beamtype,
                 d_g1xs, d_D1xs,
                 d_D1ys, d_g1ys,
                 &g1x, &D1x, &D1y, &g1y, &g2x, &D2x, &D2y, &g2y);
      }

      int out_ind = num_cross*iComponent + iBaseline;

      d_recov_g1x[out_ind] = g1x;
      d_recov_D1x[out_ind] = D1x;
      d_recov_D1y[out_ind] = D1y;
      d_recov_g1y[out_ind] = g1y;
      d_recov_g2x[out_ind] = g2x;
      d_recov_D2x[out_ind] = D2x;
      d_recov_D2y[out_ind] = D2y;
      d_recov_g2y[out_ind] = g2y;

    }
  }
}

extern "C" void test_kern_get_beam_gains(int num_freqs, int num_cross,
          int num_baselines, int num_components, int num_times, int beamtype,
          user_precision_complex_t *primay_beam_J00, user_precision_complex_t *primay_beam_J01,
          user_precision_complex_t *primay_beam_J10, user_precision_complex_t *primay_beam_J11,
          user_precision_complex_t *recover_g1x, user_precision_complex_t *recover_D1x,
          user_precision_complex_t *recover_D1y, user_precision_complex_t *recover_g1y,
          user_precision_complex_t *recover_g2x, user_precision_complex_t *recover_D2x,
          user_precision_complex_t *recover_D2y, user_precision_complex_t *recover_g2y,
          int use_twobeams, int num_ants){

  user_precision_complex_t *d_recover_g1x = NULL;
  user_precision_complex_t *d_recover_D1x = NULL;
  user_precision_complex_t *d_recover_D1y = NULL;
  user_precision_complex_t *d_recover_g1y = NULL;
  user_precision_complex_t *d_recover_g2x = NULL;
  user_precision_complex_t *d_recover_D2x = NULL;
  user_precision_complex_t *d_recover_D2y = NULL;
  user_precision_complex_t *d_recover_g2y = NULL;

  user_precision_complex_t *d_g1xs = NULL;
  user_precision_complex_t *d_D1xs = NULL;
  user_precision_complex_t *d_D1ys = NULL;
  user_precision_complex_t *d_g1ys = NULL;

  int num_recover_gains = num_components*num_cross;
  int num_input_gains = num_freqs*num_times*num_components*num_ants;

  ( gpuMalloc( (void**)&d_recover_g1x, num_recover_gains*sizeof(user_precision_complex_t) ));
  ( gpuMalloc( (void**)&d_recover_D1x, num_recover_gains*sizeof(user_precision_complex_t) ));
  ( gpuMalloc( (void**)&d_recover_D1y, num_recover_gains*sizeof(user_precision_complex_t) ));
  ( gpuMalloc( (void**)&d_recover_g1y, num_recover_gains*sizeof(user_precision_complex_t) ));
  ( gpuMalloc( (void**)&d_recover_g2x, num_recover_gains*sizeof(user_precision_complex_t) ));
  ( gpuMalloc( (void**)&d_recover_D2x, num_recover_gains*sizeof(user_precision_complex_t) ));
  ( gpuMalloc( (void**)&d_recover_D2y, num_recover_gains*sizeof(user_precision_complex_t) ));
  ( gpuMalloc( (void**)&d_recover_g2y, num_recover_gains*sizeof(user_precision_complex_t) ));

  ( gpuMalloc( (void**)&d_g1xs, num_input_gains*sizeof(user_precision_complex_t) ));
  ( gpuMalloc( (void**)&d_D1xs, num_input_gains*sizeof(user_precision_complex_t) ));
  ( gpuMalloc( (void**)&d_D1ys, num_input_gains*sizeof(user_precision_complex_t) ));
  ( gpuMalloc( (void**)&d_g1ys, num_input_gains*sizeof(user_precision_complex_t) ));

  ( gpuMemcpy(d_g1xs, primay_beam_J00, num_input_gains*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice ));
  ( gpuMemcpy(d_D1xs, primay_beam_J01, num_input_gains*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice ));
  ( gpuMemcpy(d_D1ys, primay_beam_J10, num_input_gains*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice ));
  ( gpuMemcpy(d_g1ys, primay_beam_J11, num_input_gains*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice ));

  dim3 grid, threads;

  threads.x = 128;
  grid.x = (int)ceil( (user_precision_t)num_cross / (user_precision_t)threads.x );

  //These are only needed when we're actually grabbing information for two
  //antennas instead of one
  int *d_ant1_to_baseline_map = NULL;
  int *d_ant2_to_baseline_map = NULL;

  if (use_twobeams == 1) {

    ( gpuMalloc( (void**)&d_ant1_to_baseline_map, num_baselines*sizeof(int) ));
    ( gpuMalloc( (void**)&d_ant2_to_baseline_map, num_baselines*sizeof(int) ));

    //Fill in the indexes of antenna1 and antenna2 for all cross-correlation combos
    fill_ant_to_baseline_mapping(num_ants, d_ant1_to_baseline_map,
                                           d_ant2_to_baseline_map);
  }

  gpuErrorCheckKernel("kern_get_beam_gains",
                      kern_get_beam_gains, grid, threads,
                      num_components, num_baselines,
                      num_freqs, num_cross, num_times, beamtype,
                      (gpuUserComplex *)d_g1xs,
                      (gpuUserComplex *)d_D1xs,
                      (gpuUserComplex *)d_D1ys,
                      (gpuUserComplex *)d_g1ys,
                      (gpuUserComplex *)d_recover_g1x, (gpuUserComplex *)d_recover_D1x,
                      (gpuUserComplex *)d_recover_D1y, (gpuUserComplex *)d_recover_g1y,
                      (gpuUserComplex *)d_recover_g2x, (gpuUserComplex *)d_recover_D2x,
                      (gpuUserComplex *)d_recover_D2y, (gpuUserComplex *)d_recover_g2y,
                      use_twobeams, num_ants,
                      d_ant1_to_baseline_map, d_ant2_to_baseline_map);

  ( gpuMemcpy(recover_g1x, d_recover_g1x, num_recover_gains*sizeof(user_precision_complex_t), gpuMemcpyDeviceToHost ));
  ( gpuMemcpy(recover_D1x, d_recover_D1x, num_recover_gains*sizeof(user_precision_complex_t), gpuMemcpyDeviceToHost ));
  ( gpuMemcpy(recover_D1y, d_recover_D1y, num_recover_gains*sizeof(user_precision_complex_t), gpuMemcpyDeviceToHost ));
  ( gpuMemcpy(recover_g1y, d_recover_g1y, num_recover_gains*sizeof(user_precision_complex_t), gpuMemcpyDeviceToHost ));
  ( gpuMemcpy(recover_g2x, d_recover_g2x, num_recover_gains*sizeof(user_precision_complex_t), gpuMemcpyDeviceToHost ));
  ( gpuMemcpy(recover_D2x, d_recover_D2x, num_recover_gains*sizeof(user_precision_complex_t), gpuMemcpyDeviceToHost ));
  ( gpuMemcpy(recover_D2y, d_recover_D2y, num_recover_gains*sizeof(user_precision_complex_t), gpuMemcpyDeviceToHost ));
  ( gpuMemcpy(recover_g2y, d_recover_g2y, num_recover_gains*sizeof(user_precision_complex_t), gpuMemcpyDeviceToHost ));

  ( gpuFree( d_recover_g1x ) );
  ( gpuFree( d_recover_D1x ) );
  ( gpuFree( d_recover_D1y ) );
  ( gpuFree( d_recover_g1y ) );
  ( gpuFree( d_recover_g2x ) );
  ( gpuFree( d_recover_D2x ) );
  ( gpuFree( d_recover_D2y ) );
  ( gpuFree( d_recover_g2y ) );

  ( gpuFree( d_g1xs ) );
  ( gpuFree( d_D1xs ) );
  ( gpuFree( d_D1ys ) );
  ( gpuFree( d_g1ys ) );

  if (use_twobeams == 1) {
    // free(ant1_to_baseline_map);
    // free(ant2_to_baseline_map);
    ( gpuFree( d_ant1_to_baseline_map ) );
    ( gpuFree( d_ant2_to_baseline_map ) );
  }

}

__global__ void kern_update_sum_visis_stokesIQUV(int num_freqs,
     int num_baselines, int num_components, int num_times, int beamtype,
     gpuUserComplex *d_g1xs, gpuUserComplex *d_D1xs,
     gpuUserComplex *d_D1ys, gpuUserComplex *d_g1ys,
     int *d_ant1_to_baseline_map, int *d_ant2_to_baseline_map, int use_twobeams,
     gpuUserComplex *d_visi_components,
     user_precision_t *d_flux_I, user_precision_t *d_flux_Q,
     user_precision_t *d_flux_U, user_precision_t *d_flux_V,
     user_precision_t *d_sum_visi_XX_real, user_precision_t *d_sum_visi_XX_imag,
     user_precision_t *d_sum_visi_XY_real, user_precision_t *d_sum_visi_XY_imag,
     user_precision_t *d_sum_visi_YX_real, user_precision_t *d_sum_visi_YX_imag,
     user_precision_t *d_sum_visi_YY_real, user_precision_t *d_sum_visi_YY_imag) {

  // Start by computing which baseline we're going to do
  const int iBaseline = threadIdx.x + (blockDim.x*blockIdx.x);

  if(iBaseline < num_freqs*num_baselines*num_times) {

    int time_ind = (int)floorf( (user_precision_t)iBaseline / ((user_precision_t)num_baselines * (user_precision_t)num_freqs));
    int freq_ind = (int)floorf( ((user_precision_t)iBaseline - ((user_precision_t)time_ind*(user_precision_t)num_baselines * (user_precision_t)num_freqs)) / (user_precision_t)num_baselines);

    for (int iComponent = 0; iComponent < num_components; iComponent++) {

      //There is a flux for every frequnecy and component
      int flux_ind = num_components*freq_ind + iComponent;

      update_sum_visis_stokesIQUV(iBaseline, iComponent, num_freqs,
             num_baselines, num_components, num_times, beamtype,
             d_g1xs, d_D1xs,
             d_D1ys, d_g1ys,
             d_ant1_to_baseline_map, d_ant2_to_baseline_map, use_twobeams,
             d_visi_components[iBaseline],
             d_flux_I[flux_ind], d_flux_Q[flux_ind],
             d_flux_U[flux_ind], d_flux_V[flux_ind],
             d_sum_visi_XX_real, d_sum_visi_XX_imag,
             d_sum_visi_XY_real, d_sum_visi_XY_imag,
             d_sum_visi_YX_real, d_sum_visi_YX_imag,
             d_sum_visi_YY_real, d_sum_visi_YY_imag);

    }
  }
}

extern "C" void test_kern_update_sum_visis(int num_freqs, int num_cross,
          int num_baselines, int num_components, int num_times, int beamtype,
          int use_twobeams, int num_ants,
          user_precision_complex_t *primay_beam_J00,
          user_precision_complex_t *primay_beam_J01,
          user_precision_complex_t *primay_beam_J10,
          user_precision_complex_t *primay_beam_J11,
          user_precision_complex_t *visi_components,
          user_precision_t *flux_I, user_precision_t *flux_Q,
          user_precision_t *flux_U, user_precision_t *flux_V,
          user_precision_t *sum_visi_XX_real, user_precision_t *sum_visi_XX_imag,
          user_precision_t *sum_visi_XY_real, user_precision_t *sum_visi_XY_imag,
          user_precision_t *sum_visi_YX_real, user_precision_t *sum_visi_YX_imag,
          user_precision_t *sum_visi_YY_real, user_precision_t *sum_visi_YY_imag){

  user_precision_complex_t *d_gxs = NULL;
  user_precision_complex_t *d_Dxs = NULL;
  user_precision_complex_t *d_Dys = NULL;
  user_precision_complex_t *d_gys = NULL;
  user_precision_complex_t *d_visi_components = NULL;

  //if do_ants, need to malloc more gains

  int num_input_gains = num_freqs*num_times*num_components*num_ants;

  ( gpuMalloc( (void**)&d_gxs,
                    num_input_gains*sizeof(user_precision_complex_t) ));
  ( gpuMalloc( (void**)&d_Dxs,
                    num_input_gains*sizeof(user_precision_complex_t) ));
  ( gpuMalloc( (void**)&d_Dys,
                    num_input_gains*sizeof(user_precision_complex_t) ));
  ( gpuMalloc( (void**)&d_gys,
                    num_input_gains*sizeof(user_precision_complex_t) ));
  ( gpuMalloc( (void**)&d_visi_components,
                    num_cross*sizeof(user_precision_complex_t) ));

  ( gpuMemcpy(d_gxs, primay_beam_J00,
            num_input_gains*sizeof(user_precision_complex_t),
            gpuMemcpyHostToDevice ));
  ( gpuMemcpy(d_Dxs, primay_beam_J01,
            num_input_gains*sizeof(user_precision_complex_t),
            gpuMemcpyHostToDevice ));
  ( gpuMemcpy(d_Dys, primay_beam_J10,
            num_input_gains*sizeof(user_precision_complex_t),
            gpuMemcpyHostToDevice ));
  ( gpuMemcpy(d_gys, primay_beam_J11,
            num_input_gains*sizeof(user_precision_complex_t),
            gpuMemcpyHostToDevice ));
  ( gpuMemcpy(d_visi_components, visi_components,
                                     num_cross*sizeof(user_precision_complex_t),
                                     gpuMemcpyHostToDevice ));

  user_precision_t *d_flux_I = NULL;
  user_precision_t *d_flux_Q = NULL;
  user_precision_t *d_flux_U = NULL;
  user_precision_t *d_flux_V = NULL;

  ( gpuMalloc( (void**)&d_flux_I, num_components*num_times*num_freqs*sizeof(user_precision_t) ));
  ( gpuMalloc( (void**)&d_flux_Q, num_components*num_times*num_freqs*sizeof(user_precision_t) ));
  ( gpuMalloc( (void**)&d_flux_U, num_components*num_times*num_freqs*sizeof(user_precision_t) ));
  ( gpuMalloc( (void**)&d_flux_V, num_components*num_times*num_freqs*sizeof(user_precision_t) ));

  ( gpuMemcpy(d_flux_I, flux_I,
                    num_components*num_times*num_freqs*sizeof(user_precision_t),    gpuMemcpyHostToDevice ));
  ( gpuMemcpy(d_flux_Q, flux_Q,
                    num_components*num_times*num_freqs*sizeof(user_precision_t),    gpuMemcpyHostToDevice ));
  ( gpuMemcpy(d_flux_U, flux_U,
                    num_components*num_times*num_freqs*sizeof(user_precision_t),    gpuMemcpyHostToDevice ));
  ( gpuMemcpy(d_flux_V, flux_V,
                    num_components*num_times*num_freqs*sizeof(user_precision_t),    gpuMemcpyHostToDevice ));

  user_precision_t *d_sum_visi_XX_real = NULL;
  user_precision_t *d_sum_visi_XY_real = NULL;
  user_precision_t *d_sum_visi_YX_real = NULL;
  user_precision_t *d_sum_visi_YY_real = NULL;
  user_precision_t *d_sum_visi_XX_imag = NULL;
  user_precision_t *d_sum_visi_XY_imag = NULL;
  user_precision_t *d_sum_visi_YX_imag = NULL;
  user_precision_t *d_sum_visi_YY_imag = NULL;

  ( gpuMalloc( (void**)&d_sum_visi_XX_real,
                                          num_cross*sizeof(user_precision_t) ));
  ( gpuMalloc( (void**)&d_sum_visi_XY_real,
                                          num_cross*sizeof(user_precision_t) ));
  ( gpuMalloc( (void**)&d_sum_visi_YX_real,
                                          num_cross*sizeof(user_precision_t) ));
  ( gpuMalloc( (void**)&d_sum_visi_YY_real,
                                          num_cross*sizeof(user_precision_t) ));
  ( gpuMalloc( (void**)&d_sum_visi_XX_imag,
                                          num_cross*sizeof(user_precision_t) ));
  ( gpuMalloc( (void**)&d_sum_visi_XY_imag,
                                          num_cross*sizeof(user_precision_t) ));
  ( gpuMalloc( (void**)&d_sum_visi_YX_imag,
                                          num_cross*sizeof(user_precision_t) ));
  ( gpuMalloc( (void**)&d_sum_visi_YY_imag,
                                          num_cross*sizeof(user_precision_t) ));
  gpuMemcpy(d_sum_visi_XX_real, sum_visi_XX_real,
                    num_cross*sizeof(user_precision_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_sum_visi_XY_real, sum_visi_XY_real,
                    num_cross*sizeof(user_precision_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_sum_visi_YX_real, sum_visi_YX_real,
                    num_cross*sizeof(user_precision_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_sum_visi_YY_real, sum_visi_YY_real,
                    num_cross*sizeof(user_precision_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_sum_visi_XX_imag, sum_visi_XX_imag,
                    num_cross*sizeof(user_precision_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_sum_visi_XY_imag, sum_visi_XY_imag,
                    num_cross*sizeof(user_precision_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_sum_visi_YX_imag, sum_visi_YX_imag,
                    num_cross*sizeof(user_precision_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_sum_visi_YY_imag, sum_visi_YY_imag,
                    num_cross*sizeof(user_precision_t), gpuMemcpyHostToDevice );
  


  //These are only needed when we're actually grabbing information for two
  //antennas instead of one
  int *d_ant1_to_baseline_map = NULL;
  int *d_ant2_to_baseline_map = NULL;

  if (use_twobeams == 1) {

    ( gpuMalloc( (void**)&d_ant1_to_baseline_map, num_baselines*sizeof(int) ));
    ( gpuMalloc( (void**)&d_ant2_to_baseline_map, num_baselines*sizeof(int) ));

    //Fill in the indexes of antenna1 and antenna2 for all cross-correlation combos
    fill_ant_to_baseline_mapping(num_ants, d_ant1_to_baseline_map,
                                           d_ant2_to_baseline_map);
  }

  dim3 grid, threads;

  threads.x = 128;
  grid.x = (int)ceil( (user_precision_t)num_cross / (user_precision_t)threads.x );

  gpuErrorCheckKernel("kern_update_sum_visis_stokesIQUV",
                      kern_update_sum_visis_stokesIQUV, grid, threads,
                      num_freqs, num_baselines, num_components, num_times, beamtype,
                      (gpuUserComplex *)d_gxs, (gpuUserComplex *)d_Dxs,
                      (gpuUserComplex *)d_Dys, (gpuUserComplex *)d_gys,
                      d_ant1_to_baseline_map, d_ant2_to_baseline_map, use_twobeams,
                      (gpuUserComplex *)d_visi_components,
                      d_flux_I, d_flux_Q, d_flux_U, d_flux_V,
                      d_sum_visi_XX_real, d_sum_visi_XX_imag,
                      d_sum_visi_XY_real, d_sum_visi_XY_imag,
                      d_sum_visi_YX_real, d_sum_visi_YX_imag,
                      d_sum_visi_YY_real, d_sum_visi_YY_imag );

  ( gpuMemcpy(sum_visi_XX_real, d_sum_visi_XX_real,
                  num_cross*sizeof(user_precision_t), gpuMemcpyDeviceToHost ));
  ( gpuMemcpy(sum_visi_XY_real, d_sum_visi_XY_real,
                  num_cross*sizeof(user_precision_t), gpuMemcpyDeviceToHost ));
  ( gpuMemcpy(sum_visi_YX_real, d_sum_visi_YX_real,
                  num_cross*sizeof(user_precision_t), gpuMemcpyDeviceToHost ));
  ( gpuMemcpy(sum_visi_YY_real, d_sum_visi_YY_real,
                  num_cross*sizeof(user_precision_t), gpuMemcpyDeviceToHost ));
  ( gpuMemcpy(sum_visi_XX_imag, d_sum_visi_XX_imag,
                  num_cross*sizeof(user_precision_t), gpuMemcpyDeviceToHost ));
  ( gpuMemcpy(sum_visi_XY_imag, d_sum_visi_XY_imag,
                  num_cross*sizeof(user_precision_t), gpuMemcpyDeviceToHost ));
  ( gpuMemcpy(sum_visi_YX_imag, d_sum_visi_YX_imag,
                  num_cross*sizeof(user_precision_t), gpuMemcpyDeviceToHost ));
  ( gpuMemcpy(sum_visi_YY_imag, d_sum_visi_YY_imag,
                  num_cross*sizeof(user_precision_t), gpuMemcpyDeviceToHost ));

  ( gpuFree( d_gxs ) );
  ( gpuFree( d_Dxs ) );
  ( gpuFree( d_Dys ) );
  ( gpuFree( d_gys ) );
  ( gpuFree( d_visi_components ) );
  ( gpuFree( d_flux_I ) );
  ( gpuFree( d_flux_Q ) );
  ( gpuFree( d_flux_U ) );
  ( gpuFree( d_flux_V ) );
  ( gpuFree( d_sum_visi_XX_real ) );
  ( gpuFree( d_sum_visi_XY_real ) );
  ( gpuFree( d_sum_visi_YX_real ) );
  ( gpuFree( d_sum_visi_YY_real ) );
  ( gpuFree( d_sum_visi_XX_imag ) );
  ( gpuFree( d_sum_visi_XY_imag ) );
  ( gpuFree( d_sum_visi_YX_imag ) );
  ( gpuFree( d_sum_visi_YY_imag ) );

  if (use_twobeams == 1) {
    ( gpuFree( d_ant1_to_baseline_map ) );
    ( gpuFree( d_ant2_to_baseline_map ) );
  }

}

//just make things zero pls
__global__ void kern_make_zeros(user_precision_t *array, int num_arr) {

  const int iComp = threadIdx.x + (blockDim.x*blockIdx.x);

  if (iComp < num_arr)
  {
    array[iComp] = 0.0;
  }
}


extern "C" void test_source_component_common(int num_of_each_flux_type,
           components_t components,
           double *freqs, woden_settings_t *woden_settings,
           beam_settings_t *beam_settings,
           user_precision_complex_t *gxs, user_precision_complex_t *Dxs,
           user_precision_complex_t *Dys, user_precision_complex_t *gys,
           user_precision_t *extrap_flux_I, user_precision_t *extrap_flux_Q,
           user_precision_t *extrap_flux_U, user_precision_t *extrap_flux_V,
           double *ls, double *ms, double *ns,
           e_component_type comptype){

  woden_settings->do_QUV = 1;

  source_t *chunked_source = (source_t *)malloc(sizeof(source_t));

  //TODODOD have a if (comptype == POINT) etc here so we can use same
  //componenets to test all POINT, GAUSSIAN, SHAPELET

  int NUM_FLUX_TYPES = 3;

  if (comptype == POINT) {
    chunked_source->point_components = components;
    chunked_source->n_points = NUM_FLUX_TYPES*num_of_each_flux_type;
    chunked_source->n_point_powers = num_of_each_flux_type;
    chunked_source->n_point_curves = num_of_each_flux_type;
    chunked_source->n_point_lists = num_of_each_flux_type;

    chunked_source->n_gauss = 0;
    chunked_source->n_gauss_lists = 0;
    chunked_source->n_gauss_powers = 0;
    chunked_source->n_gauss_curves = 0;
    chunked_source->n_shapes = 0;
    chunked_source->n_shape_lists = 0;
    chunked_source->n_shape_powers = 0;
    chunked_source->n_shape_curves = 0;
    chunked_source->n_shape_coeffs = 0;
  }
  else if (comptype == GAUSSIAN) {
    chunked_source->gauss_components = components;
    chunked_source->n_gauss = NUM_FLUX_TYPES*num_of_each_flux_type;
    chunked_source->n_gauss_powers = num_of_each_flux_type;
    chunked_source->n_gauss_curves = num_of_each_flux_type;
    chunked_source->n_gauss_lists = num_of_each_flux_type;

    chunked_source->n_points = 0;
    chunked_source->n_point_lists = 0;
    chunked_source->n_point_powers = 0;
    chunked_source->n_point_curves = 0;
    chunked_source->n_shapes = 0;
    chunked_source->n_shape_lists = 0;
    chunked_source->n_shape_powers = 0;
    chunked_source->n_shape_curves = 0;
    chunked_source->n_shape_coeffs = 0;
  }
  else if (comptype == SHAPELET) {
    chunked_source->shape_components = components;
    chunked_source->n_shapes = NUM_FLUX_TYPES*num_of_each_flux_type;
    chunked_source->n_shape_powers = num_of_each_flux_type;
    chunked_source->n_shape_curves = num_of_each_flux_type;
    chunked_source->n_shape_lists = num_of_each_flux_type;
    chunked_source->n_shape_coeffs = num_of_each_flux_type;

    chunked_source->n_points = 0;
    chunked_source->n_point_lists = 0;
    chunked_source->n_point_powers = 0;
    chunked_source->n_point_curves = 0;
    chunked_source->n_gauss = 0;
    chunked_source->n_gauss_lists = 0;
    chunked_source->n_gauss_powers = 0;
    chunked_source->n_gauss_curves = 0;

  }

  source_t *d_chunked_source = copy_chunked_source_to_GPU(chunked_source);

  double *d_freqs = NULL;
  ( gpuMalloc( (void**)&d_freqs,
                                     woden_settings->num_freqs*sizeof(double) ));
  ( gpuMemcpy( d_freqs, freqs,
             woden_settings->num_freqs*sizeof(double), gpuMemcpyHostToDevice) );

  d_beam_gains_t d_beam_gains;
  visibility_set_t *d_visibility_set = NULL;
  woden_settings->do_autos = 0;

  //THIS IS THE CUDA CALL ACTUALLY BEING TESTED JEEZUS--------------------------
  source_component_common(woden_settings, beam_settings, d_freqs,
       chunked_source, d_chunked_source, &d_beam_gains, comptype,
       d_visibility_set);
  //THIS IS THE CUDA CALL ACTUALLY BEING TESTED JEEZUS--------------------------

  int num_beam_values = NUM_FLUX_TYPES*num_of_each_flux_type*woden_settings->num_freqs*woden_settings->num_time_steps;

  if (woden_settings->use_dipamps == 1) {
    num_beam_values *= woden_settings->num_ants;
  }

  ( gpuMemcpy(gxs, (user_precision_complex_t*)d_beam_gains.d_gxs,
              num_beam_values*sizeof(gpuUserComplex), gpuMemcpyDeviceToHost ));

  ( gpuMemcpy(gys, (user_precision_complex_t*)d_beam_gains.d_gys,
              num_beam_values*sizeof(gpuUserComplex), gpuMemcpyDeviceToHost ));

  if (beam_settings->beamtype == FEE_BEAM || beam_settings->beamtype == FEE_BEAM_INTERP || beam_settings->beamtype == MWA_ANALY) {
    ( gpuMemcpy(Dxs, (user_precision_complex_t*)d_beam_gains.d_Dxs,
                num_beam_values*sizeof(gpuUserComplex), gpuMemcpyDeviceToHost ));
    ( gpuMemcpy(Dys, (user_precision_complex_t*)d_beam_gains.d_Dys,
                num_beam_values*sizeof(gpuUserComplex), gpuMemcpyDeviceToHost ));
  }

  //Just a little shorthand so don't have to keep writing out as much in the
  //memcpy below

  components_t d_components;

  if (comptype == POINT) {
    d_components = d_chunked_source->point_components;
  }
  else if (comptype == GAUSSIAN) {
    d_components = d_chunked_source->gauss_components;
  }
  else {
    d_components = d_chunked_source->shape_components;
  }


  ( gpuMemcpy(ls, d_components.ls,
                            NUM_FLUX_TYPES*num_of_each_flux_type*sizeof(double),
                            gpuMemcpyDeviceToHost ));
  ( gpuMemcpy(ms, d_components.ms,
                            NUM_FLUX_TYPES*num_of_each_flux_type*sizeof(double),
                            gpuMemcpyDeviceToHost ));
  ( gpuMemcpy(ns, d_components.ns,
                            NUM_FLUX_TYPES*num_of_each_flux_type*sizeof(double),
                            gpuMemcpyDeviceToHost ));

  //until we get RM synthesis working, do this for testing
  //do this because I don't want to cut out all the memcpying below, laaazy
  dim3 grid, threads;

  int num_things = NUM_FLUX_TYPES*num_of_each_flux_type*woden_settings->num_freqs;

  threads.x = 128;
  threads.y = 1;
  grid.x = (int)ceil( (float)num_things / (float)threads.x );
  grid.y = 1;

  gpuErrorCheckKernel("kern_make_zeros",
            kern_make_zeros, grid, threads,
            d_components.extrap_stokesQ, num_things);
  gpuErrorCheckKernel("kern_make_zeros",
            kern_make_zeros, grid, threads,
            d_components.extrap_stokesU, num_things);
  gpuErrorCheckKernel("kern_make_zeros",
            kern_make_zeros, grid, threads,
            d_components.extrap_stokesV, num_things);

  ( gpuMemcpy(extrap_flux_I, d_components.extrap_stokesI,
  NUM_FLUX_TYPES*num_of_each_flux_type*woden_settings->num_freqs*sizeof(user_precision_t),
                                                      gpuMemcpyDeviceToHost ));
  ( gpuMemcpy(extrap_flux_Q, d_components.extrap_stokesQ,
  NUM_FLUX_TYPES*num_of_each_flux_type*woden_settings->num_freqs*sizeof(user_precision_t),
                                                      gpuMemcpyDeviceToHost ));
  ( gpuMemcpy(extrap_flux_U, d_components.extrap_stokesU,
  NUM_FLUX_TYPES*num_of_each_flux_type*woden_settings->num_freqs*sizeof(user_precision_t),
                                                      gpuMemcpyDeviceToHost ));
  ( gpuMemcpy(extrap_flux_V, d_components.extrap_stokesV,
  NUM_FLUX_TYPES*num_of_each_flux_type*woden_settings->num_freqs*sizeof(user_precision_t),
                                                      gpuMemcpyDeviceToHost ));

  ( gpuFree( d_freqs ) );
  free_extrapolated_flux_arrays(&d_components, woden_settings->do_QUV);
  free_d_components(d_chunked_source, comptype);
  free_beam_gains(d_beam_gains, beam_settings->beamtype);
}


void malloc_lmn_arrays(source_t *d_chunked_source, components_t *components,
                        int num_components, e_component_type comptype){
  components_t *d_components=NULL;
  if (comptype == POINT) {
    d_components = &d_chunked_source->point_components;
  } else if (comptype == GAUSSIAN) {
    d_components = &d_chunked_source->gauss_components;
  // } else if (comptype == SHAPELET) {
  } else {
    d_components = &d_chunked_source->shape_components;
  }

  ( gpuMalloc( (void**)&d_components->ls,
                                          num_components*sizeof(double) ) );
  ( gpuMalloc( (void**)&d_components->ms,
                                          num_components*sizeof(double) ) );
  ( gpuMalloc( (void**)&d_components->ns,
                                          num_components*sizeof(double) ) );

  ( gpuMemcpy(d_components->ls, components->ls, num_components*sizeof(double),
                                           gpuMemcpyHostToDevice ));
  ( gpuMemcpy(d_components->ms, components->ms, num_components*sizeof(double),
                                           gpuMemcpyHostToDevice ));
  ( gpuMemcpy(d_components->ns, components->ns, num_components*sizeof(double),
                                           gpuMemcpyHostToDevice ));

}

extern "C" void test_kern_calc_visi_all(int n_powers, int n_curves, int n_lists,
          int num_baselines, int num_shape_coeffs,
          int num_freqs, int num_cross, int num_times,
          e_beamtype beamtype, e_component_type comptype,
          components_t components, double *extrap_freqs,
          user_precision_t *us, user_precision_t *vs, user_precision_t *ws,
          user_precision_t *u_shapes, user_precision_t *v_shapes,
          user_precision_t *sum_visi_XX_real, user_precision_t *sum_visi_XX_imag,
          user_precision_t *sum_visi_XY_real, user_precision_t *sum_visi_XY_imag,
          user_precision_t *sum_visi_YX_real, user_precision_t *sum_visi_YX_imag,
          user_precision_t *sum_visi_YY_real, user_precision_t *sum_visi_YY_imag,
          user_precision_t *allsteps_wavelengths, user_precision_t *sbf,
          user_precision_complex_t *gxs, user_precision_complex_t *Dxs,
          user_precision_complex_t *Dys, user_precision_complex_t *gys){

  int do_QUV = 0;
  
  int num_components = n_powers + n_curves + n_lists;

  user_precision_t *d_us = NULL;
  user_precision_t *d_vs = NULL;
  user_precision_t *d_ws = NULL;
  user_precision_t *d_allsteps_wavelengths = NULL;

  ( gpuMalloc( (void**)&d_us, num_cross*sizeof(user_precision_t) ) );
  ( gpuMalloc( (void**)&d_vs, num_cross*sizeof(user_precision_t) ) );
  ( gpuMalloc( (void**)&d_ws, num_cross*sizeof(user_precision_t) ) );
  ( gpuMalloc( (void**)&d_allsteps_wavelengths, num_cross*sizeof(user_precision_t) ) );

  ( gpuMemcpy(d_us, us,
                             num_cross*sizeof(user_precision_t), gpuMemcpyHostToDevice ));
  ( gpuMemcpy(d_vs, vs,
                             num_cross*sizeof(user_precision_t), gpuMemcpyHostToDevice ));
  ( gpuMemcpy(d_ws, ws,
                             num_cross*sizeof(user_precision_t), gpuMemcpyHostToDevice ));
  ( gpuMemcpy(d_allsteps_wavelengths, allsteps_wavelengths,
                             num_cross*sizeof(user_precision_t), gpuMemcpyHostToDevice ));

  //Here are many things that would have been done by source_component_common
  source_t *chunked_source = (source_t *)malloc(sizeof(source_t));

  double *d_extrap_freqs = NULL;
  ( gpuMalloc( (void**)&d_extrap_freqs,
                                   num_freqs*sizeof(double) ));
  ( gpuMemcpy(d_extrap_freqs, extrap_freqs,
             num_freqs*sizeof(double), gpuMemcpyHostToDevice ));

  source_t *d_chunked_source = NULL;
  components_t d_components;

  if (comptype == POINT) {

    chunked_source->point_components = components;
    chunked_source->n_points = n_powers + n_curves + n_lists;
    chunked_source->n_point_powers = n_powers;
    chunked_source->n_point_curves = n_curves;
    chunked_source->n_point_lists = n_lists;

    chunked_source->n_gauss = 0;
    chunked_source->n_gauss_lists = 0;
    chunked_source->n_gauss_powers = 0;
    chunked_source->n_gauss_curves = 0;
    chunked_source->n_shapes = 0;
    chunked_source->n_shape_lists = 0;
    chunked_source->n_shape_powers = 0;
    chunked_source->n_shape_curves = 0;
    chunked_source->n_shape_coeffs = 0;

    d_chunked_source = copy_chunked_source_to_GPU(chunked_source);
    malloc_lmn_arrays(d_chunked_source, &components, num_components, comptype);

    malloc_extrapolated_flux_arrays(&d_chunked_source->point_components,
                                    d_chunked_source->n_points,
                                    num_freqs, do_QUV);
    extrapolate_Stokes(d_chunked_source, d_extrap_freqs, num_freqs, POINT, do_QUV);
    d_components = d_chunked_source->point_components;
  }
  else if (comptype == GAUSSIAN) {

    chunked_source->gauss_components = components;
    chunked_source->n_gauss = n_powers + n_curves + n_lists;
    chunked_source->n_gauss_powers = n_powers;
    chunked_source->n_gauss_curves = n_curves;
    chunked_source->n_gauss_lists = n_lists;

    chunked_source->n_points = 0;
    chunked_source->n_point_lists = 0;
    chunked_source->n_point_powers = 0;
    chunked_source->n_point_curves = 0;
    chunked_source->n_shapes = 0;
    chunked_source->n_shape_lists = 0;
    chunked_source->n_shape_powers = 0;
    chunked_source->n_shape_curves = 0;
    chunked_source->n_shape_coeffs = 0;

    d_chunked_source = copy_chunked_source_to_GPU(chunked_source);
    malloc_lmn_arrays(d_chunked_source, &components, num_components, comptype);

    malloc_extrapolated_flux_arrays(&d_chunked_source->gauss_components,
                                    d_chunked_source->n_gauss,
                                    num_freqs, do_QUV);
    extrapolate_Stokes(d_chunked_source, d_extrap_freqs, num_freqs, GAUSSIAN, do_QUV);
    d_components = d_chunked_source->gauss_components;
  }
  else if (comptype == SHAPELET) {
    chunked_source->shape_components = components;
    chunked_source->n_shapes = n_powers + n_curves + n_lists;
    chunked_source->n_shape_powers = n_powers;
    chunked_source->n_shape_curves = n_curves;
    chunked_source->n_shape_lists = n_lists;
    chunked_source->n_shape_coeffs = num_shape_coeffs;

    chunked_source->n_points = 0;
    chunked_source->n_point_lists = 0;
    chunked_source->n_point_powers = 0;
    chunked_source->n_point_curves = 0;
    chunked_source->n_gauss = 0;
    chunked_source->n_gauss_lists = 0;
    chunked_source->n_gauss_powers = 0;
    chunked_source->n_gauss_curves = 0;

    d_chunked_source = copy_chunked_source_to_GPU(chunked_source);
    malloc_lmn_arrays(d_chunked_source, &components, num_components, comptype);

    malloc_extrapolated_flux_arrays(&d_chunked_source->shape_components,
                                    d_chunked_source->n_shapes,
                                    num_freqs, do_QUV);
    extrapolate_Stokes(d_chunked_source, d_extrap_freqs, num_freqs, SHAPELET, do_QUV);
    d_components = d_chunked_source->shape_components;
  }

  //Something to store the primary beam gains (all 4 pols) in
  d_beam_gains_t d_beam_gains;
  int num_beam_values = num_components*num_freqs*num_times;

  ( gpuMalloc( (void**)&d_beam_gains.d_gxs,
                                      num_beam_values*sizeof(gpuUserComplex) ));
  ( gpuMalloc( (void**)&d_beam_gains.d_Dxs,
                                      num_beam_values*sizeof(gpuUserComplex) ));
  ( gpuMalloc( (void**)&d_beam_gains.d_Dys,
                                      num_beam_values*sizeof(gpuUserComplex) ));
  ( gpuMalloc( (void**)&d_beam_gains.d_gys,
                                      num_beam_values*sizeof(gpuUserComplex) ));

  ( gpuMemcpy(d_beam_gains.d_gxs, (gpuUserComplex *)gxs,
              num_beam_values*sizeof(gpuUserComplex), gpuMemcpyHostToDevice ));
  ( gpuMemcpy(d_beam_gains.d_Dxs, (gpuUserComplex *)Dxs,
              num_beam_values*sizeof(gpuUserComplex), gpuMemcpyHostToDevice ));
  ( gpuMemcpy(d_beam_gains.d_Dys, (gpuUserComplex *)Dys,
              num_beam_values*sizeof(gpuUserComplex), gpuMemcpyHostToDevice ));
  ( gpuMemcpy(d_beam_gains.d_gys, (gpuUserComplex *)gys,
              num_beam_values*sizeof(gpuUserComplex), gpuMemcpyHostToDevice ));

  user_precision_t *d_sum_visi_XX_real = NULL;
  user_precision_t *d_sum_visi_XY_real = NULL;
  user_precision_t *d_sum_visi_YX_real = NULL;
  user_precision_t *d_sum_visi_YY_real = NULL;
  user_precision_t *d_sum_visi_XX_imag = NULL;
  user_precision_t *d_sum_visi_XY_imag = NULL;
  user_precision_t *d_sum_visi_YX_imag = NULL;
  user_precision_t *d_sum_visi_YY_imag = NULL;

  ( gpuMalloc( (void**)&d_sum_visi_XX_real,
                                          num_cross*sizeof(user_precision_t) ));
  ( gpuMalloc( (void**)&d_sum_visi_XY_real,
                                          num_cross*sizeof(user_precision_t) ));
  ( gpuMalloc( (void**)&d_sum_visi_YX_real,
                                          num_cross*sizeof(user_precision_t) ));
  ( gpuMalloc( (void**)&d_sum_visi_YY_real,
                                          num_cross*sizeof(user_precision_t) ));
  ( gpuMalloc( (void**)&d_sum_visi_XX_imag,
                                          num_cross*sizeof(user_precision_t) ));
  ( gpuMalloc( (void**)&d_sum_visi_XY_imag,
                                          num_cross*sizeof(user_precision_t) ));
  ( gpuMalloc( (void**)&d_sum_visi_YX_imag,
                                          num_cross*sizeof(user_precision_t) ));
  ( gpuMalloc( (void**)&d_sum_visi_YY_imag,
                                          num_cross*sizeof(user_precision_t) ));

  //Make sure the visis start at zero by copying across host versions, which
  //should be set to zero already
  ( gpuMemcpy( d_sum_visi_XX_real, sum_visi_XX_real,
    num_cross*sizeof(user_precision_t), gpuMemcpyHostToDevice ) );
  ( gpuMemcpy( d_sum_visi_XY_real, sum_visi_XY_real,
    num_cross*sizeof(user_precision_t), gpuMemcpyHostToDevice ) );
  ( gpuMemcpy( d_sum_visi_YX_real, sum_visi_YX_real,
    num_cross*sizeof(user_precision_t), gpuMemcpyHostToDevice ) );
  ( gpuMemcpy( d_sum_visi_YY_real, sum_visi_YY_real,
    num_cross*sizeof(user_precision_t), gpuMemcpyHostToDevice ) );
  ( gpuMemcpy( d_sum_visi_XX_imag, sum_visi_XX_imag,
    num_cross*sizeof(user_precision_t), gpuMemcpyHostToDevice ) );
  ( gpuMemcpy( d_sum_visi_XY_imag, sum_visi_XY_imag,
    num_cross*sizeof(user_precision_t), gpuMemcpyHostToDevice ) );
  ( gpuMemcpy( d_sum_visi_YX_imag, sum_visi_YX_imag,
    num_cross*sizeof(user_precision_t), gpuMemcpyHostToDevice ) );
  ( gpuMemcpy( d_sum_visi_YY_imag, sum_visi_YY_imag,
    num_cross*sizeof(user_precision_t), gpuMemcpyHostToDevice ) );

  dim3 grid, threads;

  threads.x = 128;
  grid.x = (int)ceil( (float)num_cross / (float)threads.x );

  //Shapelets need many many extra things

  user_precision_t *d_sbf=NULL;
  user_precision_t *d_u_shapes = NULL;
  user_precision_t *d_v_shapes = NULL;

  if (comptype == SHAPELET) {

    ( gpuMalloc( (void**)&d_u_shapes,
             num_components*num_baselines*num_times*sizeof(user_precision_t) ) );
    ( gpuMalloc( (void**)&d_v_shapes,
             num_components*num_baselines*num_times*sizeof(user_precision_t) ) );

    ( gpuMemcpy(d_u_shapes, u_shapes,
                 num_components*num_baselines*num_times*sizeof(user_precision_t),
                                                       gpuMemcpyHostToDevice ));
    ( gpuMemcpy(d_v_shapes, v_shapes,
                 num_components*num_baselines*num_times*sizeof(user_precision_t),
                                                       gpuMemcpyHostToDevice ));

    ( gpuMalloc( (void**)&d_components.shape_coeffs,
                                                num_shape_coeffs*sizeof(user_precision_t) ));
    ( gpuMalloc( (void**)&d_components.n1s,
                                                num_shape_coeffs*sizeof(user_precision_t) ));
    ( gpuMalloc( (void**)&d_components.n2s,
                                                num_shape_coeffs*sizeof(user_precision_t) ));
    ( gpuMalloc( (void**)&d_components.param_indexes,
                                                num_shape_coeffs*sizeof(user_precision_t) ));

    ( gpuMemcpy(d_components.shape_coeffs,
                          components.shape_coeffs, num_shape_coeffs*sizeof(user_precision_t),
                          gpuMemcpyHostToDevice ));
    ( gpuMemcpy(d_components.n1s,
                          components.n1s, num_shape_coeffs*sizeof(user_precision_t),
                          gpuMemcpyHostToDevice ));
    ( gpuMemcpy(d_components.n2s,
                          components.n2s, num_shape_coeffs*sizeof(user_precision_t),
                          gpuMemcpyHostToDevice ));
    ( gpuMemcpy(d_components.param_indexes,
                          components.param_indexes, num_shape_coeffs*sizeof(user_precision_t),
                          gpuMemcpyHostToDevice ));
    ( gpuMalloc( (void**)&(d_sbf), sbf_N*sbf_L*sizeof(user_precision_t) ));
    ( gpuMemcpy( d_sbf, sbf, sbf_N*sbf_L*sizeof(user_precision_t),
                        gpuMemcpyHostToDevice ));
  }

  if (comptype == POINT || comptype == GAUSSIAN ) {

    gpuErrorCheckKernel("kern_calc_visi_point_or_gauss",
                  kern_calc_visi_point_or_gauss, grid, threads,
                  d_components, d_beam_gains,
                  d_us, d_vs, d_ws,
                  d_sum_visi_XX_real, d_sum_visi_XX_imag,
                  d_sum_visi_XY_real, d_sum_visi_XY_imag,
                  d_sum_visi_YX_real, d_sum_visi_YX_imag,
                  d_sum_visi_YY_real, d_sum_visi_YY_imag,
                  num_components, num_baselines, num_freqs, num_cross,
                  num_times, beamtype, comptype, do_QUV);
  }
  else if (comptype == SHAPELET) {
    gpuErrorCheckKernel("kern_calc_visi_shapelets",
                  kern_calc_visi_shapelets, grid, threads,
                  d_components, d_beam_gains,
                  d_us, d_vs, d_ws,
                  d_allsteps_wavelengths,
                  d_u_shapes, d_v_shapes,
                  d_sum_visi_XX_real, d_sum_visi_XX_imag,
                  d_sum_visi_XY_real, d_sum_visi_XY_imag,
                  d_sum_visi_YX_real, d_sum_visi_YX_imag,
                  d_sum_visi_YY_real, d_sum_visi_YY_imag,
                  d_sbf,  num_components,
                  num_baselines, num_freqs, num_cross,
                  num_shape_coeffs, num_times, beamtype, do_QUV);
  }

  ( gpuMemcpy(sum_visi_XX_real, d_sum_visi_XX_real,
                             num_cross*sizeof(user_precision_t), gpuMemcpyDeviceToHost ));
  ( gpuMemcpy(sum_visi_XY_real, d_sum_visi_XY_real,
                             num_cross*sizeof(user_precision_t), gpuMemcpyDeviceToHost ));
  ( gpuMemcpy(sum_visi_YX_real, d_sum_visi_YX_real,
                             num_cross*sizeof(user_precision_t), gpuMemcpyDeviceToHost ));
  ( gpuMemcpy(sum_visi_YY_real, d_sum_visi_YY_real,
                             num_cross*sizeof(user_precision_t), gpuMemcpyDeviceToHost ));
  ( gpuMemcpy(sum_visi_XX_imag, d_sum_visi_XX_imag,
                             num_cross*sizeof(user_precision_t), gpuMemcpyDeviceToHost ));
  ( gpuMemcpy(sum_visi_XY_imag, d_sum_visi_XY_imag,
                             num_cross*sizeof(user_precision_t), gpuMemcpyDeviceToHost ));
  ( gpuMemcpy(sum_visi_YX_imag, d_sum_visi_YX_imag,
                             num_cross*sizeof(user_precision_t), gpuMemcpyDeviceToHost ));
  ( gpuMemcpy(sum_visi_YY_imag, d_sum_visi_YY_imag,
                             num_cross*sizeof(user_precision_t), gpuMemcpyDeviceToHost ));


  free_d_components(d_chunked_source, comptype);

  (  gpuFree( d_sum_visi_XX_real ) );
  (  gpuFree( d_sum_visi_XX_imag ) );
  (  gpuFree( d_sum_visi_XY_real ) );
  (  gpuFree( d_sum_visi_XY_imag ) );
  (  gpuFree( d_sum_visi_YX_real ) );
  (  gpuFree( d_sum_visi_YX_imag ) );
  (  gpuFree( d_sum_visi_YY_real ) );
  (  gpuFree( d_sum_visi_YY_imag ) );
  (  gpuFree( d_allsteps_wavelengths ) );



  free_beam_gains(d_beam_gains, beamtype);

  (  gpuFree( d_us ) );
  (  gpuFree( d_vs ) );
  (  gpuFree( d_ws ) );

  ( gpuFree( d_extrap_freqs ) );

  if (comptype == POINT) {
    free_extrapolated_flux_arrays(&d_chunked_source->point_components, do_QUV);
  }
  else if (comptype == GAUSSIAN) {
    free_extrapolated_flux_arrays(&d_chunked_source->gauss_components, do_QUV);
  }
  if (comptype == SHAPELET){
    free_extrapolated_flux_arrays(&d_chunked_source->shape_components, do_QUV);
    (  gpuFree( d_sbf) );
    (  gpuFree( d_u_shapes) );
    (  gpuFree( d_v_shapes) );
  }
}


extern "C" void test_kern_calc_autos(components_t *components, int beamtype,
                                     int num_components, int num_baselines,
                                     int num_freqs, int num_times, int num_ants,
                                     int num_beams,
                                     visibility_set_t *visibility_set){

  int use_twobeams = 0;
  if (num_beams > 1) {
    use_twobeams = 1;
  }

  ////malloc on device and copy extrapolated fluxes
  int num_pb_values = num_beams*num_freqs*num_times*num_components;

  int num_autos = num_ants*num_freqs*num_times;
  int num_cross = num_baselines*num_freqs*num_times;
  int num_visis = num_cross + num_autos;

  components_t *d_components = (components_t* )malloc(sizeof(components_t));
  // components_t d_components;

  ( gpuMalloc( (void**)&d_components->extrap_stokesI,
              num_components*num_freqs*sizeof(user_precision_t) ));
  ( gpuMalloc( (void**)&d_components->extrap_stokesQ,
              num_components*num_freqs*sizeof(user_precision_t) ));
  ( gpuMalloc( (void**)&d_components->extrap_stokesU,
              num_components*num_freqs*sizeof(user_precision_t) ));
  ( gpuMalloc( (void**)&d_components->extrap_stokesV,
              num_components*num_freqs*sizeof(user_precision_t) ));


  ( gpuMemcpy(d_components->extrap_stokesI,
         components->extrap_stokesI, num_components*num_freqs*sizeof(user_precision_t), gpuMemcpyHostToDevice ) );
  ( gpuMemcpy(d_components->extrap_stokesQ,
         components->extrap_stokesQ, num_components*num_freqs*sizeof(user_precision_t), gpuMemcpyHostToDevice ) );
  ( gpuMemcpy(d_components->extrap_stokesU,
         components->extrap_stokesU, num_components*num_freqs*sizeof(user_precision_t), gpuMemcpyHostToDevice ) );
  ( gpuMemcpy(d_components->extrap_stokesV,
         components->extrap_stokesV, num_components*num_freqs*sizeof(user_precision_t), gpuMemcpyHostToDevice ) );
  //
  // //malloc on device and copy beam values
  //
  d_beam_gains_t *d_component_beam_gains = (d_beam_gains_t* )malloc(sizeof(d_beam_gains_t));
  // d_beam_gains_t d_component_beam_gains;
  ( gpuMalloc( (void**)&d_component_beam_gains->d_gxs,
                                        num_pb_values*sizeof(gpuUserComplex) ));
  ( gpuMalloc( (void**)&d_component_beam_gains->d_Dxs,
                                        num_pb_values*sizeof(gpuUserComplex) ));
  ( gpuMalloc( (void**)&d_component_beam_gains->d_Dys,
                                        num_pb_values*sizeof(gpuUserComplex) ));
  ( gpuMalloc( (void**)&d_component_beam_gains->d_gys,
                                        num_pb_values*sizeof(gpuUserComplex) ));

  ( gpuMemcpy(d_component_beam_gains->d_gxs,
          (gpuUserComplex* )components->gxs, num_pb_values*sizeof(gpuUserComplex),
                                                     gpuMemcpyHostToDevice ) );
  ( gpuMemcpy(d_component_beam_gains->d_Dxs,
          (gpuUserComplex* )components->Dxs, num_pb_values*sizeof(gpuUserComplex),
                                                     gpuMemcpyHostToDevice ) );
  ( gpuMemcpy(d_component_beam_gains->d_Dys,
          (gpuUserComplex* )components->Dys, num_pb_values*sizeof(gpuUserComplex),
                                                     gpuMemcpyHostToDevice ) );
  ( gpuMemcpy(d_component_beam_gains->d_gys,
          (gpuUserComplex* )components->gys, num_pb_values*sizeof(gpuUserComplex),
                                                     gpuMemcpyHostToDevice ) );


  user_precision_t *d_sum_visi_XX_real;
  user_precision_t *d_sum_visi_XX_imag;
  user_precision_t *d_sum_visi_XY_real;
  user_precision_t *d_sum_visi_XY_imag;
  user_precision_t *d_sum_visi_YX_real;
  user_precision_t *d_sum_visi_YX_imag;
  user_precision_t *d_sum_visi_YY_real;
  user_precision_t *d_sum_visi_YY_imag;

  ( gpuMalloc( (void**)&d_sum_visi_XX_real,
                      num_visis*sizeof(user_precision_t) ) );
  ( gpuMalloc( (void**)&d_sum_visi_XX_imag,
                      num_visis*sizeof(user_precision_t) ) );
  ( gpuMalloc( (void**)&d_sum_visi_XY_real,
                      num_visis*sizeof(user_precision_t) ) );
  ( gpuMalloc( (void**)&d_sum_visi_XY_imag,
                      num_visis*sizeof(user_precision_t) ) );
  ( gpuMalloc( (void**)&d_sum_visi_YX_real,
                      num_visis*sizeof(user_precision_t) ) );
  ( gpuMalloc( (void**)&d_sum_visi_YX_imag,
                      num_visis*sizeof(user_precision_t) ) );
  ( gpuMalloc( (void**)&d_sum_visi_YY_real,
                      num_visis*sizeof(user_precision_t) ) );
  ( gpuMalloc( (void**)&d_sum_visi_YY_imag,
                      num_visis*sizeof(user_precision_t) ) );


  //ensure d_sum_visi_XX_real are set entirely to zero by copying the host
  //array values, which have been set explictly to zero during chunking
  ( gpuMemcpy(d_sum_visi_XX_real,
             visibility_set->sum_visi_XX_real,
             num_visis*sizeof(user_precision_t), gpuMemcpyHostToDevice ) );
  ( gpuMemcpy(d_sum_visi_XX_imag,
             visibility_set->sum_visi_XX_imag,
             num_visis*sizeof(user_precision_t), gpuMemcpyHostToDevice ) );
  ( gpuMemcpy(d_sum_visi_XY_real,
             visibility_set->sum_visi_XY_real,
             num_visis*sizeof(user_precision_t), gpuMemcpyHostToDevice ) );
  ( gpuMemcpy(d_sum_visi_XY_imag,
             visibility_set->sum_visi_XY_imag,
             num_visis*sizeof(user_precision_t), gpuMemcpyHostToDevice ) );
  ( gpuMemcpy(d_sum_visi_YX_real,
             visibility_set->sum_visi_YX_real,
             num_visis*sizeof(user_precision_t), gpuMemcpyHostToDevice ) );
  ( gpuMemcpy(d_sum_visi_YX_imag,
             visibility_set->sum_visi_YX_imag,
             num_visis*sizeof(user_precision_t), gpuMemcpyHostToDevice ) );
  ( gpuMemcpy(d_sum_visi_YY_real,
             visibility_set->sum_visi_YY_real,
             num_visis*sizeof(user_precision_t), gpuMemcpyHostToDevice ) );
  ( gpuMemcpy(d_sum_visi_YY_imag,
             visibility_set->sum_visi_YY_imag,
             num_visis*sizeof(user_precision_t), gpuMemcpyHostToDevice ) );

  dim3 grid, threads;

  threads.x = 64;
  threads.y = 2;
  threads.z = 1;
  grid.x = (int)ceil( (float)(num_freqs*num_times) / (float)threads.x );
  grid.y = (int)ceil( (float)(num_ants) / (float)threads.y );
  grid.z = 1;

  int do_QUV = 0;

  int *d_ant_to_auto_map = NULL;

  if (use_twobeams == 1) {

    int *ant_to_auto_map = NULL;
    ant_to_auto_map = (int *)malloc(num_ants*sizeof(int));
    for (int ant = 0; ant < num_ants; ant++){
        ant_to_auto_map[ant] = ant;
    }

    ( gpuMalloc( (void**)&d_ant_to_auto_map,
                                    num_ants*sizeof(int) ));
    ( gpuMemcpy(d_ant_to_auto_map, ant_to_auto_map,
                                    num_ants*sizeof(int), gpuMemcpyHostToDevice ));

    free(ant_to_auto_map);
  }

  gpuErrorCheckKernel("kern_calc_autos",
                kern_calc_autos, grid, threads,
                *d_components, *d_component_beam_gains,
                beamtype, num_components, num_baselines,
                num_freqs, num_times, num_ants,
                d_sum_visi_XX_real, d_sum_visi_XX_imag,
                d_sum_visi_XY_real, d_sum_visi_XY_imag,
                d_sum_visi_YX_real, d_sum_visi_YX_imag,
                d_sum_visi_YY_real, d_sum_visi_YY_imag,
                do_QUV, use_twobeams, d_ant_to_auto_map,
                d_ant_to_auto_map);

  //Copy outputs onto host so we can check our answers
  ( gpuMemcpy(visibility_set->sum_visi_XX_real,
                         d_sum_visi_XX_real, num_visis*sizeof(user_precision_t),
                                                     gpuMemcpyDeviceToHost ));
  ( gpuMemcpy(visibility_set->sum_visi_XY_real,
                         d_sum_visi_XY_real, num_visis*sizeof(user_precision_t),
                                                     gpuMemcpyDeviceToHost ));
  ( gpuMemcpy(visibility_set->sum_visi_YX_real,
                         d_sum_visi_YX_real, num_visis*sizeof(user_precision_t),
                                                     gpuMemcpyDeviceToHost ));
  ( gpuMemcpy(visibility_set->sum_visi_YY_real,
                         d_sum_visi_YY_real, num_visis*sizeof(user_precision_t),
                                                     gpuMemcpyDeviceToHost ));
  ( gpuMemcpy(visibility_set->sum_visi_XX_imag,
                         d_sum_visi_XX_imag, num_visis*sizeof(user_precision_t),
                                                     gpuMemcpyDeviceToHost ));
  ( gpuMemcpy(visibility_set->sum_visi_XY_imag,
                         d_sum_visi_XY_imag, num_visis*sizeof(user_precision_t),
                                                     gpuMemcpyDeviceToHost ));
  ( gpuMemcpy(visibility_set->sum_visi_YX_imag,
                         d_sum_visi_YX_imag, num_visis*sizeof(user_precision_t),
                                                     gpuMemcpyDeviceToHost ));
  ( gpuMemcpy(visibility_set->sum_visi_YY_imag,
                         d_sum_visi_YY_imag, num_visis*sizeof(user_precision_t),
                                                     gpuMemcpyDeviceToHost ));


  (  gpuFree( d_components->extrap_stokesI ) );
  (  gpuFree( d_components->extrap_stokesQ ) );
  (  gpuFree( d_components->extrap_stokesU ) );
  (  gpuFree( d_components->extrap_stokesV ) );
  (  gpuFree( d_component_beam_gains->d_gxs ) );
  (  gpuFree( d_component_beam_gains->d_Dxs ) );
  (  gpuFree( d_component_beam_gains->d_Dys ) );
  (  gpuFree( d_component_beam_gains->d_gys ) );


  (  gpuFree( d_sum_visi_XX_real ) );
  (  gpuFree( d_sum_visi_XX_imag ) );
  (  gpuFree( d_sum_visi_XY_real ) );
  (  gpuFree( d_sum_visi_XY_imag ) );
  (  gpuFree( d_sum_visi_YX_real ) );
  (  gpuFree( d_sum_visi_YX_imag ) );
  (  gpuFree( d_sum_visi_YY_real ) );
  (  gpuFree( d_sum_visi_YY_imag ) );

  if (use_twobeams == 1) {
    (  gpuFree( d_ant_to_auto_map ) );
  }
}