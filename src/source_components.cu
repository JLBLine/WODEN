#include <stdio.h>
#include <stdlib.h>
#include <cuComplex.h>
#include <complex.h>
#include <math.h>
#include "cudacomplex.h"
#include "fundamental_coords.h"
#include "constants.h"
#include "shapelet_basis.h"
#include "source_components.h"
#include "cudacheck.h"
#include "woden_struct_defs.h"
#include "primary_beam_cuda.h"
#include "woden_precision_defs.h"

__device__ void extrap_stokes(user_precision_t *d_allsteps_wavelengths,
           double *d_ref_freqs,
           user_precision_t *d_ref_stokesI, user_precision_t *d_ref_stokesQ,
           user_precision_t *d_ref_stokesU, user_precision_t *d_ref_stokesV,
           user_precision_t *d_SIs, int iComponent, int iBaseline,
           user_precision_t * flux_I, user_precision_t * flux_Q,
           user_precision_t * flux_U, user_precision_t * flux_V){

  double d_freq = VELC / d_allsteps_wavelengths[iBaseline];
  double d_ref_freq = d_ref_freqs[iComponent];

  user_precision_t flux_ratio = pow(d_freq / d_ref_freq, d_SIs[iComponent]);

  * flux_I = d_ref_stokesI[iComponent] * flux_ratio;
  * flux_Q = d_ref_stokesQ[iComponent] * flux_ratio;
  * flux_U = d_ref_stokesU[iComponent] * flux_ratio;
  * flux_V = d_ref_stokesV[iComponent] * flux_ratio;

}

__device__  cuUserComplex calc_measurement_equation(user_precision_t *d_us,
           user_precision_t *d_vs, user_precision_t *d_ws,
           double *d_ls, double *d_ms, double *d_ns,
           const int iBaseline, const int iComponent){

  cuUserComplex visi;

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

__device__ void apply_beam_gains(cuUserComplex g1x, cuUserComplex D1x,
          cuUserComplex D1y, cuUserComplex g1y,
          cuUserComplex g2x, cuUserComplex D2x,
          cuUserComplex D2y, cuUserComplex g2y,
          user_precision_t flux_I, user_precision_t flux_Q,
          user_precision_t flux_U, user_precision_t flux_V,
          cuUserComplex visi_component,
          cuUserComplex * visi_XX, cuUserComplex * visi_XY,
          cuUserComplex * visi_YX, cuUserComplex * visi_YY) {

  //Conjugate the second beam gains
  cuUserComplex g2x_conj = make_cuUserComplex(g2x.x,-g2x.y);
  cuUserComplex D2x_conj = make_cuUserComplex(D2x.x,-D2x.y);
  cuUserComplex D2y_conj = make_cuUserComplex(D2y.x,-D2y.y);
  cuUserComplex g2y_conj = make_cuUserComplex(g2y.x,-g2y.y);

  //Create the Stokes visibilities
  cuUserComplex visi_I = make_cuUserComplex(flux_I, 0.0)*visi_component;
  cuUserComplex visi_Q = make_cuUserComplex(flux_Q, 0.0)*visi_component;
  cuUserComplex visi_U = make_cuUserComplex(flux_U, 0.0)*visi_component;
  cuUserComplex visi_V = make_cuUserComplex(flux_V, 0.0)*visi_component;

  cuUserComplex this_XX;
  cuUserComplex this_XY;
  cuUserComplex this_YX;
  cuUserComplex this_YY;

  // this_XX = (g1x*g2x_conj + D1x*D2x_conj);
  // this_XY = (g1x*D2y_conj + D1x*g2y_conj);
  // this_YX = (D1y*g2x_conj + g1y*D2x_conj);
  // this_YY = (D1y*D2y_conj + g1y*g2y_conj);
  //
  // printf("XX %.16f %.16f\n",this_XX.x, this_XX.y );
  // printf("XY %.16f %.16f\n",this_XY.x, this_XY.y );
  // printf("YX %.16f %.16f\n",this_YX.x, this_YX.y );
  // printf("YY %.16f %.16f\n",this_YY.x, this_YY.y );

  this_XX = (g1x*g2x_conj + D1x*D2x_conj)*visi_I;
  this_XX += (g1x*g2x_conj - D1x*D2x_conj)*visi_Q;
  this_XX += (g1x*D2x_conj + D1x*g2x_conj)*visi_U;
  this_XX += (make_cuUserComplex(0.0,1.0)*visi_V)*(g1x*D2x_conj - D1x*g2x_conj);

  this_XY = (g1x*D2y_conj + D1x*g2y_conj)*visi_I;
  this_XY += (g1x*D2y_conj - D1x*g2y_conj)*visi_Q;
  this_XY += (g1x*g2y_conj + D1x*D2y_conj)*visi_U;
  this_XY += (make_cuUserComplex(0.0,1.0)*visi_V)* (g1x*g2y_conj - D1x*D2y_conj);

  this_YX = (D1y*g2x_conj + g1y*D2x_conj)*visi_I;
  this_YX += (D1y*g2x_conj - g1y*D2x_conj)*visi_Q;
  this_YX += (D1y*D2x_conj + g1y*g2x_conj)*visi_U;
  this_YX += (make_cuUserComplex(0.0,1.0)*visi_V)* (D1y*D2x_conj - g1y*g2x_conj);

  this_YY = (D1y*D2y_conj + g1y*g2y_conj)*visi_I;
  this_YY += (D1y*D2y_conj - g1y*g2y_conj)*visi_Q;
  this_YY += (D1y*g2y_conj + g1y*D2y_conj)*visi_U;
  this_YY += (make_cuUserComplex(0.0,1.0)*visi_V)* (D1y*g2y_conj - g1y*D2y_conj);

  * visi_XX = this_XX;
  * visi_XY = this_XY;
  * visi_YX = this_YX;
  * visi_YY = this_YY;

}

__device__ void get_beam_gains(int iBaseline, int iComponent, int num_freqs,
           int num_baselines, int num_components, int num_times, int beamtype,
           cuUserComplex *d_primay_beam_J00, cuUserComplex *d_primay_beam_J01,
           cuUserComplex *d_primay_beam_J10, cuUserComplex *d_primay_beam_J11,
           cuUserComplex * g1x, cuUserComplex * D1x,
           cuUserComplex * D1y, cuUserComplex * g1y,
           cuUserComplex * g2x, cuUserComplex * D2x,
           cuUserComplex * D2y, cuUserComplex * g2y){

  int beam_ind = 0;
  int time_ind = 0;
  int freq_ind = 0;

  time_ind = (int)floorf( (user_precision_t)iBaseline / ((user_precision_t)num_baselines * (user_precision_t)num_freqs));
  freq_ind = (int)floorf( ((user_precision_t)iBaseline - ((user_precision_t)time_ind*(user_precision_t)num_baselines * (user_precision_t)num_freqs)) / (user_precision_t)num_baselines);
  beam_ind = num_freqs*time_ind*num_components + (num_components*freq_ind) + iComponent;

    //Set gains to one if no beam
  if (beamtype == NO_BEAM) {
    * g1x = make_cuUserComplex(1.0, 0.0);
    * g2x = make_cuUserComplex(1.0, 0.0);
    * g1y = make_cuUserComplex(1.0, 0.0);
    * g2y = make_cuUserComplex(1.0, 0.0);
  }

  //Get gains if using a beam
  else {
    * g1x = d_primay_beam_J00[beam_ind];
    * g2x = d_primay_beam_J00[beam_ind];
    * g1y = d_primay_beam_J11[beam_ind];
    * g2y = d_primay_beam_J11[beam_ind];

  }

  //Only MWA models have leakge terms at the moment
  if (beamtype == FEE_BEAM || beamtype == FEE_BEAM_INTERP || beamtype == MWA_ANALY) {
    * D1x = d_primay_beam_J01[beam_ind];
    * D2x = d_primay_beam_J01[beam_ind];
    * D1y = d_primay_beam_J10[beam_ind];
    * D2y = d_primay_beam_J10[beam_ind];
  }
  // Set leakage to zero if no leakage
  else {
    * D1x = make_cuUserComplex(0.0, 0.0);
    * D2x = make_cuUserComplex(0.0, 0.0);
    * D1y = make_cuUserComplex(0.0, 0.0);
    * D2y = make_cuUserComplex(0.0, 0.0);
  }
} //end __device__ get_beam_gains

__device__ void update_sum_visis(int iBaseline, int iComponent, int num_freqs,
    int num_baselines, int num_components, int num_times, int beamtype,
    cuUserComplex *d_primay_beam_J00, cuUserComplex *d_primay_beam_J01,
    cuUserComplex *d_primay_beam_J10, cuUserComplex *d_primay_beam_J11,
    cuUserComplex visi_component,
    user_precision_t flux_I, user_precision_t flux_Q,
    user_precision_t flux_U, user_precision_t flux_V,
    user_precision_t *d_sum_visi_XX_real, user_precision_t *d_sum_visi_XX_imag,
    user_precision_t *d_sum_visi_XY_real, user_precision_t *d_sum_visi_XY_imag,
    user_precision_t *d_sum_visi_YX_real, user_precision_t *d_sum_visi_YX_imag,
    user_precision_t *d_sum_visi_YY_real, user_precision_t *d_sum_visi_YY_imag){

    cuUserComplex g1x;
    cuUserComplex D1x;
    cuUserComplex D1y;
    cuUserComplex g1y;
    cuUserComplex g2x;
    cuUserComplex D2x;
    cuUserComplex D2y;
    cuUserComplex g2y;

    get_beam_gains(iBaseline, iComponent, num_freqs,
               num_baselines, num_components, num_times, beamtype,
               d_primay_beam_J00, d_primay_beam_J01,
               d_primay_beam_J10, d_primay_beam_J11,
               &g1x, &D1x, &D1y, &g1y, &g2x, &D2x, &D2y, &g2y);

    cuUserComplex visi_XX;
    cuUserComplex visi_XY;
    cuUserComplex visi_YX;
    cuUserComplex visi_YY;

    apply_beam_gains(g1x, D1x, D1y, g1y, g2x, D2x, D2y, g2y,
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

extern "C" void source_component_common(int num_components,
           int num_shape_coeffs, components_t *components,
           double *d_freqs, woden_settings_t *woden_settings,
           beam_settings_t *beam_settings, e_component_type comptype,
           components_t *d_components, d_beam_gains_t *d_component_beam_gains){

  // //Copy across many things from host to Device
  cudaErrorCheckCall( cudaMalloc( (void**)&d_components->ras,
                      num_components*sizeof(double) ) );
  cudaErrorCheckCall( cudaMemcpy( d_components->ras, components->ras,
                      num_components*sizeof(double), cudaMemcpyHostToDevice ) );

  cudaErrorCheckCall( cudaMalloc( (void**)&d_components->decs,
                      num_components*sizeof(double) ) );
  cudaErrorCheckCall( cudaMemcpy( d_components->decs, components->decs,
                      num_components*sizeof(double), cudaMemcpyHostToDevice ) );

  cudaErrorCheckCall( cudaMalloc( (void**)&d_components->ref_freqs,
                      num_components*sizeof(double) ) );
  cudaErrorCheckCall( cudaMemcpy( d_components->ref_freqs, components->ref_freqs,
                      num_components*sizeof(double), cudaMemcpyHostToDevice ) );

  cudaErrorCheckCall( cudaMalloc( (void**)&d_components->ref_stokesI,
                      num_components*sizeof(user_precision_t) ) );
  cudaErrorCheckCall( cudaMemcpy( d_components->ref_stokesI, components->ref_stokesI,
                      num_components*sizeof(user_precision_t), cudaMemcpyHostToDevice ) );

  cudaErrorCheckCall( cudaMalloc( (void**)&d_components->ref_stokesQ,
                      num_components*sizeof(user_precision_t) ) );
  cudaErrorCheckCall( cudaMemcpy( d_components->ref_stokesQ, components->ref_stokesQ,
                      num_components*sizeof(user_precision_t), cudaMemcpyHostToDevice ) );

  cudaErrorCheckCall( cudaMalloc( (void**)&d_components->ref_stokesU,
                      num_components*sizeof(user_precision_t) ) );
  cudaErrorCheckCall( cudaMemcpy( d_components->ref_stokesU, components->ref_stokesU,
                      num_components*sizeof(user_precision_t), cudaMemcpyHostToDevice ) );

  cudaErrorCheckCall( cudaMalloc( (void**)&d_components->ref_stokesV,
                      num_components*sizeof(user_precision_t) ) );
  cudaErrorCheckCall( cudaMemcpy( d_components->ref_stokesV, components->ref_stokesV,
                      num_components*sizeof(user_precision_t), cudaMemcpyHostToDevice ) );

  cudaErrorCheckCall( cudaMalloc( (void**)&d_components->SIs,
                      num_components*sizeof(user_precision_t) ) );
  cudaErrorCheckCall( cudaMemcpy( d_components->SIs, components->SIs,
                      num_components*sizeof(user_precision_t), cudaMemcpyHostToDevice ) );

  if (comptype == GAUSSIAN || comptype == SHAPELET ) {
    cudaErrorCheckCall( cudaMalloc( (void**)&d_components->pas,
                        num_components*sizeof(user_precision_t) ) );
    cudaErrorCheckCall( cudaMemcpy( d_components->pas, components->pas,
                        num_components*sizeof(user_precision_t), cudaMemcpyHostToDevice ) );

    cudaErrorCheckCall( cudaMalloc( (void**)&d_components->majors,
                        num_components*sizeof(user_precision_t) ) );
    cudaErrorCheckCall( cudaMemcpy( d_components->majors, components->majors,
                        num_components*sizeof(user_precision_t), cudaMemcpyHostToDevice ) );

    cudaErrorCheckCall( cudaMalloc( (void**)&d_components->minors,
                        num_components*sizeof(user_precision_t) ) );
    cudaErrorCheckCall( cudaMemcpy( d_components->minors, components->minors,
                        num_components*sizeof(user_precision_t), cudaMemcpyHostToDevice ) );
  }

  if (comptype == SHAPELET) {
    cudaErrorCheckCall( cudaMalloc( (void**)&d_components->shape_coeffs,
                        num_shape_coeffs*sizeof(user_precision_t) ) );
    cudaErrorCheckCall( cudaMemcpy( d_components->shape_coeffs, components->shape_coeffs,
                        num_shape_coeffs*sizeof(user_precision_t),
                        cudaMemcpyHostToDevice ) );

    cudaErrorCheckCall( cudaMalloc( (void**)&d_components->n1s,
                        num_shape_coeffs*sizeof(user_precision_t) ) );
    cudaErrorCheckCall( cudaMemcpy( d_components->n1s, components->n1s,
                        num_shape_coeffs*sizeof(user_precision_t),
                        cudaMemcpyHostToDevice ) );

    cudaErrorCheckCall( cudaMalloc( (void**)&d_components->n2s,
                        num_shape_coeffs*sizeof(user_precision_t) ) );
    cudaErrorCheckCall( cudaMemcpy( d_components->n2s, components->n2s,
                        num_shape_coeffs*sizeof(user_precision_t),
                        cudaMemcpyHostToDevice ) );

    cudaErrorCheckCall( cudaMalloc( (void**)&d_components->param_indexes,
                        num_shape_coeffs*sizeof(user_precision_t) ) );
    cudaErrorCheckCall( cudaMemcpy( d_components->param_indexes, components->param_indexes,
                        num_shape_coeffs*sizeof(user_precision_t),
                        cudaMemcpyHostToDevice ) );
  }

  //Only the MWA beams currently yields cross pol values, so only malloc what
  //we need here
  //TODO in the future, this might need to be a loop over all primary beams,
  //if we have different beams for different tiles
  if (beam_settings->beamtype == FEE_BEAM || beam_settings->beamtype == MWA_ANALY || beam_settings->beamtype == FEE_BEAM_INTERP) {
    cudaErrorCheckCall( cudaMalloc( (void**)&d_component_beam_gains->d_Dxs,
                    components->num_primarybeam_values*sizeof(cuUserComplex) ));
    cudaErrorCheckCall( cudaMalloc( (void**)&d_component_beam_gains->d_Dys,
                    components->num_primarybeam_values*sizeof(cuUserComplex) ));
  }
  cudaErrorCheckCall( cudaMalloc( (void**)&d_component_beam_gains->d_gxs,
                    components->num_primarybeam_values*sizeof(cuUserComplex) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_component_beam_gains->d_gys,
                    components->num_primarybeam_values*sizeof(cuUserComplex) ));
  //
  cudaErrorCheckCall( cudaMalloc( (void**)&d_components->ls,
                                               num_components*sizeof(double) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_components->ms,
                                               num_components*sizeof(double) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_components->ns,
                                               num_components*sizeof(double) ) );


  dim3 grid, threads;

  threads.x = 128;
  threads.y = 1;
  threads.z = 1;
  grid.x = (int)ceil( (float)num_components / (float)threads.x );
  grid.y = 1;
  grid.z = 1;

  cudaErrorCheckKernel("kern_calc_lmn",
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

    int num_azza = woden_settings->num_time_steps*num_components;

    double *double_azs = (double*)malloc(num_azza*sizeof(double));
    double *double_zas = (double*)malloc(num_azza*sizeof(double));

    for (int i = 0; i < num_azza; i++) {
      double_azs[i] = (double)components->azs[i];
      double_zas[i] = (double)components->zas[i];
    }

    if (beam_settings->beamtype == FEE_BEAM_INTERP) {
      printf("\tDoing the hyperbeam (interpolated)\n");
    } else {
      printf("\tDoing the hyperbeam\n");
    }


    uint8_t parallactic = 1;
    // int num_freqs = 3;
    run_hyperbeam_cuda(num_components,
           woden_settings->num_time_steps, woden_settings->num_freqs,
           parallactic,
           beam_settings->cuda_fee_beam,
           double_azs, double_zas,
           d_component_beam_gains->d_gxs, d_component_beam_gains->d_Dxs,
           d_component_beam_gains->d_Dys, d_component_beam_gains->d_gys);

    free(double_azs);
    free(double_zas);

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

} //END source_component_common


__global__ void kern_calc_visi_point_or_gauss(components_t d_components,
           d_beam_gains_t d_component_beam_gains,
           user_precision_t *d_us, user_precision_t *d_vs, user_precision_t *d_ws,
           user_precision_t *d_sum_visi_XX_real, user_precision_t *d_sum_visi_XX_imag,
           user_precision_t *d_sum_visi_XY_real, user_precision_t *d_sum_visi_XY_imag,
           user_precision_t *d_sum_visi_YX_real, user_precision_t *d_sum_visi_YX_imag,
           user_precision_t *d_sum_visi_YY_real, user_precision_t *d_sum_visi_YY_imag,
           user_precision_t *d_allsteps_wavelengths, int num_components,
           int num_baselines, int num_freqs, int num_visis,
           int num_times, e_beamtype beamtype, e_component_type comptype) {

  // Start by computing which baseline we're going to do
  const int iBaseline = threadIdx.x + (blockDim.x*blockIdx.x);
  if(iBaseline < num_visis) {

    user_precision_t flux_I;
    user_precision_t flux_Q;
    user_precision_t flux_U;
    user_precision_t flux_V;

    cuUserComplex visi_comp;
    cuUserComplex V_envelop;

    user_precision_t pa, sinpa, cospa, u, v, x, y, invsig_x, invsig_y;

    for (size_t iComponent = 0; iComponent < num_components; iComponent++) {

      extrap_stokes(d_allsteps_wavelengths, d_components.ref_freqs,
                   d_components.ref_stokesI, d_components.ref_stokesQ,
                   d_components.ref_stokesU, d_components.ref_stokesV,
                   d_components.SIs, iComponent, iBaseline,
                   &flux_I, &flux_Q, &flux_U, &flux_V);

      visi_comp = calc_measurement_equation(d_us, d_vs, d_ws,
                             d_components.ls, d_components.ms, d_components.ns,
                             iBaseline, iComponent);

      if (comptype == GAUSSIAN) {

        V_envelop = make_cuUserComplex( 1.0, 0.0 );

        pa = d_components.pas[iComponent];
        sinpa = sin(pa);
        cospa = cos(pa);
        u = d_us[iBaseline];
        v = d_vs[iBaseline];

        x =  cospa*v + sinpa*u; // major axis
        y = -sinpa*v + cospa*u; // minor axis
        invsig_x = d_components.majors[iComponent];
        invsig_y = d_components.minors[iComponent];

        V_envelop = make_cuUserComplex( exp( -0.5 * ( x*x*invsig_x*invsig_x*M_PI_2_2_LN_2 + y*y*invsig_y*invsig_y*M_PI_2_2_LN_2 ) ), 0.0 );

        visi_comp = visi_comp*V_envelop;
      }

      update_sum_visis(iBaseline, iComponent, num_freqs,
             num_baselines, num_components, num_times, beamtype,
             d_component_beam_gains.d_gxs, d_component_beam_gains.d_Dxs,
             d_component_beam_gains.d_Dys, d_component_beam_gains.d_gys,
             visi_comp, flux_I, flux_Q, flux_U, flux_V,
             d_sum_visi_XX_real, d_sum_visi_XX_imag,
             d_sum_visi_XY_real, d_sum_visi_XY_imag,
             d_sum_visi_YX_real, d_sum_visi_YX_imag,
             d_sum_visi_YY_real, d_sum_visi_YY_imag);
    }
  }
}

__global__ void kern_calc_visi_shapelets(components_t d_components,
      d_beam_gains_t d_component_beam_gains,
      user_precision_t *d_us, user_precision_t *d_vs, user_precision_t *d_ws,
      user_precision_t *d_allsteps_wavelengths,
      user_precision_t *d_u_shapes, user_precision_t *d_v_shapes,
      user_precision_t *d_w_shapes,
      user_precision_t *d_sum_visi_XX_real, user_precision_t *d_sum_visi_XX_imag,
      user_precision_t *d_sum_visi_XY_real, user_precision_t *d_sum_visi_XY_imag,
      user_precision_t *d_sum_visi_YX_real, user_precision_t *d_sum_visi_YX_imag,
      user_precision_t *d_sum_visi_YY_real, user_precision_t *d_sum_visi_YY_imag,
      user_precision_t *d_sbf,
      int num_shapes, int num_baselines, int num_freqs, int num_visis,
      const int num_coeffs, int num_times, e_beamtype beamtype) {

  // Start by computing which baseline we're going to do
  const int iBaseline = threadIdx.x + (blockDim.x*blockIdx.x);

  if (iBaseline < num_visis) {

    user_precision_t shape_flux_I;
    user_precision_t shape_flux_Q;
    user_precision_t shape_flux_U;
    user_precision_t shape_flux_V;
    cuUserComplex visi_shape;

    for (int iCoeff = 0; iCoeff < num_coeffs; iCoeff++) {

      int iComponent = d_components.param_indexes[iCoeff];

      extrap_stokes(d_allsteps_wavelengths,  d_components.ref_freqs,
                   d_components.ref_stokesI, d_components.ref_stokesQ,
                   d_components.ref_stokesU, d_components.ref_stokesV,
                   d_components.SIs, iComponent, iBaseline,
                   &shape_flux_I, &shape_flux_Q, &shape_flux_U, &shape_flux_V);

      visi_shape = calc_measurement_equation(d_us, d_vs, d_ws,
                            d_components.ls, d_components.ms, d_components.ns,
                            iBaseline, iComponent);

      user_precision_t pa = d_components.pas[iComponent];
      user_precision_t sinpa = sin(pa);
      user_precision_t cospa = cos(pa);

      user_precision_t u_shape = d_u_shapes[iComponent*num_visis + iBaseline];
      user_precision_t v_shape = d_v_shapes[iComponent*num_visis + iBaseline];

      user_precision_t x = (cospa*v_shape + sinpa*u_shape); // major axis
      user_precision_t y = (-sinpa*v_shape + cospa*u_shape); // minor axis

      //Scales the FWHM to std to match basis functions, and account for the
      //basis functions being stored with beta = 1.0
      //Basis functions have been stored in such a way that x is in the same
      //direction as on sky, but y is opposite, so include negative here
      user_precision_t const_x = (d_components.majors[iComponent]*SQRT_M_PI_2_2_LN_2)/sbf_dx;
      user_precision_t const_y = -(d_components.minors[iComponent]*SQRT_M_PI_2_2_LN_2)/sbf_dx;

      // I^(n1+n2) = Ipow_lookup[(n1+n2) % 4]
      cuUserComplex Ipow_lookup[] = { make_cuUserComplex(  1.0,  0.0 ),
                                       make_cuUserComplex(  0.0,  1.0 ),
                                       make_cuUserComplex( -1.0,  0.0 ),
                                       make_cuUserComplex(  0.0, -1.0 ) };

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
      cuUserComplex V_envelop = make_cuUserComplex( 0.0, 0.0 );
      V_envelop = V_envelop + Ipow_lookup[(n1+n2) % 4] * f_hat * u_value*v_value;

      visi_shape = visi_shape*V_envelop;

      // printf("V_envelop %.5f %.5f\n",V_envelop.x, V_envelop.y );

      update_sum_visis(iBaseline, iComponent, num_freqs,
             num_baselines, num_shapes, num_times, beamtype,
             d_component_beam_gains.d_gxs, d_component_beam_gains.d_Dxs,
             d_component_beam_gains.d_Dys, d_component_beam_gains.d_gys,
             visi_shape,
             shape_flux_I, shape_flux_Q, shape_flux_U, shape_flux_V,
             d_sum_visi_XX_real, d_sum_visi_XX_imag,
             d_sum_visi_XY_real, d_sum_visi_XY_imag,
             d_sum_visi_YX_real, d_sum_visi_YX_imag,
             d_sum_visi_YY_real, d_sum_visi_YY_imag);
    }
  }
}



extern "C" void free_d_components(components_t d_components,
                                  e_component_type comptype){

  cudaErrorCheckCall( cudaFree( d_components.ns) );
  cudaErrorCheckCall( cudaFree( d_components.ms) );
  cudaErrorCheckCall( cudaFree( d_components.ls) );
  cudaErrorCheckCall( cudaFree( d_components.ref_freqs ) );
  cudaErrorCheckCall( cudaFree( d_components.ref_stokesI ) );
  cudaErrorCheckCall( cudaFree( d_components.ref_stokesQ ) );
  cudaErrorCheckCall( cudaFree( d_components.ref_stokesU ) );
  cudaErrorCheckCall( cudaFree( d_components.ref_stokesV ) );
  cudaErrorCheckCall( cudaFree( d_components.SIs ) );
  cudaErrorCheckCall( cudaFree( d_components.decs) );
  cudaErrorCheckCall( cudaFree( d_components.ras) );

  if (comptype == GAUSSIAN || comptype == SHAPELET) {
    cudaErrorCheckCall( cudaFree( d_components.pas ) );
    cudaErrorCheckCall( cudaFree( d_components.majors ) );
    cudaErrorCheckCall( cudaFree( d_components.minors ) );
  }

  if (comptype == SHAPELET) {
    cudaErrorCheckCall( cudaFree( d_components.shape_coeffs ) );
    cudaErrorCheckCall( cudaFree( d_components.n1s ) );
    cudaErrorCheckCall( cudaFree( d_components.n2s ) );
    cudaErrorCheckCall( cudaFree( d_components.param_indexes ) );
  }
}

extern "C" void free_beam_gains(d_beam_gains_t d_beam_gains, e_beamtype beamtype){

  cudaErrorCheckCall( cudaFree( d_beam_gains.d_gxs) );
  cudaErrorCheckCall( cudaFree( d_beam_gains.d_gys) );

  if (beamtype == FEE_BEAM){
    cudaErrorCheckCall( cudaFree( d_beam_gains.d_Dxs ) );
    cudaErrorCheckCall( cudaFree( d_beam_gains.d_Dys ) );
  }

}


/*******************************************************************************
                 Functions below to be used in unit tests
*******************************************************************************/

__global__ void kern_extrap_stokes(int num_extrap_freqs, int num_components,
           user_precision_t *d_extrap_wavelengths, double *d_ref_freqs,
           user_precision_t *d_SIs,
           user_precision_t *d_ref_stokesI, user_precision_t *d_ref_stokesQ,
           user_precision_t *d_ref_stokesU, user_precision_t *d_ref_stokesV,
           user_precision_t *d_flux_I, user_precision_t *d_flux_Q,
           user_precision_t *d_flux_U, user_precision_t *d_flux_V ) {

  // Start by computing which baseline we're going to do
  const int iComponent = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iFreq = threadIdx.y + (blockDim.y*blockIdx.y);
  // if(iBaseline < num_visis && iComponent < num_points) {
  if(iComponent < num_components && iFreq < num_extrap_freqs) {

    user_precision_t flux_I;
    user_precision_t flux_Q;
    user_precision_t flux_U;
    user_precision_t flux_V;

    extrap_stokes(d_extrap_wavelengths, d_ref_freqs,
                 d_ref_stokesI, d_ref_stokesQ,
                 d_ref_stokesU, d_ref_stokesV,
                 d_SIs, iComponent, iFreq,
                 &flux_I, &flux_Q, &flux_U, &flux_V);

    int extrap_ind = num_components*iComponent + iFreq;

    d_flux_I[extrap_ind] = flux_I;
    d_flux_Q[extrap_ind] = flux_Q;
    d_flux_U[extrap_ind] = flux_U;
    d_flux_V[extrap_ind] = flux_V;

  }
}

extern "C" void test_kern_extrap_stokes(int num_extrap_freqs, int num_components,
           user_precision_t *extrap_wavelengths, double *ref_freqs,
           user_precision_t *SIs,
           user_precision_t *ref_stokesI, user_precision_t *ref_stokesQ,
           user_precision_t *ref_stokesU, user_precision_t *ref_stokesV,
           user_precision_t *flux_I, user_precision_t *flux_Q,
           user_precision_t *flux_U, user_precision_t *flux_V){

  user_precision_t *d_extrap_wavelengths = NULL;
  double *d_ref_freqs = NULL;
  user_precision_t *d_SIs = NULL;
  user_precision_t *d_ref_stokesI = NULL;
  user_precision_t *d_ref_stokesQ = NULL;
  user_precision_t *d_ref_stokesU = NULL;
  user_precision_t *d_ref_stokesV = NULL;
  user_precision_t *d_flux_I = NULL;
  user_precision_t *d_flux_Q = NULL;
  user_precision_t *d_flux_U = NULL;
  user_precision_t *d_flux_V = NULL;

  cudaErrorCheckCall( cudaMalloc( (void**)&d_extrap_wavelengths,
                                   num_extrap_freqs*sizeof(user_precision_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_ref_freqs,
                                     num_components*sizeof(double) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_SIs,
                                     num_components*sizeof(user_precision_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_ref_stokesI,
                                     num_components*sizeof(user_precision_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_ref_stokesQ,
                                     num_components*sizeof(user_precision_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_ref_stokesU,
                                     num_components*sizeof(user_precision_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_ref_stokesV,
                                     num_components*sizeof(user_precision_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_flux_I,
                    num_components*num_extrap_freqs*sizeof(user_precision_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_flux_Q,
                    num_components*num_extrap_freqs*sizeof(user_precision_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_flux_U,
                    num_components*num_extrap_freqs*sizeof(user_precision_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_flux_V,
                    num_components*num_extrap_freqs*sizeof(user_precision_t) ));

  cudaErrorCheckCall( cudaMemcpy(d_extrap_wavelengths, extrap_wavelengths,
           num_extrap_freqs*sizeof(user_precision_t), cudaMemcpyHostToDevice ));

  cudaErrorCheckCall( cudaMemcpy(d_ref_stokesI, ref_stokesI,
             num_components*sizeof(user_precision_t), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_ref_stokesQ, ref_stokesQ,
             num_components*sizeof(user_precision_t), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_ref_stokesU, ref_stokesU,
             num_components*sizeof(user_precision_t), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_ref_stokesV, ref_stokesV,
             num_components*sizeof(user_precision_t), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_ref_freqs, ref_freqs,
             num_components*sizeof(double), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_SIs, SIs,
             num_components*sizeof(user_precision_t), cudaMemcpyHostToDevice ));

  dim3 grid, threads;

  threads.x = 16;
  threads.y = 16;
  grid.x = (int)ceilf( (float)num_components / (float)threads.x );
  grid.y = (int)ceilf( (float)num_extrap_freqs / (float)threads.y );

  cudaErrorCheckKernel("kern_extrap_stokes",
                        kern_extrap_stokes, grid, threads,
                        num_extrap_freqs, num_components,
                        d_extrap_wavelengths, d_ref_freqs, d_SIs,
                        d_ref_stokesI, d_ref_stokesQ,
                        d_ref_stokesU, d_ref_stokesV,
                        d_flux_I, d_flux_Q,
                        d_flux_U, d_flux_V);

  cudaErrorCheckCall( cudaMemcpy(flux_I, d_flux_I,
                      num_components*num_extrap_freqs*sizeof(user_precision_t),
                                                      cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(flux_Q, d_flux_Q,
                      num_components*num_extrap_freqs*sizeof(user_precision_t),
                                                      cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(flux_U, d_flux_U,
                      num_components*num_extrap_freqs*sizeof(user_precision_t),
                                                      cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(flux_V, d_flux_V,
                      num_components*num_extrap_freqs*sizeof(user_precision_t),
                                                      cudaMemcpyDeviceToHost ));

  cudaErrorCheckCall( cudaFree( d_extrap_wavelengths ) );
  cudaErrorCheckCall( cudaFree( d_ref_freqs ) );
  cudaErrorCheckCall( cudaFree( d_SIs ) );
  cudaErrorCheckCall( cudaFree( d_ref_stokesI ) );
  cudaErrorCheckCall( cudaFree( d_ref_stokesQ ) );
  cudaErrorCheckCall( cudaFree( d_ref_stokesU ) );
  cudaErrorCheckCall( cudaFree( d_ref_stokesV ) );
  cudaErrorCheckCall( cudaFree( d_flux_I ) );
  cudaErrorCheckCall( cudaFree( d_flux_Q ) );
  cudaErrorCheckCall( cudaFree( d_flux_U ) );
  cudaErrorCheckCall( cudaFree( d_flux_V ) );
}


__global__ void kern_calc_measurement_equation(int num_components, int num_baselines,
          user_precision_t *d_us, user_precision_t *d_vs, user_precision_t *d_ws,
          double *d_ls, double *d_ms, double *d_ns, cuUserComplex *d_visis) {

  // Start by computing which baseline we're going to do
  const int iBaseline = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iComponent = threadIdx.y + (blockDim.y*blockIdx.y);

  if(iComponent < num_components && iBaseline < num_baselines) {

    cuUserComplex visi;
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
  cudaErrorCheckCall( cudaMalloc( (void**)&d_us, num_baselines*sizeof(user_precision_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_vs, num_baselines*sizeof(user_precision_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_ws, num_baselines*sizeof(user_precision_t) ));
  cudaErrorCheckCall( cudaMemcpy(d_us, us, num_baselines*sizeof(user_precision_t),
                                                        cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_vs, vs, num_baselines*sizeof(user_precision_t),
                                                        cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_ws, ws, num_baselines*sizeof(user_precision_t),
                                                        cudaMemcpyHostToDevice ));

  double *d_ls = NULL;
  double *d_ms = NULL;
  double *d_ns = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_ls, num_components*sizeof(double) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_ms, num_components*sizeof(double) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_ns, num_components*sizeof(double) ));
  cudaErrorCheckCall( cudaMemcpy(d_ls, ls, num_components*sizeof(double),
                                                      cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_ms, ms, num_components*sizeof(double),
                                                      cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_ns, ns, num_components*sizeof(double),
                                                      cudaMemcpyHostToDevice ));

  user_precision_complex_t *d_visis = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_visis, num_baselines*num_components*sizeof(user_precision_complex_t) ));

  dim3 grid, threads;

  threads.x = 16;
  threads.y = 16;
  grid.x = (int)ceilf( (float)num_baselines / (float)threads.x );
  grid.y = (int)ceilf( (float)num_components / (float)threads.y );

  cudaErrorCheckKernel("kern_calc_measurement_equation",
                      kern_calc_measurement_equation, grid, threads,
                      num_components, num_baselines,
                      d_us, d_vs, d_ws,
                      d_ls, d_ms, d_ns,
                      (cuUserComplex*)d_visis );

  cudaErrorCheckCall( cudaMemcpy(visis, (user_precision_complex_t*)d_visis, num_components*num_baselines*sizeof(user_precision_complex_t),cudaMemcpyDeviceToHost ));

  cudaErrorCheckCall( cudaFree( d_us ) );
  cudaErrorCheckCall( cudaFree( d_vs ) );
  cudaErrorCheckCall( cudaFree( d_ws ) );
  cudaErrorCheckCall( cudaFree( d_ls ) );
  cudaErrorCheckCall( cudaFree( d_ms ) );
  cudaErrorCheckCall( cudaFree( d_ns ) );
  cudaErrorCheckCall( cudaFree(d_visis ) );

}

__global__ void kern_apply_beam_gains(int num_gains, cuUserComplex *d_g1xs,
          cuUserComplex *d_D1xs,
          cuUserComplex *d_D1ys, cuUserComplex *d_g1ys,
          cuUserComplex *d_g2xs, cuUserComplex *d_D2xs,
          cuUserComplex *d_D2ys, cuUserComplex *d_g2ys,
          user_precision_t *d_flux_Is, user_precision_t *d_flux_Qs,
          user_precision_t *d_flux_Us, user_precision_t *d_flux_Vs,
          cuUserComplex *d_visi_components,
          cuUserComplex *d_visi_XXs, cuUserComplex *d_visi_XYs,
          cuUserComplex *d_visi_YXs, cuUserComplex *d_visi_YYs) {

  const int iGain = threadIdx.x + (blockDim.x*blockIdx.x);
  // const int iComponent = threadIdx.y + (blockDim.y*blockIdx.y);
  // if(iBaseline < num_visis && iComponent < num_points) {
  if (iGain < num_gains) {

    cuUserComplex visi_XX;
    cuUserComplex visi_XY;
    cuUserComplex visi_YX;
    cuUserComplex visi_YY;

    apply_beam_gains(d_g1xs[iGain], d_D1xs[iGain],
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

  cudaErrorCheckCall( cudaMalloc( (void**)&d_g1xs,
                                  num_gains*sizeof(user_precision_complex_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_D1xs,
                                  num_gains*sizeof(user_precision_complex_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_D1ys,
                                  num_gains*sizeof(user_precision_complex_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_g1ys,
                                  num_gains*sizeof(user_precision_complex_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_g2xs,
                                  num_gains*sizeof(user_precision_complex_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_D2xs,
                                  num_gains*sizeof(user_precision_complex_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_D2ys,
                                  num_gains*sizeof(user_precision_complex_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_g2ys,
                                  num_gains*sizeof(user_precision_complex_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_flux_Is,
                                          num_gains*sizeof(user_precision_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_flux_Qs,
                                          num_gains*sizeof(user_precision_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_flux_Us,
                                          num_gains*sizeof(user_precision_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_flux_Vs,
                                          num_gains*sizeof(user_precision_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_visi_components,
                                  num_gains*sizeof(user_precision_complex_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_visi_XXs,
                                  num_gains*sizeof(user_precision_complex_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_visi_XYs,
                                  num_gains*sizeof(user_precision_complex_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_visi_YXs,
                                  num_gains*sizeof(user_precision_complex_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_visi_YYs,
                                  num_gains*sizeof(user_precision_complex_t) ));

  cudaErrorCheckCall( cudaMemcpy(d_g1xs, g1xs,
          num_gains*sizeof(user_precision_complex_t), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_D1xs, D1xs,
          num_gains*sizeof(user_precision_complex_t), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_D1ys, D1ys,
          num_gains*sizeof(user_precision_complex_t), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_g1ys, g1ys,
          num_gains*sizeof(user_precision_complex_t), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_g2xs, g2xs,
          num_gains*sizeof(user_precision_complex_t), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_D2xs, D2xs,
          num_gains*sizeof(user_precision_complex_t), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_D2ys, D2ys,
          num_gains*sizeof(user_precision_complex_t), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_g2ys, g2ys,
          num_gains*sizeof(user_precision_complex_t), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_visi_components, visi_components,
          num_gains*sizeof(user_precision_complex_t), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_visi_XXs, visi_XXs,
          num_gains*sizeof(user_precision_complex_t), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_visi_XYs, visi_XYs,
          num_gains*sizeof(user_precision_complex_t), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_visi_YXs, visi_YXs,
          num_gains*sizeof(user_precision_complex_t), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_visi_YYs, visi_YYs,
          num_gains*sizeof(user_precision_complex_t), cudaMemcpyHostToDevice ));

  cudaErrorCheckCall( cudaMemcpy(d_flux_Is, flux_Is,
                             num_gains*sizeof(user_precision_t), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_flux_Qs, flux_Qs,
                             num_gains*sizeof(user_precision_t), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_flux_Us, flux_Us,
                             num_gains*sizeof(user_precision_t), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_flux_Vs, flux_Vs,
                             num_gains*sizeof(user_precision_t), cudaMemcpyHostToDevice ));

  dim3 grid, threads;

  threads.x = 128;
  grid.x = (int)ceil( (user_precision_t)num_gains / (user_precision_t)threads.x );

  cudaErrorCheckKernel("kern_apply_beam_gains",
                      kern_apply_beam_gains, grid, threads,
                      num_gains,
                      (cuUserComplex *)d_g1xs, (cuUserComplex *)d_D1xs,
                      (cuUserComplex *)d_D1ys, (cuUserComplex *)d_g1ys,
                      (cuUserComplex *)d_g2xs, (cuUserComplex *)d_D2xs,
                      (cuUserComplex *)d_D2ys, (cuUserComplex *)d_g2ys,
                      d_flux_Is, d_flux_Qs,
                      d_flux_Us, d_flux_Vs,
                      (cuUserComplex *)d_visi_components,
                      (cuUserComplex *)d_visi_XXs, (cuUserComplex *)d_visi_XYs,
                      (cuUserComplex *)d_visi_YXs, (cuUserComplex *)d_visi_YYs );

  cudaErrorCheckCall( cudaMemcpy(visi_XXs, d_visi_XXs,
           num_gains*sizeof(user_precision_complex_t),cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(visi_XYs, d_visi_XYs,
           num_gains*sizeof(user_precision_complex_t),cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(visi_YXs, d_visi_YXs,
           num_gains*sizeof(user_precision_complex_t),cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(visi_YYs, d_visi_YYs,
           num_gains*sizeof(user_precision_complex_t),cudaMemcpyDeviceToHost ));

  cudaErrorCheckCall( cudaFree( d_g1xs ) );
  cudaErrorCheckCall( cudaFree( d_D1xs ) );
  cudaErrorCheckCall( cudaFree( d_D1ys ) );
  cudaErrorCheckCall( cudaFree( d_g1ys ) );
  cudaErrorCheckCall( cudaFree( d_g2xs ) );
  cudaErrorCheckCall( cudaFree( d_D2xs ) );
  cudaErrorCheckCall( cudaFree( d_D2ys ) );
  cudaErrorCheckCall( cudaFree( d_g2ys ) );
  cudaErrorCheckCall( cudaFree( d_flux_Is ) );
  cudaErrorCheckCall( cudaFree( d_flux_Qs ) );
  cudaErrorCheckCall( cudaFree( d_flux_Us ) );
  cudaErrorCheckCall( cudaFree( d_flux_Vs ) );
  cudaErrorCheckCall( cudaFree( d_visi_components ) );
  cudaErrorCheckCall( cudaFree( d_visi_XXs ) );
  cudaErrorCheckCall( cudaFree( d_visi_XYs ) );
  cudaErrorCheckCall( cudaFree( d_visi_YXs ) );
  cudaErrorCheckCall( cudaFree( d_visi_YYs ) );

}

__global__ void kern_get_beam_gains(int num_components, int num_baselines,
           int num_freqs, int num_visis, int num_times, int beamtype,
           cuUserComplex *d_primay_beam_J00, cuUserComplex *d_primay_beam_J01,
           cuUserComplex *d_primay_beam_J10, cuUserComplex *d_primay_beam_J11,
           cuUserComplex *d_recov_g1x, cuUserComplex *d_recov_D1x,
           cuUserComplex *d_recov_D1y, cuUserComplex *d_recov_g1y,
           cuUserComplex *d_recov_g2x, cuUserComplex *d_recov_D2x,
           cuUserComplex *d_recov_D2y, cuUserComplex *d_recov_g2y) {

  // Start by computing which baseline we're going to do
  const int iBaseline = threadIdx.x + (blockDim.x*blockIdx.x);
  // const int iComponent = threadIdx.y + (blockDim.y*blockIdx.y);
  // if(iBaseline < num_visis && iComponent < num_points) {
  if(iBaseline < num_visis) {

    for (int iComponent = 0; iComponent < num_components; iComponent++) {

      cuUserComplex g1x;
      cuUserComplex D1x;
      cuUserComplex D1y;
      cuUserComplex g1y;
      cuUserComplex g2x;
      cuUserComplex D2x;
      cuUserComplex D2y;
      cuUserComplex g2y;

      get_beam_gains(iBaseline, iComponent, num_freqs,
                 num_baselines, num_components, num_times, beamtype,
                 d_primay_beam_J00, d_primay_beam_J01,
                 d_primay_beam_J10, d_primay_beam_J11,
                 &g1x, &D1x, &D1y, &g1y, &g2x, &D2x, &D2y, &g2y);

      int out_ind = num_visis*iComponent + iBaseline;

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

extern "C" void test_kern_get_beam_gains(int num_freqs, int num_visis,
          int num_baselines, int num_components, int num_times, int beamtype,
          user_precision_complex_t *primay_beam_J00, user_precision_complex_t *primay_beam_J01,
          user_precision_complex_t *primay_beam_J10, user_precision_complex_t *primay_beam_J11,
          user_precision_complex_t *recover_g1x, user_precision_complex_t *recover_D1x,
          user_precision_complex_t *recover_D1y, user_precision_complex_t *recover_g1y,
          user_precision_complex_t *recover_g2x, user_precision_complex_t *recover_D2x,
          user_precision_complex_t *recover_D2y, user_precision_complex_t *recover_g2y){

  user_precision_complex_t *d_recover_g1x = NULL;
  user_precision_complex_t *d_recover_D1x = NULL;
  user_precision_complex_t *d_recover_D1y = NULL;
  user_precision_complex_t *d_recover_g1y = NULL;
  user_precision_complex_t *d_recover_g2x = NULL;
  user_precision_complex_t *d_recover_D2x = NULL;
  user_precision_complex_t *d_recover_D2y = NULL;
  user_precision_complex_t *d_recover_g2y = NULL;

  user_precision_complex_t *d_primay_beam_J00 = NULL;
  user_precision_complex_t *d_primay_beam_J01 = NULL;
  user_precision_complex_t *d_primay_beam_J10 = NULL;
  user_precision_complex_t *d_primay_beam_J11 = NULL;

  cudaErrorCheckCall( cudaMalloc( (void**)&d_recover_g1x, num_components*num_visis*sizeof(user_precision_complex_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_recover_D1x, num_components*num_visis*sizeof(user_precision_complex_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_recover_D1y, num_components*num_visis*sizeof(user_precision_complex_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_recover_g1y, num_components*num_visis*sizeof(user_precision_complex_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_recover_g2x, num_components*num_visis*sizeof(user_precision_complex_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_recover_D2x, num_components*num_visis*sizeof(user_precision_complex_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_recover_D2y, num_components*num_visis*sizeof(user_precision_complex_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_recover_g2y, num_components*num_visis*sizeof(user_precision_complex_t) ));

  cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J00, num_freqs*num_times*num_components*sizeof(user_precision_complex_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J01, num_freqs*num_times*num_components*sizeof(user_precision_complex_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J10, num_freqs*num_times*num_components*sizeof(user_precision_complex_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J11, num_freqs*num_times*num_components*sizeof(user_precision_complex_t) ));

  cudaErrorCheckCall( cudaMemcpy(d_primay_beam_J00, primay_beam_J00, num_freqs*num_times*num_components*sizeof(user_precision_complex_t), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_primay_beam_J01, primay_beam_J01, num_freqs*num_times*num_components*sizeof(user_precision_complex_t), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_primay_beam_J10, primay_beam_J10, num_freqs*num_times*num_components*sizeof(user_precision_complex_t), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_primay_beam_J11, primay_beam_J11, num_freqs*num_times*num_components*sizeof(user_precision_complex_t), cudaMemcpyHostToDevice ));

  dim3 grid, threads;

  threads.x = 128;
  grid.x = (int)ceil( (user_precision_t)num_visis / (user_precision_t)threads.x );

  cudaErrorCheckKernel("kern_get_beam_gains",
                      kern_get_beam_gains, grid, threads,
                      num_components, num_baselines,
                      num_freqs, num_visis, num_times, beamtype,
                      (cuUserComplex *)d_primay_beam_J00,
                      (cuUserComplex *)d_primay_beam_J01,
                      (cuUserComplex *)d_primay_beam_J10,
                      (cuUserComplex *)d_primay_beam_J11,
                      (cuUserComplex *)d_recover_g1x, (cuUserComplex *)d_recover_D1x,
                      (cuUserComplex *)d_recover_D1y, (cuUserComplex *)d_recover_g1y,
                      (cuUserComplex *)d_recover_g2x, (cuUserComplex *)d_recover_D2x,
                      (cuUserComplex *)d_recover_D2y, (cuUserComplex *)d_recover_g2y );

  cudaErrorCheckCall( cudaMemcpy(recover_g1x, d_recover_g1x, num_components*num_visis*sizeof(user_precision_complex_t), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(recover_D1x, d_recover_D1x, num_components*num_visis*sizeof(user_precision_complex_t), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(recover_D1y, d_recover_D1y, num_components*num_visis*sizeof(user_precision_complex_t), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(recover_g1y, d_recover_g1y, num_components*num_visis*sizeof(user_precision_complex_t), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(recover_g2x, d_recover_g2x, num_components*num_visis*sizeof(user_precision_complex_t), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(recover_D2x, d_recover_D2x, num_components*num_visis*sizeof(user_precision_complex_t), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(recover_D2y, d_recover_D2y, num_components*num_visis*sizeof(user_precision_complex_t), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(recover_g2y, d_recover_g2y, num_components*num_visis*sizeof(user_precision_complex_t), cudaMemcpyDeviceToHost ));

  cudaErrorCheckCall( cudaFree( d_recover_g1x ) );
  cudaErrorCheckCall( cudaFree( d_recover_D1x ) );
  cudaErrorCheckCall( cudaFree( d_recover_D1y ) );
  cudaErrorCheckCall( cudaFree( d_recover_g1y ) );
  cudaErrorCheckCall( cudaFree( d_recover_g2x ) );
  cudaErrorCheckCall( cudaFree( d_recover_D2x ) );
  cudaErrorCheckCall( cudaFree( d_recover_D2y ) );
  cudaErrorCheckCall( cudaFree( d_recover_g2y ) );

  cudaErrorCheckCall( cudaFree( d_primay_beam_J00 ) );
  cudaErrorCheckCall( cudaFree( d_primay_beam_J01 ) );
  cudaErrorCheckCall( cudaFree( d_primay_beam_J10 ) );
  cudaErrorCheckCall( cudaFree( d_primay_beam_J11 ) );

}

__global__ void kern_update_sum_visis(int num_freqs,
     int num_baselines, int num_components, int num_times, int beamtype,
     cuUserComplex *d_primay_beam_J00, cuUserComplex *d_primay_beam_J01,
     cuUserComplex *d_primay_beam_J10, cuUserComplex *d_primay_beam_J11,
     cuUserComplex *d_visi_components,
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

      update_sum_visis(iBaseline, iComponent, num_freqs,
             num_baselines, num_components, num_times, beamtype,
             d_primay_beam_J00, d_primay_beam_J01,
             d_primay_beam_J10, d_primay_beam_J11,
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

extern "C" void test_kern_update_sum_visis(int num_freqs, int num_visis,
          int num_baselines, int num_components, int num_times, int beamtype,
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

  user_precision_complex_t *d_primay_beam_J00 = NULL;
  user_precision_complex_t *d_primay_beam_J01 = NULL;
  user_precision_complex_t *d_primay_beam_J10 = NULL;
  user_precision_complex_t *d_primay_beam_J11 = NULL;
  user_precision_complex_t *d_visi_components = NULL;

  cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J00,
                    num_components*num_times*num_freqs*sizeof(user_precision_complex_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J01,
                    num_components*num_times*num_freqs*sizeof(user_precision_complex_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J10,
                    num_components*num_times*num_freqs*sizeof(user_precision_complex_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J11,
                    num_components*num_times*num_freqs*sizeof(user_precision_complex_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_visi_components,
                    num_visis*sizeof(user_precision_complex_t) ));

  cudaErrorCheckCall( cudaMemcpy(d_primay_beam_J00, primay_beam_J00,
            num_components*num_times*num_freqs*sizeof(user_precision_complex_t),
            cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_primay_beam_J01, primay_beam_J01,
            num_components*num_times*num_freqs*sizeof(user_precision_complex_t),
            cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_primay_beam_J10, primay_beam_J10,
            num_components*num_times*num_freqs*sizeof(user_precision_complex_t),
            cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_primay_beam_J11, primay_beam_J11,
            num_components*num_times*num_freqs*sizeof(user_precision_complex_t),
            cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_visi_components, visi_components,
                                     num_visis*sizeof(user_precision_complex_t),
                                     cudaMemcpyHostToDevice ));

  user_precision_t *d_flux_I = NULL;
  user_precision_t *d_flux_Q = NULL;
  user_precision_t *d_flux_U = NULL;
  user_precision_t *d_flux_V = NULL;

  cudaErrorCheckCall( cudaMalloc( (void**)&d_flux_I, num_components*num_times*num_freqs*sizeof(user_precision_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_flux_Q, num_components*num_times*num_freqs*sizeof(user_precision_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_flux_U, num_components*num_times*num_freqs*sizeof(user_precision_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_flux_V, num_components*num_times*num_freqs*sizeof(user_precision_t) ));

  cudaErrorCheckCall( cudaMemcpy(d_flux_I, flux_I,
                    num_components*num_times*num_freqs*sizeof(user_precision_t),    cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_flux_Q, flux_Q,
                    num_components*num_times*num_freqs*sizeof(user_precision_t),    cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_flux_U, flux_U,
                    num_components*num_times*num_freqs*sizeof(user_precision_t),    cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_flux_V, flux_V,
                    num_components*num_times*num_freqs*sizeof(user_precision_t),    cudaMemcpyHostToDevice ));

  user_precision_t *d_sum_visi_XX_real = NULL;
  user_precision_t *d_sum_visi_XY_real = NULL;
  user_precision_t *d_sum_visi_YX_real = NULL;
  user_precision_t *d_sum_visi_YY_real = NULL;
  user_precision_t *d_sum_visi_XX_imag = NULL;
  user_precision_t *d_sum_visi_XY_imag = NULL;
  user_precision_t *d_sum_visi_YX_imag = NULL;
  user_precision_t *d_sum_visi_YY_imag = NULL;

  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_XX_real,
                                          num_visis*sizeof(user_precision_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_XY_real,
                                          num_visis*sizeof(user_precision_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_YX_real,
                                          num_visis*sizeof(user_precision_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_YY_real,
                                          num_visis*sizeof(user_precision_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_XX_imag,
                                          num_visis*sizeof(user_precision_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_XY_imag,
                                          num_visis*sizeof(user_precision_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_YX_imag,
                                          num_visis*sizeof(user_precision_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_YY_imag,
                                          num_visis*sizeof(user_precision_t) ));

  dim3 grid, threads;

  threads.x = 128;
  grid.x = (int)ceil( (user_precision_t)num_visis / (user_precision_t)threads.x );

  cudaErrorCheckKernel("kern_update_sum_visis",
                      kern_update_sum_visis, grid, threads,
                      num_freqs, num_baselines, num_components, num_times, beamtype,
                      (cuUserComplex *)d_primay_beam_J00, (cuUserComplex *)d_primay_beam_J01,
                      (cuUserComplex *)d_primay_beam_J10, (cuUserComplex *)d_primay_beam_J11,
                      (cuUserComplex *)d_visi_components,
                      d_flux_I, d_flux_Q, d_flux_U, d_flux_V,
                      d_sum_visi_XX_real, d_sum_visi_XX_imag,
                      d_sum_visi_XY_real, d_sum_visi_XY_imag,
                      d_sum_visi_YX_real, d_sum_visi_YX_imag,
                      d_sum_visi_YY_real, d_sum_visi_YY_imag );

  cudaErrorCheckCall( cudaMemcpy(sum_visi_XX_real, d_sum_visi_XX_real,
                  num_visis*sizeof(user_precision_t), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(sum_visi_XY_real, d_sum_visi_XY_real,
                  num_visis*sizeof(user_precision_t), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(sum_visi_YX_real, d_sum_visi_YX_real,
                  num_visis*sizeof(user_precision_t), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(sum_visi_YY_real, d_sum_visi_YY_real,
                  num_visis*sizeof(user_precision_t), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(sum_visi_XX_imag, d_sum_visi_XX_imag,
                  num_visis*sizeof(user_precision_t), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(sum_visi_XY_imag, d_sum_visi_XY_imag,
                  num_visis*sizeof(user_precision_t), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(sum_visi_YX_imag, d_sum_visi_YX_imag,
                  num_visis*sizeof(user_precision_t), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(sum_visi_YY_imag, d_sum_visi_YY_imag,
                  num_visis*sizeof(user_precision_t), cudaMemcpyDeviceToHost ));

  cudaErrorCheckCall( cudaFree( d_primay_beam_J00 ) );
  cudaErrorCheckCall( cudaFree( d_primay_beam_J01 ) );
  cudaErrorCheckCall( cudaFree( d_primay_beam_J10 ) );
  cudaErrorCheckCall( cudaFree( d_primay_beam_J11 ) );
  cudaErrorCheckCall( cudaFree( d_visi_components ) );
  cudaErrorCheckCall( cudaFree( d_flux_I ) );
  cudaErrorCheckCall( cudaFree( d_flux_Q ) );
  cudaErrorCheckCall( cudaFree( d_flux_U ) );
  cudaErrorCheckCall( cudaFree( d_flux_V ) );
  cudaErrorCheckCall( cudaFree( d_sum_visi_XX_real ) );
  cudaErrorCheckCall( cudaFree( d_sum_visi_XY_real ) );
  cudaErrorCheckCall( cudaFree( d_sum_visi_YX_real ) );
  cudaErrorCheckCall( cudaFree( d_sum_visi_YY_real ) );
  cudaErrorCheckCall( cudaFree( d_sum_visi_XX_imag ) );
  cudaErrorCheckCall( cudaFree( d_sum_visi_XY_imag ) );
  cudaErrorCheckCall( cudaFree( d_sum_visi_YX_imag ) );
  cudaErrorCheckCall( cudaFree( d_sum_visi_YY_imag ) );

}


extern "C" void test_source_component_common(int num_components,
           int num_shape_coeffs, components_t components,
           components_t d_components,
           double *freqs, woden_settings_t *woden_settings,
           beam_settings_t *beam_settings,
           user_precision_complex_t *gxs, user_precision_complex_t *Dxs,
           user_precision_complex_t *Dys, user_precision_complex_t *gys,
           double *ls, double *ms, double *ns){

  double *d_freqs = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_freqs,
                                     woden_settings->num_freqs*sizeof(double) ));
  cudaErrorCheckCall( cudaMemcpy( d_freqs, freqs,
             woden_settings->num_freqs*sizeof(double), cudaMemcpyHostToDevice) );


  d_beam_gains_t d_beam_gains;

  source_component_common(num_components, num_shape_coeffs,
             &components, d_freqs, woden_settings, beam_settings, POINT,
             &d_components, &d_beam_gains);

  int num_beam_values = num_components*woden_settings->num_freqs*woden_settings->num_time_steps;

  cudaErrorCheckCall( cudaMemcpy(gxs, (user_precision_complex_t*)d_beam_gains.d_gxs,
              num_beam_values*sizeof(cuUserComplex), cudaMemcpyDeviceToHost ));

  cudaErrorCheckCall( cudaMemcpy(gys, (user_precision_complex_t*)d_beam_gains.d_gys,
              num_beam_values*sizeof(cuUserComplex), cudaMemcpyDeviceToHost ));

  if (beam_settings->beamtype == FEE_BEAM || beam_settings->beamtype == FEE_BEAM_INTERP || beam_settings->beamtype == MWA_ANALY) {
    cudaErrorCheckCall( cudaMemcpy(Dxs, (user_precision_complex_t*)d_beam_gains.d_Dxs,
                num_beam_values*sizeof(cuUserComplex), cudaMemcpyDeviceToHost ));
    cudaErrorCheckCall( cudaMemcpy(Dys, (user_precision_complex_t*)d_beam_gains.d_Dys,
                num_beam_values*sizeof(cuUserComplex), cudaMemcpyDeviceToHost ));
  }

  cudaErrorCheckCall( cudaMemcpy(ls, d_components.ls,
                        num_components*sizeof(double), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(ms, d_components.ms,
                        num_components*sizeof(double), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(ns, d_components.ns,
                        num_components*sizeof(double), cudaMemcpyDeviceToHost ));

  free_d_components(d_components, POINT);
  free_beam_gains(d_beam_gains, beam_settings->beamtype);
}

extern "C" void test_kern_calc_visi_all(int num_components,
          int num_baselines, int num_shape_coeffs,
          int num_freqs, int num_visis, int num_times,
          e_beamtype beamtype, e_component_type comptype,
          components_t components,
          user_precision_t *us, user_precision_t *vs, user_precision_t *ws,
          user_precision_t *u_shapes, user_precision_t *v_shapes, user_precision_t *w_shapes,
          user_precision_t *sum_visi_XX_real, user_precision_t *sum_visi_XX_imag,
          user_precision_t *sum_visi_XY_real, user_precision_t *sum_visi_XY_imag,
          user_precision_t *sum_visi_YX_real, user_precision_t *sum_visi_YX_imag,
          user_precision_t *sum_visi_YY_real, user_precision_t *sum_visi_YY_imag,
          user_precision_t *allsteps_wavelengths, user_precision_t *sbf,
          user_precision_complex_t *gxs, user_precision_complex_t *Dxs,
          user_precision_complex_t *Dys, user_precision_complex_t *gys){

  user_precision_t *d_us = NULL;
  user_precision_t *d_vs = NULL;
  user_precision_t *d_ws = NULL;
  user_precision_t *d_allsteps_wavelengths = NULL;

  cudaErrorCheckCall( cudaMalloc( (void**)&d_us, num_visis*sizeof(user_precision_t) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_vs, num_visis*sizeof(user_precision_t) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_ws, num_visis*sizeof(user_precision_t) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_allsteps_wavelengths, num_visis*sizeof(user_precision_t) ) );

  cudaErrorCheckCall( cudaMemcpy(d_us, us,
                             num_visis*sizeof(user_precision_t), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_vs, vs,
                             num_visis*sizeof(user_precision_t), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_ws, ws,
                             num_visis*sizeof(user_precision_t), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_allsteps_wavelengths, allsteps_wavelengths,
                             num_visis*sizeof(user_precision_t), cudaMemcpyHostToDevice ));

  //Something to store the device component types in
  components_t d_components;

  //just malloc these here are the overall free-ing function later needs to free them
  cudaErrorCheckCall( cudaMalloc( (void**)&d_components.ras,
                                                num_components*sizeof(double) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_components.decs,
                                                num_components*sizeof(double) ));

  cudaErrorCheckCall( cudaMalloc( (void**)&d_components.ref_stokesI,
                                               num_components*sizeof(user_precision_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_components.ref_stokesQ,
                                               num_components*sizeof(user_precision_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_components.ref_stokesU,
                                               num_components*sizeof(user_precision_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_components.ref_stokesV,
                                               num_components*sizeof(user_precision_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_components.SIs,
                                               num_components*sizeof(user_precision_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_components.ref_freqs,
                                                num_components*sizeof(double) ));

  cudaErrorCheckCall( cudaMemcpy(d_components.ref_stokesI, components.ref_stokesI,
                             num_components*sizeof(user_precision_t), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_components.ref_stokesQ, components.ref_stokesQ,
                             num_components*sizeof(user_precision_t), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_components.ref_stokesU, components.ref_stokesU,
                             num_components*sizeof(user_precision_t), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_components.ref_stokesV, components.ref_stokesV,
                             num_components*sizeof(user_precision_t), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_components.SIs, components.SIs,
                             num_components*sizeof(user_precision_t), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_components.ref_freqs, components.ref_freqs,
                        num_components*sizeof(double), cudaMemcpyHostToDevice ));

  //This would be done by source_component_common
  cudaErrorCheckCall( cudaMalloc( (void**)&d_components.ls, num_components*sizeof(double) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_components.ms, num_components*sizeof(double) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_components.ns, num_components*sizeof(double) ));

  cudaErrorCheckCall( cudaMemcpy(d_components.ls, components.ls, num_components*sizeof(double),
                                           cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_components.ms, components.ms, num_components*sizeof(double),
                                           cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_components.ns, components.ns, num_components*sizeof(double),
                                           cudaMemcpyHostToDevice ));


  //Something to store the primary beam gains (all 4 pols) in
  d_beam_gains_t d_beam_gains;
  int num_beam_values = num_components*num_freqs*num_times;

  cudaErrorCheckCall( cudaMalloc( (void**)&d_beam_gains.d_gxs,
                                      num_beam_values*sizeof(cuUserComplex) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_beam_gains.d_Dxs,
                                      num_beam_values*sizeof(cuUserComplex) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_beam_gains.d_Dys,
                                      num_beam_values*sizeof(cuUserComplex) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_beam_gains.d_gys,
                                      num_beam_values*sizeof(cuUserComplex) ));

  cudaErrorCheckCall( cudaMemcpy(d_beam_gains.d_gxs, (cuUserComplex *)gxs,
              num_beam_values*sizeof(cuUserComplex), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_beam_gains.d_Dxs, (cuUserComplex *)Dxs,
              num_beam_values*sizeof(cuUserComplex), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_beam_gains.d_Dys, (cuUserComplex *)Dys,
              num_beam_values*sizeof(cuUserComplex), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_beam_gains.d_gys, (cuUserComplex *)gys,
              num_beam_values*sizeof(cuUserComplex), cudaMemcpyHostToDevice ));

  user_precision_t *d_sum_visi_XX_real = NULL;
  user_precision_t *d_sum_visi_XY_real = NULL;
  user_precision_t *d_sum_visi_YX_real = NULL;
  user_precision_t *d_sum_visi_YY_real = NULL;
  user_precision_t *d_sum_visi_XX_imag = NULL;
  user_precision_t *d_sum_visi_XY_imag = NULL;
  user_precision_t *d_sum_visi_YX_imag = NULL;
  user_precision_t *d_sum_visi_YY_imag = NULL;

  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_XX_real,
                                          num_visis*sizeof(user_precision_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_XY_real,
                                          num_visis*sizeof(user_precision_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_YX_real,
                                          num_visis*sizeof(user_precision_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_YY_real,
                                          num_visis*sizeof(user_precision_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_XX_imag,
                                          num_visis*sizeof(user_precision_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_XY_imag,
                                          num_visis*sizeof(user_precision_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_YX_imag,
                                          num_visis*sizeof(user_precision_t) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_YY_imag,
                                          num_visis*sizeof(user_precision_t) ));

  //Make sure the visis start at zero by copying across host versions, which
  //should be set to zero already
  cudaErrorCheckCall( cudaMemcpy( d_sum_visi_XX_real, sum_visi_XX_real,
    num_visis*sizeof(user_precision_t), cudaMemcpyHostToDevice ) );
  cudaErrorCheckCall( cudaMemcpy( d_sum_visi_XY_real, sum_visi_XY_real,
    num_visis*sizeof(user_precision_t), cudaMemcpyHostToDevice ) );
  cudaErrorCheckCall( cudaMemcpy( d_sum_visi_YX_real, sum_visi_YX_real,
    num_visis*sizeof(user_precision_t), cudaMemcpyHostToDevice ) );
  cudaErrorCheckCall( cudaMemcpy( d_sum_visi_YY_real, sum_visi_YY_real,
    num_visis*sizeof(user_precision_t), cudaMemcpyHostToDevice ) );
  cudaErrorCheckCall( cudaMemcpy( d_sum_visi_XX_imag, sum_visi_XX_imag,
    num_visis*sizeof(user_precision_t), cudaMemcpyHostToDevice ) );
  cudaErrorCheckCall( cudaMemcpy( d_sum_visi_XY_imag, sum_visi_XY_imag,
    num_visis*sizeof(user_precision_t), cudaMemcpyHostToDevice ) );
  cudaErrorCheckCall( cudaMemcpy( d_sum_visi_YX_imag, sum_visi_YX_imag,
    num_visis*sizeof(user_precision_t), cudaMemcpyHostToDevice ) );
  cudaErrorCheckCall( cudaMemcpy( d_sum_visi_YY_imag, sum_visi_YY_imag,
    num_visis*sizeof(user_precision_t), cudaMemcpyHostToDevice ) );

  dim3 grid, threads;

  threads.x = 128;
  grid.x = (int)ceil( (float)num_visis / (float)threads.x );

  if (comptype == GAUSSIAN || comptype == SHAPELET ) {

    //This would be done by source_component_common
    cudaErrorCheckCall( cudaMalloc( (void**)&d_components.pas,
                                                num_components*sizeof(user_precision_t) ));
    cudaErrorCheckCall( cudaMalloc( (void**)&d_components.majors,
                                                num_components*sizeof(user_precision_t) ));
    cudaErrorCheckCall( cudaMalloc( (void**)&d_components.minors,
                                                num_components*sizeof(user_precision_t) ));

    cudaErrorCheckCall( cudaMemcpy(d_components.pas, components.pas,
                       num_components*sizeof(user_precision_t), cudaMemcpyHostToDevice ));
    cudaErrorCheckCall( cudaMemcpy(d_components.majors, components.majors,
                       num_components*sizeof(user_precision_t), cudaMemcpyHostToDevice ));
    cudaErrorCheckCall( cudaMemcpy(d_components.minors, components.minors,
                       num_components*sizeof(user_precision_t), cudaMemcpyHostToDevice ));
  }

  //Shapelets need many many extra things

  user_precision_t *d_sbf=NULL;
  user_precision_t *d_u_shapes = NULL;
  user_precision_t *d_v_shapes = NULL;
  user_precision_t *d_w_shapes = NULL;

  if (comptype == SHAPELET) {

    cudaErrorCheckCall( cudaMalloc( (void**)&d_u_shapes,
                                     num_components*num_visis*sizeof(user_precision_t) ) );
    cudaErrorCheckCall( cudaMalloc( (void**)&d_v_shapes,
                                     num_components*num_visis*sizeof(user_precision_t) ) );
    cudaErrorCheckCall( cudaMalloc( (void**)&d_w_shapes,
                                     num_components*num_visis*sizeof(user_precision_t) ) );

    cudaErrorCheckCall( cudaMemcpy(d_u_shapes, u_shapes,
               num_components*num_visis*sizeof(user_precision_t), cudaMemcpyHostToDevice ));
    cudaErrorCheckCall( cudaMemcpy(d_v_shapes, v_shapes,
               num_components*num_visis*sizeof(user_precision_t), cudaMemcpyHostToDevice ));
    cudaErrorCheckCall( cudaMemcpy(d_w_shapes, w_shapes,
               num_components*num_visis*sizeof(user_precision_t), cudaMemcpyHostToDevice ));

    cudaErrorCheckCall( cudaMalloc( (void**)&d_components.shape_coeffs,
                                                num_shape_coeffs*sizeof(user_precision_t) ));
    cudaErrorCheckCall( cudaMalloc( (void**)&d_components.n1s,
                                                num_shape_coeffs*sizeof(user_precision_t) ));
    cudaErrorCheckCall( cudaMalloc( (void**)&d_components.n2s,
                                                num_shape_coeffs*sizeof(user_precision_t) ));
    cudaErrorCheckCall( cudaMalloc( (void**)&d_components.param_indexes,
                                                num_shape_coeffs*sizeof(user_precision_t) ));

    cudaErrorCheckCall( cudaMemcpy(d_components.shape_coeffs,
                          components.shape_coeffs, num_shape_coeffs*sizeof(user_precision_t),
                          cudaMemcpyHostToDevice ));
    cudaErrorCheckCall( cudaMemcpy(d_components.n1s,
                          components.n1s, num_shape_coeffs*sizeof(user_precision_t),
                          cudaMemcpyHostToDevice ));
    cudaErrorCheckCall( cudaMemcpy(d_components.n2s,
                          components.n2s, num_shape_coeffs*sizeof(user_precision_t),
                          cudaMemcpyHostToDevice ));
    cudaErrorCheckCall( cudaMemcpy(d_components.param_indexes,
                          components.param_indexes, num_shape_coeffs*sizeof(user_precision_t),
                          cudaMemcpyHostToDevice ));
    cudaErrorCheckCall( cudaMalloc( (void**)&(d_sbf), sbf_N*sbf_L*sizeof(user_precision_t) ));
    cudaErrorCheckCall( cudaMemcpy( d_sbf, sbf, sbf_N*sbf_L*sizeof(user_precision_t),
                        cudaMemcpyHostToDevice ));
  }

  if (comptype == POINT || comptype == GAUSSIAN ) {

    cudaErrorCheckKernel("kern_calc_visi_point_or_gauss",
                  kern_calc_visi_point_or_gauss, grid, threads,
                  d_components, d_beam_gains,
                  d_us, d_vs, d_ws,
                  d_sum_visi_XX_real, d_sum_visi_XX_imag,
                  d_sum_visi_XY_real, d_sum_visi_XY_imag,
                  d_sum_visi_YX_real, d_sum_visi_YX_imag,
                  d_sum_visi_YY_real, d_sum_visi_YY_imag,
                  d_allsteps_wavelengths, num_components,
                  num_baselines, num_freqs, num_visis,
                  num_times, beamtype, comptype);
  }
  else if (comptype == SHAPELET) {
    cudaErrorCheckKernel("kern_calc_visi_shapelets",
                  kern_calc_visi_shapelets, grid, threads,
                  d_components, d_beam_gains,
                  d_us, d_vs, d_ws,
                  d_allsteps_wavelengths,
                  d_u_shapes, d_v_shapes, d_w_shapes,
                  d_sum_visi_XX_real, d_sum_visi_XX_imag,
                  d_sum_visi_XY_real, d_sum_visi_XY_imag,
                  d_sum_visi_YX_real, d_sum_visi_YX_imag,
                  d_sum_visi_YY_real, d_sum_visi_YY_imag,
                  d_sbf,  num_components,
                  num_baselines, num_freqs, num_visis,
                  num_shape_coeffs, num_times, beamtype);
  }

  cudaErrorCheckCall( cudaMemcpy(sum_visi_XX_real, d_sum_visi_XX_real,
                             num_visis*sizeof(user_precision_t), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(sum_visi_XY_real, d_sum_visi_XY_real,
                             num_visis*sizeof(user_precision_t), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(sum_visi_YX_real, d_sum_visi_YX_real,
                             num_visis*sizeof(user_precision_t), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(sum_visi_YY_real, d_sum_visi_YY_real,
                             num_visis*sizeof(user_precision_t), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(sum_visi_XX_imag, d_sum_visi_XX_imag,
                             num_visis*sizeof(user_precision_t), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(sum_visi_XY_imag, d_sum_visi_XY_imag,
                             num_visis*sizeof(user_precision_t), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(sum_visi_YX_imag, d_sum_visi_YX_imag,
                             num_visis*sizeof(user_precision_t), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(sum_visi_YY_imag, d_sum_visi_YY_imag,
                             num_visis*sizeof(user_precision_t), cudaMemcpyDeviceToHost ));


  cudaErrorCheckCall(  cudaFree( d_sum_visi_XX_real ) );
  cudaErrorCheckCall(  cudaFree( d_sum_visi_XX_imag ) );
  cudaErrorCheckCall(  cudaFree( d_sum_visi_XY_real ) );
  cudaErrorCheckCall(  cudaFree( d_sum_visi_XY_imag ) );
  cudaErrorCheckCall(  cudaFree( d_sum_visi_YX_real ) );
  cudaErrorCheckCall(  cudaFree( d_sum_visi_YX_imag ) );
  cudaErrorCheckCall(  cudaFree( d_sum_visi_YY_real ) );
  cudaErrorCheckCall(  cudaFree( d_sum_visi_YY_imag ) );
  cudaErrorCheckCall(  cudaFree( d_allsteps_wavelengths ) );

  free_d_components(d_components, comptype);
  free_beam_gains(d_beam_gains, beamtype);

  cudaErrorCheckCall(  cudaFree( d_us ) );
  cudaErrorCheckCall(  cudaFree( d_vs ) );
  cudaErrorCheckCall(  cudaFree( d_ws ) );

  if (comptype == SHAPELET){
    cudaErrorCheckCall(  cudaFree( d_sbf) );
    cudaErrorCheckCall(  cudaFree( d_u_shapes) );
    cudaErrorCheckCall(  cudaFree( d_v_shapes) );
    cudaErrorCheckCall(  cudaFree( d_w_shapes) );
  }
}
