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
#include "FEE_primary_beam_cuda.h"

__device__ void extrap_stokes(float *d_allsteps_wavelengths, float *d_ref_freqs,
           float *d_ref_stokesI, float *d_ref_stokesQ,
           float *d_ref_stokesU, float *d_ref_stokesV,
           float *d_SIs, int iComponent, int iBaseline,
           float * flux_I, float * flux_Q, float * flux_U, float * flux_V){

  float d_freq = VELC / d_allsteps_wavelengths[iBaseline];
  float d_ref_freq = d_ref_freqs[iComponent];

  float flux_ratio = powf(d_freq / d_ref_freq, d_SIs[iComponent]);

  * flux_I = d_ref_stokesI[iComponent] * flux_ratio;
  * flux_Q = d_ref_stokesQ[iComponent] * flux_ratio;
  * flux_U = d_ref_stokesU[iComponent] * flux_ratio;
  * flux_V = d_ref_stokesV[iComponent] * flux_ratio;

}

__device__  cuFloatComplex calc_measurement_equation(float *d_us,
           float *d_vs, float *d_ws, float *d_ls, float *d_ms, float *d_ns,
           const int iBaseline, const int iComponent){

  float u, v, w;
  float l, m, n;

  u = d_us[iBaseline];
  v = d_vs[iBaseline];
  w = d_ws[iBaseline];

  l = d_ls[iComponent];
  m = d_ms[iComponent];
  n = d_ns[iComponent];

  cuFloatComplex visi;

  //Not sure why, but get match with OSKAR/RTS sims, and correct location
  //on sky through WSClean, without negative infront on 2pi
  float temp = 2*M_PI*( u*l + v*m + w*(n-1) );
  // sincosf(temp, &(visi.y), &(visi.x));

  visi.y = sinf(temp);
  visi.x = cosf(temp);

  return visi;
}

__device__ void apply_beam_gains(cuFloatComplex g1x, cuFloatComplex D1x,
          cuFloatComplex D1y, cuFloatComplex g1y,
          cuFloatComplex g2x, cuFloatComplex D2x,
          cuFloatComplex D2y, cuFloatComplex g2y,
          float flux_I, float flux_Q,
          float flux_U, float flux_V,
          cuFloatComplex visi_component,
          cuFloatComplex * visi_XX, cuFloatComplex * visi_XY,
          cuFloatComplex * visi_YX, cuFloatComplex * visi_YY) {

  //Conjugate the second beam gains
  cuFloatComplex g2x_conj = make_cuFloatComplex(g2x.x,-g2x.y);
  cuFloatComplex D2x_conj = make_cuFloatComplex(D2x.x,-D2x.y);
  cuFloatComplex D2y_conj = make_cuFloatComplex(D2y.x,-D2y.y);
  cuFloatComplex g2y_conj = make_cuFloatComplex(g2y.x,-g2y.y);

  //Create the Stokes visibilities
  cuFloatComplex visi_I = cuCmulf(make_cuComplex(flux_I,0.0),visi_component );
  cuFloatComplex visi_Q = cuCmulf(make_cuComplex(flux_Q,0.0),visi_component );
  cuFloatComplex visi_U = cuCmulf(make_cuComplex(flux_U,0.0),visi_component );
  cuFloatComplex visi_V = cuCmulf(make_cuComplex(flux_V,0.0),visi_component );

  cuFloatComplex this_XX;
  cuFloatComplex this_XY;
  cuFloatComplex this_YX;
  cuFloatComplex this_YY;

  //Convert the Stokes into instrumental visibilities
  this_XX = cuCmulf(cuCmulf(g1x,g2x_conj) + cuCmulf(D1x,D2x_conj),visi_I);
  this_XX += cuCmulf(cuCmulf(g1x,g2x_conj) - cuCmulf(D1x,D2x_conj),visi_Q);
  this_XX += cuCmulf(cuCmulf(g1x,D2x_conj) + cuCmulf(D1x,g2x_conj),visi_U);
  this_XX += cuCmulf(cuCmulf(make_cuFloatComplex(0.0,1.0),visi_V), cuCmulf(g1x,D2x_conj) - cuCmulf(D1x,g2x_conj) );

  this_XY = cuCmulf(cuCmulf(g1x,D2y_conj) + cuCmulf(D1x,g2y_conj),visi_I);
  this_XY += cuCmulf(cuCmulf(g1x,D2y_conj) - cuCmulf(D1x,g2y_conj),visi_Q);
  this_XY += cuCmulf(cuCmulf(g1x,g2y_conj) + cuCmulf(D1x,D2y_conj),visi_U);
  this_XY += cuCmulf(cuCmulf(make_cuFloatComplex(0.0,1.0),visi_V), cuCmulf(g1x,g2y_conj) - cuCmulf(D1x,D2y_conj) );

  this_YX = cuCmulf(cuCmulf(D1y,g2x_conj) + cuCmulf(g1y,D2x_conj),visi_I);
  this_YX += cuCmulf(cuCmulf(D1y,g2x_conj) - cuCmulf(g1y,D2x_conj),visi_Q);
  this_YX += cuCmulf(cuCmulf(D1y,D2x_conj) + cuCmulf(g1y,g2x_conj),visi_U);
  this_YX += cuCmulf(cuCmulf(make_cuFloatComplex(0.0,1.0),visi_V), cuCmulf(D1y,D2x_conj) - cuCmulf(g1y,g2x_conj) );

  this_YY = cuCmulf(cuCmulf(D1y,D2y_conj) + cuCmulf(g1y,g2y_conj),visi_I);
  this_YY += cuCmulf(cuCmulf(D1y,D2y_conj) - cuCmulf(g1y,g2y_conj),visi_Q);
  this_YY += cuCmulf(cuCmulf(D1y,g2y_conj) + cuCmulf(g1y,D2y_conj),visi_U);
  this_YY += cuCmulf(cuCmulf(make_cuFloatComplex(0.0,1.0),visi_V), cuCmulf(D1y,g2y_conj) - cuCmulf(g1y,D2y_conj) );

  * visi_XX = this_XX;
  * visi_XY = this_XY;
  * visi_YX = this_YX;
  * visi_YY = this_YY;

}

__device__ void get_beam_gains(int iBaseline, int iComponent, int num_freqs,
           int num_baselines, int num_components, int num_times, int beamtype,
           cuFloatComplex *d_primay_beam_J00, cuFloatComplex *d_primay_beam_J01,
           cuFloatComplex *d_primay_beam_J10, cuFloatComplex *d_primay_beam_J11,
           cuFloatComplex * g1x, cuFloatComplex * D1x,
           cuFloatComplex * D1y, cuFloatComplex * g1y,
           cuFloatComplex * g2x, cuFloatComplex * D2x,
           cuFloatComplex * D2y, cuFloatComplex * g2y){

  int beam_ind = 0;
  int time_ind = 0;
  int freq_ind = 0;

  time_ind = (int)floorf( (float)iBaseline / ((float)num_baselines * (float)num_freqs));
  freq_ind = (int)floorf( ((float)iBaseline - ((float)time_ind*(float)num_baselines * (float)num_freqs)) / (float)num_baselines);
  beam_ind = num_freqs*time_ind*num_components + (num_components*freq_ind) + iComponent;

  //Get XX,YY if using a beam
  if (beamtype == FEE_BEAM || beamtype == ANALY_DIPOLE || beamtype == GAUSS_BEAM) {
    * g1x = d_primay_beam_J00[beam_ind];
    * g2x = d_primay_beam_J00[beam_ind];
    * g1y = d_primay_beam_J11[beam_ind];
    * g2y = d_primay_beam_J11[beam_ind];
  }
  else {
    * g1x = make_cuComplex(1.0, 0.0);
    * g2x = make_cuComplex(1.0, 0.0);
    * g1y = make_cuComplex(1.0, 0.0);
    * g2y = make_cuComplex(1.0, 0.0);
  }

  //Only FEE model has XY and YX at the moment
  if (beamtype == FEE_BEAM) {
    * D1x = d_primay_beam_J01[beam_ind];
    * D2x = d_primay_beam_J01[beam_ind];
    * D1y = d_primay_beam_J10[beam_ind];
    * D2y = d_primay_beam_J10[beam_ind];
  }
  else {
    * D1x = make_cuComplex(0.0, 0.0);
    * D2x = make_cuComplex(0.0, 0.0);
    * D1y = make_cuComplex(0.0, 0.0);
    * D2y = make_cuComplex(0.0, 0.0);
  }
} //end __device__ get_beam_gains

__device__ void update_sum_visis(int iBaseline, int iComponent, int num_freqs,
           int num_baselines, int num_components, int num_times, int beamtype,
           cuFloatComplex *d_primay_beam_J00, cuFloatComplex *d_primay_beam_J01,
           cuFloatComplex *d_primay_beam_J10, cuFloatComplex *d_primay_beam_J11,
           cuFloatComplex visi_component,
           float flux_I, float flux_Q, float flux_U, float flux_V,
           float *d_sum_visi_XX_real, float *d_sum_visi_XX_imag,
           float *d_sum_visi_XY_real, float *d_sum_visi_XY_imag,
           float *d_sum_visi_YX_real, float *d_sum_visi_YX_imag,
           float *d_sum_visi_YY_real, float *d_sum_visi_YY_imag){

    cuFloatComplex g1x;
    cuFloatComplex D1x;
    cuFloatComplex D1y;
    cuFloatComplex g1y;
    cuFloatComplex g2x;
    cuFloatComplex D2x;
    cuFloatComplex D2y;
    cuFloatComplex g2y;

    get_beam_gains(iBaseline, iComponent, num_freqs,
               num_baselines, num_components, num_times, beamtype,
               d_primay_beam_J00, d_primay_beam_J01,
               d_primay_beam_J10, d_primay_beam_J11,
               &g1x, &D1x, &D1y, &g1y, &g2x, &D2x, &D2y, &g2y);




    cuFloatComplex visi_XX;
    cuFloatComplex visi_XY;
    cuFloatComplex visi_YX;
    cuFloatComplex visi_YY;

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

void source_component_common(int num_components,
           cuFloatComplex *d_primay_beam_J00, cuFloatComplex *d_primay_beam_J01,
           cuFloatComplex *d_primay_beam_J10, cuFloatComplex *d_primay_beam_J11,
           float *d_freqs, float *d_ls, float *d_ms, float *d_ns,
           float *d_ras, float *d_decs, float *azs, float *zas,
           float *sin_para_angs, float *cos_para_angs,
           float *beam_has, float *beam_decs,
           woden_settings_t *woden_settings,
           beam_settings_t beam_settings,
           RTS_MWA_FEE_beam_t *FEE_beam){

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
                        d_ras, d_decs,
                        d_ls, d_ms, d_ns, num_components)

  //If using a gaussian primary beam, calculate beam values for all freqs,
  //lsts and point component locations
  if (beam_settings.beamtype == GAUSS_BEAM) {

    //TODO currently hardcoded to have beam position angle = 0.
    //Should this change with az/za?
    float cos_theta = 1.0;
    float sin_theta = 0.0;
    float sin_2theta = 0.0;
    float fwhm_lm = sinf(beam_settings.beam_FWHM_rad);

    printf("\tDoing gaussian beam tings\n");

    calculate_gaussian_beam(num_components,
         woden_settings->num_time_steps, woden_settings->num_freqs,
         beam_settings.gauss_ha, beam_settings.gauss_sdec,
         beam_settings.gauss_cdec,
         fwhm_lm, cos_theta, sin_theta, sin_2theta,
         beam_settings.beam_ref_freq, d_freqs,
         beam_has, beam_decs,
         d_primay_beam_J00, d_primay_beam_J11);

  }// end if beam == GAUSS

  else if (beam_settings.beamtype == FEE_BEAM) {

    //Rotate FEE beam by parallactic angle
    int rotation = 1;
    //Normalise FEE beam to zenith
    int scaling = 1;

    calc_CUDA_FEE_beam(azs, zas, sin_para_angs, cos_para_angs,
           num_components, woden_settings->num_time_steps, FEE_beam,
           rotation, scaling);

    threads.x = 16;
    threads.y = 16;
    grid.x = (int)ceil( ((float)num_components) / ((float)threads.x) );
    grid.y = (int)ceil( (float)woden_settings->num_time_steps / (float)threads.y );


    cudaErrorCheckKernel("kern_map_FEE_beam_gains",
              kern_map_FEE_beam_gains, grid, threads,
              (cuFloatComplex *)FEE_beam->d_FEE_beam_gain_matrices,
              d_primay_beam_J00, d_primay_beam_J01,
              d_primay_beam_J10, d_primay_beam_J11,
              woden_settings->num_freqs, num_components,
              woden_settings->num_visis, woden_settings->num_baselines,
              woden_settings->num_time_steps);
  }

  else if (beam_settings.beamtype == ANALY_DIPOLE) {
    printf("\tDoing analytic_dipole (EDA2 beam)\n");

    calculate_analytic_dipole_beam(num_components,
         woden_settings->num_time_steps, woden_settings->num_freqs,
         azs, zas, d_freqs, d_primay_beam_J00, d_primay_beam_J11);
  }
}

__global__ void kern_calc_visi_point(float *d_point_ras, float *d_point_decs,
           float *d_point_freqs, float *d_point_stokesI, float *d_point_stokesQ,
           float *d_point_stokesU, float *d_point_stokesV, float *d_point_SIs,
           float *d_us, float *d_vs, float *d_ws,
           float *d_sum_visi_XX_real, float *d_sum_visi_XX_imag,
           float *d_sum_visi_XY_real, float *d_sum_visi_XY_imag,
           float *d_sum_visi_YX_real, float *d_sum_visi_YX_imag,
           float *d_sum_visi_YY_real, float *d_sum_visi_YY_imag,
           float *d_allsteps_wavelengths,
           float *d_ls, float *d_ms, float *d_ns,
           int num_points, int num_baselines, int num_freqs, int num_visis,
           int num_times, int beamtype,
           cuFloatComplex *d_primay_beam_J00, cuFloatComplex *d_primay_beam_J01,
           cuFloatComplex *d_primay_beam_J10, cuFloatComplex *d_primay_beam_J11) {

  // Start by computing which baseline we're going to do
  const int iBaseline = threadIdx.x + (blockDim.x*blockIdx.x);
  // const int iComponent = threadIdx.y + (blockDim.y*blockIdx.y);
  // if(iBaseline < num_visis && iComponent < num_points) {
  if(iBaseline < num_visis) {

    float point_flux_I;
    float point_flux_Q;
    float point_flux_U;
    float point_flux_V;

    cuFloatComplex visi_point;

    for (size_t iComponent = 0; iComponent < num_points; iComponent++) {

      extrap_stokes(d_allsteps_wavelengths, d_point_freqs,
                   d_point_stokesI, d_point_stokesQ,
                   d_point_stokesU, d_point_stokesV,
                   d_point_SIs, iComponent, iBaseline,
                   &point_flux_I, &point_flux_Q, &point_flux_U, &point_flux_V);


      visi_point = calc_measurement_equation(d_us, d_vs, d_ws,
                             d_ls, d_ms, d_ns,
                             iBaseline, iComponent);

      update_sum_visis(iBaseline, iComponent, num_freqs,
             num_baselines, num_points, num_times, beamtype,
             d_primay_beam_J00, d_primay_beam_J01,
             d_primay_beam_J10, d_primay_beam_J11,
             visi_point,
             point_flux_I, point_flux_Q, point_flux_U, point_flux_V,
             d_sum_visi_XX_real, d_sum_visi_XX_imag,
             d_sum_visi_XY_real, d_sum_visi_XY_imag,
             d_sum_visi_YX_real, d_sum_visi_YX_imag,
             d_sum_visi_YY_real, d_sum_visi_YY_imag);
    }
  }
}

__global__ void kern_calc_visi_gaussian(float *d_gauss_ras, float *d_gauss_decs,
           float *d_gauss_freqs, float *d_gauss_stokesI, float *d_gauss_stokesQ,
           float *d_gauss_stokesU, float *d_gauss_stokesV, float *d_gauss_SIs,
           float *d_us, float *d_vs, float *d_ws,
           float *d_sum_visi_XX_real, float *d_sum_visi_XX_imag,
           float *d_sum_visi_XY_real, float *d_sum_visi_XY_imag,
           float *d_sum_visi_YX_real, float *d_sum_visi_YX_imag,
           float *d_sum_visi_YY_real, float *d_sum_visi_YY_imag,
           float *d_allsteps_wavelengths,
           float *d_ls, float *d_ms, float *d_ns,
           float *d_gauss_pas, float *d_gauss_majors, float *d_gauss_minors,
           int num_gauss, int num_baselines, int num_freqs, int num_visis,
           int num_times, int beamtype,
           cuFloatComplex *d_primay_beam_J00, cuFloatComplex *d_primay_beam_J01,
           cuFloatComplex *d_primay_beam_J10, cuFloatComplex *d_primay_beam_J11) {

  // Start by computing which baseline we're going to do
  const int iBaseline = threadIdx.x + (blockDim.x*blockIdx.x);
  // const int iComponent = threadIdx.y + (blockDim.y*blockIdx.y);

  // if(iBaseline < num_visis && iComponent < num_gauss) {
  if(iBaseline < num_visis) {

    float gauss_flux_I;
    float gauss_flux_Q;
    float gauss_flux_U;
    float gauss_flux_V;
    cuFloatComplex visi_gauss;

    for (size_t iComponent = 0; iComponent < num_gauss; iComponent++) {

      extrap_stokes(d_allsteps_wavelengths, d_gauss_freqs,
                   d_gauss_stokesI, d_gauss_stokesQ,
                   d_gauss_stokesU, d_gauss_stokesV,
                   d_gauss_SIs, iComponent, iBaseline,
                   &gauss_flux_I, &gauss_flux_Q, &gauss_flux_U, &gauss_flux_V);


      visi_gauss = calc_measurement_equation(d_us, d_vs, d_ws,
                             d_ls, d_ms, d_ns,
                             iBaseline, iComponent);

      cuFloatComplex V_envelop = make_cuFloatComplex( 1.0, 0.0 );

      float pa = d_gauss_pas[iComponent];
      float u = d_us[iBaseline];
      float v = d_vs[iBaseline];
      float sinpa = sin(pa);
      float cospa = cos(pa);

      float x =  cospa*v + sinpa*u; // major axis
      float y = -sinpa*v + cospa*u; // minor axis
      float invsig_x = d_gauss_majors[iComponent];
      float invsig_y = d_gauss_minors[iComponent];

      V_envelop = make_cuFloatComplex( exp( -0.5 * ( x*x*invsig_x*invsig_x*M_PI_2_2_LN_2 + y*y*invsig_y*invsig_y*M_PI_2_2_LN_2 ) ), 0.0 );

      visi_gauss = cuCmulf(visi_gauss, V_envelop);

      update_sum_visis(iBaseline, iComponent, num_freqs,
             num_baselines, num_gauss, num_times, beamtype,
             d_primay_beam_J00, d_primay_beam_J01,
             d_primay_beam_J10, d_primay_beam_J11,
             visi_gauss,
             gauss_flux_I, gauss_flux_Q, gauss_flux_U, gauss_flux_V,
             d_sum_visi_XX_real, d_sum_visi_XX_imag,
             d_sum_visi_XY_real, d_sum_visi_XY_imag,
             d_sum_visi_YX_real, d_sum_visi_YX_imag,
             d_sum_visi_YY_real, d_sum_visi_YY_imag);
    }
  }
}

__global__ void kern_calc_visi_shapelets(float *d_shape_freqs,
      float *d_shape_stokesI, float *d_shape_stokesQ,
      float *d_shape_stokesU, float *d_shape_stokesV, float *d_shape_SIs,
      float *d_us, float *d_vs, float *d_ws,
      float *d_allsteps_wavelengths,
      float *d_u_shapes, float *d_v_shapes, float *d_w_shapes,
      float *d_sum_visi_XX_real, float *d_sum_visi_XX_imag,
      float *d_sum_visi_XY_real, float *d_sum_visi_XY_imag,
      float *d_sum_visi_YX_real, float *d_sum_visi_YX_imag,
      float *d_sum_visi_YY_real, float *d_sum_visi_YY_imag,
      float *d_shape_pas, float *d_shape_majors,
      float *d_shape_minors,
      float *d_shape_n1s, float *d_shape_n2s, float *d_shape_coeffs,
      float *d_shape_param_indexes,
      float *d_shape_ls, float *d_shape_ms, float *d_shape_ns,
      float *d_sbf,
      int num_shapes, int num_baselines, int num_freqs, int num_visis,
      const int num_coeffs, int num_times, int beamtype,
      cuFloatComplex *d_primay_beam_J00, cuFloatComplex *d_primay_beam_J01,
      cuFloatComplex *d_primay_beam_J10, cuFloatComplex *d_primay_beam_J11) {

  // Start by computing which baseline we're going to do
  const int iBaseline = threadIdx.x + (blockDim.x*blockIdx.x);
  // const int iCoeff = threadIdx.y + (blockDim.y*blockIdx.y);

  // if (iBaseline < num_visis && iCoeff < num_coeffs) {
  if (iBaseline < num_visis) {

    float shape_flux_I;
    float shape_flux_Q;
    float shape_flux_U;
    float shape_flux_V;
    cuFloatComplex visi_shape;

    for (int iCoeff = 0; iCoeff < num_coeffs; iCoeff++) {

      int iComponent = d_shape_param_indexes[iCoeff];

      extrap_stokes(d_allsteps_wavelengths, d_shape_freqs,
                   d_shape_stokesI, d_shape_stokesQ,
                   d_shape_stokesU, d_shape_stokesV,
                   d_shape_SIs, iComponent, iBaseline,
                   &shape_flux_I, &shape_flux_Q, &shape_flux_U, &shape_flux_V);


      visi_shape = calc_measurement_equation(d_us, d_vs, d_ws,
                            d_shape_ls, d_shape_ms, d_shape_ns,
                            iBaseline, iComponent);

      float pa = d_shape_pas[iComponent];
      float sinpa = sin(pa);
      float cospa = cos(pa);

      float u_shape = d_u_shapes[iComponent*num_visis + iBaseline];
      float v_shape = d_v_shapes[iComponent*num_visis + iBaseline];

      float x = (cospa*v_shape + sinpa*u_shape); // major axis
      float y = (-sinpa*v_shape + cospa*u_shape); // minor axis

      //Scales the FWHM to std to match basis functions, and account for the
      //basis functions being stored with beta = 1.0
      //Basis functions have been stored in such a way that x is in the same
      //direction as on sky, but y is opposite, so include negative here
      float const_x = (d_shape_majors[iComponent]*SQRT_M_PI_2_2_LN_2)/sbf_dx;
      float const_y = -(d_shape_minors[iComponent]*SQRT_M_PI_2_2_LN_2)/sbf_dx;

      // I^(n1+n2) = Ipow_lookup[(n1+n2) % 4]
      cuFloatComplex Ipow_lookup[] = { make_cuFloatComplex(  1.0,  0.0 ),
                                       make_cuFloatComplex(  0.0,  1.0 ),
                                       make_cuFloatComplex( -1.0,  0.0 ),
                                       make_cuFloatComplex(  0.0, -1.0 ) };

      float xlow, xhigh, ylow, yhigh, u_value, v_value, f_hat, *sbf_n;

      // find the indices in the basis functions for u*beta_u and v*beta_v

      float xpos = x*const_x + sbf_c;
      float ypos = y*const_y + sbf_c;

      int xindex = (int)floor(xpos);
      int yindex = (int)floor(ypos);
      //
      int n1 = (int)d_shape_n1s[iCoeff];
      int n2 = (int)d_shape_n2s[iCoeff];

      // if ( n1 < 0 || n2 < 0 || n1 >= sbf_N || n2 >= sbf_N ) continue;

      f_hat = d_shape_coeffs[iCoeff];
      //
      sbf_n = &d_sbf[n1*sbf_L];
      xlow  = sbf_n[xindex];
      xhigh = sbf_n[xindex+1];
      u_value = xlow + (xhigh-xlow)*(xpos-xindex);

      sbf_n = &d_sbf[n2*sbf_L];
      ylow  = sbf_n[yindex];
      yhigh = sbf_n[yindex+1];
      v_value = ylow + (yhigh-ylow)*(ypos-yindex);

      // accumulate the intensity model for baseline pair (u,v)
      cuFloatComplex V_envelop = make_cuFloatComplex( 0.0, 0.0 );
      V_envelop = V_envelop + Ipow_lookup[(n1+n2) % 4] * f_hat * u_value*v_value;

      visi_shape = cuCmulf(visi_shape, V_envelop);

      update_sum_visis(iBaseline, iComponent, num_freqs,
             num_baselines, num_shapes, num_times, beamtype,
             d_primay_beam_J00, d_primay_beam_J01,
             d_primay_beam_J10, d_primay_beam_J11,
             visi_shape,
             shape_flux_I, shape_flux_Q, shape_flux_U, shape_flux_V,
             d_sum_visi_XX_real, d_sum_visi_XX_imag,
             d_sum_visi_XY_real, d_sum_visi_XY_imag,
             d_sum_visi_YX_real, d_sum_visi_YX_imag,
             d_sum_visi_YY_real, d_sum_visi_YY_imag);
    }
  }
}

/*******************************************************************************
                 Functions below to be used in unit tests
*******************************************************************************/

__global__ void kern_extrap_stokes(int num_extrap_freqs, int num_components,
           float *d_extrap_wavelengths, float *d_ref_freqs, float *d_SIs,
           float *d_ref_stokesI, float *d_ref_stokesQ,
           float *d_ref_stokesU, float *d_ref_stokesV,
           float *d_flux_I, float *d_flux_Q,
           float *d_flux_U, float *d_flux_V ) {

  // Start by computing which baseline we're going to do
  const int iComponent = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iFreq = threadIdx.y + (blockDim.y*blockIdx.y);
  // if(iBaseline < num_visis && iComponent < num_points) {
  if(iComponent < num_components && iFreq < num_extrap_freqs) {

    float flux_I;
    float flux_Q;
    float flux_U;
    float flux_V;

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
           float *extrap_wavelengths, float *ref_freqs, float *SIs,
           float *ref_stokesI, float *ref_stokesQ,
           float *ref_stokesU, float *ref_stokesV,
           float *flux_I, float *flux_Q,
           float *flux_U, float *flux_V){

  float *d_extrap_wavelengths = NULL;
  float *d_ref_freqs = NULL;
  float *d_SIs = NULL;
  float *d_ref_stokesI = NULL;
  float *d_ref_stokesQ = NULL;
  float *d_ref_stokesU = NULL;
  float *d_ref_stokesV = NULL;
  float *d_flux_I = NULL;
  float *d_flux_Q = NULL;
  float *d_flux_U = NULL;
  float *d_flux_V = NULL;

  cudaErrorCheckCall( cudaMalloc( (void**)&d_extrap_wavelengths, num_extrap_freqs*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_ref_freqs, num_components*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_SIs, num_components*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_ref_stokesI, num_components*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_ref_stokesQ, num_components*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_ref_stokesU, num_components*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_ref_stokesV, num_components*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_flux_I, num_components*num_extrap_freqs*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_flux_Q, num_components*num_extrap_freqs*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_flux_U, num_components*num_extrap_freqs*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_flux_V, num_components*num_extrap_freqs*sizeof(float) ));

  cudaErrorCheckCall( cudaMemcpy(d_extrap_wavelengths, extrap_wavelengths, num_extrap_freqs*sizeof(float), cudaMemcpyHostToDevice ));

  cudaErrorCheckCall( cudaMemcpy(d_ref_stokesI, ref_stokesI, num_components*sizeof(float), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_ref_stokesQ, ref_stokesQ, num_components*sizeof(float), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_ref_stokesU, ref_stokesU, num_components*sizeof(float), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_ref_stokesV, ref_stokesV, num_components*sizeof(float), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_ref_freqs, ref_freqs, num_components*sizeof(float), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_SIs, SIs, num_components*sizeof(float), cudaMemcpyHostToDevice ));

  dim3 grid, threads;

  threads.x = 16;
  threads.y = 16;
  grid.x = (int)ceil( (float)num_components / (float)threads.x );
  grid.y = (int)ceil( (float)num_extrap_freqs / (float)threads.y );

  cudaErrorCheckKernel("kern_extrap_stokes",
                        kern_extrap_stokes, grid, threads,
                        num_extrap_freqs, num_components,
                        d_extrap_wavelengths, d_ref_freqs, d_SIs,
                        d_ref_stokesI, d_ref_stokesQ,
                        d_ref_stokesU, d_ref_stokesV,
                        d_flux_I, d_flux_Q,
                        d_flux_U, d_flux_V);

  cudaErrorCheckCall( cudaMemcpy(flux_I, d_flux_I, num_components*num_extrap_freqs*sizeof(float),cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(flux_Q, d_flux_Q, num_components*num_extrap_freqs*sizeof(float),cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(flux_U, d_flux_U, num_components*num_extrap_freqs*sizeof(float),cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(flux_V, d_flux_V, num_components*num_extrap_freqs*sizeof(float),cudaMemcpyDeviceToHost ));

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
                  float *d_us, float *d_vs, float *d_ws,
                  float *d_ls, float *d_ms, float *d_ns,
                  cuFloatComplex *d_visis) {

  // Start by computing which baseline we're going to do
  const int iBaseline = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iComponent = threadIdx.y + (blockDim.y*blockIdx.y);
  // if(iBaseline < num_visis && iComponent < num_points) {
  if(iComponent < num_components && iBaseline < num_baselines) {

    cuFloatComplex visi;
    visi = calc_measurement_equation(d_us, d_vs, d_ws, d_ls, d_ms, d_ns,
                                     iBaseline, iComponent);

    int visi_ind = num_components*iBaseline + iComponent;
    d_visis[visi_ind] = visi;

  }
}

extern "C" void test_kern_calc_measurement_equation(int num_components,
                  int num_baselines,
                  float *us, float *vs, float *ws,
                  float *ls, float *ms, float *ns,
                  float _Complex *visis){

  float *d_us = NULL;
  float *d_vs = NULL;
  float *d_ws = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_us, num_baselines*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_vs, num_baselines*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_ws, num_baselines*sizeof(float) ));
  cudaErrorCheckCall( cudaMemcpy(d_us, us, num_baselines*sizeof(float), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_vs, vs, num_baselines*sizeof(float), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_ws, ws, num_baselines*sizeof(float), cudaMemcpyHostToDevice ));

  float *d_ls = NULL;
  float *d_ms = NULL;
  float *d_ns = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_ls, num_components*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_ms, num_components*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_ns, num_components*sizeof(float) ));
  cudaErrorCheckCall( cudaMemcpy(d_ls, ls, num_components*sizeof(float), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_ms, ms, num_components*sizeof(float), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_ns, ns, num_components*sizeof(float), cudaMemcpyHostToDevice ));

  float _Complex *d_visis = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_visis, num_baselines*num_components*sizeof(float _Complex) ));

  dim3 grid, threads;

  threads.x = 16;
  threads.y = 16;
  grid.x = (int)ceil( (float)num_baselines / (float)threads.x );
  grid.y = (int)ceil( (float)num_components / (float)threads.y );

  cudaErrorCheckKernel("kern_calc_measurement_equation",
                      kern_calc_measurement_equation, grid, threads,
                      num_components, num_baselines,
                      d_us, d_vs, d_ws,
                      d_ls, d_ms, d_ns,
                      (cuFloatComplex*)d_visis );

  cudaErrorCheckCall( cudaMemcpy(visis, (float _Complex*)d_visis, num_components*num_baselines*sizeof(float _Complex),cudaMemcpyDeviceToHost ));

  cudaErrorCheckCall( cudaFree( d_us ) );
  cudaErrorCheckCall( cudaFree( d_vs ) );
  cudaErrorCheckCall( cudaFree( d_ws ) );
  cudaErrorCheckCall( cudaFree( d_ls ) );
  cudaErrorCheckCall( cudaFree( d_ms ) );
  cudaErrorCheckCall( cudaFree( d_ns ) );
  cudaErrorCheckCall( cudaFree(d_visis ) );

}

__global__ void kern_apply_beam_gains(int num_gains, cuFloatComplex *d_g1xs,
          cuFloatComplex *d_D1xs,
          cuFloatComplex *d_D1ys, cuFloatComplex *d_g1ys,
          cuFloatComplex *d_g2xs, cuFloatComplex *d_D2xs,
          cuFloatComplex *d_D2ys, cuFloatComplex *d_g2ys,
          float *d_flux_Is, float *d_flux_Qs,
          float *d_flux_Us, float *d_flux_Vs,
          cuFloatComplex *d_visi_components,
          cuFloatComplex *d_visi_XXs, cuFloatComplex *d_visi_XYs,
          cuFloatComplex *d_visi_YXs, cuFloatComplex *d_visi_YYs) {

  const int iGain = threadIdx.x + (blockDim.x*blockIdx.x);
  // const int iComponent = threadIdx.y + (blockDim.y*blockIdx.y);
  // if(iBaseline < num_visis && iComponent < num_points) {
  if (iGain < num_gains) {

    cuFloatComplex visi_XX;
    cuFloatComplex visi_XY;
    cuFloatComplex visi_YX;
    cuFloatComplex visi_YY;

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

extern "C" void test_kern_apply_beam_gains(int num_gains, float _Complex *g1xs,
          float _Complex *D1xs,
          float _Complex *D1ys, float _Complex *g1ys,
          float _Complex *g2xs, float _Complex *D2xs,
          float _Complex *D2ys, float _Complex *g2ys,
          float *flux_Is, float *flux_Qs,
          float *flux_Us, float *flux_Vs,
          float _Complex *visi_components,
          float _Complex *visi_XXs, float _Complex *visi_XYs,
          float _Complex *visi_YXs, float _Complex *visi_YYs){

  float _Complex *d_g1xs = NULL;
  float _Complex *d_D1xs = NULL;
  float _Complex *d_D1ys = NULL;
  float _Complex *d_g1ys = NULL;
  float _Complex *d_g2xs = NULL;
  float _Complex *d_D2xs = NULL;
  float _Complex *d_D2ys = NULL;
  float _Complex *d_g2ys = NULL;
  float *d_flux_Is = NULL;
  float *d_flux_Qs = NULL;
  float *d_flux_Us = NULL;
  float *d_flux_Vs = NULL;
  float _Complex *d_visi_components = NULL;
  float _Complex *d_visi_XXs = NULL;
  float _Complex *d_visi_XYs = NULL;
  float _Complex *d_visi_YXs = NULL;
  float _Complex *d_visi_YYs = NULL;

  cudaErrorCheckCall( cudaMalloc( (void**)&d_g1xs, num_gains*sizeof(float _Complex) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_D1xs, num_gains*sizeof(float _Complex) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_D1ys, num_gains*sizeof(float _Complex) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_g1ys, num_gains*sizeof(float _Complex) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_g2xs, num_gains*sizeof(float _Complex) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_D2xs, num_gains*sizeof(float _Complex) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_D2ys, num_gains*sizeof(float _Complex) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_g2ys, num_gains*sizeof(float _Complex) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_flux_Is, num_gains*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_flux_Qs, num_gains*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_flux_Us, num_gains*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_flux_Vs, num_gains*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_visi_components, num_gains*sizeof(float _Complex) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_visi_XXs, num_gains*sizeof(float _Complex) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_visi_XYs, num_gains*sizeof(float _Complex) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_visi_YXs, num_gains*sizeof(float _Complex) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_visi_YYs, num_gains*sizeof(float _Complex) ));

  cudaErrorCheckCall( cudaMemcpy(d_g1xs, g1xs, num_gains*sizeof(float _Complex), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_D1xs, D1xs, num_gains*sizeof(float _Complex), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_D1ys, D1ys, num_gains*sizeof(float _Complex), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_g1ys, g1ys, num_gains*sizeof(float _Complex), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_g2xs, g2xs, num_gains*sizeof(float _Complex), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_D2xs, D2xs, num_gains*sizeof(float _Complex), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_D2ys, D2ys, num_gains*sizeof(float _Complex), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_g2ys, g2ys, num_gains*sizeof(float _Complex), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_visi_components, visi_components, num_gains*sizeof(float _Complex), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_visi_XXs, visi_XXs, num_gains*sizeof(float _Complex), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_visi_XYs, visi_XYs, num_gains*sizeof(float _Complex), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_visi_YXs, visi_YXs, num_gains*sizeof(float _Complex), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_visi_YYs, visi_YYs, num_gains*sizeof(float _Complex), cudaMemcpyHostToDevice ));

  cudaErrorCheckCall( cudaMemcpy(d_flux_Is, flux_Is, num_gains*sizeof(float), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_flux_Qs, flux_Qs, num_gains*sizeof(float), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_flux_Us, flux_Us, num_gains*sizeof(float), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_flux_Vs, flux_Vs, num_gains*sizeof(float), cudaMemcpyHostToDevice ));

  dim3 grid, threads;

  threads.x = 128;
  grid.x = (int)ceil( (float)num_gains / (float)threads.x );

  cudaErrorCheckKernel("kern_apply_beam_gains",
                      kern_apply_beam_gains, grid, threads,
                      num_gains,
                      (cuFloatComplex *)d_g1xs, (cuFloatComplex *)d_D1xs,
                      (cuFloatComplex *)d_D1ys, (cuFloatComplex *)d_g1ys,
                      (cuFloatComplex *)d_g2xs, (cuFloatComplex *)d_D2xs,
                      (cuFloatComplex *)d_D2ys, (cuFloatComplex *)d_g2ys,
                      d_flux_Is, d_flux_Qs,
                      d_flux_Us, d_flux_Vs,
                      (cuFloatComplex *)d_visi_components,
                      (cuFloatComplex *)d_visi_XXs, (cuFloatComplex *)d_visi_XYs,
                      (cuFloatComplex *)d_visi_YXs, (cuFloatComplex *)d_visi_YYs );

  cudaErrorCheckCall( cudaMemcpy(visi_XXs, d_visi_XXs,
                      num_gains*sizeof(float _Complex),cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(visi_XYs, d_visi_XYs,
                      num_gains*sizeof(float _Complex),cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(visi_YXs, d_visi_YXs,
                      num_gains*sizeof(float _Complex),cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(visi_YYs, d_visi_YYs,
                      num_gains*sizeof(float _Complex),cudaMemcpyDeviceToHost ));

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
           cuFloatComplex *d_primay_beam_J00, cuFloatComplex *d_primay_beam_J01,
           cuFloatComplex *d_primay_beam_J10, cuFloatComplex *d_primay_beam_J11,
           cuFloatComplex *d_recov_g1x, cuFloatComplex *d_recov_D1x,
           cuFloatComplex *d_recov_D1y, cuFloatComplex *d_recov_g1y,
           cuFloatComplex *d_recov_g2x, cuFloatComplex *d_recov_D2x,
           cuFloatComplex *d_recov_D2y, cuFloatComplex *d_recov_g2y) {

  // Start by computing which baseline we're going to do
  const int iBaseline = threadIdx.x + (blockDim.x*blockIdx.x);
  // const int iComponent = threadIdx.y + (blockDim.y*blockIdx.y);
  // if(iBaseline < num_visis && iComponent < num_points) {
  if(iBaseline < num_visis) {

    for (size_t iComponent = 0; iComponent < num_components; iComponent++) {

      cuFloatComplex g1x;
      cuFloatComplex D1x;
      cuFloatComplex D1y;
      cuFloatComplex g1y;
      cuFloatComplex g2x;
      cuFloatComplex D2x;
      cuFloatComplex D2y;
      cuFloatComplex g2y;

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
          float _Complex *primay_beam_J00, float _Complex *primay_beam_J01,
          float _Complex *primay_beam_J10, float _Complex *primay_beam_J11,
          float _Complex *recover_g1x, float _Complex *recover_D1x,
          float _Complex *recover_D1y, float _Complex *recover_g1y,
          float _Complex *recover_g2x, float _Complex *recover_D2x,
          float _Complex *recover_D2y, float _Complex *recover_g2y){

  float _Complex *d_recover_g1x = NULL;
  float _Complex *d_recover_D1x = NULL;
  float _Complex *d_recover_D1y = NULL;
  float _Complex *d_recover_g1y = NULL;
  float _Complex *d_recover_g2x = NULL;
  float _Complex *d_recover_D2x = NULL;
  float _Complex *d_recover_D2y = NULL;
  float _Complex *d_recover_g2y = NULL;

  float _Complex *d_primay_beam_J00 = NULL;
  float _Complex *d_primay_beam_J01 = NULL;
  float _Complex *d_primay_beam_J10 = NULL;
  float _Complex *d_primay_beam_J11 = NULL;

  cudaErrorCheckCall( cudaMalloc( (void**)&d_recover_g1x, num_components*num_visis*sizeof(float _Complex) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_recover_D1x, num_components*num_visis*sizeof(float _Complex) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_recover_D1y, num_components*num_visis*sizeof(float _Complex) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_recover_g1y, num_components*num_visis*sizeof(float _Complex) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_recover_g2x, num_components*num_visis*sizeof(float _Complex) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_recover_D2x, num_components*num_visis*sizeof(float _Complex) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_recover_D2y, num_components*num_visis*sizeof(float _Complex) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_recover_g2y, num_components*num_visis*sizeof(float _Complex) ));

  cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J00, num_freqs*num_times*num_components*sizeof(float _Complex) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J01, num_freqs*num_times*num_components*sizeof(float _Complex) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J10, num_freqs*num_times*num_components*sizeof(float _Complex) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J11, num_freqs*num_times*num_components*sizeof(float _Complex) ));

  cudaErrorCheckCall( cudaMemcpy(d_primay_beam_J00, primay_beam_J00, num_freqs*num_times*num_components*sizeof(float _Complex), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_primay_beam_J01, primay_beam_J01, num_freqs*num_times*num_components*sizeof(float _Complex), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_primay_beam_J10, primay_beam_J10, num_freqs*num_times*num_components*sizeof(float _Complex), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_primay_beam_J11, primay_beam_J11, num_freqs*num_times*num_components*sizeof(float _Complex), cudaMemcpyHostToDevice ));

  dim3 grid, threads;

  threads.x = 128;
  grid.x = (int)ceil( (float)num_visis / (float)threads.x );

  cudaErrorCheckKernel("kern_get_beam_gains",
                      kern_get_beam_gains, grid, threads,
                      num_components, num_baselines,
                      num_freqs, num_visis, num_times, beamtype,
                      (cuFloatComplex *)d_primay_beam_J00,
                      (cuFloatComplex *)d_primay_beam_J01,
                      (cuFloatComplex *)d_primay_beam_J10,
                      (cuFloatComplex *)d_primay_beam_J11,
                      (cuFloatComplex *)d_recover_g1x, (cuFloatComplex *)d_recover_D1x,
                      (cuFloatComplex *)d_recover_D1y, (cuFloatComplex *)d_recover_g1y,
                      (cuFloatComplex *)d_recover_g2x, (cuFloatComplex *)d_recover_D2x,
                      (cuFloatComplex *)d_recover_D2y, (cuFloatComplex *)d_recover_g2y );

  cudaErrorCheckCall( cudaMemcpy(recover_g1x, d_recover_g1x, num_components*num_visis*sizeof(float _Complex), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(recover_D1x, d_recover_D1x, num_components*num_visis*sizeof(float _Complex), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(recover_D1y, d_recover_D1y, num_components*num_visis*sizeof(float _Complex), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(recover_g1y, d_recover_g1y, num_components*num_visis*sizeof(float _Complex), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(recover_g2x, d_recover_g2x, num_components*num_visis*sizeof(float _Complex), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(recover_D2x, d_recover_D2x, num_components*num_visis*sizeof(float _Complex), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(recover_D2y, d_recover_D2y, num_components*num_visis*sizeof(float _Complex), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(recover_g2y, d_recover_g2y, num_components*num_visis*sizeof(float _Complex), cudaMemcpyDeviceToHost ));

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
     cuFloatComplex *d_primay_beam_J00, cuFloatComplex *d_primay_beam_J01,
     cuFloatComplex *d_primay_beam_J10, cuFloatComplex *d_primay_beam_J11,
     cuFloatComplex *d_visi_components,
     float *d_flux_I, float *d_flux_Q, float *d_flux_U, float *d_flux_V,
     float *d_sum_visi_XX_real, float *d_sum_visi_XX_imag,
     float *d_sum_visi_XY_real, float *d_sum_visi_XY_imag,
     float *d_sum_visi_YX_real, float *d_sum_visi_YX_imag,
     float *d_sum_visi_YY_real, float *d_sum_visi_YY_imag) {

  // Start by computing which baseline we're going to do
  const int iBaseline = threadIdx.x + (blockDim.x*blockIdx.x);

  if(iBaseline < num_freqs*num_baselines*num_times) {

    int time_ind = (int)floorf( (float)iBaseline / ((float)num_baselines * (float)num_freqs));
    int freq_ind = (int)floorf( ((float)iBaseline - ((float)time_ind*(float)num_baselines * (float)num_freqs)) / (float)num_baselines);

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
          float _Complex *primay_beam_J00, float _Complex *primay_beam_J01,
          float _Complex *primay_beam_J10, float _Complex *primay_beam_J11,
          float _Complex *visi_components,
          float *flux_I, float *flux_Q, float *flux_U, float *flux_V,
          float *sum_visi_XX_real, float *sum_visi_XX_imag,
          float *sum_visi_XY_real, float *sum_visi_XY_imag,
          float *sum_visi_YX_real, float *sum_visi_YX_imag,
          float *sum_visi_YY_real, float *sum_visi_YY_imag){

  float _Complex *d_primay_beam_J00 = NULL;
  float _Complex *d_primay_beam_J01 = NULL;
  float _Complex *d_primay_beam_J10 = NULL;
  float _Complex *d_primay_beam_J11 = NULL;
  float _Complex *d_visi_components = NULL;

  cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J00,
                    num_components*num_times*num_freqs*sizeof(float _Complex) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J01,
                    num_components*num_times*num_freqs*sizeof(float _Complex) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J10,
                    num_components*num_times*num_freqs*sizeof(float _Complex) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J11,
                    num_components*num_times*num_freqs*sizeof(float _Complex) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_visi_components,
                    num_visis*sizeof(float _Complex) ));

  cudaErrorCheckCall( cudaMemcpy(d_primay_beam_J00, primay_beam_J00,
                      num_components*num_times*num_freqs*sizeof(float _Complex),
                      cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_primay_beam_J01, primay_beam_J01,
                      num_components*num_times*num_freqs*sizeof(float _Complex),
                      cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_primay_beam_J10, primay_beam_J10,
                      num_components*num_times*num_freqs*sizeof(float _Complex),
                      cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_primay_beam_J11, primay_beam_J11,
                      num_components*num_times*num_freqs*sizeof(float _Complex),
                      cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_visi_components, visi_components,
                     num_visis*sizeof(float _Complex), cudaMemcpyHostToDevice ));

  float *d_flux_I = NULL;
  float *d_flux_Q = NULL;
  float *d_flux_U = NULL;
  float *d_flux_V = NULL;

  cudaErrorCheckCall( cudaMalloc( (void**)&d_flux_I, num_components*num_times*num_freqs*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_flux_Q, num_components*num_times*num_freqs*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_flux_U, num_components*num_times*num_freqs*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_flux_V, num_components*num_times*num_freqs*sizeof(float) ));

  cudaErrorCheckCall( cudaMemcpy(d_flux_I, flux_I,
                             num_components*num_times*num_freqs*sizeof(float), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_flux_Q, flux_Q,
                             num_components*num_times*num_freqs*sizeof(float), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_flux_U, flux_U,
                             num_components*num_times*num_freqs*sizeof(float), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_flux_V, flux_V,
                             num_components*num_times*num_freqs*sizeof(float), cudaMemcpyHostToDevice ));

  float *d_sum_visi_XX_real = NULL;
  float *d_sum_visi_XY_real = NULL;
  float *d_sum_visi_YX_real = NULL;
  float *d_sum_visi_YY_real = NULL;
  float *d_sum_visi_XX_imag = NULL;
  float *d_sum_visi_XY_imag = NULL;
  float *d_sum_visi_YX_imag = NULL;
  float *d_sum_visi_YY_imag = NULL;

  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_XX_real, num_visis*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_XY_real, num_visis*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_YX_real, num_visis*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_YY_real, num_visis*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_XX_imag, num_visis*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_XY_imag, num_visis*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_YX_imag, num_visis*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_YY_imag, num_visis*sizeof(float) ));

  dim3 grid, threads;

  threads.x = 128;
  grid.x = (int)ceil( (float)num_visis / (float)threads.x );

  cudaErrorCheckKernel("kern_update_sum_visis",
                      kern_update_sum_visis, grid, threads,
                      num_freqs, num_baselines, num_components, num_times, beamtype,
                      (cuFloatComplex *)d_primay_beam_J00, (cuFloatComplex *)d_primay_beam_J01,
                      (cuFloatComplex *)d_primay_beam_J10, (cuFloatComplex *)d_primay_beam_J11,
                      (cuFloatComplex *)d_visi_components,
                      d_flux_I, d_flux_Q, d_flux_U, d_flux_V,
                      d_sum_visi_XX_real, d_sum_visi_XX_imag,
                      d_sum_visi_XY_real, d_sum_visi_XY_imag,
                      d_sum_visi_YX_real, d_sum_visi_YX_imag,
                      d_sum_visi_YY_real, d_sum_visi_YY_imag );

  cudaErrorCheckCall( cudaMemcpy(sum_visi_XX_real, d_sum_visi_XX_real,
                             num_visis*sizeof(float), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(sum_visi_XY_real, d_sum_visi_XY_real,
                             num_visis*sizeof(float), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(sum_visi_YX_real, d_sum_visi_YX_real,
                             num_visis*sizeof(float), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(sum_visi_YY_real, d_sum_visi_YY_real,
                             num_visis*sizeof(float), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(sum_visi_XX_imag, d_sum_visi_XX_imag,
                             num_visis*sizeof(float), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(sum_visi_XY_imag, d_sum_visi_XY_imag,
                             num_visis*sizeof(float), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(sum_visi_YX_imag, d_sum_visi_YX_imag,
                             num_visis*sizeof(float), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(sum_visi_YY_imag, d_sum_visi_YY_imag,
                             num_visis*sizeof(float), cudaMemcpyDeviceToHost ));

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
           float _Complex *primay_beam_J00, float _Complex *primay_beam_J01,
           float _Complex *primay_beam_J10, float _Complex *primay_beam_J11,
           float *freqs, float *ls, float *ms, float *ns,
           float *ras, float *decs, float *azs, float *zas,
           float *sin_para_angs, float *cos_para_angs,
           float *beam_has, float *beam_decs,
           woden_settings_t *woden_settings,
           beam_settings_t beam_settings,
           RTS_MWA_FEE_beam_t *FEE_beam){

  float _Complex *d_primay_beam_J00 = NULL;
  float _Complex *d_primay_beam_J01 = NULL;
  float _Complex *d_primay_beam_J10 = NULL;
  float _Complex *d_primay_beam_J11 = NULL;

  int num_beam_values = num_components*woden_settings->num_freqs*woden_settings->num_time_steps;

  cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J00,
                                      num_beam_values*sizeof(float _Complex) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J01,
                                      num_beam_values*sizeof(float _Complex) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J10,
                                      num_beam_values*sizeof(float _Complex) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J11,
                                      num_beam_values*sizeof(float _Complex) ));

  float *d_ls = NULL;
  float *d_ms = NULL;
  float *d_ns = NULL;

  cudaErrorCheckCall( cudaMalloc( (void**)&d_ls, num_components*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_ms, num_components*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_ns, num_components*sizeof(float) ));

  float *d_ras = NULL;
  float *d_decs = NULL;

  cudaErrorCheckCall( cudaMalloc( (void**)&d_ras, num_components*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_decs, num_components*sizeof(float) ));

  cudaErrorCheckCall( cudaMemcpy(d_ras, ras, num_components*sizeof(float), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_decs, decs, num_components*sizeof(float), cudaMemcpyHostToDevice ));

  float *d_freqs = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_freqs, woden_settings->num_freqs*sizeof(float) ) );
  cudaErrorCheckCall( cudaMemcpy( d_freqs, freqs, woden_settings->num_freqs*sizeof(float), cudaMemcpyHostToDevice ) );

  source_component_common(num_components,
           (cuFloatComplex *)d_primay_beam_J00,
           (cuFloatComplex *)d_primay_beam_J01,
           (cuFloatComplex *)d_primay_beam_J10,
           (cuFloatComplex *)d_primay_beam_J11,
           d_freqs, d_ls, d_ms, d_ns,
           d_ras, d_decs, azs, zas,
           sin_para_angs, cos_para_angs,
           beam_has, beam_decs,
           woden_settings, beam_settings,
           FEE_beam);

  cudaErrorCheckCall( cudaMemcpy(primay_beam_J00, d_primay_beam_J00,
              num_beam_values*sizeof(float _Complex), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(primay_beam_J01, d_primay_beam_J01,
              num_beam_values*sizeof(float _Complex), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(primay_beam_J10, d_primay_beam_J10,
              num_beam_values*sizeof(float _Complex), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(primay_beam_J11, d_primay_beam_J11,
              num_beam_values*sizeof(float _Complex), cudaMemcpyDeviceToHost ));

  cudaErrorCheckCall( cudaMemcpy(ls, d_ls, num_components*sizeof(float), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(ms, d_ms, num_components*sizeof(float), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(ns, d_ns, num_components*sizeof(float), cudaMemcpyDeviceToHost ));

  cudaErrorCheckCall( cudaFree( d_primay_beam_J00 ) );
  cudaErrorCheckCall( cudaFree( d_primay_beam_J01 ) );
  cudaErrorCheckCall( cudaFree( d_primay_beam_J10 ) );
  cudaErrorCheckCall( cudaFree( d_primay_beam_J11 ) );
  cudaErrorCheckCall( cudaFree( d_ls ) );
  cudaErrorCheckCall( cudaFree( d_ms ) );
  cudaErrorCheckCall( cudaFree( d_ns ) );
}


extern "C" void test_kern_calc_visi_point(int num_components, int num_baselines,
          int num_freqs, int num_visis, int num_times, int beamtype,
          float *ras, float *decs,
          float *component_freqs,
          float *flux_I, float *flux_Q, float *flux_U, float *flux_V,
          float *SIs, float *us, float *vs, float *ws,
          float _Complex *sum_visi_XX_real, float _Complex *sum_visi_XX_imag,
          float _Complex *sum_visi_XY_real, float _Complex *sum_visi_XY_imag,
          float _Complex *sum_visi_YX_real, float _Complex *sum_visi_YX_imag,
          float _Complex *sum_visi_YY_real, float _Complex *sum_visi_YY_imag,
          float *allsteps_wavelengths,
          float *ls, float *ms, float *ns,
          float _Complex *primay_beam_J00, float _Complex *primay_beam_J01,
          float _Complex *primay_beam_J10, float _Complex *primay_beam_J11){


  float *d_us = NULL;
  float *d_vs = NULL;
  float *d_ws = NULL;
  float *d_allsteps_wavelengths = NULL;

  cudaErrorCheckCall( cudaMalloc( (void**)&d_us, num_visis*sizeof(float) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_vs, num_visis*sizeof(float) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_ws, num_visis*sizeof(float) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_allsteps_wavelengths, num_visis*sizeof(float) ) );

  cudaErrorCheckCall( cudaMemcpy(d_us, us,
                             num_visis*sizeof(float), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_vs, vs,
                             num_visis*sizeof(float), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_ws, ws,
                             num_visis*sizeof(float), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_allsteps_wavelengths, allsteps_wavelengths,
                             num_visis*sizeof(float), cudaMemcpyHostToDevice ))

  float _Complex *d_primay_beam_J00 = NULL;
  float _Complex *d_primay_beam_J01 = NULL;
  float _Complex *d_primay_beam_J10 = NULL;
  float _Complex *d_primay_beam_J11 = NULL;

  int num_beam_values = num_components*num_freqs*num_times;

  cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J00,
                                      num_beam_values*sizeof(float _Complex) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J01,
                                      num_beam_values*sizeof(float _Complex) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J10,
                                      num_beam_values*sizeof(float _Complex) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J11,
                                      num_beam_values*sizeof(float _Complex) ));

  float *d_ls = NULL;
  float *d_ms = NULL;
  float *d_ns = NULL;
  float *d_component_freqs = NULL;

  cudaErrorCheckCall( cudaMalloc( (void**)&d_ls, num_components*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_ms, num_components*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_ns, num_components*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_component_freqs, num_components*sizeof(float) ));

  cudaErrorCheckCall( cudaMemcpy(d_ls, ls, num_components*sizeof(float), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_ms, ms, num_components*sizeof(float), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_ns, ns, num_components*sizeof(float), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_component_freqs, component_freqs,
                        num_components*sizeof(float), cudaMemcpyHostToDevice ));

  float *d_ras = NULL;
  float *d_decs = NULL;

  cudaErrorCheckCall( cudaMalloc( (void**)&d_ras, num_components*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_decs, num_components*sizeof(float) ));

  cudaErrorCheckCall( cudaMemcpy(d_ras, ras, num_components*sizeof(float), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_decs, decs, num_components*sizeof(float), cudaMemcpyHostToDevice ));

  float *d_flux_I = NULL;
  float *d_flux_Q = NULL;
  float *d_flux_U = NULL;
  float *d_flux_V = NULL;
  float *d_SIs = NULL;


  cudaErrorCheckCall( cudaMalloc( (void**)&d_flux_I, num_components*num_freqs*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_flux_Q, num_components*num_freqs*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_flux_U, num_components*num_freqs*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_flux_V, num_components*num_freqs*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_SIs, num_components*num_freqs*sizeof(float) ));

  cudaErrorCheckCall( cudaMemcpy(d_flux_I, flux_I,
                             num_components*num_freqs*sizeof(float), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_flux_Q, flux_Q,
                             num_components*num_freqs*sizeof(float), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_flux_U, flux_U,
                             num_components*num_freqs*sizeof(float), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_flux_V, flux_V,
                             num_components*num_freqs*sizeof(float), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_SIs, SIs,
                             num_components*num_freqs*sizeof(float), cudaMemcpyHostToDevice ));

  float *d_sum_visi_XX_real = NULL;
  float *d_sum_visi_XY_real = NULL;
  float *d_sum_visi_YX_real = NULL;
  float *d_sum_visi_YY_real = NULL;
  float *d_sum_visi_XX_imag = NULL;
  float *d_sum_visi_XY_imag = NULL;
  float *d_sum_visi_YX_imag = NULL;
  float *d_sum_visi_YY_imag = NULL;

  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_XX_real, num_visis*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_XY_real, num_visis*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_YX_real, num_visis*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_YY_real, num_visis*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_XX_imag, num_visis*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_XY_imag, num_visis*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_YX_imag, num_visis*sizeof(float) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_YY_imag, num_visis*sizeof(float) ));

  dim3 grid, threads;

  threads.x = 128;
  grid.x = (int)ceil( (float)num_visis / (float)threads.x );

  cudaErrorCheckKernel("kern_calc_visi_point",
      kern_calc_visi_point, grid, threads,
      d_ras, d_decs,
      d_component_freqs, d_flux_I, d_flux_Q, d_flux_U, d_flux_V, d_SIs,
      d_us, d_vs, d_ws,
      d_sum_visi_XX_real, d_sum_visi_XX_imag,
      d_sum_visi_XY_real, d_sum_visi_XY_imag,
      d_sum_visi_YX_real, d_sum_visi_YX_imag,
      d_sum_visi_YY_real, d_sum_visi_YY_imag,
      d_allsteps_wavelengths,
      d_ls, d_ms, d_ns,
      num_components, num_baselines, num_freqs, num_visis,
      num_times, beamtype,
      (cuFloatComplex *)d_primay_beam_J00, (cuFloatComplex *)d_primay_beam_J01,
      (cuFloatComplex *)d_primay_beam_J10, (cuFloatComplex *)d_primay_beam_J11);

  cudaErrorCheckCall( cudaMemcpy(sum_visi_XX_real, d_sum_visi_XX_real,
                             num_visis*sizeof(float), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(sum_visi_XY_real, d_sum_visi_XY_real,
                             num_visis*sizeof(float), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(sum_visi_YX_real, d_sum_visi_YX_real,
                             num_visis*sizeof(float), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(sum_visi_YY_real, d_sum_visi_YY_real,
                             num_visis*sizeof(float), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(sum_visi_XX_imag, d_sum_visi_XX_imag,
                             num_visis*sizeof(float), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(sum_visi_XY_imag, d_sum_visi_XY_imag,
                             num_visis*sizeof(float), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(sum_visi_YX_imag, d_sum_visi_YX_imag,
                             num_visis*sizeof(float), cudaMemcpyDeviceToHost ));
  cudaErrorCheckCall( cudaMemcpy(sum_visi_YY_imag, d_sum_visi_YY_imag,
                             num_visis*sizeof(float), cudaMemcpyDeviceToHost ));


  cudaErrorCheckCall(  cudaFree( d_ras ) );
  cudaErrorCheckCall(  cudaFree( d_decs ) );
  cudaErrorCheckCall(  cudaFree( d_component_freqs ) );
  cudaErrorCheckCall(  cudaFree( d_flux_I ) );
  cudaErrorCheckCall(  cudaFree( d_flux_Q ) );
  cudaErrorCheckCall(  cudaFree( d_flux_U ) );
  cudaErrorCheckCall(  cudaFree( d_flux_V ) );
  cudaErrorCheckCall(  cudaFree( d_SIs ) );
  cudaErrorCheckCall(  cudaFree( d_us ) );
  cudaErrorCheckCall(  cudaFree( d_vs ) );
  cudaErrorCheckCall(  cudaFree( d_ws ) );
  cudaErrorCheckCall(  cudaFree( d_sum_visi_XX_real ) );
  cudaErrorCheckCall(  cudaFree( d_sum_visi_XX_imag ) );
  cudaErrorCheckCall(  cudaFree( d_sum_visi_XY_real ) );
  cudaErrorCheckCall(  cudaFree( d_sum_visi_XY_imag ) );
  cudaErrorCheckCall(  cudaFree( d_sum_visi_YX_real ) );
  cudaErrorCheckCall(  cudaFree( d_sum_visi_YX_imag ) );
  cudaErrorCheckCall(  cudaFree( d_sum_visi_YY_real ) );
  cudaErrorCheckCall(  cudaFree( d_sum_visi_YY_imag ) );
  cudaErrorCheckCall(  cudaFree( d_allsteps_wavelengths ) );

  cudaErrorCheckCall( cudaFree( d_primay_beam_J00 ) );
  cudaErrorCheckCall( cudaFree( d_primay_beam_J01 ) );
  cudaErrorCheckCall( cudaFree( d_primay_beam_J10 ) );
  cudaErrorCheckCall( cudaFree( d_primay_beam_J11 ) );
  cudaErrorCheckCall( cudaFree( d_ls ) );
  cudaErrorCheckCall( cudaFree( d_ms ) );
  cudaErrorCheckCall( cudaFree( d_ns ) );
}
