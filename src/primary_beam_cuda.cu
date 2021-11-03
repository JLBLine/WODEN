#include <stdio.h>
#include <stdlib.h>
#include <cuComplex.h>
#include <complex.h>
#include <math.h>
#include "cudacomplex.h"
#include "constants.h"
#include "fundamental_coords.h"
#include "cudacheck.h"

__host__ __device__ void twoD_Gaussian(user_precision_t x, user_precision_t y,
           user_precision_t xo, user_precision_t yo,
           user_precision_t sigma_x, user_precision_t sigma_y,
           user_precision_t cos_theta, user_precision_t sin_theta,
           user_precision_t sin_2theta,
           user_precision_t * d_beam_real, user_precision_t * d_beam_imag) {

  user_precision_t a, b, c;

  a = (cos_theta*cos_theta)/(2*sigma_x*sigma_x) + (sin_theta*sin_theta)/(2*sigma_y*sigma_y);
  b = -sin_2theta / (4*sigma_x*sigma_x) + sin_2theta / (4*sigma_y*sigma_y);
  c = (sin_theta*sin_theta)/(2*sigma_x*sigma_x) + (cos_theta*cos_theta)/(2*sigma_y*sigma_y);

  * d_beam_real = exp( -( a*(x-xo)*(x-xo) + 2*b*(x-xo)*(y-yo) + c*(y-yo)*(y-yo) ));
  * d_beam_imag = 0.0;

}

__global__ void kern_gaussian_beam(double *d_beam_ls, double *d_beam_ms,
           user_precision_t beam_ref_freq, user_precision_t *d_freqs,
           user_precision_t fwhm_lm, user_precision_t cos_theta,
           user_precision_t sin_theta, user_precision_t sin_2theta,
           int num_freqs, int num_times, int num_components,
           cuUserComplex *d_primay_beam_J00, cuUserComplex *d_primay_beam_J11) {
  // Start by computing which baseline we're going to do
  const int iLMcoord = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iFreq = threadIdx.y + (blockDim.y*blockIdx.y);

  // TODO would it be faster to have freq and time on one axis, and num
  // components on other axis?
  if (iFreq < num_freqs && iLMcoord < num_components * num_times) {

    int component = (int)floorf((float)iLMcoord / (float)num_times);
    int time_ind = iLMcoord - component*num_times;
    int beam_ind = num_freqs*time_ind*num_components + (num_components*iFreq) + component;

    user_precision_t d_beam_real, d_beam_imag;
    //Convert FWHM into standard dev, and scale for frequency
    user_precision_t std = (fwhm_lm / FWHM_FACTOR) * (beam_ref_freq / d_freqs[iFreq]);

    twoD_Gaussian((user_precision_t)d_beam_ls[iLMcoord], (user_precision_t)d_beam_ms[iLMcoord], 0, 0,
               std, std, cos_theta, sin_theta, sin_2theta,
               &d_beam_real, &d_beam_imag);

    d_primay_beam_J00[beam_ind] = make_cuUserComplex(d_beam_real, d_beam_imag);
    d_primay_beam_J11[beam_ind] = make_cuUserComplex(d_beam_real, d_beam_imag);

  }
}


extern "C" void calculate_gaussian_beam(int num_components, int num_time_steps,
           int num_freqs, user_precision_t ha0,
           user_precision_t sdec0, user_precision_t cdec0,
           user_precision_t fwhm_lm, user_precision_t cos_theta,
           user_precision_t sin_theta, user_precision_t sin_2theta,
           user_precision_t beam_ref_freq, user_precision_t *d_freqs,
           double *beam_has, double *beam_decs,
           cuUserComplex *d_primay_beam_J00, cuUserComplex *d_primay_beam_J11){

  int num_beam_hadec = num_components * num_time_steps;

  double *d_beam_has = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_beam_has, num_beam_hadec*sizeof(double)) );
  cudaErrorCheckCall( cudaMemcpy( d_beam_has, beam_has,
                      num_beam_hadec*sizeof(double), cudaMemcpyHostToDevice) );

  double *d_beam_decs = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_beam_decs, num_beam_hadec*sizeof(double)) );
  cudaErrorCheckCall( cudaMemcpy( d_beam_decs, beam_decs,
                      num_beam_hadec*sizeof(double), cudaMemcpyHostToDevice) );

  double *d_beam_ls = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_beam_ls, num_beam_hadec*sizeof(double)) );

  double *d_beam_ms = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_beam_ms, num_beam_hadec*sizeof(double)) );

  double *d_beam_ns = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_beam_ns, num_beam_hadec*sizeof(double)) );

  dim3 grid, threads;
  threads.x = 128;
  grid.x = (int)ceil( (user_precision_t)num_beam_hadec / (user_precision_t)threads.x );

  cudaErrorCheckKernel("kern_calc_lmn", kern_calc_lmn, grid, threads,
                        ha0, sdec0, cdec0,
                        d_beam_has, d_beam_decs,
                        d_beam_ls, d_beam_ms, d_beam_ns, num_beam_hadec);

  threads.y = 2;
  grid.x = (int)ceil( (float)num_time_steps*(float)num_components / (float)threads.x );
  grid.y = (int)ceil( (float)num_freqs / (float)threads.y );
  //
  cudaErrorCheckKernel("kern_gaussian_beam",kern_gaussian_beam, grid, threads,
                       d_beam_ls, d_beam_ms, beam_ref_freq, d_freqs,
                       fwhm_lm, cos_theta, sin_theta, sin_2theta,
                       num_freqs, num_time_steps, num_components,
                       d_primay_beam_J00, d_primay_beam_J11);

  cudaFree( d_beam_ns );
  cudaFree( d_beam_ms );
  cudaFree( d_beam_ls );
  cudaFree( d_beam_decs );
  cudaFree( d_beam_has );
}

__device__ void analytic_dipole(user_precision_t az, user_precision_t za,
           user_precision_t wavelength,
           cuUserComplex * d_beam_X, cuUserComplex * d_beam_Y) {

  user_precision_t dipole_height_m = 0.3;

  //Here, X means a north-south aligned dipole
  //      Y means an east-west aligned dipole

  user_precision_t theta_parallel_X = acos(sin(za)*cos(az));
  user_precision_t theta_parallel_Y = acos(sin(za)*sin(az));

  user_precision_t d_in_lambda = (2. * dipole_height_m)/wavelength;
  user_precision_t gp_effect_array = 2. * sin(M_PI*d_in_lambda*cos(za));

  user_precision_t voltage_parallel_X = sin(theta_parallel_X) * gp_effect_array;
  user_precision_t voltage_parallel_Y = sin(theta_parallel_Y) * gp_effect_array;

  cuUserComplex tempX;
  cuUserComplex tempY;

  tempX.x = voltage_parallel_X;
  tempY.x = voltage_parallel_Y;

  //No imaginary components so set to 0
  tempX.y = 0;
  tempY.y = 0;

  * d_beam_X = tempX;
  * d_beam_Y = tempY;

}

__global__ void kern_analytic_dipole_beam(user_precision_t *d_azs,
           user_precision_t *d_zas,  user_precision_t *d_freqs, int num_freqs,
           int num_times, int num_components,
           cuUserComplex *d_primay_beam_J00, cuUserComplex *d_primay_beam_J11) {
  // Start by computing which baseline we're going to do
  const int iCoord = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iFreq = threadIdx.y + (blockDim.y*blockIdx.y);

  // TODO would it be faster to have freq and time on one axis, and num
  // components on other axis?
  if (iFreq < num_freqs && iCoord < num_components * num_times) {

    int component = (int)floorf((float)iCoord / (float)num_times);
    int time_ind = iCoord - component*num_times;

    int beam_ind = num_freqs*time_ind*num_components + (num_components*iFreq) + component;

    cuUserComplex d_beam_X, d_beam_Y;
    cuUserComplex d_beam_norm_X, d_beam_norm_Y;
    user_precision_t wavelength = VELC / d_freqs[iFreq];

    analytic_dipole(d_azs[iCoord], d_zas[iCoord], wavelength,
               &d_beam_X, &d_beam_Y);

    //Should really calculate this normalisation outside of this kernel - we are
    //massively increasing the number calculations for no reason
    analytic_dipole(0.0, 0.0, wavelength,
               &d_beam_norm_X, &d_beam_norm_Y);

    cuUserComplex normed_X = d_beam_X;
    cuUserComplex normed_Y = d_beam_Y;

    //Analytic beam is entirely real, so can just normalise by real values
    normed_X.x = normed_X.x / d_beam_norm_X.x;
    normed_Y.x = normed_Y.x / d_beam_norm_Y.x;

    d_primay_beam_J00[beam_ind] = normed_X;
    d_primay_beam_J11[beam_ind] = normed_Y;

  }
}

extern "C" void calculate_analytic_dipole_beam(int num_components,
     int num_time_steps, int num_freqs,
     user_precision_t *azs, user_precision_t *zas, user_precision_t *d_freqs,
     cuUserComplex *d_primay_beam_J00, cuUserComplex *d_primay_beam_J11){

  int num_beam_azza = num_components * num_time_steps;

  user_precision_t *d_azs = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_azs, num_beam_azza*sizeof(user_precision_t)) );
  cudaErrorCheckCall( cudaMemcpy(d_azs, azs, num_beam_azza*sizeof(user_precision_t),
                      cudaMemcpyHostToDevice) );

  user_precision_t *d_zas = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_zas, num_beam_azza*sizeof(user_precision_t)) );
  cudaErrorCheckCall( cudaMemcpy(d_zas, zas, num_beam_azza*sizeof(user_precision_t),
                      cudaMemcpyHostToDevice) );

  dim3 grid, threads;
  threads.x = 128;
  threads.y = 1;

  grid.x = (int)ceil( (float)num_beam_azza  / (float)threads.x );
  grid.y = (int)ceil( (float)num_freqs / (float)threads.y );

  cudaErrorCheckKernel("kern_analytic_dipole_beam",
             kern_analytic_dipole_beam, grid, threads,
             d_azs, d_zas, d_freqs, num_freqs,
             num_time_steps, num_components,
             d_primay_beam_J00, d_primay_beam_J11);

  cudaErrorCheckCall( cudaFree(d_azs) );
  cudaErrorCheckCall( cudaFree(d_zas) );

}

/*******************************************************************************
                 Functions below to be used in unit tests
*******************************************************************************/

extern "C" void test_analytic_dipole_beam(int num_components,
     int num_time_steps, int num_freqs,
     user_precision_t *azs, user_precision_t *zas, user_precision_t *freqs,
     user_precision_complex_t *analy_beam_X,
     user_precision_complex_t *analy_beam_Y) {

  user_precision_complex_t *d_analy_beam_X = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_analy_beam_X,
    num_freqs*num_time_steps*num_components*sizeof(user_precision_complex_t)) );

  user_precision_complex_t *d_analy_beam_Y = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_analy_beam_Y,
     num_freqs*num_time_steps*num_components*sizeof(user_precision_complex_t)) );

  user_precision_t *d_freqs = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_freqs, num_freqs*sizeof(user_precision_t)) );
  cudaErrorCheckCall( cudaMemcpy(d_freqs, freqs, num_freqs*sizeof(user_precision_t),
                      cudaMemcpyHostToDevice) );

  calculate_analytic_dipole_beam(num_components,
      num_time_steps, num_freqs,
      azs, zas, d_freqs,
      (cuUserComplex *)d_analy_beam_X, (cuUserComplex *)d_analy_beam_Y);

  cudaErrorCheckCall( cudaMemcpy(analy_beam_X, d_analy_beam_X,
             num_freqs*num_time_steps*num_components*sizeof(user_precision_complex_t),
             cudaMemcpyDeviceToHost) );
  cudaErrorCheckCall( cudaMemcpy(analy_beam_Y, d_analy_beam_Y,
             num_freqs*num_time_steps*num_components*sizeof(user_precision_complex_t),
             cudaMemcpyDeviceToHost) );

  cudaErrorCheckCall( cudaFree(d_analy_beam_X) );
  cudaErrorCheckCall( cudaFree(d_analy_beam_Y) );
  cudaErrorCheckCall( cudaFree(d_freqs) );

}

extern "C" void test_kern_gaussian_beam(double *beam_ls, double *beam_ms,
           user_precision_t beam_ref_freq, user_precision_t *freqs,
           user_precision_t fwhm_lm, user_precision_t cos_theta, user_precision_t sin_theta, user_precision_t sin_2theta,
           int num_freqs, int num_time_steps, int num_components,
           user_precision_complex_t *primay_beam_J00, user_precision_complex_t *primay_beam_J11) {

  int num_beam_hadec = num_components * num_time_steps;

  double *d_beam_ls = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_beam_ls, num_beam_hadec*sizeof(double)) );
  cudaErrorCheckCall( cudaMemcpy(d_beam_ls, beam_ls,
                           num_beam_hadec*sizeof(double), cudaMemcpyHostToDevice ) );

  double *d_beam_ms = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_beam_ms, num_beam_hadec*sizeof(double)) );
  cudaErrorCheckCall( cudaMemcpy(d_beam_ms, beam_ms,
                           num_beam_hadec*sizeof(double), cudaMemcpyHostToDevice ) );

  user_precision_t *d_freqs = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_freqs, num_freqs*sizeof(user_precision_t) ) );
  cudaErrorCheckCall( cudaMemcpy(d_freqs, freqs,
                           num_freqs*sizeof(user_precision_t), cudaMemcpyHostToDevice ) );

  user_precision_complex_t *d_primay_beam_J00 = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J00,
                                   num_freqs*num_beam_hadec*sizeof(user_precision_complex_t)) );

  user_precision_complex_t *d_primay_beam_J11 = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J11,
                                   num_freqs*num_beam_hadec*sizeof(user_precision_complex_t)) );

  dim3 grid, threads;

  threads.x = 16;
  grid.x = (int)ceil( (float)num_beam_hadec / (float)threads.x );

  threads.y = 16;
  grid.y = (int)ceil( (float)num_freqs / (float)threads.y );

  cudaErrorCheckKernel("kern_gaussian_beam",
                        kern_gaussian_beam, grid, threads,
                        d_beam_ls, d_beam_ms,
                        beam_ref_freq, d_freqs,
                        fwhm_lm, cos_theta, sin_theta, sin_2theta,
                        num_freqs, num_time_steps, num_components,
                        (cuUserComplex *)d_primay_beam_J00,
                        (cuUserComplex *)d_primay_beam_J11);

  cudaErrorCheckCall( cudaMemcpy(primay_beam_J00, d_primay_beam_J00,
                 num_freqs*num_beam_hadec*sizeof(user_precision_complex_t), cudaMemcpyDeviceToHost) );
  cudaErrorCheckCall( cudaMemcpy(primay_beam_J11, d_primay_beam_J11,
                 num_freqs*num_beam_hadec*sizeof(user_precision_complex_t), cudaMemcpyDeviceToHost) );

  cudaErrorCheckCall( cudaFree(d_beam_ls ) );
  cudaErrorCheckCall( cudaFree(d_beam_ms ) );
  cudaErrorCheckCall( cudaFree(d_freqs ) );
  cudaErrorCheckCall( cudaFree(d_primay_beam_J00 ) );
  cudaErrorCheckCall( cudaFree(d_primay_beam_J11 ) );

}


extern "C" void test_calculate_gaussian_beam(int num_components, int num_time_steps,
     int num_freqs, user_precision_t ha0, user_precision_t sdec0, user_precision_t cdec0,
     user_precision_t fwhm_lm, user_precision_t cos_theta, user_precision_t sin_theta, user_precision_t sin_2theta,
     user_precision_t beam_ref_freq, user_precision_t *freqs,
     double *beam_has, double *beam_decs,
     user_precision_complex_t *primay_beam_J00, user_precision_complex_t *primay_beam_J11) {

  int num_beam_hadec = num_components * num_time_steps;

  user_precision_t *d_freqs = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_freqs, num_freqs*sizeof(user_precision_t) ) );
  cudaErrorCheckCall( cudaMemcpy(d_freqs, freqs,
                           num_freqs*sizeof(user_precision_t), cudaMemcpyHostToDevice ) );

  user_precision_complex_t *d_primay_beam_J00 = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J00,
                                   num_freqs*num_beam_hadec*sizeof(user_precision_complex_t)) );

  user_precision_complex_t *d_primay_beam_J11 = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J11,
                                   num_freqs*num_beam_hadec*sizeof(user_precision_complex_t)) );


  calculate_gaussian_beam(num_components, num_time_steps,
                         num_freqs, ha0, sdec0, cdec0,
                         fwhm_lm, cos_theta, sin_theta, sin_2theta,
                         beam_ref_freq, d_freqs,
                         beam_has, beam_decs,
                         (cuUserComplex *)d_primay_beam_J00,
                         (cuUserComplex *)d_primay_beam_J11);

  cudaErrorCheckCall( cudaMemcpy(primay_beam_J00, d_primay_beam_J00,
                 num_freqs*num_beam_hadec*sizeof(user_precision_complex_t), cudaMemcpyDeviceToHost) );
  cudaErrorCheckCall( cudaMemcpy(primay_beam_J11, d_primay_beam_J11,
                 num_freqs*num_beam_hadec*sizeof(user_precision_complex_t), cudaMemcpyDeviceToHost) );

  cudaErrorCheckCall( cudaFree(d_freqs ) );
  cudaErrorCheckCall( cudaFree(d_primay_beam_J00 ) );
  cudaErrorCheckCall( cudaFree(d_primay_beam_J11 ) );

}
