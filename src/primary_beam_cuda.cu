#include <stdio.h>
#include <stdlib.h>
#include <cuComplex.h>
#include <complex.h>
#include <math.h>
#include "cudacomplex.h"
#include "constants.h"
#include "fundamental_coords.h"
#include "cudacheck.h"


__device__ void twoD_Gaussian(float x, float y, float xo, float yo,
           float sigma_x, float sigma_y,
           float cos_theta, float sin_theta, float sin_2theta,
           float * d_beam_real, float * d_beam_imag) {

  float a, b, c;

  a = (cos_theta*cos_theta)/(2*sigma_x*sigma_x) + (sin_theta*sin_theta)/(2*sigma_y*sigma_y);
  b = -sin_2theta / (4*sigma_x*sigma_x) + sin_2theta / (4*sigma_y*sigma_y);
  c = (sin_theta*sin_theta)/(2*sigma_x*sigma_x) + (cos_theta*cos_theta)/(2*sigma_y*sigma_y);

  * d_beam_real = expf( -( a*(x-xo)*(x-xo) + 2*b*(x-xo)*(y-yo) + c*(y-yo)*(y-yo) ));
  * d_beam_imag = 0.0;

}

__global__ void kern_gaussian_beam(float *d_beam_ls, float *d_beam_ms,
           float beam_ref_freq, float *d_freqs,
           float fwhm_lm, float cos_theta, float sin_theta, float sin_2theta,
           int num_freqs, int num_times, int num_components,
           cuFloatComplex *d_primay_beam_J00, cuFloatComplex *d_primay_beam_J11) {
  // Start by computing which baseline we're going to do
  const int iLMcoord = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iFreq = threadIdx.y + (blockDim.y*blockIdx.y);

  // TODO would it be faster to have freq and time on one axis, and num
  // components on other axis?
  if (iFreq < num_freqs && iLMcoord < num_components * num_times) {

    //The l,m coords for the beam locations are for every component, every time
    //step in a 1D array. We are calculating over frequency as well and then
    //putting results into d_gauss_beam_reals, d_gauss_beam_imags, so have to do some
    //fun index maths to work out what time / component / freq we are on,
    //and where to put the output
    int time_ind = (int)floorf((float)iLMcoord / (float)num_components);
    int component = iLMcoord - time_ind*num_components;
    int beam_ind = num_freqs*time_ind*num_components + (num_components*iFreq) + component;

    float d_beam_real, d_beam_imag;
    float std = (fwhm_lm / FWHM_FACTOR) * (beam_ref_freq / d_freqs[iFreq]);

    twoD_Gaussian(d_beam_ls[iLMcoord], d_beam_ms[iLMcoord], 0, 0,
               std, std, cos_theta, sin_theta, sin_2theta,
               &d_beam_real, &d_beam_imag);

    d_primay_beam_J00[beam_ind] = make_cuComplex(d_beam_real, d_beam_imag);
    d_primay_beam_J11[beam_ind] = make_cuComplex(d_beam_real, d_beam_imag);

    // printf("%f %f %f %f %f\n",std, fwhm_lm, FWHM_FACTOR, beam_ref_freq , d_freqs[iFreq] );

    // printf("iFreq,iLMcoord,d_beam_ls[iLMcoord],d_beam_ms[iLMcoord],d_beam_real %d %d %f %f %f\n",iFreq,iLMcoord,d_beam_ls[iLMcoord],d_beam_ms[iLMcoord],d_beam_real );

  }
}


extern "C" void calculate_gaussian_beam(int num_components, int num_time_steps, int num_freqs,
     float fwhm_lm, float cos_theta, float sin_theta, float sin_2theta,
     float beam_ref_freq, float *d_freqs, float *d_beam_angles_array,
     float *beam_has, float *beam_decs,
     cuFloatComplex *d_primay_beam_J00, cuFloatComplex *d_primay_beam_J11){

  int num_beam_hadec = num_components * num_time_steps;

  float *d_beam_has = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_beam_has, num_beam_hadec*sizeof(float)) );
  cudaErrorCheckCall( cudaMemcpy( d_beam_has, beam_has,
                      num_beam_hadec*sizeof(float), cudaMemcpyHostToDevice) );

  float *d_beam_decs = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_beam_decs, num_beam_hadec*sizeof(float)) );
  cudaErrorCheckCall( cudaMemcpy( d_beam_decs, beam_decs,
                      num_beam_hadec*sizeof(float), cudaMemcpyHostToDevice) );

  float *d_beam_ls = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_beam_ls, num_beam_hadec*sizeof(float)) );

  float *d_beam_ms = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_beam_ms, num_beam_hadec*sizeof(float)) );

  float *d_beam_ns = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_beam_ns, num_beam_hadec*sizeof(float)) );

  dim3 grid, threads;
  threads.x = 128;
  grid.x = (int)ceil( (float)num_beam_hadec / (float)threads.x );

  cudaErrorCheckKernel("kern_calc_lmn", kern_calc_lmn, grid, threads,
                        d_beam_angles_array, d_beam_has, d_beam_decs,
                        d_beam_ls, d_beam_ms, d_beam_ns, num_beam_hadec);

  threads.y = 2;
  grid.x = (int)ceil( (float)num_time_steps*float(num_components) / (float)threads.x );
  grid.y = (int)ceil( (float)num_freqs / (float)threads.y );
  //
  cudaErrorCheckKernel("kern_gaussian_beam",kern_gaussian_beam, grid, threads,
                       d_beam_ls, d_beam_ms, beam_ref_freq, d_freqs,
                       fwhm_lm, cos_theta, sin_theta, sin_2theta,
                       num_freqs, num_time_steps, num_components,
                       d_primay_beam_J00, d_primay_beam_J11);

  // printf("Finished a gaussian beam kernel\n");

  cudaFree( d_beam_ns );
  cudaFree( d_beam_ms );
  cudaFree( d_beam_ls );
  cudaFree( d_beam_decs );
  cudaFree( d_beam_has );
}

// extern "C" void testing_gaussian_beam( float *beam_has, float *beam_decs,
//            float *beam_angles_array, float *beam_freqs, float *ref_freq_array,
//            float *beam_ls, float *beam_ms,
//            int num_components, int num_times, int num_freqs,
//            float *beam_reals, float *beam_imags){
//
//
//   int num_beam_hadec = num_components * num_times;
//   int num_beam_calcs = num_components * num_times * num_freqs;
//
//   float *d_beam_angles_array = NULL;
//   cudaMalloc( (void**)&d_beam_angles_array, 3*sizeof(float) );
//   cudaMemcpy( d_beam_angles_array, beam_angles_array, 3*sizeof(float), cudaMemcpyHostToDevice );
//
//   float *d_beam_has = NULL;
//   cudaMalloc( (void**)&d_beam_has, num_beam_hadec*sizeof(float) );
//   cudaMemcpy( d_beam_has, beam_has, num_beam_hadec*sizeof(float), cudaMemcpyHostToDevice );
//
//   float *d_beam_decs = NULL;
//   cudaMalloc( (void**)&d_beam_decs, num_beam_hadec*sizeof(float) );
//   cudaMemcpy( d_beam_decs, beam_decs, num_beam_hadec*sizeof(float), cudaMemcpyHostToDevice );
//
//   float *d_beam_ls = NULL;
//   cudaMalloc( (void**)&d_beam_ls, num_beam_hadec*sizeof(float) );
//
//   float *d_beam_ms = NULL;
//   cudaMalloc( (void**)&d_beam_ms, num_beam_hadec*sizeof(float) );
//
//   float *d_beam_ns = NULL;
//   cudaMalloc( (void**)&d_beam_ns, num_beam_hadec*sizeof(float) );
//
//   float *d_freqs = NULL;
//   cudaMalloc( (void**)&d_freqs, num_freqs*sizeof(float) );
//   cudaMemcpy( d_freqs, beam_freqs, num_freqs*sizeof(float), cudaMemcpyHostToDevice );
//
//   //TODO I cannot for the life of me work out how to cudaMalloc and Memcpy
//   //a single float (argh) so put the ref freq in an array (embarrassment)
//
//   float *d_beam_ref_freq = NULL;
//   cudaMalloc( (void**)&d_beam_ref_freq, sizeof(float) );
//   cudaMemcpy( d_beam_ref_freq, ref_freq_array, sizeof(float), cudaMemcpyHostToDevice );
//
//   float *d_gauss_beam_reals = NULL;
//   cudaMalloc( (void**)&d_gauss_beam_reals, num_beam_calcs*sizeof(float) );
//
//   float *d_gauss_beam_imags = NULL;
//   cudaMalloc( (void**)&d_gauss_beam_imags, num_beam_calcs*sizeof(float) );
//
//
//   dim3 grid, threads;
//
//   threads.x = 128;
//   grid.x = (int)ceil( (float)num_beam_hadec / (float)threads.x );
//
//   kern_calc_lmn<<< grid, threads >>>(d_beam_angles_array, d_beam_has, d_beam_decs,
//                 d_beam_ls, d_beam_ms, d_beam_ns, num_beam_hadec);
//
//   cudaMemcpy(beam_ls, d_beam_ls, num_beam_hadec*sizeof(float), cudaMemcpyDeviceToHost);
//   cudaMemcpy(beam_ms, d_beam_ms, num_beam_hadec*sizeof(float), cudaMemcpyDeviceToHost);
//
//   float beam_theta = 0;
//   float cos_theta = cosf(beam_theta);
//   float sin_theta = sinf(beam_theta);
//   float sin_2theta = sinf(2*beam_theta);
//
//   float fwhm_deg = 20.0;
//   float fwhm_lm = sinf(fwhm_deg * (M_PI/180.0));
//
//   threads.y = 2;
//   grid.x = (int)ceil( (float)num_times*float(num_components) / (float)threads.x );
//   grid.y = (int)ceil( (float)num_freqs / (float)threads.y );
//
//
//   kern_gaussian_beam<<< grid, threads >>>(d_beam_ls, d_beam_ms,
//              d_beam_ref_freq, d_freqs,
//              fwhm_lm, cos_theta, sin_theta, sin_2theta,
//              num_freqs, num_times, num_components,
//              d_primay_beam_J00, d_primay_beam_J11);
//
//   cudaMemcpy(beam_reals, d_gauss_beam_reals, num_beam_calcs*sizeof(float), cudaMemcpyDeviceToHost);
//   cudaMemcpy(beam_imags, d_gauss_beam_imags, num_beam_calcs*sizeof(float), cudaMemcpyDeviceToHost);
//
//   cudaFree( d_gauss_beam_imags);
//   cudaFree( d_gauss_beam_reals);
//   cudaFree( d_beam_ref_freq);
//   cudaFree( d_freqs);
//   cudaFree( d_beam_ns );
//   cudaFree( d_beam_ms );
//   cudaFree( d_beam_ls );
//   cudaFree( d_beam_decs );
//   cudaFree( d_beam_has );
//   cudaFree( d_beam_angles_array );
//
// }


__device__ void analytic_dipole(float az, float za, float wavelength,
           cuFloatComplex * d_beam_X, cuFloatComplex * d_beam_Y) {

  float dipole_height_m = 0.3;

  float theta_parallel_X = acos(sin(za)*sin(az));
  float theta_parallel_Y = acos(sin(za)*cos(az));

  float d_in_lambda = (2. * dipole_height_m)/wavelength;
  float gp_effect_array = 2. * sin(M_PI*d_in_lambda*cos(za));

  float voltage_parallel_X = sin(theta_parallel_X) * gp_effect_array;
  float voltage_parallel_Y = sin(theta_parallel_Y) * gp_effect_array;

  cuFloatComplex tempX;
  cuFloatComplex tempY;

  tempX.x = voltage_parallel_X;
  tempY.x = voltage_parallel_Y;

  // printf("%.5f %.5f %.5f %.5f %.5f\n",az,za,wavelength,tempX.x,tempY.x );

  //No imaginary components so set to 0
  tempX.y = 0;
  tempY.y = 0;

  * d_beam_X = tempX;
  * d_beam_Y = tempY;

}

__global__ void kern_analytic_dipole_beam(float *d_azs, float *d_zas,
           float *d_freqs, int num_freqs, int num_times, int num_components,
           cuFloatComplex *d_analy_beam_X, cuFloatComplex *d_analy_beam_Y) {
  // Start by computing which baseline we're going to do
  const int iCoord = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iFreq = threadIdx.y + (blockDim.y*blockIdx.y);

  // TODO would it be faster to have freq and time on one axis, and num
  // components on other axis?
  if (iFreq < num_freqs && iCoord < num_components * num_times) {

    int component = (int)floorf((float)iCoord / (float)num_times);
    int time_ind = iCoord - component*num_times;

    int beam_ind = num_freqs*time_ind*num_components + (num_components*iFreq) + component;

    cuFloatComplex d_beam_X, d_beam_Y;
    cuFloatComplex d_beam_norm_X, d_beam_norm_Y;
    float wavelength = VELC / d_freqs[iFreq];

    analytic_dipole(d_azs[iCoord], d_zas[iCoord], wavelength,
               &d_beam_X, &d_beam_Y);

    analytic_dipole(0.0, 0.0, wavelength,
               &d_beam_norm_X, &d_beam_norm_Y);

    cuFloatComplex normed_X = d_beam_X;
    cuFloatComplex normed_Y = d_beam_Y;
    normed_X.x = normed_X.x / d_beam_norm_X.x;
    normed_Y.x = normed_Y.x / d_beam_norm_Y.x;

    d_analy_beam_X[beam_ind] = normed_X;
    d_analy_beam_Y[beam_ind] = normed_Y;

  }
}


extern "C" void calculate_analytic_dipole_beam(int num_components,
     int num_time_steps, int num_freqs,
     float *azs, float *zas, float *d_freqs,
     cuFloatComplex *d_analy_beam_X, cuFloatComplex *d_analy_beam_Y){

  int num_beam_azza = num_components * num_time_steps;

  // for (size_t i = 0; i < num_components; i++) {
  //   printf("Analy beam az,za %.5f %.5f \n",azs[i],zas[i] );
  // }

  float *d_azs = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_azs, num_beam_azza*sizeof(float)) );
  cudaErrorCheckCall( cudaMemcpy(d_azs, azs, num_beam_azza*sizeof(float),
                      cudaMemcpyHostToDevice) );

  float *d_zas = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_zas, num_beam_azza*sizeof(float)) );
  cudaErrorCheckCall( cudaMemcpy(d_zas, zas, num_beam_azza*sizeof(float),
                      cudaMemcpyHostToDevice) );

  dim3 grid, threads;
  threads.x = 128;
  threads.y = 1;

  grid.x = (int)ceil( (float)num_beam_azza  / (float)threads.x );
  grid.y = (int)ceil( (float)num_freqs / (float)threads.y );

  // printf("Doing a beam gaussian beam kernel\n");
  // kern_analytic_dipole_beam<<< grid, threads >>>(d_azs, d_zas,
  //            d_freqs, num_freqs, num_time_steps, num_components,
  //            d_analy_beam_X, d_analy_beam_Y);

  cudaErrorCheckKernel("kern_analytic_dipole_beam",
             kern_analytic_dipole_beam, grid, threads,
             d_azs, d_zas, d_freqs, num_freqs,
             num_time_steps, num_components,
             d_analy_beam_X, d_analy_beam_Y);

  cudaErrorCheckCall( cudaFree(d_azs) );
  cudaErrorCheckCall( cudaFree(d_zas) );

}

extern "C" void test_analytic_dipole_beam(int num_components,
     int num_time_steps, int num_freqs,
     float *azs, float *zas, float *freqs,
     float _Complex *analy_beam_X, float _Complex *analy_beam_Y) {

  float _Complex *d_analy_beam_X = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_analy_beam_X,
              num_freqs*num_time_steps*num_components*sizeof(float _Complex)) );

  float _Complex *d_analy_beam_Y = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_analy_beam_Y,
               num_freqs*num_time_steps*num_components*sizeof(float _Complex)) );

  float *d_freqs = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_freqs, num_freqs*sizeof(float)) );
  cudaErrorCheckCall( cudaMemcpy(d_freqs, freqs, num_freqs*sizeof(float),
                      cudaMemcpyHostToDevice) );

  calculate_analytic_dipole_beam(num_components,
      num_time_steps, num_freqs,
      azs, zas, d_freqs,
      (cuFloatComplex *)d_analy_beam_X, (cuFloatComplex *)d_analy_beam_Y);

  cudaErrorCheckCall( cudaMemcpy(analy_beam_X, d_analy_beam_X,
             num_freqs*num_time_steps*num_components*sizeof(float _Complex),
             cudaMemcpyDeviceToHost) );
  cudaErrorCheckCall( cudaMemcpy(analy_beam_Y, d_analy_beam_Y,
             num_freqs*num_time_steps*num_components*sizeof(float _Complex),
             cudaMemcpyDeviceToHost) );

  cudaErrorCheckCall( cudaFree(analy_beam_X) );
  cudaErrorCheckCall( cudaFree(analy_beam_Y) );
  cudaErrorCheckCall( cudaFree(d_freqs) );

}
