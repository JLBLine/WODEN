#include <stdio.h>
#include <stdlib.h>
#include <cuComplex.h>
#include <complex.h>
#include <math.h>
#include "constants.h"
#include "fundamental_coords.h"


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
           float *ref_freq, float *d_freqs,
           float fwhm_lm, float cos_theta, float sin_theta, float sin_2theta,
           int num_freqs, int num_times, int num_components,
           float *d_beam_reals, float *d_beam_imags) {
  // Start by computing which baseline we're going to do
  const int iLMcoord = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iFreq = threadIdx.y + (blockDim.y*blockIdx.y);

  // TODO would it be faster to have freq and time on one axis, and num
  // components on other axis?
  if (iFreq < num_freqs && iLMcoord < num_components * num_times) {

    //The l,m coords for the beam locations are for every component, every time
    //step in a 1D array. We are calculating over frequency as well and then
    //putting results into d_beam_reals, d_beam_imags, so have to do some
    //fun index maths to work out what time / component / freq we are on,
    //and where to put the output
    int time_ind = (int)floorf((float)iLMcoord / (float)num_components);
    int component = iLMcoord - time_ind*num_components;
    int beam_ind = num_freqs*time_ind*num_components + (num_components*iFreq) + component;

    float d_beam_real, d_beam_imag;

    float std = (fwhm_lm / FWHM_FACTOR) * (ref_freq[0] / d_freqs[iFreq]);

    twoD_Gaussian(d_beam_ls[iLMcoord], d_beam_ms[iLMcoord], 0, 0,
               std, std, cos_theta, sin_theta, sin_2theta,
               &d_beam_real, &d_beam_imag);

    d_beam_reals[beam_ind] = d_beam_real;
    d_beam_imags[beam_ind] = d_beam_imag;

  }
}

extern "C" void testing_gaussian_beam( float *beam_has, float *beam_decs,
           float *beam_angles_array, float *beam_freqs, float *ref_freq_array,
           float *beam_ls, float *beam_ms,
           int num_components, int num_times, int num_freqs,
           float *beam_reals, float *beam_imags){


  int num_beam_hadec = num_components * num_times;
  int num_beam_calcs = num_components * num_times * num_freqs;

  float *d_beam_angles_array = NULL;
  cudaMalloc( (void**)&d_beam_angles_array, 3*sizeof(float) );
  cudaMemcpy( d_beam_angles_array, beam_angles_array, 3*sizeof(float), cudaMemcpyHostToDevice );

  float *d_beam_has = NULL;
  cudaMalloc( (void**)&d_beam_has, num_beam_hadec*sizeof(float) );
  cudaMemcpy( d_beam_has, beam_has, num_beam_hadec*sizeof(float), cudaMemcpyHostToDevice );

  float *d_beam_decs = NULL;
  cudaMalloc( (void**)&d_beam_decs, num_beam_hadec*sizeof(float) );
  cudaMemcpy( d_beam_decs, beam_decs, num_beam_hadec*sizeof(float), cudaMemcpyHostToDevice );

  float *d_beam_ls = NULL;
  cudaMalloc( (void**)&d_beam_ls, num_beam_hadec*sizeof(float) );

  float *d_beam_ms = NULL;
  cudaMalloc( (void**)&d_beam_ms, num_beam_hadec*sizeof(float) );

  float *d_beam_ns = NULL;
  cudaMalloc( (void**)&d_beam_ns, num_beam_hadec*sizeof(float) );

  float *d_freqs = NULL;
  cudaMalloc( (void**)&d_freqs, num_freqs*sizeof(float) );
  cudaMemcpy( d_freqs, beam_freqs, num_freqs*sizeof(float), cudaMemcpyHostToDevice );

  //TODO I cannot for the life of me work out how to cudaMalloc and Memcpy
  //a single float (argh) so put the ref freq in an array (embarrassment)

  float *d_beam_ref_freq = NULL;
  cudaMalloc( (void**)&d_beam_ref_freq, sizeof(float) );
  cudaMemcpy( d_beam_ref_freq, ref_freq_array, sizeof(float), cudaMemcpyHostToDevice );

  float *d_beam_reals = NULL;
  cudaMalloc( (void**)&d_beam_reals, num_beam_calcs*sizeof(float) );

  float *d_beam_imags = NULL;
  cudaMalloc( (void**)&d_beam_imags, num_beam_calcs*sizeof(float) );


  dim3 grid, threads;

  threads.x = 128;
  grid.x = (int)ceil( (float)num_beam_hadec / (float)threads.x );

  kern_calc_lmn<<< grid, threads >>>(d_beam_angles_array, d_beam_has, d_beam_decs,
                d_beam_ls, d_beam_ms, d_beam_ns, num_beam_hadec);

  cudaMemcpy(beam_ls, d_beam_ls, num_beam_hadec*sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(beam_ms, d_beam_ms, num_beam_hadec*sizeof(float), cudaMemcpyDeviceToHost);

  float beam_theta = 0;
  float cos_theta = cosf(beam_theta);
  float sin_theta = sinf(beam_theta);
  float sin_2theta = sinf(2*beam_theta);

  float fwhm_deg = 20.0;
  float fwhm_lm = sinf(fwhm_deg * (M_PI/180.0));

  threads.y = 2;
  grid.x = (int)ceil( (float)num_times*float(num_components) / (float)threads.x );
  grid.y = (int)ceil( (float)num_freqs / (float)threads.y );


  kern_gaussian_beam<<< grid, threads >>>(d_beam_ls, d_beam_ms,
             d_beam_ref_freq, d_freqs,
             fwhm_lm, cos_theta, sin_theta, sin_2theta,
             num_freqs, num_times, num_components,
             d_beam_reals, d_beam_imags);

  cudaMemcpy(beam_reals, d_beam_reals, num_beam_calcs*sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(beam_imags, d_beam_imags, num_beam_calcs*sizeof(float), cudaMemcpyDeviceToHost);

  cudaFree( d_beam_imags);
  cudaFree( d_beam_reals);
  cudaFree( d_beam_ref_freq);
  cudaFree( d_freqs);
  cudaFree( d_beam_ns );
  cudaFree( d_beam_ms );
  cudaFree( d_beam_ls );
  cudaFree( d_beam_decs );
  cudaFree( d_beam_has );
  cudaFree( d_beam_angles_array );

}
