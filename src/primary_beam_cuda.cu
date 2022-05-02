#include <stdio.h>
#include <stdlib.h>
#include <cuComplex.h>
#include <complex.h>
#include <math.h>
#include <erfa.h>
#include <mwa_hyperbeam.h>
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
           double beam_ref_freq, double *d_freqs,
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
    user_precision_t std = (fwhm_lm / FWHM_FACTOR) * (user_precision_t)(beam_ref_freq / d_freqs[iFreq]);

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
           double beam_ref_freq, double *d_freqs,
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
           user_precision_t *d_zas,  double *d_freqs, int num_freqs,
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
     user_precision_t *azs, user_precision_t *zas, double *d_freqs,
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

__device__ void RTS_MWA_beam(user_precision_t az, user_precision_t za,
           double ha, double dec,
           double wavelength, double *d_metre_delays,
           double latitude, int norm,
           cuUserComplex * gx, cuUserComplex * Dx,
           cuUserComplex * Dy, cuUserComplex * gy) {

  // set elements of the look-dir vector
  double proj_e = sin(za)*sin(az);
  double proj_n = sin(za)*cos(az);
  double proj_z = cos(za);


  int n_cols = 4;
  int n_rows = 4;

  //Used in calculating the phase later on
  double multiplier = -2 * M_PI / wavelength;
  double dipl_e, dipl_n, dipl_z;

  cuUserComplex x_dip;
  cuUserComplex y_dip;

  cuUserComplex gx_dip = {0.0, 0.0};
  cuUserComplex Dx_dip = {0.0, 0.0};
  cuUserComplex Dy_dip = {0.0, 0.0};
  cuUserComplex gy_dip = {0.0, 0.0};

  int k = 0;

  for (int i = 0; i < n_cols; i++) {
    for (int j = 0; j < n_rows; j++) {

      // set elements of the baseline vector
      dipl_e = (i - 1.5) * MWA_DIPOLE_SEP;
      dipl_n = (j - 1.5) * MWA_DIPOLE_SEP;
      dipl_z = 0.0;

      double phase = multiplier*(dipl_e*proj_e + dipl_n*proj_n + dipl_z*proj_z - d_metre_delays[k]);
      double phaseshift_re = cos(phase);
      double phaseshift_im = sin(phase);

      //TODO You could get gains here from some dipole flagging scheme in the future
      user_precision_t gainx = 1.0;
      user_precision_t gainy = 1.0;

      x_dip.x = gainx*phaseshift_re;
      x_dip.y = gainx*phaseshift_im;
      y_dip.x = gainy*phaseshift_re;
      y_dip.y = gainy*phaseshift_im;

      gx_dip += x_dip;
      Dx_dip += x_dip;
      Dy_dip += y_dip;
      gy_dip += y_dip;

      k += 1;
    }
  }

  //Calculate the effect of the ground plane
  user_precision_t ground_plane = 2.0*sin(2.0*M_PI*MWA_DIPOLE_HEIGHT/wavelength*cos(za));

  //Normalise the beam if requested
  if (norm == 1){
    ground_plane /= 2.0*sin(2.0*M_PI*MWA_DIPOLE_HEIGHT/wavelength);
  }

  //Used in some kind of parallatic rotation?
  double coslat = cos(latitude);
  double cosdec = cos(dec);
  double cosha = cos(ha);
  double sinlat = sin(latitude);
  double sindec = sin(dec);
  double sinha = sin(ha);

  //Some kind of parallatic rotation?
  user_precision_t rot0 = coslat*cosdec + sinlat*sindec*cosha;
  user_precision_t rot1 = -sinlat*sinha;
  user_precision_t rot2 = sindec*sinha;
  user_precision_t rot3 = cosha;

  //Normalise the ground plane to the number of dipoles??
  user_precision_t ground_plane_div_dipoles = ground_plane / NUM_DIPOLES;

  cuUserComplex pgx = gx_dip * rot0 * ground_plane_div_dipoles;
  cuUserComplex pDx = Dx_dip * rot1 * ground_plane_div_dipoles;
  cuUserComplex pDy = Dy_dip * rot2 * ground_plane_div_dipoles;
  cuUserComplex pgy = gy_dip * rot3 * ground_plane_div_dipoles;

  //Explicitly set the imag to zero, this beam is real
  pgx.y = 0;
  pDx.y = 0;
  pDy.y = 0;
  pgy.y = 0;

  * gx = pgx;
  * Dx = pDx;
  * Dy = pDy;
  * gy = pgy;

}

/*
d_azs              num_times*num_components
d_zas              num_times*num_components
d_has              num_times*num_components
d_decs             num_components
*/

__global__ void kern_RTS_analytic_MWA_beam(user_precision_t *d_azs,
           user_precision_t *d_zas,
           double *d_beam_has, double *d_beam_decs,
           double *d_metre_delays,
           double *d_freqs, double latitude, int norm,
           int num_freqs, int num_times, int num_components,
           cuUserComplex *d_gxs, cuUserComplex *d_Dxs,
           cuUserComplex *d_Dys, cuUserComplex *d_gys) {

  // Start by computing which baseline we're going to do
  const int iCoord = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iFreq = threadIdx.y + (blockDim.y*blockIdx.y);

  if (iFreq < num_freqs && iCoord < num_components * num_times) {

    int component = (int)floorf((float)iCoord / (float)num_times);
    int time_ind = iCoord - component*num_times;
    int beam_ind = num_freqs*time_ind*num_components + (num_components*iFreq) + component;

    cuUserComplex gx;
    cuUserComplex Dx;
    cuUserComplex Dy;
    cuUserComplex gy;

    double wavelength = VELC / d_freqs[iFreq];

    RTS_MWA_beam(d_azs[iCoord], d_zas[iCoord],
             d_beam_has[iCoord], d_beam_decs[iCoord],
             wavelength, d_metre_delays,
             latitude, norm,
             &gx, &Dx, &Dy, &gy);

    d_gxs[beam_ind] = gx;
    d_Dxs[beam_ind] = Dx;
    d_Dys[beam_ind] = Dy;
    d_gys[beam_ind] = gy;

  }
}

extern "C" void calculate_RTS_MWA_analytic_beam(int num_components,
     int num_time_steps, int num_freqs,
     user_precision_t *azs, user_precision_t *zas, user_precision_t *delays,
     double latitude, int norm,
     double *beam_has, double *beam_decs, double *d_freqs,
     cuUserComplex *d_gxs, cuUserComplex *d_Dxs,
     cuUserComplex *d_Dys, cuUserComplex *d_gys){

  int num_coords = num_components * num_time_steps;

  user_precision_t *d_azs = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_azs, num_coords*sizeof(user_precision_t)) );
  cudaErrorCheckCall( cudaMemcpy(d_azs, azs, num_coords*sizeof(user_precision_t),
                      cudaMemcpyHostToDevice) );

  user_precision_t *d_zas = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_zas, num_coords*sizeof(user_precision_t)) );
  cudaErrorCheckCall( cudaMemcpy(d_zas, zas, num_coords*sizeof(user_precision_t),
                      cudaMemcpyHostToDevice) );

  //Copy across stuff that normally gets copied by `source_component_common`
  double *d_beam_has = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_beam_has,
                      num_coords*sizeof(double)) );
  cudaErrorCheckCall( cudaMemcpy(d_beam_has, beam_has,
                      num_coords*sizeof(double),
                      cudaMemcpyHostToDevice) );

  double *d_beam_decs = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_beam_decs,
                      num_coords*sizeof(double)) );
  cudaErrorCheckCall( cudaMemcpy(d_beam_decs, beam_decs,
                      num_coords*sizeof(double),
                      cudaMemcpyHostToDevice) );

  //Apply the actual delay length added by MWA circuitry here, in metres - saves
  //on computation inside the kernel
  //Make a copy array so we don't modify the original
  //Have to reorder the delays as listed in the metafits to match what
  //the RTS analytic code
  double *metre_delays = (double*)malloc(NUM_DIPOLES*sizeof(double));

  for (int i = 0; i < 4; i++) {

      metre_delays[3-i] = delays[0+i*4];
      metre_delays[7-i] = delays[1+i*4];
      metre_delays[11-i] = delays[2+i*4];
      metre_delays[15-i] = delays[3+i*4];

  }

  //I have NO IDEA what this is doing, blindy copying the RTS code
  //One change to the RTS code is I take out the division by speed
  //of light, as later on the delays are multipled by speed of light again
  double delay_0 = 0.0;
  for (int k=0; k<NUM_DIPOLES; k++ ){
    delay_0 += metre_delays[k] * DQ;
  }

  delay_0 /= (double)NUM_DIPOLES;

  for(int k=0; k<NUM_DIPOLES; k++) {
    metre_delays[k] = metre_delays[k] * DQ - delay_0;
  }

  //Copy over to the GPU
  double *d_metre_delays = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_metre_delays, NUM_DIPOLES*sizeof(double)) );
  cudaErrorCheckCall( cudaMemcpy(d_metre_delays, metre_delays, NUM_DIPOLES*sizeof(double),
                      cudaMemcpyHostToDevice) );

  dim3 grid, threads;
  threads.x = 128;
  threads.y = 1;

  grid.x = (int)ceil( (float)num_coords  / (float)threads.x );
  grid.y = (int)ceil( (float)num_freqs / (float)threads.y );

  cudaErrorCheckKernel("kern_RTS_analytic_MWA_beam",
             kern_RTS_analytic_MWA_beam, grid, threads,
             d_azs, d_zas, d_beam_has, d_beam_decs, d_metre_delays,
             d_freqs, latitude, norm,
             num_freqs, num_time_steps, num_components,
             d_gxs, d_Dxs, d_Dys, d_gys);

  free(metre_delays);
  cudaErrorCheckCall( cudaFree(d_azs) );
  cudaErrorCheckCall( cudaFree(d_zas) );
  cudaErrorCheckCall( cudaFree(d_metre_delays) );
  cudaErrorCheckCall( cudaFree(d_beam_has) );
  cudaErrorCheckCall( cudaFree(d_beam_decs) );

}

/*******************************************************************************
                 Functions below to be used in unit tests
*******************************************************************************/

extern "C" void test_RTS_calculate_MWA_analytic_beam(int num_components,
     int num_time_steps, int num_freqs,
     user_precision_t *azs, user_precision_t *zas, user_precision_t *delays,
     double latitude, int norm,
     double *beam_has, double *beam_decs, double *freqs,
     user_precision_complex_t *gxs, user_precision_complex_t *Dxs,
     user_precision_complex_t *Dys, user_precision_complex_t *gys) {

  //Allocate space for beam gains
  user_precision_complex_t *d_gxs = NULL;
  user_precision_complex_t *d_Dxs = NULL;
  user_precision_complex_t *d_Dys = NULL;
  user_precision_complex_t *d_gys = NULL;

  cudaErrorCheckCall( cudaMalloc( (void**)&d_gxs,
    num_freqs*num_time_steps*num_components*sizeof(user_precision_complex_t)) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_Dxs,
    num_freqs*num_time_steps*num_components*sizeof(user_precision_complex_t)) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_Dys,
    num_freqs*num_time_steps*num_components*sizeof(user_precision_complex_t)) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_gys,
    num_freqs*num_time_steps*num_components*sizeof(user_precision_complex_t)) );

  double *d_freqs = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_freqs, num_freqs*sizeof(double)) );
  cudaErrorCheckCall( cudaMemcpy(d_freqs, freqs, num_freqs*sizeof(double),
                      cudaMemcpyHostToDevice) );

  //Run
  calculate_RTS_MWA_analytic_beam(num_components,
       num_time_steps, num_freqs,
       azs, zas, delays, latitude, norm,
       beam_has, beam_decs, d_freqs,
       (cuUserComplex*)d_gxs, (cuUserComplex*)d_Dxs,
       (cuUserComplex*)d_Dys, (cuUserComplex*)d_gys);

  cudaErrorCheckCall( cudaMemcpy(gxs, d_gxs,
       num_freqs*num_time_steps*num_components*sizeof(user_precision_complex_t),
       cudaMemcpyDeviceToHost) );
  cudaErrorCheckCall( cudaMemcpy(Dxs, d_Dxs,
       num_freqs*num_time_steps*num_components*sizeof(user_precision_complex_t),
       cudaMemcpyDeviceToHost) );
  cudaErrorCheckCall( cudaMemcpy(Dys, d_Dys,
       num_freqs*num_time_steps*num_components*sizeof(user_precision_complex_t),
       cudaMemcpyDeviceToHost) );
  cudaErrorCheckCall( cudaMemcpy(gys, d_gys,
       num_freqs*num_time_steps*num_components*sizeof(user_precision_complex_t),
       cudaMemcpyDeviceToHost) );

  cudaErrorCheckCall( cudaFree(d_gxs) );
  cudaErrorCheckCall( cudaFree(d_Dxs) );
  cudaErrorCheckCall( cudaFree(d_Dys) );
  cudaErrorCheckCall( cudaFree(d_gys) );

  cudaErrorCheckCall( cudaFree(d_freqs) );


}

extern "C" void test_analytic_dipole_beam(int num_components,
     int num_time_steps, int num_freqs,
     user_precision_t *azs, user_precision_t *zas, double *freqs,
     user_precision_complex_t *analy_beam_X,
     user_precision_complex_t *analy_beam_Y) {

  user_precision_complex_t *d_analy_beam_X = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_analy_beam_X,
    num_freqs*num_time_steps*num_components*sizeof(user_precision_complex_t)) );

  user_precision_complex_t *d_analy_beam_Y = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_analy_beam_Y,
     num_freqs*num_time_steps*num_components*sizeof(user_precision_complex_t)) );

  double *d_freqs = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_freqs, num_freqs*sizeof(double)) );
  cudaErrorCheckCall( cudaMemcpy(d_freqs, freqs, num_freqs*sizeof(double),
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
           double beam_ref_freq, double *freqs,
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

  double *d_freqs = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_freqs, num_freqs*sizeof(double) ) );
  cudaErrorCheckCall( cudaMemcpy(d_freqs, freqs,
                           num_freqs*sizeof(double), cudaMemcpyHostToDevice ) );

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
     double beam_ref_freq, double *freqs,
     double *beam_has, double *beam_decs,
     user_precision_complex_t *primay_beam_J00, user_precision_complex_t *primay_beam_J11) {

  int num_beam_hadec = num_components * num_time_steps;

  double *d_freqs = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_freqs, num_freqs*sizeof(double) ) );
  cudaErrorCheckCall( cudaMemcpy(d_freqs, freqs,
                           num_freqs*sizeof(double), cudaMemcpyHostToDevice ) );

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

// we have 4 dimensions to loop over here:
//  - component
//  - freqs
//  - times
//  - tiles (at some point in future)
//
// it's likely that time will usually be the smallest, so I'm going to
// loop over that. Time will tell whether that's a good idea or not
__global__ void kern_map_hyperbeam_gains(int num_components,
           int num_times, int num_freqs, int iTime, int num_unique_fee_freqs,
           double *d_jones, const int *d_tile_map, const int *d_freq_map,
           cuUserComplex *d_primay_beam_J00,
           cuUserComplex *d_primay_beam_J01,
           cuUserComplex *d_primay_beam_J10,
           cuUserComplex *d_primay_beam_J11) {

  //All baselines at all freqs and all times
  int iComponent = threadIdx.x + (blockDim.x*blockIdx.x);

  //The tile to map - currently unused but hopefully in the future tis possible
  int iTile = threadIdx.y + (blockDim.y*blockIdx.y);

  //freq step
  int iFreq = threadIdx.z + (blockDim.z*blockIdx.z);

  if(iComponent < num_components && iTime < num_times && iFreq < num_freqs ) {

    cuUserComplex d_beam_J00;
    cuUserComplex d_beam_J01;
    cuUserComplex d_beam_J10;
    cuUserComplex d_beam_J11;

    // For *this tile* and *this frequency*, access the de-duplicated beam
    // response.
    int i_row = d_tile_map[iTile];
    int i_col = d_freq_map[iFreq];

    int num_directions = num_times * num_components;

    int current_ind = ((num_directions * num_unique_fee_freqs * i_row) + num_directions * i_col) + iComponent*num_times + iTime;

    //Get an index to split into the WODEN style containers
    int new_ind = num_freqs*iTime*num_components + (num_components*iFreq) + iComponent;

    d_beam_J00.x = (user_precision_t)d_jones[2*MAX_POLS*current_ind + 0];
    d_beam_J00.y = (user_precision_t)d_jones[2*MAX_POLS*current_ind + 1];
    d_beam_J01.x = (user_precision_t)d_jones[2*MAX_POLS*current_ind + 2];
    d_beam_J01.y = (user_precision_t)d_jones[2*MAX_POLS*current_ind + 3];
    d_beam_J10.x = (user_precision_t)d_jones[2*MAX_POLS*current_ind + 4];
    d_beam_J10.y = (user_precision_t)d_jones[2*MAX_POLS*current_ind + 5];
    d_beam_J11.x = (user_precision_t)d_jones[2*MAX_POLS*current_ind + 6];
    d_beam_J11.y = (user_precision_t)d_jones[2*MAX_POLS*current_ind + 7];

    d_primay_beam_J00[new_ind] = d_beam_J00;
    d_primay_beam_J01[new_ind] = d_beam_J01;
    d_primay_beam_J10[new_ind] = d_beam_J10;
    d_primay_beam_J11[new_ind] = d_beam_J11;

  }
}

extern "C" void run_hyperbeam_cuda(int num_components,
           int num_time_steps, int num_freqs,
           uint8_t parallactic,
           struct FEEBeamCUDA *cuda_fee_beam,
           double *azs, double *zas,
           cuUserComplex *d_primay_beam_J00,
           cuUserComplex *d_primay_beam_J01,
           cuUserComplex *d_primay_beam_J10,
           cuUserComplex *d_primay_beam_J11){


  int num_azza = num_components * num_time_steps;
  int num_beam_values = num_azza * num_freqs;

  double *d_jones = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&(d_jones),
                      2*MAX_POLS*num_beam_values*sizeof(double)) );

  //Call hyperbeam and feed it a character array to return an error string
  //if something goes wrong
  int32_t status;
  char error_str[100];
  status = calc_jones_cuda_device(cuda_fee_beam,
                  (uint32_t)num_azza,
                  azs, zas,
                  parallactic,
                  &(d_jones),
                  error_str);

  //
  if (status != 0) {
    printf("hyperbeam error %d %s\n",status,error_str );
  }

  dim3 grid, threads;

  //For now, only using one beam
  grid.y = threads.y = 1;
  //In the future, different tiles can have different beams, so something like
  // threads.x = 32;
  // threads.y = 2;
  // grid.y = (int)ceil( (float)num_tiles / (float)threads.y );


  threads.x = 64;
  threads.z = 2;

  grid.x = (int)ceil( (float)num_components / (float)threads.x );
  grid.z = (int)ceil( (float)num_freqs / (float)threads.z );


  int num_unique_fee_freqs = get_num_unique_fee_freqs(cuda_fee_beam);
  const int *d_tile_map = get_tile_map(cuda_fee_beam);
  const int *d_freq_map = get_freq_map(cuda_fee_beam);

  for (int iTime = 0; iTime < num_time_steps; iTime++) {

    cudaErrorCheckKernel("kern_map_hyperbeam_gains",
                          kern_map_hyperbeam_gains, grid, threads,
                          num_components, num_time_steps, num_freqs,
                          iTime, num_unique_fee_freqs,
                          d_jones, d_tile_map, d_freq_map,
                          d_primay_beam_J00,
                          d_primay_beam_J01,
                          d_primay_beam_J10,
                          d_primay_beam_J11);
  }

  cudaErrorCheckCall( cudaFree(d_jones) );

}

extern "C" void test_run_hyperbeam_cuda(int num_components,
           int num_time_steps, int num_freqs,
           uint8_t parallatic,
           struct FEEBeamCUDA *cuda_fee_beam,
           double *azs, double *zas,
           user_precision_complex_t *primay_beam_J00,
           user_precision_complex_t *primay_beam_J01,
           user_precision_complex_t *primay_beam_J10,
           user_precision_complex_t *primay_beam_J11){
  uint32_t num_azza = num_components * num_time_steps;
  int num_beam_values = num_azza * num_freqs;

  user_precision_complex_t *d_primay_beam_J00 = NULL;
  user_precision_complex_t *d_primay_beam_J01 = NULL;
  user_precision_complex_t *d_primay_beam_J10 = NULL;
  user_precision_complex_t *d_primay_beam_J11 = NULL;

  cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J00,
                      num_beam_values*sizeof(user_precision_complex_t)) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J01,
                      num_beam_values*sizeof(user_precision_complex_t)) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J10,
                      num_beam_values*sizeof(user_precision_complex_t)) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J11,
                      num_beam_values*sizeof(user_precision_complex_t)) );

  run_hyperbeam_cuda(num_components,
             num_time_steps, num_freqs,
             parallatic,
             cuda_fee_beam,
             azs, zas,
             (cuUserComplex *)d_primay_beam_J00,
             (cuUserComplex *)d_primay_beam_J01,
             (cuUserComplex *)d_primay_beam_J10,
             (cuUserComplex *)d_primay_beam_J11);

  cudaErrorCheckCall( cudaMemcpy( primay_beam_J00, d_primay_beam_J00,
                      num_beam_values*sizeof(user_precision_complex_t),
                      cudaMemcpyDeviceToHost) );

  cudaErrorCheckCall( cudaMemcpy( primay_beam_J01, d_primay_beam_J01,
                      num_beam_values*sizeof(user_precision_complex_t),
                      cudaMemcpyDeviceToHost) );

  cudaErrorCheckCall( cudaMemcpy( primay_beam_J10, d_primay_beam_J10,
                      num_beam_values*sizeof(user_precision_complex_t),
                      cudaMemcpyDeviceToHost) );

  cudaErrorCheckCall( cudaMemcpy( primay_beam_J11, d_primay_beam_J11,
                      num_beam_values*sizeof(user_precision_complex_t),
                      cudaMemcpyDeviceToHost) );

  cudaErrorCheckCall( cudaFree(d_primay_beam_J00) );
  cudaErrorCheckCall( cudaFree(d_primay_beam_J01) );
  cudaErrorCheckCall( cudaFree(d_primay_beam_J10) );
  cudaErrorCheckCall( cudaFree(d_primay_beam_J11) );

}
