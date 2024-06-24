#include <stdio.h>
#include <stdlib.h>
#include <cuComplex.h>
#include <complex.h>
#include <math.h>
#include "constants.h"
#include "fundamental_coords.h"
#include "cudacheck.h"
#include "primary_beam_cuda.h"

//__host__??
__device__ void twoD_Gaussian(user_precision_t x, user_precision_t y,
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
           cuUserComplex *d_g1xs, cuUserComplex *d_g1ys) {
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

    d_g1xs[beam_ind] = make_cuUserComplex(d_beam_real, d_beam_imag);
    d_g1ys[beam_ind] = make_cuUserComplex(d_beam_real, d_beam_imag);

  }
}


extern "C" void calculate_gaussian_beam(int num_components, int num_time_steps,
           int num_freqs, user_precision_t ha0,
           user_precision_t sdec0, user_precision_t cdec0,
           user_precision_t fwhm_lm, user_precision_t cos_theta,
           user_precision_t sin_theta, user_precision_t sin_2theta,
           double beam_ref_freq, double *d_freqs,
           double *beam_has, double *beam_decs,
           cuUserComplex *d_g1xs, cuUserComplex *d_g1ys){

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
                       d_g1xs, d_g1ys);

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
           cuUserComplex *d_g1xs, cuUserComplex *d_g1ys) {
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

    d_g1xs[beam_ind] = normed_X;
    d_g1ys[beam_ind] = normed_Y;

  }
}

extern "C" void calculate_analytic_dipole_beam(int num_components,
     int num_time_steps, int num_freqs,
     user_precision_t *azs, user_precision_t *zas, double *d_freqs,
     cuUserComplex *d_g1xs, cuUserComplex *d_g1ys){

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
             d_g1xs, d_g1ys);

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

      metre_delays[3-i] = (double)delays[0+i*4];
      metre_delays[7-i] = (double)delays[1+i*4];
      metre_delays[11-i] = (double)delays[2+i*4];
      metre_delays[15-i] = (double)delays[3+i*4];

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


// we have 3 dimensions to loop over here:
//  - component
//  - freqs
//  - times
//
// it's likely that time will usually be the smallest, so I'm going to
// loop over that. Time will tell whether that's a good idea or not LOL
__global__ void kern_map_hyperbeam_gains(int num_components,
           int num_times, int num_freqs, int iTime, int num_unique_fee_freqs,
           double *d_jones, const int *d_tile_map, const int *d_freq_map,
           int parallactic,
           cuUserComplex *d_gxs,
           cuUserComplex *d_Dxs,
           cuUserComplex *d_Dys,
           cuUserComplex *d_gys) {

  //All baselines at all freqs and all times
  int iComponent = threadIdx.x + (blockDim.x*blockIdx.x);

  //The tile to map - currently unused but hopefully in the future tis possible
  int iTile = threadIdx.y + (blockDim.y*blockIdx.y);

  //freq step
  int iFreq = threadIdx.z + (blockDim.z*blockIdx.z);

  if (iComponent < num_components && iTime < num_times && iFreq < num_freqs ) {

    cuUserComplex d_beam_J00;
    cuUserComplex d_beam_J01;
    cuUserComplex d_beam_J10;
    cuUserComplex d_beam_J11;

    // For *this tile* and *this frequency*, access the de-duplicated beam
    // response.
    int i_row = d_tile_map[iTile];
    int i_col = d_freq_map[iFreq];

    int hyper_ind, current_ind;

    //If we are doing parallactic rotation corrections, the latitude needs
    //to be fed in. As this change change with time, we call hyperbeam for
    //each time step, and use pointer arithmatic when storing the outputs
    //This means the striping for retrieving outputs is different.
    if (parallactic == 1)
    {
      //Where the desired data sits within a single time step hyperdrive
      //output
      hyper_ind = ((num_components * num_unique_fee_freqs * i_row) + num_components * i_col);

      //Add on how many previous time-steps-worth of outputs we've
      //already remapped
      current_ind = hyper_ind + iTime*num_components*num_freqs + iComponent;
    } else {
      hyper_ind = ((num_components*num_times * num_unique_fee_freqs * i_row) + num_components * i_col);

      current_ind = hyper_ind + iTime*num_components + iComponent;
    }

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

    d_gxs[new_ind] = d_beam_J00;
    d_Dxs[new_ind] = d_beam_J01;
    d_Dys[new_ind] = d_beam_J10;
    d_gys[new_ind] = d_beam_J11;

  }
}



// we have 4 dimensions to loop over here:
//  - component
//  - freqs
//  - times
//  - tiles
//
// it's likely that time will usually be the smallest, so I'm going to
// loop over that. Time will tell whether that's a good idea or not LOL
__global__ void kern_map_hyperbeam_gains_multi_antennas(int num_components,
           int num_times, int num_freqs, int num_tiles, int iTime, int num_unique_fee_freqs,
           double *d_jones, const int *d_tile_map, const int *d_freq_map,
           int parallactic,
           cuUserComplex *d_gxs_ants,
           cuUserComplex *d_Dxs_ants,
           cuUserComplex *d_Dys_ants,
           cuUserComplex *d_gys_ants) {

  //All baselines at all freqs and all times
  int iComponent = threadIdx.x + (blockDim.x*blockIdx.x);

  //The tile to map - currently unused but hopefully in the future tis possible
  int iTile = threadIdx.y + (blockDim.y*blockIdx.y);

  //freq step
  int iFreq = threadIdx.z + (blockDim.z*blockIdx.z);

  if(iComponent < num_components && iTime < num_times && iFreq < num_freqs && iTile < num_tiles) {

    cuUserComplex d_beam_J00;
    cuUserComplex d_beam_J01;
    cuUserComplex d_beam_J10;
    cuUserComplex d_beam_J11;

    // For *this tile* and *this frequency*, access the de-duplicated beam
    // response.
    int i_row = d_tile_map[iTile];
    int i_col = d_freq_map[iFreq];

    int hyper_ind, current_ind;

    //If we are doing parallactic rotation corrections, the latitude needs
    //to be fed in. As this change change with time, we call hyperbeam for
    //each time step, and use pointer arithmatic when storing the outputs
    //This means the striping for retrieving outputs is different.
    if (parallactic == 1)
    {
      //Where the desired data sits within a single time step hyperdrive
      //output
      hyper_ind = ((num_components*num_unique_fee_freqs*i_row) + num_components * i_col);

      //Add on how many previous time-steps-worth of outputs we've
      //already remapped
      current_ind = hyper_ind + iTime*num_components*num_freqs*num_tiles + iComponent;
    } else {
      hyper_ind = ((num_components*num_times*num_unique_fee_freqs*i_row) + num_components * i_col);

      current_ind = hyper_ind + iTime*num_components*num_tiles + iComponent;
    }

    d_beam_J00.x = (user_precision_t)d_jones[2*MAX_POLS*current_ind + 0];
    d_beam_J00.y = (user_precision_t)d_jones[2*MAX_POLS*current_ind + 1];
    d_beam_J01.x = (user_precision_t)d_jones[2*MAX_POLS*current_ind + 2];
    d_beam_J01.y = (user_precision_t)d_jones[2*MAX_POLS*current_ind + 3];
    d_beam_J10.x = (user_precision_t)d_jones[2*MAX_POLS*current_ind + 4];
    d_beam_J10.y = (user_precision_t)d_jones[2*MAX_POLS*current_ind + 5];
    d_beam_J11.x = (user_precision_t)d_jones[2*MAX_POLS*current_ind + 6];
    d_beam_J11.y = (user_precision_t)d_jones[2*MAX_POLS*current_ind + 7];


    //Get an index to split into the WODEN style containers
    int new_ind =  num_times*num_freqs*num_components*iTile + num_freqs*iTime*num_components + (num_components*iFreq) + iComponent;

    d_gxs_ants[new_ind] = d_beam_J00;
    d_Dxs_ants[new_ind] = d_beam_J01;
    d_Dys_ants[new_ind] = d_beam_J10;
    d_gys_ants[new_ind] = d_beam_J11;

  }
}


__global__ void fill_with_ones(int num_azza, double *d_jones) {

  //All baselines at all freqs and all times
  int iAzza = threadIdx.x + (blockDim.x*blockIdx.x);

  if (iAzza < num_azza) {
    d_jones[iAzza] = 1;
  }
}


//TODO need this function to take in a second set of jones matrices
extern "C" void run_hyperbeam_cuda(int num_components,
           int num_time_steps, int num_freqs,
           uint8_t parallactic,
           struct FEEBeamGpu *cuda_fee_beam,
           double *azs, double *zas,
           double *latitudes,
           cuUserComplex *d_gxs,
           cuUserComplex *d_Dxs,
           cuUserComplex *d_Dys,
           cuUserComplex *d_gys){


  int num_azza = num_components * num_time_steps;
  int num_beam_values = num_azza * num_freqs;

  //TODO if not doing parallactic, we don't have to malloc this much memory;
  //only enough for ONE time step. Once we get into a beam for each tile,
  //this will be a useful memory save. For now, not too expensive, so leave
  //as is
  double *d_jones = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&(d_jones),
                      2*MAX_POLS*num_beam_values*sizeof(double)) );

  int32_t status = 0;
  // The beam-response Jones matrix be in the IAU polarisation order
  int iau_order = 1;

  //This is used for the remapping below
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
  int num_unique_fee_tiles = get_num_unique_fee_tiles(cuda_fee_beam);

  int num_tiles = 1;

  //Get host pointers to the tile and freq maps
  const int *tile_map;
  const int *freq_map;

  tile_map = get_fee_tile_map(cuda_fee_beam);
  freq_map = get_fee_freq_map(cuda_fee_beam);

  //Copy the tile and freq maps to the GPU
  int32_t *d_tile_map = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&(d_tile_map),
                      num_tiles*sizeof(int32_t) ) );
  int32_t *d_freq_map = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&(d_freq_map),
                      num_freqs*sizeof(int32_t) ) );

  cudaErrorCheckCall( cudaMemcpy(d_tile_map, tile_map,
            num_tiles*sizeof(int32_t), cudaMemcpyHostToDevice ) );
  cudaErrorCheckCall( cudaMemcpy(d_freq_map, freq_map,
            num_freqs*sizeof(int32_t), cudaMemcpyHostToDevice ) );

  if (parallactic) {

    //The latitude of the array should change with time if we are precessing
    //it back to J2000. This means we have to call hyperbeam for as many time
    //steps as we have. Supply chunks of az, za, and d_jones via pointer
    //arithmatic
    int increment;
    for (int time_ind = 0; time_ind < num_time_steps; time_ind++) {

      //The size of the increment depends on the number of unique
      //frequencies - 
      increment = time_ind*num_components;
      //You get real + imag for MAX_POLS polarisations hence the larger
      //increment on the d_jones array below
      status = fee_calc_jones_gpu_device(cuda_fee_beam,
                      (uint32_t)num_components,
                      azs + increment, zas + increment,
                      &latitudes[time_ind],
                      iau_order,
                      (double *)d_jones + 2*MAX_POLS*num_freqs*increment);
      //TODO this needs to take in a BASELINE mapping function, so we can
      //get the correct tile going into d_gx_1 and d_gx_2
      cudaErrorCheckKernel("kern_map_hyperbeam_gains",
                            kern_map_hyperbeam_gains, grid, threads,
                            num_components, num_time_steps, num_freqs,
                            time_ind, num_unique_fee_freqs,
                            d_jones, d_tile_map, d_freq_map,
                            parallactic,
                            d_gxs,
                            d_Dxs,
                            d_Dys,
                            d_gys);
    }

  }
  else {
    status = fee_calc_jones_gpu_device(cuda_fee_beam,
                    (uint32_t)num_azza,
                    azs, zas,
                    NULL,
                    iau_order,
                    (double *)d_jones);

    for (int iTime = 0; iTime < num_time_steps; iTime++) {

      cudaErrorCheckKernel("kern_map_hyperbeam_gains",
                            kern_map_hyperbeam_gains, grid, threads,
                            num_components, num_time_steps, num_freqs,
                            iTime, num_unique_fee_freqs,
                            d_jones, d_tile_map, d_freq_map,
                            parallactic,
                            d_gxs,
                            d_Dxs,
                            d_Dys,
                            d_gys);
    }

  }

  if (status != 0) {
    const char *func_name = "fee_calc_jones_gpu_device";
    handle_hyperbeam_error(__FILE__, __LINE__, func_name);
  }

  cudaErrorCheckCall( cudaFree(d_jones) );
  cudaErrorCheckCall( cudaFree(d_tile_map) );
  cudaErrorCheckCall( cudaFree(d_freq_map) );

}




//TODO need this function to take in a second set of jones matrices
extern "C" void run_hyperbeam_cuda_multi_ants(int num_components,
           int num_time_steps, int num_freqs,
           int num_ants,
           uint8_t parallactic,
           struct FEEBeamGpu *cuda_fee_beam,
           double *azs, double *zas,
           double *latitudes,
           cuUserComplex *d_gxs_ants,
           cuUserComplex *d_Dxs_ants,
           cuUserComplex *d_Dys_ants,
           cuUserComplex *d_gys_ants){


  int num_azza = num_components * num_time_steps;
  int num_beam_values = num_azza * num_freqs * num_ants;

  //TODO if not doing parallactic, we don't have to malloc this much memory;
  //only enough for ONE time step. Once we get into a beam for each tile,
  //this will be a useful memory save. For now, not too expensive, so leave
  //as is
  double *d_jones = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&(d_jones),
                      2*MAX_POLS*num_beam_values*sizeof(double)) );

  int32_t status = 0;
  // Should the beam-response Jones matrix be in the IAU polarisation order?
  int iau_order;

  //This is used for the remapping below
  dim3 grid, threads;

  //In the future, different tiles can have different beams, so something like
  threads.x = 32;
  threads.y = 2;
  threads.z = 2;

  grid.x = (int)ceil( (float)num_components / (float)threads.x );
  grid.y = (int)ceil( (float)num_ants / (float)threads.y );
  grid.z = (int)ceil( (float)num_freqs / (float)threads.z );

  int num_unique_fee_freqs;
  const int *d_tile_map;
  const int *d_freq_map;

  if (parallactic) {
    iau_order = 1;

    //The latitude of the array should change with time if we are precessing
    //it back to J2000. This means we have to call hyperbeam for as many time
    //steps as we have. Supply chunks of az, za, and d_jones via pointer
    //arithmatic
    int increment;
    for (int time_ind = 0; time_ind < num_time_steps; time_ind++) {

      //The size of the increment depends on the number of unique
      //frequencies - 
      increment = time_ind*num_components;
      //You get real + imag for MAX_POLS polarisations hence the larger
      //increment on the d_jones array below
      status = fee_calc_jones_gpu_device(cuda_fee_beam,
                      (uint32_t)num_components,
                      azs + increment, zas + increment,
                      &latitudes[time_ind],
                      iau_order,
                      (double *)d_jones + 2*MAX_POLS*num_freqs*num_ants*increment);

      num_unique_fee_freqs = get_num_unique_fee_freqs(cuda_fee_beam);
      d_tile_map = get_fee_tile_map(cuda_fee_beam);
      d_freq_map = get_fee_freq_map(cuda_fee_beam);


      //TODO this needs to take in a BASELINE mapping function, so we can
      //get the correct tile going into d_gx_1 and d_gx_2
      cudaErrorCheckKernel("kern_map_hyperbeam_gains_multi_antennas",
                            kern_map_hyperbeam_gains_multi_antennas, grid, threads,
                            num_components, num_time_steps, num_freqs, num_ants,
                            time_ind, num_unique_fee_freqs,
                            d_jones, d_tile_map, d_freq_map,
                            parallactic,
                            d_gxs_ants,
                            d_Dxs_ants,
                            d_Dys_ants,
                            d_gys_ants);
    }

  }
  else {
    int iau_order = 1;
    status = fee_calc_jones_gpu_device(cuda_fee_beam,
                    (uint32_t)num_azza,
                    azs, zas,
                    NULL,
                    iau_order,
                    (double *)d_jones);

    num_unique_fee_freqs = get_num_unique_fee_freqs(cuda_fee_beam);
    d_tile_map = get_fee_tile_map(cuda_fee_beam);
    d_freq_map = get_fee_freq_map(cuda_fee_beam);

    // for (int iTime = 0; iTime < num_time_steps; iTime++) {

    //   cudaErrorCheckKernel("kern_map_hyperbeam_gains",
    //                         kern_map_hyperbeam_gains, grid, threads,
    //                         num_components, num_time_steps, num_freqs,
    //                         iTime, num_unique_fee_freqs,
    //                         d_jones, d_tile_map, d_freq_map,
    //                         parallactic,
    //                         d_gxs,
    //                         d_Dxs,
    //                         d_Dys,
    //                         d_gys);
    // }

  }

  if (status != 0) {
    const char *func_name = "fee_calc_jones_gpu_device";
    handle_hyperbeam_error(__FILE__, __LINE__, func_name);
  }

  cudaErrorCheckCall( cudaFree(d_jones) );

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

  user_precision_complex_t *d_g1xs = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_g1xs,
                                   num_freqs*num_beam_hadec*sizeof(user_precision_complex_t)) );

  user_precision_complex_t *d_g1ys = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_g1ys,
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
                        (cuUserComplex *)d_g1xs,
                        (cuUserComplex *)d_g1ys);

  cudaErrorCheckCall( cudaMemcpy(primay_beam_J00, d_g1xs,
                 num_freqs*num_beam_hadec*sizeof(user_precision_complex_t), cudaMemcpyDeviceToHost) );
  cudaErrorCheckCall( cudaMemcpy(primay_beam_J11, d_g1ys,
                 num_freqs*num_beam_hadec*sizeof(user_precision_complex_t), cudaMemcpyDeviceToHost) );

  cudaErrorCheckCall( cudaFree(d_beam_ls ) );
  cudaErrorCheckCall( cudaFree(d_beam_ms ) );
  cudaErrorCheckCall( cudaFree(d_freqs ) );
  cudaErrorCheckCall( cudaFree(d_g1xs ) );
  cudaErrorCheckCall( cudaFree(d_g1ys ) );

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

  user_precision_complex_t *d_g1xs = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_g1xs,
                                   num_freqs*num_beam_hadec*sizeof(user_precision_complex_t)) );

  user_precision_complex_t *d_g1ys = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_g1ys,
                                   num_freqs*num_beam_hadec*sizeof(user_precision_complex_t)) );


  calculate_gaussian_beam(num_components, num_time_steps,
                         num_freqs, ha0, sdec0, cdec0,
                         fwhm_lm, cos_theta, sin_theta, sin_2theta,
                         beam_ref_freq, d_freqs,
                         beam_has, beam_decs,
                         (cuUserComplex *)d_g1xs,
                         (cuUserComplex *)d_g1ys);

  cudaErrorCheckCall( cudaMemcpy(primay_beam_J00, d_g1xs,
                 num_freqs*num_beam_hadec*sizeof(user_precision_complex_t), cudaMemcpyDeviceToHost) );
  cudaErrorCheckCall( cudaMemcpy(primay_beam_J11, d_g1ys,
                 num_freqs*num_beam_hadec*sizeof(user_precision_complex_t), cudaMemcpyDeviceToHost) );

  cudaErrorCheckCall( cudaFree(d_freqs ) );
  cudaErrorCheckCall( cudaFree(d_g1xs ) );
  cudaErrorCheckCall( cudaFree(d_g1ys ) );

}

extern "C" void test_run_hyperbeam_cuda(int num_components,
           int num_time_steps, int num_freqs,
           uint8_t parallatic,
           struct FEEBeamGpu *cuda_fee_beam,
           double *azs, double *zas,
           double *latitudes,
           user_precision_complex_t *primay_beam_J00,
           user_precision_complex_t *primay_beam_J01,
           user_precision_complex_t *primay_beam_J10,
           user_precision_complex_t *primay_beam_J11){
  uint32_t num_azza = num_components * num_time_steps;
  int num_beam_values = num_azza * num_freqs;

  user_precision_complex_t *d_gxs = NULL;
  user_precision_complex_t *d_Dxs = NULL;
  user_precision_complex_t *d_Dys = NULL;
  user_precision_complex_t *d_gys = NULL;

  cudaErrorCheckCall( cudaMalloc( (void**)&d_gxs,
                      num_beam_values*sizeof(user_precision_complex_t)) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_Dxs,
                      num_beam_values*sizeof(user_precision_complex_t)) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_Dys,
                      num_beam_values*sizeof(user_precision_complex_t)) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_gys,
                      num_beam_values*sizeof(user_precision_complex_t)) );

  double *reordered_azs = (double *)malloc(num_azza*sizeof(double));
  double *reordered_zas = (double *)malloc(num_azza*sizeof(double));

  int stripe_new, stripe_old;

  for (int time_ind = 0; time_ind < num_time_steps; time_ind++) {
    for (int comp_ind = 0; comp_ind < num_components; comp_ind++) {
      stripe_new = time_ind*num_components + comp_ind;
      stripe_old = comp_ind*num_time_steps + time_ind;
      reordered_azs[stripe_new] = azs[stripe_old];
      reordered_zas[stripe_new] = zas[stripe_old];
    }
  }

  run_hyperbeam_cuda(num_components,
             num_time_steps, num_freqs,
             parallatic,
             cuda_fee_beam,
             reordered_azs, reordered_zas,
             latitudes,
             (cuUserComplex *)d_gxs,
             (cuUserComplex *)d_Dxs,
             (cuUserComplex *)d_Dys,
             (cuUserComplex *)d_gys);

  cudaErrorCheckCall( cudaMemcpy( primay_beam_J00, d_gxs,
                      num_beam_values*sizeof(user_precision_complex_t),
                      cudaMemcpyDeviceToHost) );

  cudaErrorCheckCall( cudaMemcpy( primay_beam_J01, d_Dxs,
                      num_beam_values*sizeof(user_precision_complex_t),
                      cudaMemcpyDeviceToHost) );

  cudaErrorCheckCall( cudaMemcpy( primay_beam_J10, d_Dys,
                      num_beam_values*sizeof(user_precision_complex_t),
                      cudaMemcpyDeviceToHost) );

  cudaErrorCheckCall( cudaMemcpy( primay_beam_J11, d_gys,
                      num_beam_values*sizeof(user_precision_complex_t),
                      cudaMemcpyDeviceToHost) );

  cudaErrorCheckCall( cudaFree(d_gxs) );
  cudaErrorCheckCall( cudaFree(d_Dxs) );
  cudaErrorCheckCall( cudaFree(d_Dys) );
  cudaErrorCheckCall( cudaFree(d_gys) );

  free(reordered_azs);
  free(reordered_zas);

}





extern "C" void test_run_hyperbeam_cuda_multi_ants(int num_components,
           int num_time_steps, int num_freqs, int num_ants,
           uint8_t parallatic, 
           struct FEEBeamGpu *cuda_fee_beam,
           double *azs, double *zas,
           double *latitudes,
           user_precision_complex_t *primay_beam_J00,
           user_precision_complex_t *primay_beam_J01,
           user_precision_complex_t *primay_beam_J10,
           user_precision_complex_t *primay_beam_J11){
  uint32_t num_azza = num_components * num_time_steps;
  int num_beam_values = num_azza * num_freqs * num_ants;

  user_precision_complex_t *d_gxs_ants = NULL;
  user_precision_complex_t *d_Dxs_ants = NULL;
  user_precision_complex_t *d_Dys_ants = NULL;
  user_precision_complex_t *d_gys_ants = NULL;

  cudaErrorCheckCall( cudaMalloc( (void**)&d_gxs_ants,
                      num_beam_values*sizeof(user_precision_complex_t)) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_Dxs_ants,
                      num_beam_values*sizeof(user_precision_complex_t)) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_Dys_ants,
                      num_beam_values*sizeof(user_precision_complex_t)) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_gys_ants,
                      num_beam_values*sizeof(user_precision_complex_t)) );

  double *reordered_azs = (double *)malloc(num_azza*sizeof(double));
  double *reordered_zas = (double *)malloc(num_azza*sizeof(double));

  int stripe_new, stripe_old;

  for (int time_ind = 0; time_ind < num_time_steps; time_ind++) {
    for (int comp_ind = 0; comp_ind < num_components; comp_ind++) {
      stripe_new = time_ind*num_components + comp_ind;
      stripe_old = comp_ind*num_time_steps + time_ind;
      reordered_azs[stripe_new] = azs[stripe_old];
      reordered_zas[stripe_new] = zas[stripe_old];
    }
  }

  run_hyperbeam_cuda_multi_ants(num_components,
             num_time_steps, num_freqs, num_ants,
             parallatic,
             cuda_fee_beam,
             reordered_azs, reordered_zas,
             latitudes,
             (cuUserComplex *)d_gxs_ants,
             (cuUserComplex *)d_Dxs_ants,
             (cuUserComplex *)d_Dys_ants,
             (cuUserComplex *)d_gys_ants);

  cudaErrorCheckCall( cudaMemcpy( primay_beam_J00, d_gxs_ants,
                      num_beam_values*sizeof(user_precision_complex_t),
                      cudaMemcpyDeviceToHost) );

  cudaErrorCheckCall( cudaMemcpy( primay_beam_J01, d_Dxs_ants,
                      num_beam_values*sizeof(user_precision_complex_t),
                      cudaMemcpyDeviceToHost) );

  cudaErrorCheckCall( cudaMemcpy( primay_beam_J10, d_Dys_ants,
                      num_beam_values*sizeof(user_precision_complex_t),
                      cudaMemcpyDeviceToHost) );

  cudaErrorCheckCall( cudaMemcpy( primay_beam_J11, d_gys_ants,
                      num_beam_values*sizeof(user_precision_complex_t),
                      cudaMemcpyDeviceToHost) );

  cudaErrorCheckCall( cudaFree(d_gxs_ants) );
  cudaErrorCheckCall( cudaFree(d_Dxs_ants) );
  cudaErrorCheckCall( cudaFree(d_Dys_ants) );
  cudaErrorCheckCall( cudaFree(d_gys_ants) );

  free(reordered_azs);
  free(reordered_zas);

}