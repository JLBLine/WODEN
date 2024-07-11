#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include "constants.h"
#include "fundamental_coords.h"
#include "primary_beam_gpu.h"
#include "gpu_macros.h"

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
           gpuUserComplex *d_g1xs, gpuUserComplex *d_g1ys) {
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

    d_g1xs[beam_ind] = make_gpuUserComplex(d_beam_real, d_beam_imag);
    d_g1ys[beam_ind] = make_gpuUserComplex(d_beam_real, d_beam_imag);

  }
}


extern "C" void calculate_gaussian_beam(int num_components, int num_time_steps,
           int num_freqs, user_precision_t ha0,
           user_precision_t sdec0, user_precision_t cdec0,
           user_precision_t fwhm_lm, user_precision_t cos_theta,
           user_precision_t sin_theta, user_precision_t sin_2theta,
           double beam_ref_freq, double *d_freqs,
           double *beam_has, double *beam_decs,
           gpuUserComplex *d_g1xs, gpuUserComplex *d_g1ys){

  int num_beam_hadec = num_components * num_time_steps;

  double *d_beam_has = NULL;
  ( gpuMalloc( (void**)&d_beam_has, num_beam_hadec*sizeof(double)) );
  ( gpuMemcpy( d_beam_has, beam_has,
                      num_beam_hadec*sizeof(double), gpuMemcpyHostToDevice) );

  double *d_beam_decs = NULL;
  ( gpuMalloc( (void**)&d_beam_decs, num_beam_hadec*sizeof(double)) );
  ( gpuMemcpy( d_beam_decs, beam_decs,
                      num_beam_hadec*sizeof(double), gpuMemcpyHostToDevice) );

  double *d_beam_ls = NULL;
  ( gpuMalloc( (void**)&d_beam_ls, num_beam_hadec*sizeof(double)) );

  double *d_beam_ms = NULL;
  ( gpuMalloc( (void**)&d_beam_ms, num_beam_hadec*sizeof(double)) );

  double *d_beam_ns = NULL;
  ( gpuMalloc( (void**)&d_beam_ns, num_beam_hadec*sizeof(double)) );

  dim3 grid, threads;
  threads.x = 128;
  grid.x = (int)ceil( (user_precision_t)num_beam_hadec / (user_precision_t)threads.x );

  gpuErrorCheckKernel("kern_calc_lmn", kern_calc_lmn, grid, threads,
                        ha0, sdec0, cdec0,
                        d_beam_has, d_beam_decs,
                        d_beam_ls, d_beam_ms, d_beam_ns, num_beam_hadec);

  threads.y = 2;
  grid.x = (int)ceil( (float)num_time_steps*(float)num_components / (float)threads.x );
  grid.y = (int)ceil( (float)num_freqs / (float)threads.y );
  //
  gpuErrorCheckKernel("kern_gaussian_beam",kern_gaussian_beam, grid, threads,
                       d_beam_ls, d_beam_ms, beam_ref_freq, d_freqs,
                       fwhm_lm, cos_theta, sin_theta, sin_2theta,
                       num_freqs, num_time_steps, num_components,
                       d_g1xs, d_g1ys);

  gpuFree( d_beam_ns );
  gpuFree( d_beam_ms );
  gpuFree( d_beam_ls );
  gpuFree( d_beam_decs );
  gpuFree( d_beam_has );
}

__device__ void analytic_dipole(user_precision_t az, user_precision_t za,
           user_precision_t wavelength,
           gpuUserComplex * d_beam_X, gpuUserComplex * d_beam_Y) {

  user_precision_t dipole_height_m = 0.3;

  //Here, X means a north-south aligned dipole
  //      Y means an east-west aligned dipole

  user_precision_t theta_parallel_X = acos(sin(za)*cos(az));
  user_precision_t theta_parallel_Y = acos(sin(za)*sin(az));

  user_precision_t d_in_lambda = (2. * dipole_height_m)/wavelength;
  user_precision_t gp_effect_array = 2. * sin(M_PI*d_in_lambda*cos(za));

  user_precision_t voltage_parallel_X = sin(theta_parallel_X) * gp_effect_array;
  user_precision_t voltage_parallel_Y = sin(theta_parallel_Y) * gp_effect_array;

  gpuUserComplex tempX;
  gpuUserComplex tempY;

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
           gpuUserComplex *d_g1xs, gpuUserComplex *d_g1ys) {
  // Start by computing which baseline we're going to do
  const int iCoord = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iFreq = threadIdx.y + (blockDim.y*blockIdx.y);

  // TODO would it be faster to have freq and time on one axis, and num
  // components on other axis?
  if (iFreq < num_freqs && iCoord < num_components * num_times) {

    int component = (int)floorf((float)iCoord / (float)num_times);
    int time_ind = iCoord - component*num_times;

    int beam_ind = num_freqs*time_ind*num_components + (num_components*iFreq) + component;

    gpuUserComplex d_beam_X, d_beam_Y;
    gpuUserComplex d_beam_norm_X, d_beam_norm_Y;
    user_precision_t wavelength = VELC / d_freqs[iFreq];

    analytic_dipole(d_azs[iCoord], d_zas[iCoord], wavelength,
               &d_beam_X, &d_beam_Y);

    //Should really calculate this normalisation outside of this kernel - we are
    //massively increasing the number calculations for no reason
    analytic_dipole(0.0, 0.0, wavelength,
               &d_beam_norm_X, &d_beam_norm_Y);

    gpuUserComplex normed_X = d_beam_X;
    gpuUserComplex normed_Y = d_beam_Y;

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
     gpuUserComplex *d_g1xs, gpuUserComplex *d_g1ys){

  int num_beam_azza = num_components * num_time_steps;

  user_precision_t *d_azs = NULL;
  ( gpuMalloc( (void**)&d_azs, num_beam_azza*sizeof(user_precision_t)) );
  ( gpuMemcpy(d_azs, azs, num_beam_azza*sizeof(user_precision_t),
                      gpuMemcpyHostToDevice) );

  user_precision_t *d_zas = NULL;
  ( gpuMalloc( (void**)&d_zas, num_beam_azza*sizeof(user_precision_t)) );
  ( gpuMemcpy(d_zas, zas, num_beam_azza*sizeof(user_precision_t),
                      gpuMemcpyHostToDevice) );

  dim3 grid, threads;
  threads.x = 128;
  threads.y = 1;

  grid.x = (int)ceil( (float)num_beam_azza  / (float)threads.x );
  grid.y = (int)ceil( (float)num_freqs / (float)threads.y );

  gpuErrorCheckKernel("kern_analytic_dipole_beam",
             kern_analytic_dipole_beam, grid, threads,
             d_azs, d_zas, d_freqs, num_freqs,
             num_time_steps, num_components,
             d_g1xs, d_g1ys);

  ( gpuFree(d_azs) );
  ( gpuFree(d_zas) );

}

__device__ void RTS_MWA_beam(user_precision_t az, user_precision_t za,
           double ha, double dec,
           double wavelength, double *d_metre_delays,
           double latitude, int norm,
           gpuUserComplex * gx, gpuUserComplex * Dx,
           gpuUserComplex * Dy, gpuUserComplex * gy) {

  // set elements of the look-dir vector
  double proj_e = sin(za)*sin(az);
  double proj_n = sin(za)*cos(az);
  double proj_z = cos(za);

  int n_cols = 4;
  int n_rows = 4;

  //Used in calculating the phase later on
  double multiplier = -2 * M_PI / wavelength;
  double dipl_e, dipl_n, dipl_z;

  gpuUserComplex x_dip;
  gpuUserComplex y_dip;

  gpuUserComplex gx_dip = {0.0, 0.0};
  gpuUserComplex Dx_dip = {0.0, 0.0};
  gpuUserComplex Dy_dip = {0.0, 0.0};
  gpuUserComplex gy_dip = {0.0, 0.0};

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

  gpuUserComplex pgx = gx_dip * rot0 * ground_plane_div_dipoles;
  gpuUserComplex pDx = Dx_dip * rot1 * ground_plane_div_dipoles;
  gpuUserComplex pDy = Dy_dip * rot2 * ground_plane_div_dipoles;
  gpuUserComplex pgy = gy_dip * rot3 * ground_plane_div_dipoles;

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
           gpuUserComplex *d_gxs, gpuUserComplex *d_Dxs,
           gpuUserComplex *d_Dys, gpuUserComplex *d_gys) {

  // Start by computing which baseline we're going to do
  const int iCoord = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iFreq = threadIdx.y + (blockDim.y*blockIdx.y);

  if (iFreq < num_freqs && iCoord < num_components * num_times) {

    int component = (int)floorf((float)iCoord / (float)num_times);
    int time_ind = iCoord - component*num_times;
    int beam_ind = num_freqs*time_ind*num_components + (num_components*iFreq) + component;

    gpuUserComplex gx;
    gpuUserComplex Dx;
    gpuUserComplex Dy;
    gpuUserComplex gy;

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
     user_precision_t *azs, user_precision_t *zas, int *delays,
     double latitude, int norm,
     double *beam_has, double *beam_decs, double *d_freqs,
     gpuUserComplex *d_gxs, gpuUserComplex *d_Dxs,
     gpuUserComplex *d_Dys, gpuUserComplex *d_gys){

  int num_coords = num_components * num_time_steps;

  user_precision_t *d_azs = NULL;
  ( gpuMalloc( (void**)&d_azs, num_coords*sizeof(user_precision_t)) );
  ( gpuMemcpy(d_azs, azs, num_coords*sizeof(user_precision_t),
                      gpuMemcpyHostToDevice) );

  user_precision_t *d_zas = NULL;
  ( gpuMalloc( (void**)&d_zas, num_coords*sizeof(user_precision_t)) );
  ( gpuMemcpy(d_zas, zas, num_coords*sizeof(user_precision_t),
                      gpuMemcpyHostToDevice) );

  //Copy across stuff that normally gets copied by `source_component_common`
  double *d_beam_has = NULL;
  ( gpuMalloc( (void**)&d_beam_has,
                      num_coords*sizeof(double)) );
  ( gpuMemcpy(d_beam_has, beam_has,
                      num_coords*sizeof(double),
                      gpuMemcpyHostToDevice) );

  double *d_beam_decs = NULL;
  ( gpuMalloc( (void**)&d_beam_decs,
                      num_coords*sizeof(double)) );
  ( gpuMemcpy(d_beam_decs, beam_decs,
                      num_coords*sizeof(double),
                      gpuMemcpyHostToDevice) );

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
  ( gpuMalloc( (void**)&d_metre_delays, NUM_DIPOLES*sizeof(double)) );
  ( gpuMemcpy(d_metre_delays, metre_delays, NUM_DIPOLES*sizeof(double),
                      gpuMemcpyHostToDevice) );

  dim3 grid, threads;
  threads.x = 128;
  threads.y = 1;

  grid.x = (int)ceil( (float)num_coords  / (float)threads.x );
  grid.y = (int)ceil( (float)num_freqs / (float)threads.y );

  gpuErrorCheckKernel("kern_RTS_analytic_MWA_beam",
             kern_RTS_analytic_MWA_beam, grid, threads,
             d_azs, d_zas, d_beam_has, d_beam_decs, d_metre_delays,
             d_freqs, latitude, norm,
             num_freqs, num_time_steps, num_components,
             d_gxs, d_Dxs, d_Dys, d_gys);

  free(metre_delays);
  ( gpuFree(d_azs) );
  ( gpuFree(d_zas) );
  ( gpuFree(d_metre_delays) );
  ( gpuFree(d_beam_has) );
  ( gpuFree(d_beam_decs) );

}


// we have 4 dimensions to loop over here:
//  - component
//  - freqs
//  - times
//  - tiles
//
// it's likely that time will usually be the smallest, so I'm going to
// loop over that. Time will tell whether that's a good idea or not LOL
__global__ void kern_map_hyperbeam_gains(int num_components,
           int num_times, int num_freqs, int num_tiles, int iTime, int num_unique_fee_freqs,
           double *d_jones, const int *d_tile_map, const int *d_freq_map,
           int parallactic,
           gpuUserComplex *d_gxs,
           gpuUserComplex *d_Dxs,
           gpuUserComplex *d_Dys,
           gpuUserComplex *d_gys) {

  //All baselines at all freqs and all times
  int iComponent = threadIdx.x + (blockDim.x*blockIdx.x);

  //The tile to map - currently unused but hopefully in the future tis possible
  int iTile = threadIdx.y + (blockDim.y*blockIdx.y);

  //freq step
  int iFreq = threadIdx.z + (blockDim.z*blockIdx.z);

  if(iComponent < num_components && iTime < num_times && iFreq < num_freqs && iTile < num_tiles) {

    gpuUserComplex d_beam_J00;
    gpuUserComplex d_beam_J01;
    gpuUserComplex d_beam_J10;
    gpuUserComplex d_beam_J11;

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
      current_ind = hyper_ind + iComponent;
    } else {
      hyper_ind = ((num_components*num_times*num_unique_fee_freqs*i_row) + num_components * i_col);

      //There is no frequency used in first chunk of this remapping, as that's handled
      //by the hyper_ind
      current_ind = iTime*num_components*num_tiles + hyper_ind + iComponent;
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

    d_gxs[new_ind] = d_beam_J00;
    d_Dxs[new_ind] = d_beam_J01;
    d_Dys[new_ind] = d_beam_J10;
    d_gys[new_ind] = d_beam_J11;

  }
}


__global__ void fill_with_ones(int num_azza, double *d_jones) {

  //All baselines at all freqs and all times
  int iAzza = threadIdx.x + (blockDim.x*blockIdx.x);

  if (iAzza < num_azza) {
    d_jones[iAzza] = 1;
  }
}


extern "C" void run_hyperbeam_gpu(int num_components,
           int num_time_steps, int num_freqs,
           int num_beams, uint8_t parallactic,
           struct FEEBeamGpu *gpu_fee_beam,
           double *azs, double *zas,
           double *latitudes,
           gpuUserComplex *d_gxs,
           gpuUserComplex *d_Dxs,
           gpuUserComplex *d_Dys,
           gpuUserComplex *d_gys){

  //Always do things in IAU order
  int iau_order = 1;

  int num_azza = num_components * num_time_steps;
  int num_beam_values = num_azza * num_freqs * num_beams;

  //If doing parallactic rotation, call hyperbeam for each time step, meaning
  //we need less memory
  double *d_jones = NULL;

  if (parallactic) {
    ( gpuMalloc( (void**)&(d_jones),
                      2*MAX_POLS*num_beam_values*sizeof(double)) );
  } else {
    ( gpuMalloc( (void**)&(d_jones),
                      2*MAX_POLS*num_beam_values*num_time_steps*sizeof(double)) );
  }

  int32_t status = 0;

  //This is used for the remapping below
  dim3 grid, threads;

  threads.x = 32;
  threads.z = 2;

  grid.x = (int)ceil( (float)num_components / (float)threads.x );
  grid.z = (int)ceil( (float)num_freqs / (float)threads.z );

  if (num_beams > 1) {
    threads.y = 2;
    grid.y = (int)ceil( (float)num_beams / (float)threads.y );
  } else {
    grid.y = threads.y = 1;
  }

  // int num_tiles = num_beams;
  int num_unique_fee_freqs = get_num_unique_fee_freqs(gpu_fee_beam);
  // int num_unique_fee_tiles = get_num_unique_fee_tiles(gpu_fee_beam);

  //Get host pointers to the tile and freq maps
  const int *tile_map;
  const int *freq_map;

  tile_map = get_fee_tile_map(gpu_fee_beam);
  freq_map = get_fee_freq_map(gpu_fee_beam);

  //Copy the tile and freq maps to the GPU
  int32_t *d_tile_map = NULL;
  ( gpuMalloc( (void**)&(d_tile_map),
                      num_beams*sizeof(int32_t) ) );
  int32_t *d_freq_map = NULL;
  ( gpuMalloc( (void**)&(d_freq_map),
                      num_freqs*sizeof(int32_t) ) );

  ( gpuMemcpy(d_tile_map, tile_map,
            num_beams*sizeof(int32_t), gpuMemcpyHostToDevice ) );
  ( gpuMemcpy(d_freq_map, freq_map,
            num_freqs*sizeof(int32_t), gpuMemcpyHostToDevice ) );

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

      status = fee_calc_jones_gpu_device(gpu_fee_beam,
                      (uint32_t)num_components,
                      azs + increment, zas + increment,
                      &latitudes[time_ind],
                      iau_order,
                      (double *)d_jones);

      gpuErrorCheckKernel("kern_map_hyperbeam_gains",
                            kern_map_hyperbeam_gains, grid, threads,
                            num_components, num_time_steps, num_freqs, num_beams,
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
    status = fee_calc_jones_gpu_device(gpu_fee_beam,
                    (uint32_t)num_azza,
                    azs, zas,
                    NULL,
                    iau_order,
                    (double *)d_jones);

    for (int iTime = 0; iTime < num_time_steps; iTime++) {

      gpuErrorCheckKernel("kern_map_hyperbeam_gains",
                            kern_map_hyperbeam_gains, grid, threads,
                            num_components, num_time_steps, num_freqs, num_beams,
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
    handle_hyperbeam_error(__FILE__, __LINE__, "fee_calc_jones_gpu_device");
    // printf("Something went wrong running fee_calc_jones_gpu_device\n");
  }

  ( gpuFree(d_jones) );

}

/*******************************************************************************
                 Functions below to be used in unit tests
*******************************************************************************/

extern "C" void test_RTS_calculate_MWA_analytic_beam(int num_components,
     int num_time_steps, int num_freqs,
     user_precision_t *azs, user_precision_t *zas, int *delays,
     double latitude, int norm,
     double *beam_has, double *beam_decs, double *freqs,
     user_precision_complex_t *gxs, user_precision_complex_t *Dxs,
     user_precision_complex_t *Dys, user_precision_complex_t *gys) {

  //Allocate space for beam gains
  user_precision_complex_t *d_gxs = NULL;
  user_precision_complex_t *d_Dxs = NULL;
  user_precision_complex_t *d_Dys = NULL;
  user_precision_complex_t *d_gys = NULL;

  ( gpuMalloc( (void**)&d_gxs,
    num_freqs*num_time_steps*num_components*sizeof(user_precision_complex_t)) );
  ( gpuMalloc( (void**)&d_Dxs,
    num_freqs*num_time_steps*num_components*sizeof(user_precision_complex_t)) );
  ( gpuMalloc( (void**)&d_Dys,
    num_freqs*num_time_steps*num_components*sizeof(user_precision_complex_t)) );
  ( gpuMalloc( (void**)&d_gys,
    num_freqs*num_time_steps*num_components*sizeof(user_precision_complex_t)) );

  double *d_freqs = NULL;
  ( gpuMalloc( (void**)&d_freqs, num_freqs*sizeof(double)) );
  ( gpuMemcpy(d_freqs, freqs, num_freqs*sizeof(double),
                      gpuMemcpyHostToDevice) );

  //Run
  calculate_RTS_MWA_analytic_beam(num_components,
       num_time_steps, num_freqs,
       azs, zas, delays, latitude, norm,
       beam_has, beam_decs, d_freqs,
       (gpuUserComplex*)d_gxs, (gpuUserComplex*)d_Dxs,
       (gpuUserComplex*)d_Dys, (gpuUserComplex*)d_gys);

  ( gpuMemcpy(gxs, d_gxs,
       num_freqs*num_time_steps*num_components*sizeof(user_precision_complex_t),
       gpuMemcpyDeviceToHost) );
  ( gpuMemcpy(Dxs, d_Dxs,
       num_freqs*num_time_steps*num_components*sizeof(user_precision_complex_t),
       gpuMemcpyDeviceToHost) );
  ( gpuMemcpy(Dys, d_Dys,
       num_freqs*num_time_steps*num_components*sizeof(user_precision_complex_t),
       gpuMemcpyDeviceToHost) );
  ( gpuMemcpy(gys, d_gys,
       num_freqs*num_time_steps*num_components*sizeof(user_precision_complex_t),
       gpuMemcpyDeviceToHost) );

  ( gpuFree(d_gxs) );
  ( gpuFree(d_Dxs) );
  ( gpuFree(d_Dys) );
  ( gpuFree(d_gys) );

  ( gpuFree(d_freqs) );


}

extern "C" void test_analytic_dipole_beam(int num_components,
     int num_time_steps, int num_freqs,
     user_precision_t *azs, user_precision_t *zas, double *freqs,
     user_precision_complex_t *analy_beam_X,
     user_precision_complex_t *analy_beam_Y) {

  user_precision_complex_t *d_analy_beam_X = NULL;
  ( gpuMalloc( (void**)&d_analy_beam_X,
    num_freqs*num_time_steps*num_components*sizeof(user_precision_complex_t)) );

  user_precision_complex_t *d_analy_beam_Y = NULL;
  ( gpuMalloc( (void**)&d_analy_beam_Y,
     num_freqs*num_time_steps*num_components*sizeof(user_precision_complex_t)) );

  double *d_freqs = NULL;
  ( gpuMalloc( (void**)&d_freqs, num_freqs*sizeof(double)) );
  ( gpuMemcpy(d_freqs, freqs, num_freqs*sizeof(double),
                      gpuMemcpyHostToDevice) );

  calculate_analytic_dipole_beam(num_components,
      num_time_steps, num_freqs,
      azs, zas, d_freqs,
      (gpuUserComplex *)d_analy_beam_X, (gpuUserComplex *)d_analy_beam_Y);

  ( gpuMemcpy(analy_beam_X, d_analy_beam_X,
             num_freqs*num_time_steps*num_components*sizeof(user_precision_complex_t),
             gpuMemcpyDeviceToHost) );
  ( gpuMemcpy(analy_beam_Y, d_analy_beam_Y,
             num_freqs*num_time_steps*num_components*sizeof(user_precision_complex_t),
             gpuMemcpyDeviceToHost) );

  ( gpuFree(d_analy_beam_X) );
  ( gpuFree(d_analy_beam_Y) );
  ( gpuFree(d_freqs) );

}

extern "C" void test_kern_gaussian_beam(double *beam_ls, double *beam_ms,
           double beam_ref_freq, double *freqs,
           user_precision_t fwhm_lm, user_precision_t cos_theta, user_precision_t sin_theta, user_precision_t sin_2theta,
           int num_freqs, int num_time_steps, int num_components,
           user_precision_complex_t *primay_beam_J00, user_precision_complex_t *primay_beam_J11) {

  int num_beam_hadec = num_components * num_time_steps;

  double *d_beam_ls = NULL;
  ( gpuMalloc( (void**)&d_beam_ls, num_beam_hadec*sizeof(double)) );
  ( gpuMemcpy(d_beam_ls, beam_ls,
                           num_beam_hadec*sizeof(double), gpuMemcpyHostToDevice ) );

  double *d_beam_ms = NULL;
  ( gpuMalloc( (void**)&d_beam_ms, num_beam_hadec*sizeof(double)) );
  ( gpuMemcpy(d_beam_ms, beam_ms,
                           num_beam_hadec*sizeof(double), gpuMemcpyHostToDevice ) );

  double *d_freqs = NULL;
  ( gpuMalloc( (void**)&d_freqs, num_freqs*sizeof(double) ) );
  ( gpuMemcpy(d_freqs, freqs,
                           num_freqs*sizeof(double), gpuMemcpyHostToDevice ) );

  user_precision_complex_t *d_g1xs = NULL;
  ( gpuMalloc( (void**)&d_g1xs,
                                   num_freqs*num_beam_hadec*sizeof(user_precision_complex_t)) );

  user_precision_complex_t *d_g1ys = NULL;
  ( gpuMalloc( (void**)&d_g1ys,
                                   num_freqs*num_beam_hadec*sizeof(user_precision_complex_t)) );

  dim3 grid, threads;

  threads.x = 16;
  grid.x = (int)ceil( (float)num_beam_hadec / (float)threads.x );

  threads.y = 16;
  grid.y = (int)ceil( (float)num_freqs / (float)threads.y );

  gpuErrorCheckKernel("kern_gaussian_beam",
                        kern_gaussian_beam, grid, threads,
                        d_beam_ls, d_beam_ms,
                        beam_ref_freq, d_freqs,
                        fwhm_lm, cos_theta, sin_theta, sin_2theta,
                        num_freqs, num_time_steps, num_components,
                        (gpuUserComplex *)d_g1xs,
                        (gpuUserComplex *)d_g1ys);

  ( gpuMemcpy(primay_beam_J00, d_g1xs,
                 num_freqs*num_beam_hadec*sizeof(user_precision_complex_t), gpuMemcpyDeviceToHost) );
  ( gpuMemcpy(primay_beam_J11, d_g1ys,
                 num_freqs*num_beam_hadec*sizeof(user_precision_complex_t), gpuMemcpyDeviceToHost) );

  ( gpuFree(d_beam_ls ) );
  ( gpuFree(d_beam_ms ) );
  ( gpuFree(d_freqs ) );
  ( gpuFree(d_g1xs ) );
  ( gpuFree(d_g1ys ) );

}


extern "C" void test_calculate_gaussian_beam(int num_components, int num_time_steps,
     int num_freqs, user_precision_t ha0, user_precision_t sdec0, user_precision_t cdec0,
     user_precision_t fwhm_lm, user_precision_t cos_theta, user_precision_t sin_theta, user_precision_t sin_2theta,
     double beam_ref_freq, double *freqs,
     double *beam_has, double *beam_decs,
     user_precision_complex_t *primay_beam_J00, user_precision_complex_t *primay_beam_J11) {

  int num_beam_hadec = num_components * num_time_steps;

  double *d_freqs = NULL;
  ( gpuMalloc( (void**)&d_freqs, num_freqs*sizeof(double) ) );
  ( gpuMemcpy(d_freqs, freqs,
                           num_freqs*sizeof(double), gpuMemcpyHostToDevice ) );

  user_precision_complex_t *d_g1xs = NULL;
  ( gpuMalloc( (void**)&d_g1xs,
                                   num_freqs*num_beam_hadec*sizeof(user_precision_complex_t)) );

  user_precision_complex_t *d_g1ys = NULL;
  ( gpuMalloc( (void**)&d_g1ys,
                                   num_freqs*num_beam_hadec*sizeof(user_precision_complex_t)) );


  calculate_gaussian_beam(num_components, num_time_steps,
                         num_freqs, ha0, sdec0, cdec0,
                         fwhm_lm, cos_theta, sin_theta, sin_2theta,
                         beam_ref_freq, d_freqs,
                         beam_has, beam_decs,
                         (gpuUserComplex *)d_g1xs,
                         (gpuUserComplex *)d_g1ys);

  ( gpuMemcpy(primay_beam_J00, d_g1xs,
                 num_freqs*num_beam_hadec*sizeof(user_precision_complex_t), gpuMemcpyDeviceToHost) );
  ( gpuMemcpy(primay_beam_J11, d_g1ys,
                 num_freqs*num_beam_hadec*sizeof(user_precision_complex_t), gpuMemcpyDeviceToHost) );

  ( gpuFree(d_freqs ) );
  ( gpuFree(d_g1xs ) );
  ( gpuFree(d_g1ys ) );

}

extern "C" void test_run_hyperbeam_gpu(int num_components,
           int num_time_steps, int num_freqs, int num_ants,
           uint8_t parallatic,
           struct FEEBeamGpu *gpu_fee_beam,
           double *azs, double *zas,
           double *latitudes,
           user_precision_complex_t *primay_beam_J00,
           user_precision_complex_t *primay_beam_J01,
           user_precision_complex_t *primay_beam_J10,
           user_precision_complex_t *primay_beam_J11){
  uint32_t num_azza = num_components * num_time_steps;
  int num_beam_values = num_azza * num_freqs * num_ants;

  user_precision_complex_t *d_gxs = NULL;
  user_precision_complex_t *d_Dxs = NULL;
  user_precision_complex_t *d_Dys = NULL;
  user_precision_complex_t *d_gys = NULL;

  ( gpuMalloc( (void**)&d_gxs,
                      num_beam_values*sizeof(user_precision_complex_t)) );
  ( gpuMalloc( (void**)&d_Dxs,
                      num_beam_values*sizeof(user_precision_complex_t)) );
  ( gpuMalloc( (void**)&d_Dys,
                      num_beam_values*sizeof(user_precision_complex_t)) );
  ( gpuMalloc( (void**)&d_gys,
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

  run_hyperbeam_gpu(num_components,
             num_time_steps, num_freqs, num_ants,
             parallatic,
             gpu_fee_beam,
             reordered_azs, reordered_zas,
             latitudes,
             (gpuUserComplex *)d_gxs,
             (gpuUserComplex *)d_Dxs,
             (gpuUserComplex *)d_Dys,
             (gpuUserComplex *)d_gys);

  ( gpuMemcpy( primay_beam_J00, d_gxs,
                      num_beam_values*sizeof(user_precision_complex_t),
                      gpuMemcpyDeviceToHost) );

  ( gpuMemcpy( primay_beam_J01, d_Dxs,
                      num_beam_values*sizeof(user_precision_complex_t),
                      gpuMemcpyDeviceToHost) );

  ( gpuMemcpy( primay_beam_J10, d_Dys,
                      num_beam_values*sizeof(user_precision_complex_t),
                      gpuMemcpyDeviceToHost) );

  ( gpuMemcpy( primay_beam_J11, d_gys,
                      num_beam_values*sizeof(user_precision_complex_t),
                      gpuMemcpyDeviceToHost) );

  ( gpuFree(d_gxs) );
  ( gpuFree(d_Dxs) );
  ( gpuFree(d_Dys) );
  ( gpuFree(d_gys) );

  free(reordered_azs);
  free(reordered_zas);

}