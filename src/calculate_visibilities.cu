#include <stdio.h>
#include <stdlib.h>
#include <cuComplex.h>
#include <complex.h>
#include <math.h>
#include "woden.h"
#include "calculate_visibilities.h"
#include "shapelet_basis.h"
#include "cudacomplex.h"
#include "fundamental_coords.h"
#include "constants.h"
#include "source_components.h"
#include "primary_beam.h"
#include "FEE_primary_beam_cuda.h"


extern "C" void calculate_visibilities(float *X_diff_metres, float *Y_diff_metres, float *Z_diff_metres,
                    catsource_t catsource, float *angles_array, beam_settings_t beam_settings,
                    const int num_baselines, const int num_time_steps, const int num_visis,
                    const int num_freqs,
                    visibility_set_t *visibility_set,
                    float *sbf2) {

  // printf("CUDA error 0: %s\n", cudaGetErrorString( cudaGetLastError() ) );

  /*START We should be able to do all this outside of this function and transfer in--------------*/
  //==========================================================================================
  float *d_X_diff = NULL;
  float *d_Y_diff = NULL;
  float *d_Z_diff = NULL;

  cudaMalloc( (void**)&d_X_diff, num_baselines*sizeof(float) );
  cudaMemcpy( d_X_diff, X_diff_metres, num_baselines*sizeof(float), cudaMemcpyHostToDevice );
  cudaMalloc( (void**)&d_Y_diff, num_baselines*sizeof(float) );
  cudaMemcpy( d_Y_diff, Y_diff_metres, num_baselines*sizeof(float), cudaMemcpyHostToDevice );
  cudaMalloc( (void**)&d_Z_diff, num_baselines*sizeof(float) );
  cudaMemcpy( d_Z_diff, Z_diff_metres, num_baselines*sizeof(float), cudaMemcpyHostToDevice );

  float *d_angles_array = NULL;
  cudaMalloc( (void**)&d_angles_array, 3*sizeof(float) );
  cudaMemcpy( d_angles_array, angles_array, 3*sizeof(float), cudaMemcpyHostToDevice );

  float *d_sha0s = NULL;
  float *d_cha0s = NULL;
  float *d_wavelengths = NULL;
  cudaMalloc( (void**)&d_sha0s, num_visis*sizeof(float) );
  cudaMemcpy( d_sha0s, visibility_set->sha0s, num_visis*sizeof(float), cudaMemcpyHostToDevice );
  cudaMalloc( (void**)&d_cha0s, num_visis*sizeof(float) );
  cudaMemcpy( d_cha0s, visibility_set->cha0s, num_visis*sizeof(float), cudaMemcpyHostToDevice );
  cudaMalloc( (void**)&d_wavelengths, num_visis*sizeof(float) );
  cudaMemcpy( d_wavelengths, visibility_set->wavelengths, num_visis*sizeof(float), cudaMemcpyHostToDevice );
  //
  /* END We should be able to do all this outside of this function and transfer in--------------*/

  // float *d_sum_visi_real = NULL;
  // float *d_sum_visi_imag = NULL;
  float *d_u_metres = NULL;
  float *d_v_metres = NULL;
  float *d_w_metres = NULL;
  float *d_us = NULL;
  float *d_vs = NULL;
  float *d_ws = NULL;

  cudaMalloc( (void**)&d_u_metres, num_visis*sizeof(float) );
  cudaMalloc( (void**)&d_v_metres, num_visis*sizeof(float) );
  cudaMalloc( (void**)&d_w_metres, num_visis*sizeof(float) );
  cudaMalloc( (void**)&d_us, num_visis*sizeof(float) );
  cudaMalloc( (void**)&d_vs, num_visis*sizeof(float) );
  cudaMalloc( (void**)&d_ws, num_visis*sizeof(float) );

  // printf("CUDA error 1: %s\n", cudaGetErrorString( cudaGetLastError() ) );



  // printf("CUDA error 3: %s\n", cudaGetErrorString( cudaGetLastError() ) );

  float *d_sum_visi_XX_real;
  float *d_sum_visi_XX_imag;
  float *d_sum_visi_XY_real;
  float *d_sum_visi_XY_imag;
  float *d_sum_visi_YX_real;
  float *d_sum_visi_YX_imag;
  float *d_sum_visi_YY_real;
  float *d_sum_visi_YY_imag;

  cudaMalloc( (void**)&d_sum_visi_XX_real, num_visis*sizeof(float) );
  cudaMalloc( (void**)&d_sum_visi_XX_imag, num_visis*sizeof(float) );
  cudaMalloc( (void**)&d_sum_visi_XY_real, num_visis*sizeof(float) );
  cudaMalloc( (void**)&d_sum_visi_XY_imag, num_visis*sizeof(float) );
  cudaMalloc( (void**)&d_sum_visi_YX_real, num_visis*sizeof(float) );
  cudaMalloc( (void**)&d_sum_visi_YX_imag, num_visis*sizeof(float) );
  cudaMalloc( (void**)&d_sum_visi_YY_real, num_visis*sizeof(float) );
  cudaMalloc( (void**)&d_sum_visi_YY_imag, num_visis*sizeof(float) );

  cudaMemcpy(d_sum_visi_XX_real, visibility_set->sum_visi_XX_real, num_visis*sizeof(float), cudaMemcpyHostToDevice );
  cudaMemcpy(d_sum_visi_XX_imag, visibility_set->sum_visi_XX_imag, num_visis*sizeof(float), cudaMemcpyHostToDevice );
  cudaMemcpy(d_sum_visi_XY_real, visibility_set->sum_visi_XY_real, num_visis*sizeof(float), cudaMemcpyHostToDevice );
  cudaMemcpy(d_sum_visi_XY_imag, visibility_set->sum_visi_XY_imag, num_visis*sizeof(float), cudaMemcpyHostToDevice );
  cudaMemcpy(d_sum_visi_YX_real, visibility_set->sum_visi_YX_real, num_visis*sizeof(float), cudaMemcpyHostToDevice );
  cudaMemcpy(d_sum_visi_YX_imag, visibility_set->sum_visi_YX_imag, num_visis*sizeof(float), cudaMemcpyHostToDevice );
  cudaMemcpy(d_sum_visi_YY_real, visibility_set->sum_visi_YY_real, num_visis*sizeof(float), cudaMemcpyHostToDevice );
  cudaMemcpy(d_sum_visi_YY_imag, visibility_set->sum_visi_YY_imag, num_visis*sizeof(float), cudaMemcpyHostToDevice );

  dim3 grid, threads;

  threads.x = 128;
  threads.y = 1;
  grid.x = (int)ceil( (float)num_visis / (float)threads.x );
  grid.y = 1;

  kern_calc_uvw<<< grid, threads >>>( d_X_diff,
          d_Y_diff, d_Z_diff,
          d_u_metres, d_v_metres, d_w_metres,
          d_us, d_vs, d_ws, d_wavelengths,
          d_angles_array, d_cha0s, d_sha0s,
          num_visis, num_baselines);

  int num_points = catsource.n_points;
  int num_gauss = catsource.n_gauss;
  int num_shapes = catsource.n_shapes;

  //Setup some beam related arrays
  float *d_beam_angles_array = NULL;
  cudaMalloc( (void**)&d_beam_angles_array, 3*sizeof(float) );
  cudaMemcpy( d_beam_angles_array, beam_settings.beam_angles_array, 3*sizeof(float), cudaMemcpyHostToDevice );

  float *d_freqs = NULL;
  cudaMalloc( (void**)&d_freqs, num_freqs*sizeof(float) );
  cudaMemcpy( d_freqs, visibility_set->channel_frequencies, num_freqs*sizeof(float), cudaMemcpyHostToDevice );

  //TODO currently hardcoded to have beam position angle = 0. Should this change with az/za?
  float cos_theta = 1.0;
  float sin_theta = 0.0;
  float sin_2theta = 0.0;
  float fwhm_lm; //= 20.0 * D2R;

  float *d_beam_ref_freq = NULL;
  if (beam_settings.beamtype == GAUSS_BEAM) {

    cudaMalloc( (void**)&d_beam_ref_freq, sizeof(float) );
    cudaMemcpy( d_beam_ref_freq, beam_settings.beam_ref_freq_array, sizeof(float), cudaMemcpyHostToDevice );

    fwhm_lm = sinf(beam_settings.beam_FWHM_rad);

  }

  float *d_gauss_beam_reals = NULL;
  float *d_gauss_beam_imags = NULL;

  cuFloatComplex *d_analy_beam_X = NULL;
  cuFloatComplex *d_analy_beam_Y = NULL;

  if (num_points > 0) {
    printf("\tDoing point components\n");

    float *d_point_ras=NULL;
    float *d_point_decs=NULL;
    float *d_point_freqs=NULL;
    float *d_point_stokesI=NULL;
    float *d_point_stokesQ=NULL;
    float *d_point_stokesU=NULL;
    float *d_point_stokesV=NULL;
    float *d_point_SIs=NULL;

    cudaMalloc( (void**)&(d_point_ras), num_points*sizeof(float) );
    cudaMemcpy( d_point_ras, catsource.point_ras, num_points*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_point_decs), num_points*sizeof(float) );
    cudaMemcpy( d_point_decs, catsource.point_decs, num_points*sizeof(float), cudaMemcpyHostToDevice );

    // cudaMalloc( (void**)&(d_point_fluxes), num_points*sizeof(float) );
    // cudaMemcpy( d_point_fluxes, catsource.point_fluxes, num_points*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_point_freqs), num_points*sizeof(float) );
    cudaMemcpy( d_point_freqs, catsource.point_ref_freqs, num_points*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_point_stokesI), num_points*sizeof(float) );
    cudaMemcpy( d_point_stokesI, catsource.point_ref_stokesI, num_points*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_point_stokesQ), num_points*sizeof(float) );
    cudaMemcpy( d_point_stokesQ, catsource.point_ref_stokesQ, num_points*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_point_stokesU), num_points*sizeof(float) );
    cudaMemcpy( d_point_stokesU, catsource.point_ref_stokesU, num_points*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_point_stokesV), num_points*sizeof(float) );
    cudaMemcpy( d_point_stokesV, catsource.point_ref_stokesV, num_points*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_point_SIs), num_points*sizeof(float) );
    cudaMemcpy( d_point_SIs, catsource.point_SIs, num_points*sizeof(float), cudaMemcpyHostToDevice );

    float *d_ls=NULL;
    float *d_ms=NULL;
    float *d_ns=NULL;

    cudaMalloc( (void**)&d_ls, num_points*sizeof(float) );
    cudaMalloc( (void**)&d_ms, num_points*sizeof(float) );
    cudaMalloc( (void**)&d_ns, num_points*sizeof(float) );

    threads.x = 128;
    threads.y = 1;
    threads.z = 1;
    grid.x = (int)ceil( (float)num_points / (float)threads.x );
    grid.y = 1;
    grid.z = 1;

    kern_calc_lmn<<< grid, threads >>>(d_angles_array, d_point_ras, d_point_decs,
                               d_ls, d_ms, d_ns, num_points);

    //If using a gaussian primary beam, calculate beam values for all freqs,
    //lsts and point component locations
    if (beam_settings.beamtype == GAUSS_BEAM) {
      cudaMalloc( (void**)&d_gauss_beam_reals, beam_settings.num_point_beam_values*sizeof(float) );
      cudaMalloc( (void**)&d_gauss_beam_imags, beam_settings.num_point_beam_values*sizeof(float) );

      printf("\tDoing gaussian beam tings\n");

      calculate_gaussian_beam(num_points, num_time_steps, num_freqs,
           fwhm_lm, cos_theta, sin_theta, sin_2theta,
           d_beam_ref_freq, d_freqs, d_beam_angles_array,
           beam_settings.beam_point_has, beam_settings.beam_point_decs,
           d_gauss_beam_reals, d_gauss_beam_imags);

    }// end if beam == GAUSS

    if (beam_settings.beamtype == FEE_BEAM) {

      cudaMalloc( (void**)&d_gauss_beam_reals, num_time_steps*num_points*sizeof(float) );
      cudaMalloc( (void**)&d_gauss_beam_imags, num_time_steps*num_points*sizeof(float) );

      calc_CUDA_FEE_beam(d_gauss_beam_reals, d_gauss_beam_imags,
             catsource.point_azs, catsource.point_zas,
             catsource.sin_point_para_angs, catsource.cos_point_para_angs,
             num_points, num_time_steps, beam_settings.FEE_beam);

      // printf("calculate_visibilities error 2: %s\n", cudaGetErrorString( cudaGetLastError() ) );
    }

    if (beam_settings.beamtype == ANALY_DIPOLE) {
      printf("\tTrying to do analytic_dipole\n");

      cudaMalloc( (void**)&d_gauss_beam_reals, sizeof(float) );
      cudaMalloc( (void**)&d_gauss_beam_imags, sizeof(float) );

      cudaMalloc( (void**)&d_analy_beam_X, num_time_steps*num_points*sizeof(cuFloatComplex) );
      cudaMalloc( (void**)&d_analy_beam_Y, num_time_steps*num_points*sizeof(cuFloatComplex) );

      calculate_analytic_dipole_beam(num_points, num_time_steps, num_freqs,
           catsource.point_azs, catsource.point_zas, d_freqs,
           d_analy_beam_X, d_analy_beam_Y);

      // printf("calculate_visibilities error 2: %s\n", cudaGetErrorString( cudaGetLastError() ) );
    }

    if (num_points == 1) {
      threads.x = 128;
      threads.y = 1;
      grid.x = grid.x = (int)ceil( (float)num_visis / (float)threads.x );
      grid.y = 1;
    }
    else {
      threads.x = 64;
      threads.y = 4;
      grid.x = (int)ceil( (float)num_visis / (float)threads.x );
      grid.y = (int)ceil( ((float)num_points) / ((float)threads.y) );
    }

    kern_calc_visi_point<<<grid , threads>>>(d_point_ras, d_point_decs,
            d_point_freqs, d_point_stokesI, d_point_stokesQ,
            d_point_stokesU, d_point_stokesV, d_point_SIs,
            d_us, d_vs, d_ws,
            d_sum_visi_XX_real, d_sum_visi_XX_imag,
            d_sum_visi_XY_real, d_sum_visi_XY_imag,
            d_sum_visi_YX_real, d_sum_visi_YX_imag,
            d_sum_visi_YY_real, d_sum_visi_YY_imag,
            d_angles_array, d_wavelengths,
            d_ls, d_ms, d_ns,
            num_points, num_baselines, num_freqs, num_visis,
            num_time_steps,
            d_gauss_beam_reals, d_gauss_beam_imags,  beam_settings.beamtype,
            (cuFloatComplex *)beam_settings.FEE_beam->d_FEE_beam_gain_matrices,
            d_analy_beam_X, d_analy_beam_Y);

    cudaFree( beam_settings.FEE_beam->d_FEE_beam_gain_matrices);
    cudaFree( d_gauss_beam_imags);
    cudaFree( d_gauss_beam_reals);
    cudaFree( d_ns);
    cudaFree( d_ms);
    cudaFree( d_ls);
    cudaFree( d_point_freqs );
    // cudaFree( d_point_fluxes );
    cudaFree( d_point_stokesI );
    cudaFree( d_point_stokesQ );
    cudaFree( d_point_stokesU );
    cudaFree( d_point_stokesV );
    cudaFree( d_point_SIs );
    cudaFree( d_point_decs);
    cudaFree( d_point_ras);

    cudaFree( d_analy_beam_X );
    cudaFree( d_analy_beam_Y );

  }//if point sources

  if (num_gauss > 0) {
    printf("\tDoing gaussian components\n");

    float *d_gauss_ras=NULL;
    float *d_gauss_decs=NULL;
    float *d_gauss_pas=NULL;
    float *d_gauss_majors=NULL;
    float *d_gauss_minors=NULL;

    float *d_gauss_freqs=NULL;
    float *d_gauss_stokesI=NULL;
    float *d_gauss_stokesQ=NULL;
    float *d_gauss_stokesU=NULL;
    float *d_gauss_stokesV=NULL;
    float *d_gauss_SIs=NULL;

    cudaMalloc( (void**)&(d_gauss_freqs), num_gauss*sizeof(float) );
    cudaMemcpy( d_gauss_freqs, catsource.gauss_ref_freqs, num_gauss*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_gauss_stokesI), num_gauss*sizeof(float) );
    cudaMemcpy( d_gauss_stokesI, catsource.gauss_ref_stokesI, num_gauss*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_gauss_stokesQ), num_gauss*sizeof(float) );
    cudaMemcpy( d_gauss_stokesQ, catsource.gauss_ref_stokesQ, num_gauss*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_gauss_stokesU), num_gauss*sizeof(float) );
    cudaMemcpy( d_gauss_stokesU, catsource.gauss_ref_stokesU, num_gauss*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_gauss_stokesV), num_gauss*sizeof(float) );
    cudaMemcpy( d_gauss_stokesV, catsource.gauss_ref_stokesV, num_gauss*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_gauss_SIs), num_gauss*sizeof(float) );
    cudaMemcpy( d_gauss_SIs, catsource.gauss_SIs, num_gauss*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_gauss_ras), num_gauss*sizeof(float) );
    cudaMemcpy( d_gauss_ras, catsource.gauss_ras, num_gauss*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_gauss_decs), num_gauss*sizeof(float) );
    cudaMemcpy( d_gauss_decs, catsource.gauss_decs, num_gauss*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_gauss_pas), num_gauss*sizeof(float) );
    cudaMemcpy( d_gauss_pas, catsource.gauss_pas, num_gauss*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_gauss_majors), num_gauss*sizeof(float) );
    cudaMemcpy( d_gauss_majors, catsource.gauss_majors, num_gauss*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_gauss_minors), num_gauss*sizeof(float) );
    cudaMemcpy( d_gauss_minors, catsource.gauss_minors, num_gauss*sizeof(float), cudaMemcpyHostToDevice );

    float *d_ls=NULL;
    float *d_ms=NULL;
    float *d_ns=NULL;

    cudaMalloc( (void**)&d_ls, num_gauss*sizeof(float) );
    cudaMalloc( (void**)&d_ms, num_gauss*sizeof(float) );
    cudaMalloc( (void**)&d_ns, num_gauss*sizeof(float) );

    threads.x = 128;
    threads.y = 1;
    threads.z = 1;
    grid.x = (int)ceil( (float)num_gauss / (float)threads.x );
    grid.y = 1;
    grid.z = 1;

    kern_calc_lmn<<< grid , threads >>>(d_angles_array, d_gauss_ras, d_gauss_decs,
                               d_ls, d_ms, d_ns, num_gauss);


    //TODO make these 2 by 2 (or have 4 arrays) to get instrumental pol going when
    //we implement the FEE beam
    //Make some empty beam stuff
    // float *d_gauss_beam_reals = NULL;
    // float *d_gauss_beam_imags = NULL;

    if (beam_settings.beamtype == GAUSS_BEAM) {

      cudaMalloc( (void**)&d_gauss_beam_imags, beam_settings.num_gausscomp_beam_values*sizeof(float) );
      cudaMalloc( (void**)&d_gauss_beam_reals, beam_settings.num_gausscomp_beam_values*sizeof(float) );

      calculate_gaussian_beam(num_gauss, num_time_steps, num_freqs,
           fwhm_lm, cos_theta, sin_theta, sin_2theta,
           d_beam_ref_freq, d_freqs, d_beam_angles_array,
           beam_settings.beam_gausscomp_has, beam_settings.beam_gausscomp_decs,
           d_gauss_beam_reals, d_gauss_beam_imags);

    }// end if beam == GAUSS

    if (beam_settings.beamtype == FEE_BEAM) {

      cudaMalloc( (void**)&d_gauss_beam_reals, num_time_steps*num_gauss*sizeof(float) );
      cudaMalloc( (void**)&d_gauss_beam_imags, num_time_steps*num_gauss*sizeof(float) );

      calc_CUDA_FEE_beam(d_gauss_beam_reals, d_gauss_beam_imags,
             catsource.gauss_azs, catsource.gauss_zas,
             catsource.sin_gauss_para_angs, catsource.cos_gauss_para_angs,
             num_gauss, num_time_steps, beam_settings.FEE_beam);
    }

    if (num_gauss == 1) {
     threads.x = 128;
     threads.y = 1;
     grid.x = grid.x = (int)ceil( (float)num_visis / (float)threads.x );
     grid.y = 1;
    }
    else {
     threads.x = 64;
     threads.y = 4;
     grid.x = (int)ceil( (float)num_visis / (float)threads.x );
     grid.y = (int)ceil( ((float)num_gauss) / ((float)threads.y) );
    }

    kern_calc_visi_gaussian<<<grid , threads>>>(d_gauss_ras, d_gauss_decs,
            d_gauss_freqs, d_gauss_stokesI, d_gauss_stokesQ,
            d_gauss_stokesU, d_gauss_stokesV, d_gauss_SIs,
            d_us, d_vs, d_ws,
            d_sum_visi_XX_real, d_sum_visi_XX_imag,
            d_sum_visi_XY_real, d_sum_visi_XY_imag,
            d_sum_visi_YX_real, d_sum_visi_YX_imag,
            d_sum_visi_YY_real, d_sum_visi_YY_imag,
            d_angles_array, d_wavelengths,
            d_ls, d_ms, d_ns,
            d_gauss_pas, d_gauss_majors, d_gauss_minors,
            num_gauss, num_baselines, num_freqs, num_visis, num_time_steps,
            d_gauss_beam_reals, d_gauss_beam_imags,  beam_settings.beamtype,
            (cuFloatComplex *)beam_settings.FEE_beam->d_FEE_beam_gain_matrices);

    cudaFree( beam_settings.FEE_beam->d_FEE_beam_gain_matrices);
    cudaFree( d_gauss_beam_imags);
    cudaFree( d_gauss_beam_reals);
    cudaFree( d_ns);
    cudaFree( d_ms);
    cudaFree( d_ls);
    cudaFree( d_gauss_minors );
    cudaFree( d_gauss_majors);
    cudaFree( d_gauss_pas);
    cudaFree( d_gauss_decs);
    cudaFree( d_gauss_ras);
    cudaFree( d_gauss_freqs );
    cudaFree( d_gauss_stokesI );
    cudaFree( d_gauss_stokesQ );
    cudaFree( d_gauss_stokesU );
    cudaFree( d_gauss_stokesV );
    cudaFree( d_gauss_SIs );

  }//if gauss sources

  if (num_shapes > 0) {
    printf("\tDoing shapelet components\n");

    float *d_shape_ras=NULL;
    float *d_shape_decs=NULL;
    float *d_shape_fluxes=NULL;
    float *d_shape_freqs=NULL;
    float *d_shape_pas=NULL;
    float *d_shape_majors=NULL;
    float *d_shape_minors=NULL;

    float *d_shape_coeffs=NULL;
    float *d_shape_n1s=NULL;
    float *d_shape_n2s=NULL;
    float *d_shape_param_indexes=NULL;

    float *d_sbf=NULL;
    float *d_lsts=NULL;

    float *d_u_s_metres = NULL;
    float *d_v_s_metres = NULL;
    float *d_w_s_metres = NULL;

    float *d_shape_ls=NULL;
    float *d_shape_ms=NULL;
    float *d_shape_ns=NULL;

    //Who likes cudaMalloc cudaMallocs? We like cudaMalloc cudaMallocs
    cudaMalloc( (void**)&(d_shape_ras), num_shapes*sizeof(float) );
    cudaMemcpy( d_shape_ras, catsource.shape_ras, num_shapes*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_shape_decs), num_shapes*sizeof(float) );
    cudaMemcpy( d_shape_decs, catsource.shape_decs, num_shapes*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_shape_fluxes), num_shapes*sizeof(float) );
    cudaMemcpy( d_shape_fluxes, catsource.shape_fluxes, num_shapes*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_shape_pas), num_shapes*sizeof(float) );
    cudaMemcpy( d_shape_pas, catsource.shape_pas, num_shapes*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_shape_majors), num_shapes*sizeof(float) );
    cudaMemcpy( d_shape_majors, catsource.shape_majors, num_shapes*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_shape_minors), num_shapes*sizeof(float) );
    cudaMemcpy( d_shape_minors, catsource.shape_minors, num_shapes*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_shape_freqs), num_shapes*sizeof(float) );
    cudaMemcpy( d_shape_freqs, catsource.shape_freqs, num_shapes*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_shape_coeffs), catsource.n_shape_coeffs*sizeof(float) );
    cudaMemcpy( d_shape_coeffs, catsource.shape_coeffs, catsource.n_shape_coeffs*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_shape_n1s), catsource.n_shape_coeffs*sizeof(float) );
    cudaMemcpy( d_shape_n1s, catsource.shape_n1s, catsource.n_shape_coeffs*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_shape_n2s), catsource.n_shape_coeffs*sizeof(float) );
    cudaMemcpy( d_shape_n2s, catsource.shape_n2s, catsource.n_shape_coeffs*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_shape_param_indexes), catsource.n_shape_coeffs*sizeof(float) );
    cudaMemcpy( d_shape_param_indexes, catsource.shape_param_indexes, catsource.n_shape_coeffs*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_sbf), sbf_N*sbf_L*sizeof(float) );
    cudaMemcpy( d_sbf, sbf2, sbf_N*sbf_L*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_shape_ls), num_shapes*sizeof(float) );
    cudaMalloc( (void**)&(d_shape_ms), num_shapes*sizeof(float) );
    cudaMalloc( (void**)&(d_shape_ns), num_shapes*sizeof(float) );

    cudaMalloc( (void**)&(d_lsts), num_visis*sizeof(float) );
    cudaMemcpy( d_lsts, visibility_set->lsts, num_visis*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&d_u_s_metres, num_shapes*num_visis*sizeof(float) );
    cudaMalloc( (void**)&d_v_s_metres, num_shapes*num_visis*sizeof(float) );
    cudaMalloc( (void**)&d_w_s_metres, num_shapes*num_visis*sizeof(float) );


    threads.x = 128;
    threads.y = 1;
    grid.x = (int)ceil( ((float)num_shapes / (float)threads.x) );
    grid.y = 1;

    kern_calc_lmn<<< grid , threads, 0 >>>(d_angles_array, d_shape_ras, d_shape_decs,
                             d_shape_ls, d_shape_ms, d_shape_ns, num_shapes);


    if (num_shapes == 1) {
      threads.x = 128;
      threads.y = 1;
      grid.x = (int)ceil( (float)num_visis / (float)threads.x );
      grid.y = 1;
    }
    else {
      threads.x = 64;
      threads.y = 2;
      grid.x = (int)ceil( (float)num_visis / (float)threads.x );
      grid.y = (int)ceil( ((float)num_shapes) / ((float)threads.y) );
    }

    kern_calc_uvw_shapelet<<< grid, threads >>>(d_X_diff, d_Y_diff, d_Z_diff,
          d_u_s_metres, d_v_s_metres, d_w_s_metres,
          d_lsts, d_shape_ras, d_shape_decs,
          num_baselines, num_visis, num_shapes);

    //TODO make these 2 by 2 (or have 4 arrays) to get instrumental pol going when
    //we implement the FEE beam
    //Make some empty beam stuff
    // float *d_gauss_beam_reals = NULL;
    // float *d_gauss_beam_imags = NULL;

    if (beam_settings.beamtype == GAUSS_BEAM) {
      cudaMalloc( (void**)&d_gauss_beam_reals, beam_settings.num_shape_beam_values*sizeof(float) );
      cudaMalloc( (void**)&d_gauss_beam_imags, beam_settings.num_shape_beam_values*sizeof(float) );

      calculate_gaussian_beam(num_shapes, num_time_steps, num_freqs,
           fwhm_lm, cos_theta, sin_theta, sin_2theta,
           d_beam_ref_freq, d_freqs, d_beam_angles_array,
           beam_settings.beam_shape_has, beam_settings.beam_shape_decs,
           d_gauss_beam_reals, d_gauss_beam_imags);

    }// end if beam == GAUSS

    if (beam_settings.beamtype == FEE_BEAM) {

      cudaMalloc( (void**)&d_gauss_beam_reals, num_time_steps*num_shapes*sizeof(float) );
      cudaMalloc( (void**)&d_gauss_beam_imags, num_time_steps*num_shapes*sizeof(float) );

      calc_CUDA_FEE_beam(d_gauss_beam_reals, d_gauss_beam_imags,
                         catsource.shape_azs, catsource.shape_zas,
                         catsource.sin_shape_para_angs, catsource.cos_shape_para_angs,
                         num_shapes, num_time_steps, beam_settings.FEE_beam);
    }

    if (catsource.n_shape_coeffs == 1) {
      threads.x = 64;
      threads.y = 1;
      grid.x = (int)ceil( (float)num_visis / (float)threads.x );
      grid.y = 1;
    }
    else {
      threads.x = 64;
      threads.y = 2;
      grid.x = (int)ceil( (float)num_visis / (float)threads.x );
      grid.y = (int)ceil( ((float)catsource.n_shape_coeffs) / ((float)threads.y) );
    }

    kern_calc_visi_shapelets<<< grid, threads >>>(d_shape_ras,
            d_shape_decs, d_shape_fluxes, d_shape_freqs,
            d_us, d_vs, d_ws, d_wavelengths,
            d_u_s_metres, d_v_s_metres, d_w_s_metres,
            d_sum_visi_XX_real, d_sum_visi_XX_imag,
            d_sum_visi_XY_real, d_sum_visi_XY_imag,
            d_sum_visi_YX_real, d_sum_visi_YX_imag,
            d_sum_visi_YY_real, d_sum_visi_YY_imag,
            d_angles_array, d_shape_pas, d_shape_majors, d_shape_minors,
            d_shape_n1s, d_shape_n2s, d_shape_coeffs, d_shape_param_indexes,
            d_shape_ls, d_shape_ms, d_shape_ns,
            d_sbf,
            num_shapes, num_baselines, num_freqs, num_visis,
            catsource.n_shape_coeffs, num_time_steps,
            d_gauss_beam_reals, d_gauss_beam_imags, beam_settings.beamtype,
            (cuFloatComplex *)beam_settings.FEE_beam->d_FEE_beam_gain_matrices);

    cudaFree( beam_settings.FEE_beam->d_FEE_beam_gain_matrices);

    cudaFree( d_gauss_beam_imags);
    cudaFree( d_gauss_beam_reals);

    cudaFree( d_shape_ns );
    cudaFree( d_shape_ms );
    cudaFree( d_shape_ls );
    cudaFree(d_w_s_metres);
    cudaFree(d_v_s_metres);
    cudaFree(d_u_s_metres);
    cudaFree(d_lsts);
    cudaFree(d_sbf);
    cudaFree(d_shape_param_indexes);
    cudaFree(d_shape_n2s);
    cudaFree(d_shape_n1s);
    cudaFree(d_shape_coeffs);
    cudaFree( d_shape_minors );
    cudaFree( d_shape_majors);
    cudaFree( d_shape_pas);
    cudaFree( d_shape_freqs );
    cudaFree( d_shape_fluxes );
    cudaFree( d_shape_decs);
    cudaFree( d_shape_ras);

  }//if shapelet

  //Get the results into host memory
  // cudaMemcpy(visibility_set->sum_visi_real,d_sum_visi_real,num_visis*sizeof(float),cudaMemcpyDeviceToHost);
  // cudaMemcpy(visibility_set->sum_visi_imag,d_sum_visi_imag,num_visis*sizeof(float),cudaMemcpyDeviceToHost);

  cudaMemcpy(visibility_set->sum_visi_XX_real,d_sum_visi_XX_real,num_visis*sizeof(float),cudaMemcpyDeviceToHost);
  cudaMemcpy(visibility_set->sum_visi_XX_imag,d_sum_visi_XX_imag,num_visis*sizeof(float),cudaMemcpyDeviceToHost);

  cudaMemcpy(visibility_set->sum_visi_XY_real,d_sum_visi_XY_real,num_visis*sizeof(float),cudaMemcpyDeviceToHost);
  cudaMemcpy(visibility_set->sum_visi_XY_imag,d_sum_visi_XY_imag,num_visis*sizeof(float),cudaMemcpyDeviceToHost);

  cudaMemcpy(visibility_set->sum_visi_YX_real,d_sum_visi_YX_real,num_visis*sizeof(float),cudaMemcpyDeviceToHost);
  cudaMemcpy(visibility_set->sum_visi_YX_imag,d_sum_visi_YX_imag,num_visis*sizeof(float),cudaMemcpyDeviceToHost);

  cudaMemcpy(visibility_set->sum_visi_YY_real,d_sum_visi_YY_real,num_visis*sizeof(float),cudaMemcpyDeviceToHost);
  cudaMemcpy(visibility_set->sum_visi_YY_imag,d_sum_visi_YY_imag,num_visis*sizeof(float),cudaMemcpyDeviceToHost);

  cudaMemcpy(visibility_set->us_metres,d_u_metres,num_visis*sizeof(float),cudaMemcpyDeviceToHost);
  cudaMemcpy(visibility_set->vs_metres,d_v_metres,num_visis*sizeof(float),cudaMemcpyDeviceToHost);
  cudaMemcpy(visibility_set->ws_metres,d_w_metres,num_visis*sizeof(float),cudaMemcpyDeviceToHost);

  //Free up the GPU memory

  cudaFree( d_freqs );
  cudaFree( d_beam_angles_array );
  cudaFree( d_beam_ref_freq );

  cudaFree( d_ws );
  cudaFree( d_vs );
  cudaFree( d_us );
  cudaFree( d_w_metres );
  cudaFree( d_v_metres );
  cudaFree( d_u_metres );
  // cudaFree( d_sum_visi_imag );
  // cudaFree( d_sum_visi_real );
  cudaFree( d_wavelengths );
  cudaFree( d_cha0s );
  cudaFree( d_sha0s );
  cudaFree( d_angles_array );
  cudaFree( d_Z_diff );
  cudaFree( d_Y_diff );
  cudaFree( d_X_diff );


  cudaFree( d_sum_visi_XX_imag );
  cudaFree( d_sum_visi_XX_real );
  cudaFree( d_sum_visi_XY_imag );
  cudaFree( d_sum_visi_XY_real );
  cudaFree( d_sum_visi_YX_imag );
  cudaFree( d_sum_visi_YX_real );
  cudaFree( d_sum_visi_YY_imag );
  cudaFree( d_sum_visi_YY_real );

}
