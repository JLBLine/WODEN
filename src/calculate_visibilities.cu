#include <stdio.h>
#include <stdlib.h>
#include <cuComplex.h>
#include <complex.h>
#include <math.h>

#include "woden_precision_defs.h"
#include "cudacomplex.h"

#include "calculate_visibilities.h"
#include "shapelet_basis.h"
#include "fundamental_coords.h"
#include "constants.h"
#include "source_components.h"
#include "primary_beam_cuda.h"
#include "cudacheck.h"
#include "visibility_set.h"

extern "C" void calculate_visibilities(array_layout_t *array_layout,
  source_catalogue_t *cropped_sky_models, beam_settings_t *beam_settings,
  woden_settings_t *woden_settings, visibility_set_t *visibility_set,
  user_precision_t *sbf) {

  const int num_baselines = woden_settings->num_baselines;
  const int num_time_steps = woden_settings->num_time_steps;
  const int num_visis = woden_settings->num_visis;
  const int num_freqs = woden_settings->num_freqs;

  //Setup chunk_visibility_set to hold the visibility outputs of each
  visibility_set_t *chunk_visibility_set = setup_visibility_set(num_visis);

  //Copy across settings for simulation - we always simulate all time steps
  //and frequncies for each chunk
  chunk_visibility_set->allsteps_sha0s = visibility_set->allsteps_sha0s;
  chunk_visibility_set->allsteps_cha0s = visibility_set->allsteps_cha0s;
  chunk_visibility_set->allsteps_lsts = visibility_set->allsteps_lsts;
  chunk_visibility_set->allsteps_wavelengths = visibility_set->allsteps_wavelengths;
  chunk_visibility_set->channel_frequencies = visibility_set->channel_frequencies;

  double *d_X_diff = NULL;
  double *d_Y_diff = NULL;
  double *d_Z_diff = NULL;

  cudaErrorCheckCall( cudaMalloc( (void**)&d_X_diff, num_baselines*sizeof(double) ) );
  cudaErrorCheckCall( cudaMemcpy( d_X_diff, array_layout->X_diff_metres,
                      num_baselines*sizeof(double), cudaMemcpyHostToDevice ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_Y_diff, num_baselines*sizeof(double) ) );
  cudaErrorCheckCall( cudaMemcpy( d_Y_diff, array_layout->Y_diff_metres,
                      num_baselines*sizeof(double), cudaMemcpyHostToDevice ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_Z_diff, num_baselines*sizeof(double) ) );
  cudaErrorCheckCall( cudaMemcpy( d_Z_diff, array_layout->Z_diff_metres,
                      num_baselines*sizeof(double), cudaMemcpyHostToDevice ) );

  double *d_allsteps_sha0s = NULL;
  double *d_allsteps_cha0s = NULL;
  user_precision_t *d_allsteps_wavelengths = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_allsteps_sha0s, num_visis*sizeof(double) ) );
  cudaErrorCheckCall( cudaMemcpy( d_allsteps_sha0s, visibility_set->allsteps_sha0s,
                      num_visis*sizeof(double), cudaMemcpyHostToDevice ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_allsteps_cha0s, num_visis*sizeof(double) ) );
  cudaErrorCheckCall( cudaMemcpy( d_allsteps_cha0s, visibility_set->allsteps_cha0s,
                      num_visis*sizeof(double), cudaMemcpyHostToDevice ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_allsteps_wavelengths,
                                         num_visis*sizeof(user_precision_t) ) );
  cudaErrorCheckCall( cudaMemcpy( d_allsteps_wavelengths, visibility_set->allsteps_wavelengths,
                      num_visis*sizeof(user_precision_t), cudaMemcpyHostToDevice ) );

  user_precision_t *d_u_metres = NULL;
  user_precision_t *d_v_metres = NULL;
  user_precision_t *d_w_metres = NULL;
  user_precision_t *d_us = NULL;
  user_precision_t *d_vs = NULL;
  user_precision_t *d_ws = NULL;

  cudaErrorCheckCall( cudaMalloc( (void**)&d_u_metres, num_visis*sizeof(user_precision_t) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_v_metres, num_visis*sizeof(user_precision_t) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_w_metres, num_visis*sizeof(user_precision_t) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_us, num_visis*sizeof(user_precision_t) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_vs, num_visis*sizeof(user_precision_t) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_ws, num_visis*sizeof(user_precision_t) ) );

  user_precision_t *d_sum_visi_XX_real;
  user_precision_t *d_sum_visi_XX_imag;
  user_precision_t *d_sum_visi_XY_real;
  user_precision_t *d_sum_visi_XY_imag;
  user_precision_t *d_sum_visi_YX_real;
  user_precision_t *d_sum_visi_YX_imag;
  user_precision_t *d_sum_visi_YY_real;
  user_precision_t *d_sum_visi_YY_imag;

  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_XX_real,
                      num_visis*sizeof(user_precision_t) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_XX_imag,
                      num_visis*sizeof(user_precision_t) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_XY_real,
                      num_visis*sizeof(user_precision_t) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_XY_imag,
                      num_visis*sizeof(user_precision_t) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_YX_real,
                      num_visis*sizeof(user_precision_t) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_YX_imag,
                      num_visis*sizeof(user_precision_t) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_YY_real,
                      num_visis*sizeof(user_precision_t) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_YY_imag,
                      num_visis*sizeof(user_precision_t) ) );

  double *d_freqs = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_freqs, num_freqs*sizeof(double) ) );
  cudaErrorCheckCall( cudaMemcpy( d_freqs, visibility_set->channel_frequencies,
                      num_freqs*sizeof(double), cudaMemcpyHostToDevice ) );

  //if we have shapelets in our sky model, copy the shapelet basis functions
  //into GPU memory
  user_precision_t *d_sbf=NULL;
  if (cropped_sky_models->num_shapelets > 0) {
    cudaErrorCheckCall( cudaMalloc( (void**)&(d_sbf), sbf_N*sbf_L*sizeof(user_precision_t) ));
    cudaErrorCheckCall( cudaMemcpy( d_sbf, sbf, sbf_N*sbf_L*sizeof(user_precision_t),
                        cudaMemcpyHostToDevice ));
  }

  if (beam_settings->beamtype == FEE_BEAM || beam_settings->beamtype == FEE_BEAM_INTERP){

    //Convert some WODEN precision stuff into hyperbeam precision stuff
    uint32_t *freqs_hz = (uint32_t*)malloc(woden_settings->num_freqs*sizeof(uint32_t));

    for (int freq_ind = 0; freq_ind < woden_settings->num_freqs; freq_ind++) {
        freqs_hz[freq_ind] = (uint32_t)visibility_set->channel_frequencies[freq_ind];
        // freqs_hz[freq_ind] = (uint32_t)beam_settings->base_middle_freq;
    }

    beam_settings->hyper_delays = (uint32_t*)malloc(16*sizeof(uint32_t));

    for (int delay = 0; delay < 16; delay++) {
      beam_settings->hyper_delays[delay] = (uint32_t)woden_settings->FEE_ideal_delays[delay];
    }

    double amps[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    // uint32_t num_freqs_hyper = ();
    uint32_t num_tiles = 1;
    uint32_t num_amps = 16;
    uint8_t norm_to_zenith = 1;

    int32_t status = new_cuda_fee_beam(beam_settings->fee_beam,
                            freqs_hz,
                            beam_settings->hyper_delays,
                            amps,
                            woden_settings->num_freqs,
                            num_tiles,
                            num_amps,
                            norm_to_zenith,
                            &beam_settings->cuda_fee_beam,
                            beam_settings->hyper_error_str);

    if (status != 0) {
      printf("hyperbeam error %d %s\n", status, beam_settings->hyper_error_str );
    }
  }

  //Iterate through all sky model chunks, calculated visibilities are
  //added to chunk_visibility_set, and then summed onto visibility_set
  for (int chunk = 0; chunk < cropped_sky_models->num_sources; chunk++) {

    catsource_t catsource;
    catsource = cropped_sky_models->catsources[chunk];

    //Make sure the temp visis are 0 at the start of each chunk
    for (int visi = 0; visi < num_visis; visi++) {
      chunk_visibility_set->sum_visi_XX_real[visi] = 0;
      chunk_visibility_set->sum_visi_XX_imag[visi] = 0;
      chunk_visibility_set->sum_visi_XY_real[visi] = 0;
      chunk_visibility_set->sum_visi_XY_imag[visi] = 0;
      chunk_visibility_set->sum_visi_YX_real[visi] = 0;
      chunk_visibility_set->sum_visi_YX_imag[visi] = 0;
      chunk_visibility_set->sum_visi_YY_real[visi] = 0;
      chunk_visibility_set->sum_visi_YY_imag[visi] = 0;
    }

    printf("Processing chunk %d\n", chunk);
    printf("\tNumber of components in chunk are: P %d G %d S_coeffs %d\n",
              catsource.n_points,
              catsource.n_gauss,
              catsource.n_shape_coeffs );

    //ensure d_sum_visi_XX_real are set entirely to zero by copying the host
    //array values, which have been set explictly to zero during chunking
    cudaErrorCheckCall( cudaMemcpy(d_sum_visi_XX_real,
               chunk_visibility_set->sum_visi_XX_real,
               num_visis*sizeof(user_precision_t), cudaMemcpyHostToDevice ) );
    cudaErrorCheckCall( cudaMemcpy(d_sum_visi_XX_imag,
               chunk_visibility_set->sum_visi_XX_imag,
               num_visis*sizeof(user_precision_t), cudaMemcpyHostToDevice ) );
    cudaErrorCheckCall( cudaMemcpy(d_sum_visi_XY_real,
               chunk_visibility_set->sum_visi_XY_real,
               num_visis*sizeof(user_precision_t), cudaMemcpyHostToDevice ) );
    cudaErrorCheckCall( cudaMemcpy(d_sum_visi_XY_imag,
               chunk_visibility_set->sum_visi_XY_imag,
               num_visis*sizeof(user_precision_t), cudaMemcpyHostToDevice ) );
    cudaErrorCheckCall( cudaMemcpy(d_sum_visi_YX_real,
               chunk_visibility_set->sum_visi_YX_real,
               num_visis*sizeof(user_precision_t), cudaMemcpyHostToDevice ) );
    cudaErrorCheckCall( cudaMemcpy(d_sum_visi_YX_imag,
               chunk_visibility_set->sum_visi_YX_imag,
               num_visis*sizeof(user_precision_t), cudaMemcpyHostToDevice ) );
    cudaErrorCheckCall( cudaMemcpy(d_sum_visi_YY_real,
               chunk_visibility_set->sum_visi_YY_real,
               num_visis*sizeof(user_precision_t), cudaMemcpyHostToDevice ) );
    cudaErrorCheckCall( cudaMemcpy(d_sum_visi_YY_imag,
               chunk_visibility_set->sum_visi_YY_imag,
               num_visis*sizeof(user_precision_t), cudaMemcpyHostToDevice ) );

    dim3 grid, threads;

    threads.x = 128;
    threads.y = 1;
    grid.x = (int)ceil( (float)num_visis / (float)threads.x );
    grid.y = 1;

    cudaErrorCheckKernel("kern_calc_uvw",
            kern_calc_uvw, grid, threads,
            d_X_diff, d_Y_diff, d_Z_diff,
            d_u_metres, d_v_metres, d_w_metres,
            d_us, d_vs, d_ws, d_allsteps_wavelengths,
            woden_settings->sdec0,
            woden_settings->cdec0,
            d_allsteps_cha0s, d_allsteps_sha0s,
            num_visis, num_baselines);

    int num_points = catsource.n_points;
    int num_gauss = catsource.n_gauss;
    int num_shapes = catsource.n_shapes;

    //All components types needs beam Jones matrices, so declare here
    cuUserComplex *d_primay_beam_J00 = NULL;
    cuUserComplex *d_primay_beam_J01 = NULL;
    cuUserComplex *d_primay_beam_J10 = NULL;
    cuUserComplex *d_primay_beam_J11 = NULL;
    //Same goes for l,m,n coords
    double *d_ls=NULL;
    double *d_ms=NULL;
    double *d_ns=NULL;

    if (num_points > 0) {
      printf("\tDoing point components\n");

      double *d_point_ras=NULL;
      double *d_point_decs=NULL;
      double *d_point_freqs=NULL;
      user_precision_t *d_point_stokesI=NULL;
      user_precision_t *d_point_stokesQ=NULL;
      user_precision_t *d_point_stokesU=NULL;
      user_precision_t *d_point_stokesV=NULL;
      user_precision_t *d_point_SIs=NULL;

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_point_ras),
                          num_points*sizeof(double) ) );
      cudaErrorCheckCall( cudaMemcpy( d_point_ras, catsource.point_ras,
                          num_points*sizeof(double), cudaMemcpyHostToDevice ) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_point_decs),
                          num_points*sizeof(double) ) );
      cudaErrorCheckCall( cudaMemcpy( d_point_decs, catsource.point_decs,
                          num_points*sizeof(double), cudaMemcpyHostToDevice ) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_point_freqs),
                          num_points*sizeof(double) ) );
      cudaErrorCheckCall( cudaMemcpy( d_point_freqs, catsource.point_ref_freqs,
                          num_points*sizeof(double),
                          cudaMemcpyHostToDevice ) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_point_stokesI),
                          num_points*sizeof(user_precision_t) ) );
      cudaErrorCheckCall( cudaMemcpy( d_point_stokesI, catsource.point_ref_stokesI,
                          num_points*sizeof(user_precision_t),
                          cudaMemcpyHostToDevice ) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_point_stokesQ),
                          num_points*sizeof(user_precision_t) ) );
      cudaErrorCheckCall( cudaMemcpy( d_point_stokesQ, catsource.point_ref_stokesQ,
                          num_points*sizeof(user_precision_t),
                          cudaMemcpyHostToDevice ) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_point_stokesU),
                          num_points*sizeof(user_precision_t) ) );
      cudaErrorCheckCall( cudaMemcpy( d_point_stokesU, catsource.point_ref_stokesU,
                          num_points*sizeof(user_precision_t),
                          cudaMemcpyHostToDevice ) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_point_stokesV),
                          num_points*sizeof(user_precision_t) ) );
      cudaErrorCheckCall( cudaMemcpy( d_point_stokesV, catsource.point_ref_stokesV,
                          num_points*sizeof(user_precision_t),
                          cudaMemcpyHostToDevice ) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_point_SIs),
                          num_points*sizeof(user_precision_t) ) );
      cudaErrorCheckCall( cudaMemcpy( d_point_SIs, catsource.point_SIs,
                          num_points*sizeof(user_precision_t),
                          cudaMemcpyHostToDevice ) );

      //Only the FEE beam currently yields cross pol values, so only malloc what
      //we need here
      if (beam_settings->beamtype == FEE_BEAM || beam_settings->beamtype == FEE_BEAM_INTERP || beam_settings->beamtype == MWA_ANALY) {
        cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J01,
                  catsource.num_point_primarybeam_values*sizeof(cuUserComplex) ));
        cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J10,
                  catsource.num_point_primarybeam_values*sizeof(cuUserComplex) ));
      }

      cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J00,
                  catsource.num_point_primarybeam_values*sizeof(cuUserComplex) ));
      cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J11,
                  catsource.num_point_primarybeam_values*sizeof(cuUserComplex) ));
      //
      cudaErrorCheckCall( cudaMalloc( (void**)&d_ls, num_points*sizeof(double) ) );
      cudaErrorCheckCall( cudaMalloc( (void**)&d_ms, num_points*sizeof(double) ) );
      cudaErrorCheckCall( cudaMalloc( (void**)&d_ns, num_points*sizeof(double) ) );

      //This calculates l,m,n for all components, and any beam calculations
      //that are needed - uses the same methods for all component types
      source_component_common(num_points,
                 d_primay_beam_J00, d_primay_beam_J01,
                 d_primay_beam_J10, d_primay_beam_J11,
                 d_freqs, d_ls, d_ms, d_ns,
                 d_point_ras, d_point_decs,
                 catsource.point_azs, catsource.point_zas,
                 catsource.sin_point_para_angs, catsource.cos_point_para_angs,
                 catsource.point_gaussbeam_has, catsource.point_gaussbeam_decs,
                 woden_settings, beam_settings);

      threads.x = 128;
      threads.y = 1;
      grid.x = (int)ceil( (float)num_visis / (float)threads.x );
      grid.y = 1;

      cudaErrorCheckKernel("kern_calc_visi_point",
                           kern_calc_visi_point, grid, threads,
                           d_point_freqs, d_point_stokesI, d_point_stokesQ,
                           d_point_stokesU, d_point_stokesV, d_point_SIs,
                           d_us, d_vs, d_ws,
                           d_sum_visi_XX_real, d_sum_visi_XX_imag,
                           d_sum_visi_XY_real, d_sum_visi_XY_imag,
                           d_sum_visi_YX_real, d_sum_visi_YX_imag,
                           d_sum_visi_YY_real, d_sum_visi_YY_imag,
                           d_allsteps_wavelengths,
                           d_ls, d_ms, d_ns,
                           num_points, num_baselines, num_freqs, num_visis,
                           num_time_steps, beam_settings->beamtype,
                           d_primay_beam_J00, d_primay_beam_J01,
                           d_primay_beam_J10, d_primay_beam_J11);


      cudaErrorCheckCall( cudaFree( d_ns) );
      cudaErrorCheckCall( cudaFree( d_ms) );
      cudaErrorCheckCall( cudaFree( d_ls) );
      cudaErrorCheckCall( cudaFree( d_point_freqs ) );
      cudaErrorCheckCall( cudaFree( d_point_stokesI ) );
      cudaErrorCheckCall( cudaFree( d_point_stokesQ ) );
      cudaErrorCheckCall( cudaFree( d_point_stokesU ) );
      cudaErrorCheckCall( cudaFree( d_point_stokesV ) );
      cudaErrorCheckCall( cudaFree( d_point_SIs ) );
      cudaErrorCheckCall( cudaFree( d_point_decs) );
      cudaErrorCheckCall( cudaFree( d_point_ras) );

      cudaErrorCheckCall( cudaFree( d_primay_beam_J00 ) );
      cudaErrorCheckCall( cudaFree( d_primay_beam_J11 ) );

      if (beam_settings->beamtype == FEE_BEAM || beam_settings->beamtype == FEE_BEAM_INTERP || beam_settings->beamtype == MWA_ANALY) {
        cudaErrorCheckCall( cudaFree( d_primay_beam_J01 ) );
        cudaErrorCheckCall( cudaFree( d_primay_beam_J10 ) );
      }

    }//if point sources

    if (num_gauss > 0) {
      printf("\tDoing gaussian components\n");

      double *d_gauss_ras=NULL;
      double *d_gauss_decs=NULL;
      user_precision_t *d_gauss_pas=NULL;
      user_precision_t *d_gauss_majors=NULL;
      user_precision_t *d_gauss_minors=NULL;

      double *d_gauss_freqs=NULL;
      user_precision_t *d_gauss_stokesI=NULL;
      user_precision_t *d_gauss_stokesQ=NULL;
      user_precision_t *d_gauss_stokesU=NULL;
      user_precision_t *d_gauss_stokesV=NULL;
      user_precision_t *d_gauss_SIs=NULL;

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_gauss_freqs),
                          num_gauss*sizeof(double)) );
      cudaErrorCheckCall( cudaMemcpy( d_gauss_freqs, catsource.gauss_ref_freqs,
                          num_gauss*sizeof(double),
                          cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_gauss_stokesI),
                          num_gauss*sizeof(user_precision_t)) );
      cudaErrorCheckCall( cudaMemcpy( d_gauss_stokesI, catsource.gauss_ref_stokesI,
                          num_gauss*sizeof(user_precision_t),
                          cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_gauss_stokesQ),
                          num_gauss*sizeof(user_precision_t)) );
      cudaErrorCheckCall( cudaMemcpy( d_gauss_stokesQ, catsource.gauss_ref_stokesQ,
                          num_gauss*sizeof(user_precision_t),
                          cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_gauss_stokesU),
                          num_gauss*sizeof(user_precision_t)) );
      cudaErrorCheckCall( cudaMemcpy( d_gauss_stokesU, catsource.gauss_ref_stokesU,
                          num_gauss*sizeof(user_precision_t),
                          cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_gauss_stokesV),
                          num_gauss*sizeof(user_precision_t)) );
      cudaErrorCheckCall( cudaMemcpy( d_gauss_stokesV, catsource.gauss_ref_stokesV,
                          num_gauss*sizeof(user_precision_t),
                          cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_gauss_SIs),
                          num_gauss*sizeof(user_precision_t)) );
      cudaErrorCheckCall( cudaMemcpy( d_gauss_SIs, catsource.gauss_SIs,
                          num_gauss*sizeof(user_precision_t),
                          cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_gauss_ras),
                          num_gauss*sizeof(double)) );
      cudaErrorCheckCall( cudaMemcpy( d_gauss_ras, catsource.gauss_ras,
                          num_gauss*sizeof(double), cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_gauss_decs),
                          num_gauss*sizeof(double)) );
      cudaErrorCheckCall( cudaMemcpy( d_gauss_decs, catsource.gauss_decs,
                          num_gauss*sizeof(double), cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_gauss_pas),
                          num_gauss*sizeof(user_precision_t)) );
      cudaErrorCheckCall( cudaMemcpy( d_gauss_pas, catsource.gauss_pas,
                          num_gauss*sizeof(user_precision_t),
                          cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_gauss_majors),
                          num_gauss*sizeof(user_precision_t)) );
      cudaErrorCheckCall( cudaMemcpy( d_gauss_majors, catsource.gauss_majors,
                          num_gauss*sizeof(user_precision_t),
                          cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_gauss_minors),
                          num_gauss*sizeof(user_precision_t)) );
      cudaErrorCheckCall( cudaMemcpy( d_gauss_minors, catsource.gauss_minors,
                          num_gauss*sizeof(user_precision_t),
                          cudaMemcpyHostToDevice) );

      //Only the FEE beam currently yields cross pol values, so only malloc what
      //we need here
      if (beam_settings->beamtype == FEE_BEAM || beam_settings->beamtype == FEE_BEAM_INTERP || beam_settings->beamtype == MWA_ANALY) {
        cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J01,
              catsource.num_gauss_primarybeam_values*sizeof(cuUserComplex)) );
        cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J10,
              catsource.num_gauss_primarybeam_values*sizeof(cuUserComplex)) );
      }

      cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J00,
            catsource.num_gauss_primarybeam_values*sizeof(cuUserComplex)) );
      cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J11,
            catsource.num_gauss_primarybeam_values*sizeof(cuUserComplex)) );

      cudaErrorCheckCall( cudaMalloc( (void**)&d_ls, num_gauss*sizeof(double)) );
      cudaErrorCheckCall( cudaMalloc( (void**)&d_ms, num_gauss*sizeof(double)) );
      cudaErrorCheckCall( cudaMalloc( (void**)&d_ns, num_gauss*sizeof(double)) );

      source_component_common(num_gauss,
                 d_primay_beam_J00, d_primay_beam_J01,
                 d_primay_beam_J10, d_primay_beam_J11,
                 d_freqs, d_ls, d_ms, d_ns,
                 d_gauss_ras, d_gauss_decs,
                 catsource.gauss_azs, catsource.gauss_zas,
                 catsource.sin_gauss_para_angs, catsource.cos_gauss_para_angs,
                 catsource.gauss_gaussbeam_has,
                 catsource.gauss_gaussbeam_decs,
                 woden_settings, beam_settings);

      threads.x = 128;
      threads.y = 1;
      grid.x = (int)ceil( (float)num_visis / (float)threads.x );
      grid.y = 1;

      cudaErrorCheckKernel("kern_calc_visi_gaussian",
              kern_calc_visi_gaussian, grid, threads,
              d_gauss_freqs, d_gauss_stokesI, d_gauss_stokesQ,
              d_gauss_stokesU, d_gauss_stokesV, d_gauss_SIs,
              d_us, d_vs, d_ws,
              d_sum_visi_XX_real, d_sum_visi_XX_imag,
              d_sum_visi_XY_real, d_sum_visi_XY_imag,
              d_sum_visi_YX_real, d_sum_visi_YX_imag,
              d_sum_visi_YY_real, d_sum_visi_YY_imag,
              d_allsteps_wavelengths,
              d_ls, d_ms, d_ns,
              d_gauss_pas, d_gauss_majors, d_gauss_minors,
              num_gauss, num_baselines, num_freqs, num_visis,
              num_time_steps, beam_settings->beamtype,
              d_primay_beam_J00, d_primay_beam_J01,
              d_primay_beam_J10, d_primay_beam_J11);

      cudaErrorCheckCall( cudaFree(d_primay_beam_J00) );
      cudaErrorCheckCall( cudaFree(d_primay_beam_J11) );

      if (beam_settings->beamtype == FEE_BEAM || beam_settings->beamtype == FEE_BEAM_INTERP || beam_settings->beamtype == MWA_ANALY) {
        cudaErrorCheckCall( cudaFree( d_primay_beam_J01 ) );
        cudaErrorCheckCall( cudaFree( d_primay_beam_J10 ) );
      }

      cudaErrorCheckCall( cudaFree(d_ns) );
      cudaErrorCheckCall( cudaFree(d_ms) );
      cudaErrorCheckCall( cudaFree(d_ls) );
      cudaErrorCheckCall( cudaFree(d_gauss_minors ) );
      cudaErrorCheckCall( cudaFree(d_gauss_majors) );
      cudaErrorCheckCall( cudaFree(d_gauss_pas) );
      cudaErrorCheckCall( cudaFree(d_gauss_decs) );
      cudaErrorCheckCall( cudaFree(d_gauss_ras) );
      cudaErrorCheckCall( cudaFree(d_gauss_freqs ) );
      cudaErrorCheckCall( cudaFree(d_gauss_stokesI ) );
      cudaErrorCheckCall( cudaFree(d_gauss_stokesQ ) );
      cudaErrorCheckCall( cudaFree(d_gauss_stokesU ) );
      cudaErrorCheckCall( cudaFree(d_gauss_stokesV ) );
      cudaErrorCheckCall( cudaFree(d_gauss_SIs ) );

    }//if gauss sources

    if (num_shapes > 0) {
      printf("\tDoing shapelet components\n");

      double *d_shape_ras=NULL;
      double *d_shape_decs=NULL;
      double *d_shape_freqs=NULL;
      user_precision_t *d_shape_stokesI=NULL;
      user_precision_t *d_shape_stokesQ=NULL;
      user_precision_t *d_shape_stokesU=NULL;
      user_precision_t *d_shape_stokesV=NULL;
      user_precision_t *d_shape_SIs=NULL;
      user_precision_t *d_shape_pas=NULL;
      user_precision_t *d_shape_majors=NULL;
      user_precision_t *d_shape_minors=NULL;

      user_precision_t *d_shape_coeffs=NULL;
      user_precision_t *d_shape_n1s=NULL;
      user_precision_t *d_shape_n2s=NULL;
      user_precision_t *d_shape_param_indexes=NULL;

      double *d_allsteps_lsts=NULL;

      user_precision_t *d_u_shapes = NULL;
      user_precision_t *d_v_shapes = NULL;
      user_precision_t *d_w_shapes = NULL;

      //Who likes cudaMalloc cudaMallocs? We like cudaMalloc cudaMallocs
      cudaErrorCheckCall( cudaMalloc( (void**)&(d_shape_ras),
                          num_shapes*sizeof(double)) );
      cudaErrorCheckCall( cudaMemcpy( d_shape_ras, catsource.shape_ras,
              num_shapes*sizeof(double), cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_shape_decs),
                          num_shapes*sizeof(double)) );
      cudaErrorCheckCall( cudaMemcpy( d_shape_decs, catsource.shape_decs,
              num_shapes*sizeof(double), cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_shape_freqs),
                          num_shapes*sizeof(double)) );
      cudaErrorCheckCall( cudaMemcpy( d_shape_freqs, catsource.shape_ref_freqs,
                          num_shapes*sizeof(double),
                          cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_shape_stokesI),
                          num_shapes*sizeof(user_precision_t)) );
      cudaErrorCheckCall( cudaMemcpy( d_shape_stokesI, catsource.shape_ref_stokesI,
                          num_shapes*sizeof(user_precision_t),
                          cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_shape_stokesQ),
                          num_shapes*sizeof(user_precision_t)) );
      cudaErrorCheckCall( cudaMemcpy( d_shape_stokesQ, catsource.shape_ref_stokesQ,
                          num_shapes*sizeof(user_precision_t),
                          cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_shape_stokesU),
                          num_shapes*sizeof(user_precision_t)) );
      cudaErrorCheckCall( cudaMemcpy( d_shape_stokesU, catsource.shape_ref_stokesU,
                          num_shapes*sizeof(user_precision_t),
                          cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_shape_stokesV),
                          num_shapes*sizeof(user_precision_t)) );
      cudaErrorCheckCall( cudaMemcpy( d_shape_stokesV, catsource.shape_ref_stokesV,
                          num_shapes*sizeof(user_precision_t),
                          cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_shape_SIs),
                          num_shapes*sizeof(user_precision_t)) );
      cudaErrorCheckCall( cudaMemcpy( d_shape_SIs, catsource.shape_SIs,
                          num_shapes*sizeof(user_precision_t),
                          cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_shape_pas),
                          num_shapes*sizeof(user_precision_t)) );
      cudaErrorCheckCall( cudaMemcpy( d_shape_pas, catsource.shape_pas,
              num_shapes*sizeof(user_precision_t), cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_shape_majors),
                          num_shapes*sizeof(user_precision_t)) );
      cudaErrorCheckCall( cudaMemcpy( d_shape_majors, catsource.shape_majors,
              num_shapes*sizeof(user_precision_t), cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_shape_minors),
                          num_shapes*sizeof(user_precision_t)) );
      cudaErrorCheckCall( cudaMemcpy( d_shape_minors, catsource.shape_minors,
              num_shapes*sizeof(user_precision_t), cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_shape_coeffs),
                          catsource.n_shape_coeffs*sizeof(user_precision_t)) );
      cudaErrorCheckCall( cudaMemcpy( d_shape_coeffs, catsource.shape_coeffs,
              catsource.n_shape_coeffs*sizeof(user_precision_t),
              cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_shape_n1s),
                          catsource.n_shape_coeffs*sizeof(user_precision_t)) );
      cudaErrorCheckCall( cudaMemcpy( d_shape_n1s, catsource.shape_n1s,
              catsource.n_shape_coeffs*sizeof(user_precision_t),
              cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_shape_n2s),
                          catsource.n_shape_coeffs*sizeof(user_precision_t)) );
      cudaErrorCheckCall( cudaMemcpy( d_shape_n2s, catsource.shape_n2s,
              catsource.n_shape_coeffs*sizeof(user_precision_t),
              cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_shape_param_indexes),
                          catsource.n_shape_coeffs*sizeof(user_precision_t)) );
      cudaErrorCheckCall( cudaMemcpy( d_shape_param_indexes,
             catsource.shape_param_indexes,
             catsource.n_shape_coeffs*sizeof(user_precision_t),
             cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_ls), num_shapes*sizeof(double)) );
      cudaErrorCheckCall( cudaMalloc( (void**)&(d_ms), num_shapes*sizeof(double)) );
      cudaErrorCheckCall( cudaMalloc( (void**)&(d_ns), num_shapes*sizeof(double)) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_allsteps_lsts),
                                          num_visis*sizeof(double)) );
      cudaErrorCheckCall( cudaMemcpy( d_allsteps_lsts, visibility_set->allsteps_lsts,
                             num_visis*sizeof(double),
                             cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&d_u_shapes,
                          num_shapes*num_visis*sizeof(user_precision_t)) );
      cudaErrorCheckCall( cudaMalloc( (void**)&d_v_shapes,
                          num_shapes*num_visis*sizeof(user_precision_t)) );
      cudaErrorCheckCall( cudaMalloc( (void**)&d_w_shapes,
                          num_shapes*num_visis*sizeof(user_precision_t)) );

      //Only the FEE beam currently yields cross pol values, so only malloc what
      //we need here
      if (beam_settings->beamtype == FEE_BEAM || beam_settings->beamtype == FEE_BEAM_INTERP || beam_settings->beamtype == MWA_ANALY) {
        cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J01,
                  catsource.num_shape_primarybeam_values*sizeof(cuUserComplex)) );
        cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J10,
                  catsource.num_shape_primarybeam_values*sizeof(cuUserComplex)) );
      }

      cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J00,
                catsource.num_shape_primarybeam_values*sizeof(cuUserComplex)) );
      cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J11,
                catsource.num_shape_primarybeam_values*sizeof(cuUserComplex)) );


      source_component_common(num_shapes,
           d_primay_beam_J00, d_primay_beam_J01,
           d_primay_beam_J10, d_primay_beam_J11,
           d_freqs, d_ls, d_ms, d_ns,
           d_shape_ras, d_shape_decs,
           catsource.shape_azs, catsource.shape_zas,
           catsource.sin_shape_para_angs, catsource.cos_shape_para_angs,
           catsource.shape_gaussbeam_has,
           catsource.shape_gaussbeam_decs,
           woden_settings, beam_settings);


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

      cudaErrorCheckKernel("kern_calc_uvw_shapelet",
                            kern_calc_uvw_shapelet, grid, threads,
                            d_X_diff, d_Y_diff, d_Z_diff,
                            d_u_shapes, d_v_shapes, d_w_shapes,
                            d_allsteps_wavelengths,
                            d_allsteps_lsts, d_shape_ras, d_shape_decs,
                            num_baselines, num_visis, num_shapes);

      //Splitting over visibilities, but looping over shapelet coeffs inside
      //the kernel
      threads.x = 128;
      threads.y = 1;
      grid.x = (int)ceil( (float)num_visis / (float)threads.x );
      grid.y = 1;

      cudaErrorCheckKernel("kern_calc_visi_shapelets",
              kern_calc_visi_shapelets, grid, threads,
              d_shape_freqs, d_shape_stokesI, d_shape_stokesQ,
              d_shape_stokesU, d_shape_stokesV, d_shape_SIs,
              d_us, d_vs, d_ws, d_allsteps_wavelengths,
              d_u_shapes, d_v_shapes, d_w_shapes,
              d_sum_visi_XX_real, d_sum_visi_XX_imag,
              d_sum_visi_XY_real, d_sum_visi_XY_imag,
              d_sum_visi_YX_real, d_sum_visi_YX_imag,
              d_sum_visi_YY_real, d_sum_visi_YY_imag,
              d_shape_pas, d_shape_majors, d_shape_minors,
              d_shape_n1s, d_shape_n2s, d_shape_coeffs, d_shape_param_indexes,
              d_ls, d_ms, d_ns,
              d_sbf,
              num_shapes, num_baselines, num_freqs, num_visis,
              catsource.n_shape_coeffs, num_time_steps, beam_settings->beamtype,
              d_primay_beam_J00, d_primay_beam_J01,
              d_primay_beam_J10, d_primay_beam_J11);

      cudaErrorCheckCall( cudaFree(d_primay_beam_J00) );
      cudaErrorCheckCall( cudaFree(d_primay_beam_J11) );

      if (beam_settings->beamtype == FEE_BEAM || beam_settings->beamtype == FEE_BEAM_INTERP || beam_settings->beamtype == MWA_ANALY) {
        cudaErrorCheckCall( cudaFree( d_primay_beam_J01 ) );
        cudaErrorCheckCall( cudaFree( d_primay_beam_J10 ) );
      }

      cudaErrorCheckCall( cudaFree(d_ns) );
      cudaErrorCheckCall( cudaFree(d_ms) );
      cudaErrorCheckCall( cudaFree(d_ls) );
      cudaErrorCheckCall( cudaFree(d_w_shapes) );
      cudaErrorCheckCall( cudaFree(d_v_shapes) );
      cudaErrorCheckCall( cudaFree(d_u_shapes) );
      cudaErrorCheckCall( cudaFree(d_allsteps_lsts) );
      cudaErrorCheckCall( cudaFree(d_shape_param_indexes) );
      cudaErrorCheckCall( cudaFree(d_shape_n2s) );
      cudaErrorCheckCall( cudaFree(d_shape_n1s) );
      cudaErrorCheckCall( cudaFree(d_shape_coeffs) );
      cudaErrorCheckCall( cudaFree(d_shape_minors ) );
      cudaErrorCheckCall( cudaFree(d_shape_majors) );
      cudaErrorCheckCall( cudaFree(d_shape_pas) );
      cudaErrorCheckCall( cudaFree(d_shape_decs) );
      cudaErrorCheckCall( cudaFree(d_shape_ras) );
      cudaErrorCheckCall( cudaFree(d_shape_freqs ) );
      cudaErrorCheckCall( cudaFree(d_shape_stokesI ) );
      cudaErrorCheckCall( cudaFree(d_shape_stokesQ ) );
      cudaErrorCheckCall( cudaFree(d_shape_stokesU ) );
      cudaErrorCheckCall( cudaFree(d_shape_stokesV ) );
      cudaErrorCheckCall( cudaFree(d_shape_SIs ) );

    }//if shapelet

    cudaErrorCheckCall( cudaMemcpy(chunk_visibility_set->sum_visi_XX_real,
                       d_sum_visi_XX_real,num_visis*sizeof(user_precision_t),
                                                     cudaMemcpyDeviceToHost) );
    cudaErrorCheckCall( cudaMemcpy(chunk_visibility_set->sum_visi_XX_imag,
                       d_sum_visi_XX_imag,num_visis*sizeof(user_precision_t),
                                                     cudaMemcpyDeviceToHost) );

    cudaErrorCheckCall( cudaMemcpy(chunk_visibility_set->sum_visi_XY_real,
                       d_sum_visi_XY_real,num_visis*sizeof(user_precision_t),
                                                     cudaMemcpyDeviceToHost) );
    cudaErrorCheckCall( cudaMemcpy(chunk_visibility_set->sum_visi_XY_imag,
                       d_sum_visi_XY_imag,num_visis*sizeof(user_precision_t),
                                                     cudaMemcpyDeviceToHost) );

    cudaErrorCheckCall( cudaMemcpy(chunk_visibility_set->sum_visi_YX_real,
                       d_sum_visi_YX_real,num_visis*sizeof(user_precision_t),
                                                     cudaMemcpyDeviceToHost) );
    cudaErrorCheckCall( cudaMemcpy(chunk_visibility_set->sum_visi_YX_imag,
                       d_sum_visi_YX_imag,num_visis*sizeof(user_precision_t),
                                                     cudaMemcpyDeviceToHost) );

    cudaErrorCheckCall( cudaMemcpy(chunk_visibility_set->sum_visi_YY_real,
                       d_sum_visi_YY_real,num_visis*sizeof(user_precision_t),
                                                     cudaMemcpyDeviceToHost) );
    cudaErrorCheckCall( cudaMemcpy(chunk_visibility_set->sum_visi_YY_imag,
                       d_sum_visi_YY_imag,num_visis*sizeof(user_precision_t),
                                                     cudaMemcpyDeviceToHost) );

    cudaErrorCheckCall( cudaMemcpy(chunk_visibility_set->us_metres,
                                  d_u_metres,num_visis*sizeof(user_precision_t),
                                                      cudaMemcpyDeviceToHost) );
    cudaErrorCheckCall( cudaMemcpy(chunk_visibility_set->vs_metres,
                                  d_v_metres,num_visis*sizeof(user_precision_t),
                                                      cudaMemcpyDeviceToHost) );
    cudaErrorCheckCall( cudaMemcpy(chunk_visibility_set->ws_metres,
                                  d_w_metres,num_visis*sizeof(user_precision_t),
                                                      cudaMemcpyDeviceToHost) );

    //add to visiblity_set
    for (int visi = 0; visi < num_visis; visi++) {
      //if the first chunk then initialise our values, and copy across
      //the u,v,w coords
      if (chunk == 0) {
        //ensure temp visi's are 0.0
        visibility_set->sum_visi_XX_real[visi] = 0;
        visibility_set->sum_visi_XX_imag[visi] = 0;
        visibility_set->sum_visi_XY_real[visi] = 0;
        visibility_set->sum_visi_XY_imag[visi] = 0;
        visibility_set->sum_visi_YX_real[visi] = 0;
        visibility_set->sum_visi_YX_imag[visi] = 0;
        visibility_set->sum_visi_YY_real[visi] = 0;
        visibility_set->sum_visi_YY_imag[visi] = 0;

        visibility_set->us_metres[visi] = chunk_visibility_set->us_metres[visi];
        visibility_set->vs_metres[visi] = chunk_visibility_set->vs_metres[visi];
        visibility_set->ws_metres[visi] = chunk_visibility_set->ws_metres[visi];
      }

      //add each chunk of components to visibility set
      visibility_set->sum_visi_XX_real[visi] += chunk_visibility_set->sum_visi_XX_real[visi];
      visibility_set->sum_visi_XX_imag[visi] += chunk_visibility_set->sum_visi_XX_imag[visi];
      visibility_set->sum_visi_XY_real[visi] += chunk_visibility_set->sum_visi_XY_real[visi];
      visibility_set->sum_visi_XY_imag[visi] += chunk_visibility_set->sum_visi_XY_imag[visi];
      visibility_set->sum_visi_YX_real[visi] += chunk_visibility_set->sum_visi_YX_real[visi];
      visibility_set->sum_visi_YX_imag[visi] += chunk_visibility_set->sum_visi_YX_imag[visi];
      visibility_set->sum_visi_YY_real[visi] += chunk_visibility_set->sum_visi_YY_real[visi];
      visibility_set->sum_visi_YY_imag[visi] += chunk_visibility_set->sum_visi_YY_imag[visi];

    }//visi loop

  } //chunk loop

  //Free the chunk_visibility_set
  free_visi_set_outputs(chunk_visibility_set);
  free( chunk_visibility_set );
  //We don't call free_visi_set_inputs as arrays free-ed by that function are
  //just pointers to the original `visiblity_set` in this instance. Only call
  //that in woden::main

  //Free up the GPU memory
  cudaErrorCheckCall( cudaFree(d_freqs) );

  cudaErrorCheckCall( cudaFree(d_ws) );
  cudaErrorCheckCall( cudaFree(d_vs) );
  cudaErrorCheckCall( cudaFree(d_us) );
  cudaErrorCheckCall( cudaFree(d_w_metres) );
  cudaErrorCheckCall( cudaFree(d_v_metres) );
  cudaErrorCheckCall( cudaFree(d_u_metres) );
  cudaErrorCheckCall( cudaFree(d_allsteps_wavelengths) );
  cudaErrorCheckCall( cudaFree(d_allsteps_cha0s) );
  cudaErrorCheckCall( cudaFree(d_allsteps_sha0s) );
  cudaErrorCheckCall( cudaFree(d_Z_diff) );
  cudaErrorCheckCall( cudaFree(d_Y_diff) );
  cudaErrorCheckCall( cudaFree(d_X_diff) );

  if (cropped_sky_models->num_shapelets > 0) {
    cudaErrorCheckCall( cudaFree( d_sbf ) );
  }

  if (beam_settings->beamtype == FEE_BEAM || beam_settings->beamtype == FEE_BEAM_INTERP) {
    free_cuda_fee_beam(beam_settings->cuda_fee_beam);
    free(beam_settings->hyper_delays);
  }

  cudaErrorCheckCall( cudaFree(d_sum_visi_XX_imag) );
  cudaErrorCheckCall( cudaFree(d_sum_visi_XX_real) );
  cudaErrorCheckCall( cudaFree(d_sum_visi_XY_imag) );
  cudaErrorCheckCall( cudaFree(d_sum_visi_XY_real) );
  cudaErrorCheckCall( cudaFree(d_sum_visi_YX_imag) );
  cudaErrorCheckCall( cudaFree(d_sum_visi_YX_real) );
  cudaErrorCheckCall( cudaFree(d_sum_visi_YY_imag) );
  cudaErrorCheckCall( cudaFree(d_sum_visi_YY_real) );

}
