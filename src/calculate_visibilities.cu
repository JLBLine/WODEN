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
#include "primary_beam_cuda.h"
#include "FEE_primary_beam_cuda.h"
#include "cudacheck.h"

extern "C" void calculate_visibilities(array_layout_t * array_layout,
  source_catalogue_t *cropped_sky_models,
  woden_settings_t *woden_settings, visibility_set_t *visibility_set,
  float *sbf) {

  const int num_baselines = woden_settings->num_baselines;
  const int num_time_steps = woden_settings->num_time_steps;
  const int num_visis = woden_settings->num_visis;
  const int num_freqs = woden_settings->num_freqs;

  //Setup chunk_visibility_set to hold the visibility outputs of each
  //sky chunk simulation
  visibility_set_t *chunk_visibility_set = (visibility_set_t *)malloc(sizeof(visibility_set_t));

  chunk_visibility_set->us_metres = (float *)malloc( num_visis * sizeof(float) );
  chunk_visibility_set->vs_metres = (float *)malloc( num_visis * sizeof(float) );
  chunk_visibility_set->ws_metres = (float *)malloc( num_visis * sizeof(float) );

  chunk_visibility_set->allsteps_sha0s = visibility_set->allsteps_sha0s;
  chunk_visibility_set->allsteps_cha0s = visibility_set->allsteps_cha0s;
  chunk_visibility_set->allsteps_lsts = visibility_set->allsteps_lsts;
  chunk_visibility_set->allsteps_wavelengths = visibility_set->allsteps_wavelengths;
  chunk_visibility_set->channel_frequencies = visibility_set->channel_frequencies;

  chunk_visibility_set->sum_visi_XX_real = (float *)malloc( num_visis * sizeof(float) );
  chunk_visibility_set->sum_visi_XX_imag = (float *)malloc( num_visis * sizeof(float) );
  chunk_visibility_set->sum_visi_XY_real = (float *)malloc( num_visis * sizeof(float) );
  chunk_visibility_set->sum_visi_XY_imag = (float *)malloc( num_visis * sizeof(float) );
  chunk_visibility_set->sum_visi_YX_real = (float *)malloc( num_visis * sizeof(float) );
  chunk_visibility_set->sum_visi_YX_imag = (float *)malloc( num_visis * sizeof(float) );
  chunk_visibility_set->sum_visi_YY_real = (float *)malloc( num_visis * sizeof(float) );
  chunk_visibility_set->sum_visi_YY_imag = (float *)malloc( num_visis * sizeof(float) );

  float *d_X_diff = NULL;
  float *d_Y_diff = NULL;
  float *d_Z_diff = NULL;

  cudaErrorCheckCall( cudaMalloc( (void**)&d_X_diff, num_baselines*sizeof(float) ) );
  cudaErrorCheckCall( cudaMemcpy( d_X_diff, array_layout->X_diff_metres,
                      num_baselines*sizeof(float), cudaMemcpyHostToDevice ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_Y_diff, num_baselines*sizeof(float) ) );
  cudaErrorCheckCall( cudaMemcpy( d_Y_diff, array_layout->Y_diff_metres,
                      num_baselines*sizeof(float), cudaMemcpyHostToDevice ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_Z_diff, num_baselines*sizeof(float) ) );
  cudaErrorCheckCall( cudaMemcpy( d_Z_diff, array_layout->Z_diff_metres,
                      num_baselines*sizeof(float), cudaMemcpyHostToDevice ) );

  float *d_allsteps_sha0s = NULL;
  float *d_allsteps_cha0s = NULL;
  float *d_allsteps_wavelengths = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_allsteps_sha0s, num_visis*sizeof(float) ) );
  cudaErrorCheckCall( cudaMemcpy( d_allsteps_sha0s, visibility_set->allsteps_sha0s,
                      num_visis*sizeof(float), cudaMemcpyHostToDevice ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_allsteps_cha0s, num_visis*sizeof(float) ) );
  cudaErrorCheckCall( cudaMemcpy( d_allsteps_cha0s, visibility_set->allsteps_cha0s,
                      num_visis*sizeof(float), cudaMemcpyHostToDevice ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_allsteps_wavelengths, num_visis*sizeof(float) ) );
  cudaErrorCheckCall( cudaMemcpy( d_allsteps_wavelengths, visibility_set->allsteps_wavelengths,
                      num_visis*sizeof(float), cudaMemcpyHostToDevice ) );

  float *d_u_metres = NULL;
  float *d_v_metres = NULL;
  float *d_w_metres = NULL;
  float *d_us = NULL;
  float *d_vs = NULL;
  float *d_ws = NULL;

  cudaErrorCheckCall( cudaMalloc( (void**)&d_u_metres, num_visis*sizeof(float) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_v_metres, num_visis*sizeof(float) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_w_metres, num_visis*sizeof(float) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_us, num_visis*sizeof(float) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_vs, num_visis*sizeof(float) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_ws, num_visis*sizeof(float) ) );

  float *d_sum_visi_XX_real;
  float *d_sum_visi_XX_imag;
  float *d_sum_visi_XY_real;
  float *d_sum_visi_XY_imag;
  float *d_sum_visi_YX_real;
  float *d_sum_visi_YX_imag;
  float *d_sum_visi_YY_real;
  float *d_sum_visi_YY_imag;

  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_XX_real,
                      num_visis*sizeof(float) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_XX_imag,
                      num_visis*sizeof(float) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_XY_real,
                      num_visis*sizeof(float) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_XY_imag,
                      num_visis*sizeof(float) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_YX_real,
                      num_visis*sizeof(float) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_YX_imag,
                      num_visis*sizeof(float) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_YY_real,
                      num_visis*sizeof(float) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_sum_visi_YY_imag,
                      num_visis*sizeof(float) ) );

  float *d_freqs = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_freqs, num_freqs*sizeof(float) ) );
  cudaErrorCheckCall( cudaMemcpy( d_freqs, visibility_set->channel_frequencies,
                      num_freqs*sizeof(float), cudaMemcpyHostToDevice ) );

  //if we have shapelets in our sky model, copy the shapelet basis functions
  //into GPU memory
  float *d_sbf=NULL;
  if (cropped_sky_models->num_shapelets > 0) {
    cudaErrorCheckCall( cudaMalloc( (void**)&(d_sbf), sbf_N*sbf_L*sizeof(float) ));
    cudaErrorCheckCall( cudaMemcpy( d_sbf, sbf, sbf_N*sbf_L*sizeof(float),
                        cudaMemcpyHostToDevice ));
  }

  RTS_MWA_FEE_beam_t *FEE_beam;
  RTS_MWA_FEE_beam_t *FEE_beam_zenith;
  beam_settings_t base_beam_settings;

  if (woden_settings->beamtype == FEE_BEAM){

    base_beam_settings = cropped_sky_models->beam_settings[0];
    FEE_beam = base_beam_settings.FEE_beam;
    FEE_beam_zenith = base_beam_settings.FEE_beam_zenith;

    //Only needs to be done once per frequency as zenith never moves, so do
    //this outside the sky model chunk loop
    printf("Getting FEE beam normalisation...\n");
    get_HDFBeam_normalisation(FEE_beam_zenith, FEE_beam);
    //Free the zenith pointing as done with it now
    free_FEE_primary_beam_from_GPU(base_beam_settings.FEE_beam_zenith);
    printf(" done.\n");

    printf("Copying the FEE beam across to the GPU...");
    copy_FEE_primary_beam_to_GPU(FEE_beam);
    printf(" done.\n");
  }

  //Iterate through all sky model chunks, calculated visibilities are
  //added to chunk_visibility_set, and then summed onto visibility_set
  for (int chunk = 0; chunk < cropped_sky_models->num_sources; chunk++) {

    catsource_t catsource;
    catsource = cropped_sky_models->catsources[chunk];

    beam_settings_t beam_settings;
    beam_settings = cropped_sky_models->beam_settings[chunk];

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
              cropped_sky_models->catsources[chunk].n_points,
              cropped_sky_models->catsources[chunk].n_gauss,
              cropped_sky_models->catsources[chunk].n_shape_coeffs );

    //ensure d_sum_visi_XX_real are set entirely to zero by copying the host
    //array values, which have been set explictly to zero during chunking
    cudaErrorCheckCall( cudaMemcpy(d_sum_visi_XX_real,
               chunk_visibility_set->sum_visi_XX_real,
               num_visis*sizeof(float), cudaMemcpyHostToDevice ) );
    cudaErrorCheckCall( cudaMemcpy(d_sum_visi_XX_imag,
               chunk_visibility_set->sum_visi_XX_imag,
               num_visis*sizeof(float), cudaMemcpyHostToDevice ) );
    cudaErrorCheckCall( cudaMemcpy(d_sum_visi_XY_real,
               chunk_visibility_set->sum_visi_XY_real,
               num_visis*sizeof(float), cudaMemcpyHostToDevice ) );
    cudaErrorCheckCall( cudaMemcpy(d_sum_visi_XY_imag,
               chunk_visibility_set->sum_visi_XY_imag,
               num_visis*sizeof(float), cudaMemcpyHostToDevice ) );
    cudaErrorCheckCall( cudaMemcpy(d_sum_visi_YX_real,
               chunk_visibility_set->sum_visi_YX_real,
               num_visis*sizeof(float), cudaMemcpyHostToDevice ) );
    cudaErrorCheckCall( cudaMemcpy(d_sum_visi_YX_imag,
               chunk_visibility_set->sum_visi_YX_imag,
               num_visis*sizeof(float), cudaMemcpyHostToDevice ) );
    cudaErrorCheckCall( cudaMemcpy(d_sum_visi_YY_real,
               chunk_visibility_set->sum_visi_YY_real,
               num_visis*sizeof(float), cudaMemcpyHostToDevice ) );
    cudaErrorCheckCall( cudaMemcpy(d_sum_visi_YY_imag,
               chunk_visibility_set->sum_visi_YY_imag,
               num_visis*sizeof(float), cudaMemcpyHostToDevice ) );

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
            woden_settings->sdec0, woden_settings->cdec0,
            d_allsteps_cha0s, d_allsteps_sha0s,
            num_visis, num_baselines);

    int num_points = catsource.n_points;
    int num_gauss = catsource.n_gauss;
    int num_shapes = catsource.n_shapes;

    //All components types needs beam Jones matrices, so declare here
    cuFloatComplex *d_primay_beam_J00 = NULL;
    cuFloatComplex *d_primay_beam_J01 = NULL;
    cuFloatComplex *d_primay_beam_J10 = NULL;
    cuFloatComplex *d_primay_beam_J11 = NULL;
    //Same goes for l,m,n coords
    float *d_ls=NULL;
    float *d_ms=NULL;
    float *d_ns=NULL;

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

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_point_ras),
                          num_points*sizeof(float) ) );
      cudaErrorCheckCall( cudaMemcpy( d_point_ras, catsource.point_ras,
                          num_points*sizeof(float), cudaMemcpyHostToDevice ) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_point_decs),
                          num_points*sizeof(float) ) );
      cudaErrorCheckCall( cudaMemcpy( d_point_decs, catsource.point_decs,
                          num_points*sizeof(float), cudaMemcpyHostToDevice ) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_point_freqs),
                          num_points*sizeof(float) ) );
      cudaErrorCheckCall( cudaMemcpy( d_point_freqs, catsource.point_ref_freqs,
                          num_points*sizeof(float), cudaMemcpyHostToDevice ) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_point_stokesI),
                          num_points*sizeof(float) ) );
      cudaErrorCheckCall( cudaMemcpy( d_point_stokesI, catsource.point_ref_stokesI,
                          num_points*sizeof(float), cudaMemcpyHostToDevice ) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_point_stokesQ),
                          num_points*sizeof(float) ) );
      cudaErrorCheckCall( cudaMemcpy( d_point_stokesQ, catsource.point_ref_stokesQ,
                          num_points*sizeof(float), cudaMemcpyHostToDevice ) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_point_stokesU),
                          num_points*sizeof(float) ) );
      cudaErrorCheckCall( cudaMemcpy( d_point_stokesU, catsource.point_ref_stokesU,
                          num_points*sizeof(float), cudaMemcpyHostToDevice ) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_point_stokesV),
                          num_points*sizeof(float) ) );
      cudaErrorCheckCall( cudaMemcpy( d_point_stokesV, catsource.point_ref_stokesV,
                          num_points*sizeof(float), cudaMemcpyHostToDevice ) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_point_SIs),
                          num_points*sizeof(float) ) );
      cudaErrorCheckCall( cudaMemcpy( d_point_SIs, catsource.point_SIs,
                          num_points*sizeof(float), cudaMemcpyHostToDevice ) );
      //
      //Only the FEE beam currently yields cross pol values, so only malloc what
      //we need here
      if (beam_settings.beamtype == FEE_BEAM) {
        cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J01,
                  beam_settings.num_point_beam_values*sizeof(cuFloatComplex) ));
        cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J10,
                  beam_settings.num_point_beam_values*sizeof(cuFloatComplex) ));
      }

      cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J00,
                  beam_settings.num_point_beam_values*sizeof(cuFloatComplex) ));
      cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J11,
                  beam_settings.num_point_beam_values*sizeof(cuFloatComplex) ));
      //
      cudaErrorCheckCall( cudaMalloc( (void**)&d_ls, num_points*sizeof(float) ) );
      cudaErrorCheckCall( cudaMalloc( (void**)&d_ms, num_points*sizeof(float) ) );
      cudaErrorCheckCall( cudaMalloc( (void**)&d_ns, num_points*sizeof(float) ) );

      //This calculates l,m,n for all components, and any beam calculations
      //that are needed - uses the same methods for all component types
      source_component_common(num_points,
                 d_primay_beam_J00, d_primay_beam_J01,
                 d_primay_beam_J10, d_primay_beam_J11,
                 d_freqs, d_ls, d_ms, d_ns,
                 d_point_ras, d_point_decs,
                 catsource.point_azs, catsource.point_zas,
                 catsource.sin_point_para_angs, catsource.cos_point_para_angs,
                 beam_settings.beam_point_has, beam_settings.beam_point_decs,
                 woden_settings, beam_settings, FEE_beam);


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

      cudaErrorCheckKernel("kern_calc_visi_point",
                           kern_calc_visi_point, grid, threads,
                           d_point_ras, d_point_decs,
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
                           num_time_steps, beam_settings.beamtype,
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

      if (beam_settings.beamtype == FEE_BEAM){
        cudaErrorCheckCall( cudaFree( FEE_beam->d_FEE_beam_gain_matrices) );
        cudaErrorCheckCall( cudaFree( d_primay_beam_J01 ) );
        cudaErrorCheckCall( cudaFree( d_primay_beam_J10 ) );
      }

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

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_gauss_freqs),
                          num_gauss*sizeof(float)) );
      cudaErrorCheckCall( cudaMemcpy( d_gauss_freqs, catsource.gauss_ref_freqs,
                          num_gauss*sizeof(float), cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_gauss_stokesI),
                          num_gauss*sizeof(float)) );
      cudaErrorCheckCall( cudaMemcpy( d_gauss_stokesI, catsource.gauss_ref_stokesI,
                          num_gauss*sizeof(float), cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_gauss_stokesQ),
                          num_gauss*sizeof(float)) );
      cudaErrorCheckCall( cudaMemcpy( d_gauss_stokesQ, catsource.gauss_ref_stokesQ,
                          num_gauss*sizeof(float), cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_gauss_stokesU),
                          num_gauss*sizeof(float)) );
      cudaErrorCheckCall( cudaMemcpy( d_gauss_stokesU, catsource.gauss_ref_stokesU,
                          num_gauss*sizeof(float), cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_gauss_stokesV),
                          num_gauss*sizeof(float)) );
      cudaErrorCheckCall( cudaMemcpy( d_gauss_stokesV, catsource.gauss_ref_stokesV,
                          num_gauss*sizeof(float), cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_gauss_SIs),
                          num_gauss*sizeof(float)) );
      cudaErrorCheckCall( cudaMemcpy( d_gauss_SIs, catsource.gauss_SIs,
                          num_gauss*sizeof(float), cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_gauss_ras),
                          num_gauss*sizeof(float)) );
      cudaErrorCheckCall( cudaMemcpy( d_gauss_ras, catsource.gauss_ras,
                          num_gauss*sizeof(float), cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_gauss_decs),
                          num_gauss*sizeof(float)) );
      cudaErrorCheckCall( cudaMemcpy( d_gauss_decs, catsource.gauss_decs,
                          num_gauss*sizeof(float), cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_gauss_pas),
                          num_gauss*sizeof(float)) );
      cudaErrorCheckCall( cudaMemcpy( d_gauss_pas, catsource.gauss_pas,
                          num_gauss*sizeof(float), cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_gauss_majors),
                          num_gauss*sizeof(float)) );
      cudaErrorCheckCall( cudaMemcpy( d_gauss_majors, catsource.gauss_majors,
                          num_gauss*sizeof(float), cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_gauss_minors),
                          num_gauss*sizeof(float)) );
      cudaErrorCheckCall( cudaMemcpy( d_gauss_minors, catsource.gauss_minors,
                          num_gauss*sizeof(float), cudaMemcpyHostToDevice) );

      //Only the FEE beam currently yields cross pol values, so only malloc what
      //we need here
      if (beam_settings.beamtype == FEE_BEAM) {
        cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J01,
              beam_settings.num_gausscomp_beam_values*sizeof(cuFloatComplex)) );
        cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J10,
              beam_settings.num_gausscomp_beam_values*sizeof(cuFloatComplex)) );
      }

      cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J00,
            beam_settings.num_gausscomp_beam_values*sizeof(cuFloatComplex)) );
      cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J11,
            beam_settings.num_gausscomp_beam_values*sizeof(cuFloatComplex)) );

      cudaErrorCheckCall( cudaMalloc( (void**)&d_ls, num_gauss*sizeof(float)) );
      cudaErrorCheckCall( cudaMalloc( (void**)&d_ms, num_gauss*sizeof(float)) );
      cudaErrorCheckCall( cudaMalloc( (void**)&d_ns, num_gauss*sizeof(float)) );

      source_component_common(num_gauss,
                 d_primay_beam_J00, d_primay_beam_J01,
                 d_primay_beam_J10, d_primay_beam_J11,
                 d_freqs, d_ls, d_ms, d_ns,
                 d_gauss_ras, d_gauss_decs,
                 catsource.gauss_azs, catsource.gauss_zas,
                 catsource.sin_gauss_para_angs, catsource.cos_gauss_para_angs,
                 beam_settings.beam_gausscomp_has,
                 beam_settings.beam_gausscomp_decs,
                 woden_settings, beam_settings, FEE_beam);

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

      cudaErrorCheckKernel("kern_calc_visi_gaussian",
              kern_calc_visi_gaussian, grid, threads,
              d_gauss_ras, d_gauss_decs,
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
              num_time_steps, beam_settings.beamtype,
              d_primay_beam_J00, d_primay_beam_J01,
              d_primay_beam_J10, d_primay_beam_J11);

      cudaErrorCheckCall( cudaFree(d_primay_beam_J00) );
      cudaErrorCheckCall( cudaFree(d_primay_beam_J11) );

      if (beam_settings.beamtype == FEE_BEAM){
        cudaErrorCheckCall( cudaFree(FEE_beam->d_FEE_beam_gain_matrices) );
        cudaErrorCheckCall( cudaFree(d_primay_beam_J01) );
        cudaErrorCheckCall( cudaFree(d_primay_beam_J10) );
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

      float *d_shape_ras=NULL;
      float *d_shape_decs=NULL;
      float *d_shape_freqs=NULL;
      float *d_shape_stokesI=NULL;
      float *d_shape_stokesQ=NULL;
      float *d_shape_stokesU=NULL;
      float *d_shape_stokesV=NULL;
      float *d_shape_SIs=NULL;
      float *d_shape_pas=NULL;
      float *d_shape_majors=NULL;
      float *d_shape_minors=NULL;

      float *d_shape_coeffs=NULL;
      float *d_shape_n1s=NULL;
      float *d_shape_n2s=NULL;
      float *d_shape_param_indexes=NULL;

      float *d_allsteps_lsts=NULL;

      float *d_u_s_metres = NULL;
      float *d_v_s_metres = NULL;
      float *d_w_s_metres = NULL;

      //Who likes cudaMalloc cudaMallocs? We like cudaMalloc cudaMallocs
      cudaErrorCheckCall( cudaMalloc( (void**)&(d_shape_ras),
                          num_shapes*sizeof(float)) );
      cudaErrorCheckCall( cudaMemcpy( d_shape_ras, catsource.shape_ras,
              num_shapes*sizeof(float), cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_shape_decs),
                          num_shapes*sizeof(float)) );
      cudaErrorCheckCall( cudaMemcpy( d_shape_decs, catsource.shape_decs,
              num_shapes*sizeof(float), cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_shape_freqs),
                          num_shapes*sizeof(float)) );
      cudaErrorCheckCall( cudaMemcpy( d_shape_freqs, catsource.shape_ref_freqs,
                          num_shapes*sizeof(float), cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_shape_stokesI),
                          num_shapes*sizeof(float)) );
      cudaErrorCheckCall( cudaMemcpy( d_shape_stokesI, catsource.shape_ref_stokesI,
                          num_shapes*sizeof(float), cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_shape_stokesQ),
                          num_shapes*sizeof(float)) );
      cudaErrorCheckCall( cudaMemcpy( d_shape_stokesQ, catsource.shape_ref_stokesQ,
                          num_shapes*sizeof(float), cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_shape_stokesU),
                          num_shapes*sizeof(float)) );
      cudaErrorCheckCall( cudaMemcpy( d_shape_stokesU, catsource.shape_ref_stokesU,
                          num_shapes*sizeof(float), cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_shape_stokesV),
                          num_shapes*sizeof(float)) );
      cudaErrorCheckCall( cudaMemcpy( d_shape_stokesV, catsource.shape_ref_stokesV,
                          num_shapes*sizeof(float), cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_shape_SIs),
                          num_shapes*sizeof(float)) );
      cudaErrorCheckCall( cudaMemcpy( d_shape_SIs, catsource.shape_SIs,
                          num_shapes*sizeof(float), cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_shape_pas),
                          num_shapes*sizeof(float)) );
      cudaErrorCheckCall( cudaMemcpy( d_shape_pas, catsource.shape_pas,
              num_shapes*sizeof(float), cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_shape_majors),
                          num_shapes*sizeof(float)) );
      cudaErrorCheckCall( cudaMemcpy( d_shape_majors, catsource.shape_majors,
              num_shapes*sizeof(float), cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_shape_minors),
                          num_shapes*sizeof(float)) );
      cudaErrorCheckCall( cudaMemcpy( d_shape_minors, catsource.shape_minors,
              num_shapes*sizeof(float), cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_shape_coeffs),
                          catsource.n_shape_coeffs*sizeof(float)) );
      cudaErrorCheckCall( cudaMemcpy( d_shape_coeffs, catsource.shape_coeffs,
              catsource.n_shape_coeffs*sizeof(float), cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_shape_n1s),
                          catsource.n_shape_coeffs*sizeof(float)) );
      cudaErrorCheckCall( cudaMemcpy( d_shape_n1s, catsource.shape_n1s,
              catsource.n_shape_coeffs*sizeof(float), cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_shape_n2s),
                          catsource.n_shape_coeffs*sizeof(float)) );
      cudaErrorCheckCall( cudaMemcpy( d_shape_n2s, catsource.shape_n2s,
              catsource.n_shape_coeffs*sizeof(float), cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_shape_param_indexes),
                          catsource.n_shape_coeffs*sizeof(float)) );
      cudaErrorCheckCall( cudaMemcpy( d_shape_param_indexes,
             catsource.shape_param_indexes, catsource.n_shape_coeffs*sizeof(float), cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_ls), num_shapes*sizeof(float)) );
      cudaErrorCheckCall( cudaMalloc( (void**)&(d_ms), num_shapes*sizeof(float)) );
      cudaErrorCheckCall( cudaMalloc( (void**)&(d_ns), num_shapes*sizeof(float)) );

      cudaErrorCheckCall( cudaMalloc( (void**)&(d_allsteps_lsts), num_visis*sizeof(float)) );
      cudaErrorCheckCall( cudaMemcpy( d_allsteps_lsts, visibility_set->allsteps_lsts,
                             num_visis*sizeof(float), cudaMemcpyHostToDevice) );

      cudaErrorCheckCall( cudaMalloc( (void**)&d_u_s_metres,
                          num_shapes*num_visis*sizeof(float)) );
      cudaErrorCheckCall( cudaMalloc( (void**)&d_v_s_metres,
                          num_shapes*num_visis*sizeof(float)) );
      cudaErrorCheckCall( cudaMalloc( (void**)&d_w_s_metres,
                          num_shapes*num_visis*sizeof(float)) );

      //Only the FEE beam currently yields cross pol values, so only malloc what
      //we need here
      if (beam_settings.beamtype == FEE_BEAM) {
        cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J01,
                  beam_settings.num_shape_beam_values*sizeof(cuFloatComplex)) );
        cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J10,
                  beam_settings.num_shape_beam_values*sizeof(cuFloatComplex)) );
      }

      cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J00,
                beam_settings.num_shape_beam_values*sizeof(cuFloatComplex)) );
      cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J11,
                beam_settings.num_shape_beam_values*sizeof(cuFloatComplex)) );


      source_component_common(num_shapes,
           d_primay_beam_J00, d_primay_beam_J01,
           d_primay_beam_J10, d_primay_beam_J11,
           d_freqs, d_ls, d_ms, d_ns,
           d_shape_ras, d_shape_decs,
           catsource.shape_azs, catsource.shape_zas,
           catsource.sin_shape_para_angs, catsource.cos_shape_para_angs,
           beam_settings.beam_shape_has,
           beam_settings.beam_shape_decs,
           woden_settings, beam_settings, FEE_beam);


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
                            d_u_s_metres, d_v_s_metres, d_w_s_metres,
                            d_allsteps_lsts, d_shape_ras, d_shape_decs,
                            num_baselines, num_visis, num_shapes);

      if (catsource.n_shape_coeffs == 1) {
        threads.x = 64;
        threads.y = 1;
        grid.x = (int)ceil( (float)num_visis / (float)threads.x );
        grid.y = 1;
      }
      else {
        threads.x = 64;
        threads.y = 4;
        grid.x = (int)ceil( (float)num_visis / (float)threads.x );
        grid.y = (int)ceil( ((float)catsource.n_shape_coeffs) / ((float)threads.y) );
      }

      cudaErrorCheckKernel("kern_calc_visi_shapelets",
              kern_calc_visi_shapelets, grid, threads,
              d_shape_freqs, d_shape_stokesI, d_shape_stokesQ,
              d_shape_stokesU, d_shape_stokesV, d_shape_SIs,
              d_us, d_vs, d_ws, d_allsteps_wavelengths,
              d_u_s_metres, d_v_s_metres, d_w_s_metres,
              d_sum_visi_XX_real, d_sum_visi_XX_imag,
              d_sum_visi_XY_real, d_sum_visi_XY_imag,
              d_sum_visi_YX_real, d_sum_visi_YX_imag,
              d_sum_visi_YY_real, d_sum_visi_YY_imag,
              d_shape_pas, d_shape_majors, d_shape_minors,
              d_shape_n1s, d_shape_n2s, d_shape_coeffs, d_shape_param_indexes,
              d_ls, d_ms, d_ns,
              d_sbf,
              num_shapes, num_baselines, num_freqs, num_visis,
              catsource.n_shape_coeffs, num_time_steps, beam_settings.beamtype,
              d_primay_beam_J00, d_primay_beam_J01,
              d_primay_beam_J10, d_primay_beam_J11);

      cudaErrorCheckCall( cudaFree(d_primay_beam_J00) );
      cudaErrorCheckCall( cudaFree(d_primay_beam_J11) );

      if (beam_settings.beamtype == FEE_BEAM){
        cudaErrorCheckCall( cudaFree( FEE_beam->d_FEE_beam_gain_matrices) );
        cudaErrorCheckCall( cudaFree(d_primay_beam_J01) );
        cudaErrorCheckCall( cudaFree(d_primay_beam_J10) );
      }

      cudaErrorCheckCall( cudaFree(d_ns) );
      cudaErrorCheckCall( cudaFree(d_ms) );
      cudaErrorCheckCall( cudaFree(d_ls) );
      cudaErrorCheckCall( cudaFree(d_w_s_metres) );
      cudaErrorCheckCall( cudaFree(d_v_s_metres) );
      cudaErrorCheckCall( cudaFree(d_u_s_metres) );
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
           d_sum_visi_XX_real,num_visis*sizeof(float),cudaMemcpyDeviceToHost) );
    cudaErrorCheckCall( cudaMemcpy(chunk_visibility_set->sum_visi_XX_imag,
           d_sum_visi_XX_imag,num_visis*sizeof(float),cudaMemcpyDeviceToHost) );

    cudaErrorCheckCall( cudaMemcpy(chunk_visibility_set->sum_visi_XY_real,
           d_sum_visi_XY_real,num_visis*sizeof(float),cudaMemcpyDeviceToHost) );
    cudaErrorCheckCall( cudaMemcpy(chunk_visibility_set->sum_visi_XY_imag,
           d_sum_visi_XY_imag,num_visis*sizeof(float),cudaMemcpyDeviceToHost) );

    cudaErrorCheckCall( cudaMemcpy(chunk_visibility_set->sum_visi_YX_real,
           d_sum_visi_YX_real,num_visis*sizeof(float),cudaMemcpyDeviceToHost) );
    cudaErrorCheckCall( cudaMemcpy(chunk_visibility_set->sum_visi_YX_imag,
           d_sum_visi_YX_imag,num_visis*sizeof(float),cudaMemcpyDeviceToHost) );

    cudaErrorCheckCall( cudaMemcpy(chunk_visibility_set->sum_visi_YY_real,
           d_sum_visi_YY_real,num_visis*sizeof(float),cudaMemcpyDeviceToHost) );
    cudaErrorCheckCall( cudaMemcpy(chunk_visibility_set->sum_visi_YY_imag,
           d_sum_visi_YY_imag,num_visis*sizeof(float),cudaMemcpyDeviceToHost) );

    cudaErrorCheckCall( cudaMemcpy(chunk_visibility_set->us_metres,
              d_u_metres,num_visis*sizeof(float),cudaMemcpyDeviceToHost) );
    cudaErrorCheckCall( cudaMemcpy(chunk_visibility_set->vs_metres,
              d_v_metres,num_visis*sizeof(float),cudaMemcpyDeviceToHost) );
    cudaErrorCheckCall( cudaMemcpy(chunk_visibility_set->ws_metres,
              d_w_metres,num_visis*sizeof(float),cudaMemcpyDeviceToHost) );

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
  free( chunk_visibility_set->us_metres );
  free( chunk_visibility_set->vs_metres );
  free( chunk_visibility_set->ws_metres );

  free(chunk_visibility_set->sum_visi_XX_real);
  free(chunk_visibility_set->sum_visi_XX_imag);
  free(chunk_visibility_set->sum_visi_XY_real);
  free(chunk_visibility_set->sum_visi_XY_imag);
  free(chunk_visibility_set->sum_visi_YX_real);
  free(chunk_visibility_set->sum_visi_YX_imag);
  free(chunk_visibility_set->sum_visi_YY_real);
  free(chunk_visibility_set->sum_visi_YY_imag);
  //
  free( chunk_visibility_set );

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

  if (woden_settings->beamtype == FEE_BEAM) {
    free_FEE_primary_beam_from_GPU(base_beam_settings.FEE_beam);
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
