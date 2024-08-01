#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

#include "woden_precision_defs.h"
#include "gpucomplex.h"

#include "calculate_visibilities.h"
#include "fundamental_coords.h"
#include "constants.h"
#include "source_components.h"
#include "primary_beam_gpu.h"
#include "hyperbeam_error.h"

//Helpful C code we are also using
#include "visibility_set.h"

#include "gpu_macros.h"

extern "C" void calculate_visibilities(array_layout_t *array_layout,
  source_catalogue_t *cropped_sky_models, beam_settings_t *beam_settings,
  woden_settings_t *woden_settings, visibility_set_t *visibility_set,
  user_precision_t *sbf) {

  //Boolean for if we are using two beams per visibility or assuming all
  //beams are the same
  int use_twobeams = 0;
  int num_beams = 1;
  //Currently only do this is we have the use_dipamps flag; could be expanded
  if (woden_settings->use_dipamps == 1) {
    use_twobeams = 1;
    num_beams = woden_settings->num_ants;
  }
  // else {
  //   use_twobeams = 0;
  //   num_beams = 1;
  // }

  //TODO - once rotation measure has been implemented, this should be set
  //only if we are using a rotation measure

  const int num_baselines = woden_settings->num_baselines;
  const int num_time_steps = woden_settings->num_time_steps;
  const int num_visis = woden_settings->num_visis;
  const int num_autos = woden_settings->num_autos;
  const int num_cross = woden_settings->num_cross;
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

   gpuMalloc( (void**)&d_X_diff,
                                 num_time_steps*num_baselines*sizeof(double) );
   gpuMemcpy( d_X_diff, array_layout->X_diff_metres,
        num_time_steps*num_baselines*sizeof(double), gpuMemcpyHostToDevice );
   gpuMalloc( (void**)&d_Y_diff,
                                 num_time_steps*num_baselines*sizeof(double) );
   gpuMemcpy( d_Y_diff, array_layout->Y_diff_metres,
        num_time_steps*num_baselines*sizeof(double), gpuMemcpyHostToDevice );
   gpuMalloc( (void**)&d_Z_diff,
                                 num_time_steps*num_baselines*sizeof(double) );
   gpuMemcpy( d_Z_diff, array_layout->Z_diff_metres,
        num_time_steps*num_baselines*sizeof(double), gpuMemcpyHostToDevice );

  double *d_allsteps_sha0s = NULL;
  double *d_allsteps_cha0s = NULL;
  user_precision_t *d_allsteps_wavelengths = NULL;
  gpuMalloc( (void**)&d_allsteps_sha0s, num_cross*sizeof(double) );
  gpuMemcpy( d_allsteps_sha0s, visibility_set->allsteps_sha0s,
                      num_cross*sizeof(double), gpuMemcpyHostToDevice );
  gpuMalloc( (void**)&d_allsteps_cha0s, num_cross*sizeof(double) );
  gpuMemcpy( d_allsteps_cha0s, visibility_set->allsteps_cha0s,
                      num_cross*sizeof(double), gpuMemcpyHostToDevice );
  gpuMalloc( (void**)&d_allsteps_wavelengths,
                                         num_cross*sizeof(user_precision_t) );
  gpuMemcpy( d_allsteps_wavelengths, visibility_set->allsteps_wavelengths,
                      num_cross*sizeof(user_precision_t), gpuMemcpyHostToDevice );

  user_precision_t *d_u_metres = NULL;
  user_precision_t *d_v_metres = NULL;
  user_precision_t *d_w_metres = NULL;
  user_precision_t *d_us = NULL;
  user_precision_t *d_vs = NULL;
  user_precision_t *d_ws = NULL;

  ( gpuMalloc( (void**)&d_u_metres, num_visis*sizeof(user_precision_t) ) );
  ( gpuMalloc( (void**)&d_v_metres, num_visis*sizeof(user_precision_t) ) );
  ( gpuMalloc( (void**)&d_w_metres, num_visis*sizeof(user_precision_t) ) );
  ( gpuMalloc( (void**)&d_us, num_visis*sizeof(user_precision_t) ) );
  ( gpuMalloc( (void**)&d_vs, num_visis*sizeof(user_precision_t) ) );
  ( gpuMalloc( (void**)&d_ws, num_visis*sizeof(user_precision_t) ) );

  visibility_set_t *d_visibility_set =  (visibility_set_t* )malloc(sizeof(visibility_set_t));

  ( gpuMalloc( (void**)&d_visibility_set->sum_visi_XX_real,
                      num_visis*sizeof(user_precision_t) ) );
  ( gpuMalloc( (void**)&d_visibility_set->sum_visi_XX_imag,
                      num_visis*sizeof(user_precision_t) ) );
  ( gpuMalloc( (void**)&d_visibility_set->sum_visi_XY_real,
                      num_visis*sizeof(user_precision_t) ) );
  ( gpuMalloc( (void**)&d_visibility_set->sum_visi_XY_imag,
                      num_visis*sizeof(user_precision_t) ) );
  ( gpuMalloc( (void**)&d_visibility_set->sum_visi_YX_real,
                      num_visis*sizeof(user_precision_t) ) );
  ( gpuMalloc( (void**)&d_visibility_set->sum_visi_YX_imag,
                      num_visis*sizeof(user_precision_t) ) );
  ( gpuMalloc( (void**)&d_visibility_set->sum_visi_YY_real,
                      num_visis*sizeof(user_precision_t) ) );
  ( gpuMalloc( (void**)&d_visibility_set->sum_visi_YY_imag,
                      num_visis*sizeof(user_precision_t) ) );

  double *d_freqs = NULL;
  ( gpuMalloc( (void**)&d_freqs, num_freqs*sizeof(double) ) );
  ( gpuMemcpy( d_freqs, visibility_set->channel_frequencies,
                      num_freqs*sizeof(double), gpuMemcpyHostToDevice ) );

  //if we have shapelets in our sky model, copy the shapelet basis functions
  //into GPU memory
  user_precision_t *d_sbf=NULL;
  if (cropped_sky_models->num_shapelets > 0) {
    ( gpuMalloc( (void**)&(d_sbf), sbf_N*sbf_L*sizeof(user_precision_t) ));
    ( gpuMemcpy( d_sbf, sbf, sbf_N*sbf_L*sizeof(user_precision_t),
                        gpuMemcpyHostToDevice ));
  }

  //TODO - this could be a function in some other file
  if (beam_settings->beamtype == FEE_BEAM || beam_settings->beamtype == FEE_BEAM_INTERP){

    //Convert some WODEN precision stuff into hyperbeam precision stuff
    uint32_t *freqs_hz = (uint32_t*)malloc(woden_settings->num_freqs*sizeof(uint32_t));

    for (int freq_ind = 0; freq_ind < woden_settings->num_freqs; freq_ind++) {
        freqs_hz[freq_ind] = (uint32_t)visibility_set->channel_frequencies[freq_ind];
    }

    beam_settings->hyper_delays = (uint32_t*)malloc(16*num_beams*sizeof(uint32_t));

    for (int delay = 0; delay < 16*num_beams; delay++) {
      beam_settings->hyper_delays[delay] = (uint32_t)woden_settings->FEE_ideal_delays[delay];
    }

    uint32_t num_amps;

    //32 means we have amplitudes for both X and Y
    if (use_twobeams == 0) {
      num_amps = 16;
      double mwa_dipole_amps[16] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
                                    1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
      woden_settings->mwa_dipole_amps = mwa_dipole_amps;
      printf("HIP needs this printf otherwise it doesnt work\n");

    } else {
      num_amps = 32;
      //woden_settings->mwa_dipole_amps should already have been set so
      //leave it alone
      //TODO set this outside of calculate_visibilities
    }

    uint8_t norm_to_zenith = 1;


    int32_t status = new_gpu_fee_beam(beam_settings->fee_beam,
                            freqs_hz,
                            beam_settings->hyper_delays,
                            woden_settings->mwa_dipole_amps,
                            woden_settings->num_freqs,
                            (uint32_t)num_beams,
                            num_amps,
                            norm_to_zenith,
                            &beam_settings->gpu_fee_beam);

    if (status != 0) {
      // handle_hyperbeam_error(__FILE__, __LINE__, "new_gpu_fee_beam");
      printf("Something went wrong launching new_gpu_fee_beam\n");
    }
  }

  //If we a different primary beam for each antenna, setup the baseline
  //to anetenna mapping arrays
  int *d_ant1_to_baseline_map = NULL;
  int *d_ant2_to_baseline_map = NULL;

  if (use_twobeams == 1) {

    // printf("WHAT YA DO? woden_settings->num_baselines: %d\n", woden_settings->num_baselines);

    gpuMalloc( (void**)&d_ant1_to_baseline_map,
                                woden_settings->num_baselines*sizeof(int) );

    gpuMalloc( (void**)&d_ant2_to_baseline_map,
                                woden_settings->num_baselines*sizeof(int) );

    fill_ant_to_baseline_mapping(woden_settings->num_ants, d_ant1_to_baseline_map,
                                                           d_ant2_to_baseline_map);

  }


  //Iterate through all sky model chunks, calculated visibilities are
  //added to chunk_visibility_set, and then summed onto visibility_set
  for (int chunk = 0; chunk < cropped_sky_models->num_sources; chunk++) {

    // source_t *source = (source_t *)malloc(sizeof(source_t));

    source_t *source = &cropped_sky_models->sources[chunk];

    // printf("\tsource->n_comps %d\n", source->n_comps);
    // printf("\tsource->n_points %d\n", source->n_points);
    // printf("\tsource->n_point_lists %d\n", source->n_point_lists);
    // printf("\tsource->n_point_powers %d\n", source->n_point_powers);
    // printf("\tsource->n_point_curves %d\n", source->n_point_curves);
    // printf("\tsource->n_gauss %d\n", source->n_gauss);
    // printf("\tsource->n_gauss_lists %d\n", source->n_gauss_lists);
    // printf("\tsource->n_gauss_powers %d\n", source->n_gauss_powers);
    // printf("\tsource->n_gauss_curves %d\n", source->n_gauss_curves);
    // printf("\tsource->n_shapes %d\n", source->n_shapes);
    // printf("\tsource->n_shape_lists %d\n", source->n_shape_lists);
    // printf("\tsource->n_shape_powers %d\n", source->n_shape_powers);
    // printf("\tsource->n_shape_curves %d\n", source->n_shape_curves);
    // printf("\tsource->n_shape_coeffs %d\n", source->n_shape_coeffs);

    printf("About to copy the chunked source to the GPU\n");

    source_t *d_chunked_source = copy_chunked_source_to_GPU(source);
    printf("Have copied across the chunk to the GPU\n");

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

    //TODO setup some kind of logging and have an option to silence all this
    //printed output
    printf("Processing chunk %d\n", chunk);
    printf("\tNumber of components in chunk are: P %d G %d S_coeffs %d\n",
              source->n_points,
              source->n_gauss,
              source->n_shape_coeffs );

    //ensure d_sum_visi_XX_real are set entirely to zero by copying the host
    //array values, which have been set explictly to zero during chunking
    gpuMemcpy(d_visibility_set->sum_visi_XX_real,
               chunk_visibility_set->sum_visi_XX_real,
               num_visis*sizeof(user_precision_t), gpuMemcpyHostToDevice );
    gpuMemcpy(d_visibility_set->sum_visi_XX_imag,
               chunk_visibility_set->sum_visi_XX_imag,
               num_visis*sizeof(user_precision_t), gpuMemcpyHostToDevice );
    gpuMemcpy(d_visibility_set->sum_visi_XY_real,
               chunk_visibility_set->sum_visi_XY_real,
               num_visis*sizeof(user_precision_t), gpuMemcpyHostToDevice );
    gpuMemcpy(d_visibility_set->sum_visi_XY_imag,
               chunk_visibility_set->sum_visi_XY_imag,
               num_visis*sizeof(user_precision_t), gpuMemcpyHostToDevice );
    gpuMemcpy(d_visibility_set->sum_visi_YX_real,
               chunk_visibility_set->sum_visi_YX_real,
               num_visis*sizeof(user_precision_t), gpuMemcpyHostToDevice );
    gpuMemcpy(d_visibility_set->sum_visi_YX_imag,
               chunk_visibility_set->sum_visi_YX_imag,
               num_visis*sizeof(user_precision_t), gpuMemcpyHostToDevice );
    gpuMemcpy(d_visibility_set->sum_visi_YY_real,
               chunk_visibility_set->sum_visi_YY_real,
               num_visis*sizeof(user_precision_t), gpuMemcpyHostToDevice );
    gpuMemcpy(d_visibility_set->sum_visi_YY_imag,
               chunk_visibility_set->sum_visi_YY_imag,
               num_visis*sizeof(user_precision_t), gpuMemcpyHostToDevice );

    dim3 grid, threads;

    threads.x = 128;
    threads.y = 1;
    grid.x = (int)ceil( (float)num_cross / (float)threads.x );
    grid.y = 1;

    gpuErrorCheckKernel("kern_calc_uvw",
            kern_calc_uvw, grid, threads,
            d_X_diff, d_Y_diff, d_Z_diff,
            d_u_metres, d_v_metres, d_w_metres,
            d_us, d_vs, d_ws, d_allsteps_wavelengths,
            woden_settings->sdec0,
            woden_settings->cdec0,
            d_allsteps_cha0s, d_allsteps_sha0s,
            num_cross, num_baselines, num_time_steps, num_freqs);

    //Set auto-correlation uvw to zero if doing autos
    if (woden_settings->do_autos == 1){
      grid.x = (int)ceil( (float)num_autos / (float)threads.x );
      grid.y = threads.y = 1;

      gpuErrorCheckKernel("set_auto_uvw_to_zero",
              set_auto_uvw_to_zero, grid, threads,
              woden_settings->num_cross, woden_settings->num_visis,
              d_us, d_vs, d_ws);

      gpuErrorCheckKernel("set_auto_uvw_to_zero",
              set_auto_uvw_to_zero, grid, threads,
              woden_settings->num_cross, woden_settings->num_visis,
              d_u_metres, d_v_metres, d_w_metres);
    }

    int num_points = source->n_points;
    int num_gauss = source->n_gauss;
    int num_shapes = source->n_shapes;

    if (num_points > 0) {
      printf("\tDoing point components\n");

      //Something to store the primary beam gains (all 4 pols) in
      d_beam_gains_t d_point_beam_gains;
      if (use_twobeams == 1) {
        d_point_beam_gains.d_ant1_to_baseline_map = d_ant1_to_baseline_map;
        d_point_beam_gains.d_ant2_to_baseline_map = d_ant2_to_baseline_map;
        d_point_beam_gains.use_twobeams = 1;
      } else {
        d_point_beam_gains.use_twobeams = 0;
      }

      printf("\tExtrapolating fluxes and beams...\n");
      source_component_common(woden_settings, beam_settings, d_freqs,
                              source, d_chunked_source,
                              &d_point_beam_gains, POINT,
                              d_visibility_set);
      printf("\tExtrapolating fluxes and beams done.\n");

      threads.x = 128;
      threads.y = 1;
      grid.x = (int)ceil( (float)num_cross / (float)threads.x );
      grid.y = 1;

      printf("\tDoing visi kernel...\n");
      gpuErrorCheckKernel("kern_calc_visi_point_or_gauss",
                           kern_calc_visi_point_or_gauss, grid, threads,
                           d_chunked_source->point_components, d_point_beam_gains,
                           d_us, d_vs, d_ws,
                           d_visibility_set->sum_visi_XX_real,
                           d_visibility_set->sum_visi_XX_imag,
                           d_visibility_set->sum_visi_XY_real,
                           d_visibility_set->sum_visi_XY_imag,
                           d_visibility_set->sum_visi_YX_real,
                           d_visibility_set->sum_visi_YX_imag,
                           d_visibility_set->sum_visi_YY_real,
                           d_visibility_set->sum_visi_YY_imag,
                           num_points, num_baselines, num_freqs, num_cross,
                           num_time_steps, beam_settings->beamtype, POINT);
      printf("\tVisi kernel done\n");

      free_d_components(d_chunked_source, POINT);
      free_extrapolated_flux_arrays(&d_chunked_source->point_components);
      free_beam_gains(d_point_beam_gains, beam_settings->beamtype);

    }//if point sources

    if (num_gauss > 0) {
      printf("\tDoing gaussian components\n");


      //Something to store the primary beam gains (all 4 pols) in
      d_beam_gains_t d_gauss_beam_gains;
      if (use_twobeams == 1) {
        d_gauss_beam_gains.d_ant1_to_baseline_map = d_ant1_to_baseline_map;
        d_gauss_beam_gains.d_ant2_to_baseline_map = d_ant2_to_baseline_map;
        d_gauss_beam_gains.use_twobeams = 1;
      } else {
        d_gauss_beam_gains.use_twobeams = 0;
      }

      source_component_common(woden_settings, beam_settings, d_freqs,
                              source, d_chunked_source,
                              &d_gauss_beam_gains, GAUSSIAN,
                              d_visibility_set);

      threads.x = 128;
      threads.y = 1;
      grid.x = (int)ceil( (float)num_cross / (float)threads.x );
      grid.y = 1;

      gpuErrorCheckKernel("kern_calc_visi_point_or_gauss",
                           kern_calc_visi_point_or_gauss, grid, threads,
                           d_chunked_source->gauss_components, d_gauss_beam_gains,
                           d_us, d_vs, d_ws,
                           d_visibility_set->sum_visi_XX_real,
                           d_visibility_set->sum_visi_XX_imag,
                           d_visibility_set->sum_visi_XY_real,
                           d_visibility_set->sum_visi_XY_imag,
                           d_visibility_set->sum_visi_YX_real,
                           d_visibility_set->sum_visi_YX_imag,
                           d_visibility_set->sum_visi_YY_real,
                           d_visibility_set->sum_visi_YY_imag,
                           num_gauss, num_baselines, num_freqs, num_cross,
                           num_time_steps, beam_settings->beamtype, GAUSSIAN);

      free_d_components(d_chunked_source, GAUSSIAN);
      free_extrapolated_flux_arrays(&d_chunked_source->gauss_components);
      free_beam_gains(d_gauss_beam_gains, beam_settings->beamtype);

    }//if gauss sources

    if (num_shapes > 0) {
      printf("\tDoing shapelet components\n");

      double *d_lsts=NULL;
      user_precision_t *d_u_shapes = NULL;
      user_precision_t *d_v_shapes = NULL;

      ( gpuMalloc( (void**)&(d_lsts),
                                          woden_settings->num_time_steps*sizeof(double)) );
      ( gpuMemcpy( d_lsts, woden_settings->lsts,
                             woden_settings->num_time_steps*sizeof(double),
                             gpuMemcpyHostToDevice) );

      ( gpuMalloc( (void**)&d_u_shapes,
            num_shapes*num_baselines*num_time_steps*sizeof(user_precision_t)) );
      ( gpuMalloc( (void**)&d_v_shapes,
            num_shapes*num_baselines*num_time_steps*sizeof(user_precision_t)) );

      //Something to store the primary beam gains (all 4 pols) in
      d_beam_gains_t d_shape_beam_gains;
      if (use_twobeams == 1) {
        d_shape_beam_gains.d_ant1_to_baseline_map = d_ant1_to_baseline_map;
        d_shape_beam_gains.d_ant2_to_baseline_map = d_ant2_to_baseline_map;
        d_shape_beam_gains.use_twobeams = 1;
      } else {
        d_shape_beam_gains.use_twobeams = 0;
      }

      source_component_common(woden_settings, beam_settings, d_freqs,
                              source, d_chunked_source,
                              &d_shape_beam_gains, SHAPELET,
                              d_visibility_set);


      if (num_shapes == 1) {
        threads.x = 128;
        threads.y = 1;
        grid.x = (int)ceil( (float)(num_baselines*num_time_steps) / (float)threads.x );
        grid.y = 1;
      }
      else {
        threads.x = 64;
        threads.y = 2;
        grid.x = (int)ceil( (float)(num_baselines*num_time_steps) / (float)threads.x );
        grid.y = (int)ceil( ((float)num_shapes) / ((float)threads.y) );
      }

      //Extra set of coords, centred on sky location of each shapelet source
      gpuErrorCheckKernel("kern_calc_uv_shapelet",
                            kern_calc_uv_shapelet, grid, threads,
                            d_X_diff, d_Y_diff, d_Z_diff,
                            d_u_shapes, d_v_shapes, d_lsts,
                            d_chunked_source->shape_components.ras,
                            d_chunked_source->shape_components.decs,
                            num_baselines, num_time_steps, num_shapes);

      //Splitting over visibilities, but looping over shapelet coeffs inside
      //the kernel
      //Splitting over visibilities, but looping over shapelet coeffs inside
      //the kernel
      threads.x = 128;
      threads.y = 1;
      grid.x = (int)ceil( (float)num_cross / (float)threads.x );
      grid.y = 1;

      gpuErrorCheckKernel("kern_calc_visi_shapelets",
            kern_calc_visi_shapelets, grid, threads,
            d_chunked_source->shape_components, d_shape_beam_gains,
            d_us, d_vs, d_ws,
            d_allsteps_wavelengths,
            d_u_shapes, d_v_shapes,
            d_visibility_set->sum_visi_XX_real,
            d_visibility_set->sum_visi_XX_imag,
            d_visibility_set->sum_visi_XY_real,
            d_visibility_set->sum_visi_XY_imag,
            d_visibility_set->sum_visi_YX_real,
            d_visibility_set->sum_visi_YX_imag,
            d_visibility_set->sum_visi_YY_real,
            d_visibility_set->sum_visi_YY_imag,
            d_sbf,
            num_shapes, num_baselines, num_freqs, num_cross,
            d_chunked_source->n_shape_coeffs, num_time_steps,
            beam_settings->beamtype);

      ( gpuFree(d_v_shapes) );
      ( gpuFree(d_u_shapes) );
      ( gpuFree(d_lsts) );

      printf("Making it to this call here\n");
      free_d_components(d_chunked_source, SHAPELET);
      free_extrapolated_flux_arrays(&d_chunked_source->shape_components);
      free_beam_gains(d_shape_beam_gains, beam_settings->beamtype);

    }//if shapelet

    ( gpuMemcpy(chunk_visibility_set->sum_visi_XX_real,
        d_visibility_set->sum_visi_XX_real, num_visis*sizeof(user_precision_t),
                                                     gpuMemcpyDeviceToHost) );
    ( gpuMemcpy(chunk_visibility_set->sum_visi_XX_imag,
        d_visibility_set->sum_visi_XX_imag, num_visis*sizeof(user_precision_t),
                                                     gpuMemcpyDeviceToHost) );

    ( gpuMemcpy(chunk_visibility_set->sum_visi_XY_real,
        d_visibility_set->sum_visi_XY_real, num_visis*sizeof(user_precision_t),
                                                     gpuMemcpyDeviceToHost) );
    ( gpuMemcpy(chunk_visibility_set->sum_visi_XY_imag,
        d_visibility_set->sum_visi_XY_imag, num_visis*sizeof(user_precision_t),
                                                     gpuMemcpyDeviceToHost) );

    ( gpuMemcpy(chunk_visibility_set->sum_visi_YX_real,
        d_visibility_set->sum_visi_YX_real, num_visis*sizeof(user_precision_t),
                                                     gpuMemcpyDeviceToHost) );
    ( gpuMemcpy(chunk_visibility_set->sum_visi_YX_imag,
        d_visibility_set->sum_visi_YX_imag, num_visis*sizeof(user_precision_t),
                                                     gpuMemcpyDeviceToHost) );

    ( gpuMemcpy(chunk_visibility_set->sum_visi_YY_real,
        d_visibility_set->sum_visi_YY_real, num_visis*sizeof(user_precision_t),
                                                     gpuMemcpyDeviceToHost) );
    ( gpuMemcpy(chunk_visibility_set->sum_visi_YY_imag,
        d_visibility_set->sum_visi_YY_imag, num_visis*sizeof(user_precision_t),
                                                     gpuMemcpyDeviceToHost) );

    ( gpuMemcpy(chunk_visibility_set->us_metres,
                                  d_u_metres,num_visis*sizeof(user_precision_t),
                                                      gpuMemcpyDeviceToHost) );
    ( gpuMemcpy(chunk_visibility_set->vs_metres,
                                  d_v_metres,num_visis*sizeof(user_precision_t),
                                                      gpuMemcpyDeviceToHost) );
    ( gpuMemcpy(chunk_visibility_set->ws_metres,
                                  d_w_metres,num_visis*sizeof(user_precision_t),
                                                      gpuMemcpyDeviceToHost) );

    //add to visiblity_set
    for (int visi = 0; visi < num_visis; visi++) {
      //if the first chunk then initialise our values, and copy across
      //the u,v,w coords
      if (chunk == 0) {
        //ensure temp visi's are 0.0
        // visibility_set->sum_visi_XX_real[visi] = 0;
        // visibility_set->sum_visi_XX_imag[visi] = 0;
        // visibility_set->sum_visi_XY_real[visi] = 0;
        // visibility_set->sum_visi_XY_imag[visi] = 0;
        // visibility_set->sum_visi_YX_real[visi] = 0;
        // visibility_set->sum_visi_YX_imag[visi] = 0;
        // visibility_set->sum_visi_YY_real[visi] = 0;
        // visibility_set->sum_visi_YY_imag[visi] = 0;

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
  ( gpuFree(d_freqs) );

  ( gpuFree(d_ws) );
  ( gpuFree(d_vs) );
  ( gpuFree(d_us) );
  ( gpuFree(d_w_metres) );
  ( gpuFree(d_v_metres) );
  ( gpuFree(d_u_metres) );
  ( gpuFree(d_allsteps_wavelengths) );
  ( gpuFree(d_allsteps_cha0s) );
  ( gpuFree(d_allsteps_sha0s) );
  ( gpuFree(d_Z_diff) );
  ( gpuFree(d_Y_diff) );
  ( gpuFree(d_X_diff) );

  if (cropped_sky_models->num_shapelets > 0) {
    ( gpuFree( d_sbf ) );
  }

  if (use_twobeams == 1) {
    ( gpuFree( d_ant1_to_baseline_map ) );
    ( gpuFree( d_ant2_to_baseline_map ) );
  }

  if (beam_settings->beamtype == FEE_BEAM || beam_settings->beamtype == FEE_BEAM_INTERP) {
    free_gpu_fee_beam(beam_settings->gpu_fee_beam);
    free(beam_settings->hyper_delays);
  }

  ( gpuFree(d_visibility_set->sum_visi_XX_imag) );
  ( gpuFree(d_visibility_set->sum_visi_XX_real) );
  ( gpuFree(d_visibility_set->sum_visi_XY_imag) );
  ( gpuFree(d_visibility_set->sum_visi_XY_real) );
  ( gpuFree(d_visibility_set->sum_visi_YX_imag) );
  ( gpuFree(d_visibility_set->sum_visi_YX_real) );
  ( gpuFree(d_visibility_set->sum_visi_YY_imag) );
  ( gpuFree(d_visibility_set->sum_visi_YY_real) );

}
