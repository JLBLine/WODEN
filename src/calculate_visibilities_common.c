#include "calculate_visibilities_common.h"


void calculate_component_visis(e_component_type comptype,
                               calc_visi_inouts_t *mem_calc_visi_inouts,
                               woden_settings_t *woden_settings,
                               beam_settings_t *beam_settings,
                               source_t *source, source_t *mem_chunked_source,
                               visibility_set_t *mem_visibility_set,
                               int num_beams, int use_twobeams,
                               int do_gpu) {

  int verbose = 0;

  int num_components;
  components_t mem_components;
  if (comptype == POINT){
    mem_components = mem_chunked_source->point_components;
    num_components = mem_chunked_source->n_points;
  } else if (comptype == GAUSSIAN){
    mem_components = mem_chunked_source->gauss_components;
    num_components = mem_chunked_source->n_gauss;
  } else {
    mem_components = mem_chunked_source->shape_components;
    num_components = mem_chunked_source->n_shapes;
  }

  //Something to store the primary beam gains (all 4 pols) in
  beam_gains_t *mem_beam_gains = malloc(sizeof(beam_gains_t));
  // beam_gains_t beam_gains;

  if (use_twobeams == 1) {
    mem_beam_gains->ant1_to_baseline_map = mem_calc_visi_inouts->ant1_to_baseline_map;
    mem_beam_gains->ant2_to_baseline_map = mem_calc_visi_inouts->ant2_to_baseline_map;
    mem_beam_gains->use_twobeams = 1;
  } else {
    mem_beam_gains->use_twobeams = 0;
  }
  //If an everybeam model, already calculated beam gains on the CPU
  //So just copy them across
  if (beam_settings->beamtype == EB_LOFAR || beam_settings->beamtype == EB_OSKAR  || beam_settings->beamtype == EB_MWA) {
    if (do_gpu == 1){
      copy_CPU_beam_gains_to_GPU(&mem_components, mem_beam_gains,
        num_components*woden_settings->num_freqs*woden_settings->num_time_steps*num_beams);
    } else {
      mem_beam_gains->gxs = mem_components.gxs;
      mem_beam_gains->Dxs = mem_components.Dxs;
      mem_beam_gains->Dys = mem_components.Dys;
      mem_beam_gains->gys = mem_components.gys;
    }
    
  }
  if (verbose == 1){
    printf("\tExtrapolating fluxes and beams...\n");
  }
  source_component_common(woden_settings, beam_settings,
                          mem_calc_visi_inouts->freqs,
                          source, mem_chunked_source,
                          mem_beam_gains, comptype,
                          mem_visibility_set);

  if (comptype == POINT){
    mem_components = mem_chunked_source->point_components;
  } else if (comptype == GAUSSIAN){
    mem_components = mem_chunked_source->gauss_components;
  } else {
    mem_components = mem_chunked_source->shape_components;
  }

  if (verbose == 1){
    printf("\tExtrapolating fluxes and beams done.\n");
    printf("\tDoing visi kernel...\n");
  }
  
  if (comptype == POINT || comptype == GAUSSIAN){
    if (do_gpu == 1){
      calc_visi_point_or_gauss_gpu(mem_components, *mem_beam_gains,
                                mem_calc_visi_inouts, mem_visibility_set,
                                num_components, beam_settings->beamtype,
                                comptype, woden_settings);
    } else {
      calc_visi_point_or_gauss_cpu(mem_components, *mem_beam_gains,
                                mem_calc_visi_inouts, mem_visibility_set,
                                num_components, beam_settings->beamtype,
                                comptype, woden_settings);
      // printf("CPU inside %.3f\n", mem_visibility_set->sum_visi_XX_real[0]);

    }
  } else {
    if (do_gpu == 1){
      calc_visi_shapelets_gpu(mem_components, *mem_beam_gains,
                              mem_calc_visi_inouts, mem_visibility_set,
                              num_components,
                              mem_chunked_source->n_shape_coeffs,
                              beam_settings->beamtype, woden_settings);
    } else {
      calc_visi_shapelets_cpu(mem_components, *mem_beam_gains,
                              mem_calc_visi_inouts, mem_visibility_set,
                              num_components,
                              mem_chunked_source->n_shape_coeffs,
                              beam_settings->beamtype, woden_settings);
    }
  } //end else it's a shapelet
  if (verbose == 1){
      printf("\tVisi kernel done\n");
    }
 
  if (do_gpu == 1){
    free_extrapolated_flux_arrays_gpu(&mem_components);
    free_components_gpu(mem_chunked_source, comptype);
    free_beam_gains_gpu(mem_beam_gains, beam_settings->beamtype);
  } else {

    free_components_cpu(&mem_components);
    free_beam_gains_cpu(mem_beam_gains, beam_settings->beamtype);
    free(mem_beam_gains);
  }
}


void calculate_visibilities(array_layout_t *array_layout,
  source_catalogue_t *cropped_sky_models, beam_settings_t *beam_settings,
  woden_settings_t *woden_settings, visibility_set_t *visibility_set,
  user_precision_t *sbf) {

    // printf("WE IS DOING THE THING\n");

  int verbose = 0;

  //Boolean for if we are using two beams per visibility or assuming all
  //beams are the same
  int use_twobeams = 0;
  int num_beams = 1;
  //Currently only do this is we have the use_dipamps flag; could be expanded
  if (woden_settings->use_dipamps == 1) {
    use_twobeams = 1;
    num_beams = woden_settings->num_ants;
  }

  int do_gpu = woden_settings->do_gpu;

  //Default behaviour with everybeam is to use a different beam for each station
  if (beam_settings->beamtype == EB_LOFAR || beam_settings->beamtype == EB_OSKAR  || beam_settings->beamtype == EB_MWA) {
    use_twobeams = 1;
    num_beams = woden_settings->num_ants;

    //However if single_everybeam_station is set, we're only using one beam
    //for all stations
    if (woden_settings->single_everybeam_station == 1) {
      use_twobeams = 0;
      num_beams = 1;
    }
  }

  const int num_visis = woden_settings->num_visis;
  const int num_cross = woden_settings->num_cross;

  //Setup chunk_visibility_set to hold the visibility outputs of each
  visibility_set_t *chunk_visibility_set = setup_visibility_set(num_visis);

  //Copy across settings for simulation - we always simulate all time steps
  //and frequncies for each chunk
  chunk_visibility_set->allsteps_sha0s = visibility_set->allsteps_sha0s;
  chunk_visibility_set->allsteps_cha0s = visibility_set->allsteps_cha0s;
  chunk_visibility_set->allsteps_lsts = visibility_set->allsteps_lsts;
  chunk_visibility_set->allsteps_wavelengths = visibility_set->allsteps_wavelengths;
  chunk_visibility_set->channel_frequencies = visibility_set->channel_frequencies;

  // Containers for GPU inputs and outputs. Putting everything in a struct
  // makes it easier to allocated and free GPU memory
  visibility_set_t *mem_visibility_set = malloc(sizeof(visibility_set_t));
  calc_visi_inouts_t *mem_calc_visi_inouts = malloc(sizeof(calc_visi_inouts_t));
  
  if (do_gpu == 1) {
    
    mem_calc_visi_inouts = create_calc_visi_inouts_gpu(array_layout,
                visibility_set, mem_visibility_set, sbf, woden_settings,
                cropped_sky_models->num_shapelets, use_twobeams);
  } else {
    mem_visibility_set = setup_visibility_set(num_visis);
    mem_calc_visi_inouts = create_calc_visi_inouts_cpu(array_layout,
                visibility_set, sbf, woden_settings,
                cropped_sky_models->num_shapelets, use_twobeams);
  }
  
  //TODO - this could be a function in some other file?
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

    if (do_gpu == 1) {
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
  }

  //Iterate through all sky model chunks, calculated visibilities are
  //added to chunk_visibility_set, and then summed onto visibility_set
  for (int chunk = 0; chunk < cropped_sky_models->num_sources; chunk++) {
    source_t *source = &cropped_sky_models->sources[chunk];

  //   // printf("\tsource->n_comps %d\n", source->n_comps);
  //   // printf("\tsource->n_points %d\n", source->n_points);
  //   // printf("\tsource->n_point_lists %d\n", source->n_point_lists);
  //   // printf("\tsource->n_point_powers %d\n", source->n_point_powers);
  //   // printf("\tsource->n_point_curves %d\n", source->n_point_curves);
  //   // printf("\tsource->n_gauss %d\n", source->n_gauss);
  //   // printf("\tsource->n_gauss_lists %d\n", source->n_gauss_lists);
  //   // printf("\tsource->n_gauss_powers %d\n", source->n_gauss_powers);
  //   // printf("\tsource->n_gauss_curves %d\n", source->n_gauss_curves);
  //   // printf("\tsource->n_shapes %d\n", source->n_shapes);
  //   // printf("\tsource->n_shape_lists %d\n", source->n_shape_lists);
  //   // printf("\tsource->n_shape_powers %d\n", source->n_shape_powers);
  //   // printf("\tsource->n_shape_curves %d\n", source->n_shape_curves);
  //   // printf("\tsource->n_shape_coeffs %d\n", source->n_shape_coeffs);

    source_t *mem_chunked_source;
    //If doing GPU, we have to copy contents across to GPU mem
    if (do_gpu == 1) {
      if (verbose == 1){
        printf("About to copy the chunked source to the GPU\n");
      }
      mem_chunked_source = copy_chunked_source_to_GPU(source);
      if (verbose == 1){
        printf("Have copied across the chunk to the GPU\n");
      }
    //If doing CPU, we can just use the source directly
    } else {
      mem_chunked_source = source;
    }

    //Make sure the temp visis are 0 at the start of each chunk
    // for (int visi = 0; visi < num_visis; visi++) {
    //   chunk_visibility_set->sum_visi_XX_real[visi] = 0;
    //   chunk_visibility_set->sum_visi_XX_imag[visi] = 0;
    //   chunk_visibility_set->sum_visi_XY_real[visi] = 0;
    //   chunk_visibility_set->sum_visi_XY_imag[visi] = 0;
    //   chunk_visibility_set->sum_visi_YX_real[visi] = 0;
    //   chunk_visibility_set->sum_visi_YX_imag[visi] = 0;
    //   chunk_visibility_set->sum_visi_YY_real[visi] = 0;
    //   chunk_visibility_set->sum_visi_YY_imag[visi] = 0;
    // }

    set_visi_set_to_zero_cpu(chunk_visibility_set, num_visis);

    //TODO setup some kind of logging and have an option to silence all this
    //printed output
    if (verbose == 1){
      printf("Processing chunk %d\n", chunk);
      printf("\tNumber of components in chunk are: P %d G %d S_coeffs %d\n",
                source->n_points,
                source->n_gauss,
                source->n_shape_coeffs );
    }

    if (do_gpu == 1) {
      set_visibilities_to_zero_gpu(mem_visibility_set, chunk_visibility_set, num_visis);

      calc_uvw_gpu(mem_calc_visi_inouts->X_diff, mem_calc_visi_inouts->Y_diff,
                   mem_calc_visi_inouts->Z_diff, mem_calc_visi_inouts->u_metres,
                   mem_calc_visi_inouts->v_metres, mem_calc_visi_inouts->w_metres,
                   mem_calc_visi_inouts->us, mem_calc_visi_inouts->vs,
                   mem_calc_visi_inouts->ws, mem_calc_visi_inouts->allsteps_wavelengths,
                   mem_calc_visi_inouts->allsteps_cha0s,
                   mem_calc_visi_inouts->allsteps_sha0s, woden_settings);

    } else {
      set_visi_set_to_zero_cpu(mem_visibility_set, num_visis);

      calc_uvw_cpu(array_layout->X_diff_metres, array_layout->Y_diff_metres,
                   array_layout->Z_diff_metres,
                   mem_calc_visi_inouts->u_metres, mem_calc_visi_inouts->v_metres,
                   mem_calc_visi_inouts->w_metres, mem_calc_visi_inouts->us,
                   mem_calc_visi_inouts->vs, mem_calc_visi_inouts->ws,
                   visibility_set->allsteps_wavelengths,
                   woden_settings->sdec0, woden_settings->cdec0,
                   visibility_set->allsteps_cha0s, visibility_set->allsteps_sha0s,
                   woden_settings->num_cross, woden_settings->num_baselines,
                   woden_settings->num_time_steps, woden_settings->num_freqs);

      // printf("WE IS DOING THE THING %.5f\n", mem_calc_visi_inouts->u_metres[num_visis-1]);
    }

    // printf("do_autos num_visis num_cross %d %d %d\n", woden_settings->do_autos, woden_settings->num_visis, woden_settings->num_cross);

    //Set auto-correlation uvw to zero if doing autos
    if (woden_settings->do_autos == 1){
      if (do_gpu == 1) {

        set_auto_uvw_to_zero_gpu(num_cross, num_visis,
                              mem_calc_visi_inouts->us, mem_calc_visi_inouts->vs,
                              mem_calc_visi_inouts->ws);

        set_auto_uvw_to_zero_gpu(num_cross, num_visis,
                              mem_calc_visi_inouts->u_metres,
                              mem_calc_visi_inouts->v_metres,
                              mem_calc_visi_inouts->w_metres);
      } else {
        for (int auto_ind = num_cross; auto_ind < num_visis; auto_ind++) {
          // printf("setting %d to zero out of %d %d\n", auto_ind, num_visis, num_cross);
          mem_calc_visi_inouts->us[auto_ind] = 0.0;
          mem_calc_visi_inouts->vs[auto_ind] = 0.0;
          mem_calc_visi_inouts->ws[auto_ind] = 0.0;

          mem_calc_visi_inouts->u_metres[auto_ind] = 0.0;
          mem_calc_visi_inouts->v_metres[auto_ind] = 0.0;
          mem_calc_visi_inouts->w_metres[auto_ind] = 0.0;
        }
      }
    }

    int num_points = source->n_points;
    int num_gauss = source->n_gauss;
    int num_shapes = source->n_shapes;

    if (num_points > 0) {
      if (verbose == 1){
        printf("\tDoing point components\n");
      }
      //The zero in this input is number of shaelet coefficients, as obviously
      //there are none for point sources
      calculate_component_visis(POINT, mem_calc_visi_inouts, woden_settings,
                               beam_settings, source, mem_chunked_source,
                               mem_visibility_set, num_beams, use_twobeams, do_gpu);

      // printf("CPU inside %.3f\n", mem_visibility_set->sum_visi_XX_real[0]);
      // printf("CPU inside %.3f\n", mem_visibility_set->sum_visi_XX_imag[0]);
      // printf("CPU inside %.3f\n", mem_visibility_set->sum_visi_XY_real[0]);
      // printf("CPU inside %.3f\n", mem_visibility_set->sum_visi_XY_imag[0]);
      // printf("CPU inside %.3f\n", mem_visibility_set->sum_visi_YX_real[0]);
      // printf("CPU inside %.3f\n", mem_visibility_set->sum_visi_YX_imag[0]);
      // printf("CPU inside %.3f\n", mem_visibility_set->sum_visi_YY_real[0]);
      // printf("CPU inside %.3f\n", mem_visibility_set->sum_visi_YY_imag[0]);
    }//if point sources

    if (num_gauss > 0) {
      if (verbose == 1){
        printf("\tDoing Gaussian components\n");
      }
      calculate_component_visis(GAUSSIAN, mem_calc_visi_inouts, woden_settings,
                               beam_settings, source, mem_chunked_source,
                               mem_visibility_set, num_beams, use_twobeams, do_gpu);
    }//if gauss sources

    if (num_shapes > 0) {
      if (verbose == 1){
        printf("\tDoing shapelet components\n");
      }

      if (do_gpu == 1) {
        calc_uv_shapelet_gpu(mem_calc_visi_inouts->u_shapes, mem_calc_visi_inouts->v_shapes,
                           num_shapes, mem_calc_visi_inouts->X_diff,
                           mem_calc_visi_inouts->Y_diff, mem_calc_visi_inouts->Z_diff,
                           mem_chunked_source->shape_components.ras,
                           mem_chunked_source->shape_components.decs,
                           woden_settings);
      } else {
        calc_uv_shapelet_cpu(array_layout->X_diff_metres, array_layout->Y_diff_metres,
                             array_layout->Z_diff_metres,
                             mem_calc_visi_inouts->u_shapes,
                             mem_calc_visi_inouts->v_shapes,
                             woden_settings->lsts,
                             mem_chunked_source->shape_components.ras,
                             mem_chunked_source->shape_components.decs,
                             woden_settings->num_baselines,
                             woden_settings->num_time_steps, num_shapes);
      }

      calculate_component_visis(SHAPELET, mem_calc_visi_inouts, woden_settings,
                                 beam_settings, source, mem_chunked_source,
                                 mem_visibility_set, num_beams, use_twobeams, do_gpu);
    }//if shapelet

    if (do_gpu == 1){
      copy_gpu_visi_set_to_host(mem_visibility_set, mem_calc_visi_inouts,
                                chunk_visibility_set, num_visis);
    } else {
      for (int visi = 0; visi < num_visis; visi++) {
        // // printf("wut %.1f\n", mem_calc_visi_inouts->w_metres[visi]);
        // printf("wut2 %.1f\n", chunk_visibility_set->ws_metres[visi]);
        chunk_visibility_set->us_metres[visi] = mem_calc_visi_inouts->u_metres[visi];
        chunk_visibility_set->vs_metres[visi] = mem_calc_visi_inouts->v_metres[visi];
        chunk_visibility_set->ws_metres[visi] = mem_calc_visi_inouts->w_metres[visi];
        chunk_visibility_set->sum_visi_XX_real[visi] = mem_visibility_set->sum_visi_XX_real[visi];
        chunk_visibility_set->sum_visi_XX_imag[visi] = mem_visibility_set->sum_visi_XX_imag[visi];
        chunk_visibility_set->sum_visi_XY_real[visi] = mem_visibility_set->sum_visi_XY_real[visi];
        chunk_visibility_set->sum_visi_XY_imag[visi] = mem_visibility_set->sum_visi_XY_imag[visi];
        chunk_visibility_set->sum_visi_YX_real[visi] = mem_visibility_set->sum_visi_YX_real[visi];
        chunk_visibility_set->sum_visi_YX_imag[visi] = mem_visibility_set->sum_visi_YX_imag[visi];
        chunk_visibility_set->sum_visi_YY_real[visi] = mem_visibility_set->sum_visi_YY_real[visi];
        chunk_visibility_set->sum_visi_YY_imag[visi] = mem_visibility_set->sum_visi_YY_imag[visi];
      }


      // memcpy(chunk_visibility_set->us_metres, mem_calc_visi_inouts->u_metres,
      //        num_visis*sizeof(user_precision_t));
      // memcpy(chunk_visibility_set->vs_metres, mem_calc_visi_inouts->v_metres,
      //         num_visis*sizeof(user_precision_t));
      // memcpy(chunk_visibility_set->ws_metres, mem_calc_visi_inouts->w_metres,
      //         num_visis*sizeof(user_precision_t));

      // chunk_visibility_set->us_metres = mem_calc_visi_inouts->u_metres;
      // chunk_visibility_set->vs_metres = mem_calc_visi_inouts->v_metres;
      // chunk_visibility_set->ws_metres = mem_calc_visi_inouts->w_metres;
    }

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

        // visibility_set->us_metres[visi] = mem_calc_visi_inouts->u_metres[visi];
        // visibility_set->vs_metres[visi] = mem_calc_visi_inouts->v_metres[visi];
        // visibility_set->ws_metres[visi] = mem_calc_visi_inouts->w_metres[visi];

        // printf("WE IS DOING THE THING %.5f\n", chunk_visibility_set->us_metres[visi]);
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
  if (do_gpu == 1) {
    free_calc_visi_inouts_gpu(mem_calc_visi_inouts, mem_visibility_set,
                              cropped_sky_models->num_shapelets, use_twobeams);
  } else {
    free_visi_set_outputs(mem_visibility_set);
    free( mem_visibility_set );
    free_calc_visi_inouts_cpu(mem_calc_visi_inouts, cropped_sky_models->num_shapelets,
                             use_twobeams);
  }
}
