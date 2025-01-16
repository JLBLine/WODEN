#include "calculate_visibilities_common.h"


void calculate_component_visis(e_component_type comptype,
                               calc_visi_inouts_t *mem_calc_visi_inouts,
                               woden_settings_t *woden_settings,
                               beam_settings_t *beam_settings,
                               source_t *source, source_t *mem_chunked_source,
                               visibility_set_t *mem_visibility_set,
                               int num_beams, int use_twobeams,
                               int do_gpu) {

  // char log_buffer[128];
  // int log_len = sizeof log_buffer;

  if (woden_settings->verbose == 1){
    log_message("Starting calculate_component_visis");
  }

  int num_components;
  components_t mem_components;
  // components_t *components;
  if (comptype == POINT){
    mem_components = mem_chunked_source->point_components;
    // components = &source->point_components;
    num_components = mem_chunked_source->n_points;
  } else if (comptype == GAUSSIAN){
    mem_components = mem_chunked_source->gauss_components;
    // components = &source->gauss_components;
    num_components = mem_chunked_source->n_gauss;
  } else {
    mem_components = mem_chunked_source->shape_components;
    // components = &source->shape_components;
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

  if (woden_settings->verbose == 1){
    log_message("\tExtrapolating fluxes and beams...");
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

  if (woden_settings->verbose == 1){
    log_message("\tExtrapolating fluxes and beams done.");
    log_message("\tDoing visi kernel...");
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
  if (woden_settings->verbose == 1){
      log_message("\tVisi kernel done");
    }
 
  if (do_gpu == 1){
    free_extrapolated_flux_arrays_gpu(&mem_components);
    free_components_gpu(mem_chunked_source, comptype);
    free_beam_gains_gpu(mem_beam_gains, beam_settings->beamtype);
    free(mem_beam_gains);
  } else {

    free_components_cpu(&mem_components);
    free_beam_gains_cpu(mem_beam_gains, beam_settings->beamtype);
    free(mem_beam_gains);
  }
  if (woden_settings->verbose == 1){
    log_message("Finished calculate_component_visis");
  }
}


void calculate_visibilities(array_layout_t *array_layout,
  source_catalogue_t *cropped_sky_models, beam_settings_t *beam_settings,
  woden_settings_t *woden_settings, visibility_set_t *visibility_set,
  user_precision_t *sbf) {

  char log_buffer[128];
  int log_len = sizeof log_buffer;

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

    if (use_twobeams == 0) {
      num_amps = 16;

    } else {
      num_amps = 32;
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
        log_message("WARNING - Something went wrong launching new_gpu_fee_beam");
      }
    }
  }

  //Iterate through all sky model chunks, calculated visibilities are
  //added to chunk_visibility_set, and then summed onto visibility_set
  for (int chunk = 0; chunk < cropped_sky_models->num_sources; chunk++) {
    source_t *source = &cropped_sky_models->sources[chunk];

    source_t *mem_chunked_source;
    //If doing GPU, we have to copy contents across to GPU mem
    if (do_gpu == 1) {
      if (woden_settings->verbose == 1){
        log_message("About to copy the chunked source to the GPU");
      }
      mem_chunked_source = copy_chunked_source_to_GPU(source);
      if (woden_settings->verbose == 1){
        log_message("Have copied across the chunk to the GPU");
      }
    //If doing CPU, we can just use the source directly
    } else {
      mem_chunked_source = source;
    }

    set_visi_set_to_zero_cpu(chunk_visibility_set, num_visis);

    //TODO setup some kind of logging and have an option to silence all this
    //printed output
    if (woden_settings->verbose == 1){
      snprintf(log_buffer, log_len, "Processing chunk %d", chunk);
      log_message(log_buffer);

      snprintf(log_buffer, log_len,
              "\tNumber of components in chunk are: P %d G %d S_coeffs %d",
              source->n_points,
              source->n_gauss,
              source->n_shape_coeffs );
      log_message(log_buffer);
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
      if (woden_settings->verbose == 1){
        log_message("\tDoing point components");
      }
      //The zero in this input is number of shaelet coefficients, as obviously
      //there are none for point sources
      calculate_component_visis(POINT, mem_calc_visi_inouts, woden_settings,
                               beam_settings, source, mem_chunked_source,
                               mem_visibility_set, num_beams, use_twobeams, do_gpu);

    }//if point sources

    if (num_gauss > 0) {
      if (woden_settings->verbose == 1){
        log_message("\tDoing Gaussian components");
      }
      calculate_component_visis(GAUSSIAN, mem_calc_visi_inouts, woden_settings,
                               beam_settings, source, mem_chunked_source,
                               mem_visibility_set, num_beams, use_twobeams, do_gpu);
    }//if gauss sources

    if (num_shapes > 0) {
      if (woden_settings->verbose == 1){
        log_message("\tDoing shapelet components");
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
    }

    //add to visiblity_set
    for (int visi = 0; visi < num_visis; visi++) {
      //if the first chunk then initialise our values, and copy across
      //the u,v,w coords
      if (chunk == 0) {
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
