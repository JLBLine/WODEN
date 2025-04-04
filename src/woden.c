#include "woden.h"

int run_woden(woden_settings_t *woden_settings, visibility_set_t *visibility_sets,
             source_catalogue_t *cropped_sky_models, array_layout_t * array_layout,
             user_precision_t *sbf) {
  

  // if (woden_settings->do_gpu) {
  //   log_message("Will be using the GPU");
  // } else {
  //   log_message("Will be using only the CPU");
  // }

  //Setup a buffer on the stack to hold log messages
  //Things made on the stack don't have to be freed
  char log_buffer[128];
  int log_len = sizeof log_buffer;

  if (woden_settings->verbose == 1){
    #ifdef DOUBLE_PRECISION
    // printf("WODEN is using DOUBLE precision\n");
    log_message("WODEN is using DOUBLE precision");
    #else
    // printf("WODEN is using FLOAT precision\n");
    log_message("WODEN is using FLOAT precision");
    #endif
  }
  

  //Is everything OK integer
  int status=0;

  //Setup some beam settings given user chose parameters
  beam_settings_t *beam_settings = fill_primary_beam_settings(woden_settings,
                                                              woden_settings->lsts);

  //MWA correlator data is split into 24 'coarse' bands of 1.28MHz bandwidth,
  //which is typically split into 10, 20, or 40kHz fine channels
  //User can change these settings using run_woden.py / in the json
  //Loop through each coarse frequency band, run the simulation and dump to
  //a binary file
  for (int band = 0; band < woden_settings->num_bands; band++) {
    //Set the lower frequency edge for this coarse band
    int band_num = woden_settings->band_nums[band];
    double base_band_freq = ((band_num - 1)*woden_settings->coarse_band_width) + woden_settings->base_low_freq;

    if (woden_settings->verbose == 1){
      snprintf(log_buffer, log_len,
              "Simulating band %02d with bottom freq %.8e", band_num, base_band_freq);
      log_message(log_buffer);
    }

    woden_settings->base_band_freq = base_band_freq;

    //Fill in the time/freq/baseline settings in `visiblity_set` needed by
    //calculate_visibilities
    fill_timefreq_visibility_set(&visibility_sets[band], woden_settings,
                                  base_band_freq, woden_settings->lsts);

    //The intial setup of the FEE beam is done on the CPU, so call it here
    if (woden_settings->beamtype == FEE_BEAM || woden_settings->beamtype == FEE_BEAM_INTERP){
      double base_middle_freq = base_band_freq + (woden_settings->coarse_band_width/2.0);
      beam_settings->base_middle_freq = base_middle_freq;

      if (woden_settings->verbose == 1){
        snprintf(log_buffer, log_len,
                "Middle freq for FEE beam is %.8e",base_middle_freq);
        log_message(log_buffer);
      }

      status = new_fee_beam(woden_settings->hdf5_beam_path,
                            &beam_settings->fee_beam);
      if (status != 0) {
        handle_hyperbeam_error(__FILE__, __LINE__, "new_fee_beam");
        // printf("Something went wrong launching new_fee_beam\n");
      }
    }

    //Launch the GPU or CPU code
    calculate_visibilities(array_layout, cropped_sky_models, beam_settings,
                  woden_settings, &visibility_sets[band], sbf);

    if (woden_settings->verbose == 1){
      snprintf(log_buffer, log_len, "Calls for band %d finished",band_num );
        log_message(log_buffer);
    }

    //Release the CPU MWA FEE beam if required
    if (woden_settings->beamtype == FEE_BEAM || woden_settings->beamtype == FEE_BEAM_INTERP){
      free_fee_beam(beam_settings->fee_beam);
    }
  }//band loop
  // log_message("run_woden.c is ending");
  return status;
}//run_woden
