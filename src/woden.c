#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <erfa.h>
#include <complex.h>

#include "woden_precision_defs.h"
#include "constants.h"
#include "woden_struct_defs.h"

#include "create_sky_model.h"
#include "shapelet_basis.h"
#include "chunk_sky_model.h"
#include "print_help.h"
#include "primary_beam.h"
#include "FEE_primary_beam.h"
#include "woden_settings.h"
#include "visibility_set.h"
#include "array_layout.h"

//Main CUDA executable to link in
extern void calculate_visibilities(array_layout_t * array_layout,
  source_catalogue_t *cropped_sky_models, beam_settings_t *beam_settings,
  woden_settings_t *woden_settings,  visibility_set_t *visibility_set,
  user_precision_t *sbf);

int main(int argc, char **argv) {

  #ifdef DOUBLE_PRECISION
  printf("WODEN is using DOUBLE precision\n");
  #else
  printf("WODEN is using FLOAT precision\n");
  #endif

  //If not enough arguments, print help
  if (argc < 2) {
    print_cmdline_help();
    exit(1);
  }

  //If --help is passed, print help
  if (strcmp("--help", argv[1]) == 0) {
    print_cmdline_help();
    exit(1);

  }

  //If -h is passed, print help
  if (strcmp("-h", argv[1]) == 0) {
    print_cmdline_help();
    exit(1);

  }

  //Is everything OK integer
  int status=0;

  //Create the shapelet basis function array
  user_precision_t *sbf = malloc( sbf_N * sbf_L * sizeof(user_precision_t) );
  sbf = create_sbf(sbf);

  //Read in the settings from the controlling json file
  woden_settings_t *woden_settings = malloc( sizeof(woden_settings_t) );
  status = read_json_settings(argv[1], woden_settings);

  if (status == 1) {
    printf("read_json_settings failed. Exiting now\n");
    exit(1);
  }

  //Create the array layout in instrument-centric X,Y,Z using positions
  //Rotate back to J2000 if necessary
  array_layout_t * array_layout;
  array_layout = calc_XYZ_diffs(woden_settings, woden_settings->do_precession);

  //Setup all LSTs array for all time steps in this simulation
  double *lsts = setup_lsts_and_phase_centre(woden_settings);

  //Read in the source catalogue
  source_catalogue_t *raw_srccat = malloc( sizeof(source_catalogue_t) );
  status = read_source_catalogue(woden_settings->cat_filename, raw_srccat);

  if (status == 1) {
    printf("read_source_catalogue failed. Exiting now\n");
    exit(1);
  }

  //Crop emission below the horizon, and collapse all SOURCES from raw_srccat
  //into one single SOURCE
  printf("Horizon cropping sky model and calculating az/za for all components\nfor observation\n");
  // catsource_t *cropped_src;
  catsource_t *cropped_src = crop_sky_model(raw_srccat, lsts, woden_settings->latitude,
                               woden_settings->num_time_steps, woden_settings->sky_crop_type);

  printf("Finished cropping and calculating az/za\n");

  //Setup some beam settings given user chose parameters
  beam_settings_t *beam_settings = fill_primary_beam_settings(woden_settings,
                                                             cropped_src, lsts);

  // Chunk the sky models into smaller pieces that fit onto the GPU
  source_catalogue_t *cropped_sky_models = create_chunked_sky_models(cropped_src,
                                                                     woden_settings);


  //MWA correlator data is split into 24 'coarse' bands of 1.28MHz bandwidth,
  //which is typically split into 10, 20, or 40kHz fine channels
  //User can change these settings using run_woden.py / in the json
  //Loop through each coarse frequency band, run the simulation and dump to
  //a binary file
  for (size_t band = 0; band < woden_settings->num_bands; band++) {
    //Set the lower frequency edge for this coarse band
    int band_num = woden_settings->band_nums[band];
    user_precision_t base_band_freq = ((band_num - 1)*woden_settings->coarse_band_width) + woden_settings->base_low_freq;
    printf("Simulating band %02d with bottom freq %.8e\n",band_num,base_band_freq);

    woden_settings->base_band_freq = base_band_freq;

    beam_settings->FEE_beam = malloc(sizeof(RTS_MWA_FEE_beam_t));
    //We need the zenith beam to get the normalisation
    beam_settings->FEE_beam_zenith = malloc(sizeof(RTS_MWA_FEE_beam_t));

    //The intial setup of the FEE beam is done on the CPU, so call it here
    if (woden_settings->beamtype == FEE_BEAM){
      user_precision_t base_middle_freq = base_band_freq + (woden_settings->coarse_band_width/2.0);
    //
      printf("Middle freq is %.8e \n",base_middle_freq );
    //
      user_precision_t user_precision_t_zenith_delays[16] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    //
      printf("Setting up the zenith FEE beam...");
      status = RTS_MWAFEEInit(woden_settings->hdf5_beam_path, base_middle_freq, beam_settings->FEE_beam_zenith, user_precision_t_zenith_delays);
      printf(" done.\n");

      printf("Setting up the FEE beam...");
      status = RTS_MWAFEEInit(woden_settings->hdf5_beam_path, base_middle_freq,
            beam_settings->FEE_beam, woden_settings->FEE_ideal_delays);
      printf(" done.\n");

    }

    //Setup the visibility container
    visibility_set_t *visibility_set = setup_visibility_set(woden_settings->num_visis);

    //Fill in the time/freq/baseline settings in `visiblity_set` needed by
    //calculate_visibilities
    fill_timefreq_visibility_set(visibility_set, woden_settings,
                                 base_band_freq, lsts);

    //Launch the CUDA code
    calculate_visibilities(array_layout, cropped_sky_models, beam_settings,
                  woden_settings, visibility_set, sbf);

    printf("GPU calls for band %d finished\n",band_num );

    //Write out binary file for python code to read and convert to uvfits
    write_visi_set_binary(visibility_set, band_num, woden_settings->num_visis);

    //Writes out a text file with u,v,w XX_re, XX_im. Useful for bug hunting
    //in desperation
    // write_visi_set_text(visibility_set, band_num, woden_settings);

    //Free up that memory
    free_visi_set_inputs(visibility_set);
    free_visi_set_outputs(visibility_set);
    free( visibility_set );

    //Release the CPU MWA FEE beam if required
    if (woden_settings->beamtype == FEE_BEAM){
      RTS_freeHDFBeam(beam_settings->FEE_beam);
      RTS_freeHDFBeam(beam_settings->FEE_beam_zenith);
    }

  }//band loop
  printf("WODEN is done\n");
}//main
