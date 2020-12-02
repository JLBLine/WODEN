#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fitsio.h>
#include <erfa.h>
#include <complex.h>

#include "create_sky_model.h"
#include "shapelet_basis.h"
#include "woden.h"
#include "constants.h"
#include "chunk_source.h"
#include "print_help.h"
#include "primary_beam.h"
#include "FEE_primary_beam.h"

// #include "FEE_primary_beam_cuda.h"

int main(int argc, char **argv) {

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

  //Create the shapelet basis function array
  float *sbf;
  sbf = NULL;
  sbf = malloc( sbf_N * sbf_L * sizeof(float) );
  sbf = create_sbf(sbf);

  //Read in the settings from the controlling json file
  woden_settings_t *woden_settings;
  woden_settings = read_json_settings(argv[1]);

  if (woden_settings->chunking_size > MAX_CHUNKING_SIZE) {
    printf("Current maximum allowable chunk size is %d.  Defaulting to this value.", MAX_CHUNKING_SIZE);
    woden_settings->chunking_size = MAX_CHUNKING_SIZE;
  }
  else if (woden_settings->chunking_size < 1 ) {
    woden_settings->chunking_size = MAX_CHUNKING_SIZE;
  }

  //Create the array layout in instrument-centric X,Y,Z using positions
  //Rotate back to J2000 if necessary
  array_layout_t * array_layout;
  array_layout = calc_XYZ_diffs(woden_settings);
  woden_settings->num_baselines = array_layout->num_baselines;

  //Set some constants based on the settings
  float ha0, sha0, cha0;
  float wavelength;
  float frequency;
  const int num_visis = woden_settings->num_baselines * woden_settings->num_time_steps * woden_settings->num_freqs;

  woden_settings->num_visis = num_visis;

  float sdec0,cdec0;
  sdec0 = sin(woden_settings->dec0); cdec0=cos(woden_settings->dec0);

  printf("Setting phase centre RA,DEC %.5fdeg %.5fdeg\n",woden_settings->ra0/DD2R, woden_settings->dec0/DD2R);

  //Used for calculating l,m,n for components
  float angles_array[3] = {sdec0, cdec0, woden_settings->ra0};
  int num_time_steps = woden_settings->num_time_steps;

  //Calculate all lsts for this observation
  float lsts[num_time_steps];

  for ( int time_step = 0; time_step < woden_settings->num_time_steps; time_step++ ) {
    float lst = woden_settings->lst_base + time_step*woden_settings->time_res*SOLAR2SIDEREAL*DS2R;

    //TODO add half a time step is good? Add time decorrelation?
    lst += 0.5*woden_settings->time_res*SOLAR2SIDEREAL*DS2R;
    lsts[time_step] = lst;
  }

  //Read in the source catalogue
  source_catalogue_t *raw_srccat;
  raw_srccat = read_source_catalogue(woden_settings->cat_filename);

  //Crop emission below the horizon, and collapse all SOURCES from raw_srccat
  //into one single SOURCE
  printf("Horizon cropping sky model and calculating az/za for all components for observation\n");
  catsource_t *cropped_src;
  cropped_src = crop_sky_model(raw_srccat, lsts, num_time_steps, woden_settings->sky_crop_type);

  printf("Finished cropping and calculating az/za\n");

  //Setup some beam settings given user chose parameters
  beam_settings_t beam_settings;
  beam_settings = fill_primary_beam_settings(woden_settings, cropped_src,
                                            lsts, num_time_steps);

  //TODO allow this frequency resolution to over-written by user commands
  //MWA correlator data is split into 24 'coarse' bands of 1.28MHz bandwidth,
  //which is typically split into 10, 20, or 40kHz fine channels
  //Loop through each coarse frequency band, run the simulation and dump to
  //a binary file
  for (size_t band = 0; band < woden_settings->num_bands; band++) {
    //Set the lower frequency edge for this coarse band
    int band_num = woden_settings->band_nums[band];
    float base_band_freq = ((band_num - 1)*woden_settings->coarse_band_width) + woden_settings->base_low_freq;
    printf("Simulating band %02d with bottom freq %.8e\n",band_num,base_band_freq);

    woden_settings->base_band_freq = base_band_freq;

    //TODO - add half a freq resolution in here? minus? Leave as is?
    // base_band_freq += woden_settings->frequency_resolution/2.0;



      beam_settings.FEE_beam = malloc(sizeof(copy_primary_beam_t));
      // //We need the zenith beam to get the normalisation
      beam_settings.FEE_beam_zenith = malloc(sizeof(copy_primary_beam_t));

      //The intial setup of the FEE beam is done on the CPU, so call it here
      if (woden_settings->beamtype == FEE_BEAM){
        float base_middle_freq = base_band_freq + woden_settings->coarse_band_width/2.0;
      //
        // Just use one single tile beam for all for now - will need a certain
        // number in the future to include dipole flagging
        int st = 0;
        printf("Middle freq is %f\n",base_middle_freq );
      //
        float float_zenith_delays[16] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
      //
        printf("Setting up the zenith FEE beam...");
        RTS_HDFBeamInit(woden_settings->hdf5_beam_path, base_middle_freq, beam_settings.FEE_beam_zenith, float_zenith_delays, st);
        printf(" done.\n");

        printf("Setting up the FEE beam...");
        RTS_HDFBeamInit(woden_settings->hdf5_beam_path, base_middle_freq,
              beam_settings.FEE_beam, woden_settings->FEE_ideal_delays, st);
        printf(" done.\n");

      }

    visibility_set_t *visibility_set = malloc(sizeof(visibility_set_t));
    visibility_set->sum_visi_real = malloc( num_visis * sizeof(float) );
    visibility_set->sum_visi_imag = malloc( num_visis * sizeof(float) );
    visibility_set->us_metres = malloc( num_visis * sizeof(float) );
    visibility_set->vs_metres = malloc( num_visis * sizeof(float) );
    visibility_set->ws_metres = malloc( num_visis * sizeof(float) );
    visibility_set->sha0s = malloc( num_visis * sizeof(float) );
    visibility_set->cha0s = malloc( num_visis * sizeof(float) );
    visibility_set->lsts = malloc( num_visis * sizeof(float) );
    visibility_set->wavelengths = malloc( num_visis * sizeof(float) );
    visibility_set->channel_frequencies = malloc( (int)woden_settings->num_freqs * sizeof(float) );

    visibility_set->sum_visi_XX_real = malloc( num_visis * sizeof(float) );
    visibility_set->sum_visi_XX_imag = malloc( num_visis * sizeof(float) );
    visibility_set->sum_visi_XY_real = malloc( num_visis * sizeof(float) );
    visibility_set->sum_visi_XY_imag = malloc( num_visis * sizeof(float) );
    visibility_set->sum_visi_YX_real = malloc( num_visis * sizeof(float) );
    visibility_set->sum_visi_YX_imag = malloc( num_visis * sizeof(float) );
    visibility_set->sum_visi_YY_real = malloc( num_visis * sizeof(float) );
    visibility_set->sum_visi_YY_imag = malloc( num_visis * sizeof(float) );

    visibility_set->sum_visi_real = malloc( num_visis * sizeof(float) );
    visibility_set->sum_visi_imag = malloc( num_visis * sizeof(float) );

    //Fill in the fine channel frequencies
    for (int freq_step = 0; freq_step < woden_settings->num_freqs; freq_step++) {
      frequency = base_band_freq + (woden_settings->frequency_resolution*freq_step);
      visibility_set->channel_frequencies[freq_step] = frequency;
    }

    //Fill in visibility settings in order of baseline,freq,time
    //Order matches that of a uvfits file (I live in the past)

    for ( int time_step = 0; time_step < woden_settings->num_time_steps; time_step++ ) {
      ha0 = lsts[time_step] - woden_settings->ra0;
      sha0 = sin(ha0); cha0=cos(ha0);

      for (int freq_step = 0; freq_step < woden_settings->num_freqs; freq_step++) {
        frequency = base_band_freq + (woden_settings->frequency_resolution*freq_step);
        wavelength = VELC / frequency;
        int step = woden_settings->num_baselines*(time_step*woden_settings->num_freqs + freq_step);

        for (int baseline = 0; baseline < woden_settings->num_baselines; baseline++) {
          visibility_set->cha0s[step + baseline] = cha0;
          visibility_set->sha0s[step + baseline] = sha0;
          visibility_set->lsts[step + baseline] = lsts[time_step];
          visibility_set->wavelengths[step + baseline] = wavelength;
        }//baseline loop
      }//freq loop
    }//time loop

    //Calculating a single shapelet coeff is equivalent to a point/gauss so treat as a
    //component here
    int num_components = cropped_src->n_points + cropped_src->n_gauss + cropped_src->n_shape_coeffs;
    int num_chunks;
    //TODO should we chunk outside the band for-loop so that we can reuse the chunks for each band (should be the same)
    if (num_components > woden_settings->chunking_size) {
      printf("Chunking sky model\n");
      num_chunks = num_components / woden_settings->chunking_size;

      if (num_components % woden_settings->chunking_size != 0) {
        num_chunks = num_chunks + 1;
      }
      printf("Number of chunks required is %d\n",  num_chunks);
    } else {
      num_chunks = 1;
    }

    //Make a struct to contain sky as many sky models and beam settings as
    //needed
    source_catalogue_t *cropped_sky_models;
    cropped_sky_models = malloc(sizeof(source_catalogue_t));

    cropped_sky_models->num_sources = num_chunks;
    cropped_sky_models->num_shapelets = 0;
    cropped_sky_models->catsources = malloc(num_chunks*sizeof(catsource_t));
    cropped_sky_models->beam_settings = malloc(num_chunks*sizeof(beam_settings_t));

    catsource_t *temp_cropped_src = malloc(sizeof(catsource_t));
    int point_iter = 0;
    int gauss_iter = 0;
    int shape_iter = 0;

    printf("Chunking sky model..\n");
    //For each chunk, calculate the visibilities for those components
    for (int chunk = 0; chunk < num_chunks; chunk++) {

      fill_chunk_src(temp_cropped_src, cropped_src, num_chunks, chunk,
                     woden_settings->chunking_size, woden_settings->num_time_steps,
                     &point_iter, &gauss_iter, &shape_iter);

      //Add the number of shapelets onto the full source catalogue value
      //so we know if we need to setup shapelet basis functions in GPU memory
      //or not

      cropped_sky_models->num_shapelets += temp_cropped_src->n_shapes;

      beam_settings_t beam_settings_chunk;
      beam_settings_chunk = make_beam_settings_chunk(beam_settings, temp_cropped_src,
                            cropped_src, woden_settings, point_iter, gauss_iter, shape_iter);

      // printf("Managed to chunk the beam %d\n",chunk );

      cropped_sky_models->catsources[chunk] = *temp_cropped_src;
      cropped_sky_models->beam_settings[chunk] = beam_settings_chunk;

    }

    printf("Sky model chunked.\n");

    //setup a temporary visibility set that calculate_visibilities will populate
    visibility_set_t *chunk_visibility_set = malloc(sizeof(visibility_set_t));

    chunk_visibility_set->us_metres = malloc( num_visis * sizeof(float) );
    chunk_visibility_set->vs_metres = malloc( num_visis * sizeof(float) );
    chunk_visibility_set->ws_metres = malloc( num_visis * sizeof(float) );

    chunk_visibility_set->sha0s = visibility_set->sha0s;
    chunk_visibility_set->cha0s = visibility_set->cha0s;
    chunk_visibility_set->lsts = visibility_set->lsts;
    chunk_visibility_set->wavelengths = visibility_set->wavelengths;
    chunk_visibility_set->channel_frequencies = visibility_set->channel_frequencies;

    chunk_visibility_set->sum_visi_XX_real = malloc( num_visis * sizeof(float) );
    chunk_visibility_set->sum_visi_XX_imag = malloc( num_visis * sizeof(float) );
    chunk_visibility_set->sum_visi_XY_real = malloc( num_visis * sizeof(float) );
    chunk_visibility_set->sum_visi_XY_imag = malloc( num_visis * sizeof(float) );
    chunk_visibility_set->sum_visi_YX_real = malloc( num_visis * sizeof(float) );
    chunk_visibility_set->sum_visi_YX_imag = malloc( num_visis * sizeof(float) );
    chunk_visibility_set->sum_visi_YY_real = malloc( num_visis * sizeof(float) );
    chunk_visibility_set->sum_visi_YY_imag = malloc( num_visis * sizeof(float) );

    calculate_visibilities(array_layout, cropped_sky_models,
                  angles_array, woden_settings,
                  visibility_set, chunk_visibility_set, sbf,
                  num_chunks);

    printf("GPU calls for band %d finished\n",band_num );

    free(temp_cropped_src);

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

    // if (woden_settings->beamtype == FEE_BEAM) {
    //   free_FEE_primary_beam_from_GPU(beam_settings.FEE_beam);
    // }

    //Dumps u,v,w (metres), Re(vis), Im(vis) to a binary file
    FILE *output_visi;
    char buf[0x100];
    snprintf(buf, sizeof(buf), "output_visi_band%02d.dat", band_num);

    output_visi = fopen(buf,"wb");

    if(output_visi == NULL)
    {
        printf("Could not open output_visi_band%02d.dat - exiting", band_num);
        exit(1);
    }

    fwrite(visibility_set->us_metres, num_visis*sizeof(float), 1, output_visi);
    fwrite(visibility_set->vs_metres, num_visis*sizeof(float), 1, output_visi);
    fwrite(visibility_set->ws_metres, num_visis*sizeof(float), 1, output_visi);
    // fwrite(visibility_set->sum_visi_real, num_visis*sizeof(float), 1, output_visi);
    // fwrite(visibility_set->sum_visi_imag, num_visis*sizeof(float), 1, output_visi);

    fwrite(visibility_set->sum_visi_XX_real, num_visis*sizeof(float), 1, output_visi);
    fwrite(visibility_set->sum_visi_XX_imag, num_visis*sizeof(float), 1, output_visi);
    fwrite(visibility_set->sum_visi_XY_real, num_visis*sizeof(float), 1, output_visi);
    fwrite(visibility_set->sum_visi_XY_imag, num_visis*sizeof(float), 1, output_visi);
    fwrite(visibility_set->sum_visi_YX_real, num_visis*sizeof(float), 1, output_visi);
    fwrite(visibility_set->sum_visi_YX_imag, num_visis*sizeof(float), 1, output_visi);
    fwrite(visibility_set->sum_visi_YY_real, num_visis*sizeof(float), 1, output_visi);
    fwrite(visibility_set->sum_visi_YY_imag, num_visis*sizeof(float), 1, output_visi);

    fflush(output_visi);
    fclose(output_visi);


    // // Dumps u,v,w (metres), Re(vis), Im(vis) directly to text file - useful for
    // // bug hunting with small outputs
    // FILE *output_visi_text;
    // char buff[0x100];
    // snprintf(buff, sizeof(buff), "output_visi_band%02d.txt", band_num);
    // output_visi_text = fopen(buff,"w");
    // for ( int time_step = 0; time_step < woden_settings->num_time_steps; time_step++ ) {
    //   for ( int freq_step = 0; freq_step < woden_settings->num_freqs; freq_step++ ) {
    //     for (int baseline = 0; baseline < woden_settings->num_baselines; baseline++) {
    //       int step = woden_settings->num_baselines*(time_step*woden_settings->num_freqs + freq_step);
    //       fprintf(output_visi_text,"%f %f %f %f %f\n",visibility_set->us_metres[step + baseline],
    //               visibility_set->vs_metres[step + baseline],visibility_set->ws_metres[step + baseline],
    //               visibility_set->sum_visi_XX_real[step + baseline],visibility_set->sum_visi_XX_imag[step + baseline]);
    //     }
    //   }
    // }
    // fflush(output_visi_text);
    // fclose(output_visi_text);
    //

    //Free up that memory
    free( visibility_set->us_metres );
    free( visibility_set->vs_metres );
    free( visibility_set->ws_metres );
    free( visibility_set->sha0s );
    free( visibility_set->cha0s );
    free( visibility_set->lsts );
    free( visibility_set->wavelengths );

    free(visibility_set->sum_visi_XX_real);
    free(visibility_set->sum_visi_XX_imag);
    free(visibility_set->sum_visi_XY_real);
    free(visibility_set->sum_visi_XY_imag);
    free(visibility_set->sum_visi_YX_real);
    free(visibility_set->sum_visi_YX_imag);
    free(visibility_set->sum_visi_YY_real);
    free(visibility_set->sum_visi_YY_imag);

    free( visibility_set );

  }//band loop
  printf("WODEN is done\n");
}//main
