#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fitsio.h>
#include <erfa.h>

#include "create_sky_model.h"
#include "shapelet_basis.h"
#include "woden.h"
#include "constants.h"
#include "chunk_source.h"

int main(int argc, char **argv) {

  if (argc < 2) {
    printf("Must input a .json settings file to run woden. This file must include:\n");
    printf("\tra0: ra phase centre (float in degrees)\n");
    printf("\tdec0: dec phase centre (float in degrees)\n");
    printf("\tnum_freqs: number of fine frequency channels to simulate (int)\n");
    printf("\tnum_time_steps: number of time steps to simulate (int)\n");
    printf("\tcat_filename: path to and name of WODEN-style srclist (string)\n");
    printf("\tmetafits_filename: path to MWA and name of metafits file to base\n\t\tsimulation on (string)\n");
    printf("\tchunking_size: size of point source chunks to process (int)\n");
    printf("\n");
    printf("Optionally, the .json can include:\n");
    printf("\tsky_crop_components=True: WODEN crops sources with any component\n\t\tbelow the horizon. Add this arg to include all components\n\t\tabove horizon, regardless of which source they belong to\n");
    exit(1);

  }

  if (strcmp("--help", argv[1]) == 0) {
    printf("Must input a .json settings file to run woden. This file must include:\n");
    printf("\tra0: ra phase centre (float in degrees)\n");
    printf("\tdec0: dec phase centre (float in degrees)\n");
    printf("\tnum_freqs: number of fine frequency channels to simulate (int)\n");
    printf("\tnum_time_steps: number of time steps to simulate (int)\n");
    printf("\tcat_filename: path to and name of WODEN-style srclist (string)\n");
    printf("\tmetafits_filename: path to MWA and name of metafits file to base\n\t\tsimulation on (string)\n");
    printf("\tchunking_size: size of point source chunks to process (int)\n");
    printf("\n");
    printf("Optionally, the .json can include:\n");
    printf("\tsky_crop_components=True: WODEN crops sources with any component\n\t\tbelow the horizon. Add this arg to include all components\n\t\tabove horizon, regardless of which source they belong to\n");
    exit(1);

  }

  float *sbf;
  sbf = NULL;
  sbf = malloc( sbf_N * sbf_L * sizeof(float) );
  sbf = create_sbf(sbf);

  woden_settings_t *woden_settings;
  woden_settings = read_json_settings(argv[1]);

  if (woden_settings->chunking_size > MAX_CHUNKING_SIZE) {
    printf("Current maximum allowable chunk size is %d.  Defaulting to this value.", MAX_CHUNKING_SIZE);
    woden_settings->chunking_size = MAX_CHUNKING_SIZE;
  }
  else if (woden_settings->chunking_size < 1 ) {
    woden_settings->chunking_size = MAX_CHUNKING_SIZE;
  }

  int status=0;
  static fitsfile *metaf_file=NULL;
  MetaFfile_t metafits;
  fits_open_file(&metaf_file, woden_settings->metafits_filename, READONLY, &status);
  status = init_meta_file(metaf_file, &metafits, woden_settings->metafits_filename);

  array_layout_t *array_layout;
  array_layout = calc_XYZ_diffs(&metafits, metafits.num_tiles);

  woden_settings->lst_base = metafits.lst_base;
  woden_settings->base_low_freq = metafits.base_low_freq;
  woden_settings->num_baselines = array_layout->num_baselines;

  float ha0, sha0, cha0;
  float wavelength;
  float frequency;
  const int num_visis = woden_settings->num_baselines * woden_settings->num_time_steps * woden_settings->num_freqs;

  float sdec0,cdec0;
  sdec0 = sin(woden_settings->dec0); cdec0=cos(woden_settings->dec0);

  printf("Setting phase centre (rad) to %f %f\n",woden_settings->ra0, woden_settings->dec0);
  printf("Obs pointing centre (rad) is %f %f\n",metafits.ra_point, metafits.dec_point);

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

  printf("Horizon cropping sky model and calculating az/za for all components for observation\n");

  catsource_t *cropped_src;
  cropped_src = crop_sky_model(raw_srccat, lsts, num_time_steps, woden_settings->sky_crop_type);

  printf("Finished cropping and calculating az/za\n");

  //Setup beam settings for observation
  beam_settings_t beam_settings;
  //Angles used in calculating beam style l,m,ns
  beam_settings.beam_angles_array = malloc(3*sizeof(float));
  beam_settings.beam_angles_array[0] = sinf(metafits.dec_point);
  beam_settings.beam_angles_array[1] = cosf(metafits.dec_point);
  beam_settings.beam_angles_array[2] = woden_settings->lst_base - metafits.ra_point;

  //Number of beam calculations needed for point components
  beam_settings.num_point_beam_values = cropped_src->n_points * woden_settings->num_time_steps * woden_settings->num_freqs;

  if (woden_settings->beamtype == GAUSS_BEAM) {
    beam_settings.beamtype = GAUSSIAN;

    printf("Setting up Gaussian primary beam settings\n");
    printf("   setting beam FWHM to %.5fdeg and ref freq to %.3fMHz\n",woden_settings->gauss_beam_FWHM,woden_settings->gauss_beam_ref_freq / 1e+6  );

    beam_settings.beam_FWHM_rad = woden_settings->gauss_beam_FWHM * D2R;
    //TODO I cannot for the life of me work out how to cudaMalloc and Memcpy
    //a single float (argh) so put the ref freq in an array (embarrassment)
    float beam_ref_freq_array[1] = {woden_settings->gauss_beam_ref_freq};
    beam_settings.beam_ref_freq_array = malloc(sizeof(float));
    beam_settings.beam_ref_freq_array = beam_ref_freq_array;

    //Store all ha (which change with lst) that the beam needs to be calculated at.
    beam_settings.beam_point_has = malloc(woden_settings->num_time_steps * cropped_src->n_points * sizeof(float));
    beam_settings.beam_point_decs = malloc(woden_settings->num_time_steps * cropped_src->n_points * sizeof(float));

    for ( int time_step = 0; time_step < woden_settings->num_time_steps; time_step++ ) {
      float lst = woden_settings->lst_base + time_step*woden_settings->time_res*SOLAR2SIDEREAL*DS2R;
      //TODO add half a time step is good? Add time decorrelation?
      lst += 0.5*woden_settings->time_res*SOLAR2SIDEREAL*DS2R;

      for (int component = 0; component < cropped_src->n_points; component++) {
        int step = cropped_src->n_points*time_step + component;

        beam_settings.beam_point_has[step] = lst - cropped_src->point_ras[component];
        beam_settings.beam_point_decs[step] = cropped_src->point_decs[component];

      }//point loop

    //TODO add in gaussian components
    //TODO add in shapelet components

    }//gaussian beam time loop
  } // End if (woden_settings->gaussian_beam)



  for (size_t band = 0; band < woden_settings->num_bands; band++) {
    int band_num = woden_settings->band_nums[band];
    float base_band_freq = ((band_num - 1)*(metafits.bandwidth/24.0)) + woden_settings->base_low_freq;
    printf("Simulating band %02d with bottom freq %.8e\n",band_num,base_band_freq);

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
    visibility_set->channel_frequencies = malloc( woden_settings->num_freqs * sizeof(float) );

    // //Useful for testing beam things
    // visibility_set->beam_has = malloc(woden_settings->num_time_steps * cropped_src->n_points * sizeof(float) );
    // visibility_set->beam_decs = malloc(woden_settings->num_time_steps * cropped_src->n_points * sizeof(float) );
    // visibility_set->beam_ls = malloc(woden_settings->num_time_steps * cropped_src->n_points * sizeof(float) );
    // visibility_set->beam_ms = malloc(woden_settings->num_time_steps * cropped_src->n_points * sizeof(float) );
    // visibility_set->beam_reals = malloc(woden_settings->num_time_steps * cropped_src->n_points * woden_settings->num_freqs * sizeof(float) );
    // visibility_set->beam_imags = malloc(woden_settings->num_time_steps * cropped_src->n_points * woden_settings->num_freqs * sizeof(float) );

    //Fill in the channel frequencies
    for (int freq_step = 0; freq_step < woden_settings->num_freqs; freq_step++) {
      frequency = base_band_freq + (woden_settings->frequency_resolution*freq_step);
      visibility_set->channel_frequencies[freq_step] = frequency;
    }

    //Fill in visibility settings in order of baseline,freq,time
    //Order matches that of a uvfits file (I live in the past)

    for ( int time_step = 0; time_step < woden_settings->num_time_steps; time_step++ ) {
      float lst = woden_settings->lst_base + time_step*woden_settings->time_res*SOLAR2SIDEREAL*DS2R;
      //TODO add half a time step is good? Add time decorrelation?
      lst += 0.5*woden_settings->time_res*SOLAR2SIDEREAL*DS2R;
      ha0 = lst - woden_settings->ra0;
      sha0 = sin(ha0); cha0=cos(ha0);

      for (int freq_step = 0; freq_step < woden_settings->num_freqs; freq_step++) {
        frequency = base_band_freq + (woden_settings->frequency_resolution*freq_step);
        wavelength = VELC / frequency;

        int step = woden_settings->num_baselines*(time_step*woden_settings->num_freqs + freq_step);

        for (int baseline = 0; baseline < woden_settings->num_baselines; baseline++) {

          visibility_set->cha0s[step + baseline] = cha0;
          visibility_set->sha0s[step + baseline] = sha0;
          visibility_set->lsts[step + baseline] = lst;
          visibility_set->wavelengths[step + baseline] = wavelength;
        }//baseline loop
      }//freq loop
    }//time loop

    //Calculating a single shapelet coeff is equivalent to a point/gauss so treat as a
    //component here
    int num_components = cropped_src->n_points + cropped_src->n_gauss + cropped_src->n_shape_coeffs;

    //TODO should we chunk outside the band for-loop so that we can reuse the chunks for each band (should be the same)
    if (num_components > woden_settings->chunking_size) {
      printf("Chunking sky model\n");
      int num_chunks = num_components / woden_settings->chunking_size;

      if (num_components % woden_settings->chunking_size != 0) {
        num_chunks = num_chunks + 1;
      }

      printf("Number of chunks required is %d\n",  num_chunks);

      //setup a temporary visibility set that calculate_visibilities will populate
      visibility_set_t *temp_visibility_set = malloc(sizeof(visibility_set_t));
      temp_visibility_set->sum_visi_real = malloc( num_visis * sizeof(float) );
      temp_visibility_set->sum_visi_imag = malloc( num_visis * sizeof(float) );
      temp_visibility_set->us_metres = malloc( num_visis * sizeof(float) );
      temp_visibility_set->vs_metres = malloc( num_visis * sizeof(float) );
      temp_visibility_set->ws_metres = malloc( num_visis * sizeof(float) );

      temp_visibility_set->sha0s = visibility_set->sha0s;
      temp_visibility_set->cha0s = visibility_set->cha0s;
      temp_visibility_set->lsts = visibility_set->lsts;
      temp_visibility_set->wavelengths = visibility_set->wavelengths;

      catsource_t *temp_cropped_src = malloc(sizeof(catsource_t));

      //For each chunk, calculate the visibilities for those components
      for (int chunk = 0; chunk < num_chunks; chunk++) {
        printf("Processing chunk %d\n", chunk);

        fill_chunk_src(temp_cropped_src, cropped_src, num_chunks, chunk,
                       woden_settings->chunking_size, woden_settings->num_time_steps );

        printf("\tNumber of components in chunk are: P %d G %d S_coeffs %d\n",temp_cropped_src->n_points,temp_cropped_src->n_gauss,temp_cropped_src->n_shape_coeffs );

        calculate_visibilities(array_layout->X_diff_metres, array_layout->Y_diff_metres, array_layout->Z_diff_metres,
                      *temp_cropped_src, angles_array, beam_settings,
                      woden_settings->num_baselines, woden_settings->num_time_steps,
                      num_visis, woden_settings->num_freqs, visibility_set,
                      sbf);

        // printf("Adding temporary visibility set.\n");
	      //add to visiblity_set
        for (int visi = 0; visi < num_visis; visi++) {
          //if the first chunk then initialise our values, and copy across
          //the u,v,w coords
          if (chunk == 0) {
            visibility_set->sum_visi_real[visi] = 0;
            visibility_set->sum_visi_imag[visi] = 0;

            visibility_set->us_metres[visi] = temp_visibility_set->us_metres[visi];
            visibility_set->vs_metres[visi] = temp_visibility_set->vs_metres[visi];
            visibility_set->ws_metres[visi] = temp_visibility_set->ws_metres[visi];
          }

          //add each chunk of components to visibility set
          visibility_set->sum_visi_real[visi] += temp_visibility_set->sum_visi_real[visi];
          visibility_set->sum_visi_imag[visi] += temp_visibility_set->sum_visi_imag[visi];
        }//visi loop
      }//chunk loop

      free( temp_visibility_set->sum_visi_real );
      free( temp_visibility_set->sum_visi_imag );
      free( temp_visibility_set->us_metres );
      free( temp_visibility_set->vs_metres );
      free( temp_visibility_set->ws_metres );

      free( temp_visibility_set );

    }
    //If not chunking the components, just simulate all in one go
    else {
      calculate_visibilities(array_layout->X_diff_metres, array_layout->Y_diff_metres, array_layout->Z_diff_metres,
                    *cropped_src, angles_array, beam_settings,
                    woden_settings->num_baselines, woden_settings->num_time_steps,
                    num_visis, woden_settings->num_freqs, visibility_set,
                    sbf);

    }

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
    fwrite(visibility_set->sum_visi_real, num_visis*sizeof(float), 1, output_visi);
    fwrite(visibility_set->sum_visi_imag, num_visis*sizeof(float), 1, output_visi);

    fflush(output_visi);
    fclose(output_visi);

    // // // Dumps u,v,w (metres), Re(vis), Im(vis) directly to text file - useful for
    // // // bug hunting with small outputs
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
    //               visibility_set->sum_visi_real[step + baseline],visibility_set->sum_visi_imag[step + baseline]);
    //     }
    //   }
    // }
    // fflush(output_visi_text);
    // fclose(output_visi_text);
    //
    // //// Beam testing text file
    // FILE *output_beamcoords;
    // char bufff[0x100];
    // snprintf(bufff, sizeof(bufff), "output_beam_coords%02d.txt", band_num);
    // output_beamcoords = fopen(bufff,"w");
    // for ( int time_step = 0; time_step < woden_settings->num_time_steps; time_step++ ) {
    //     for (int component = 0; component < cropped_src->n_points; component++) {
    //       int step = cropped_src->n_points*time_step;
    //       fprintf(output_beamcoords,"%f %f %f %f %f %f\n",
    //               beam_settings.beam_point_has[step + component], beam_settings.beam_point_decs[step + component],
    //               visibility_set->beam_has[step + component], visibility_set->beam_decs[step + component],
    //               visibility_set->beam_ls[step + component], visibility_set->beam_ms[step + component]);
    //     }
    // }
    // fflush(output_beamcoords);
    // fclose(output_beamcoords);
    // //
    // //
    // //// Beam testing text file
    // FILE *output_beamvals;
    // char buffff[0x100];
    // snprintf(buffff, sizeof(buffff), "output_beam_values%02d.txt", band_num);
    // output_beamvals = fopen(buffff,"w");
    // for ( int time_step = 0; time_step < woden_settings->num_time_steps; time_step++ ) {
    //   for ( int freq_step = 0; freq_step < woden_settings->num_freqs; freq_step++ ) {
    //     for (int component = 0; component < cropped_src->n_points; component++) {
    //       int step = cropped_src->n_points*(time_step*woden_settings->num_freqs + freq_step);
    //       fprintf(output_beamvals,"%f %f\n",
    //               visibility_set->beam_reals[step + component],visibility_set->beam_imags[step + component]);
    //     }
    //   }
    // }
    // fflush(output_beamvals);
    // fclose(output_beamvals);

    // free( visibility_set->beam_imags );
    // free( visibility_set->beam_reals );
    // free( visibility_set->beam_ls );
    // free( visibility_set->beam_ms );
    // free( visibility_set->beam_has );
    // free( visibility_set->beam_decs );

    free( visibility_set->sum_visi_real );
    free( visibility_set->sum_visi_imag );
    free( visibility_set->us_metres );
    free( visibility_set->vs_metres );
    free( visibility_set->ws_metres );
    free( visibility_set->sha0s );
    free( visibility_set->cha0s );
    free( visibility_set->lsts );
    free( visibility_set->wavelengths );

    free( visibility_set );

  }//band loop
}//main
