#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fitsio.h>
#include <erfa.h>
// #include "read_and_write.h"
#include "create_sky_model.h"
#include "shapelet_basis.h"
#include "woden.h"
#include "constants.h"

int main(int argc, char **argv) {

  if (argc < 2) {
    printf("Must input a .json settings file to run woden. This file must include:\n");
    printf("\tra0: ra phase centre (float in degrees)\n");
    printf("\tdec0: dec phase centre (float in degrees)\n");
    printf("\tnum_freqs: number of fine frequency channels to simulate (int)\n");
    printf("\tnum_time_steps: number of time steps to simulate (int)\n");
    printf("\tcat_filename: path to and name of RTS-new-style srclist (string)\n");
    printf("\tmetafits_filename: path to MWA and name of metafits file to base\n\t\tsimulation on (string)\n");
    exit(1);

  }

  if (strcmp("--help", argv[1]) == 0) {
    printf("wodan needs a .json settings file to run. This file must include:\n");
    printf("\tra0: ra phase centre (float in degrees)\n");
    printf("\tdec0: dec phase centre (float in degrees)\n");
    printf("\tnum_freqs: number of fine frequency channels to simulate (int)\n");
    printf("\tnum_time_steps: number of time steps to simulate (int)\n");
    printf("\tcat_filename: path to and name of RTS-new-style srclist (string)\n");
    printf("\tmetafits_filename: path to MWA and name of metafits file to base\n\t\tsimulation on (string)\n");
    exit(1);

  }

  float *sbf;
  sbf = NULL;
  sbf = malloc( sbf_N * sbf_L * sizeof(float) );
  sbf = create_sbf(sbf);

  woden_settings_t * woden_settings;
  woden_settings = read_json_settings(argv[1]);

  int status=0;
  // array_layout_t * array_layout;
  static fitsfile *metaf_file=NULL;
  MetaFfile_t metafits;
  fits_open_file(&metaf_file, woden_settings->metafits_filename, READONLY, &status);
  status = init_meta_file(metaf_file, &metafits, woden_settings->metafits_filename);

  array_layout_t * array_layout;
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

    calculate_visibilities(array_layout->X_diff_metres, array_layout->Y_diff_metres, array_layout->Z_diff_metres,
                    *cropped_src, angles_array,
                    woden_settings->num_baselines, num_visis, visibility_set,
                    sbf);

    // calculate_visibilities(array_layout->X_diff_metres, array_layout->Y_diff_metres, array_layout->Z_diff_metres,
    //                 raw_srccat->catsources[0], angles_array,
    //                 woden_settings->num_baselines, num_visis, visibility_set,
    //                 sbf);

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

    // // Dumps u,v,w (metres), Re(vis), Im(vis) directly to text file - useful for
    // // bug hunting with small outputs
    // char buff[0x100];
    // snprintf(buff, sizeof(buff), "output_visi_band%02d.txt", band_num);
    // output_visi = fopen(buff,"w");
    // for ( int time_step = 0; time_step < woden_settings->num_time_steps; time_step++ ) {
    //   for ( int freq_step = 0; freq_step < woden_settings->num_freqs; freq_step++ ) {
    //     for (int baseline = 0; baseline < woden_settings->num_baselines; baseline++) {
    //       int step = woden_settings->num_baselines*(time_step*woden_settings->num_freqs + freq_step);
    //       fprintf(output_visi,"%f %f %f %f %f\n",visibility_set->us_metres[step + baseline],
    //               visibility_set->vs_metres[step + baseline],visibility_set->ws_metres[step + baseline],
    //               visibility_set->sum_visi_real[step + baseline],visibility_set->sum_visi_imag[step + baseline]);
    //     }
    //   }
    // }
    // fclose(output_visi);

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
