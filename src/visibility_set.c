#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "woden_struct_defs.h"
#include "woden_precision_defs.h"


visibility_set_t * setup_visibility_set(int num_visis) {
  visibility_set_t *visibility_set = (visibility_set_t *)malloc(sizeof(visibility_set_t));
  visibility_set->us_metres = (user_precision_t *)malloc( num_visis * sizeof(user_precision_t) );
  visibility_set->vs_metres = (user_precision_t *)malloc( num_visis * sizeof(user_precision_t) );
  visibility_set->ws_metres = (user_precision_t *)malloc( num_visis * sizeof(user_precision_t) );

  visibility_set->sum_visi_XX_real = (user_precision_t *)malloc( num_visis * sizeof(user_precision_t) );
  visibility_set->sum_visi_XX_imag = (user_precision_t *)malloc( num_visis * sizeof(user_precision_t) );
  visibility_set->sum_visi_XY_real = (user_precision_t *)malloc( num_visis * sizeof(user_precision_t) );
  visibility_set->sum_visi_XY_imag = (user_precision_t *)malloc( num_visis * sizeof(user_precision_t) );
  visibility_set->sum_visi_YX_real = (user_precision_t *)malloc( num_visis * sizeof(user_precision_t) );
  visibility_set->sum_visi_YX_imag = (user_precision_t *)malloc( num_visis * sizeof(user_precision_t) );
  visibility_set->sum_visi_YY_real = (user_precision_t *)malloc( num_visis * sizeof(user_precision_t) );
  visibility_set->sum_visi_YY_imag = (user_precision_t *)malloc( num_visis * sizeof(user_precision_t) );

  return visibility_set;
}

void fill_timefreq_visibility_set(visibility_set_t *visibility_set,
                                  woden_settings_t *woden_settings,
                                  double base_band_freq,
                                  double *lsts) {
  //For easy indexing when running on GPUs, make 4 arrays that match
  //the settings for every baseline, frequency, and time step in the
  //simulation.
  //Fill in visibility settings in order of baseline,freq,time
  //Order matches that of a uvfits file (I live in the past)

  int num_visis = woden_settings->num_visis;
  visibility_set->allsteps_sha0s = malloc( num_visis * sizeof(double) );
  visibility_set->allsteps_cha0s = malloc( num_visis * sizeof(double) );
  visibility_set->allsteps_lsts = malloc( num_visis * sizeof(double) );
  visibility_set->allsteps_wavelengths = malloc( num_visis * sizeof(user_precision_t) );
  visibility_set->channel_frequencies = malloc( (int)woden_settings->num_freqs * sizeof(double) );

  user_precision_t wavelength;
  double frequency;
  //Fill in the fine channel frequencies
  for (int freq_step = 0; freq_step < woden_settings->num_freqs; freq_step++) {
    frequency = base_band_freq + (woden_settings->frequency_resolution*freq_step);
    visibility_set->channel_frequencies[freq_step] = frequency;
  }

  //Fill the time/frequency settings for all baseline/freq/time
  double ha0, sha0, cha0;

  for ( int time_step = 0; time_step < woden_settings->num_time_steps; time_step++ ) {
    ha0 = lsts[time_step] - woden_settings->ra0;
    sha0 = sin(ha0); cha0=cos(ha0);

    for (int freq_step = 0; freq_step < woden_settings->num_freqs; freq_step++) {
      frequency = base_band_freq + (woden_settings->frequency_resolution*freq_step);
      wavelength = (VELC / frequency);
      int step = woden_settings->num_baselines*(time_step*woden_settings->num_freqs + freq_step);

      for (int baseline = 0; baseline < woden_settings->num_baselines; baseline++) {
        visibility_set->allsteps_cha0s[step + baseline] = cha0;
        visibility_set->allsteps_sha0s[step + baseline] = sha0;
        visibility_set->allsteps_lsts[step + baseline] = lsts[time_step];
        visibility_set->allsteps_wavelengths[step + baseline] = wavelength;
      }//baseline loop
    }//freq loop
  }//time loop
}

void write_visi_set_binary(visibility_set_t *visibility_set,
                           int band_num, int num_visis) {
  //Dumps u,v,w (metres), Re(vis), Im(vis) to a binary file
  //Output order is by baseline (fastest changing), frequency, time (slowest
  //changing)
  //This means the us_metres, vs_metres, ws_metres are repeated over frequency,
  //but keeps the dimensions of the output sane
  FILE *output_visi;
  char buf[0x100];
  snprintf(buf, sizeof(buf), "output_visi_band%02d.dat", band_num);

  output_visi = fopen(buf,"wb");

  if(output_visi == NULL)
  {
      printf("Could not open output_visi_band%02d.dat - exiting", band_num);
      exit(1);
  }

  fwrite(visibility_set->us_metres, num_visis*sizeof(user_precision_t), 1, output_visi);
  fwrite(visibility_set->vs_metres, num_visis*sizeof(user_precision_t), 1, output_visi);
  fwrite(visibility_set->ws_metres, num_visis*sizeof(user_precision_t), 1, output_visi);

  fwrite(visibility_set->sum_visi_XX_real, num_visis*sizeof(user_precision_t), 1, output_visi);
  fwrite(visibility_set->sum_visi_XX_imag, num_visis*sizeof(user_precision_t), 1, output_visi);
  fwrite(visibility_set->sum_visi_XY_real, num_visis*sizeof(user_precision_t), 1, output_visi);
  fwrite(visibility_set->sum_visi_XY_imag, num_visis*sizeof(user_precision_t), 1, output_visi);
  fwrite(visibility_set->sum_visi_YX_real, num_visis*sizeof(user_precision_t), 1, output_visi);
  fwrite(visibility_set->sum_visi_YX_imag, num_visis*sizeof(user_precision_t), 1, output_visi);
  fwrite(visibility_set->sum_visi_YY_real, num_visis*sizeof(user_precision_t), 1, output_visi);
  fwrite(visibility_set->sum_visi_YY_imag, num_visis*sizeof(user_precision_t), 1, output_visi);

  fflush(output_visi);
  fclose(output_visi);
}

void write_visi_set_text(visibility_set_t *visibility_set, int band_num,
                         woden_settings_t *woden_settings) {

  // Dumps u,v,w (metres), XX Re(vis), XX Im(vis) directly to text file - useful for
  // bug hunting with small outputs
  FILE *output_visi_text;
  char buff[0x100];
  snprintf(buff, sizeof(buff), "output_visi_band%02d.txt", band_num);
  output_visi_text = fopen(buff,"w");
  for ( int time_step = 0; time_step < woden_settings->num_time_steps; time_step++ ) {
    for ( int freq_step = 0; freq_step < woden_settings->num_freqs; freq_step++ ) {
      for (int baseline = 0; baseline < woden_settings->num_baselines; baseline++) {
        int step = woden_settings->num_baselines*(time_step*woden_settings->num_freqs + freq_step);

        #ifdef DOUBLE_PRECISION
          fprintf(output_visi_text,"%lf %lf %lf %lf %lf\n",
                  visibility_set->us_metres[step + baseline],
                  visibility_set->vs_metres[step + baseline],
                  visibility_set->ws_metres[step + baseline],
                  visibility_set->sum_visi_XX_real[step + baseline],
                  visibility_set->sum_visi_XX_imag[step + baseline]);
        #else
          fprintf(output_visi_text,"%f %f %f %f %f\n",
                  visibility_set->us_metres[step + baseline],
                  visibility_set->vs_metres[step + baseline],
                  visibility_set->ws_metres[step + baseline],
                  visibility_set->sum_visi_XX_real[step + baseline],
                  visibility_set->sum_visi_XX_imag[step + baseline]);
        #endif
      }
    }
  }
  fflush(output_visi_text);
  fclose(output_visi_text);
}

void free_visi_set_inputs(visibility_set_t *visibility_set) {
  free(visibility_set->allsteps_sha0s);
  free(visibility_set->allsteps_cha0s);
  free(visibility_set->allsteps_lsts);
  free(visibility_set->allsteps_wavelengths);
  free(visibility_set->channel_frequencies);
}

void free_visi_set_outputs(visibility_set_t *visibility_set) {
  free(visibility_set->us_metres);
  free(visibility_set->vs_metres);
  free(visibility_set->ws_metres);

  free(visibility_set->sum_visi_XX_real);
  free(visibility_set->sum_visi_XX_imag);
  free(visibility_set->sum_visi_XY_real);
  free(visibility_set->sum_visi_XY_imag);
  free(visibility_set->sum_visi_YX_real);
  free(visibility_set->sum_visi_YX_imag);
  free(visibility_set->sum_visi_YY_real);
  free(visibility_set->sum_visi_YY_imag);
}
