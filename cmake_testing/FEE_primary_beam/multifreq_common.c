#include "multifreq_common.h"


//Different delays settings, which control the pointing of the MWA beam
user_precision_t delays1[16] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

user_precision_t delays2[16] = {0.0, 4.0, 8.0, 12.0, 0.0, 4.0, 8.0, 12.0,
                     0.0, 4.0, 8.0, 12.0, 0.0, 4.0, 8.0, 12.0};

user_precision_t delays3[16] = {0.0, 2.0, 4.0, 8.0, 2.0, 4.0, 8.0, 12.0,
                     4.0, 8.0, 12.0, 16.0, 8.0, 12.0, 16.0, 20.0};


void create_metafits_and_beam_freqs(woden_settings_t *woden_settings,
                                    double *beam_freqs, double base_low_freq,
                                    int num_freqs, double freq_res,
                                    user_precision_t* delays) {

  woden_settings->base_low_freq = base_low_freq;
  woden_settings->num_freqs = num_freqs;

  for(int i=0;i<16;i++) {
      woden_settings->FEE_ideal_delays[i] = delays[i];
  }

  for (int find = 0; find < num_freqs; find++) {
    beam_freqs[find] = base_low_freq + find*freq_res;
  }

}

void create_metafits_and_beam_freqs_delays1_freqs1(woden_settings_t *woden_settings,
                                                   double *beam_freqs){

  double freq_res = 40e+3;
  double base_low_freq = 167.035e+6;
  int num_freqs = NUM_FREQS1;

  create_metafits_and_beam_freqs(woden_settings, beam_freqs,
                                base_low_freq, num_freqs,
                                freq_res, delays1);

}


void create_metafits_and_beam_freqs_delays2_freqs2(woden_settings_t *woden_settings,
                                                   double *beam_freqs){

  double freq_res = 1.28e+6;
  double base_low_freq = 167.035e+6;
  int num_freqs = NUM_FREQS2;

  create_metafits_and_beam_freqs(woden_settings, beam_freqs,
                                base_low_freq, num_freqs,
                                freq_res, delays2);

}

void create_metafits_and_beam_freqs_delays3_freqs3(woden_settings_t *woden_settings,
                                                   double *beam_freqs){

  double freq_res = 0.64e+3;
  double base_low_freq = 180.0e+6;
  int num_freqs = NUM_FREQS3;

  create_metafits_and_beam_freqs(woden_settings, beam_freqs,
                                base_low_freq, num_freqs,
                                freq_res, delays3);

}
