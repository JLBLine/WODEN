#pragma once

#include <unity.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "woden_struct_defs.h"

#define NUM_FREQS1 32
#define NUM_FREQS2 5
#define NUM_FREQS3 5

//Different delays settings, which control the pointing of the MWA beam
user_precision_t delays1[16];

user_precision_t delays2[16];

user_precision_t delays3[16];

void create_metafits_and_beam_freqs(woden_settings_t *woden_settings,
                                    double *beam_freqs, double base_low_freq,
                                    int num_freqs, double freq_res,
                                    user_precision_t* delays);

void create_metafits_and_beam_freqs_delays1_freqs1(woden_settings_t *woden_settings,
                                                   double *beam_freqs);

void create_metafits_and_beam_freqs_delays2_freqs2(woden_settings_t *woden_settings,
                                                   double *beam_freqs);

void create_metafits_and_beam_freqs_delays3_freqs3(woden_settings_t *woden_settings,
                                                   double *beam_freqs);
