#pragma once

#include "constants.h"
#include "woden_precision_defs.h"

const int NUM_COMPS = 5;

const double FREQ1 =  167e+6;
const double FREQ_RES1 = 80e+3;
const int NUM_FREQS1 = 32;

const double FREQ2 =  167e+6;
const double FREQ_RES2 = 1.28e+6;
const int NUM_FREQS2 = 12;

const double FREQ3 =  190e+6;
const double FREQ_RES3 = 320e+3;
const int NUM_FREQS3 = 12;

double azs[] = {270.0*DD2R, 270.0*DD2R, 0.0*DD2R, 90.0*DD2R, 90.0*DD2R};
double zas[] = {60.0*DD2R, 30.0*DD2R, 0.0*DD2R, 30.0*DD2R, 60.0*DD2R};

//Different delays settings, which control the pointing of the MWA beam
user_precision_t delays1[16] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

user_precision_t delays2[16] = {0.0, 4.0, 8.0, 12.0, 0.0, 4.0, 8.0, 12.0,
                     0.0, 4.0, 8.0, 12.0, 0.0, 4.0, 8.0, 12.0};

user_precision_t delays3[16] = {0.0, 2.0, 4.0, 8.0, 2.0, 4.0, 8.0, 12.0,
                     4.0, 8.0, 12.0, 16.0, 8.0, 12.0, 16.0, 20.0};

// user_precision_t azs[] = {260.0*DD2R, 260.0*DD2R, 0.0*DD2R, 100.0*DD2R, 100.0*DD2R};
// user_precision_t zas[] = {60.0*DD2R, 30.0*DD2R, 0.0*DD2R, 30.0*DD2R, 60.0*DD2R};
