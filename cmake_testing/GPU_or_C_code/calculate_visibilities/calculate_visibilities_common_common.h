#pragma once

#include <math.h>
#include <unity.h>
#include <stdlib.h>
#include <complex.h>

#include "constants.h"
#include "woden_struct_defs.h"
#include "shapelet_basis.h"
#include "visibility_set.h"
#include "shapelet_basis.h"
#include "hyperbeam_error.h"
// #include "calculate_visibilities_common.h"

#define NUM_ANTS 3
#define NUM_BASELINES 3
#define NUM_FREQS 3
#define NUM_TIME_STEPS 2
#define NUM_CROSS NUM_BASELINES*NUM_FREQS*NUM_TIME_STEPS
#define NUM_VISI NUM_CROSS
#define RA0 0.0
#define BASE_BAND_FREQ 120000000.0
#define STOKESI 0.3333333333333333
// #define STOKESI 1.0
#define POL_FRAC 1.0

source_catalogue_t * make_cropped_sky_models(double ra0, double dec0,
                                             int n_points, int n_gauss,
                                             int n_shapes,
                                             int num_sources);

void free_sky_model(source_catalogue_t *cropped_sky_models);

woden_settings_t * make_woden_settings(double ra0, double dec0);

void test_uvw(visibility_set_t *visibility_set,  double *lsts,
              double ra0, double dec0, woden_settings_t *woden_settings);

void fill_array_layout(array_layout_t *array_layout);

visibility_set_t * test_calculate_visibilities(source_catalogue_t *cropped_sky_models,
                                 beam_settings_t *beam_settings,
                                 woden_settings_t *woden_settings,
                                 double ra0, double dec0,
                                 int beamtype);

//Due to the way I've set up the sky model, the I,Q,U,V fluxes should all be
//some multiple of the number of components. Everything is also stuck at phase
//centre, so visis should be the same for all baselines. Give given the complex
//beam values, we can predict the visibilities
void predict_inst_stokes(int num_comps, double _Complex g1x, double _Complex D1x,
                                        double _Complex D1y, double _Complex g1y,
                                        double _Complex g2x, double _Complex D2x,
                                        double _Complex D2y, double _Complex g2y,
                                        double _Complex * xx, double _Complex * xy,
                                        double _Complex * yx, double _Complex * yy);

void test_comp_phase_centre_twogains(visibility_set_t *visibility_set,
                                     int num_comps,
                                     double _Complex gain1x, double _Complex gain1y,
                                     double _Complex gain2x, double _Complex gain2y,
                                     woden_settings_t *woden_settings);

void test_comp_phase_centre_allgains(visibility_set_t *visibility_set,
                              int num_comps,
                              double _Complex gain1x, double _Complex leak1x,
                              double _Complex leak1y, double _Complex gain1y,
                              double _Complex gain2x, double _Complex leak2x,
                              double _Complex leak2y, double _Complex gain2y,
                              woden_settings_t *woden_settings, double tol);

//Aight this one tests things when we have all gains varying, and all antennas
//have different primary beams
void test_comp_phase_centre_allgains_multiants(visibility_set_t *visibility_set,
                                int num_comps, 
                                double _Complex gain1x, double _Complex leak1x,
                                double _Complex leak1y, double _Complex gain1y,
                                double _Complex gain2x, double _Complex leak2x,
                                double _Complex leak2y, double _Complex gain2y,
                                double *antx_mult, double *anty_mult, int num_ants,
                                woden_settings_t *woden_settings,
                                double tol);

void test_comp_phase_centre_allgains_diffants(visibility_set_t *visibility_set,
                                int num_comps, 
                                double _Complex *expec_gainx, double _Complex *expec_leakx,
                                double _Complex *expec_leaky, double _Complex *expec_gainy,
                                int num_ants,
                                woden_settings_t *woden_settings,
                                double tol);