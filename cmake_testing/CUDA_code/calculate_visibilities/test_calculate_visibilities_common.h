#pragma once

#include <math.h>
#include <unity.h>
#include <stdlib.h>
#include <complex.h>

#include "constants.h"
#include "woden_struct_defs.h"
#include "shapelet_basis.h"
#include "woden_settings.h"
#include "visibility_set.h"
#include "shapelet_basis.h"
// #include "hdf5.h"
// #include "hdf5_hl.h"

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

//External CUDA code we're linking in
extern void calculate_visibilities(array_layout_t *array_layout,
  source_catalogue_t *cropped_sky_models, beam_settings_t *beam_settings,
  woden_settings_t *woden_settings, visibility_set_t *visibility_set,
  user_precision_t *sbf);

source_catalogue_t * make_cropped_sky_models(double ra0, double dec0,
                                             int n_points, int n_gauss,
                                             int n_shapes,
                                             int num_sources);

void free_sky_model(source_catalogue_t *cropped_sky_models);

woden_settings_t * make_woden_settings(double ra0, double dec0);

void test_uvw(visibility_set_t *visibility_set,  double *lsts,
              double ra0, double dec0, woden_settings_t *woden_settings);

visibility_set_t * test_calculate_visibilities(source_catalogue_t *cropped_sky_models,
                                 beam_settings_t *beam_settings,
                                 woden_settings_t *woden_settings,
                                 double ra0, double dec0,
                                 int beamtype);

void test_comp_phase_centre_twogains(visibility_set_t *visibility_set,
                                     double gain1xx, double gain1yy,
                                     double gain2xx, double gain2yy,
                                     woden_settings_t *woden_settings);

void test_comp_phase_centre_allgains(visibility_set_t *visibility_set,
                                     double gain1xx_re, double gain1xx_im,
                                     double gain1xy_re, double gain1xy_im,
                                     double gain1yx_re, double gain1yx_im,
                                     double gain1yy_re, double gain1yy_im,
                                     double gain2xx_re, double gain2xx_im,
                                     double gain2xy_re, double gain2xy_im,
                                     double gain2yx_re, double gain2yx_im,
                                     double gain2yy_re, double gain2yy_im,
                                     woden_settings_t *woden_settings,
                                     double tol);
