#pragma once

#include <math.h>
#include <unity.h>
#include <stdlib.h>
#include <complex.h>

#include "constants.h"
#include "woden_struct_defs.h"
#include "FEE_primary_beam.h"
#include "shapelet_basis.h"
#include "woden_settings.h"
#include "visibility_set.h"
#include "shapelet_basis.h"
#include "hdf5.h"
#include "hdf5_hl.h"

#define NUM_BASELINES 3
#define NUM_FREQS 3
#define NUM_TIME_STEPS 2
#define NUM_VISI NUM_BASELINES*NUM_FREQS*NUM_TIME_STEPS
#define RA0 0.0
#define BASE_BAND_FREQ 120000000.0

source_catalogue_t * make_cropped_sky_models(float ra0, float dec0,
                                             int n_points, int n_gauss,
                                             int n_shapes,
                                             int num_sources);

void free_sky_model(source_catalogue_t *cropped_sky_models);

woden_settings_t * make_woden_settings(float ra0, float dec0);

void test_uvw(visibility_set_t *visibility_set,  float *lsts,
              float ra0, float dec0);

visibility_set_t * test_calculate_visibilities(source_catalogue_t *cropped_sky_models,
                                 beam_settings_t *beam_settings,
                                 woden_settings_t *woden_settings,
                                 float ra0, float dec0,
                                 int beamtype);

void test_comp_phase_centre_twogains(visibility_set_t *visibility_set,
                                     float gain1xx, float gain1yy,
                                     float gain2xx, float gain2yy);
