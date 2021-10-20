#include <math.h>
#include <unity.h>
#include <stdlib.h>
#include <complex.h>
#include <erfa.h>

#include "constants.h"
#include "woden_struct_defs.h"
#include "FEE_primary_beam.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

void sincosf(float x, float *sin, float *cos);

//External CUDA code we're linking in
extern void test_source_component_common(int num_components,
           float _Complex *primay_beam_J00, float _Complex *primay_beam_J01,
           float _Complex *primay_beam_J10, float _Complex *primay_beam_J11,
           float *freqs, float *ls, float *ms, float *ns,
           float *ras, float *decs, float *azs, float *zas,
           float *sin_para_angs, float *cos_para_angs,
           float *beam_has, float *beam_decs,
           woden_settings_t *woden_settings,
           beam_settings_t *beam_settings);

extern void get_HDFBeam_normalisation(RTS_MWA_FEE_beam_t *FEE_beam_zenith,
                RTS_MWA_FEE_beam_t *FEE_beam);

extern void free_FEE_primary_beam_from_GPU(RTS_MWA_FEE_beam_t *primary_beam);

extern void copy_FEE_primary_beam_to_GPU(RTS_MWA_FEE_beam_t *FEE_beam);

#define UNITY_INCLUDE_FLOAT

/*******************************************************************************
Here are a bunch of predefined arrays to be used
Some are inputs for tests, others are expected outputs
*******************************************************************************/

float azs[] = {4.71238899, 4.71238899, 4.71238899, 4.71238899, 4.71238899,
               4.71238899, 4.71238899, 4.71238899, 4.71238899, 4.71238899,
               4.71238899, 4.71238899, 0.00000000, 0.00000000, 0.00000000,
               1.57079637, 1.57079637, 1.57079637, 1.57079637, 1.57079637,
               1.57079637, 1.57079637, 1.57079637, 1.57079637, 1.57079637,
               1.57079637, 1.57079637};
float zas[] = {1.57079631, 1.57079631, 1.57079631, 1.04719766, 1.04719766,
               1.04719766, 0.78539832, 0.78539832, 0.78539832, 0.52359898,
               0.52359898, 0.52359898, 0.00000000, 0.00000000, 0.00000000,
               0.52359875, 0.52359875, 0.52359875, 0.78539820, 0.78539820,
               0.78539820, 1.04719760, 1.04719760, 1.04719760, 1.57079637,
               1.57079637, 1.57079637};

float sin_para_angs[] = {1.00000000, 1.00000000, 1.00000000, 1.00000000,
                         1.00000000, 1.00000000, 1.00000000, 1.00000000,
                         1.00000000, 1.00000000, 1.00000000, 1.00000000,
                         0.00000000, 0.00000000, 0.00000000, -1.00000000,
                         -1.00000000, -1.00000000, -1.00000000, -1.00000000,
                         -1.00000000, -1.00000000, -1.00000000, -1.00000000,
                         -1.00000000, -1.00000000, -1.00000000};

float cos_para_angs[] = {0.00000000, 0.00000000, 0.00000000, 0.00000000,
                         0.00000000, 0.00000000, 0.00000000, 0.00000000,
                         0.00000000, 0.00000000, 0.00000000, 0.00000000,
                         1.00000000, 1.00000000, 1.00000000, 0.00000000,
                         0.00000000, 0.00000000, 0.00000000, 0.00000000,
                         0.00000000, 0.00000000, 0.00000000, 0.00000000,
                         0.00000000, 0.00000000, 0.00000000};

float gauss_expected[] = { 0.2806705, 0.3856093, 0.5297833, 0.7278620,
                           1.0000000, 0.7278622, 0.5297834, 0.3856093,
                           0.2806705, 0.0062056, 0.0221101, 0.0787758,
                           0.2806702, 1.0000000, 0.2806705, 0.0787759,
                           0.0221101, 0.0062056,
                           0.2806705, 0.3856093, 0.5297833, 0.7278620,
                           1.0000000, 0.7278622, 0.5297834, 0.3856093,
                           0.2806705, 0.0062056, 0.0221101, 0.0787758,
                           0.2806702, 1.0000000, 0.2806705, 0.0787759,
                           0.0221101, 0.0062056,
                           0.2806705, 0.3856093, 0.5297833, 0.7278620,
                           1.0000000, 0.7278622, 0.5297834, 0.3856093,
                           0.2806705, 0.0062056, 0.0221101, 0.0787758,
                           0.2806702, 1.0000000, 0.2806705, 0.0787759,
                           0.0221101, 0.0062056};

float analy_expec_J00[] = { 0.0000000, 0.5257682, 0.7312802, 0.8807548, 1.0000000,
                    0.8807548, 0.7312803, 0.5257683, 0.0000000, 0.0000000,
                    0.6182293, 0.8162959, 0.9315211, 1.0000000, 0.9315211,
                    0.8162960, 0.6182294, 0.0000000, 0.0000001, 0.5257682,
                    0.7312802, 0.8807548, 1.0000000, 0.8807548, 0.7312803,
                    0.5257683, 0.0000000, 0.0000000, 0.6182293, 0.8162959,
                    0.9315211, 1.0000000, 0.9315211, 0.8162960, 0.6182294,
                    0.0000000, 0.0000000, 0.5257682, 0.7312802, 0.8807548,
                    1.0000000, 0.8807548, 0.7312803, 0.5257683, 0.0000000,
                    0.0000001, 0.6182293, 0.8162959, 0.9315211, 1.0000000,
                    0.9315211, 0.8162960, 0.6182294, 0.0000000 };

float analy_expec_J11[] = { 0.0000000, 0.2628839, 0.5170931, 0.7627559, 1.0000000,
                    0.7627561, 0.5170933, 0.2628841, 0.0000000, 0.0000000,
                    0.3091144, 0.5772082, 0.8067207, 1.0000000, 0.8067210,
                    0.5772084, 0.3091146, 0.0000000, 0.0000000, 0.2628839,
                    0.5170931, 0.7627559, 1.0000000, 0.7627561, 0.5170933,
                    0.2628841, 0.0000000, 0.0000000, 0.3091144, 0.5772082,
                    0.8067207, 1.0000000, 0.8067210, 0.5772084, 0.3091146,
                    0.0000000, 0.0000000, 0.2628839, 0.5170931, 0.7627559,
                    1.0000000, 0.7627561, 0.5170933, 0.2628841, 0.0000000,
                    0.0000000, 0.3091144, 0.5772082, 0.8067207, 1.0000000,
                    0.8067210, 0.5772084, 0.3091146, 0.0000000 };

float fee_expec_J00_re[] = { -0.0000036, -0.0000707, 0.0000195, 0.0000860,
                  0.0002512, 0.0000867, 0.0000197, -0.0000710, -0.0000035,
                  -0.0000036, -0.0000707, 0.0000195, 0.0000860, 0.0002512,
                  0.0000867, 0.0000197, -0.0000710, -0.0000035, -0.0000036,
                  -0.0000707, 0.0000195, 0.0000860, 0.0002512, 0.0000867,
                  0.0000197, -0.0000710, -0.0000035, -0.0000036, -0.0000707,
                  0.0000195, 0.0000860, 0.0002512, 0.0000867, 0.0000197,
                  -0.0000710, -0.0000035, -0.0000036, -0.0000707, 0.0000195,
                  0.0000860, 0.0002512, 0.0000867, 0.0000197, -0.0000710,
                  -0.0000035, -0.0000036, -0.0000707, 0.0000195, 0.0000860,
                  0.0002512, 0.0000867, 0.0000197, -0.0000710, -0.0000035 };

float fee_expec_J00_im[] = { 0.0000008, -0.0000328, -0.0000167, 0.0000288,
                  0.0001243, 0.0000291, -0.0000164, -0.0000324, 0.0000006,
                  0.0000008, -0.0000328, -0.0000167, 0.0000288, 0.0001243,
                  0.0000291, -0.0000164, -0.0000324, 0.0000006, 0.0000008,
                  -0.0000328, -0.0000167, 0.0000288, 0.0001243, 0.0000291,
                  -0.0000164, -0.0000324, 0.0000006, 0.0000008, -0.0000328,
                  -0.0000167, 0.0000288, 0.0001243, 0.0000291, -0.0000164,
                  -0.0000324, 0.0000006, 0.0000008, -0.0000328, -0.0000167,
                  0.0000288, 0.0001243, 0.0000291, -0.0000164, -0.0000324,
                  0.0000006, 0.0000008, -0.0000328, -0.0000167, 0.0000288,
                  0.0001243, 0.0000291, -0.0000164, -0.0000324, 0.0000006 };

float fee_expec_J01_re[] = { -0.0029124, 0.0510724, 0.2057645, 0.0950377,
                   0.9703817, 0.0950371, 0.2057633, 0.0510722, -0.0029126,
                   -0.0029124, 0.0510724, 0.2057645, 0.0950377, 0.9703817,
                   0.0950371, 0.2057633, 0.0510722, -0.0029126, -0.0029124,
                   0.0510724, 0.2057645, 0.0950377, 0.9703817, 0.0950371,
                   0.2057633, 0.0510722, -0.0029126, -0.0029124, 0.0510724,
                   0.2057645, 0.0950377, 0.9703817, 0.0950371, 0.2057633,
                   0.0510722, -0.0029126, -0.0029124, 0.0510724, 0.2057645,
                   0.0950377, 0.9703817, 0.0950371, 0.2057633, 0.0510722,
                   -0.0029126, -0.0029124, 0.0510724, 0.2057645, 0.0950377,
                   0.9703817, 0.0950371, 0.2057633, 0.0510722, -0.0029126 };

float fee_expec_J01_im[] = { -0.0000438, 0.0094490, 0.0335542, 0.0188814,
                  0.2415767, 0.0188819, 0.0335540, 0.0094488, -0.0000439,
                  -0.0000438, 0.0094490, 0.0335542, 0.0188814, 0.2415767,
                  0.0188819, 0.0335540, 0.0094488, -0.0000439, -0.0000438,
                  0.0094490, 0.0335542, 0.0188814, 0.2415767, 0.0188819,
                  0.0335540, 0.0094488, -0.0000439, -0.0000438, 0.0094490,
                  0.0335542, 0.0188814, 0.2415767, 0.0188819, 0.0335540,
                  0.0094488, -0.0000439, -0.0000438, 0.0094490, 0.0335542,
                  0.0188814, 0.2415767, 0.0188819, 0.0335540, 0.0094488,
                  -0.0000439, -0.0000438, 0.0094490, 0.0335542, 0.0188814,
                  0.2415767, 0.0188819, 0.0335540, 0.0094488, -0.0000439 };

float fee_expec_J10_re[] = { 0.0016467, -0.0212285, -0.1320742, -0.0687671,
                -0.9701791, -0.0687667, -0.1320757, -0.0212311, 0.0016451,
                0.0016467, -0.0212285, -0.1320742, -0.0687671, -0.9701791,
                -0.0687667, -0.1320757, -0.0212311, 0.0016451, 0.0016467,
                -0.0212285, -0.1320742, -0.0687671, -0.9701791, -0.0687667,
                -0.1320757, -0.0212311, 0.0016451, 0.0016467, -0.0212285,
                -0.1320742, -0.0687671, -0.9701791, -0.0687667, -0.1320757,
                -0.0212311, 0.0016451, 0.0016467, -0.0212285, -0.1320742,
                -0.0687671, -0.9701791, -0.0687667, -0.1320757, -0.0212311,
                0.0016451, 0.0016467, -0.0212285, -0.1320742, -0.0687671,
                -0.9701791, -0.0687667, -0.1320757, -0.0212311, 0.0016451 };

float fee_expec_J10_im[] = { 0.0008441, 0.0186237, -0.0165366, -0.0381620,
                -0.2423892, -0.0381634, -0.0165374, 0.0186229, 0.0008454,
                0.0008441, 0.0186237, -0.0165366, -0.0381620, -0.2423892,
                -0.0381634, -0.0165374, 0.0186229, 0.0008454, 0.0008441,
                0.0186237, -0.0165366, -0.0381620, -0.2423892, -0.0381634,
                -0.0165374, 0.0186229, 0.0008454, 0.0008441, 0.0186237,
                -0.0165366, -0.0381620, -0.2423892, -0.0381634, -0.0165374,
                0.0186229, 0.0008454, 0.0008441, 0.0186237, -0.0165366,
                -0.0381620, -0.2423892, -0.0381634, -0.0165374, 0.0186229,
                0.0008454, 0.0008441, 0.0186237, -0.0165366, -0.0381620,
                -0.2423892, -0.0381634, -0.0165374, 0.0186229, 0.0008454 };

float fee_expec_J11_re[] = { -0.0000014, -0.0000088, -0.0000159, -0.0000640,
                 -0.0002477, -0.0000634, -0.0000149, -0.0000075, -0.0000018,
                 -0.0000014, -0.0000088, -0.0000159, -0.0000640, -0.0002477,
                 -0.0000634, -0.0000149, -0.0000075, -0.0000018, -0.0000014,
                 -0.0000088, -0.0000159, -0.0000640, -0.0002477, -0.0000634,
                 -0.0000149, -0.0000075, -0.0000018, -0.0000014, -0.0000088,
                 -0.0000159, -0.0000640, -0.0002477, -0.0000634, -0.0000149,
                 -0.0000075, -0.0000018, -0.0000014, -0.0000088, -0.0000159,
                 -0.0000640, -0.0002477, -0.0000634, -0.0000149, -0.0000075,
                 -0.0000018, -0.0000014, -0.0000088, -0.0000159, -0.0000640,
                 -0.0002477, -0.0000634, -0.0000149, -0.0000075, -0.0000018 };

float fee_expec_J11_im[] = { 0.0000041, 0.0000744, 0.0000205, -0.0001054,
                -0.0001180, -0.0001041, 0.0000226, 0.0000758, 0.0000037,
                0.0000041, 0.0000744, 0.0000205, -0.0001054, -0.0001180,
                -0.0001041, 0.0000226, 0.0000758, 0.0000037, 0.0000041,
                0.0000744, 0.0000205, -0.0001054, -0.0001180, -0.0001041,
                0.0000226, 0.0000758, 0.0000037, 0.0000041, 0.0000744,
                0.0000205, -0.0001054, -0.0001180, -0.0001041, 0.0000226,
                0.0000758, 0.0000037, 0.0000041, 0.0000744, 0.0000205,
                -0.0001054, -0.0001180, -0.0001041, 0.0000226, 0.0000758,
                0.0000037, 0.0000041, 0.0000744, 0.0000205, -0.0001054,
                -0.0001180, -0.0001041, 0.0000226, 0.0000758, 0.0000037 };

/*
Test that l,m,n and beam values are calculated correctly by `source_component_common`
for a constant Dec=0, latitude=0, and all beam types
There are other tests for the l,m,n coords and beam functions, so no need to
test millions of scenarios here, so stick with Dec=0
*/
void test_source_component_common_ConstantDecChooseBeams(int beamtype, char* mwa_fee_hdf5) {

  //Set up some test condition inputs
  int num_times = 3;
  int num_freqs = 2;
  int num_components = 9;
  int num_baselines = 5;

  int num_beam_values = num_times*num_freqs*num_components;

  float ra0 = 0.0*DD2R;
  float dec0 = 0.0*DD2R;

  float *decs = malloc(num_components*sizeof(float));
  float *zeroes = calloc(num_components, sizeof(float));

  //Keep RA between 0 and 2*pi here but enter RAs that should return
  //negative l values
  float ras[9] = {(3*M_PI)/2, (5*M_PI)/3, (7*M_PI)/4, (11*M_PI)/6,
                   0.0, M_PI/6, M_PI/4, M_PI/3, M_PI/2};

  for (int i = 0; i < num_components; i++) {
    decs[i] = dec0;
  }

  //Get the settings into a woden_settings_t struct
  woden_settings_t *woden_settings = malloc(sizeof(woden_settings_t));
  woden_settings->num_freqs = num_freqs;
  woden_settings->num_time_steps = num_times;
  woden_settings->ra0 = ra0;
  woden_settings->sdec0 = sinf(dec0);
  woden_settings->cdec0 = cosf(dec0);

  beam_settings_t *beam_settings = malloc(sizeof(beam_settings_t));
  beam_settings->beamtype = beamtype;

  /*********************************************************************
  Code used to generate the az / za and parallactic angles is below
  Stick into permanent arrays so test doesn't rely on erfa.
  **********************************************************************

  float *zas = malloc(num_components*num_times*sizeof(float));
  float *azs = malloc(num_components*num_times*sizeof(float));
  float *sin_para_angs = malloc(num_components*num_times*sizeof(float));
  float *cos_para_angs = malloc(num_components*num_times*sizeof(float));

  //Let's say the sources are stationary on the sky to keep the maths in the
  //test to a minimum
  double erfa_az, el, ha, para_angle;
  double lst = 0.0;
  double latitude = 0.0;

  for (size_t component = 0; component < num_components; component++) {
    for (size_t time_ind = 0; time_ind < num_times; time_ind++) {

      double ha = lst - ras[component];

      eraHd2ae(ha, (double)decs[component], latitude, &erfa_az, &el );

      azs[component*num_times + time_ind] = (float)erfa_az;
      zas[component*num_times + time_ind] = M_PI / 2. - (float)el;

      para_angle = eraHd2pa(ha, (double)decs[component], latitude);

      sin_para_angs[component*num_times + time_ind] = (float)sin(para_angle);
      cos_para_angs[component*num_times + time_ind] = (float)cos(para_angle);

      printf("%d %.8f %.8f %.8f %.8f\n",component*num_times + time_ind,
                                 (float)erfa_az, M_PI / 2. - (float)el,
                                 (float)sin(para_angle), (float)cos(para_angle));
    }
  }
  */

  float freqs[] = {100e+6, 200e+6};

  float *beam_has = malloc(num_components*num_times*sizeof(float));
  float *beam_decs = malloc(num_components*num_times*sizeof(float));

  //Make ha/dec coords if using the Gaussian beam
  if (beamtype == GAUSS_BEAM) {
    //Stick the Gaussian beam pointed at zenith
    //We're testing at latitude=zero
    beam_settings->gauss_ha = 0.0;
    beam_settings->gauss_sdec = 0.0;
    beam_settings->gauss_cdec = 1.0;
    beam_settings->beam_FWHM_rad = 80.0*DD2R;
    beam_settings->beam_ref_freq = 150e+6;

    for (size_t component = 0; component < num_components; component++) {
      for (size_t time_ind = 0; time_ind < num_times; time_ind++) {

        beam_has[component*num_times + time_ind] = ras[component];
        beam_decs[component*num_times + time_ind] = decs[component];

      }
    }
  }

  beam_settings->FEE_beam = malloc(sizeof(RTS_MWA_FEE_beam_t));
  //If FEE_BEAM, call the C code to interrogate the hdf5 file and set beam
  //things up
  if (beamtype == FEE_BEAM) {

    //Get a zenith pointing beam for normalisation purposes
    RTS_MWA_FEE_beam_t *FEE_beam_zenith = malloc(sizeof(RTS_MWA_FEE_beam_t));

    float base_middle_freq = 150e+6;
//
    float float_zenith_delays[16] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    //
    printf("\n\tSetting up the zenith FEE beam...");
    RTS_MWAFEEInit(mwa_fee_hdf5, base_middle_freq, FEE_beam_zenith, float_zenith_delays);
    printf(" done.\n");

    printf("\tSetting up the FEE beam...");
    RTS_MWAFEEInit(mwa_fee_hdf5, base_middle_freq, beam_settings->FEE_beam, float_zenith_delays);
    printf(" done.\n");

    printf("\tGetting FEE beam normalisation...");
    get_HDFBeam_normalisation(FEE_beam_zenith, beam_settings->FEE_beam);
    //Free the zenith pointing as done with it now
    free_FEE_primary_beam_from_GPU(FEE_beam_zenith);
    printf(" done.\n");

    printf("\tCopying the FEE beam across to the GPU...");
    copy_FEE_primary_beam_to_GPU(beam_settings->FEE_beam);
    printf(" done.\n");

  }

  //Output arrays
  float _Complex *primay_beam_J00 = calloc(num_beam_values, sizeof(float _Complex));
  float _Complex *primay_beam_J01 = calloc(num_beam_values, sizeof(float _Complex));
  float _Complex *primay_beam_J10 = calloc(num_beam_values, sizeof(float _Complex));
  float _Complex *primay_beam_J11 = calloc(num_beam_values, sizeof(float _Complex));

  float *ls = malloc(num_components*sizeof(float));
  float *ms = malloc(num_components*sizeof(float));
  float *ns = malloc(num_components*sizeof(float));

  //Run the CUDA code
  test_source_component_common(num_components,
             primay_beam_J00, primay_beam_J01,
             primay_beam_J10, primay_beam_J11,
             freqs, ls, ms, ns,
             ras, decs, azs, zas,
             sin_para_angs, cos_para_angs,
             beam_has, beam_decs,
             woden_settings,
             beam_settings);

  float l_expected[9] = {-1.0, -sqrt(3)/2.0, -sqrt(2)/2.0, -0.5,
                          0.0, 0.5, sqrt(2)/2.0, sqrt(3)/2.0, 1.0};
  float n_expected[9] = {0.0, 0.5, sqrt(2)/2.0, sqrt(3)/2.0,
                         1.0, sqrt(3)/2.0, sqrt(2)/2.0, 0.5, 0.0};

  //Check the l values match expectations
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(l_expected, ls, num_components);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(zeroes, ms, num_components);
  for (int i = 0; i < num_components; i++) {
    TEST_ASSERT_FLOAT_WITHIN(2e-7, n_expected[i], ns[i]);
  }

  //Depending on beamtype, check results match expectations
  if (beamtype == GAUSS_BEAM) {
    for (size_t output = 0; output < num_beam_values; output++) {
      // printf("%.7f %.5f %.5f %.5f %.5f %.5f %.7f %.5f\n",
          // creal(primay_beam_J00[output]), cimag(primay_beam_J00[output]),
          // creal(primay_beam_J01[output]), cimag(primay_beam_J01[output]),
          // creal(primay_beam_J10[output]), cimag(primay_beam_J10[output]),
          // creal(primay_beam_J11[output]), cimag(primay_beam_J11[output]));

          TEST_ASSERT_EQUAL_FLOAT(gauss_expected[output], creal(primay_beam_J00[output]));
          TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(primay_beam_J00[output]));
          TEST_ASSERT_EQUAL_FLOAT(0.0, creal(primay_beam_J01[output]));
          TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(primay_beam_J01[output]));
          TEST_ASSERT_EQUAL_FLOAT(0.0, creal(primay_beam_J10[output]));
          TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(primay_beam_J10[output]));
          TEST_ASSERT_EQUAL_FLOAT(gauss_expected[output], creal(primay_beam_J11[output]));
          TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(primay_beam_J11[output]));

    }
  }
  else if (beamtype == ANALY_DIPOLE) {
    for (size_t output = 0; output < num_beam_values; output++) {
          TEST_ASSERT_FLOAT_WITHIN(1e-7, analy_expec_J00[output], creal(primay_beam_J00[output]));
          TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(primay_beam_J00[output]));
          TEST_ASSERT_EQUAL_FLOAT(0.0, creal(primay_beam_J01[output]));
          TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(primay_beam_J01[output]));
          TEST_ASSERT_EQUAL_FLOAT(0.0, creal(primay_beam_J10[output]));
          TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(primay_beam_J10[output]));
          TEST_ASSERT_FLOAT_WITHIN(1e-7, analy_expec_J11[output], creal(primay_beam_J11[output]));
          TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(primay_beam_J11[output]));
    }
  }
  else if (beamtype == FEE_BEAM) {
    for (size_t output = 0; output < num_beam_values; output++) {

      TEST_ASSERT_FLOAT_WITHIN(1e-7, fee_expec_J00_re[output], creal(primay_beam_J00[output]));
      TEST_ASSERT_FLOAT_WITHIN(1e-7, fee_expec_J00_im[output], cimag(primay_beam_J00[output]));
      TEST_ASSERT_FLOAT_WITHIN(1e-7, fee_expec_J01_re[output], creal(primay_beam_J01[output]));
      TEST_ASSERT_FLOAT_WITHIN(1e-7, fee_expec_J01_im[output], cimag(primay_beam_J01[output]));
      TEST_ASSERT_FLOAT_WITHIN(1e-7, fee_expec_J10_re[output], creal(primay_beam_J10[output]));
      TEST_ASSERT_FLOAT_WITHIN(1e-7, fee_expec_J10_im[output], cimag(primay_beam_J10[output]));
      TEST_ASSERT_FLOAT_WITHIN(1e-7, fee_expec_J11_re[output], creal(primay_beam_J11[output]));
      TEST_ASSERT_FLOAT_WITHIN(1e-7, fee_expec_J11_im[output], cimag(primay_beam_J11[output]));
    }
  }
  else {
    //Don't need to calculate beam values for NO_BEAM, so these values
    //should still be zero, as we initialised the array with calloc
    for (size_t output = 0; output < num_beam_values; output++) {

          TEST_ASSERT_EQUAL_FLOAT(0.0, creal(primay_beam_J00[output]));
          TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(primay_beam_J00[output]));
          TEST_ASSERT_EQUAL_FLOAT(0.0, creal(primay_beam_J01[output]));
          TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(primay_beam_J01[output]));
          TEST_ASSERT_EQUAL_FLOAT(0.0, creal(primay_beam_J10[output]));
          TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(primay_beam_J10[output]));
          TEST_ASSERT_EQUAL_FLOAT(0.0, creal(primay_beam_J11[output]));
          TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(primay_beam_J11[output]));

    }
  }

  //Be free my beauties
  free(zeroes);
  free(decs);
  free(primay_beam_J00);
  free(primay_beam_J01);
  free(primay_beam_J10);
  free(primay_beam_J11);
  free(ls);
  free(ms);
  free(ns);

  // free(zas);
  // free(azs);
  // free(sin_para_angs);
  // free(cos_para_angs);

}

/*
This test checks source_component_common with beamtype=FEE_BEAM
*/
void test_source_component_common_ConstantDecFEEBeam(void) {
  //Look for environment variable MWA_FEE_HDF5, which should point
  //towards mwa_full_embedded_element_pattern.h5. If we can't find it,
  //can't run the tests, so don't
  char* mwa_fee_hdf5 = getenv("MWA_FEE_HDF5");

  if (mwa_fee_hdf5) {
    printf("MWA_FEE_HDF5: %s", mwa_fee_hdf5 );
    test_source_component_common_ConstantDecChooseBeams(FEE_BEAM, mwa_fee_hdf5);
  }
  else {
    printf("MWA_FEE_HDF5 not found - not running MWA_FEE beam test");
  }
}

/*
This test checks source_component_common with beamtype=ANALY_DIPOLE
*/
void test_source_component_common_ConstantDecAnalyBeam(void) {
  test_source_component_common_ConstantDecChooseBeams(ANALY_DIPOLE, " ");
}

/*
This test checks source_component_common with beamtype=GAUSS_BEAM
*/
void test_source_component_common_ConstantDecGaussBeam(void) {
  test_source_component_common_ConstantDecChooseBeams(GAUSS_BEAM, " ");
}

/*
This test checks source_component_common with beamtype=NO_BEAM
*/
void test_source_component_common_ConstantDecNoBeam(void) {
  test_source_component_common_ConstantDecChooseBeams(NO_BEAM, " ");
}

//Run the test with unity
int main(void)
{
    UNITY_BEGIN();
    RUN_TEST(test_source_component_common_ConstantDecNoBeam);
    RUN_TEST(test_source_component_common_ConstantDecGaussBeam);
    RUN_TEST(test_source_component_common_ConstantDecAnalyBeam);
    RUN_TEST(test_source_component_common_ConstantDecFEEBeam);

    return UNITY_END();
}
