#include <math.h>
#include <unity.h>
#include <stdlib.h>
#include <complex.h>
#include <erfa.h>

#include "constants.h"
#include "woden_struct_defs.h"
#include "woden_precision_defs.h"
#include "FEE_primary_beam.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

// void sincos(user_precision_t x, user_precision_t *sin, user_precision_t *cos);

//External CUDA code we're linking in
extern void test_source_component_common(int num_components,
           user_precision_complex_t *primay_beam_J00, user_precision_complex_t *primay_beam_J01,
           user_precision_complex_t *primay_beam_J10, user_precision_complex_t *primay_beam_J11,
           double *freqs, double *ls, double *ms, double *ns,
           double *ras, double *decs, user_precision_t *azs, user_precision_t *zas,
           user_precision_t *sin_para_angs, user_precision_t *cos_para_angs,
           double *beam_has, double *beam_decs,
           woden_settings_t *woden_settings,
           beam_settings_t *beam_settings);

extern void get_HDFBeam_normalisation(RTS_MWA_FEE_beam_t *FEE_beam_zenith,
                RTS_MWA_FEE_beam_t *FEE_beam);

extern void free_FEE_primary_beam_from_GPU(RTS_MWA_FEE_beam_t *primary_beam);

extern void copy_FEE_primary_beam_to_GPU(RTS_MWA_FEE_beam_t *FEE_beam);


/*******************************************************************************
Here are a bunch of predefined arrays to be used
Some are inputs for tests, others are expected outputs
*******************************************************************************/

user_precision_t azs[] = {4.71238899, 4.71238899, 4.71238899, 4.71238899, 4.71238899,
               4.71238899, 4.71238899, 4.71238899, 4.71238899, 4.71238899,
               4.71238899, 4.71238899, 0.00000000, 0.00000000, 0.00000000,
               1.57079637, 1.57079637, 1.57079637, 1.57079637, 1.57079637,
               1.57079637, 1.57079637, 1.57079637, 1.57079637, 1.57079637,
               1.57079637, 1.57079637};
user_precision_t zas[] = {1.57079631, 1.57079631, 1.57079631, 1.04719766, 1.04719766,
               1.04719766, 0.78539832, 0.78539832, 0.78539832, 0.52359898,
               0.52359898, 0.52359898, 0.00000000, 0.00000000, 0.00000000,
               0.52359875, 0.52359875, 0.52359875, 0.78539820, 0.78539820,
               0.78539820, 1.04719760, 1.04719760, 1.04719760, 1.57079637,
               1.57079637, 1.57079637};

user_precision_t sin_para_angs[] = {1.00000000, 1.00000000, 1.00000000, 1.00000000,
                         1.00000000, 1.00000000, 1.00000000, 1.00000000,
                         1.00000000, 1.00000000, 1.00000000, 1.00000000,
                         0.00000000, 0.00000000, 0.00000000, -1.00000000,
                         -1.00000000, -1.00000000, -1.00000000, -1.00000000,
                         -1.00000000, -1.00000000, -1.00000000, -1.00000000,
                         -1.00000000, -1.00000000, -1.00000000};

user_precision_t cos_para_angs[] = {0.00000000, 0.00000000, 0.00000000, 0.00000000,
                         0.00000000, 0.00000000, 0.00000000, 0.00000000,
                         0.00000000, 0.00000000, 0.00000000, 0.00000000,
                         1.00000000, 1.00000000, 1.00000000, 0.00000000,
                         0.00000000, 0.00000000, 0.00000000, 0.00000000,
                         0.00000000, 0.00000000, 0.00000000, 0.00000000,
                         0.00000000, 0.00000000, 0.00000000};

double analy_expec_J00[] = { 0.00000002, 0.52576818, 0.73128019,
  0.88075472, 1.00000000, 0.88075483, 0.73128027, 0.52576824, -0.00000005,
  0.00000002, 0.61822932, 0.81629589, 0.93152102, 1.00000000, 0.93152109,
  0.81629596, 0.61822938, -0.00000006, 0.00000002, 0.52576818, 0.73128019,
  0.88075472, 1.00000000, 0.88075483, 0.73128027, 0.52576824, -0.00000005,
  0.00000002, 0.61822932, 0.81629589, 0.93152102, 1.00000000, 0.93152109,
  0.81629596, 0.61822938, -0.00000006, 0.00000002, 0.52576818, 0.73128019,
  0.88075472, 1.00000000, 0.88075483, 0.73128027, 0.52576824, -0.00000005,
  0.00000002, 0.61822932, 0.81629589, 0.93152102, 1.00000000, 0.93152109,
  0.81629596, 0.61822938, -0.00000006 };

double analy_expec_J11[] = { 0.00000000, 0.26288404, 0.51709310,
  0.76275588, 1.00000000, 0.76275607, 0.51709322, 0.26288410, -0.00000000,
  0.00000000, 0.30911460, 0.57720827, 0.80672078, 1.00000000, 0.80672094,
  0.57720839, 0.30911466, -0.00000000, 0.00000000, 0.26288404, 0.51709310,
  0.76275588, 1.00000000, 0.76275607, 0.51709322, 0.26288410, -0.00000000,
  0.00000000, 0.30911460, 0.57720827, 0.80672078, 1.00000000, 0.80672094,
  0.57720839, 0.30911466, -0.00000000, 0.00000000, 0.26288404, 0.51709310,
  0.76275588, 1.00000000, 0.76275607, 0.51709322, 0.26288410, -0.00000000,
  0.00000000, 0.30911460, 0.57720827, 0.80672078, 1.00000000, 0.80672094,
  0.57720839, 0.30911466, -0.00000000 };

double fee_expec_J00_re[] = { -0.00000375, -0.00007259, 0.00001992, 0.00008827,
  0.00025053, 0.00008904, 0.00002027, -0.00007286, -0.00000362, -0.00000375, -0.00007259,
  0.00001992, 0.00008827, 0.00025053, 0.00008904, 0.00002027, -0.00007286, -0.00000362,
  -0.00000375, -0.00007259, 0.00001992, 0.00008827, 0.00025053, 0.00008904, 0.00002027,
  -0.00007286, -0.00000362, -0.00000375, -0.00007259, 0.00001992, 0.00008827,
  0.00025053, 0.00008904, 0.00002027, -0.00007286, -0.00000362, -0.00000375,
  -0.00007259, 0.00001992, 0.00008827, 0.00025053, 0.00008904, 0.00002027, -0.00007286,
  -0.00000362, -0.00000375, -0.00007259, 0.00001992, 0.00008827, 0.00025053,
  0.00008904, 0.00002027, -0.00007286, -0.00000362 };

double fee_expec_J00_im[] = { 0.00000086, -0.00003368, -0.00001715,
  0.00002948, 0.00012354, 0.00002986, -0.00001680, -0.00003321, 0.00000057,
  0.00000086, -0.00003368, -0.00001715, 0.00002948, 0.00012354, 0.00002986,
  -0.00001680, -0.00003321, 0.00000057, 0.00000086, -0.00003368, -0.00001715,
  0.00002948, 0.00012354, 0.00002986, -0.00001680, -0.00003321, 0.00000057,
  0.00000086, -0.00003368, -0.00001715, 0.00002948, 0.00012354, 0.00002986,
  -0.00001680, -0.00003321, 0.00000057, 0.00000086, -0.00003368, -0.00001715,
  0.00002948, 0.00012354, 0.00002986, -0.00001680, -0.00003321, 0.00000057,
  0.00000086, -0.00003368, -0.00001715, 0.00002948, 0.00012354, 0.00002986,
  -0.00001680, -0.00003321, 0.00000057 };

double fee_expec_J01_re[] = { -0.00299196, 0.05242207, 0.21119083,
  0.09754611, 0.97041277, 0.09754551, 0.21118965, 0.05242191, -0.00299213,
  -0.00299196, 0.05242207, 0.21119083, 0.09754611, 0.97041277, 0.09754551,
  0.21118965, 0.05242191, -0.00299213, -0.00299196, 0.05242207, 0.21119083,
  0.09754611, 0.97041277, 0.09754551, 0.21118965, 0.05242191, -0.00299213,
  -0.00299196, 0.05242207, 0.21119083, 0.09754611, 0.97041277, 0.09754551,
  0.21118965, 0.05242191, -0.00299213, -0.00299196, 0.05242207, 0.21119083,
  0.09754611, 0.97041277, 0.09754551, 0.21118965, 0.05242191, -0.00299213,
  -0.00299196, 0.05242207, 0.21119083, 0.09754611, 0.97041277, 0.09754551,
  0.21118965, 0.05242191, -0.00299213 };

double fee_expec_J01_im[] = { -0.00004542, 0.00969854, 0.03443863,
  0.01937939, 0.24145197, 0.01937985, 0.03443851, 0.00969827, -0.00004550,
  -0.00004542, 0.00969854, 0.03443863, 0.01937939, 0.24145197, 0.01937985,
   0.03443851, 0.00969827, -0.00004550, -0.00004542, 0.00969854, 0.03443863,
   0.01937939, 0.24145197, 0.01937985, 0.03443851, 0.00969827, -0.00004550,
   -0.00004542, 0.00969854, 0.03443863, 0.01937939, 0.24145197, 0.01937985,
   0.03443851, 0.00969827, -0.00004550, -0.00004542, 0.00969854, 0.03443863,
    0.01937939, 0.24145197, 0.01937985, 0.03443851, 0.00969827, -0.00004550,
    -0.00004542, 0.00969854, 0.03443863, 0.01937939, 0.24145197, 0.01937985,
     0.03443851, 0.00969827, -0.00004550 };

double fee_expec_J10_re[] = {0.00169129, -0.02178988, -0.13555926,
  -0.07058298, -0.97021026, -0.07058250, -0.13556084, -0.02179252, 0.00168970,
  0.00169129, -0.02178988, -0.13555926, -0.07058298, -0.97021026, -0.07058250,
  -0.13556084, -0.02179252, 0.00168970, 0.00169129, -0.02178988, -0.13555926,
  -0.07058298, -0.97021026, -0.07058250, -0.13556084, -0.02179252, 0.00168970,
  0.00169129, -0.02178988, -0.13555926, -0.07058298, -0.97021026, -0.07058250,
  -0.13556084, -0.02179252, 0.00168970, 0.00169129, -0.02178988, -0.13555926,
  -0.07058298, -0.97021026, -0.07058250, -0.13556084, -0.02179252, 0.00168970,
  0.00169129, -0.02178988, -0.13555926, -0.07058298, -0.97021026, -0.07058250,
  -0.13556084, -0.02179252, 0.00168970  };

double fee_expec_J10_im[] = { 0.00086604, 0.01911517, -0.01697326,
  -0.03916865, -0.24226442, -0.03917003, -0.01697404, 0.01911435, 0.00086739,
  0.00086604, 0.01911517, -0.01697326, -0.03916865, -0.24226442, -0.03917003,
  -0.01697404, 0.01911435, 0.00086739, 0.00086604, 0.01911517, -0.01697326,
  -0.03916865, -0.24226442, -0.03917003, -0.01697404, 0.01911435, 0.00086739,
  0.00086604, 0.01911517, -0.01697326, -0.03916865, -0.24226442, -0.03917003,
  -0.01697404, 0.01911435, 0.00086739, 0.00086604, 0.01911517, -0.01697326,
  -0.03916865, -0.24226442, -0.03917003, -0.01697404, 0.01911435, 0.00086739,
   0.00086604, 0.01911517, -0.01697326, -0.03916865, -0.24226442, -0.03917003,
   -0.01697404, 0.01911435, 0.00086739 };

double fee_expec_J11_re[] = { -0.00000147, -0.00000903, -0.00001642,
  -0.00006570, -0.00024691, -0.00006504, -0.00001520, -0.00000766, -0.00000182,
  -0.00000147, -0.00000903, -0.00001642, -0.00006570, -0.00024691, -0.00006504,
  -0.00001520, -0.00000766, -0.00000182, -0.00000147, -0.00000903, -0.00001642,
  -0.00006570, -0.00024691, -0.00006504, -0.00001520, -0.00000766, -0.00000182,
  -0.00000147, -0.00000903, -0.00001642, -0.00006570, -0.00024691, -0.00006504,
  -0.00001520, -0.00000766, -0.00000182, -0.00000147, -0.00000903, -0.00001642,
  -0.00006570, -0.00024691, -0.00006504, -0.00001520, -0.00000766, -0.00000182,
  -0.00000147, -0.00000903, -0.00001642, -0.00006570, -0.00024691, -0.00006504,
  -0.00001520, -0.00000766, -0.00000182 };

double fee_expec_J11_im[] = { 0.00000424, 0.00007635, 0.00002101,
  -0.00010817, -0.00011722, -0.00010677, 0.00002321, 0.00007781, 0.00000384,
  0.00000424, 0.00007635, 0.00002101, -0.00010817, -0.00011722, -0.00010677,
   0.00002321, 0.00007781, 0.00000384, 0.00000424, 0.00007635, 0.00002101,
   -0.00010817, -0.00011722, -0.00010677, 0.00002321, 0.00007781, 0.00000384,
    0.00000424, 0.00007635, 0.00002101, -0.00010817, -0.00011722, -0.00010677,
    0.00002321, 0.00007781, 0.00000384, 0.00000424, 0.00007635, 0.00002101,
    -0.00010817, -0.00011722, -0.00010677, 0.00002321, 0.00007781, 0.00000384,
    0.00000424, 0.00007635, 0.00002101, -0.00010817, -0.00011722, -0.00010677,
    0.00002321, 0.00007781, 0.00000384 };

double TOL;


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
  // int num_baselines = 5;

  int num_beam_values = num_times*num_freqs*num_components;

  double ra0 = 0.0*DD2R;
  double dec0 = 0.0*DD2R;

  double *decs = malloc(num_components*sizeof(double));
  double *zeroes = calloc(num_components, sizeof(double));

  //Keep RA between 0 and 2*pi here but enter RAs that should return
  //negative l values
  double ras[9] = {(3*M_PI)/2, (5*M_PI)/3, (7*M_PI)/4, (11*M_PI)/6,
                   0.0, M_PI/6, M_PI/4, M_PI/3, M_PI/2};

  for (int i = 0; i < num_components; i++) {
    decs[i] = dec0;
  }

  //Get the settings into a woden_settings_t struct
  woden_settings_t *woden_settings = malloc(sizeof(woden_settings_t));
  woden_settings->num_freqs = num_freqs;
  woden_settings->num_time_steps = num_times;
  woden_settings->ra0 = ra0;
  woden_settings->sdec0 = sin(dec0);
  woden_settings->cdec0 = cos(dec0);

  beam_settings_t *beam_settings = malloc(sizeof(beam_settings_t));
  beam_settings->beamtype = beamtype;

  /*********************************************************************
  Code used to generate the az / za and parallactic angles is below
  Stick into permanent arrays so test doesn't rely on erfa.
  **********************************************************************

  user_precision_t *zas = malloc(num_components*num_times*sizeof(user_precision_t));
  user_precision_t *azs = malloc(num_components*num_times*sizeof(user_precision_t));
  user_precision_t *sin_para_angs = malloc(num_components*num_times*sizeof(user_precision_t));
  user_precision_t *cos_para_angs = malloc(num_components*num_times*sizeof(user_precision_t));

  //Let's say the sources are stationary on the sky to keep the maths in the
  //test to a minimum
  double erfa_az, el, ha, para_angle;
  double lst = 0.0;
  double latitude = 0.0;

  for (int component = 0; component < num_components; component++) {
    for (int time_ind = 0; time_ind < num_times; time_ind++) {

      double ha = lst - ras[component];

      eraHd2ae(ha, (double)decs[component], latitude, &erfa_az, &el );

      azs[component*num_times + time_ind] = (user_precision_t)erfa_az;
      zas[component*num_times + time_ind] = M_PI / 2. - (user_precision_t)el;

      para_angle = eraHd2pa(ha, (double)decs[component], latitude);

      sin_para_angs[component*num_times + time_ind] = (user_precision_t)sin(para_angle);
      cos_para_angs[component*num_times + time_ind] = (user_precision_t)cos(para_angle);

      printf("%d %.8f %.8f %.8f %.8f\n",component*num_times + time_ind,
                                 (user_precision_t)erfa_az, M_PI / 2. - (user_precision_t)el,
                                 (user_precision_t)sin(para_angle), (user_precision_t)cos(para_angle));
    }
  }
  */

  double freqs[] = {100e+6, 200e+6};

  double *beam_has = malloc(num_components*num_times*sizeof(double));
  double *beam_decs = malloc(num_components*num_times*sizeof(double));

  //Make ha/dec coords if using the Gaussian beam
  if (beamtype == GAUSS_BEAM) {
    //Stick the Gaussian beam pointed at zenith
    //We're testing at latitude=zero
    beam_settings->gauss_ha = 0.0;
    beam_settings->gauss_sdec = 0.0;
    beam_settings->gauss_cdec = 1.0;
    beam_settings->beam_FWHM_rad = 80.0*DD2R;
    beam_settings->beam_ref_freq = 150e+6;

    for (int component = 0; component < num_components; component++) {
      for (int time_ind = 0; time_ind < num_times; time_ind++) {

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

    double base_middle_freq = 150e+6;
//
    user_precision_t user_precision_t_zenith_delays[16] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    //
    printf("\n\tSetting up the zenith FEE beam...");
    RTS_MWAFEEInit(mwa_fee_hdf5, base_middle_freq, FEE_beam_zenith, user_precision_t_zenith_delays);
    printf(" done.\n");

    printf("\tSetting up the FEE beam...");
    RTS_MWAFEEInit(mwa_fee_hdf5, base_middle_freq, beam_settings->FEE_beam, user_precision_t_zenith_delays);
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

  // //Output arrays
  user_precision_complex_t *primay_beam_J00 = calloc(num_beam_values, sizeof(user_precision_complex_t));
  user_precision_complex_t *primay_beam_J01 = calloc(num_beam_values, sizeof(user_precision_complex_t));
  user_precision_complex_t *primay_beam_J10 = calloc(num_beam_values, sizeof(user_precision_complex_t));
  user_precision_complex_t *primay_beam_J11 = calloc(num_beam_values, sizeof(user_precision_complex_t));

  double *ls = malloc(num_components*sizeof(double));
  double *ms = malloc(num_components*sizeof(double));
  double *ns = malloc(num_components*sizeof(double));

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

  double l_expected[9] = {-1.0, -sqrt(3)/2.0, -sqrt(2)/2.0, -0.5,
                          0.0, 0.5, sqrt(2)/2.0, sqrt(3)/2.0, 1.0};
  double n_expected[9] = {0.0, 0.5, sqrt(2)/2.0, sqrt(3)/2.0,
                         1.0, sqrt(3)/2.0, sqrt(2)/2.0, 0.5, 0.0};

  //Check the l values match expectations
  //Both FLOAT and DOUBLE use 64 bit here so same tolerance
  TOL = 1e-15;

  for (int i = 0; i < num_components; i++) {
    TEST_ASSERT_DOUBLE_WITHIN(TOL, l_expected[i], ls[i]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, ms[i]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, n_expected[i], ns[i]);
  }

  //Depending on beamtype, check results match expectations

  //For the Gaussian beam, the way the l,m coords are set up measns we can
  //analytically predict the values, so go through results and calcluate
  //expected values, and compare to what we got
  int beam_ind = 0;
  if (beamtype == GAUSS_BEAM) {

    #ifdef DOUBLE_PRECISION
      TOL = 1e-12;
    #else
      TOL = 1e-7;
    #endif

    double fwhm_lm = sin(beam_settings->beam_FWHM_rad);
    double beam_ref_freq = 150e+6;

    for (int time_ind = 0; time_ind < num_times; time_ind++) {
      for (int freq_ind = 0; freq_ind < num_freqs; freq_ind++) {
        for (int comp_ind = 0; comp_ind < num_components; comp_ind++) {

          double std = (fwhm_lm / FWHM_FACTOR) * (beam_ref_freq / freqs[freq_ind]);
          double exp_inside = (double)ls[comp_ind] / std;
          double estimate = exp(-0.5*exp_inside*exp_inside);

          //Only real gains have value in GAUSSIAN beam model
          TEST_ASSERT_DOUBLE_WITHIN(TOL, estimate,
                                   creal(primay_beam_J00[beam_ind]));

          TEST_ASSERT_DOUBLE_WITHIN(TOL, estimate,
                                   creal(primay_beam_J11[beam_ind]));

          //Everything else should be zero
          TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(primay_beam_J10[beam_ind]));
          TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(primay_beam_J11[beam_ind]));
          TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(primay_beam_J00[beam_ind]));
          TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, creal(primay_beam_J01[beam_ind]));
          TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(primay_beam_J01[beam_ind]));
          TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, creal(primay_beam_J10[beam_ind]));

          beam_ind += 1;

        }
      }
    }
  }
  else if (beamtype == ANALY_DIPOLE) {

    #ifdef DOUBLE_PRECISION
      TOL = 1e-7;
    #else
      TOL = 1e-6;
    #endif

    for (int output = 0; output < num_beam_values; output++) {

      TEST_ASSERT_DOUBLE_WITHIN(TOL, analy_expec_J00[output], creal(primay_beam_J00[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, analy_expec_J11[output], creal(primay_beam_J11[output]));

      TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(primay_beam_J00[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, creal(primay_beam_J01[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(primay_beam_J01[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, creal(primay_beam_J10[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(primay_beam_J10[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(primay_beam_J11[output]));
    }
  }
  else if (beamtype == FEE_BEAM) {

    #ifdef DOUBLE_PRECISION
      TOL = 1e-7;
    #else
      TOL = 3e-2;
    #endif

    for (int output = 0; output < num_beam_values; output++) {
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_J00_re[output], creal(primay_beam_J00[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_J00_im[output], cimag(primay_beam_J00[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_J01_re[output], creal(primay_beam_J01[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_J01_im[output], cimag(primay_beam_J01[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_J10_re[output], creal(primay_beam_J10[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_J10_im[output], cimag(primay_beam_J10[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_J11_re[output], creal(primay_beam_J11[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_J11_im[output], cimag(primay_beam_J11[output]));
    }
  }
  if (beamtype == NO_BEAM) {
    //Don't need to calculate beam values for NO_BEAM, these values
    //should still be zero, as we initialised the array with calloc
    //The function `source_components::get_beam_gains` which is called later
    //by the COMPONENT kernels assigns gains on one and leakage of zero
    //in the NO_BEAM case

    TOL = 1e-15;

    for (int output = 0; output < num_beam_values; output++) {

      TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, creal(primay_beam_J00[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(primay_beam_J00[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, creal(primay_beam_J01[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(primay_beam_J01[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, creal(primay_beam_J10[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(primay_beam_J10[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, creal(primay_beam_J11[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(primay_beam_J11[output]));

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
  free(beam_has);
  free(beam_decs);

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
