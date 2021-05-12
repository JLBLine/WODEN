#include <math.h>
#include <unity.h>
#include <stdlib.h>
#include <complex.h>
#include <erfa.h>

#include "constants.h"
#include "woden_struct_defs.h"

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
           beam_settings_t beam_settings,
           RTS_MWA_FEE_beam_t *FEE_beam);

#define UNITY_INCLUDE_FLOAT

/*
Test that l,m,n and beam values are calculated correctly by `source_component_common`
for a constant Dec=0, latitude=0, and all beam types
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


  woden_settings_t *woden_settings = malloc(sizeof(woden_settings_t));
  woden_settings->num_freqs = num_freqs;
  woden_settings->num_time_steps = num_times;
  woden_settings->ra0 = ra0;
  woden_settings->sdec0 = sinf(dec0);
  woden_settings->cdec0 = cosf(dec0);

  beam_settings_t beam_settings;
  beam_settings.beamtype = beamtype;


  RTS_MWA_FEE_beam_t *FEE_beam = malloc(sizeof(RTS_MWA_FEE_beam_t));

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

  float _Complex *primay_beam_J00 = calloc(num_beam_values, sizeof(float _Complex));
  float _Complex *primay_beam_J01 = calloc(num_beam_values, sizeof(float _Complex));
  float _Complex *primay_beam_J10 = calloc(num_beam_values, sizeof(float _Complex));
  float _Complex *primay_beam_J11 = calloc(num_beam_values, sizeof(float _Complex));

  float *ls = malloc(num_components*sizeof(float));
  float *ms = malloc(num_components*sizeof(float));
  float *ns = malloc(num_components*sizeof(float));

  float l_expected[9] = {-1.0, -sqrt(3)/2.0, -sqrt(2)/2.0, -0.5,
                          0.0, 0.5, sqrt(2)/2.0, sqrt(3)/2.0, 1.0};
  float n_expected[9] = {0.0, 0.5, sqrt(2)/2.0, sqrt(3)/2.0,
                         1.0, sqrt(3)/2.0, sqrt(2)/2.0, 0.5, 0.0};

  float freqs[] = {100e+6, 200e+6};

  float *beam_has = malloc(num_components*num_times*sizeof(float));
  float *beam_decs = malloc(num_components*num_times*sizeof(float));

  if (beamtype == GAUSS_BEAM) {
    //Stick the Gaussian beam pointed at zenith
    //We're testing at latitude=zero
    beam_settings.gauss_ha = 0.0;
    beam_settings.gauss_sdec = 0.0;
    beam_settings.gauss_cdec = 1.0;
    beam_settings.beam_FWHM_rad = 80.0*DD2R;
    beam_settings.beam_ref_freq = 150e+6;

    for (size_t component = 0; component < num_components; component++) {
      for (size_t time_ind = 0; time_ind < num_times; time_ind++) {

        beam_has[component*num_times + time_ind] = ras[component];
        beam_decs[component*num_times + time_ind] = decs[component];

      }
    }
  }

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
  // //Run the CUDA code
  test_source_component_common(num_components,
             primay_beam_J00, primay_beam_J01,
             primay_beam_J10, primay_beam_J11,
             freqs, ls, ms, ns,
             ras, decs, azs, zas,
             sin_para_angs, cos_para_angs,
             beam_has, beam_decs,
             woden_settings,
             beam_settings,
             FEE_beam);

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


  //Check the l values match expectations
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(l_expected, ls, num_components);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(zeroes, ms, num_components);
  for (int i = 0; i < num_components; i++) {
    TEST_ASSERT_FLOAT_WITHIN(2e-7, n_expected[i], ns[i]);
  }

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
      // printf("%.7f %.5f %.5f %.5f %.5f %.5f %.7f %.5f\n",
      //     creal(primay_beam_J00[output]), cimag(primay_beam_J00[output]),
      //     creal(primay_beam_J01[output]), cimag(primay_beam_J01[output]),
      //     creal(primay_beam_J10[output]), cimag(primay_beam_J10[output]),
      //     creal(primay_beam_J11[output]), cimag(primay_beam_J11[output]));

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
  //No cross-pols so multiply expected values by 2.0
  else if (beamtype == FEE_BEAM) {
    for (size_t output = 0; output < num_beam_values; output++) {
      // printf("%.7f %.5f %.5f %.5f %.5f %.5f %.7f %.5f\n",
      //     creal(primay_beam_J00[output]), cimag(primay_beam_J00[output]),
      //     creal(primay_beam_J01[output]), cimag(primay_beam_J01[output]),
      //     creal(primay_beam_J10[output]), cimag(primay_beam_J10[output]),
      //     creal(primay_beam_J11[output]), cimag(primay_beam_J11[output]));
    }
  }
  else {
    //Don't need to calculat beam values for NO_BEAM, so these values
    //should still be zero, as we initialised the array with calloc
    for (size_t output = 0; output < num_beam_values; output++) {
      // printf("%.1f %.1f %.1f %.1f %.1f %.1f %.1f\n",
      //     creal(primay_beam_J00[output]), cimag(primay_beam_J00[output]),
      //     creal(primay_beam_J01[output]), cimag(primay_beam_J01[output]),
      //     creal(primay_beam_J10[output]), cimag(primay_beam_J10[output]),
      //     creal(primay_beam_J11[output]), cimag(primay_beam_J11[output]));

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
This test checks varying the measurement equation with beamtype=FEE_BEAM
*/
void test_source_component_common_ConstantDecFEEBeam(void) {
  char* mwa_fee_hdf5 = getenv("MWA_FEE_HDF5");

  if (mwa_fee_hdf5) {
    printf("MWA_FEE_HDF5: %s", mwa_fee_hdf5 );
  }
  else {
    printf("MWA_FEE_HDF5 not found - not running MWA_FEE beam test");
  }


  test_source_component_common_ConstantDecChooseBeams(FEE_BEAM, mwa_fee_hdf5);
}

/*
This test checks varying the measurement equation with beamtype=ANALY_DIPOLE
*/
void test_source_component_common_ConstantDecAnalyBeam(void) {
  test_source_component_common_ConstantDecChooseBeams(ANALY_DIPOLE, " ");
}

/*
This test checks varying the measurement equation with beamtype=GAUSS_BEAM
*/
void test_source_component_common_ConstantDecGaussBeam(void) {
  test_source_component_common_ConstantDecChooseBeams(GAUSS_BEAM, " ");
}

/*
This test checks varying the measurement equation with beamtype=NO_BEAM
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
