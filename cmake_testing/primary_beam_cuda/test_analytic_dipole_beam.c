#include <unity.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "constants.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

//External CUDA code we're linking in
extern void test_analytic_dipole_beam(int num_components,
             int num_time_steps, int num_freqs,
             float *azs, float *zas, float *freqs,
             float _Complex *analy_beam_X, float _Complex *analy_beam_Y);

#define UNITY_INCLUDE_FLOAT

/*
Test that the analytic dipole beam code returns the correct values, in the
correct order, for two time and frequency steps, with 25 directions on the sky
*/
void test_analytic_dipole_beam_GivesCorrectValues(void) {
  int num_freqs = 2;
  int num_times = 2;
  int num_az = 5;
  int num_za = 5;

  int num_components = num_az*num_za;
  int num_azza = num_az*num_za*num_times;
  int num_beam_values = num_azza*num_freqs;

  float az_inc = 2*M_PI / num_az;
  float za_inc = (M_PI/2) / num_za;

  float *azs = malloc(num_azza*sizeof(float));
  float *zas = malloc(num_azza*sizeof(float));

  //Function expects all az/za to include all time steps, and all
  //components. For each component, put az/za as they change with time
  //Here, we'll just stick each component to stay stationary on the sky,
  //so input the same az/za twice, one for each time step
  int azza_ind = 0;
  for (size_t az_ind = 0; az_ind < num_az; az_ind++) {
    for (size_t za_ind = 0; za_ind < num_za; za_ind++) {
      azs[azza_ind] = 0 + az_ind*az_inc;
      zas[azza_ind] = 0 + za_ind*za_inc;
      azza_ind += 1;

      azs[azza_ind] = 0 + az_ind*az_inc;
      zas[azza_ind] = 0 + za_ind*za_inc;
      azza_ind += 1;
    }
  }

  //Just test two frequencies
  float freqs[2] = {100e+6, 200e+6};

  float _Complex *analy_beam_X = malloc(num_beam_values*sizeof(float _Complex));
  float _Complex *analy_beam_Y = malloc(num_beam_values*sizeof(float _Complex));

  //Run the CUDA code
  test_analytic_dipole_beam(num_components,
               num_times, num_freqs,
               azs, zas, freqs,
               analy_beam_X, analy_beam_Y);

  //Expected real values for the X (north-south) dipole
  float expected_X[100] = { 1.0000000, 0.9103666, 0.6699219, 0.3610003,
      0.1014450, 1.0000000, 0.9528418, 0.8142948, 0.5946680, 0.3137854,
      1.0000000, 0.9268205, 0.7284526, 0.4643463, 0.2096881, 1.0000000,
      0.9268205, 0.7284527, 0.4643464, 0.2096881, 1.0000000, 0.9528418,
      0.8142948, 0.5946679, 0.3137854, 1.0000000, 0.9303051, 0.7234575,
      0.4162246, 0.1230724, 1.0000000, 0.9737105, 0.8793676, 0.6856377,
      0.3806824, 1.0000000, 0.9471194, 0.7866654, 0.5353799, 0.2543922,
      1.0000000, 0.9471194, 0.7866656, 0.5353801, 0.2543922, 1.0000000,
      0.9737105, 0.8793676, 0.6856377, 0.3806824, 1.0000000, 0.9103666,
      0.6699219, 0.3610003, 0.1014450, 1.0000000, 0.9528418, 0.8142948,
      0.5946680, 0.3137854, 1.0000000, 0.9268205, 0.7284526, 0.4643463,
      0.2096881, 1.0000000, 0.9268205, 0.7284527, 0.4643464, 0.2096881,
      1.0000000, 0.9528418, 0.8142948, 0.5946679, 0.3137854, 1.0000000,
      0.9303051, 0.7234575, 0.4162246, 0.1230724, 1.0000000, 0.9737105,
      0.8793676, 0.6856377, 0.3806824, 1.0000000, 0.9471194, 0.7866654,
      0.5353799, 0.2543922, 1.0000000, 0.9471194, 0.7866656, 0.5353801,
      0.2543922, 1.0000000, 0.9737105, 0.8793676, 0.6856377, 0.3806824 };

  //Expected real values for the Y (east-west) dipole
  float expected_Y[100] = { 1.0000000, 0.9572161, 0.8280690, 0.6141704,
      0.3282829, 1.0000000, 0.9149439, 0.6865985, 0.3922965, 0.1399981,
      1.0000000, 0.9412937, 0.7770780, 0.5402860, 0.2721978, 1.0000000,
      0.9412937, 0.7770780, 0.5402859, 0.2721978, 1.0000000, 0.9149439,
      0.6865985, 0.3922965, 0.1399981, 1.0000000, 0.9781806, 0.8942425,
      0.7081236, 0.3982707, 1.0000000, 0.9349827, 0.7414668, 0.4523083,
      0.1698447, 1.0000000, 0.9619095, 0.8391768, 0.6229365, 0.3302286,
      1.0000000, 0.9619095, 0.8391768, 0.6229365, 0.3302286, 1.0000000,
      0.9349827, 0.7414667, 0.4523083, 0.1698447, 1.0000000, 0.9572161,
      0.8280690, 0.6141704, 0.3282829, 1.0000000, 0.9149439, 0.6865985,
      0.3922965, 0.1399981, 1.0000000, 0.9412937, 0.7770780, 0.5402860,
      0.2721978, 1.0000000, 0.9412937, 0.7770780, 0.5402859, 0.2721978,
      1.0000000, 0.9149439, 0.6865985, 0.3922965, 0.1399981, 1.0000000,
      0.9781806, 0.8942425, 0.7081236, 0.3982707, 1.0000000, 0.9349827,
      0.7414668, 0.4523083, 0.1698447, 1.0000000, 0.9619095, 0.8391768,
      0.6229365, 0.3302286, 1.0000000, 0.9619095, 0.8391768, 0.6229365,
      0.3302286, 1.0000000, 0.9349827, 0.7414667, 0.4523083, 0.1698447 };

  //Separate real and imaginary components
  float *for_testing_X_re = malloc(num_beam_values*sizeof(float));
  float *for_testing_Y_re = malloc(num_beam_values*sizeof(float));
  float *for_testing_X_im = malloc(num_beam_values*sizeof(float));
  float *for_testing_Y_im = malloc(num_beam_values*sizeof(float));
  float *zeroes = (float*)calloc(num_beam_values, sizeof(float));

  for (size_t i = 0; i < num_beam_values; i++) {
    for_testing_X_re[i] = creal(analy_beam_X[i]);
    for_testing_X_im[i] = cimag(analy_beam_X[i]);
    for_testing_Y_re[i] = creal(analy_beam_Y[i]);
    for_testing_Y_im[i] = cimag(analy_beam_Y[i]);
  }

  //Test the arrays make sense. Analytic dipole beam is purely real, so
  //imaginary should be zero.
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(expected_X, for_testing_X_re, num_beam_values);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(zeroes, for_testing_X_im, num_beam_values);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(expected_Y, for_testing_Y_re, num_beam_values);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(zeroes, for_testing_Y_im, num_beam_values);

  //****************************************************************************
  //UNCOMMENT TO PRINT STUFF OUT AND PLOT USING PYTHON

  // float *for_printing_X = malloc(num_beam_values*sizeof(float));
  // float *for_printing_Y = malloc(num_beam_values*sizeof(float));
  // float *for_printing_az = malloc(num_beam_values*sizeof(float));
  // float *for_printing_za = malloc(num_beam_values*sizeof(float));

  // for (size_t iCoord = 0; iCoord < num_azza; iCoord++) {
  //   for (size_t freq = 0; freq < num_freqs; freq++) {
  //
  //     int component = (int)floorf((float)iCoord / (float)num_times);
  //     int time_ind = iCoord - component*num_times;
  //     int beam_ind = num_freqs*time_ind*num_components + (num_components*freq) + component;
  //
  //     for_printing_X[beam_ind] = creal(analy_beam_X[beam_ind]);
  //     for_printing_Y[beam_ind] = creal(analy_beam_Y[beam_ind]);
  //     for_printing_az[beam_ind] = azs[iCoord];
  //     for_printing_za[beam_ind] = zas[iCoord];
  //   }
  // }
  //
  // for (size_t i = 0; i < num_beam_values; i++) {
  //   printf("%.7f %.7f %.7f 0 %.7f 0\n",for_printing_az[i], for_printing_za[i],
  //                                      for_printing_X[i], for_printing_Y[i] );
  // }
  //****************************************************************************
}

//Run the test with unity
int main(void)
{
    UNITY_BEGIN();
    RUN_TEST(test_analytic_dipole_beam_GivesCorrectValues);
    return UNITY_END();
}
