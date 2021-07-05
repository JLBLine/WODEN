#include <unity.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "constants.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

//External CUDA code we're linking in
extern void test_kern_extrap_stokes(int num_extrap_freqs, int num_components,
           float *extrap_wavelengths, float *ref_freqs, float *SIs,
           float *ref_stokesI, float *ref_stokesQ,
           float *ref_stokesU, float *ref_stokesV,
           float *flux_I, float *flux_Q,
           float *flux_U, float *flux_V);

#define UNITY_INCLUDE_FLOAT

/*
Test that the linear SI flux extrapolation code works correctly
*/
void test_kern_extrap_stokes_GivesCorrectValues(void) {

  //Set up some test condition inputs
  int num_extrap_freqs = 5;
  int num_components = 5;

  float extrap_wavelengths[5] = {VELC/50e+6, VELC/100e+6, VELC/150e+6,
                                 VELC/200e+6, VELC/250e+6};
  float ref_freqs[5] = {50e+6, 100e+6, 150e+6, 200e+6, 250e+6};
  float SIs[5] = {0.0, -0.8, 0.5, -0.5, 1.0};

  float ref_stokesI[5] = {1.0, 1.0, 1.0, 1.0, 1.0};
  float ref_stokesQ[5] = {0.0, 0.0, 1.0, 0.0, 0.0};
  float ref_stokesU[5] = {0.0, 0.0, 0.0, 1.0, 0.0};
  float ref_stokesV[5] = {0.0, 0.0, 0.0, 0.0, 1.0};

  //Space for outputs
  float *flux_I = malloc(num_extrap_freqs*num_components*sizeof(float));
  float *flux_Q = malloc(num_extrap_freqs*num_components*sizeof(float));
  float *flux_U = malloc(num_extrap_freqs*num_components*sizeof(float));
  float *flux_V = malloc(num_extrap_freqs*num_components*sizeof(float));

  //Run the CUDA code
  test_kern_extrap_stokes(num_extrap_freqs, num_components,
                           extrap_wavelengths, ref_freqs, SIs,
                           ref_stokesI, ref_stokesQ,
                           ref_stokesU, ref_stokesV,
                           flux_I, flux_Q,
                           flux_U, flux_V);

  //Make some expected value arrays
  float *expec_flux_I = malloc(num_extrap_freqs*num_components*sizeof(float));
  float *expec_flux_Q = malloc(num_extrap_freqs*num_components*sizeof(float));
  float *expec_flux_U = malloc(num_extrap_freqs*num_components*sizeof(float));
  float *expec_flux_V = malloc(num_extrap_freqs*num_components*sizeof(float));

  //Fill values with what should have been found
  int ind = 0;
  for (size_t comp = 0; comp < num_components; comp++) {
    for (size_t extrap = 0; extrap < num_extrap_freqs; extrap++) {

      float extrap_freq = VELC / extrap_wavelengths[extrap];
      float flux_ratio = powf(extrap_freq / ref_freqs[comp], SIs[comp]);

      expec_flux_I[ind] = ref_stokesI[comp] * flux_ratio;
      expec_flux_Q[ind] = ref_stokesQ[comp] * flux_ratio;
      expec_flux_U[ind] = ref_stokesU[comp] * flux_ratio;
      expec_flux_V[ind] = ref_stokesV[comp] * flux_ratio;

      ind ++;
    }
  }

  //Check the two are equal
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_flux_I, flux_I, num_extrap_freqs*num_components);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_flux_Q, flux_Q, num_extrap_freqs*num_components);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_flux_U, flux_U, num_extrap_freqs*num_components);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_flux_V, flux_V, num_extrap_freqs*num_components);

  //Be free my beauties
  free(flux_I);
  free(flux_Q);
  free(flux_U);
  free(flux_V);
  free(expec_flux_I);
  free(expec_flux_Q);
  free(expec_flux_U);
  free(expec_flux_V);

}

//Run the test with unity
int main(void)
{
    UNITY_BEGIN();
    RUN_TEST(test_kern_extrap_stokes_GivesCorrectValues);
    return UNITY_END();
}
