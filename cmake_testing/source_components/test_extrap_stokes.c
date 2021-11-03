#include <unity.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "constants.h"
#include "woden_precision_defs.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

//External CUDA code we're linking in
extern void test_kern_extrap_stokes(int num_extrap_freqs, int num_components,
           user_precision_t *extrap_wavelengths, user_precision_t *ref_freqs, user_precision_t *SIs,
           user_precision_t *ref_stokesI, user_precision_t *ref_stokesQ,
           user_precision_t *ref_stokesU, user_precision_t *ref_stokesV,
           user_precision_t *flux_I, user_precision_t *flux_Q,
           user_precision_t *flux_U, user_precision_t *flux_V);

#define UNITY_INCLUDE_FLOAT

/*
Test that the linear SI flux extrapolation code works correctly
*/
void test_kern_extrap_stokes_GivesCorrectValues(void) {

  //Set up some test condition inputs
  int num_extrap_freqs = 5;
  int num_components = 5;

  user_precision_t extrap_wavelengths[5] = {VELC/50e+6, VELC/100e+6, VELC/150e+6,
                                 VELC/200e+6, VELC/250e+6};
  user_precision_t ref_freqs[5] = {50e+6, 100e+6, 150e+6, 200e+6, 250e+6};
  user_precision_t SIs[5] = {0.0, -0.8, 0.5, -0.5, 1.0};

  user_precision_t ref_stokesI[5] = {1.0, 1.0, 1.0, 1.0, 1.0};
  user_precision_t ref_stokesQ[5] = {0.0, 0.0, 1.0, 0.0, 0.0};
  user_precision_t ref_stokesU[5] = {0.0, 0.0, 0.0, 1.0, 0.0};
  user_precision_t ref_stokesV[5] = {0.0, 0.0, 0.0, 0.0, 1.0};

  //Space for outputs
  user_precision_t *flux_I = malloc(num_extrap_freqs*num_components*sizeof(user_precision_t));
  user_precision_t *flux_Q = malloc(num_extrap_freqs*num_components*sizeof(user_precision_t));
  user_precision_t *flux_U = malloc(num_extrap_freqs*num_components*sizeof(user_precision_t));
  user_precision_t *flux_V = malloc(num_extrap_freqs*num_components*sizeof(user_precision_t));

  //Run the CUDA code
  test_kern_extrap_stokes(num_extrap_freqs, num_components,
                           extrap_wavelengths, ref_freqs, SIs,
                           ref_stokesI, ref_stokesQ,
                           ref_stokesU, ref_stokesV,
                           flux_I, flux_Q,
                           flux_U, flux_V);

  //Make some expected value arrays
  user_precision_t *expec_flux_I = malloc(num_extrap_freqs*num_components*sizeof(user_precision_t));
  user_precision_t *expec_flux_Q = malloc(num_extrap_freqs*num_components*sizeof(user_precision_t));
  user_precision_t *expec_flux_U = malloc(num_extrap_freqs*num_components*sizeof(user_precision_t));
  user_precision_t *expec_flux_V = malloc(num_extrap_freqs*num_components*sizeof(user_precision_t));

  //Fill values with what should have been found
  int ind = 0;
  for (int comp = 0; comp < num_components; comp++) {
    for (int extrap = 0; extrap < num_extrap_freqs; extrap++) {

      user_precision_t extrap_freq = VELC / extrap_wavelengths[extrap];
      user_precision_t flux_ratio = pow(extrap_freq / ref_freqs[comp], SIs[comp]);

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
