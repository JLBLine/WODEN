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


void test_analytic_dipole_beam_GivesCorrectValues(void) {
  int num_freqs = 1;
  int num_times = 1;
  int num_az = 5;
  int num_za = 5;

  int num_azza = num_az*num_za;

  int num_beam_values = num_azza*num_freqs*num_times;

  float az_inc = 2*M_PI / num_az;
  float za_inc = M_PI / num_za;

  float *azs = malloc(num_beam_values*sizeof(float));
  float *zas = malloc(num_beam_values*sizeof(float));

  int azza_ind = 0;
  for (size_t az_ind = 0; az_ind < num_az; az_ind++) {
    for (size_t za_ind = 0; za_ind < num_za; za_ind++) {
      azs[azza_ind] = 0 + az_ind*az_inc;
      zas[azza_ind] = 0 + za_ind*za_inc;

      azza_ind += 1;

    }
  }

  float start_freq = 100e+6;
  float freq_inc = 25e+6;

  float *freqs = malloc(num_freqs*sizeof(float));

  for (size_t i = 0; i < num_freqs; i++) {
    freqs[i] = start_freq + i*freq_inc;
  }

  float _Complex *analy_beam_X = malloc(num_beam_values*sizeof(float _Complex));
  float _Complex *analy_beam_Y = malloc(num_beam_values*sizeof(float _Complex));

  test_analytic_dipole_beam(num_azza,
               num_times, num_freqs,
               azs, zas, freqs,
               analy_beam_X, analy_beam_Y);

  for (size_t i = 0; i < num_beam_values; i++) {
    printf("%.5f %.5f %.5f %.5f %.5f %.5f\n",azs[i], zas[i],
                                   creal(analy_beam_X[i]), cimag(analy_beam_X[i]),
                                   creal(analy_beam_Y[i]), cimag(analy_beam_Y[i]) );
  }

  //TODO figure out what expected values are, do comparison
  //Run the same __host__ function as a comparison??
  //Or just compare to a list of numbers?
  //
  // for (size_t i = 0; i < num_components; i++) {
  //
  //   float std = (fwhm_lm / FWHM_FACTOR);
  //   float exp_inside = beam_ls[i] / std;
  //   float estimate = expf(-0.5*exp_inside*exp_inside);
  //
  //   TEST_ASSERT_FLOAT_WITHIN(1e-7, estimate, creal(primay_beam_J00[i]));
  //   TEST_ASSERT_FLOAT_WITHIN(1e-7, estimate, creal(primay_beam_J11[i]));
  //   TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(primay_beam_J00[i]));
  //   TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(primay_beam_J11[i]));

}


int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_analytic_dipole_beam_GivesCorrectValues);

    return UNITY_END();
}
