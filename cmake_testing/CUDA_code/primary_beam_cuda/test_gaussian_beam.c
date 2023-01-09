#include <unity.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "constants.h"
#include "woden_precision_defs.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

//External CUDA code we're linking in
extern void test_kern_gaussian_beam(double *beam_ls, double *beam_ms,
           double beam_ref_freq, double *freqs,
           user_precision_t fwhm_lm, user_precision_t cos_theta, user_precision_t sin_theta, user_precision_t sin_2theta,
           int num_freqs, int num_time_steps, int num_components,
           user_precision_complex_t *primay_beam_J00, user_precision_complex_t *primay_beam_J11);

#ifdef DOUBLE_PRECISION
  double TOL = 1e-16;
#else
  double TOL = 1e-10;
#endif

/*
This function calls test_kern_gaussian_beam. It calculates a set of l,m input
coords, either varying l or m, keeping the other coord set to zero. This way
it tests a strip in l or m
*/
void get_1D_gaussian_values(double *beam_ls, double *beam_ms,
                            user_precision_complex_t *primay_beam_J00,
                            user_precision_complex_t *primay_beam_J11,
                            double *freqs, double beam_ref_freq,
                            int vary_l, user_precision_t fwhm_lm,
                            int num_freqs, int num_times, int num_components){

  // user_precision_t beam_ref_freq = 150e+6;

  int l_range = 2;
  double l_inc = l_range / num_components;

  //TODO currently hardcoded to have beam position angle = 0.
  //Should this change with az/za?
  user_precision_t cos_theta = 1.0;
  user_precision_t sin_theta = 0.0;
  user_precision_t sin_2theta = 0.0;

  if (vary_l == 0) {
    for (int i = 0; i < num_components; i++) {
      beam_ls[i] = 0.0;
      beam_ms[i] = -1 + i*l_inc;
    }
  }
  else {
    for (int i = 0; i < num_components; i++) {
      beam_ls[i] = -1 + i*l_inc;
      beam_ms[i] = 0.0;
    }
  }

  double freq_inc = 25e+5;

  for (int i = 0; i < num_freqs; i++) {
    freqs[i] = beam_ref_freq + i*freq_inc;
  }

  test_kern_gaussian_beam(beam_ls, beam_ms,
                          beam_ref_freq, freqs,
                          fwhm_lm, cos_theta, sin_theta, sin_2theta,
                          num_freqs, num_times, num_components,
                          primay_beam_J00, primay_beam_J11);
}

void test_GaussBeam_GivesCorrectlValues(void) {
  int num_freqs = 1;
  int num_times = 1;
  int num_components = 100;
  int num_beam_values = num_freqs*num_components*num_times;

  //Set full width half max to 20 deg, convert to l,m coord
  user_precision_t fwhm = 20.0;
  user_precision_t fwhm_lm = sin(fwhm*DD2R);
  double beam_ref_freq = 150e+6;

  double *beam_ls = malloc(num_beam_values*sizeof(double));
  double *beam_ms = malloc(num_beam_values*sizeof(double));

  user_precision_complex_t *primay_beam_J00 = malloc(num_beam_values*sizeof(user_precision_complex_t));
  user_precision_complex_t *primay_beam_J11 = malloc(num_beam_values*sizeof(user_precision_complex_t));

  double *freqs = malloc(num_freqs*sizeof(double));

  int vary_l = 1;
  get_1D_gaussian_values(beam_ls, beam_ms,
                         primay_beam_J00, primay_beam_J11,
                         freqs, beam_ref_freq,
                         vary_l, fwhm_lm,
                         num_freqs, num_times, num_components);

  for (int i = 0; i < num_components; i++) {

    double std = (fwhm_lm / FWHM_FACTOR);
    double exp_inside = beam_ls[i] / std;
    double estimate = exp(-0.5*exp_inside*exp_inside);

    TEST_ASSERT_DOUBLE_WITHIN(TOL, estimate, creal(primay_beam_J00[i]));
    TEST_ASSERT_DOUBLE_WITHIN(TOL, estimate, creal(primay_beam_J11[i]));
    TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(primay_beam_J00[i]));
    TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(primay_beam_J11[i]));

  }

  free(beam_ls);
  free(beam_ms);
  free(primay_beam_J00);
  free(primay_beam_J11);
  free(freqs);
}

void test_GaussBeam_GivesCorrectmValues(void) {
  int num_freqs = 1;
  int num_times = 1;
  int num_components = 100;
  int num_beam_values = num_freqs*num_components*num_times;

  //Set full width half max to 20 deg, convert to l,m coord
  user_precision_t fwhm = 20.0;
  user_precision_t fwhm_lm = sin(fwhm*DD2R);
  double beam_ref_freq = 150e+6;

  double *beam_ls = malloc(num_beam_values*sizeof(double));
  double *beam_ms = malloc(num_beam_values*sizeof(double));

  user_precision_complex_t *primay_beam_J00 = malloc(num_beam_values*sizeof(user_precision_complex_t));
  user_precision_complex_t *primay_beam_J11 = malloc(num_beam_values*sizeof(user_precision_complex_t));

  double *freqs = malloc(num_freqs*sizeof(double));

  int vary_l = 0;
  get_1D_gaussian_values(beam_ls, beam_ms,
                         primay_beam_J00, primay_beam_J11,
                         freqs, beam_ref_freq,
                         vary_l, fwhm_lm,
                         num_freqs, num_times, num_components);

  for (int i = 0; i < num_components; i++) {

    double std = (fwhm_lm / FWHM_FACTOR);
    double exp_inside = beam_ms[i] / std;
    double estimate = exp(-0.5*exp_inside*exp_inside);

    TEST_ASSERT_DOUBLE_WITHIN(TOL, estimate, creal(primay_beam_J00[i]));
    TEST_ASSERT_DOUBLE_WITHIN(TOL, estimate, creal(primay_beam_J11[i]));
    TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(primay_beam_J00[i]));
    TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(primay_beam_J11[i]));

  }

  free(beam_ls);
  free(beam_ms);
  free(primay_beam_J00);
  free(primay_beam_J11);
  free(freqs);
}


void test_GaussBeam_GivesCorrectlValuesByFreq(void) {
  int num_freqs = 5;
  int num_times = 1;
  int num_components = 100;
  int num_beam_values = num_freqs*num_components*num_times;

  //Set full width half max to 20 deg, convert to l,m coord
  user_precision_t fwhm = 20.0;
  user_precision_t fwhm_lm = sin(fwhm*DD2R);
  double beam_ref_freq = 150e+6;

  double *beam_ls = malloc(num_beam_values*sizeof(double));
  double *beam_ms = malloc(num_beam_values*sizeof(double));

  user_precision_complex_t *primay_beam_J00 = malloc(num_beam_values*sizeof(user_precision_complex_t));
  user_precision_complex_t *primay_beam_J11 = malloc(num_beam_values*sizeof(user_precision_complex_t));

  double *freqs = malloc(num_freqs*sizeof(double));

  int vary_l = 1;
  get_1D_gaussian_values(beam_ls, beam_ms,
                         primay_beam_J00, primay_beam_J11,
                         freqs, beam_ref_freq,
                         vary_l, fwhm_lm,
                         num_freqs, num_times, num_components);

  int beam_ind = 0;
  for (int freq_ind = 0; freq_ind < num_freqs; freq_ind++) {
    for (int comp_ind = 0; comp_ind < num_components; comp_ind++) {

      double std = (fwhm_lm / FWHM_FACTOR) * (beam_ref_freq / freqs[freq_ind]);
      double exp_inside = beam_ls[comp_ind] / std;
      double estimate = exp(-0.5*exp_inside*exp_inside);

      TEST_ASSERT_DOUBLE_WITHIN(TOL, estimate, creal(primay_beam_J00[beam_ind]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, estimate, creal(primay_beam_J11[beam_ind]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(primay_beam_J00[beam_ind]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(primay_beam_J11[beam_ind]));

      beam_ind += 1;

    }
  }

  free(beam_ls);
  free(beam_ms);
  free(primay_beam_J00);
  free(primay_beam_J11);
  free(freqs);
}

void test_GaussBeam_GivesCorrectmValuesByFreq(void) {
  int num_freqs = 5;
  int num_times = 1;
  int num_components = 100;
  int num_beam_values = num_freqs*num_components*num_times;

  //Set full width half max to 20 deg, convert to l,m coord
  user_precision_t fwhm = 20.0;
  user_precision_t fwhm_lm = sin(fwhm*DD2R);
  double beam_ref_freq = 150e+6;

  double *beam_ls = malloc(num_beam_values*sizeof(double));
  double *beam_ms = malloc(num_beam_values*sizeof(double));

  user_precision_complex_t *primay_beam_J00 = malloc(num_beam_values*sizeof(user_precision_complex_t));
  user_precision_complex_t *primay_beam_J11 = malloc(num_beam_values*sizeof(user_precision_complex_t));

  double *freqs = malloc(num_freqs*sizeof(double));

  int vary_l = 0;
  get_1D_gaussian_values(beam_ls, beam_ms,
                         primay_beam_J00, primay_beam_J11,
                         freqs, beam_ref_freq,
                         vary_l, fwhm_lm,
                         num_freqs, num_times, num_components);

  int beam_ind = 0;
  for (int freq_ind = 0; freq_ind < num_freqs; freq_ind++) {
    for (int comp_ind = 0; comp_ind < num_components; comp_ind++) {

      double std = (fwhm_lm / FWHM_FACTOR) * (beam_ref_freq / freqs[freq_ind]);
      double exp_inside = beam_ms[comp_ind] / std;
      double estimate = expf(-0.5*exp_inside*exp_inside);

      TEST_ASSERT_DOUBLE_WITHIN(TOL, estimate, creal(primay_beam_J00[beam_ind]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, estimate, creal(primay_beam_J11[beam_ind]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(primay_beam_J00[beam_ind]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(primay_beam_J11[beam_ind]));

      beam_ind += 1;

    }
  }

  free(beam_ls);
  free(beam_ms);
  free(primay_beam_J00);
  free(primay_beam_J11);
  free(freqs);
}

int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_GaussBeam_GivesCorrectlValues);
    RUN_TEST(test_GaussBeam_GivesCorrectmValues);
    RUN_TEST(test_GaussBeam_GivesCorrectlValuesByFreq);
    RUN_TEST(test_GaussBeam_GivesCorrectmValuesByFreq);

    return UNITY_END();
}
