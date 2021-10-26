#include <unity.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "constants.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

//External CUDA code we're linking in
extern void test_kern_gaussian_beam(double *beam_ls, double *beam_ms,
           float beam_ref_freq, float *freqs,
           float fwhm_lm, float cos_theta, float sin_theta, float sin_2theta,
           int num_freqs, int num_time_steps, int num_components,
           float _Complex *primay_beam_J00, float _Complex *primay_beam_J11);

#define UNITY_INCLUDE_FLOAT

/*
This function calls test_kern_gaussian_beam. It calculates a set of l,m input
coords, either varying l or m, keeping the other coord set to zero. This way
it tests a strip in l or m
*/
void get_1D_gaussian_values(double *beam_ls, double *beam_ms,
                            float _Complex *primay_beam_J00,
                            float _Complex *primay_beam_J11,
                            float *freqs, float beam_ref_freq,
                            int vary_l, float fwhm_lm,
                            int num_freqs, int num_times, int num_components){

  // float beam_ref_freq = 150e+6;

  int l_range = 2;
  double l_inc = l_range / num_components;

  int num_beam_values = num_freqs*num_components*num_times;

  //TODO currently hardcoded to have beam position angle = 0.
  //Should this change with az/za?
  float cos_theta = 1.0;
  float sin_theta = 0.0;
  float sin_2theta = 0.0;

  if (vary_l == 0) {
    for (size_t i = 0; i < num_components; i++) {
      beam_ls[i] = 0.0;
      beam_ms[i] = -1 + i*l_inc;
    }
  }
  else {
    for (size_t i = 0; i < num_components; i++) {
      beam_ls[i] = -1 + i*l_inc;
      beam_ms[i] = 0.0;
    }
  }

  float freq_inc = 25e+5;

  for (size_t i = 0; i < num_freqs; i++) {
    freqs[i] = beam_ref_freq + i*freq_inc;
  }

  test_kern_gaussian_beam(beam_ls, beam_ms,
                          beam_ref_freq, freqs,
                          fwhm_lm, cos_theta, sin_theta, sin_2theta,
                          num_freqs, num_times, num_components,
                          primay_beam_J00, primay_beam_J11);
}

void test_analytic_dipole_beam_GivesCorrectlValues(void) {
  int num_freqs = 1;
  int num_times = 1;
  int num_components = 100;
  int num_beam_values = num_freqs*num_components*num_times;

  //Set full width half max to 20 deg, convert to l,m coord
  float fwhm = 20.0;
  float fwhm_lm = sinf(fwhm*DD2R);
  float beam_ref_freq = 150e+6;

  double *beam_ls = malloc(num_beam_values*sizeof(double));
  double *beam_ms = malloc(num_beam_values*sizeof(double));

  float _Complex *primay_beam_J00 = malloc(num_beam_values*sizeof(float _Complex));
  float _Complex *primay_beam_J11 = malloc(num_beam_values*sizeof(float _Complex));

  float *freqs = malloc(num_freqs*sizeof(float));

  int vary_l = 1;
  get_1D_gaussian_values(beam_ls, beam_ms,
                         primay_beam_J00, primay_beam_J11,
                         freqs, beam_ref_freq,
                         vary_l, fwhm_lm,
                         num_freqs, num_times, num_components);

  for (size_t i = 0; i < num_components; i++) {

    float std = (fwhm_lm / FWHM_FACTOR);
    float exp_inside = (float)beam_ls[i] / std;
    float estimate = expf(-0.5*exp_inside*exp_inside);

    TEST_ASSERT_FLOAT_WITHIN(1e-7, estimate, creal(primay_beam_J00[i]));
    TEST_ASSERT_FLOAT_WITHIN(1e-7, estimate, creal(primay_beam_J11[i]));
    TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(primay_beam_J00[i]));
    TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(primay_beam_J11[i]));

  }
}

void test_analytic_dipole_beam_GivesCorrectmValues(void) {
  int num_freqs = 1;
  int num_times = 1;
  int num_components = 100;
  int num_beam_values = num_freqs*num_components*num_times;

  //Set full width half max to 20 deg, convert to l,m coord
  float fwhm = 20.0;
  float fwhm_lm = sinf(fwhm*DD2R);
  float beam_ref_freq = 150e+6;

  double *beam_ls = malloc(num_beam_values*sizeof(double));
  double *beam_ms = malloc(num_beam_values*sizeof(double));

  float _Complex *primay_beam_J00 = malloc(num_beam_values*sizeof(float _Complex));
  float _Complex *primay_beam_J11 = malloc(num_beam_values*sizeof(float _Complex));

  float *freqs = malloc(num_freqs*sizeof(float));

  int vary_l = 0;
  get_1D_gaussian_values(beam_ls, beam_ms,
                         primay_beam_J00, primay_beam_J11,
                         freqs, beam_ref_freq,
                         vary_l, fwhm_lm,
                         num_freqs, num_times, num_components);

  for (size_t i = 0; i < num_components; i++) {

    float std = (fwhm_lm / FWHM_FACTOR);
    float exp_inside = (float)beam_ms[i] / std;
    float estimate = expf(-0.5*exp_inside*exp_inside);

    TEST_ASSERT_FLOAT_WITHIN(1e-7, estimate, creal(primay_beam_J00[i]));
    TEST_ASSERT_FLOAT_WITHIN(1e-7, estimate, creal(primay_beam_J11[i]));
    TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(primay_beam_J00[i]));
    TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(primay_beam_J11[i]));

  }
}


void test_analytic_dipole_beam_GivesCorrectlValuesByFreq(void) {
  int num_freqs = 5;
  int num_times = 1;
  int num_components = 100;
  int num_beam_values = num_freqs*num_components*num_times;

  //Set full width half max to 20 deg, convert to l,m coord
  float fwhm = 20.0;
  float fwhm_lm = sinf(fwhm*DD2R);
  float beam_ref_freq = 150e+6;

  double *beam_ls = malloc(num_beam_values*sizeof(double));
  double *beam_ms = malloc(num_beam_values*sizeof(double));

  float _Complex *primay_beam_J00 = malloc(num_beam_values*sizeof(float _Complex));
  float _Complex *primay_beam_J11 = malloc(num_beam_values*sizeof(float _Complex));

  float *freqs = malloc(num_freqs*sizeof(float));

  int vary_l = 1;
  get_1D_gaussian_values(beam_ls, beam_ms,
                         primay_beam_J00, primay_beam_J11,
                         freqs, beam_ref_freq,
                         vary_l, fwhm_lm,
                         num_freqs, num_times, num_components);

  int beam_ind = 0;
  for (size_t freq_ind = 0; freq_ind < num_freqs; freq_ind++) {
    for (size_t comp_ind = 0; comp_ind < num_components; comp_ind++) {

      float std = (fwhm_lm / FWHM_FACTOR) * (beam_ref_freq / freqs[freq_ind]);
      float exp_inside = (float)beam_ls[comp_ind] / std;
      float estimate = expf(-0.5*exp_inside*exp_inside);

      TEST_ASSERT_FLOAT_WITHIN(1e-7, estimate, creal(primay_beam_J00[beam_ind]));
      TEST_ASSERT_FLOAT_WITHIN(1e-7, estimate, creal(primay_beam_J11[beam_ind]));
      TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(primay_beam_J00[beam_ind]));
      TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(primay_beam_J11[beam_ind]));

      beam_ind += 1;

    }
  }
}

void test_analytic_dipole_beam_GivesCorrectmValuesByFreq(void) {
  int num_freqs = 5;
  int num_times = 1;
  int num_components = 100;
  int num_beam_values = num_freqs*num_components*num_times;

  //Set full width half max to 20 deg, convert to l,m coord
  float fwhm = 20.0;
  float fwhm_lm = sinf(fwhm*DD2R);
  float beam_ref_freq = 150e+6;

  double *beam_ls = malloc(num_beam_values*sizeof(double));
  double *beam_ms = malloc(num_beam_values*sizeof(double));

  float _Complex *primay_beam_J00 = malloc(num_beam_values*sizeof(float _Complex));
  float _Complex *primay_beam_J11 = malloc(num_beam_values*sizeof(float _Complex));

  float *freqs = malloc(num_freqs*sizeof(float));

  int vary_l = 0;
  get_1D_gaussian_values(beam_ls, beam_ms,
                         primay_beam_J00, primay_beam_J11,
                         freqs, beam_ref_freq,
                         vary_l, fwhm_lm,
                         num_freqs, num_times, num_components);

  int beam_ind = 0;
  for (size_t freq_ind = 0; freq_ind < num_freqs; freq_ind++) {
    for (size_t comp_ind = 0; comp_ind < num_components; comp_ind++) {

      float std = (fwhm_lm / FWHM_FACTOR) * (beam_ref_freq / freqs[freq_ind]);
      float exp_inside = (float)beam_ms[comp_ind] / std;
      float estimate = expf(-0.5*exp_inside*exp_inside);

      TEST_ASSERT_FLOAT_WITHIN(1e-7, estimate, creal(primay_beam_J00[beam_ind]));
      TEST_ASSERT_FLOAT_WITHIN(1e-7, estimate, creal(primay_beam_J11[beam_ind]));
      TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(primay_beam_J00[beam_ind]));
      TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(primay_beam_J11[beam_ind]));

      beam_ind += 1;

    }
  }
}

int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_analytic_dipole_beam_GivesCorrectlValues);
    RUN_TEST(test_analytic_dipole_beam_GivesCorrectmValues);
    RUN_TEST(test_analytic_dipole_beam_GivesCorrectlValuesByFreq);
    // RUN_TEST(test_analytic_dipole_beam_GivesCorrectmValuesByFreq);

    return UNITY_END();
}
