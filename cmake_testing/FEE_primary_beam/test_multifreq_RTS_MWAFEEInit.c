#include <unity.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "constants.h"
#include "FEE_primary_beam.h"
#include "woden_struct_defs.h"
#include "multifreq_common.h"
#include "test_multifreq_RTS_MWAFEEInit.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

/*
Test the initial setting up of the spherical harmonic values for the FEE beam
works for multiple frequencies when using the interpolated FEE beam values

Test the first and second entry of the Q1 and Q2 arrays
*/
void test_multifreq_RTS_MWAFEEInit(beam_settings_t *beam_settings,
                             woden_settings_t *woden_settings,
                             double *beam_freqs, double *expec_outputs,
                             double *expec_outputs_zen){


  multifreq_RTS_MWAFEEInit(beam_settings, woden_settings, beam_freqs);

  // double *expected_outputs = malloc(beam_settings->num_MWAFEE*8*sizeof(double));
  // double *expected_outputs_zenith = malloc(beam_settings->num_MWAFEE*8*sizeof(double));

  for (int freq_ind = 0; freq_ind < beam_settings->num_MWAFEE; freq_ind++) {

    // printf("%.1f\n",beam_freqs[freq_ind] );

    RTS_MWA_FEE_beam_t *FEE_beam = &beam_settings->FEE_beams[freq_ind];
    RTS_MWA_FEE_beam_t *FEE_beam_zenith = &beam_settings->FEE_beam_zeniths[freq_ind];

    for (int index = 0; index < 2; index++) {

      TEST_ASSERT_DOUBLE_WITHIN(1e-7, expec_outputs[freq_ind*8 + index*4 + 0],
                                      creal(FEE_beam->Q1[0][index]));
      TEST_ASSERT_DOUBLE_WITHIN(1e-7, expec_outputs[freq_ind*8 + index*4 + 1],
                                      creal(FEE_beam->Q2[0][index]));
      TEST_ASSERT_DOUBLE_WITHIN(1e-7, expec_outputs[freq_ind*8 + index*4 + 2],
                                      creal(FEE_beam->Q1[1][index]));
      TEST_ASSERT_DOUBLE_WITHIN(1e-7, expec_outputs[freq_ind*8 + index*4 + 3],
                                      creal(FEE_beam->Q2[1][index]));
      //
      TEST_ASSERT_DOUBLE_WITHIN(1e-7, expec_outputs_zen[freq_ind*8 + index*4 + 0],
                                      creal(FEE_beam_zenith->Q1[0][index]));
      TEST_ASSERT_DOUBLE_WITHIN(1e-7, expec_outputs_zen[freq_ind*8 + index*4 + 1],
                                      creal(FEE_beam_zenith->Q2[0][index]));
      TEST_ASSERT_DOUBLE_WITHIN(1e-7, expec_outputs_zen[freq_ind*8 + index*4 + 2],
                                      creal(FEE_beam_zenith->Q1[1][index]));
      TEST_ASSERT_DOUBLE_WITHIN(1e-7, expec_outputs_zen[freq_ind*8 + index*4 + 3],
                                      creal(FEE_beam_zenith->Q2[1][index]));

      // expected_outputs[freq_ind*8 + index*4 + 0] = creal(FEE_beam->Q1[0][index]);
      // expected_outputs[freq_ind*8 + index*4 + 1] = creal(FEE_beam->Q2[0][index]);
      // expected_outputs[freq_ind*8 + index*4 + 2] = creal(FEE_beam->Q1[1][index]);
      // expected_outputs[freq_ind*8 + index*4 + 3] = creal(FEE_beam->Q2[1][index]);
      //
      // expected_outputs_zenith[freq_ind*8 + index*4 + 0] = creal(FEE_beam_zenith->Q1[0][index]);
      // expected_outputs_zenith[freq_ind*8 + index*4 + 1] = creal(FEE_beam_zenith->Q2[0][index]);
      // expected_outputs_zenith[freq_ind*8 + index*4 + 2] = creal(FEE_beam_zenith->Q1[1][index]);
      // expected_outputs_zenith[freq_ind*8 + index*4 + 3] = creal(FEE_beam_zenith->Q2[1][index]);
    }
    // printf("--------------------------------------------------\n");

    RTS_freeHDFBeam(FEE_beam);
    RTS_freeHDFBeam(FEE_beam_zenith);

  }

  // for (int i = 0; i < beam_settings->num_MWAFEE*8; i++) {
  //   printf("%.8f, ",expected_outputs[i]);
  // }
  // printf("\n");
  // for (int i = 0; i < beam_settings->num_MWAFEE*8; i++) {
  //   printf("%.8f, ",expected_outputs_zenith[i]);
  // }
  // printf("\n");

  // free(expected_outputs);
  // free(expected_outputs_zenith);

}

void check_for_env_and_run_finetest(woden_settings_t *woden_settings,
                                beam_settings_t *beam_settings,
                                double *beam_freqs, double *expec_outputs,
                                double *expec_outputs_zen) {

  char* mwa_fee_hdf5 = getenv("MWA_FEE_HDF5_INTERP");

  if (mwa_fee_hdf5) {
    printf("MWA_FEE_HDF5_INTERP: %s\n", mwa_fee_hdf5 );

    woden_settings->hdf5_beam_path = mwa_fee_hdf5;

    test_multifreq_RTS_MWAFEEInit(beam_settings, woden_settings, beam_freqs,
                                  expec_outputs, expec_outputs_zen);

  }
  else {
    printf("MWA_FEE_HDF5 not found - not running test_RTS_FEE_beam test");
  }
}

void test_multifreq_RTS_MWAFEEInit_delays1_freqs1(){

  woden_settings_t *woden_settings = malloc(sizeof(woden_settings_t));
  beam_settings_t *beam_settings = malloc(sizeof(beam_settings_t));
  double *beam_freqs = malloc(NUM_FREQS1*sizeof(double));

  create_metafits_and_beam_freqs_delays1_freqs1(woden_settings, beam_freqs);
  check_for_env_and_run_finetest(woden_settings, beam_settings, beam_freqs,
                                 expec_freq1delay1_zen, expec_freq1delay1_zen);

}

void test_multifreq_RTS_MWAFEEInit_delays2_freqs2(){

  woden_settings_t *woden_settings = malloc(sizeof(woden_settings_t));
  beam_settings_t *beam_settings = malloc(sizeof(beam_settings_t));
  double *beam_freqs = malloc(NUM_FREQS2*sizeof(double));

  create_metafits_and_beam_freqs_delays2_freqs2(woden_settings, beam_freqs);
  check_for_env_and_run_finetest(woden_settings, beam_settings, beam_freqs,
                                 expec_freq2delay2, expec_freq2delay2_zen);

}

void test_multifreq_RTS_MWAFEEInit_delays3_freqs3(){

  woden_settings_t *woden_settings = malloc(sizeof(woden_settings_t));
  beam_settings_t *beam_settings = malloc(sizeof(beam_settings_t));
  double *beam_freqs = malloc(NUM_FREQS3*sizeof(double));

  create_metafits_and_beam_freqs_delays3_freqs3(woden_settings, beam_freqs);
  check_for_env_and_run_finetest(woden_settings, beam_settings, beam_freqs,
                                 expec_freq3delay3, expec_freq3delay3_zen);

}


//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_multifreq_RTS_MWAFEEInit_delays1_freqs1);
    RUN_TEST(test_multifreq_RTS_MWAFEEInit_delays2_freqs2);
    RUN_TEST(test_multifreq_RTS_MWAFEEInit_delays3_freqs3);

    return UNITY_END();
}
