#include <math.h>
#include <unity.h>
#include <stdlib.h>
#include <stdint.h>
#include <complex.h>

#include "constants.h"
#include "woden_struct_defs.h"
// #include "test_run_hyperbeam.h"
// #include "test_RTS_FEE_beam.h"
#include <mwa_hyperbeam.h>
// #include "azza_radec_nside101.h"
#include "test_run_hyperbeam_interp.h"
#include "test_interp_hyper_expected.h"

// #define ROTATION 0

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

// //External CUDA code we're linking in
extern void test_run_hyperbeam_cuda(int num_components,
           int num_times, int num_freqs,
           uint8_t parallatic,
           struct FEEBeamCUDA *cuda_fee_beam,
           double *azs, double *zas,
           user_precision_complex_t *primay_beam_J00,
           user_precision_complex_t *primay_beam_J01,
           user_precision_complex_t *primay_beam_J10,
           user_precision_complex_t *primay_beam_J11);

//Different delays settings, which control the pointing of the MWA beam
user_precision_t zenith_delays[16] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                           0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

user_precision_t off_zenith1_delays[16] = {0.0, 4.0, 8.0, 12.0, 0.0, 4.0, 8.0, 12.0,
                                0.0, 4.0, 8.0, 12.0, 0.0, 4.0, 8.0, 12.0};

user_precision_t off_zenith2_delays[16] = {0.0, 2.0, 4.0, 8.0, 2.0, 4.0, 8.0, 12.0,
                                4.0, 8.0, 12.0, 16.0, 8.0, 12.0, 16.0, 20.0};

void test_hyperbeam_interp(int freq_int,
                  int delay_int, user_precision_t *delays,
                  double base_low_freq, double freq_res, int num_freqs,
                  woden_settings_t *woden_settings,
                  double *expected_values) {

  // beam_settings_t *beam_settings = malloc(sizeof(beam_settings_t));

  woden_settings->base_low_freq = base_low_freq;
  woden_settings->num_freqs = num_freqs;

  struct FEEBeam *fee_beam;
  char error_str[100];

  int32_t status = 0;

  status =  new_fee_beam(woden_settings->hdf5_beam_path, &fee_beam, error_str);

  if (status != 0) {
    printf("hyperbeam error %d %s\n",status,error_str );
  }

  uint32_t *hyper_delays = malloc(16*sizeof(uint32_t));

  for (int delay = 0; delay < 16; delay++) {
    hyper_delays[delay] = (uint32_t)delays[delay];
  }

  double amps[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

  uint8_t num_times = 1;
  uint32_t num_tiles = 1;
  uint32_t num_amps = 16;
  uint8_t norm_to_zenith = 1;
  uint32_t num_freqs_hyper = woden_settings->num_freqs;

  int num_azza = NUM_COMPS*num_times;
  int num_beam_values = num_azza*num_freqs;

  //Convert some WODEN precision stuff into hyperbeam precision stuff
  uint32_t *freqs_hz = (uint32_t*)malloc(woden_settings->num_freqs*sizeof(uint32_t));

  for (int find = 0; find < woden_settings->num_freqs; find++) {
      // freqs_hz[freq_ind] = (uint32_t)visibility_set->channel_frequencies[freq_ind];
      freqs_hz[find] = (uint32_t)(base_low_freq + find*freq_res);

  }

  struct FEEBeamCUDA *cuda_fee_beam;

  status = new_cuda_fee_beam(fee_beam,
                          freqs_hz,
                          hyper_delays,
                          amps,
                          num_freqs_hyper,
                          num_tiles,
                          num_amps,
                          norm_to_zenith,
                          &cuda_fee_beam,
                          error_str);

  if (status != 0) {
    printf("hyperbeam error %d %s\n",status,error_str );
  }

  user_precision_complex_t *primay_beam_J00 = malloc(num_beam_values*sizeof(user_precision_complex_t));
  user_precision_complex_t *primay_beam_J01 = malloc(num_beam_values*sizeof(user_precision_complex_t));
  user_precision_complex_t *primay_beam_J10 = malloc(num_beam_values*sizeof(user_precision_complex_t));
  user_precision_complex_t *primay_beam_J11 = malloc(num_beam_values*sizeof(user_precision_complex_t));


  uint8_t parallatic = 1;

  // printf("ROTATES %d %d\n", parallatic, rotate );

  test_run_hyperbeam_cuda(NUM_COMPS,
             num_times, num_freqs,
             parallatic,
             cuda_fee_beam,
             azs, zas,
             primay_beam_J00,
             primay_beam_J01,
             primay_beam_J10,
             primay_beam_J11);

  free_cuda_fee_beam(cuda_fee_beam);
  free_fee_beam(fee_beam);
  // free(jones);
  #ifdef DOUBLE_PRECISION
    double TOL = 1e-10;
  #else
    double TOL = 1e-7;
  #endif

  //Because I'm using ctest to make some outputs for plotting, only do the ctest
  //for things I'm actually trying to test (I'm running some of the same
  //frequencies through the coarse beam model here to plot the difference)

  //   //Check answers are as expected
  int beam_ind, expec_ind;
  //

  if (freq_int != 4) {

    for (int freq = 0; freq < num_freqs; freq++) {
      for (int comp = 0; comp < NUM_COMPS; comp++) {

          // printf("%.1f %.1f %.1f %.1f %.1f %.1f %.1f %.1f\n",
          //          creal(primay_beam_J00[beam_ind]), cimag(primay_beam_J00[beam_ind]),
          //          creal(primay_beam_J01[beam_ind]), cimag(primay_beam_J01[beam_ind]),
          //          creal(primay_beam_J10[beam_ind]), cimag(primay_beam_J10[beam_ind]),
          //          creal(primay_beam_J11[beam_ind]), cimag(primay_beam_J11[beam_ind]));

          beam_ind = NUM_COMPS*freq + comp;
          //All pols and reals/imags are stored in one expected array so index
          //accordingly
          expec_ind = 2*MAX_POLS*beam_ind;
          TEST_ASSERT_DOUBLE_WITHIN(TOL, expected_values[expec_ind+0],
                                    creal(primay_beam_J00[beam_ind]));
          TEST_ASSERT_DOUBLE_WITHIN(TOL, expected_values[expec_ind+1],
                                    cimag(primay_beam_J00[beam_ind]));
          TEST_ASSERT_DOUBLE_WITHIN(TOL, expected_values[expec_ind+2],
                                    creal(primay_beam_J01[beam_ind]));
          TEST_ASSERT_DOUBLE_WITHIN(TOL, expected_values[expec_ind+3],
                                    cimag(primay_beam_J01[beam_ind]));
          TEST_ASSERT_DOUBLE_WITHIN(TOL, expected_values[expec_ind+4],
                                    creal(primay_beam_J10[beam_ind]));
          TEST_ASSERT_DOUBLE_WITHIN(TOL, expected_values[expec_ind+5],
                                    cimag(primay_beam_J10[beam_ind]));
          TEST_ASSERT_DOUBLE_WITHIN(TOL, expected_values[expec_ind+6],
                                    creal(primay_beam_J11[beam_ind]));
          TEST_ASSERT_DOUBLE_WITHIN(TOL, expected_values[expec_ind+7],
                                    cimag(primay_beam_J11[beam_ind]));
      }
    }

  }

  FILE *output_text;

  char filename[100];
  sprintf(filename, "hyperbeam_interp_delays%01d_freqs%01d.txt", delay_int, freq_int);

  output_text = fopen(filename,"w");

  int count = 0;


  for (int freq = 0; freq < num_freqs; freq ++) {
    for (int comp = 0; comp < NUM_COMPS; comp ++) {

      int beam_ind = NUM_COMPS*freq + comp;

      fprintf(output_text,"%.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.1d\n",
           creal(primay_beam_J00[beam_ind]), cimag(primay_beam_J00[beam_ind]),
           creal(primay_beam_J01[beam_ind]), cimag(primay_beam_J01[beam_ind]),
           creal(primay_beam_J10[beam_ind]), cimag(primay_beam_J10[beam_ind]),
           creal(primay_beam_J11[beam_ind]), cimag(primay_beam_J11[beam_ind]),
       freqs_hz[freq] );

       count ++;
    }
  }

  fflush(output_text);
  fclose(output_text);

  free(primay_beam_J00);
  free(primay_beam_J01);
  free(primay_beam_J10);
  free(primay_beam_J11);


}

/*
Check whether the environment variable for the FEE hdf5 beam exists, don't run
the test if it's missing
*/
void check_for_env_and_run_finetest(int freq_int, int delay_int,
                  user_precision_t *delays, double base_low_freq,
                  double freq_res, int num_freqs,
                  double *expected_values) {

  char* mwa_fee_hdf5 = getenv("MWA_FEE_HDF5_INTERP");
  woden_settings_t *woden_settings = malloc(sizeof(woden_settings_t));

  if (mwa_fee_hdf5) {
    printf("MWA_FEE_HDF5_INTERP: %s\n", mwa_fee_hdf5 );

    woden_settings->hdf5_beam_path = mwa_fee_hdf5;

    test_hyperbeam_interp(freq_int, delay_int,
                      delays, base_low_freq,
                      freq_res, num_freqs, woden_settings,
                      expected_values);
  }
  else {
    printf("MWA_FEE_HDF5_INTERP not found - not running test_RTS_FEE_beam test");
  }
}

/*
Check whether the environment variable for the FEE hdf5 beam exists, don't run
the test if it's missing
*/
void check_for_env_and_run_coarsetest(int freq_int, int delay_int,
                  user_precision_t *delays, double base_low_freq,
                  double freq_res, int num_freqs,
                  double *expected_values) {

  char* mwa_fee_hdf5 = getenv("MWA_FEE_HDF5");
  woden_settings_t *woden_settings = malloc(sizeof(woden_settings_t));

  if (mwa_fee_hdf5) {
    printf("MWA_FEE_HDF5: %s\n", mwa_fee_hdf5 );

    woden_settings->hdf5_beam_path = mwa_fee_hdf5;

    test_hyperbeam_interp(freq_int, delay_int,
                      delays, base_low_freq,
                      freq_res, num_freqs, woden_settings,
                      expected_values);
  }
  else {
    printf("MWA_FEE_HDF5 not found - not running test_RTS_FEE_beam test");
  }
}

void test_delays1_freqs1(){
  check_for_env_and_run_finetest(1, 1, delays1, FREQ1, FREQ_RES1, NUM_FREQS1,
                                 hyper_f1d1);
}

void test_delays2_freqs2(){
  check_for_env_and_run_finetest(2, 2, delays2, FREQ2, FREQ_RES2, NUM_FREQS2,
                                 hyper_f2d2);
}

void test_delays3_freqs3(){
  check_for_env_and_run_finetest(3, 3, delays3, FREQ3, FREQ_RES3, NUM_FREQS3,
                                 hyper_f3d3);
}

void test_delays4_freqs4(){
  check_for_env_and_run_coarsetest(4, 4, delays1, FREQ1, FREQ_RES1, NUM_FREQS1,
                                 hyper_f3d3);
}

//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_delays1_freqs1);
    RUN_TEST(test_delays2_freqs2);
    RUN_TEST(test_delays3_freqs3);
    RUN_TEST(test_delays4_freqs4);

    return UNITY_END();
}
