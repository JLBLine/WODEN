#include <math.h>
#include <unity.h>
#include <stdlib.h>
#include <complex.h>

#include "constants.h"
#include "woden_struct_defs.h"
#include "FEE_primary_beam.h"
#include "test_RTS_FEE_beam.h"
// #include "test_kern_calc_visi_common.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

//External CUDA code we're linking in
extern void test_RTS_CUDA_FEE_beam(int num_components,
           float *azs, float *zas, float latitude,
           RTS_MWA_FEE_beam_t *FEE_beam_zenith,
           RTS_MWA_FEE_beam_t *FEE_beam,
           int rotation, int scaling,
           float _Complex *FEE_beam_gains);

#define UNITY_INCLUDE_FLOAT

//Different delays settings, which control the pointing of the MWA beam
float zenith_delays[16] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                           0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

float off_zenith1_delays[16] = {0.0, 4.0, 8.0, 12.0, 0.0, 4.0, 8.0, 12.0,
                                0.0, 4.0, 8.0, 12.0, 0.0, 4.0, 8.0, 12.0};

float off_zenith2_delays[16] = {0.0, 2.0, 4.0, 8.0, 2.0, 4.0, 8.0, 12.0,
                                4.0, 8.0, 12.0, 16.0, 8.0, 12.0, 16.0, 20.0};

void test_RTS_CUDA_FEE_beam_VaryFreqVaryPointing(float freq, float *delays,
                                                 char* mwa_fee_hdf5,
                                                 float *expected,
                                                 char *outname) {

  //Call the C code to interrogate the hdf5 file and set beam things up
  RTS_MWA_FEE_beam_t *FEE_beam = malloc(sizeof(RTS_MWA_FEE_beam_t));

  //Get a zenith pointing beam for normalisation purposes
  RTS_MWA_FEE_beam_t *FEE_beam_zenith = malloc(sizeof(RTS_MWA_FEE_beam_t));

  printf("\n\tSetting up the zenith FEE beam...");
  RTS_MWAFEEInit(mwa_fee_hdf5, freq, FEE_beam_zenith, zenith_delays);
  printf(" done.\n");

  printf("\tSetting up the FEE beam...");
  RTS_MWAFEEInit(mwa_fee_hdf5, freq, FEE_beam, delays);
  printf(" done.\n");

  //Set up a bunch of az/za
  int num_za = 36;
  int num_az = 144;

  float za_res = (90.0 / num_za)*DD2R;
  float az_res = (360.0 / num_az)*DD2R;

  int num_azza = num_za*num_az;

  float *zas = malloc(num_azza*sizeof(float));
  float *azs = malloc(num_azza*sizeof(float));

  int count = 0;
  for (int az = 0; az < num_az; az++) {
    for (int za = 0; za < num_za; za++) {
      zas[count] = za*za_res;
      azs[count] = az*az_res;
      count ++;
    }
  }

  //Rotate by parallactic angles
  int rotation = 1;
  //Scale to zenith
  int scaling = 1;

  int num_pols = 4;
  float _Complex *FEE_beam_gains = malloc(num_pols*num_azza*sizeof(float _Complex));

  //Run the CUDA code
  test_RTS_CUDA_FEE_beam(num_azza,
             azs, zas, MWA_LAT,
             FEE_beam_zenith,
             FEE_beam,
             rotation, scaling,
             FEE_beam_gains);

  //Check the values are within tolerance
  float tol = 1e-6;
  for (size_t comp = 0; comp < num_azza; comp++) {
    TEST_ASSERT_FLOAT_WITHIN(tol, expected[2*num_pols*comp+0], creal(FEE_beam_gains[num_pols*comp+0]) );
    TEST_ASSERT_FLOAT_WITHIN(tol, expected[2*num_pols*comp+1], cimag(FEE_beam_gains[num_pols*comp+0]) );
    TEST_ASSERT_FLOAT_WITHIN(tol, expected[2*num_pols*comp+2], creal(FEE_beam_gains[num_pols*comp+1]) );
    TEST_ASSERT_FLOAT_WITHIN(tol, expected[2*num_pols*comp+3], cimag(FEE_beam_gains[num_pols*comp+1]) );
    TEST_ASSERT_FLOAT_WITHIN(tol, expected[2*num_pols*comp+4], creal(FEE_beam_gains[num_pols*comp+2]) );
    TEST_ASSERT_FLOAT_WITHIN(tol, expected[2*num_pols*comp+5], cimag(FEE_beam_gains[num_pols*comp+2]) );
    TEST_ASSERT_FLOAT_WITHIN(tol, expected[2*num_pols*comp+6], creal(FEE_beam_gains[num_pols*comp+3]) );
    TEST_ASSERT_FLOAT_WITHIN(tol, expected[2*num_pols*comp+7], cimag(FEE_beam_gains[num_pols*comp+3]) );
  }

  // FILE *beam_values_out;
  // char buff[0x100];
  // snprintf(buff, sizeof(buff), outname);
  // beam_values_out = fopen(buff,"w");
  //
  // for (size_t comp = 0; comp < num_azza; comp++) {
  //   fprintf(beam_values_out, "%.5f %.5f %.7f %.7f %.7f %.7f %.7f %.7f %.7f %.7f\n",
  //   azs[comp], zas[comp],
  //   creal(FEE_beam_gains[num_pols*comp+0]), cimag(FEE_beam_gains[num_pols*comp+0]),
  //   creal(FEE_beam_gains[num_pols*comp+1]), cimag(FEE_beam_gains[num_pols*comp+1]),
  //   creal(FEE_beam_gains[num_pols*comp+2]), cimag(FEE_beam_gains[num_pols*comp+2]),
  //   creal(FEE_beam_gains[num_pols*comp+3]), cimag(FEE_beam_gains[num_pols*comp+3]));
  // }

}

/*
Check whether the environment variable for the FEE hdf5 beam exists, don't run
the test if it's missing
*/
void check_for_env_and_run_test(float freq, float *delays, float *expected,
                                char *outname) {
  char* mwa_fee_hdf5 = getenv("MWA_FEE_HDF5");

  if (mwa_fee_hdf5) {
    printf("MWA_FEE_HDF5: %s", mwa_fee_hdf5 );
    test_RTS_CUDA_FEE_beam_VaryFreqVaryPointing(freq, delays, mwa_fee_hdf5, expected, outname);
  }
  else {
    printf("MWA_FEE_HDF5 not found - not running test_RTS_FEE_beam test");
  }
}

/*
Run the test but vary the frequency and pointings. Compare to pre-calculated
values that are stored in test_RTS_FEE_beam.h
*/
void test_RTS_CUDA_FEE_beam_100MHz_zenith(void) {
  check_for_env_and_run_test(100e+6, zenith_delays, zenith_100, "zenith_100.txt");
}

void test_RTS_CUDA_FEE_beam_150MHz_zenith(void) {
  check_for_env_and_run_test(150e+6, zenith_delays, zenith_150, "zenith_150.txt");
}

void test_RTS_CUDA_FEE_beam_200MHz_zenith(void) {
  check_for_env_and_run_test(200e+6, zenith_delays, zenith_200, "zenith_200.txt");
}

void test_RTS_CUDA_FEE_beam_100MHz_off_zenith1(void) {
  check_for_env_and_run_test(100e+6, off_zenith1_delays, offzen1_100, "offzen1_100.txt");
}

void test_RTS_CUDA_FEE_beam_150MHz_off_zenith1(void) {
  check_for_env_and_run_test(150e+6, off_zenith1_delays, offzen1_150, "offzen1_150.txt");
}

void test_RTS_CUDA_FEE_beam_200MHz_off_zenith1(void) {
  check_for_env_and_run_test(200e+6, off_zenith1_delays, offzen1_200, "offzen1_200.txt");
}

void test_RTS_CUDA_FEE_beam_100MHz_off_zenith2(void) {
  check_for_env_and_run_test(100e+6, off_zenith2_delays, offzen2_100, "offzen2_100.txt");
}

void test_RTS_CUDA_FEE_beam_150MHz_off_zenith2(void) {
  check_for_env_and_run_test(150e+6, off_zenith2_delays, offzen2_150, "offzen2_150.txt");
}

void test_RTS_CUDA_FEE_beam_200MHz_off_zenith2(void) {
  check_for_env_and_run_test(200e+6, off_zenith2_delays, offzen2_200, "offzen2_200.txt");
}


//Run the test with unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_RTS_CUDA_FEE_beam_100MHz_zenith);
    RUN_TEST(test_RTS_CUDA_FEE_beam_150MHz_zenith);
    RUN_TEST(test_RTS_CUDA_FEE_beam_200MHz_zenith);

    RUN_TEST(test_RTS_CUDA_FEE_beam_100MHz_off_zenith1);
    RUN_TEST(test_RTS_CUDA_FEE_beam_150MHz_off_zenith1);
    RUN_TEST(test_RTS_CUDA_FEE_beam_200MHz_off_zenith1);

    RUN_TEST(test_RTS_CUDA_FEE_beam_100MHz_off_zenith2);
    RUN_TEST(test_RTS_CUDA_FEE_beam_150MHz_off_zenith2);
    RUN_TEST(test_RTS_CUDA_FEE_beam_200MHz_off_zenith2);

    return UNITY_END();
}
