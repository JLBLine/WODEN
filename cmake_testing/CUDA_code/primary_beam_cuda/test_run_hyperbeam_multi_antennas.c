#include <math.h>
#include <unity.h>
#include <stdlib.h>
#include <stdint.h>
#include <complex.h>

#include "constants.h"
#include "woden_struct_defs.h"
#include "test_run_hyperbeam.h"
#include "hyperbeam_error.h"
// #include <mwa_hyperbeam.h>
#include "azza_radec_nside051.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

// //External CUDA code we're linking in
extern void test_run_hyperbeam_cuda(int num_components,
           int num_time_steps, int num_freqs, int num_ants,
           uint8_t parallatic, 
           struct FEEBeamGpu *cuda_fee_beam,
           double *azs, double *zas,
           double *latitudes, 
           user_precision_complex_t *primay_beam_J00,
           user_precision_complex_t *primay_beam_J01,
           user_precision_complex_t *primay_beam_J10,
           user_precision_complex_t *primay_beam_J11);

//Different delays settings, which control the pointing of the MWA beam
user_precision_t zenith_delays[16] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

// user_precision_t zenith_delays[32] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                      // 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                      // 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                      // 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,};

user_precision_t off_zenith1_delays[16] = {0.0, 4.0, 8.0, 12.0, 0.0, 4.0, 8.0, 12.0,
                                0.0, 4.0, 8.0, 12.0, 0.0, 4.0, 8.0, 12.0};

user_precision_t off_zenith2_delays[16] = {0.0, 2.0, 4.0, 8.0, 2.0, 4.0, 8.0, 12.0,
                                4.0, 8.0, 12.0, 16.0, 8.0, 12.0, 16.0, 20.0};

void test_hyperbeam_VaryFreqVaryPointing(double freq,
                                         user_precision_t *delays,
                                         char* mwa_fee_hdf5,
                                         double *expected,
                                         char *outname, int rotate) {

  int nside = 51;

  int num_times = 2;

  int num_components = nside*nside;
  int num_azza = num_components*num_times;

  double *azs = malloc(num_azza*sizeof(double));
  double *zas = malloc(num_azza*sizeof(double));

  for (int coord = 0; coord < num_components; coord++) {
    for (int time = 0; time < num_times; time++) {

      azs[coord*num_times + time] = nside051_azs[coord];
      zas[coord*num_times + time] = nside051_zas[coord];
    }
  }

  struct FEEBeam *fee_beam;
  // char error_str[100];

  int32_t status = 0;
  //
  // status =  new_fee_beam(mwa_fee_hdf5, &fee_beam, error_str);

  status = new_fee_beam(mwa_fee_hdf5, &fee_beam);
  if (status != 0) {
    handle_hyperbeam_error(__FILE__, __LINE__, "new_fee_beam");
  }

  uint32_t num_tiles = 3;

  //These are the amplitudes for the dipoles, as read in from metafits
  //I believe that they have X - east-west, Y - north-south
  double amps[96] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                     0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                     0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4,
                     0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6,
                     0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8,
                     1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};


  int num_delays_per_tile = 16;

  uint32_t *hyper_delays = malloc(num_tiles*num_delays_per_tile*sizeof(uint32_t));

  for (int delay = 0; delay < num_delays_per_tile; delay++) {
    for (int tile = 0; tile < num_tiles; tile++) {
      hyper_delays[tile*num_delays_per_tile + delay] = (uint32_t)delays[delay];
      hyper_delays[tile*num_delays_per_tile + delay] = (uint32_t)delays[delay];
    }
  }

  //MAKE a 2D array of amps, should be num_tiles * 32 (16 for X 16 for Y)
  
  //This num_amps is either 16 or 32, meaning either same amps for X,Y or
  //unique amps for X,Y
  uint32_t num_amps = 32;
  uint8_t norm_to_zenith = 1;
  uint32_t num_freqs = 3;
  int num_beam_values = num_azza*num_freqs*num_tiles;

  //Check that it runs with all channel freqs, but returns the same thing
  //for all frequencies
  uint32_t freqs_hz[3] = {freq - 40e+3, freq, freq + 40e+3};
  // uint32_t freqs_hz[1] = {freq};

  struct FEEBeamGpu *cuda_fee_beam;

  status = new_gpu_fee_beam(fee_beam,
                             freqs_hz,
                             hyper_delays,
                             amps,
                             num_freqs,
                             num_tiles,
                             num_amps,
                             norm_to_zenith,
                             &cuda_fee_beam);

  if (status != 0) {
    handle_hyperbeam_error(__FILE__, __LINE__, "new_gpu_fee_beam");
  }

  user_precision_complex_t *primay_beam_J00 = malloc(num_beam_values*sizeof(user_precision_complex_t));
  user_precision_complex_t *primay_beam_J01 = malloc(num_beam_values*sizeof(user_precision_complex_t));
  user_precision_complex_t *primay_beam_J10 = malloc(num_beam_values*sizeof(user_precision_complex_t));
  user_precision_complex_t *primay_beam_J11 = malloc(num_beam_values*sizeof(user_precision_complex_t));


  uint8_t parallatic = (uint8_t)rotate;

  double latitudes[] = {-0.4660608448386394, -0.498};
  test_run_hyperbeam_cuda(num_components,
             num_times, num_freqs, (int)num_tiles,
             parallatic,
             cuda_fee_beam,
             azs, zas,
             latitudes,
             primay_beam_J00,
             primay_beam_J01,
             primay_beam_J10,
             primay_beam_J11);

  free_gpu_fee_beam(cuda_fee_beam);
  free_fee_beam(fee_beam);

  double TOL = 1e-6;

  // int num_big_diffs = 0;

  // // // Check the values are within TOLerance
  // for (int comp = 0; comp < num_azza; comp++) {
  //
  //   // printf("%d %d %.16f %.16f\n",comp, num_azza, expected[2*MAX_POLS*comp+0], creal(primay_beam_J00[comp]) );
  //
  //   // printf("%.2e\n",fabs(expected[2*MAX_POLS*comp+0] - creal(primay_beam_J00[comp])) );
  //
  //   // if (fabs(expected[2*MAX_POLS*comp+0] - creal(primay_beam_J00[comp])) > TOL) {
  //   //   num_big_diffs += 1;
  //   //   printf("%d %.5f %.12f %.16f %.16f\n",comp, azs[comp], zas[comp], expected[2*MAX_POLS*comp+0], creal(primay_beam_J00[comp]) );
  //   // }
  //

  //Given the dip amps we set earlier, we can multiply the expected values by
  //one of these constants as appropriate

  //OKOK so when hyperdrive reads in amps, X = east-west, Y = north-south
  //I use it with iau_order = 1;, which switches Y = east-west, X = north-south
  //This means the expected values should be as below
  double antx_mult[3] = {0.2, 0.6, 1.0};
  double anty_mult[3] = {0.0, 0.4, 0.8};

  for (int ant = 0; ant < num_tiles; ant ++) {
    for (int time = 0; time < num_times; time ++) {
      for (int freq = 0; freq < num_freqs; freq ++) {
        for (int comp = 0; comp < num_components; comp ++) {

          int beam_ind = ant*num_freqs*num_times*num_components + num_freqs*time*num_components + num_components*freq + comp;

          int expected_base = 2*MAX_POLS*comp + 2*MAX_POLS*time*num_components;

          TEST_ASSERT_DOUBLE_WITHIN(TOL, antx_mult[ant]*expected[expected_base+0],
                                    creal(primay_beam_J00[beam_ind]) );
          TEST_ASSERT_DOUBLE_WITHIN(TOL, antx_mult[ant]*expected[expected_base+1],
                                    cimag(primay_beam_J00[beam_ind]) );
          TEST_ASSERT_DOUBLE_WITHIN(TOL, antx_mult[ant]*expected[expected_base+2],
                                    creal(primay_beam_J01[beam_ind]) );
          TEST_ASSERT_DOUBLE_WITHIN(TOL, antx_mult[ant]*expected[expected_base+3],
                                    cimag(primay_beam_J01[beam_ind]) );
          TEST_ASSERT_DOUBLE_WITHIN(TOL, anty_mult[ant]*expected[expected_base+4],
                                    creal(primay_beam_J10[beam_ind]) );
          TEST_ASSERT_DOUBLE_WITHIN(TOL, anty_mult[ant]*expected[expected_base+5],
                                    cimag(primay_beam_J10[beam_ind]) );
          TEST_ASSERT_DOUBLE_WITHIN(TOL, anty_mult[ant]*expected[expected_base+6],
                                    creal(primay_beam_J11[beam_ind]) );
          TEST_ASSERT_DOUBLE_WITHIN(TOL, anty_mult[ant]*expected[expected_base+7],
                                    cimag(primay_beam_J11[beam_ind]) );

        }
      }
    }
  }

  FILE *beam_values_out;
  char buff[0x100];

  #ifdef DOUBLE_PRECISION
  if (rotate == 1) {
    snprintf(buff, sizeof(buff), "%s_rot_double.txt", outname);
  } else {
    snprintf(buff, sizeof(buff), "%s_double.txt", outname);
  }

  #else
      snprintf(buff, sizeof(buff), "%s_float.txt", outname);
  #endif
  //
  beam_values_out = fopen(buff,"w");

  // for (int comp = 0; comp < num_azza; comp++) {
  //   fprintf(beam_values_out, "%.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.6e\n",
  //   azs[comp], zas[comp],
  //   creal(primay_beam_J00[comp]), cimag(primay_beam_J00[comp]),
  //   creal(primay_beam_J01[comp]), cimag(primay_beam_J01[comp]),
  //   creal(primay_beam_J10[comp]), cimag(primay_beam_J10[comp]),
  //   creal(primay_beam_J11[comp]), cimag(primay_beam_J11[comp]),
  //   freq);
  // }

  for (int ant = 0; ant < num_tiles; ant ++) {
    for (int time = 0; time < num_times; time ++) {
      for (int freq = 0; freq < num_freqs; freq ++) {
        for (int comp = 0; comp < num_components; comp ++) {
  
          int beam_ind = ant*num_freqs*num_times*num_components + num_freqs*time*num_components + num_components*freq + comp;
          // int beam_ind = num_freqs*time*num_components + num_components*freq + comp;
          int coord_ind = comp*num_times + time;
  
          fprintf(beam_values_out,"%.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.1d\n",
          azs[coord_ind], zas[coord_ind],
          creal(primay_beam_J00[beam_ind]), cimag(primay_beam_J00[beam_ind]),
          creal(primay_beam_J01[beam_ind]), cimag(primay_beam_J01[beam_ind]),
          creal(primay_beam_J10[beam_ind]), cimag(primay_beam_J10[beam_ind]),
          creal(primay_beam_J11[beam_ind]), cimag(primay_beam_J11[beam_ind]),
          freqs_hz[freq] );
  
        }
      }
    }
  }

  free(azs);
  free(zas);

  fflush(beam_values_out);
  fclose(beam_values_out);

  free(primay_beam_J00);
  free(primay_beam_J01);
  free(primay_beam_J10);
  free(primay_beam_J11);


}

/*
Check whether the environment variable for the FEE hdf5 beam exists, don't run
the test if it's missing
*/
void check_for_env_and_run_test(double freq, user_precision_t *delays,
                                double *expected, double *expected_rot,
                                char *outname) {
  char* mwa_fee_hdf5 = getenv("MWA_FEE_HDF5");

  #ifdef DOUBLE_PRECISION
  printf("WODEN is using DOUBLE precision\n");
  #else
  printf("WODEN is using FLOAT precision\n");
  #endif

  if (mwa_fee_hdf5) {
    printf("MWA_FEE_HDF5: %s\n", mwa_fee_hdf5 );

    int rotate;

    // // Without rotation by parallactic angle
    // rotate = 0;
    // test_hyperbeam_VaryFreqVaryPointing(freq, delays, mwa_fee_hdf5, expected, outname, rotate);

    // //With rotation by parallactic angle
    rotate = 1;
    test_hyperbeam_VaryFreqVaryPointing(freq, delays, mwa_fee_hdf5, expected_rot, outname, rotate);
  }
  else {
    printf("MWA_FEE_HDF5 not found - not running test_RTS_FEE_beam test\n");
  }
}

/*
Run the test but vary the frequency and pointings. Compare to pre-calculated
values that are stored in test_RTS_FEE_beam.h
*/
void test_hyperbeam_100MHz_zenith(void) {
  check_for_env_and_run_test(100e+6, zenith_delays, zenith_100,
                             zenith_100_rot, "hyperbeam_zenith_100");
}

void test_hyperbeam_150MHz_zenith(void) {
  check_for_env_and_run_test(150e+6, zenith_delays, zenith_150,
                             zenith_150_rot, "hyperbeam_zenith_150");
}

void test_hyperbeam_200MHz_zenith(void) {
  check_for_env_and_run_test(200e+6, zenith_delays, zenith_200,
                             zenith_200_rot, "hyperbeam_zenith_200_two_ants");
}

void test_hyperbeam_100MHz_off_zenith1(void) {
  check_for_env_and_run_test(100e+6, off_zenith1_delays, offzen1_100,
                             offzen1_100_rot, "hyperbeam_offzen1_100");
}

void test_hyperbeam_150MHz_off_zenith1(void) {
  check_for_env_and_run_test(150e+6, off_zenith1_delays, offzen1_150,
                             offzen1_150_rot, "hyperbeam_offzen1_150");
}

void test_hyperbeam_200MHz_off_zenith1(void) {
  check_for_env_and_run_test(200e+6, off_zenith1_delays, offzen1_200,
                             offzen1_200_rot, "hyperbeam_offzen1_200");
}

void test_hyperbeam_100MHz_off_zenith2(void) {
  check_for_env_and_run_test(100e+6, off_zenith2_delays, offzen2_100,
                             offzen2_100_rot, "hyperbeam_offzen2_100");
}

void test_hyperbeam_150MHz_off_zenith2(void) {
  check_for_env_and_run_test(150e+6, off_zenith2_delays, offzen2_150,
                             offzen2_150_rot, "hyperbeam_offzen2_150");
}

void test_hyperbeam_200MHz_off_zenith2(void) {
  check_for_env_and_run_test(200e+6, off_zenith2_delays, offzen2_200,
                             offzen2_200_rot, "hyperbeam_offzen2_200");
}


//Run the test with unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_hyperbeam_100MHz_zenith);
    RUN_TEST(test_hyperbeam_150MHz_zenith);
    RUN_TEST(test_hyperbeam_200MHz_zenith);
    
    RUN_TEST(test_hyperbeam_100MHz_off_zenith1);
    RUN_TEST(test_hyperbeam_150MHz_off_zenith1);
    RUN_TEST(test_hyperbeam_200MHz_off_zenith1);
    
    RUN_TEST(test_hyperbeam_100MHz_off_zenith2);
    RUN_TEST(test_hyperbeam_150MHz_off_zenith2);
    RUN_TEST(test_hyperbeam_200MHz_off_zenith2);

    return UNITY_END();
}
