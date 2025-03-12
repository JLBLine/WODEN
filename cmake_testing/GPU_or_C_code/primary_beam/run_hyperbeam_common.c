#include "run_hyperbeam_common.h"
// #include <mwa_hyperbeam.h>
#include "azza_radec_nside051.h"
#include "test_run_hyperbeam.h"

//Different delays settings, which control the pointing of the MWA beam
user_precision_t zenith_delays[16] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                           0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

user_precision_t off_zenith1_delays[16] = {0.0, 4.0, 8.0, 12.0, 0.0, 4.0, 8.0, 12.0,
                                0.0, 4.0, 8.0, 12.0, 0.0, 4.0, 8.0, 12.0};

user_precision_t off_zenith2_delays[16] = {0.0, 2.0, 4.0, 8.0, 2.0, 4.0, 8.0, 12.0,
                                4.0, 8.0, 12.0, 16.0, 8.0, 12.0, 16.0, 20.0};

void test_hyperbeam_VaryFreqVaryPointing(double freq, user_precision_t *delays,
                                         char* mwa_fee_hdf5, double *expected,
                                         char *outname, int rotate, int do_gpu) {

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
    // printf("There was an error calling new_fee_beam\n");
  }

  uint32_t *hyper_delays = malloc(16*sizeof(uint32_t));

  for (int delay = 0; delay < 16; delay++) {
    hyper_delays[delay] = (uint32_t)delays[delay];
  }

  double amps[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

  uint32_t num_tiles = 1;
  uint32_t num_amps = 16;
  uint8_t norm_to_zenith = 1;
  uint32_t num_freqs = 3;
  int num_beam_values = num_azza*num_freqs;

  //Check that it runs with all channel freqs, but returns the same thing
  //for all frequencies
  uint32_t freqs_hz[3] = {freq - 40e+3, freq, freq + 40e+3};
  double freqs[3] = {freq - 40e+3, freq, freq + 40e+3};

  // double freqs[3] = {freq - 1.28e+6, freq, freq + 1.28e+6};

  // uint32_t freqs_hz[1] = {freq};

  // if (do_gpu == 1){
    
  // }

  user_precision_complex_t *primay_beam_J00 = malloc(num_beam_values*sizeof(user_precision_complex_t));
  user_precision_complex_t *primay_beam_J01 = malloc(num_beam_values*sizeof(user_precision_complex_t));
  user_precision_complex_t *primay_beam_J10 = malloc(num_beam_values*sizeof(user_precision_complex_t));
  user_precision_complex_t *primay_beam_J11 = malloc(num_beam_values*sizeof(user_precision_complex_t));


  uint8_t parallactic = (uint8_t)rotate;
  int num_ants = 1;

  double latitudes[] = {-0.4660608448386394, -0.498};


  if (do_gpu == 1){

    struct FEEBeamGpu *gpu_fee_beam;
  
    status = new_gpu_fee_beam(fee_beam,
                               freqs_hz,
                               hyper_delays,
                               amps,
                               num_freqs,
                               num_tiles,
                               num_amps,
                               norm_to_zenith,
                               &gpu_fee_beam);
    if (status != 0) {
    handle_hyperbeam_error(__FILE__, __LINE__, "new_gpu_fee_beam");
    // printf("There was an error calling new_fee_beam\n");
    }

    test_run_hyperbeam_gpu(num_components,
               num_times, num_freqs, num_ants,
               parallactic,
               gpu_fee_beam,
               azs, zas,
               latitudes,
               primay_beam_J00,
               primay_beam_J01,
               primay_beam_J10,
               primay_beam_J11);
  
    free_gpu_fee_beam(gpu_fee_beam);
    
  } else {

    //When doing parallactic angle rotation, feed a different latitude for each
    //time step. az,za are stored normally by component, then time. Can do an
    //easy pointer arithmatic below to iterate over chunks of time if we reorder
    //them here (this is done internally to test_run_hyperbeam_gpu above)
    double *reordered_azs = malloc(num_azza*sizeof(double));
    double *reordered_zas = malloc(num_azza*sizeof(double));

    int stripe_new, stripe_old;

    for (int time_ind = 0; time_ind < num_times; time_ind++) {
      for (int comp_ind = 0; comp_ind < num_components; comp_ind++) {
        stripe_new = time_ind*num_components + comp_ind;
        stripe_old = comp_ind*num_times + time_ind;
        reordered_azs[stripe_new] = azs[stripe_old];
        reordered_zas[stripe_new] = zas[stripe_old];
      }
    }

    run_hyperbeam_cpu(num_components, num_times, num_freqs, num_ants,
                      parallactic, freqs, fee_beam,
                      hyper_delays, num_amps, amps,
                      reordered_azs, reordered_zas, latitudes,
                      primay_beam_J00, primay_beam_J01,
                      primay_beam_J10, primay_beam_J11);

    free(reordered_azs);
    free(reordered_zas);

    // test_run_hyperbeam_cpu(num_components, num_times, num_freqs, num_ants,
    //                   parallactic, freqs, fee_beam,
    //                   hyper_delays, num_amps, amps,
    //                   azs, zas, latitudes,
    //                   primay_beam_J00, primay_beam_J01,
    //                   primay_beam_J10, primay_beam_J11);
  }

  free_fee_beam(fee_beam);

  double TOL = 1e-6;

  for (int time = 0; time < num_times; time ++) {
    for (int freq = 0; freq < num_freqs; freq ++) {
      for (int comp = 0; comp < num_components; comp ++) {

        int beam_ind = num_freqs*time*num_components + num_components*freq + comp;

        int expected_base = 2*MAX_POLS*comp + 2*MAX_POLS*time*num_components;

        // printf(" %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n",
        //  creal(primay_beam_J00[beam_ind]), cimag(primay_beam_J00[beam_ind]),
        //  creal(primay_beam_J01[beam_ind]), cimag(primay_beam_J01[beam_ind]),
        //  creal(primay_beam_J10[beam_ind]), cimag(primay_beam_J10[beam_ind]),
        //  creal(primay_beam_J11[beam_ind]), cimag(primay_beam_J11[beam_ind]));


        // printf("want %.4f %.4f %.4f %.4f\n", expected[expected_base+0], expected[expected_base+1],
        //                                      expected[expected_base+6], expected[expected_base+7]);
        
        // printf(" got %.4f %.4f %.4f %.4f\n", creal(primay_beam_J00[beam_ind]),
        //                                      cimag(primay_beam_J00[beam_ind]),
        //                                      creal(primay_beam_J11[beam_ind]),
        //                                      cimag(primay_beam_J11[beam_ind]) );

        TEST_ASSERT_DOUBLE_WITHIN(TOL, expected[expected_base+0],
                                  creal(primay_beam_J00[beam_ind]) );
        TEST_ASSERT_DOUBLE_WITHIN(TOL, expected[expected_base+1],
                                  cimag(primay_beam_J00[beam_ind]) );
        TEST_ASSERT_DOUBLE_WITHIN(TOL, expected[expected_base+2],
                                  creal(primay_beam_J01[beam_ind]) );
        TEST_ASSERT_DOUBLE_WITHIN(TOL, expected[expected_base+3],
                                  cimag(primay_beam_J01[beam_ind]) );
        TEST_ASSERT_DOUBLE_WITHIN(TOL, expected[expected_base+4],
                                  creal(primay_beam_J10[beam_ind]) );
        TEST_ASSERT_DOUBLE_WITHIN(TOL, expected[expected_base+5],
                                  cimag(primay_beam_J10[beam_ind]) );
        TEST_ASSERT_DOUBLE_WITHIN(TOL, expected[expected_base+6],
                                  creal(primay_beam_J11[beam_ind]) );
        TEST_ASSERT_DOUBLE_WITHIN(TOL, expected[expected_base+7],
                                  cimag(primay_beam_J11[beam_ind]) );

      }
    }
  }

  if (do_gpu){
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


    for (int time = 0; time < num_times; time ++) {
      for (int freq = 0; freq < num_freqs; freq ++) {
        for (int comp = 0; comp < num_components; comp ++) {

          int beam_ind = num_freqs*time*num_components + num_components*freq + comp;
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
    fflush(beam_values_out);
    fclose(beam_values_out);
  }

  free(azs);
  free(zas);

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
                                char *outname, int do_gpu) {
  char* mwa_fee_hdf5 = getenv("MWA_FEE_HDF5");

  // #ifdef DOUBLE_PRECISION
  // printf("WODEN is using DOUBLE precision\n");
  // #else
  // printf("WODEN is using FLOAT precision\n");
  // #endif

  if (mwa_fee_hdf5) {
    printf("MWA_FEE_HDF5: %s\n", mwa_fee_hdf5 );

    int rotate;

    // Without rotation by parallactic angle
    rotate = 0;
    test_hyperbeam_VaryFreqVaryPointing(freq, delays, mwa_fee_hdf5,
                                         expected, outname, rotate, do_gpu);

    // //With rotation by parallactic angle
    rotate = 1;
    test_hyperbeam_VaryFreqVaryPointing(freq, delays, mwa_fee_hdf5,
                                         expected_rot, outname, rotate, do_gpu);
  }
  else {
    printf("MWA_FEE_HDF5 not found - not running test_RTS_FEE_beam test\n");
  }
}

/*
Run the test but vary the frequency and pointings. Compare to pre-calculated
values that are stored in test_RTS_FEE_beam.h
*/
void test_hyperbeam_100MHz_zenith(int do_gpu) {
  check_for_env_and_run_test(100e+6, zenith_delays, zenith_100,
                             zenith_100_rot, "hyperbeam_zenith_100", do_gpu);
}

void test_hyperbeam_150MHz_zenith(int do_gpu) {
  check_for_env_and_run_test(150e+6, zenith_delays, zenith_150,
                             zenith_150_rot, "hyperbeam_zenith_150", do_gpu);
}

void test_hyperbeam_200MHz_zenith(int do_gpu) {
  check_for_env_and_run_test(200e+6, zenith_delays, zenith_200,
                             zenith_200_rot, "hyperbeam_zenith_200", do_gpu);
}

void test_hyperbeam_100MHz_off_zenith1(int do_gpu) {
  check_for_env_and_run_test(100e+6, off_zenith1_delays, offzen1_100,
                             offzen1_100_rot, "hyperbeam_offzen1_100", do_gpu);
}

void test_hyperbeam_150MHz_off_zenith1(int do_gpu) {
  check_for_env_and_run_test(150e+6, off_zenith1_delays, offzen1_150,
                             offzen1_150_rot, "hyperbeam_offzen1_150", do_gpu);
}

void test_hyperbeam_200MHz_off_zenith1(int do_gpu) {
  check_for_env_and_run_test(200e+6, off_zenith1_delays, offzen1_200,
                             offzen1_200_rot, "hyperbeam_offzen1_200", do_gpu);
}

void test_hyperbeam_100MHz_off_zenith2(int do_gpu) {
  check_for_env_and_run_test(100e+6, off_zenith2_delays, offzen2_100,
                             offzen2_100_rot, "hyperbeam_offzen2_100", do_gpu);
}

void test_hyperbeam_150MHz_off_zenith2(int do_gpu) {
  check_for_env_and_run_test(150e+6, off_zenith2_delays, offzen2_150,
                             offzen2_150_rot, "hyperbeam_offzen2_150", do_gpu);
}

void test_hyperbeam_200MHz_off_zenith2(int do_gpu) {
  check_for_env_and_run_test(200e+6, off_zenith2_delays, offzen2_200,
                             offzen2_200_rot, "hyperbeam_offzen2_200", do_gpu);
}