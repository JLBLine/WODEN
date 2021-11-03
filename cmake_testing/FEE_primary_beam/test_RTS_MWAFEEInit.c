#include <unity.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "FEE_primary_beam.h"
#include "woden_struct_defs.h"
#include "test_RTS_MWAFEEInit.h"
#include <complex.h>

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

#define UNITY_INCLUDE_FLOAT

//Different delays settings, which control the pointing of the MWA beam
user_precision_t delays1[16] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

user_precision_t delays2[16] = {0.0, 4.0, 8.0, 12.0, 0.0, 4.0, 8.0, 12.0,
                     0.0, 4.0, 8.0, 12.0, 0.0, 4.0, 8.0, 12.0};

user_precision_t delays3[16] = {0.0, 2.0, 4.0, 8.0, 2.0, 4.0, 8.0, 12.0,
                     4.0, 8.0, 12.0, 16.0, 8.0, 12.0, 16.0, 20.0};


/*
FEE_primary_beam::RTS_MWAFEEInit grabs data out of an HDF5 file and selects
the correct subset based on the given delays and frequency. Check it works
for a number of delay settings and frequencies by checking a few values in the
6 output arrays
*/



/*
Run RTS_MWAFEEInit, and check that populated arrays match exected values
(which are defined in test_RTS_MWAFEEInit.h)
*/
void test_RTS_MWAFEEInit(user_precision_t freq, user_precision_t *delays,
                         char* mwa_fee_hdf5,
                         double *expec_Q1_real, double *expec_Q1_imag,
                         double *expec_Q2_real, double *expec_Q2_imag){

  //Call the C code to interrogate the hdf5 file and set beam things up
  RTS_MWA_FEE_beam_t *FEE_beam = malloc(sizeof(RTS_MWA_FEE_beam_t));

  RTS_MWAFEEInit(mwa_fee_hdf5, freq, FEE_beam, delays);

  // printf("nmax, nMN %d %d\n", FEE_beam->nmax, FEE_beam->nMN );

  //Loop through indexes and check results match expectations
  int test_index;
  int loop_index = 0;
  //There are two pols to test
  for (int pol = 0; pol < 2; pol++) {
    for (int index = 0; index < NUM_TEST_INDEXES; index++) {
      test_index = indexes_to_test[index];
      TEST_ASSERT_FLOAT_WITHIN(1e-7, (user_precision_t)expec_Q1_real[loop_index],creal(FEE_beam->Q1[pol][test_index]));
      TEST_ASSERT_FLOAT_WITHIN(1e-7, (user_precision_t)expec_Q1_imag[loop_index],cimag(FEE_beam->Q1[pol][test_index]));
      TEST_ASSERT_FLOAT_WITHIN(1e-7, (user_precision_t)expec_Q2_real[loop_index],creal(FEE_beam->Q2[pol][test_index]));
      TEST_ASSERT_FLOAT_WITHIN(1e-7, (user_precision_t)expec_Q2_imag[loop_index],cimag(FEE_beam->Q2[pol][test_index]));
      TEST_ASSERT_FLOAT_WITHIN(1e-7, (user_precision_t)expec_M[loop_index],FEE_beam->M[pol][test_index]);
      TEST_ASSERT_FLOAT_WITHIN(1e-7, (user_precision_t)expec_N[loop_index],FEE_beam->N[pol][test_index]);

      loop_index += 1;
    }
  }

  // printf("%.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f\n",
  //        creal(FEE_beam->Q1[0][indexes_to_test[0]]), creal(FEE_beam->Q1[0][indexes_to_test[1]]),
  //        creal(FEE_beam->Q1[0][indexes_to_test[2]]), creal(FEE_beam->Q1[0][indexes_to_test[3]]),
  //        creal(FEE_beam->Q1[0][indexes_to_test[4]]), creal(FEE_beam->Q1[0][indexes_to_test[5]]),
  //        creal(FEE_beam->Q1[1][indexes_to_test[0]]), creal(FEE_beam->Q1[1][indexes_to_test[1]]),
  //        creal(FEE_beam->Q1[1][indexes_to_test[2]]), creal(FEE_beam->Q1[1][indexes_to_test[3]]),
  //        creal(FEE_beam->Q1[1][indexes_to_test[4]]), creal(FEE_beam->Q1[1][indexes_to_test[5]]) );
  //
  // printf("%.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f\n",
  //        cimag(FEE_beam->Q1[0][indexes_to_test[0]]), cimag(FEE_beam->Q1[0][indexes_to_test[1]]),
  //        cimag(FEE_beam->Q1[0][indexes_to_test[2]]), cimag(FEE_beam->Q1[0][indexes_to_test[3]]),
  //        cimag(FEE_beam->Q1[0][indexes_to_test[4]]), cimag(FEE_beam->Q1[0][indexes_to_test[5]]),
  //        cimag(FEE_beam->Q1[1][indexes_to_test[0]]), cimag(FEE_beam->Q1[1][indexes_to_test[1]]),
  //        cimag(FEE_beam->Q1[1][indexes_to_test[2]]), cimag(FEE_beam->Q1[1][indexes_to_test[3]]),
  //        cimag(FEE_beam->Q1[1][indexes_to_test[4]]), cimag(FEE_beam->Q1[1][indexes_to_test[5]]) );
  //
  // printf("%.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f\n",
  //        creal(FEE_beam->Q2[0][indexes_to_test[0]]), creal(FEE_beam->Q2[0][indexes_to_test[1]]),
  //        creal(FEE_beam->Q2[0][indexes_to_test[2]]), creal(FEE_beam->Q2[0][indexes_to_test[3]]),
  //        creal(FEE_beam->Q2[0][indexes_to_test[4]]), creal(FEE_beam->Q2[0][indexes_to_test[5]]),
  //        creal(FEE_beam->Q2[1][indexes_to_test[0]]), creal(FEE_beam->Q2[1][indexes_to_test[1]]),
  //        creal(FEE_beam->Q2[1][indexes_to_test[2]]), creal(FEE_beam->Q2[1][indexes_to_test[3]]),
  //        creal(FEE_beam->Q2[1][indexes_to_test[4]]), creal(FEE_beam->Q2[1][indexes_to_test[5]]) );
  //
  // printf("%.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f\n",
  //        cimag(FEE_beam->Q2[0][indexes_to_test[0]]), cimag(FEE_beam->Q2[0][indexes_to_test[1]]),
  //        cimag(FEE_beam->Q2[0][indexes_to_test[2]]), cimag(FEE_beam->Q2[0][indexes_to_test[3]]),
  //        cimag(FEE_beam->Q2[0][indexes_to_test[4]]), cimag(FEE_beam->Q2[0][indexes_to_test[5]]),
  //        cimag(FEE_beam->Q2[1][indexes_to_test[0]]), cimag(FEE_beam->Q2[1][indexes_to_test[1]]),
  //        cimag(FEE_beam->Q2[1][indexes_to_test[2]]), cimag(FEE_beam->Q2[1][indexes_to_test[3]]),
  //        cimag(FEE_beam->Q2[1][indexes_to_test[4]]), cimag(FEE_beam->Q2[1][indexes_to_test[5]]) );
  //
  // printf("%.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f\n",
  //        FEE_beam->M[0][indexes_to_test[0]], FEE_beam->M[0][indexes_to_test[1]],
  //        FEE_beam->M[0][indexes_to_test[2]], FEE_beam->M[0][indexes_to_test[3]],
  //        FEE_beam->M[0][indexes_to_test[4]], FEE_beam->M[0][indexes_to_test[5]],
  //        FEE_beam->M[1][indexes_to_test[0]], FEE_beam->M[1][indexes_to_test[1]],
  //        FEE_beam->M[1][indexes_to_test[2]], FEE_beam->M[1][indexes_to_test[3]],
  //        FEE_beam->M[1][indexes_to_test[4]], FEE_beam->M[1][indexes_to_test[5]] );
  //
  // printf("%.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f\n",
  //        FEE_beam->N[0][indexes_to_test[0]], FEE_beam->N[0][indexes_to_test[1]],
  //        FEE_beam->N[0][indexes_to_test[2]], FEE_beam->N[0][indexes_to_test[3]],
  //        FEE_beam->N[0][indexes_to_test[4]], FEE_beam->N[0][indexes_to_test[5]],
  //        FEE_beam->N[1][indexes_to_test[0]], FEE_beam->N[1][indexes_to_test[1]],
  //        FEE_beam->N[1][indexes_to_test[2]], FEE_beam->N[1][indexes_to_test[3]],
  //        FEE_beam->N[1][indexes_to_test[4]], FEE_beam->N[1][indexes_to_test[5]] );

  RTS_freeHDFBeam(FEE_beam);

}

/*
Check whether the environment variable for the FEE hdf5 beam exists, don't run
the test if it's missing
*/
void check_for_env_and_run_test(user_precision_t freq, user_precision_t *delays,
                                double *expec_Q1_real, double *expec_Q1_imag,
                                double *expec_Q2_real, double *expec_Q2_imag) {

  char* mwa_fee_hdf5 = getenv("MWA_FEE_HDF5");

  if (mwa_fee_hdf5) {
    printf("MWA_FEE_HDF5: %s\n", mwa_fee_hdf5 );
    test_RTS_MWAFEEInit(freq, delays, mwa_fee_hdf5,
                        expec_Q1_real, expec_Q1_imag,
                        expec_Q2_real, expec_Q2_imag );

  }
  else {
    printf("MWA_FEE_HDF5 not found - not running test_RTS_FEE_beam test");
  }
}

/*
Check for first set of delays at 50MHz
*/
void test_horizon_test_Freq050MHzDelays1(void){

  check_for_env_and_run_test(50e+6, delays1,
                             expec_50MHz_d1_Q1_real, expec_50MHz_d1_Q1_imag,
                             expec_50MHz_d1_Q2_real, expec_50MHz_d1_Q2_imag);
}

/*
Check for first set of delays at 150MHz
*/
void test_horizon_test_Freq150MHzDelays1(void){

  check_for_env_and_run_test(150e+6, delays1,
                             expec_150MHz_d1_Q1_real, expec_150MHz_d1_Q1_imag,
                             expec_150MHz_d1_Q2_real, expec_150MHz_d1_Q2_imag);
}

/*
Check for first set of delays at 250MHz
*/
void test_horizon_test_Freq250MHzDelays1(void){

  check_for_env_and_run_test(250e+6, delays1,
                             expec_250MHz_d1_Q1_real, expec_250MHz_d1_Q1_imag,
                             expec_250MHz_d1_Q2_real, expec_250MHz_d1_Q2_imag);
}

/*
Check for second set of delays at 50MHz
*/
void test_horizon_test_Freq050MHzDelays2(void){

  check_for_env_and_run_test(50e+6, delays2,
                             expec_50MHz_d2_Q1_real, expec_50MHz_d2_Q1_imag,
                             expec_50MHz_d2_Q2_real, expec_50MHz_d2_Q2_imag);
}

/*
Check for second set of delays at 150MHz
*/
void test_horizon_test_Freq150MHzDelays2(void){

  check_for_env_and_run_test(150e+6, delays2,
                             expec_150MHz_d2_Q1_real, expec_150MHz_d2_Q1_imag,
                             expec_150MHz_d2_Q2_real, expec_150MHz_d2_Q2_imag);
}

/*
Check for second set of delays at 250MHz
*/
void test_horizon_test_Freq250MHzDelays2(void){

  check_for_env_and_run_test(250e+6, delays2,
                             expec_250MHz_d2_Q1_real, expec_250MHz_d2_Q1_imag,
                             expec_250MHz_d2_Q2_real, expec_250MHz_d2_Q2_imag);
}

/*
Check for third set of delays at 50MHz
*/
void test_horizon_test_Freq050MHzDelays3(void){

  check_for_env_and_run_test(50e+6, delays3,
                             expec_50MHz_d3_Q1_real, expec_50MHz_d3_Q1_imag,
                             expec_50MHz_d3_Q2_real, expec_50MHz_d3_Q2_imag);
}

/*
Check for third set of delays at 150MHz
*/
void test_horizon_test_Freq150MHzDelays3(void){

  check_for_env_and_run_test(150e+6, delays3,
                             expec_150MHz_d3_Q1_real, expec_150MHz_d3_Q1_imag,
                             expec_150MHz_d3_Q2_real, expec_150MHz_d3_Q2_imag);
}

/*
Check for third set of delays at 250MHz
*/
void test_horizon_test_Freq250MHzDelays3(void){

  check_for_env_and_run_test(250e+6, delays3,
                             expec_250MHz_d3_Q1_real, expec_250MHz_d3_Q1_imag,
                             expec_250MHz_d3_Q2_real, expec_250MHz_d3_Q2_imag);
}

//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_horizon_test_Freq050MHzDelays1);
    RUN_TEST(test_horizon_test_Freq150MHzDelays1);
    RUN_TEST(test_horizon_test_Freq250MHzDelays1);

    RUN_TEST(test_horizon_test_Freq050MHzDelays2);
    RUN_TEST(test_horizon_test_Freq150MHzDelays2);
    RUN_TEST(test_horizon_test_Freq250MHzDelays2);

    RUN_TEST(test_horizon_test_Freq050MHzDelays3);
    RUN_TEST(test_horizon_test_Freq150MHzDelays3);
    RUN_TEST(test_horizon_test_Freq250MHzDelays3);

    return UNITY_END();
}
