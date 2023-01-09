#include <math.h>
#include <unity.h>
#include <stdlib.h>
#include <complex.h>

#include "constants.h"
#include "woden_struct_defs.h"
#include "woden_precision_defs.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

//External CUDA code we're linking in
extern void test_kern_get_beam_gains(int num_freqs, int num_visis,
          int num_baselines, int num_components, int num_times, int beamtype,
          user_precision_complex_t *primay_beam_J00, user_precision_complex_t *primay_beam_J01,
          user_precision_complex_t *primay_beam_J10, user_precision_complex_t *primay_beam_J11,
          user_precision_complex_t *recover_g1x, user_precision_complex_t *recover_D1x,
          user_precision_complex_t *recover_D1y, user_precision_complex_t *recover_g1y,
          user_precision_complex_t *recover_g2x, user_precision_complex_t *recover_D2x,
          user_precision_complex_t *recover_D2y, user_precision_complex_t *recover_g2y);

#define UNITY_INCLUDE_FLOAT

/*
Test the __device__ code that grabs the beam gains works correctly
for all beam types
*/
void test_kern_get_beam_gains_ChooseBeams(int beamtype) {

  #ifdef DOUBLE_PRECISION
  printf("WODEN is using DOUBLE precision\n");
  #else
  printf("WODEN is using FLOAT precision\n");
  #endif

  int num_baselines = 3.0;
  int num_times = 2.0;
  int num_freqs = 2.0;

  int num_visis = num_baselines*num_times*num_freqs;

  int num_components = 4.0;

  user_precision_complex_t *primay_beam_J00 = malloc(num_freqs*num_times*num_components*sizeof(user_precision_complex_t));
  user_precision_complex_t *primay_beam_J01 = malloc(num_freqs*num_times*num_components*sizeof(user_precision_complex_t));
  user_precision_complex_t *primay_beam_J10 = malloc(num_freqs*num_times*num_components*sizeof(user_precision_complex_t));
  user_precision_complex_t *primay_beam_J11 = malloc(num_freqs*num_times*num_components*sizeof(user_precision_complex_t));

  //All models apart from NO_BEAM should have gains set other than 1.0
  int count = 0;
  if (beamtype != NO_BEAM) {
    for (int visi = 0; visi < num_components*num_times*num_freqs; visi++) {
      primay_beam_J00[visi] = count + I*0.0;
      primay_beam_J11[visi] = count + I*0.0;
      count ++ ;
    }
  }

  //Only MWA beams have cross-pols
  count = 0;
  if (beamtype == FEE_BEAM || beamtype == FEE_BEAM_INTERP || beamtype == MWA_ANALY) {
    for (int visi = 0; visi < num_components*num_times*num_freqs; visi++) {
      primay_beam_J01[visi] = count + I*0.0;
      primay_beam_J10[visi] = count + I*0.0;
      count ++ ;
    }
  }

  user_precision_complex_t *recover_g1x = malloc(num_visis*num_components*sizeof(user_precision_complex_t));
  user_precision_complex_t *recover_D1x = malloc(num_visis*num_components*sizeof(user_precision_complex_t));
  user_precision_complex_t *recover_D1y = malloc(num_visis*num_components*sizeof(user_precision_complex_t));
  user_precision_complex_t *recover_g1y = malloc(num_visis*num_components*sizeof(user_precision_complex_t));
  user_precision_complex_t *recover_g2x = malloc(num_visis*num_components*sizeof(user_precision_complex_t));
  user_precision_complex_t *recover_D2x = malloc(num_visis*num_components*sizeof(user_precision_complex_t));
  user_precision_complex_t *recover_D2y = malloc(num_visis*num_components*sizeof(user_precision_complex_t));
  user_precision_complex_t *recover_g2y = malloc(num_visis*num_components*sizeof(user_precision_complex_t));

  //Run the CUDA code
  test_kern_get_beam_gains(num_freqs, num_visis,
          num_baselines, num_components, num_times, beamtype,
          primay_beam_J00, primay_beam_J01,
          primay_beam_J10, primay_beam_J11,
          recover_g1x, recover_D1x,
          recover_D1y, recover_g1y,
          recover_g2x, recover_D2x,
          recover_D2y, recover_g2y);

  user_precision_t expected_order[] = { 0.00, 0.00, 0.00, 4.00, 4.00, 4.00, 8.00,
                             8.00, 8.00, 12.00, 12.00, 12.00, 1.00, 1.00,
                             1.00, 5.00, 5.00, 5.00, 9.00, 9.00, 9.00,
                             13.00, 13.00, 13.00, 2.00, 2.00, 2.00, 6.00,
                             6.00, 6.00, 10.00, 10.00, 10.00, 14.00, 14.00,
                             14.00, 3.00, 3.00, 3.00, 7.00, 7.00, 7.00,
                             11.00, 11.00, 11.00, 15.00, 15.00, 15.00 };

  //These models are real only, and have no leakage terms
  if (beamtype == ANALY_DIPOLE || beamtype == GAUSS_BEAM) {
    for (int output = 0; output < num_visis*num_components; output++) {

      TEST_ASSERT_EQUAL_FLOAT(expected_order[output], creal(recover_g1x[output]));
      TEST_ASSERT_EQUAL_FLOAT(0.0, creal(recover_D1x[output]));
      TEST_ASSERT_EQUAL_FLOAT(0.0, creal(recover_D1y[output]));
      TEST_ASSERT_EQUAL_FLOAT(expected_order[output], creal(recover_g1y[output]));
      TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(recover_g1x[output]));
      TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(recover_D1x[output]));
      TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(recover_D1y[output]));
      TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(recover_g1y[output]));

      TEST_ASSERT_EQUAL_FLOAT(expected_order[output], creal(recover_g2x[output]));
      TEST_ASSERT_EQUAL_FLOAT(0.0, creal(recover_D2x[output]));
      TEST_ASSERT_EQUAL_FLOAT(0.0, creal(recover_D2y[output]));
      TEST_ASSERT_EQUAL_FLOAT(expected_order[output], creal(recover_g2y[output]));
      TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(recover_g2x[output]));
      TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(recover_D2x[output]));
      TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(recover_D2y[output]));
      TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(recover_g2y[output]));

    }
  }
  else if (beamtype == FEE_BEAM || beamtype == FEE_BEAM_INTERP || beamtype == MWA_ANALY) {
    for (int output = 0; output < num_visis*num_components; output++) {

      TEST_ASSERT_EQUAL_FLOAT(expected_order[output], creal(recover_g1x[output]));
      TEST_ASSERT_EQUAL_FLOAT(expected_order[output], creal(recover_D1x[output]));
      TEST_ASSERT_EQUAL_FLOAT(expected_order[output], creal(recover_D1y[output]));
      TEST_ASSERT_EQUAL_FLOAT(expected_order[output], creal(recover_g1y[output]));
      TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(recover_g1x[output]));
      TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(recover_D1x[output]));
      TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(recover_D1y[output]));
      TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(recover_g1y[output]));

      TEST_ASSERT_EQUAL_FLOAT(expected_order[output], creal(recover_g2x[output]));
      TEST_ASSERT_EQUAL_FLOAT(expected_order[output], creal(recover_D2x[output]));
      TEST_ASSERT_EQUAL_FLOAT(expected_order[output], creal(recover_D2y[output]));
      TEST_ASSERT_EQUAL_FLOAT(expected_order[output], creal(recover_g2y[output]));
      TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(recover_g2x[output]));
      TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(recover_D2x[output]));
      TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(recover_D2y[output]));
      TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(recover_g2y[output]));

    }
  }
  //Should only be here if NO_BEAM, which just has gains of 1.0
  else {
    for (int output = 0; output < num_visis*num_components; output++) {
      TEST_ASSERT_EQUAL_FLOAT(1.0, creal(recover_g1x[output]));
      TEST_ASSERT_EQUAL_FLOAT(0.0, creal(recover_D1x[output]));
      TEST_ASSERT_EQUAL_FLOAT(0.0, creal(recover_D1y[output]));
      TEST_ASSERT_EQUAL_FLOAT(1.0, creal(recover_g1y[output]));
      TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(recover_g1x[output]));
      TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(recover_D1x[output]));
      TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(recover_D1y[output]));
      TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(recover_g1y[output]));

      TEST_ASSERT_EQUAL_FLOAT(1.0, creal(recover_g2x[output]));
      TEST_ASSERT_EQUAL_FLOAT(0.0, creal(recover_D2x[output]));
      TEST_ASSERT_EQUAL_FLOAT(0.0, creal(recover_D2y[output]));
      TEST_ASSERT_EQUAL_FLOAT(1.0, creal(recover_g2y[output]));
      TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(recover_g2x[output]));
      TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(recover_D2x[output]));
      TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(recover_D2y[output]));
      TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(recover_g2y[output]));
    }
  }

  //Be free my beauties
  free(primay_beam_J00);
  free(primay_beam_J01);
  free(primay_beam_J10);
  free(primay_beam_J11);
  free(recover_g1x);
  free(recover_D1x);
  free(recover_D1y);
  free(recover_g1y);
  free(recover_g2x);
  free(recover_D2x);
  free(recover_D2y);
  free(recover_g2y);

}


/*
This test checks the correct behaviour when simuating with beamtype=FEE_BEAM
*/
void test_kern_get_beam_gains_FEEBeam(void) {
  test_kern_get_beam_gains_ChooseBeams(FEE_BEAM);
}

/*
This test checks the correct behaviour when simuating with beamtype=ANALY_DIPOLE
*/
void test_kern_get_beam_gains_AnalyDipoleBeam(void) {
  test_kern_get_beam_gains_ChooseBeams(ANALY_DIPOLE);
}

/*
This test checks the correct behaviour when simuating with beamtype=GAUSS_BEAM
*/
void test_kern_get_beam_gains_GaussBeam(void) {
  test_kern_get_beam_gains_ChooseBeams(GAUSS_BEAM);
}

/*
This test checks the correct behaviour when simuating with beamtype=NO_BEAM
*/
void test_kern_get_beam_gains_NoBeam(void) {
  test_kern_get_beam_gains_ChooseBeams(NO_BEAM);
}

/*
This test checks the correct behaviour when simuating with beamtype=FEE_BEAM_INTERP
*/
void test_kern_get_beam_gains_FEEBeamInterp(void) {
  test_kern_get_beam_gains_ChooseBeams(FEE_BEAM_INTERP);
}

/*
This test checks the correct behaviour when simuating with beamtype=MWA_ANALY
*/
void test_kern_get_beam_gains_MWAAnaly(void) {
  test_kern_get_beam_gains_ChooseBeams(MWA_ANALY);
}


//Run the test with unity
int main(void)
{
    UNITY_BEGIN();
    RUN_TEST(test_kern_get_beam_gains_FEEBeam);
    RUN_TEST(test_kern_get_beam_gains_AnalyDipoleBeam);
    RUN_TEST(test_kern_get_beam_gains_GaussBeam);
    RUN_TEST(test_kern_get_beam_gains_NoBeam);
    RUN_TEST(test_kern_get_beam_gains_FEEBeamInterp);
    RUN_TEST(test_kern_get_beam_gains_MWAAnaly);


    return UNITY_END();
}
