#include <math.h>
#include <unity.h>
#include <stdlib.h>
#include <complex.h>

#include "constants.h"
#include "woden_struct_defs.h"
#include "woden_precision_defs.h"

void setUp (void) {} /* Is run before eVary test, put unit init calls here. */
void tearDown (void) {} /* Is run after eVary test, put unit clean-up calls here. */

//External CUDA code we're linking in
extern void test_kern_update_sum_visis(int num_freqs, int num_cross,
          int num_baselines, int num_components, int num_times, int beamtype,
          int use_twoants, int num_ants,
          user_precision_complex_t *primay_beam_J00,
          user_precision_complex_t *primay_beam_J01,
          user_precision_complex_t *primay_beam_J10,
          user_precision_complex_t *primay_beam_J11,
          user_precision_complex_t *visi_components,
          user_precision_t *flux_I, user_precision_t *flux_Q,
          user_precision_t *flux_U, user_precision_t *flux_V,
          user_precision_t *sum_visi_XX_real, user_precision_t *sum_visi_XX_imag,
          user_precision_t *sum_visi_XY_real, user_precision_t *sum_visi_XY_imag,
          user_precision_t *sum_visi_YX_real, user_precision_t *sum_visi_YX_imag,
          user_precision_t *sum_visi_YY_real, user_precision_t *sum_visi_YY_imag);

/*
Test the __device__ code that updates the summed visibilities by grabbing the
correct beam gain and mesurement equation, multiplying and summing onto the visi
Here we keep the component visibilities and fluxes constant and vary the beam gains
Test works for all primary beam types
*/
void test_kern_update_sum_visis_VaryGainChooseBeams(int beamtype) {


  int num_ants = 3;
  int use_twoants = 1;
  int num_baselines = ((num_ants - 1)*num_ants) / 2;
  int num_times = 2;
  int num_freqs = 2;

  int num_visis = num_baselines*num_times*num_freqs;

  int num_components = 3;

  int num_input_gains = num_freqs*num_times*num_components*num_ants;

  user_precision_complex_t *primay_beam_J00 = malloc(num_input_gains*sizeof(user_precision_complex_t));
  user_precision_complex_t *primay_beam_J01 = malloc(num_input_gains*sizeof(user_precision_complex_t));
  user_precision_complex_t *primay_beam_J10 = malloc(num_input_gains*sizeof(user_precision_complex_t));
  user_precision_complex_t *primay_beam_J11 = malloc(num_input_gains*sizeof(user_precision_complex_t));
  user_precision_complex_t *visi_components = malloc(num_visis*sizeof(user_precision_complex_t));

  //Just stick base measurement equation to 1.0
  for (int visi = 0; visi < num_visis; visi++) {
    visi_components[visi] = 1.0 + I*0.0;
  }

  //Fill in some varying beam gains based on beam type
  int count = 0;
  if (beamtype != NO_BEAM) {
    for (int visi = 0; visi < num_input_gains; visi++) {
      primay_beam_J00[visi] = count + I*0.0;
      primay_beam_J11[visi] = count + I*0.0;
      count ++ ;
    }
  }

  //Only FEE_BEAM has cross-pols
  count = 0;
  if (beamtype == FEE_BEAM || beamtype == FEE_BEAM_INTERP || beamtype == MWA_ANALY) {
    for (int visi = 0; visi < num_input_gains; visi++) {
      primay_beam_J01[visi] = count + I*0.0;
      primay_beam_J10[visi] = count + I*0.0;
      // primay_beam_J01[visi] = 0.0 + I*0.0;
      // primay_beam_J10[visi] = 0.0 + I*0.0;
      count ++ ;
    }
  }

  //Set all the Stokes I to 1.0 and other stokes to zero
  user_precision_t *flux_I = malloc(num_freqs*num_components*sizeof(user_precision_t));
  user_precision_t *flux_Q = calloc(num_freqs*num_components, sizeof(user_precision_t));
  user_precision_t *flux_U = calloc(num_freqs*num_components, sizeof(user_precision_t));
  user_precision_t *flux_V = calloc(num_freqs*num_components, sizeof(user_precision_t));

  //Just stick base measurement equation to 1.0
  for (int visi = 0; visi < num_freqs*num_components; visi++) {
    flux_I[visi] = 1.0;
  }

  //Make sure arrays to hold summed visis are initialised to zero
  user_precision_t *sum_visi_XX_real = calloc(num_visis, sizeof(user_precision_t));
  user_precision_t *sum_visi_XX_imag = calloc(num_visis, sizeof(user_precision_t));
  user_precision_t *sum_visi_XY_real = calloc(num_visis, sizeof(user_precision_t));
  user_precision_t *sum_visi_XY_imag = calloc(num_visis, sizeof(user_precision_t));
  user_precision_t *sum_visi_YX_real = calloc(num_visis, sizeof(user_precision_t));
  user_precision_t *sum_visi_YX_imag = calloc(num_visis, sizeof(user_precision_t));
  user_precision_t *sum_visi_YY_real = calloc(num_visis, sizeof(user_precision_t));
  user_precision_t *sum_visi_YY_imag = calloc(num_visis, sizeof(user_precision_t));

  

  //Run the CUDA code
  test_kern_update_sum_visis(num_freqs, num_visis,
          num_baselines, num_components, num_times, beamtype,
          use_twoants, num_ants,
          primay_beam_J00, primay_beam_J01,
          primay_beam_J10, primay_beam_J11,
          visi_components,
          flux_I, flux_Q, flux_U, flux_V,
          sum_visi_XX_real, sum_visi_XX_imag,
          sum_visi_XY_real, sum_visi_XY_imag,
          sum_visi_YX_real, sum_visi_YX_imag,
          sum_visi_YY_real, sum_visi_YY_imag);

  //Do a bunch of interesting maths to predict what the outcomes should be
  int *ant1_indexes = malloc(num_ants*sizeof(int));
  int *ant2_indexes = malloc(num_ants*sizeof(int));

  int cross_index = 0;
  for (int ant1 = 0; ant1 < num_ants-1; ant1++)
  {
    for (int ant2 = ant1 + 1; ant2 < num_ants; ant2++)
    {
      ant1_indexes[cross_index] = ant1;
      ant2_indexes[cross_index] = ant2;

      cross_index += 1;
    }
  }

  int iBaseline = 0; 
  int time_ind = 0;
  int freq_ind = 0;
  int baseline_ind = 0;

  double expected_gain1;
  double expected_gain2;
  double expected_visi;

  double *expected_visis = malloc(num_visis*sizeof(double));

  int num_coords = num_components*num_times*num_freqs;

  //This does some reverse engineering to figure out for a give iBaseline (which
  //is a combined index over all times, freqs, cross-correlations) what the
  //gains that would go into that visibility would be
  for (int visi = 0; visi < num_visis; visi++) {
    // comp_ind = (int)floorf(visi / num_visis);
    iBaseline = visi;

    time_ind = (int)floorf( iBaseline / (num_baselines * num_freqs));
    freq_ind = (int)floorf( (iBaseline - (time_ind*num_baselines * num_freqs)) / num_baselines);
    baseline_ind = iBaseline - (time_ind*num_baselines * num_freqs) - (freq_ind*num_baselines);

    //OK so the visibilities are a sum over all components, so loop over
    //number of components add both gains

    expected_visi = 0;

    for (int comp_ind = 0; comp_ind < num_components; comp_ind++)
    {
      expected_gain1 = ant1_indexes[baseline_ind]*num_coords + time_ind*num_components*num_freqs + freq_ind*num_components + comp_ind; 
      expected_gain2 = ant2_indexes[baseline_ind]*num_coords + time_ind*num_components*num_freqs + freq_ind*num_components + comp_ind;

      expected_visi += expected_gain1*expected_gain2;
    }
    
    expected_visis[visi] = expected_visi;

  }

  if (beamtype == FEE_BEAM || beamtype == FEE_BEAM_INTERP || beamtype == MWA_ANALY) {
    for (int output = 0; output < num_visis; output++) {

      // printf("%.1f %.1f %.1f %.1f %.1f %.1f %.1f %.1f\n",
      //     sum_visi_XX_real[output], sum_visi_XX_imag[output],
      //     sum_visi_XY_real[output], sum_visi_XY_imag[output],
      //     sum_visi_YX_real[output], sum_visi_YX_imag[output],
      //     sum_visi_YY_real[output], sum_visi_YY_imag[output]);

      //We do a times two here because e.g. XX is g1x*g2x + D1y*D2y and
      //we've set the gains and leakages to be the same

      TEST_ASSERT_EQUAL_DOUBLE(2*expected_visis[output], sum_visi_XX_real[output]);
      TEST_ASSERT_EQUAL_DOUBLE(2*expected_visis[output], sum_visi_XY_real[output]);
      TEST_ASSERT_EQUAL_DOUBLE(2*expected_visis[output], sum_visi_YX_real[output]);
      TEST_ASSERT_EQUAL_DOUBLE(2*expected_visis[output], sum_visi_YY_real[output]);
      TEST_ASSERT_EQUAL_DOUBLE(0.0, sum_visi_XX_imag[output]);
      TEST_ASSERT_EQUAL_DOUBLE(0.0, sum_visi_XY_imag[output]);
      TEST_ASSERT_EQUAL_DOUBLE(0.0, sum_visi_YX_imag[output]);
      TEST_ASSERT_EQUAL_DOUBLE(0.0, sum_visi_YY_imag[output]);

    }
  }

  //This could be useful for future expansion
  // else if (beamtype == ANALY_DIPOLE || beamtype == GAUSS_BEAM) {
  //   for (int output = 0; output < num_visis; output++) {

  //     TEST_ASSERT_EQUAL_DOUBLE(expected_order[output], sum_visi_XX_real[output]);
  //     TEST_ASSERT_EQUAL_DOUBLE(0.0, sum_visi_XY_real[output]);
  //     TEST_ASSERT_EQUAL_DOUBLE(0.0, sum_visi_YX_real[output]);
  //     TEST_ASSERT_EQUAL_DOUBLE(expected_order[output], sum_visi_YY_real[output]);
  //     TEST_ASSERT_EQUAL_DOUBLE(0.0, sum_visi_XX_imag[output]);
  //     TEST_ASSERT_EQUAL_DOUBLE(0.0, sum_visi_XY_imag[output]);
  //     TEST_ASSERT_EQUAL_DOUBLE(0.0, sum_visi_YX_imag[output]);
  //     TEST_ASSERT_EQUAL_DOUBLE(0.0, sum_visi_YY_imag[output]);

  //   }
  // }
  // else {
  //   //If there is no beam, should just end up with the number of components
  //   //as we apply a gain of 1.0
  //   for (int output = 0; output < num_visis; output++) {

  //     TEST_ASSERT_EQUAL_DOUBLE(num_components, sum_visi_XX_real[output]);
  //     TEST_ASSERT_EQUAL_DOUBLE(0.0, sum_visi_XY_real[output]);
  //     TEST_ASSERT_EQUAL_DOUBLE(0.0, sum_visi_YX_real[output]);
  //     TEST_ASSERT_EQUAL_DOUBLE(num_components, sum_visi_YY_real[output]);
  //     TEST_ASSERT_EQUAL_DOUBLE(0.0, sum_visi_XX_imag[output]);
  //     TEST_ASSERT_EQUAL_DOUBLE(0.0, sum_visi_XY_imag[output]);
  //     TEST_ASSERT_EQUAL_DOUBLE(0.0, sum_visi_YX_imag[output]);
  //     TEST_ASSERT_EQUAL_DOUBLE(0.0, sum_visi_YY_imag[output]);

  //   }
  // }
  //
  // //Be free my beauties
  free( primay_beam_J00 );
  free( primay_beam_J01 );
  free( primay_beam_J10 );
  free( primay_beam_J11 );
  free( visi_components );
  free( flux_I );
  free( flux_Q );
  free( flux_U );
  free( flux_V );
  free( sum_visi_XX_real );
  free( sum_visi_XX_imag );
  free( sum_visi_XY_real );
  free( sum_visi_XY_imag );
  free( sum_visi_YX_real );
  free( sum_visi_YX_imag );
  free( sum_visi_YY_real );
  free( sum_visi_YY_imag );

}

/*
This test checks varying the gain with beamtype=FEE_BEAM
*/
void test_kern_update_sum_visis_VaryGainFEEBeam(void) {
  test_kern_update_sum_visis_VaryGainChooseBeams(FEE_BEAM);
}

/*
This test checks varying the measurement equation with beamtype=FEE_BEAM_INTERP
*/
void test_kern_update_sum_visis_VaryGainFEEInterpBeam(void) {
  test_kern_update_sum_visis_VaryGainChooseBeams(FEE_BEAM_INTERP);
}

// /*
// This test checks varying the gain with beamtype=ANALY_DIPOLE
// */
// void test_kern_update_sum_visis_VaryGainAnalyBeam(void) {
//   test_kern_update_sum_visis_VaryGainChooseBeams(ANALY_DIPOLE);
// }

// /*
// This test checks varying the gain with beamtype=GAUSS_BEAM
// */
// void test_kern_update_sum_visis_VaryGainGaussBeam(void) {
//   test_kern_update_sum_visis_VaryGainChooseBeams(GAUSS_BEAM);
// }

// /*
// This test checks varying the gain with beamtype=NO_BEAM
// */
// void test_kern_update_sum_visis_VaryGainNoBeam(void) {
//   test_kern_update_sum_visis_VaryGainChooseBeams(NO_BEAM);
// }



// /*
// This test checks varying the measurement equation with beamtype=MWA_ANALY
// */
// void test_kern_update_sum_visis_VaryGainMWAAnaly(void) {
//   test_kern_update_sum_visis_VaryGainChooseBeams(MWA_ANALY);
// }

//Run the test with unity
int main(void)
{
    UNITY_BEGIN();
    //Test while varying beam gain for all beam types
    RUN_TEST(test_kern_update_sum_visis_VaryGainFEEBeam);
    RUN_TEST(test_kern_update_sum_visis_VaryGainFEEInterpBeam);
    // RUN_TEST(test_kern_update_sum_visis_VaryGainAnalyBeam);
    // RUN_TEST(test_kern_update_sum_visis_VaryGainGaussBeam);
    // RUN_TEST(test_kern_update_sum_visis_VaryGainNoBeam);
    // RUN_TEST(test_kern_update_sum_visis_VaryGainMWAAnaly);

    return UNITY_END();
}
