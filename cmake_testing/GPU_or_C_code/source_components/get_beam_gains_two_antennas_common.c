#include "get_beam_gains_two_antennas_common.h"

#define UNITY_INCLUDE_FLOAT

/*
Test the __device__ code that grabs the beam gains works correctly
for all beam types
*/
void test_get_beam_gains_ChooseBeams(int beamtype, int do_gpu) {

  // #ifdef DOUBLE_PRECISION
  // printf("WODEN is using DOUBLE precision\n");
  // #else
  // printf("WODEN is using FLOAT precision\n");
  // #endif


  int num_ants = 3;
  int use_twoants = 1;

  int num_baselines = ((num_ants - 1)*num_ants) / 2;
  int num_times = 2;
  int num_freqs = 2;
  

  int num_visis = num_baselines*num_times*num_freqs;

  int num_components = 4;

  int num_input_gains = num_ants*num_components*num_times*num_freqs;

  user_precision_complex_t *primay_beam_J00 = malloc(num_input_gains*sizeof(user_precision_complex_t));
  user_precision_complex_t *primay_beam_J01 = malloc(num_input_gains*sizeof(user_precision_complex_t));
  user_precision_complex_t *primay_beam_J10 = malloc(num_input_gains*sizeof(user_precision_complex_t));
  user_precision_complex_t *primay_beam_J11 = malloc(num_input_gains*sizeof(user_precision_complex_t));

  //All models apart from NO_BEAM should have gains set other than 1.0
  int count = 0;
  if (beamtype != NO_BEAM) {
    for (int visi = 0; visi < num_ants*num_components*num_times*num_freqs; visi++) {
      primay_beam_J00[visi] = count + I*0.0;
      primay_beam_J11[visi] = count + I*0.0;
      count ++ ;
    }
  }

  //Only MWA beams have cross-pols
  count = 0;
  if (beamtype == FEE_BEAM || beamtype == FEE_BEAM_INTERP || beamtype == MWA_ANALY) {
    for (int visi = 0; visi < num_ants*num_components*num_times*num_freqs; visi++) {
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

  if (do_gpu == 1) {
    test_kern_get_beam_gains(num_freqs, num_visis,
            num_baselines, num_components, num_times, beamtype,
            primay_beam_J00, primay_beam_J01,
            primay_beam_J10, primay_beam_J11,
            recover_g1x, recover_D1x,
            recover_D1y, recover_g1y,
            recover_g2x, recover_D2x,
            recover_D2y, recover_g2y,
            use_twoants, num_ants);
  }

  //OK we're setting things up to have num_ants slowest changing, then num_times, num_freqs, num_components
  //`kern_get_beam_gains`, which is called by `test_kern_get_beam_gains`, is supposed
  //to output results in order of components slowest change, baselines fastest.
  //a.k.a. the first num_visis are all for the first component

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

  //iBaseline is all baselines all freqs all times; goes time, freq, baselines
  //in slowest to fastest changing order
  
  int comp_ind = 0;
  int iBaseline = 0; 
  int time_ind = 0;
  int freq_ind = 0;
  int baseline_ind = 0;
  int output_ind = 0;

  double *expected_gain1 = malloc(num_visis*num_components*sizeof(double));
  double *expected_gain2 = malloc(num_visis*num_components*sizeof(double));

  int num_coords = num_components*num_times*num_freqs;

  for (int visi = 0; visi < num_ants*num_components*num_times*num_freqs; visi++) {
    comp_ind = (int)floorf(visi / num_visis);
    iBaseline = (visi - (comp_ind*num_visis));

    time_ind = (int)floorf( iBaseline / (num_baselines * num_freqs));
    freq_ind = (int)floorf( (iBaseline - (time_ind*num_baselines * num_freqs)) / num_baselines);
    baseline_ind = iBaseline - (time_ind*num_baselines * num_freqs) - (freq_ind*num_baselines);

    output_ind = num_visis*comp_ind + iBaseline;

    //Given we just incremented gains, as we go slow->fast with antennas,time,freq,comp
    //this is what we would expect
    expected_gain1[output_ind] = ant1_indexes[baseline_ind]*num_coords + time_ind*num_components*num_freqs + freq_ind*num_components + comp_ind; 
    expected_gain2[output_ind] = ant2_indexes[baseline_ind]*num_coords + time_ind*num_components*num_freqs + freq_ind*num_components + comp_ind;

  }


  if (beamtype == FEE_BEAM || beamtype == FEE_BEAM_INTERP) {
    for (int output = 0; output < num_visis*num_components; output++) {

      // printf("output1, found1, expec1 %d %.1f %.1f\n", output,creal(recover_g1x[output]), expected_gain1[output]);
      // printf("output, found, expec %d %.1f %.1f\n", output,creal(recover_D1x[output]), expected_gain1[output]);
      // printf("output, found, expec %d %.1f %.1f\n", output,creal(recover_D1y[output]), expected_gain1[output]);
      // printf("output, found, expec %d %.1f %.1f\n", output,creal(recover_g1y[output]), expected_gain1[output]);

      // printf("output2, found2, expec2 %d %.1f %.1f\n", output,creal(recover_g2x[output]), expected_gain2[output]);
      // printf("output, found, expec %d %.1f %.1f\n", output,creal(recover_D2x[output]), expected_gain2[output]);
      // printf("output, found, expec %d %.1f %.1f\n", output,creal(recover_D2y[output]), expected_gain2[output]);
      // printf("output, found, expec %d %.1f %.1f\n", output,creal(recover_g2y[output]), expected_gain2[output]);

      TEST_ASSERT_EQUAL_FLOAT(expected_gain1[output], creal(recover_g1x[output]));
      TEST_ASSERT_EQUAL_FLOAT(expected_gain1[output], creal(recover_D1x[output]));
      TEST_ASSERT_EQUAL_FLOAT(expected_gain1[output], creal(recover_D1y[output]));
      TEST_ASSERT_EQUAL_FLOAT(expected_gain1[output], creal(recover_g1y[output]));
      //We explicitly set the imaginary part to 0
      TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(recover_g1x[output]));
      TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(recover_D1x[output]));
      TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(recover_D1y[output]));
      TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(recover_g1y[output]));

      TEST_ASSERT_EQUAL_FLOAT(expected_gain2[output], creal(recover_g2x[output]));
      TEST_ASSERT_EQUAL_FLOAT(expected_gain2[output], creal(recover_D2x[output]));
      TEST_ASSERT_EQUAL_FLOAT(expected_gain2[output], creal(recover_D2y[output]));
      TEST_ASSERT_EQUAL_FLOAT(expected_gain2[output], creal(recover_g2y[output]));
      //We explicitly set the imaginary part to 0
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
  free(expected_gain1);
  free(expected_gain2);

}