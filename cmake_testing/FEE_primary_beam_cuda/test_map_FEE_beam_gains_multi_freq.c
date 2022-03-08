#include <unity.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "constants.h"
#include "FEE_primary_beam.h"
// #include "FEE_primary_beam_cuda.h"
#include "woden_struct_defs.h"
// #include "multifreq_common.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

#define UNITY_INCLUDE_FLOAT

//CUDA code we are testing
extern void test_map_FEE_beam_gains_multi_freq(beam_settings_t *beam_settings,
                              user_precision_complex_t *primay_beam_J00,
                              user_precision_complex_t *primay_beam_J01,
                              user_precision_complex_t *primay_beam_J10,
                              user_precision_complex_t *primay_beam_J11,
                              int num_freqs, int num_components, int num_times);

/*
Test that the code that grabs the FEE beam gains out of the `num_freqs`
length array of `RTS_MWA_FEE_beam_t` when running multiple frequencies of the
MWA FEE beam code does the correct mapping into `primay_beam_J0*`
*/
void run_test_map_FEE_beam_gains_multi_freq(int num_freqs, int num_components,
                                            int num_times){

  //Initialise beam_settings
  beam_settings_t *beam_settings = malloc(sizeof(beam_settings_t));
  //This will hold the gain values on the GPU - gain values are assigned
  //by `test_map_FEE_beam_gains_multi_freq` on the device
  beam_settings->FEE_beams = malloc(num_freqs*sizeof(RTS_MWA_FEE_beam_t));

  //Containers for outputs
  user_precision_complex_t *primay_beam_J00 = malloc(num_freqs*num_times*num_components*sizeof(user_precision_complex_t));
  user_precision_complex_t *primay_beam_J01 = malloc(num_freqs*num_times*num_components*sizeof(user_precision_complex_t));
  user_precision_complex_t *primay_beam_J10 = malloc(num_freqs*num_times*num_components*sizeof(user_precision_complex_t));
  user_precision_complex_t *primay_beam_J11 = malloc(num_freqs*num_times*num_components*sizeof(user_precision_complex_t));

  //Test the CUDA code
  test_map_FEE_beam_gains_multi_freq(beam_settings,
    primay_beam_J00, primay_beam_J01,
    primay_beam_J10, primay_beam_J11,
    num_freqs, num_components, num_times);

  //Check answers are as expected
  int beam_ind;
  double expec_gain;

  for (int freq = 0; freq < num_freqs; freq++) {
    for (int time = 0; time < num_times; time++) {
      for (int comp = 0; comp < num_components; comp++) {
        beam_ind = num_freqs*num_components*time + num_components*freq + comp;
        //This is the way the gain is set by `test_map_FEE_beam_gains_multi_freq`
        expec_gain = time*num_freqs + freq;

        TEST_ASSERT_EQUAL_DOUBLE(expec_gain, creal(primay_beam_J00[beam_ind]));
        TEST_ASSERT_EQUAL_DOUBLE(expec_gain, cimag(primay_beam_J00[beam_ind]));
        TEST_ASSERT_EQUAL_DOUBLE(expec_gain, creal(primay_beam_J01[beam_ind]));
        TEST_ASSERT_EQUAL_DOUBLE(expec_gain, cimag(primay_beam_J01[beam_ind]));
        TEST_ASSERT_EQUAL_DOUBLE(expec_gain, creal(primay_beam_J10[beam_ind]));
        TEST_ASSERT_EQUAL_DOUBLE(expec_gain, cimag(primay_beam_J10[beam_ind]));
        TEST_ASSERT_EQUAL_DOUBLE(expec_gain, creal(primay_beam_J11[beam_ind]));
        TEST_ASSERT_EQUAL_DOUBLE(expec_gain, cimag(primay_beam_J11[beam_ind]));

      }
    }
  }


  free(beam_settings->FEE_beams);
  free(primay_beam_J00);
  free(primay_beam_J01);
  free(primay_beam_J10);
  free(primay_beam_J11);

}

/*
Run with one set of params
*/

void run_test_map_FEE_beam_gains_multi_freq_test1(){
  int num_freqs = 3;
  int num_components = 4;
  int num_times = 2;

  run_test_map_FEE_beam_gains_multi_freq(num_freqs, num_components, num_times);
}

/*
Run with another set of params
*/

void run_test_map_FEE_beam_gains_multi_freq_test2(){
  int num_freqs = 4;
  int num_components = 3;
  int num_times = 10;

  run_test_map_FEE_beam_gains_multi_freq(num_freqs, num_components, num_times);
}


void run_test_map_FEE_beam_gains_multi_freq_test3(){
  int num_freqs = 32;
  int num_components = 5000;
  int num_times = 56;

  run_test_map_FEE_beam_gains_multi_freq(num_freqs, num_components, num_times);
}


//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(run_test_map_FEE_beam_gains_multi_freq_test1);
    RUN_TEST(run_test_map_FEE_beam_gains_multi_freq_test2);
    RUN_TEST(run_test_map_FEE_beam_gains_multi_freq_test3);

    return UNITY_END();
}
