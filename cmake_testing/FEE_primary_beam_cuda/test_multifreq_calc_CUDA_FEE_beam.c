#include <unity.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <erfa.h>
#include <string.h>

#include "constants.h"
#include "FEE_primary_beam.h"
#include "woden_struct_defs.h"
#include "test_run_and_map_multifreq_calc_CUDA_FEE_beam.h"
#include "test_multifreq_RTS_FEE_beam.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

#define UNITY_INCLUDE_FLOAT

//CUDA code we are testing------------------------------------------------------
extern void multifreq_get_MWAFEE_normalisation(beam_settings_t *beam_settings);

extern void free_FEE_primary_beam_from_GPU(RTS_MWA_FEE_beam_t *primary_beam);

extern void test_multifreq_calc_CUDA_FEE_beam(beam_settings_t *beam_settings,
                            user_precision_t *azs, user_precision_t *zas,
                            int num_time_steps, double latitude,
                            int num_components, int rotation, int scaling,
                            user_precision_complex_t *all_FEE_beam_gains);


#ifdef DOUBLE_PRECISION
  double TOL = 1e-13;
#else
  double TOL = 3e-2;
#endif


//Testing code------------------------------------------------------------------
void run_test_multifreq_calc_CUDA_FEE_beam(user_precision_t *delays,
                  double base_low_freq, double freq_res, int num_freqs,
                  woden_settings_t *woden_settings,
                  double *expected_values){

  // //Initialise beam_settings
  beam_settings_t *beam_settings = malloc(sizeof(beam_settings_t));
  double *beam_freqs = malloc(num_freqs*sizeof(double));

  woden_settings->base_low_freq = base_low_freq;
  woden_settings->num_freqs = num_freqs;

  for(int i=0;i<16;i++) {
      woden_settings->FEE_ideal_delays[i] = delays[i];
  }

  for (int find = 0; find < num_freqs; find++) {
    beam_freqs[find] = base_low_freq + find*freq_res;
  }

  //Setup up the MWA FEE beams on the CPU
  multifreq_RTS_MWAFEEInit(beam_settings, woden_settings, beam_freqs);

  //Send them to GPU and calculate normalisations
  multifreq_get_MWAFEE_normalisation(beam_settings);

  int num_times = 1;
  //Rotate by parallatic angle and scale to zenith
  int rotation = 1;
  int scaling = 1;
  //Store the outputs from all frequencies, for all az, za, for all polarisations
  int num_gains = num_times*NUM_COMPS*MAX_POLS*beam_settings->num_MWAFEE;
  user_precision_complex_t *all_FEE_beam_gains = malloc(num_gains*sizeof(user_precision_complex_t));

  test_multifreq_calc_CUDA_FEE_beam(beam_settings, azs, zas, num_times,
                            MWA_LAT_RAD, NUM_COMPS, rotation, scaling,
                            all_FEE_beam_gains);

  //Check answers are as expected
  for (int freq_ind = 0; freq_ind < beam_settings->num_MWAFEE; freq_ind++) {
    int beam_ind, expec_ind;

    for (int comp = 0; comp < NUM_COMPS; comp++) {

        //Beam indexes are stored as complex
        beam_ind = MAX_POLS*(NUM_COMPS*freq_ind + comp);
        //Expected values stored as real/imag so double the index
        expec_ind = 2*beam_ind;

        // printf("%.4f %.4f\n", expected_values[expec_ind+0], creal(all_FEE_beam_gains[beam_ind+0]));

        TEST_ASSERT_DOUBLE_WITHIN(TOL, expected_values[expec_ind+0],
                                  creal(all_FEE_beam_gains[beam_ind+0]));
        TEST_ASSERT_DOUBLE_WITHIN(TOL, expected_values[expec_ind+1],
                                  cimag(all_FEE_beam_gains[beam_ind+0]));
        TEST_ASSERT_DOUBLE_WITHIN(TOL, expected_values[expec_ind+2],
                                  creal(all_FEE_beam_gains[beam_ind+1]));
        TEST_ASSERT_DOUBLE_WITHIN(TOL, expected_values[expec_ind+3],
                                  cimag(all_FEE_beam_gains[beam_ind+1]));
        TEST_ASSERT_DOUBLE_WITHIN(TOL, expected_values[expec_ind+4],
                                  creal(all_FEE_beam_gains[beam_ind+2]));
        TEST_ASSERT_DOUBLE_WITHIN(TOL, expected_values[expec_ind+5],
                                  cimag(all_FEE_beam_gains[beam_ind+2]));
        TEST_ASSERT_DOUBLE_WITHIN(TOL, expected_values[expec_ind+6],
                                  creal(all_FEE_beam_gains[beam_ind+3]));
        TEST_ASSERT_DOUBLE_WITHIN(TOL, expected_values[expec_ind+7],
                                  cimag(all_FEE_beam_gains[beam_ind+3]));
    }

    //Grab a beam at a frequency and free it
    RTS_MWA_FEE_beam_t *FEE_beam = &beam_settings->FEE_beams[freq_ind];
    RTS_MWA_FEE_beam_t *FEE_beam_zenith = &beam_settings->FEE_beam_zeniths[freq_ind];

    free_FEE_primary_beam_from_GPU(FEE_beam);

    RTS_freeHDFBeam(FEE_beam);
    RTS_freeHDFBeam(FEE_beam_zenith);

  }

  free(all_FEE_beam_gains);

}


void check_for_env_and_run_finetest(user_precision_t *delays,
                  double base_low_freq,
                  double freq_res, int num_freqs,
                  double *expected_values) {

  char* mwa_fee_hdf5 = getenv("MWA_FEE_HDF5_INTERP");
  woden_settings_t *woden_settings = malloc(sizeof(woden_settings_t));

  if (mwa_fee_hdf5) {
    printf("MWA_FEE_HDF5_INTERP: %s\n", mwa_fee_hdf5 );

    woden_settings->hdf5_beam_path = mwa_fee_hdf5;

    run_test_multifreq_calc_CUDA_FEE_beam(delays, base_low_freq,
                      freq_res, num_freqs, woden_settings,
                      expected_values);
  }
  else {
    printf("MWA_FEE_HDF5_INTERP not found - not running test_RTS_FEE_beam test");
  }
}

void test_delays1_freqs1(){
  check_for_env_and_run_finetest(delays1, FREQ1, FREQ_RES1, NUM_FREQS1,
                                 hyper_f1d1);
}

void test_delays2_freqs2(){
  check_for_env_and_run_finetest(delays2, FREQ2, FREQ_RES2, NUM_FREQS2,
                                 hyper_f2d2);
}


//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_delays1_freqs1);
    RUN_TEST(test_delays2_freqs2);

    return UNITY_END();
}
