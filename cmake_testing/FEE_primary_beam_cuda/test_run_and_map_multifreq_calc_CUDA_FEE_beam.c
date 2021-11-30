#include <unity.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

// #include "constants.h"
#include "FEE_primary_beam.h"
// #include "FEE_primary_beam_cuda.h"
// #include "woden_struct_defs.h"
// #include "multifreq_common.h"

#include "test_run_and_map_multifreq_calc_CUDA_FEE_beam.h"
#include "test_multifreq_RTS_FEE_beam.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

//CUDA code we are testing
extern void multifreq_get_MWAFEE_normalisation(beam_settings_t *beam_settings);

extern void free_FEE_primary_beam_from_GPU(RTS_MWA_FEE_beam_t *primary_beam);

//CUDA code we are testing
extern void test_run_and_map_multifreq_calc_CUDA_FEE_beam(beam_settings_t *beam_settings,
    user_precision_t *azs, user_precision_t *zas, double latitude,
    user_precision_complex_t *primay_beam_J00,
    user_precision_complex_t *primay_beam_J01,
    user_precision_complex_t *primay_beam_J10,
    user_precision_complex_t *primay_beam_J11,
    int num_freqs, int NUM_COMPS, int num_times,
    int rotation, int scaling);

#ifdef DOUBLE_PRECISION
  double TOL = 1e-13;
#else
  double TOL = 3e-2;
#endif

/*
Run the multiple frequency MWA FEE beam from start to finish, where start is
the az,za directions for all time steps, we calculate the beam in these
directions for all requested frequencies and time, then map them into
single polarisation arrays, and grab off the GPU.
*/
void run_test_run_and_map_multifreq_calc_CUDA_FEE_beam(int freq_int,
                  int delay_int, user_precision_t *delays,
                  double base_low_freq, double freq_res, int num_freqs,
                  woden_settings_t *woden_settings,
                  double *expected_values){

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
  multifreq_RTS_MWAFEEInit(beam_settings,  woden_settings, beam_freqs);

  //Send them to GPU and calculate normalisations
  multifreq_get_MWAFEE_normalisation(beam_settings);

  //Test is setup for a single time step
  int num_times = 1;

  //Containers for outputs
  user_precision_complex_t *primay_beam_J00 = malloc(num_freqs*NUM_COMPS*sizeof(user_precision_complex_t));
  user_precision_complex_t *primay_beam_J01 = malloc(num_freqs*NUM_COMPS*sizeof(user_precision_complex_t));
  user_precision_complex_t *primay_beam_J10 = malloc(num_freqs*NUM_COMPS*sizeof(user_precision_complex_t));
  user_precision_complex_t *primay_beam_J11 = malloc(num_freqs*NUM_COMPS*sizeof(user_precision_complex_t));

  int rotation = 1;
  int scaling = 1;

  //Test the CUDA code
  test_run_and_map_multifreq_calc_CUDA_FEE_beam(beam_settings,
    azs, zas, MWA_LAT_RAD,
    primay_beam_J00, primay_beam_J01,
    primay_beam_J10, primay_beam_J11,
    num_freqs, NUM_COMPS, num_times,
    rotation, scaling);

  //Check answers are as expected
  int beam_ind, expec_ind;

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

  FILE *output_text;

  char filename[100];
  sprintf(filename, "MWA_FEE_multifreq_gains_delays%01d_freqs%01d.txt", delay_int, freq_int);

  output_text = fopen(filename,"w");

  int count = 0;


  for (int freq = 0; freq < num_freqs; freq ++) {
    for (int comp = 0; comp < NUM_COMPS; comp ++) {

      int beam_ind = NUM_COMPS*freq + comp;

      fprintf(output_text,"%.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.1f\n",
           creal(primay_beam_J00[beam_ind]), cimag(primay_beam_J00[beam_ind]),
           creal(primay_beam_J01[beam_ind]), cimag(primay_beam_J01[beam_ind]),
           creal(primay_beam_J10[beam_ind]), cimag(primay_beam_J10[beam_ind]),
           creal(primay_beam_J11[beam_ind]), cimag(primay_beam_J11[beam_ind]),
       beam_freqs[freq] );

       count ++;
    }
  }

  fflush(output_text);
  fclose(output_text);

  for (int freq_ind = 0; freq_ind < beam_settings->num_MWAFEE; freq_ind++) {

    RTS_MWA_FEE_beam_t *FEE_beam = &beam_settings->FEE_beams[freq_ind];
    RTS_MWA_FEE_beam_t *FEE_beam_zenith = &beam_settings->FEE_beam_zeniths[freq_ind];

    free_FEE_primary_beam_from_GPU(FEE_beam);

    RTS_freeHDFBeam(FEE_beam);
    RTS_freeHDFBeam(FEE_beam_zenith);

  }

  free(primay_beam_J00);
  free(primay_beam_J01);
  free(primay_beam_J10);
  free(primay_beam_J11);

}

void check_for_env_and_run_finetest(int freq_int, int delay_int,
                  user_precision_t *delays, double base_low_freq,
                  double freq_res, int num_freqs,
                  double *expected_values) {

  char* mwa_fee_hdf5 = getenv("MWA_FEE_HDF5_INTERP");
  woden_settings_t *woden_settings = malloc(sizeof(woden_settings_t));

  if (mwa_fee_hdf5) {
    printf("MWA_FEE_HDF5_INTERP: %s\n", mwa_fee_hdf5 );

    woden_settings->hdf5_beam_path = mwa_fee_hdf5;

    run_test_run_and_map_multifreq_calc_CUDA_FEE_beam(freq_int, delay_int,
                      delays, base_low_freq,
                      freq_res, num_freqs, woden_settings,
                      expected_values);
  }
  else {
    printf("MWA_FEE_HDF5_INTERP not found - not running test_RTS_FEE_beam test");
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

//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_delays1_freqs1);
    RUN_TEST(test_delays2_freqs2);
    RUN_TEST(test_delays3_freqs3);

    return UNITY_END();
}
