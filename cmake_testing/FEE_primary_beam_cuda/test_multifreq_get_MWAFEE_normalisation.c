#include <unity.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "constants.h"
#include "FEE_primary_beam.h"
// #include "FEE_primary_beam_cuda.h"
#include "woden_struct_defs.h"
#include "test_run_and_map_multifreq_calc_CUDA_FEE_beam.h"
#include "test_multifreq_get_MWAFEE_normalisation.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

//CUDA code we are testing
extern void multifreq_get_MWAFEE_normalisation(beam_settings_t *beam_settings);

extern void free_FEE_primary_beam_from_GPU(RTS_MWA_FEE_beam_t *primary_beam);

#ifdef DOUBLE_PRECISION
  double TOL = 1e-11;
#else
  double TOL = 3e-2;
#endif

/*
Run the code that gets the normalisation for the interpolated beam for
multiple frequencies
*/
void test_multifreq_get_MWAFEE_normalisation(user_precision_t *delays,
                  double base_low_freq, double freq_res, int num_freqs,
                  woden_settings_t *woden_settings,
                  user_precision_t *expec_gx, user_precision_t *expec_Dx,
                  user_precision_t *expec_Dy, user_precision_t *expec_gy){

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

  for (int freq_ind = 0; freq_ind < beam_settings->num_MWAFEE; freq_ind++) {

    RTS_MWA_FEE_beam_t *FEE_beam = &beam_settings->FEE_beams[freq_ind];
    RTS_MWA_FEE_beam_t *FEE_beam_zenith = &beam_settings->FEE_beam_zeniths[freq_ind];

    // printf("MULTI %.12f %.12f %.12f %.12f\n",creal(FEE_beam->norm_fac[0]), creal(FEE_beam->norm_fac[1]), creal(FEE_beam->norm_fac[2]),creal(FEE_beam->norm_fac[3]) );

    //Just test the real values, full integration test later will tell
    //us more, real is enough for testing here
    TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_gx[freq_ind],
                              creal(FEE_beam->norm_fac[0]));
    TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_Dx[freq_ind],
                              creal(FEE_beam->norm_fac[1]));
    TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_Dy[freq_ind],
                              creal(FEE_beam->norm_fac[2]));
    TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_gy[freq_ind],
                              creal(FEE_beam->norm_fac[3]));

    //Free all the things
    free_FEE_primary_beam_from_GPU(FEE_beam);

    RTS_freeHDFBeam(FEE_beam);
    RTS_freeHDFBeam(FEE_beam_zenith);

  }
}

void check_for_env_and_run_finetest(user_precision_t *delays,
                  double base_low_freq,
                  double freq_res, int num_freqs,
                  user_precision_t *expec_gx, user_precision_t *expec_Dx,
                  user_precision_t *expec_Dy, user_precision_t *expec_gy) {

  char* mwa_fee_hdf5 = getenv("MWA_FEE_HDF5_INTERP");

  woden_settings_t *woden_settings = malloc(sizeof(woden_settings_t));

  if (mwa_fee_hdf5) {
    printf("MWA_FEE_HDF5_INTERP: %s\n", mwa_fee_hdf5 );

    woden_settings->hdf5_beam_path = mwa_fee_hdf5;

    test_multifreq_get_MWAFEE_normalisation(delays, base_low_freq,
                      freq_res, num_freqs, woden_settings,
                      expec_gx, expec_Dx,
                      expec_Dy, expec_gy);

  }
  else {
    printf("MWA_FEE_HDF5 not found - not running test_RTS_FEE_beam test");
  }
}

void test_multifreq_get_MWAFEE_normalisation_delays1_freqs1(){
  check_for_env_and_run_finetest(delays1, FREQ1, FREQ_RES1, NUM_FREQS1,
                                 expec_gx_1, expec_Dx_1,
                                 expec_Dy_1, expec_gy_1);
}

void test_multifreq_get_MWAFEE_normalisation_delays2_freqs2(){
  check_for_env_and_run_finetest(delays2, FREQ2, FREQ_RES2, NUM_FREQS2,
                                 expec_gx_2, expec_Dx_2,
                                 expec_Dy_2, expec_gy_2);
}


//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_multifreq_get_MWAFEE_normalisation_delays1_freqs1);
    RUN_TEST(test_multifreq_get_MWAFEE_normalisation_delays2_freqs2);

    return UNITY_END();
}
