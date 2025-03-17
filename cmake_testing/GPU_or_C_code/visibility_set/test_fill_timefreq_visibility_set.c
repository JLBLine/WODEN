#include <unity.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "woden_struct_defs.h"
#include "visibility_set.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

// #define UNITY_INCLUDE_FLOAT
// #define UNITY_INCLUDE_DOUBLE

/*
Check that this function sets a number of simulation settings in
`visibility_set` correctly
*/
void test_fill_timefreq_visibility() {

  //Dummy simulation settings
  woden_settings_t *woden_settings = malloc(sizeof(woden_settings_t));
  woden_settings->num_freqs = 3;
  woden_settings->frequency_resolution = (VELC / 4.0);
  woden_settings->num_time_steps = 3;
  woden_settings->ra0 = 0.0;
  woden_settings->num_baselines = 2;
  woden_settings->num_visis = woden_settings->num_freqs*woden_settings->num_time_steps*woden_settings->num_baselines;

  double base_band_freq = 149896229.0;

  double lsts[] = {0.0, M_PI/6.0, M_PI/4.0};

  double expec_lsts[] = {0, 0, 0, 0, 0, 0,
                         M_PI/6.0,  M_PI/6.0,  M_PI/6.0,  M_PI/6.0,  M_PI/6.0,  M_PI/6.0,
                         M_PI/4.0, M_PI/4.0, M_PI/4.0, M_PI/4.0, M_PI/4.0, M_PI/4.0};

  double expec_sha0s[] = {0, 0, 0, 0, 0, 0,
                          0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                         sqrt(2)/2.0, sqrt(2)/2.0, sqrt(2)/2.0, sqrt(2)/2.0, sqrt(2)/2.0, sqrt(2)/2.0 };

  double expec_cha0s[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                          sqrt(3.0)/2.0, sqrt(3.0)/2.0, sqrt(3.0)/2.0, sqrt(3.0)/2.0, sqrt(3.0)/2.0, sqrt(3.0)/2.0,
                          sqrt(2)/2.0, sqrt(2)/2.0, sqrt(2)/2.0, sqrt(2)/2.0, sqrt(2)/2.0, sqrt(2)/2.0 };

  double expec_wavelengths[] = {2.0, 2.0, 1.3333333333333333, 1.3333333333333333, 1.0, 1.0,
                               2.0, 2.0, 1.3333333333333333, 1.3333333333333333, 1.0, 1.0,
                               2.0, 2.0, 1.3333333333333333, 1.3333333333333333, 1.0, 1.0};

  //Output container
  visibility_set_t *visibility_set = malloc(sizeof(visibility_set_t));

  //Code to be tested
  fill_timefreq_visibility_set(visibility_set, woden_settings,
                              base_band_freq, lsts);
  #ifdef DOUBLE_PRECISION
    double TOL = 1e-15;
  #else
    double TOL = 1e-7;
  #endif

  //Check results are good
  for (size_t visi = 0; visi < woden_settings->num_visis; visi++) {
    // printf("%.8f %.8f %.8f %.8f\n", visibility_set->allsteps_cha0s[visi],
    //                                 visibility_set->allsteps_sha0s[visi],
    //                                 visibility_set->allsteps_lsts[visi],
    //                                 visibility_set->allsteps_wavelengths[visi]);

    // printf("%.16f %.16f\n", expec_wavelengths[visi], visibility_set->allsteps_wavelengths[visi] );

    TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_wavelengths[visi], visibility_set->allsteps_wavelengths[visi]);

    TEST_ASSERT_DOUBLE_WITHIN(1e-15, expec_sha0s[visi], visibility_set->allsteps_sha0s[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(1e-15, expec_cha0s[visi], visibility_set->allsteps_cha0s[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(1e-15, expec_lsts[visi], visibility_set->allsteps_lsts[visi]);

  }

  double expec_freqs[] = {VELC / 2, (3*VELC) / 4, VELC};

  for (int freq = 0; freq < woden_settings->num_freqs; freq++) {
    TEST_ASSERT_DOUBLE_WITHIN(1e-15, expec_freqs[freq], visibility_set->channel_frequencies[freq]);
  }

}

//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_fill_timefreq_visibility);

    return UNITY_END();
}
