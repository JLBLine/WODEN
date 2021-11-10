#include <unity.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "woden_struct_defs.h"
#include "visibility_set.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

#define UNITY_INCLUDE_FLOAT

/*
Check that this function sets a number of simulation settings in
`visibility_set` correctly
*/
void test_fill_timefreq_visibility() {

  //Dummy simulation settings
  woden_settings_t *woden_settings = malloc(sizeof(woden_settings_t));
  woden_settings->num_freqs = 3;
  woden_settings->frequency_resolution = 50e6;
  woden_settings->num_time_steps = 3;
  woden_settings->ra0 = M_PI;
  woden_settings->num_baselines = 2;
  woden_settings->num_visis = woden_settings->num_freqs*woden_settings->num_time_steps*woden_settings->num_baselines;

  user_precision_t base_band_freq = VELC / 2;

  double lsts[] = {0.0, M_PI/2, 2*M_PI};

  //Output container
  visibility_set_t *visibility_set = malloc(sizeof(visibility_set_t));

  //Code to be tested
  fill_timefreq_visibility_set(visibility_set, woden_settings,
                              base_band_freq, lsts);

  //Expected results
  double expec_lsts[] = {0, 0, 0, 0, 0, 0,
                        M_PI/2, M_PI/2, M_PI/2, M_PI/2, M_PI/2, M_PI/2,
                        2*M_PI, 2*M_PI, 2*M_PI, 2*M_PI, 2*M_PI, 2*M_PI};

  double expec_cha0s[] = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
                          0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                         -1.0, -1.0, -1.0, -1.0, -1.0, -1.0 };

  double expec_sha0s[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                        -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
                         0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  user_precision_t expec_wavelengths[] = {2.0, 2.0, 1.49974048, 1.49974048, 1.19966781, 1.19966781,
                               2.0, 2.0, 1.49974048, 1.49974048, 1.19966781, 1.19966781,
                               2.0, 2.0, 1.49974048, 1.49974048, 1.19966781, 1.19966781};

  //Check results are good
  for (size_t visi = 0; visi < woden_settings->num_visis; visi++) {
    // printf("%.8f %.8f %.8f %.8f\n", visibility_set->allsteps_cha0s[visi],
    //                                 visibility_set->allsteps_sha0s[visi],
    //                                 visibility_set->allsteps_lsts[visi],
    //                                 visibility_set->allsteps_wavelengths[visi]);

    TEST_ASSERT_FLOAT_WITHIN(1e-7, expec_sha0s[visi], visibility_set->allsteps_sha0s[visi]);
    TEST_ASSERT_FLOAT_WITHIN(1e-7, expec_cha0s[visi], visibility_set->allsteps_cha0s[visi]);
    TEST_ASSERT_FLOAT_WITHIN(1e-7, expec_lsts[visi], visibility_set->allsteps_lsts[visi]);
    TEST_ASSERT_FLOAT_WITHIN(1e-7, expec_wavelengths[visi], visibility_set->allsteps_wavelengths[visi]);
  }

  user_precision_t expec_freqs[] = {VELC / 2, VELC / 2 + 50e+6, VELC / 2 + 100e+6};

  TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_freqs, visibility_set->channel_frequencies, woden_settings->num_freqs);








}

//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_fill_timefreq_visibility);

    return UNITY_END();
}
