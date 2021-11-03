#include <unity.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "woden_settings.h"
// #include "create_sky_model.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

#define UNITY_INCLUDE_FLOAT

/*
Check the function that uses `woden_settings` to create the LSTs for the
simulation works as expected
*/
void check_setup_lsts_and_phase_centre(woden_settings_t *woden_settings) {

  //Function to be tested
  user_precision_t *lsts = setup_lsts_and_phase_centre(woden_settings);

  user_precision_t *expected_lsts = malloc(woden_settings->num_time_steps*sizeof(user_precision_t));

  for (int time_step = 0; time_step < woden_settings->num_time_steps ; time_step++) {
    expected_lsts[time_step] = woden_settings->lst_base + (time_step + 0.5)*woden_settings->time_res*SOLAR2SIDEREAL*DS2R;
  }

  int num_visis = woden_settings->num_baselines*woden_settings->num_time_steps*woden_settings->num_freqs;
  TEST_ASSERT_EQUAL_INT(num_visis, woden_settings->num_visis);

  TEST_ASSERT_EQUAL_FLOAT(sin(woden_settings->dec0), woden_settings->sdec0);
  TEST_ASSERT_EQUAL_FLOAT(cos(woden_settings->dec0), woden_settings->cdec0);

  TEST_ASSERT_EQUAL_FLOAT_ARRAY(expected_lsts, lsts, woden_settings->num_time_steps);

  free(lsts);
}

/*
Check things work with an LST of zero
*/
void check_setup_lsts_and_phase_centre_set1() {

  woden_settings_t *woden_settings = malloc(sizeof(woden_settings_t));

  woden_settings->num_baselines = 1;
  woden_settings->num_time_steps = 3;
  woden_settings->num_freqs = 2;
  woden_settings->time_res = 1.0;

  woden_settings->lst_base = 0.0;
  woden_settings->ra0 = 0.0;
  woden_settings->dec0 = M_PI/2;

  check_setup_lsts_and_phase_centre(woden_settings);

  free(woden_settings);


}

/*
Check things work with a non-zero LST
*/
void check_setup_lsts_and_phase_centre_set2() {

  woden_settings_t *woden_settings = malloc(sizeof(woden_settings_t));

  woden_settings->num_baselines = 1;
  woden_settings->num_time_steps = 8;
  woden_settings->num_freqs = 2;
  woden_settings->time_res = 1.0;

  woden_settings->lst_base = 1.2345;
  woden_settings->ra0 = 0.0;
  woden_settings->dec0 = M_PI/2;

  check_setup_lsts_and_phase_centre(woden_settings);

  free(woden_settings);

}

/*
Different non-zero LST, different dec0
*/
void check_setup_lsts_and_phase_centre_set3() {

  woden_settings_t *woden_settings = malloc(sizeof(woden_settings_t));

  woden_settings->num_baselines = 1;
  woden_settings->num_time_steps = 8;
  woden_settings->num_freqs = 2;
  woden_settings->time_res = 1.0;

  woden_settings->lst_base = 1.2345;
  woden_settings->ra0 = 0.0;
  woden_settings->dec0 = MWA_LAT_RAD;

  check_setup_lsts_and_phase_centre(woden_settings);

  free(woden_settings);

}

//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(check_setup_lsts_and_phase_centre_set1);
    RUN_TEST(check_setup_lsts_and_phase_centre_set2);
    RUN_TEST(check_setup_lsts_and_phase_centre_set3);

    return UNITY_END();
}
