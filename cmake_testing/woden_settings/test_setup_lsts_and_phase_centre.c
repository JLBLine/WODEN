#include <unity.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "woden_settings.h"
// #include "create_sky_model.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

#ifdef DOUBLE_PRECISION
  double TOL = 1e-15;
#else
  double TOL = 1e-7;
#endif

/*
Check the function that uses `woden_settings` to create the LSTs for the
simulation works as expected
*/
void check_setup_lsts_and_phase_centre(woden_settings_t *woden_settings,
                                       double *expec_prec_lsts) {

  //Do it without precession
  woden_settings->do_precession = 0;

  //Function to be tested
  double *lsts = setup_lsts_and_phase_centre(woden_settings);

  double *expected_lsts = malloc(woden_settings->num_time_steps*sizeof(double));

  for (int time_step = 0; time_step < woden_settings->num_time_steps ; time_step++) {
    expected_lsts[time_step] = woden_settings->lst_base + (time_step + 0.5)*woden_settings->time_res*SOLAR2SIDEREAL*DS2R;
  }

  TEST_ASSERT_DOUBLE_WITHIN(TOL, sin(woden_settings->dec0), woden_settings->sdec0);
  TEST_ASSERT_DOUBLE_WITHIN(TOL, cos(woden_settings->dec0), woden_settings->cdec0);

  TEST_ASSERT_EQUAL_DOUBLE_ARRAY(expected_lsts, lsts, woden_settings->num_time_steps);

  free(lsts);

  //Do it all again but with precession
  woden_settings->do_precession = 1;

  //Function to be tested
  double *lsts_precess = setup_lsts_and_phase_centre(woden_settings);

  TEST_ASSERT_DOUBLE_WITHIN(TOL, sin(woden_settings->dec0), woden_settings->sdec0);
  TEST_ASSERT_DOUBLE_WITHIN(TOL, cos(woden_settings->dec0), woden_settings->cdec0);

  TEST_ASSERT_EQUAL_DOUBLE_ARRAY(expec_prec_lsts, lsts_precess, woden_settings->num_time_steps);

  free(lsts_precess);

  double *expected_mjd = malloc(woden_settings->num_time_steps*sizeof(double));
  double sec_to_d = 1 / (24.0*60.0*60.0);
  double mjd = woden_settings->jd_date - 2400000.5;

  for (int time_step = 0; time_step < woden_settings->num_time_steps ; time_step++) {
    expected_mjd[time_step] = mjd + (time_step + 0.5)*woden_settings->time_res*sec_to_d;
    // printf("%.8f %.8f\n",expected_mjd[time_step], woden_settings->mjds[time_step] );
  }

  TEST_ASSERT_EQUAL_DOUBLE_ARRAY(expected_mjd, woden_settings->mjds, woden_settings->num_time_steps);

  free(expected_mjd);

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
  woden_settings->latitude = MWA_LAT_RAD;
  woden_settings->latitude_obs_epoch_base = MWA_LAT_RAD;

  woden_settings->lst_base = 0.0;
  woden_settings->lst_obs_epoch_base = 0.0;
  woden_settings->ra0 = 0.0;
  woden_settings->dec0 = M_PI/2;
  woden_settings->jd_date = 2457278.201145833;

  double expec_prec_lsts[3] = { 6.2797277312015503, 6.2798007086408179,
                                6.2798736860799780 };

  check_setup_lsts_and_phase_centre(woden_settings, expec_prec_lsts);

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
  woden_settings->latitude = MWA_LAT_RAD;
  woden_settings->latitude_obs_epoch_base = MWA_LAT_RAD;

  woden_settings->lst_base = 1.2345;
  woden_settings->lst_obs_epoch_base = 1.2345;
  woden_settings->ra0 = 0.0;
  woden_settings->dec0 = M_PI/2;
  woden_settings->jd_date = 2457278.201145833;

  double expec_prec_lsts[8] = { 1.2317545128009186, 1.2318274511154514,
                                1.2319003894260832, 1.2319733277328140,
                                1.2320462660356444, 1.2321192043345734,
                                1.2321921426296016, 1.2322650809207287 };

  check_setup_lsts_and_phase_centre(woden_settings, expec_prec_lsts);

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
  woden_settings->latitude = MWA_LAT_RAD;
  woden_settings->latitude_obs_epoch_base = MWA_LAT_RAD;

  woden_settings->lst_base = 1.2345;
  woden_settings->lst_obs_epoch_base = 1.2345;
  woden_settings->ra0 = 0.0;
  woden_settings->dec0 = MWA_LAT_RAD;
  woden_settings->jd_date = 2457278.201145833;

  double expec_prec_lsts[8] = { 1.2317545128009186, 1.2318274511154514,
                                1.2319003894260832, 1.2319733277328140,
                                1.2320462660356444, 1.2321192043345734,
                                1.2321921426296016, 1.2322650809207287 };

  check_setup_lsts_and_phase_centre(woden_settings, expec_prec_lsts);

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
