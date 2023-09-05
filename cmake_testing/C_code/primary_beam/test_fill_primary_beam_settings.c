#include <unity.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "primary_beam.h"
#include "woden_struct_defs.h"
#include "woden_precision_defs.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

#ifdef DOUBLE_PRECISION
  double TOL = 1e-15;
#else
  double TOL = 1e-7;
#endif

double lsts[] = {1*DD2R, 120*DD2R, 240*DD2R};

/*
Fill in some example simulation settings
*/
woden_settings_t * make_woden_settings(){
  woden_settings_t *woden_settings = malloc(sizeof(woden_settings_t));

  woden_settings->gauss_dec_point = -30*DD2R;
  woden_settings->lst_base = 15*DD2R;
  woden_settings->gauss_ra_point = 30*DD2R;
  woden_settings->gauss_beam_FWHM = 20*DD2R;
  woden_settings->gauss_beam_ref_freq = 100e+6;
  woden_settings->num_time_steps = 3;
  woden_settings->latitude = MWA_LAT_RAD;

  return woden_settings;
}

/*
Call `fill_primary_beam_settings` with the given `woden_settings`
and sky model `src`
*/
void test_fill_primary_beam_settings(woden_settings_t *woden_settings) {

  //Function to be tested
  beam_settings_t *beam_settings = fill_primary_beam_settings(woden_settings,
                                                              lsts);

  //Both the Gaussian Beam or analytic MWA Beam need the HA/Dec coords for
  //all directions, so test both of those here
  if (woden_settings->beamtype == GAUSS_BEAM || woden_settings->beamtype == MWA_ANALY) {

    //Only GAUSS_BEAM should have these things set
    if (woden_settings->beamtype == GAUSS_BEAM) {
      TEST_ASSERT_EQUAL_INT(GAUSS_BEAM, beam_settings->beamtype );

      //Check major settings are copied across / calculated
      TEST_ASSERT_EQUAL_FLOAT(woden_settings->gauss_beam_FWHM * DD2R,
                              beam_settings->beam_FWHM_rad );
      TEST_ASSERT_EQUAL_DOUBLE(woden_settings->gauss_beam_ref_freq,
                               beam_settings->beam_ref_freq );

      TEST_ASSERT_DOUBLE_WITHIN(TOL, sin(woden_settings->gauss_dec_point),
                              beam_settings->gauss_sdec );
      TEST_ASSERT_DOUBLE_WITHIN(TOL, cos(woden_settings->gauss_dec_point),
                              beam_settings->gauss_cdec );
      TEST_ASSERT_EQUAL_DOUBLE(woden_settings->lst_base - woden_settings->gauss_ra_point,
                              beam_settings->gauss_ha );
    } else {
      TEST_ASSERT_EQUAL_INT(MWA_ANALY, beam_settings->beamtype );
    }

  }
  //For everything else, just check that the beamtype has been copied from
  //woden_settings to beam_settings
  else if (woden_settings->beamtype == FEE_BEAM){
    TEST_ASSERT_EQUAL_INT(FEE_BEAM, beam_settings->beamtype );
  }

  else if (woden_settings->beamtype == FEE_BEAM_INTERP ) {
    TEST_ASSERT_EQUAL_INT(FEE_BEAM_INTERP, beam_settings->beamtype );
  }

  else if (woden_settings->beamtype == ANALY_DIPOLE) {
    TEST_ASSERT_EQUAL_INT(ANALY_DIPOLE, beam_settings->beamtype );
  }

  else {
    TEST_ASSERT_EQUAL_INT(NO_BEAM, beam_settings->beamtype );
  }
  free(woden_settings);
}


/*
Test `fill_primary_beam_settings` when `beamtype = GAUSS_BEAM`
*/
void test_fill_primary_beam_settingsGaussBeam(void) {

  woden_settings_t *woden_settings = make_woden_settings();
  woden_settings->beamtype = GAUSS_BEAM;

  test_fill_primary_beam_settings(woden_settings);

}

/*
Test `fill_primary_beam_settings` when `beamtype = FEE_BEAM`
*/
void test_fill_primary_beam_settingsMWAFEEBeam(void) {

  woden_settings_t *woden_settings = make_woden_settings();
  woden_settings->beamtype = FEE_BEAM;

  test_fill_primary_beam_settings(woden_settings);

}

/*
Test `fill_primary_beam_settings` when `beamtype = ANALY_DIPOLE`
*/
void test_fill_primary_beam_settingsEDA2Beam(void) {

  woden_settings_t *woden_settings = make_woden_settings();
  woden_settings->beamtype = ANALY_DIPOLE;

  test_fill_primary_beam_settings(woden_settings);

}

/*
Test `fill_primary_beam_settings` when `beamtype = NO_BEAM`
*/
void test_fill_primary_beam_settingsNoBeam(void) {

  woden_settings_t *woden_settings = make_woden_settings();
  woden_settings->beamtype = NO_BEAM;

  test_fill_primary_beam_settings(woden_settings);

}

/*
Test `fill_primary_beam_settings` when `beamtype = FEE_BEAM_INTERP`
*/
void test_fill_primary_beam_settingsMWAFEEInterpBeam(void) {

  woden_settings_t *woden_settings = make_woden_settings();
  woden_settings->beamtype = FEE_BEAM_INTERP;

  test_fill_primary_beam_settings(woden_settings);

}

/*
Test `fill_primary_beam_settings` when `beamtype = MWA_ANALY`
*/
void test_fill_primary_beam_settingsMWAAnalyBeam(void) {

  woden_settings_t *woden_settings = make_woden_settings();
  woden_settings->beamtype = MWA_ANALY;

  test_fill_primary_beam_settings(woden_settings);

}


//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_fill_primary_beam_settingsGaussBeam);
    RUN_TEST(test_fill_primary_beam_settingsMWAAnalyBeam);
    RUN_TEST(test_fill_primary_beam_settingsMWAFEEBeam);
    RUN_TEST(test_fill_primary_beam_settingsEDA2Beam);
    RUN_TEST(test_fill_primary_beam_settingsNoBeam);
    RUN_TEST(test_fill_primary_beam_settingsMWAFEEInterpBeam);

    return UNITY_END();
}
