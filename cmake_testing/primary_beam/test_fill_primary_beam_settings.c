#include <unity.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "primary_beam.h"
#include "woden_struct_defs.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

#define UNITY_INCLUDE_FLOAT

/*
DESCRIPTION HERE
*/

float point_ras[] = {0 , 5*DD2R, 10*DD2R};
float point_decs[] = {MWA_LAT_RAD, MWA_LAT_RAD, MWA_LAT_RAD};

float gauss_ras[] = {115*DD2R , 125*DD2R, 130*DD2R};
float gauss_decs[] = {-15*DD2R , -20*DD2R, -25*DD2R};

float shape_ras[] = {235*DD2R , 250*DD2R, 265*DD2R};
float shape_decs[] = {-35*DD2R , -40*DD2R, -45*DD2R};

float expec_point_sin_para[] = {-0.0039216, -0.6142135, -0.6142136, -0.0156904,
                                -0.5764040, -0.6534451, -0.0353442, -0.5401087,
                                -0.6939089 };

float expec_point_cos_para[] = {-0.9999923, -0.7891399, 0.7891398, 0.9998769,
                                -0.8171648, 0.7569739, 0.9993752, -0.8415953,
                                 0.7200627 };

float expec_gauss_sin_para[] = {-0.5432732, -0.9341043, -0.6122640, -0.6251023,
                                -0.8344397, -0.5628954, -0.6805615, -0.2228374,
                                -0.5384376 };

float expec_gauss_cos_para[] = {0.8395559, -0.3570003, -0.7906534, 0.7805428,
                                0.5510992, -0.8265281, 0.7326910, 0.9748556,
                                -0.8426654 };

float expec_shape_sin_para[] = {-0.6794566, -0.5854309, 0.8773373, -0.5505476,
                                -0.7216252, 0.8188159, -0.3965520, -0.8523725,
                                 0.5593383 };

float expec_shape_cos_para[] = {-0.7337157, 0.8107223, -0.4798742, -0.8348038,
                                 0.6922840, 0.5740561, -0.9180123, 0.5229351,
                                 0.8289395 };
/*
Make the polpulated catsource_t struct. Stick in some necessary values
*/
catsource_t * make_sky_model(void) {

  catsource_t *src = malloc(sizeof(catsource_t));

  src->n_comps = 9;
  src->n_points = 3;
  src->n_gauss = 3;
  src->n_shapes = 3;

  src->point_ras = point_ras;
  src->point_decs = point_decs;

  src->gauss_ras = gauss_ras;
  src->gauss_decs = gauss_decs;

  src->shape_ras = shape_ras;
  src->shape_decs = shape_decs;

  return src;
}

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
  //Make sky model
  catsource_t *src = make_sky_model();

  float lsts[] = {1*DD2R, 120*DD2R, 240*DD2R};

  //Function to be tested
  beam_settings_t *beam_settings = fill_primary_beam_settings(woden_settings,
                                                              src, lsts);

  //If running a Gaussian Beam simulation
  if (woden_settings->beamtype == GAUSS_BEAM) {
    TEST_ASSERT_EQUAL_INT(GAUSS_BEAM, beam_settings->beamtype );

    //Check major settings are copied across / calculated
    TEST_ASSERT_EQUAL_FLOAT(woden_settings->gauss_beam_FWHM * DD2R,
                            beam_settings->beam_FWHM_rad );
    TEST_ASSERT_EQUAL_FLOAT(woden_settings->gauss_beam_ref_freq,
                            beam_settings->beam_ref_freq );

    TEST_ASSERT_EQUAL_FLOAT(sinf(woden_settings->gauss_dec_point),
                            beam_settings->gauss_sdec );
    TEST_ASSERT_EQUAL_FLOAT(cosf(woden_settings->gauss_dec_point),
                            beam_settings->gauss_cdec );
    TEST_ASSERT_EQUAL_FLOAT(woden_settings->lst_base - woden_settings->gauss_ra_point,
                            beam_settings->gauss_ha );

    //Test the hour angle / decs are calculated correctly
    //Loop over all time and point components and calculate ha
    for (int component = 0; component < src->n_points; component++) {
      for ( int time_step = 0; time_step < woden_settings->num_time_steps; time_step++ ) {
        int step = component*woden_settings->num_time_steps + time_step;

        TEST_ASSERT_EQUAL_FLOAT(src->point_decs[component],
                                src->point_gaussbeam_decs[step]);
        TEST_ASSERT_EQUAL_FLOAT(lsts[time_step] - src->point_ras[component],
                                src->point_gaussbeam_has[step]);

        TEST_ASSERT_EQUAL_FLOAT(src->gauss_decs[component],
                                src->gauss_gaussbeam_decs[step]);
        TEST_ASSERT_EQUAL_FLOAT(lsts[time_step] - src->gauss_ras[component],
                                src->gauss_gaussbeam_has[step]);

        TEST_ASSERT_EQUAL_FLOAT(src->shape_decs[component],
                                src->shape_gaussbeam_decs[step]);
        TEST_ASSERT_EQUAL_FLOAT(lsts[time_step] - src->shape_ras[component],
                                src->shape_gaussbeam_has[step]);
      }
    }
  } else if (woden_settings->beamtype == FEE_BEAM) {
    // printf("HERE %d %d\n",woden_settings->beamtype, beam_settings->beamtype );
    TEST_ASSERT_EQUAL_INT(FEE_BEAM, beam_settings->beamtype);

    TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_point_sin_para,
                                  src->sin_point_para_angs, 9);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_point_cos_para,
                                  src->cos_point_para_angs, 9);

    TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_gauss_sin_para,
                                  src->sin_gauss_para_angs, 9);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_gauss_cos_para,
                                  src->cos_gauss_para_angs, 9);

    TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_shape_sin_para,
                                  src->sin_shape_para_angs, 9);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_shape_cos_para,
                                  src->cos_shape_para_angs, 9);

    // for (int para = 0; para < 9; para++) {
    //   printf("%.7f %.7f\n",src->sin_shape_para_angs[para],
    //                        src->cos_shape_para_angs[para] );
    // }
  }

  else if (woden_settings->beamtype == ANALY_DIPOLE) {
    TEST_ASSERT_EQUAL_INT(ANALY_DIPOLE, beam_settings->beamtype );
  }

  else {
    TEST_ASSERT_EQUAL_INT(NO_BEAM, beam_settings->beamtype );
  }
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



//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    // RUN_TEST(test_fill_primary_beam_settingsGaussBeam);
    RUN_TEST(test_fill_primary_beam_settingsMWAFEEBeam);

    return UNITY_END();
}