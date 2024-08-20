/*
`calculate_visibilities::calculate_visibilities` is the gateway function
to all CUDA functionality in WODEN. We'll test here in one baseline, frequency,
and time configuration. We'll vary the sky model and the primary beam. By
sticking all COMPONENTs at phase centre, we can just sum the expected fluxes
in XX / YY real to check things are being lanuched.

More variations like different phase centres / array configs etc are tested
in different test suites, so really just test that the correct CUDA functions
are launched by calculate_visibilities::calculate_visibilities`
*/

#include <math.h>
#include <unity.h>
#include <stdlib.h>
#include <complex.h>

#include "test_calculate_visibilities_common.h"

#ifdef DOUBLE_PRECISION
  double TOL = 1e-8;
#else
  double TOL = 5e-5;
#endif

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */


/*
For a numbre of COMPONENTs at phase centre, with no beam model
should just have a single gain in the XX and YY
real visis
*/
void test_comp_phase_centre_nobeam(int num_comps,
                                   visibility_set_t *visibility_set,
                                   woden_settings_t *woden_settings) {

  int num_cross = woden_settings->num_cross;

  double _Complex gain = 1.0 + I*0.0;
  double _Complex leak = 0.0 + I*0.0;

  double _Complex xx, xy, yx, yy;

  predict_inst_stokes(num_comps, gain, leak, leak, gain, gain, leak, leak, gain,
                        &xx, &xy, &yx, &yy);

  // printf("%.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n", creal(xx), cimag(xx),
  //          creal(xy), cimag(xy), creal(yx), cimag(yx), creal(yy), cimag(yy));

  // //First up check that the cross-correlations are as expected
  for (int visi = 0; visi < num_cross; visi++) {

  // printf("%d %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n",
  //           visi,
  //           visibility_set->sum_visi_XX_real[visi],
  //           visibility_set->sum_visi_XX_imag[visi],
  //           visibility_set->sum_visi_XY_real[visi],
  //           visibility_set->sum_visi_XY_imag[visi],
  //           visibility_set->sum_visi_YX_real[visi],
  //           visibility_set->sum_visi_YX_imag[visi],
  //           visibility_set->sum_visi_YY_imag[visi],
  //           visibility_set->sum_visi_YY_real[visi] );

    TEST_ASSERT_DOUBLE_WITHIN(TOL, creal(xx), visibility_set->sum_visi_XX_real[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, cimag(xx), visibility_set->sum_visi_XX_imag[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, creal(xy), visibility_set->sum_visi_XY_real[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, cimag(xy), visibility_set->sum_visi_XY_imag[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, creal(yx), visibility_set->sum_visi_YX_real[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, cimag(yx), visibility_set->sum_visi_YX_imag[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, creal(yy), visibility_set->sum_visi_YY_real[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, cimag(yy), visibility_set->sum_visi_YY_imag[visi]);
  }

  //Next, if there are auto-correlations, check they are correct
  //For the no beam case with sources all at phase centre they are the same
  //as the cross-correlations
  for (int visi = 0; visi < woden_settings->num_autos; visi++) {

    TEST_ASSERT_DOUBLE_WITHIN(TOL, creal(xx),
                             visibility_set->sum_visi_XX_real[num_cross+visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, cimag(xx),
                             visibility_set->sum_visi_XX_imag[num_cross+visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, creal(xy),
                             visibility_set->sum_visi_XY_real[num_cross+visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, cimag(xy),
                             visibility_set->sum_visi_XY_imag[num_cross+visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, creal(yx),
                             visibility_set->sum_visi_YX_real[num_cross+visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, cimag(yx),
                             visibility_set->sum_visi_YX_imag[num_cross+visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, creal(yy),
                             visibility_set->sum_visi_YY_imag[num_cross+visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, cimag(yy),
                             visibility_set->sum_visi_YY_real[num_cross+visi]);
  }
}



// #define UNITY_INCLUDE_FLOAT

void test_calculate_visibilities_NoBeam(int n_points, int n_gauss, int n_shapes,
                                        int num_sources) {

  source_catalogue_t *cropped_sky_models = make_cropped_sky_models(RA0, -0.46606083776035967,
                                                    n_points, n_gauss, n_shapes,
                                                    num_sources);

  beam_settings_t *beam_settings = malloc(sizeof(beam_settings_t));
  beam_settings->beamtype = NO_BEAM;

  woden_settings_t *woden_settings = make_woden_settings(RA0, -0.46606083776035967);
  woden_settings->beamtype = NO_BEAM;

  // woden_settings->do_autos = 1;
  // woden_settings->num_autos = NUM_CROSS;
  // woden_settings->num_visis = woden_settings->num_cross + woden_settings->num_autos;

  printf("We have this many visis %d %d %d\n",woden_settings->num_visis,woden_settings->num_autos,woden_settings->num_cross );

  visibility_set_t *visibility_set = test_calculate_visibilities(cropped_sky_models,
                                          beam_settings, woden_settings, RA0, -0.46606083776035967,
                                          beam_settings->beamtype);

  int num_comps = (n_points + n_gauss + n_shapes)*num_sources;
  test_comp_phase_centre_nobeam(num_comps, visibility_set, woden_settings);

  free_visi_set_inputs(visibility_set);
  free_visi_set_outputs(visibility_set);

  woden_settings->do_autos = 1;
  woden_settings->num_autos = NUM_CROSS;
  woden_settings->num_visis = woden_settings->num_cross + woden_settings->num_autos;

  cropped_sky_models = make_cropped_sky_models(RA0, -0.46606083776035967,
                                                    n_points, n_gauss, n_shapes,
                                                    num_sources);

  printf("We have this many visis %d %d %d\n",woden_settings->num_visis,woden_settings->num_autos,woden_settings->num_cross );
  visibility_set = test_calculate_visibilities(cropped_sky_models,
                                          beam_settings, woden_settings, RA0, -0.46606083776035967,
                                          beam_settings->beamtype);
  test_comp_phase_centre_nobeam(num_comps, visibility_set, woden_settings);

  free_visi_set_inputs(visibility_set);
  free_visi_set_outputs(visibility_set);

}

//Test with a single SOURCE, single COMPONENT
void test_calculate_visibilities_NoBeam_OneSource_SinglePoint(void) {
  int n_points = 1;
  int n_gauss = 0;
  int n_shapes = 0;
  int num_sources = 1;
  test_calculate_visibilities_NoBeam(n_points, n_gauss, n_shapes, num_sources);
}

void test_calculate_visibilities_NoBeam_OneSource_SingleGauss(void) {
  int n_points = 0;
  int n_gauss = 1;
  int n_shapes = 0;
  int num_sources = 1;
  test_calculate_visibilities_NoBeam(n_points, n_gauss, n_shapes, num_sources);

}

void test_calculate_visibilities_NoBeam_OneSource_SingleShape(void) {
  int n_points = 0;
  int n_gauss = 0;
  int n_shapes = 1;
  int num_sources = 1;
  test_calculate_visibilities_NoBeam(n_points, n_gauss, n_shapes, num_sources);
}

void test_calculate_visibilities_NoBeam_OneSource_SingleAll(void) {
  int n_points = 1;
  int n_gauss = 1;
  int n_shapes = 1;
  int num_sources = 1;
  test_calculate_visibilities_NoBeam(n_points, n_gauss, n_shapes, num_sources);
}


//Test with a three SOURCEs, single COMPONENT
void test_calculate_visibilities_NoBeam_ThreeSource_SinglePoint(void) {
  int n_points = 1;
  int n_gauss = 0;
  int n_shapes = 0;
  int num_sources = 3;
  test_calculate_visibilities_NoBeam(n_points, n_gauss, n_shapes, num_sources);

}

void test_calculate_visibilities_NoBeam_ThreeSource_SingleGauss(void) {
  int n_points = 0;
  int n_gauss = 1;
  int n_shapes = 0;
  int num_sources = 3;
  test_calculate_visibilities_NoBeam(n_points, n_gauss, n_shapes, num_sources);
}

void test_calculate_visibilities_NoBeam_ThreeSource_SingleShape(void) {
  int n_points = 0;
  int n_gauss = 0;
  int n_shapes = 1;
  int num_sources = 3;
  test_calculate_visibilities_NoBeam(n_points, n_gauss, n_shapes, num_sources);
}

void test_calculate_visibilities_NoBeam_ThreeSource_SingleAll(void) {
  int n_points = 1;
  int n_gauss = 1;
  int n_shapes = 1;
  int num_sources = 3;
  test_calculate_visibilities_NoBeam(n_points, n_gauss, n_shapes, num_sources);
}

//Test with three SOURCEs, three COPMONENTs
void test_calculate_visibilities_NoBeam_ThreeSource_FivePoint(void) {
  int n_points = 5;
  int n_gauss = 0;
  int n_shapes = 0;
  int num_sources = 3;
  test_calculate_visibilities_NoBeam(n_points, n_gauss, n_shapes, num_sources);

}

void test_calculate_visibilities_NoBeam_ThreeSource_FiveGauss(void) {
  int n_points = 0;
  int n_gauss = 5;
  int n_shapes = 0;
  int num_sources = 3;
  test_calculate_visibilities_NoBeam(n_points, n_gauss, n_shapes, num_sources);
}

void test_calculate_visibilities_NoBeam_ThreeSource_FiveShape(void) {
  int n_points = 0;
  int n_gauss = 0;
  int n_shapes = 5;
  int num_sources = 3;
  test_calculate_visibilities_NoBeam(n_points, n_gauss, n_shapes, num_sources);
}

void test_calculate_visibilities_NoBeam_ThreeSource_FiveAll(void) {
  int n_points = 5;
  int n_gauss = 5;
  int n_shapes = 5;
  int num_sources = 3;
  test_calculate_visibilities_NoBeam(n_points, n_gauss, n_shapes, num_sources);
}


//Run the test with unity
int main(void)
{
    UNITY_BEGIN();
    // //Test with a single SOURCE, single COMPONENT
    RUN_TEST(test_calculate_visibilities_NoBeam_OneSource_SinglePoint);
    RUN_TEST(test_calculate_visibilities_NoBeam_OneSource_SingleGauss);
    RUN_TEST(test_calculate_visibilities_NoBeam_OneSource_SingleShape);
    RUN_TEST(test_calculate_visibilities_NoBeam_OneSource_SingleAll);

    // //Test with three SOURCEs, single COPMONENT
    RUN_TEST(test_calculate_visibilities_NoBeam_ThreeSource_SinglePoint);
    RUN_TEST(test_calculate_visibilities_NoBeam_ThreeSource_SingleGauss);
    RUN_TEST(test_calculate_visibilities_NoBeam_ThreeSource_SingleShape);
    RUN_TEST(test_calculate_visibilities_NoBeam_ThreeSource_SingleAll);

    // //Test with three SOURCEs, five COPMONENTs
    RUN_TEST(test_calculate_visibilities_NoBeam_ThreeSource_FivePoint);
    RUN_TEST(test_calculate_visibilities_NoBeam_ThreeSource_FiveGauss);
    RUN_TEST(test_calculate_visibilities_NoBeam_ThreeSource_FiveShape);
    RUN_TEST(test_calculate_visibilities_NoBeam_ThreeSource_FiveAll);

    return UNITY_END();
}
