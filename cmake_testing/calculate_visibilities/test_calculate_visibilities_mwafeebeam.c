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

//External CUDA code we're linking in
extern void calculate_visibilities(array_layout_t *array_layout,
  source_catalogue_t *cropped_sky_models, beam_settings_t *beam_settings,
  woden_settings_t *woden_settings, visibility_set_t *visibility_set,
  user_precision_t *sbf);

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

#define UNITY_INCLUDE_FLOAT

void test_comp_phase_centre_allgains(visibility_set_t *visibility_set,
                                     user_precision_t gain1xx_re, user_precision_t gain1xx_im,
                                     user_precision_t gain1xy_re, user_precision_t gain1xy_im,
                                     user_precision_t gain1yx_re, user_precision_t gain1yx_im,
                                     user_precision_t gain1yy_re, user_precision_t gain1yy_im,
                                     user_precision_t gain2xx_re, user_precision_t gain2xx_im,
                                     user_precision_t gain2xy_re, user_precision_t gain2xy_im,
                                     user_precision_t gain2yx_re, user_precision_t gain2yx_im,
                                     user_precision_t gain2yy_re, user_precision_t gain2yy_im) {

  //Small user_precision_t errors in the measurement equation mean that when at phase
  //centre, altough imaginary should be zero, you get fluctuations dependent
  //on u,v,w. So need a different tolerance here
  //The accuracy we get depends on whether we are float or double
  #ifdef DOUBLE_PRECISION
    user_precision_t main_pol_im_atol = 1e-10;
    user_precision_t rtol = 1e-10;
    user_precision_t atol = 1e-10;
  #else
    user_precision_t main_pol_im_atol = 1e-5;
    user_precision_t rtol = 1e-6;
    user_precision_t atol = 1e-6;
  #endif

  for (size_t visi = 0; visi < NUM_VISI; visi++) {
    if (visi < NUM_VISI / 2) {
      TEST_ASSERT_FLOAT_WITHIN(gain1xx_re*rtol + atol, gain1xx_re,
                                        visibility_set->sum_visi_XX_real[visi]);
      TEST_ASSERT_FLOAT_WITHIN(main_pol_im_atol + atol, gain1xx_im,
                                        visibility_set->sum_visi_XX_imag[visi]);
      TEST_ASSERT_FLOAT_WITHIN(gain1xy_re*rtol + atol, gain1xy_re,
                                        visibility_set->sum_visi_XY_real[visi]);
      TEST_ASSERT_FLOAT_WITHIN(gain1xy_im*rtol + atol, gain1xy_im,
                                        visibility_set->sum_visi_XY_imag[visi]);
      TEST_ASSERT_FLOAT_WITHIN(gain1yx_re*rtol + atol, gain1yx_re,
                                        visibility_set->sum_visi_YX_real[visi]);
      TEST_ASSERT_FLOAT_WITHIN(gain1yx_im*rtol + atol, gain1yx_im,
                                        visibility_set->sum_visi_YX_imag[visi]);
      TEST_ASSERT_FLOAT_WITHIN(gain1yy_re*rtol + atol, gain1yy_re,
                                        visibility_set->sum_visi_YY_real[visi]);
      TEST_ASSERT_FLOAT_WITHIN(main_pol_im_atol + atol, gain1yy_im,
                                        visibility_set->sum_visi_YY_imag[visi]);
    }
    else {
      TEST_ASSERT_FLOAT_WITHIN(gain2xx_re*rtol + atol, gain2xx_re,
                                        visibility_set->sum_visi_XX_real[visi]);
      TEST_ASSERT_FLOAT_WITHIN(main_pol_im_atol + atol, gain2xx_im,
                                        visibility_set->sum_visi_XX_imag[visi]);
      TEST_ASSERT_FLOAT_WITHIN(gain2xy_re*rtol + atol, gain2xy_re,
                                        visibility_set->sum_visi_XY_real[visi]);
      TEST_ASSERT_FLOAT_WITHIN(gain2xy_im*rtol + atol, gain2xy_im,
                                        visibility_set->sum_visi_XY_imag[visi]);
      TEST_ASSERT_FLOAT_WITHIN(gain2yx_re*rtol + atol, gain2yx_re,
                                        visibility_set->sum_visi_YX_real[visi]);
      TEST_ASSERT_FLOAT_WITHIN(gain2yx_im*rtol + atol, gain2yx_im,
                                        visibility_set->sum_visi_YX_imag[visi]);
      TEST_ASSERT_FLOAT_WITHIN(gain2yy_re*rtol + atol, gain2yy_re,
                                        visibility_set->sum_visi_YY_real[visi]);
      TEST_ASSERT_FLOAT_WITHIN(main_pol_im_atol + atol, gain2yy_im,
                                        visibility_set->sum_visi_YY_imag[visi]);
    }
  }
}

void test_calculate_visibilities_MWAFEEBeam(int n_points, int n_gauss, int n_shapes,
                                           int num_sources) {

  source_catalogue_t *cropped_sky_models = make_cropped_sky_models(RA0, MWA_LAT_RAD,
                                                    n_points, n_gauss, n_shapes,
                                                    num_sources);

  woden_settings_t *woden_settings = make_woden_settings(RA0, MWA_LAT_RAD);
  woden_settings->beamtype = FEE_BEAM;

  beam_settings_t *beam_settings = malloc(sizeof(beam_settings_t));
  beam_settings->beamtype = FEE_BEAM;

  visibility_set_t *visibility_set = test_calculate_visibilities(cropped_sky_models,
                                          beam_settings, woden_settings, RA0, MWA_LAT_RAD,
                                          beam_settings->beamtype);

  printf("WODEN settings %d %d %d\n",woden_settings->num_time_steps,
                                     woden_settings->num_freqs,
                                     woden_settings->num_baselines );
  // for (size_t visi = 0; visi < NUM_VISI; visi++) {
  //   printf("%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n",
  //           visibility_set->sum_visi_XX_real[visi],
  //           visibility_set->sum_visi_XX_imag[visi],
  //           visibility_set->sum_visi_XY_real[visi],
  //           visibility_set->sum_visi_XY_imag[visi],
  //           visibility_set->sum_visi_YX_real[visi],
  //           visibility_set->sum_visi_YX_imag[visi],
  //           visibility_set->sum_visi_YY_real[visi],
  //           visibility_set->sum_visi_YY_imag[visi]);
  // }

  int multiplier = (n_points + n_gauss + n_shapes)*num_sources;

  user_precision_t gain1xx_re = 1.0000000388 * multiplier;
  user_precision_t gain1xx_im = -0.0000000000 * multiplier;
  user_precision_t gain1xy_re = -0.0003180012 * multiplier;
  user_precision_t gain1xy_im = -0.0000022794 * multiplier;
  user_precision_t gain1yx_re = -0.0003180012 * multiplier;
  user_precision_t gain1yx_im = 0.0000022794 * multiplier;
  user_precision_t gain1yy_re = 1.0000000416 * multiplier;
  user_precision_t gain1yy_im = -0.0000000000 * multiplier;

  //The different precisions gives different answers, so switch expectations
  //depending on what precision we are using
  #ifdef DOUBLE_PRECISION
  printf("WODEN is using DOUBLE precision\n");
    user_precision_t gain2xx_re = 0.0069048407 * multiplier;
    user_precision_t gain2xx_im = -0.0000000000 * multiplier;
    user_precision_t gain2xy_re = -0.0004663870 * multiplier;
    user_precision_t gain2xy_im = 0.0000452961 * multiplier;
    user_precision_t gain2yx_re = -0.0004663870 * multiplier;
    user_precision_t gain2yx_im = -0.0000452961 * multiplier;
    user_precision_t gain2yy_re = 0.0038896733 * multiplier;
    user_precision_t gain2yy_im = -0.0000000000 * multiplier;
  #else
    printf("WODEN is using FLOAT precision\n");
    user_precision_t gain2xx_re = 0.0065037487 * multiplier;
    user_precision_t gain2xx_im = -0.0000000000 * multiplier;
    user_precision_t gain2xy_re = -0.0004392957 * multiplier;
    user_precision_t gain2xy_im = 0.0000426647 * multiplier;
    user_precision_t gain2yx_re = -0.0004392957 * multiplier;
    user_precision_t gain2yx_im = -0.0000426647 * multiplier;
    user_precision_t gain2yy_re = 0.0036637327 * multiplier;
    user_precision_t gain2yy_im = -0.0000000000 * multiplier;
  #endif

  test_comp_phase_centre_allgains(visibility_set,
                                  gain1xx_re, gain1xx_im,
                                  gain1xy_re, gain1xy_im,
                                  gain1yx_re, gain1yx_im,
                                  gain1yy_re, gain1yy_im,
                                  gain2xx_re, gain2xx_im,
                                  gain2xy_re, gain2xy_im,
                                  gain2yx_re, gain2yx_im,
                                  gain2yy_re, gain2yy_im);

  free(beam_settings);
  free_visi_set_inputs(visibility_set);
  free_visi_set_outputs(visibility_set);

  // free(visibility_set);

}

//Test with a single SOURCE, single COMPONENT
void test_calculate_visibilities_MWAFEEBeam_OneSource_SinglePoint(void) {
  int n_points = 1;
  int n_gauss = 0;
  int n_shapes = 0;
  int num_sources = 1;
  test_calculate_visibilities_MWAFEEBeam(n_points, n_gauss, n_shapes, num_sources);
}

void test_calculate_visibilities_MWAFEEBeam_OneSource_SingleGauss(void) {
  int n_points = 0;
  int n_gauss = 1;
  int n_shapes = 0;
  int num_sources = 1;
  test_calculate_visibilities_MWAFEEBeam(n_points, n_gauss, n_shapes, num_sources);

}

void test_calculate_visibilities_MWAFEEBeam_OneSource_SingleShape(void) {
  int n_points = 0;
  int n_gauss = 0;
  int n_shapes = 1;
  int num_sources = 1;
  test_calculate_visibilities_MWAFEEBeam(n_points, n_gauss, n_shapes, num_sources);
}

void test_calculate_visibilities_MWAFEEBeam_OneSource_SingleAll(void) {
  int n_points = 1;
  int n_gauss = 1;
  int n_shapes = 1;
  int num_sources = 1;
  test_calculate_visibilities_MWAFEEBeam(n_points, n_gauss, n_shapes, num_sources);
}


//Test with a three SOURCEs, single COMPONENT
void test_calculate_visibilities_MWAFEEBeam_ThreeSource_SinglePoint(void) {
  int n_points = 1;
  int n_gauss = 0;
  int n_shapes = 0;
  int num_sources = 3;
  test_calculate_visibilities_MWAFEEBeam(n_points, n_gauss, n_shapes, num_sources);

}

void test_calculate_visibilities_MWAFEEBeam_ThreeSource_SingleGauss(void) {
  int n_points = 0;
  int n_gauss = 1;
  int n_shapes = 0;
  int num_sources = 3;
  test_calculate_visibilities_MWAFEEBeam(n_points, n_gauss, n_shapes, num_sources);
}

void test_calculate_visibilities_MWAFEEBeam_ThreeSource_SingleShape(void) {
  int n_points = 0;
  int n_gauss = 0;
  int n_shapes = 1;
  int num_sources = 3;
  test_calculate_visibilities_MWAFEEBeam(n_points, n_gauss, n_shapes, num_sources);
}

void test_calculate_visibilities_MWAFEEBeam_ThreeSource_SingleAll(void) {
  int n_points = 1;
  int n_gauss = 1;
  int n_shapes = 1;
  int num_sources = 3;
  test_calculate_visibilities_MWAFEEBeam(n_points, n_gauss, n_shapes, num_sources);
}


//Test with three SOURCEs, three COPMONENTs
void test_calculate_visibilities_MWAFEEBeam_ThreeSource_ThreePoint(void) {
  int n_points = 3;
  int n_gauss = 0;
  int n_shapes = 0;
  int num_sources = 3;
  test_calculate_visibilities_MWAFEEBeam(n_points, n_gauss, n_shapes, num_sources);

}

void test_calculate_visibilities_MWAFEEBeam_ThreeSource_ThreeMWAFEE(void) {
  int n_points = 0;
  int n_gauss = 3;
  int n_shapes = 0;
  int num_sources = 3;
  test_calculate_visibilities_MWAFEEBeam(n_points, n_gauss, n_shapes, num_sources);
}

void test_calculate_visibilities_MWAFEEBeam_ThreeSource_ThreeShape(void) {
  int n_points = 0;
  int n_gauss = 0;
  int n_shapes = 3;
  int num_sources = 3;
  test_calculate_visibilities_MWAFEEBeam(n_points, n_gauss, n_shapes, num_sources);
}

void test_calculate_visibilities_MWAFEEBeam_ThreeSource_ThreeAll(void) {
  int n_points = 3;
  int n_gauss = 3;
  int n_shapes = 3;
  int num_sources = 3;
  test_calculate_visibilities_MWAFEEBeam(n_points, n_gauss, n_shapes, num_sources);
}



// Run the test with unity
int main(void)
{
    UNITY_BEGIN();
    //Test with a single SOURCE, single COMPONENT
    RUN_TEST(test_calculate_visibilities_MWAFEEBeam_OneSource_SinglePoint);
    // RUN_TEST(test_calculate_visibilities_MWAFEEBeam_OneSource_SingleGauss);
    // RUN_TEST(test_calculate_visibilities_MWAFEEBeam_OneSource_SingleShape);
    // RUN_TEST(test_calculate_visibilities_MWAFEEBeam_OneSource_SingleAll);

    // //Test with three SOURCEs, single COPMONENT
    // RUN_TEST(test_calculate_visibilities_MWAFEEBeam_ThreeSource_SinglePoint);
    // RUN_TEST(test_calculate_visibilities_MWAFEEBeam_ThreeSource_SingleGauss);
    // RUN_TEST(test_calculate_visibilities_MWAFEEBeam_ThreeSource_SingleShape);
    // RUN_TEST(test_calculate_visibilities_MWAFEEBeam_ThreeSource_SingleAll);
    //
    // //Test with three SOURCEs, three COPMONENTs
    // RUN_TEST(test_calculate_visibilities_MWAFEEBeam_ThreeSource_ThreePoint);
    // RUN_TEST(test_calculate_visibilities_MWAFEEBeam_ThreeSource_ThreeMWAFEE);
    // RUN_TEST(test_calculate_visibilities_MWAFEEBeam_ThreeSource_ThreeShape);
    // RUN_TEST(test_calculate_visibilities_MWAFEEBeam_ThreeSource_ThreeAll);

    return UNITY_END();
}
