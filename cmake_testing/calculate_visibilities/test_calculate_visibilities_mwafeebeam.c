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
  float *sbf);

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

#define UNITY_INCLUDE_FLOAT

void test_comp_phase_centre_allgains(visibility_set_t *visibility_set,
                                     float gain1xx_re, float gain1xx_im,
                                     float gain1xy_re, float gain1xy_im,
                                     float gain1yx_re, float gain1yx_im,
                                     float gain1yy_re, float gain1yy_im,
                                     float gain2xx_re, float gain2xx_im,
                                     float gain2xy_re, float gain2xy_im,
                                     float gain2yx_re, float gain2yx_im,
                                     float gain2yy_re, float gain2yy_im) {

  float rtol = 1e-4;
  float atol = 1e-5;

  //Small float errors in the measurement equation mean that when at phase
  //centre, altough imaginary should be zero, you get fluctuations dependent
  //on u,v,w. So need a different tolerance here
  //This is trade off we make for speed of float over using double
  float main_pol_im_atol = 1e-3;

  for (size_t visi = 0; visi < NUM_VISI; visi++) {
    if (visi < NUM_VISI / 2) {
      TEST_ASSERT_FLOAT_WITHIN(gain1xx_re*rtol + atol, gain1xx_re, visibility_set->sum_visi_XX_real[visi]);
      TEST_ASSERT_FLOAT_WITHIN(main_pol_im_atol + atol, gain1xx_im, visibility_set->sum_visi_XX_imag[visi]);
      TEST_ASSERT_FLOAT_WITHIN(gain1xy_re*rtol + atol, gain1xy_re, visibility_set->sum_visi_XY_real[visi]);
      TEST_ASSERT_FLOAT_WITHIN(gain1xy_im*rtol + atol, gain1xy_im, visibility_set->sum_visi_XY_imag[visi]);
      TEST_ASSERT_FLOAT_WITHIN(gain1yx_re*rtol + atol, gain1yx_re, visibility_set->sum_visi_YX_real[visi]);
      TEST_ASSERT_FLOAT_WITHIN(gain1yx_im*rtol + atol, gain1yx_im, visibility_set->sum_visi_YX_imag[visi]);
      TEST_ASSERT_FLOAT_WITHIN(gain1yy_re*rtol + atol, gain1yy_re, visibility_set->sum_visi_YY_real[visi]);
      TEST_ASSERT_FLOAT_WITHIN(main_pol_im_atol + atol, gain1yy_im, visibility_set->sum_visi_YY_imag[visi]);
    }
    // } else {
    //   TEST_ASSERT_FLOAT_WITHIN(gainxx*1e-4, gainxx, visibility_set->sum_visi_XX_real[visi]);
    //   TEST_ASSERT_FLOAT_WITHIN(gainxx*1e-4, 0.0, visibility_set->sum_visi_XX_imag[visi]);
    //   TEST_ASSERT_FLOAT_WITHIN(gainxx*1e-4, 0.0, visibility_set->sum_visi_XY_real[visi]);
    //   TEST_ASSERT_FLOAT_WITHIN(gainxx*1e-4, 0.0, visibility_set->sum_visi_XY_imag[visi]);
    //   TEST_ASSERT_FLOAT_WITHIN(gainxx*1e-4, 0.0, visibility_set->sum_visi_YX_real[visi]);
    //   TEST_ASSERT_FLOAT_WITHIN(gainxx*1e-4, 0.0, visibility_set->sum_visi_YX_imag[visi]);
    //   TEST_ASSERT_FLOAT_WITHIN(gainyy*1e-4, gainyy, visibility_set->sum_visi_YY_real[visi]);
    //   TEST_ASSERT_FLOAT_WITHIN(gainyy*1e-4, 0.0, visibility_set->sum_visi_YY_imag[visi]);
    // }
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

  float gain1xx_re = 0.99999422 * (n_points + n_gauss + n_shapes)*num_sources;
  float gain1xx_im = -0.00000337 * (n_points + n_gauss + n_shapes)*num_sources;
  float gain1xy_re = -0.00031884 * (n_points + n_gauss + n_shapes)*num_sources;
  float gain1xy_im = -0.00000219 * (n_points + n_gauss + n_shapes)*num_sources;
  float gain1yx_re = -0.00031884 * (n_points + n_gauss + n_shapes)*num_sources;
  float gain1yx_im = 0.00000219 * (n_points + n_gauss + n_shapes)*num_sources;
  float gain1yy_re = 0.99999422 * (n_points + n_gauss + n_shapes)*num_sources;
  float gain1yy_im = -0.00000337 * (n_points + n_gauss + n_shapes)*num_sources;

  float gain2xx_re = 0.00650369 * (n_points + n_gauss + n_shapes)*num_sources;
  float gain2xx_im = 0.00000000 * (n_points + n_gauss + n_shapes)*num_sources;
  float gain2xy_re = -0.00043929 * (n_points + n_gauss + n_shapes)*num_sources;
  float gain2xy_im = 0.00004266 * (n_points + n_gauss + n_shapes)*num_sources;
  float gain2yx_re = -0.00043929 * (n_points + n_gauss + n_shapes)*num_sources;
  float gain2yx_im = -0.00004266 * (n_points + n_gauss + n_shapes)*num_sources;
  float gain2yy_re = 0.00366370 * (n_points + n_gauss + n_shapes)*num_sources;
  float gain2yy_im = 0.00000000 * (n_points + n_gauss + n_shapes)*num_sources;

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
    RUN_TEST(test_calculate_visibilities_MWAFEEBeam_OneSource_SingleGauss);
    RUN_TEST(test_calculate_visibilities_MWAFEEBeam_OneSource_SingleShape);
    RUN_TEST(test_calculate_visibilities_MWAFEEBeam_OneSource_SingleAll);

    //Test with three SOURCEs, single COPMONENT
    RUN_TEST(test_calculate_visibilities_MWAFEEBeam_ThreeSource_SinglePoint);
    RUN_TEST(test_calculate_visibilities_MWAFEEBeam_ThreeSource_SingleGauss);
    RUN_TEST(test_calculate_visibilities_MWAFEEBeam_ThreeSource_SingleShape);
    RUN_TEST(test_calculate_visibilities_MWAFEEBeam_ThreeSource_SingleAll);

    //Test with three SOURCEs, three COPMONENTs
    RUN_TEST(test_calculate_visibilities_MWAFEEBeam_ThreeSource_ThreePoint);
    RUN_TEST(test_calculate_visibilities_MWAFEEBeam_ThreeSource_ThreeMWAFEE);
    RUN_TEST(test_calculate_visibilities_MWAFEEBeam_ThreeSource_ThreeShape);
    RUN_TEST(test_calculate_visibilities_MWAFEEBeam_ThreeSource_ThreeAll);

    return UNITY_END();
}
