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

void test_comp_phase_centre_allgains(visibility_set_t *visibility_set,
                                     double gain1xx_re, double gain1xx_im,
                                     double gain1xy_re, double gain1xy_im,
                                     double gain1yx_re, double gain1yx_im,
                                     double gain1yy_re, double gain1yy_im,
                                     double gain2xx_re, double gain2xx_im,
                                     double gain2xy_re, double gain2xy_im,
                                     double gain2yx_re, double gain2yx_im,
                                     double gain2yy_re, double gain2yy_im) {

  //Small user_precision_t errors in the measurement equation mean that when at phase
  //centre, altough imaginary should be zero, you get fluctuations dependent
  //on u,v,w. So need a different tolerance here
  //The accuracy we get depends on whether we are float or double
  #ifdef DOUBLE_PRECISION
    double TOL = 1e-7;
  #else
    double TOL = 4e-3;
  #endif

  for (int visi = 0; visi < NUM_VISI; visi++) {
    if (visi < NUM_VISI / 2) {
      TEST_ASSERT_DOUBLE_WITHIN(TOL, gain1xx_re,
                                        visibility_set->sum_visi_XX_real[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(TOL, gain1xx_im,
                                        visibility_set->sum_visi_XX_imag[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(TOL, gain1xy_re,
                                        visibility_set->sum_visi_XY_real[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(TOL, gain1xy_im,
                                        visibility_set->sum_visi_XY_imag[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(TOL, gain1yx_re,
                                        visibility_set->sum_visi_YX_real[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(TOL, gain1yx_im,
                                        visibility_set->sum_visi_YX_imag[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(TOL, gain1yy_re,
                                        visibility_set->sum_visi_YY_real[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(TOL, gain1yy_im,
                                        visibility_set->sum_visi_YY_imag[visi]);
    }
    else {
      TEST_ASSERT_DOUBLE_WITHIN(TOL, gain2xx_re,
                                        visibility_set->sum_visi_XX_real[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(TOL, gain2xx_im,
                                        visibility_set->sum_visi_XX_imag[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(TOL, gain2xy_re,
                                        visibility_set->sum_visi_XY_real[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(TOL, gain2xy_im,
                                        visibility_set->sum_visi_XY_imag[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(TOL, gain2yx_re,
                                        visibility_set->sum_visi_YX_real[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(TOL, gain2yx_im,
                                        visibility_set->sum_visi_YX_imag[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(TOL, gain2yy_re,
                                        visibility_set->sum_visi_YY_real[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(TOL, gain2yy_im,
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

  double multiplier = (n_points + n_gauss + n_shapes)*num_sources*STOKESI;


  //These values are taken from the double precision version of the MWA FEE
  //beam code
  double gain1xx_re = 1.0000000388296133 * multiplier;
  double gain1xx_im = 0.0 * multiplier;
  double gain1xy_re = -0.0003180012451422 * multiplier;
  double gain1xy_im = -0.0000022794060068 * multiplier;
  double gain1yx_re = -0.0003180012451422 * multiplier;
  double gain1yx_im = 0.0000022794060068 * multiplier;
  double gain1yy_re = 1.0000000415988461 * multiplier;
  double gain1yy_im = 0.0 * multiplier;

  double gain2xx_re = 0.0069048406570405 * multiplier;
  double gain2xx_im = 0.0 * multiplier;
  double gain2xy_re = -0.0004663870439568 * multiplier;
  double gain2xy_im = 0.0000452960552731 * multiplier;
  double gain2yx_re = -0.0004663870439568 * multiplier;
  double gain2yx_im = -0.0000452960552731 * multiplier;
  double gain2yy_re = 0.0038896732780955* multiplier;
  double gain2yy_im = 0.0 * multiplier;

  test_comp_phase_centre_allgains(visibility_set,
                                  gain1xx_re, gain1xx_im,
                                  gain1xy_re, gain1xy_im,
                                  gain1yx_re, gain1yx_im,
                                  gain1yy_re, gain1yy_im,
                                  gain2xx_re, gain2xx_im,
                                  gain2xy_re, gain2xy_im,
                                  gain2yx_re, gain2yx_im,
                                  gain2yy_re, gain2yy_im);

  RTS_freeHDFBeam(beam_settings->FEE_beam);
  RTS_freeHDFBeam(beam_settings->FEE_beam_zenith);

  free(beam_settings);
  free_visi_set_inputs(visibility_set);
  free_visi_set_outputs(visibility_set);
  free(woden_settings);

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
