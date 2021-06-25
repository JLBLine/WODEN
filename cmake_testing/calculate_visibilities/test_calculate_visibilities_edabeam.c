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

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

#define UNITY_INCLUDE_FLOAT

void test_calculate_visibilities_EDA2Beam(int n_points, int n_gauss, int n_shapes,
                                           int num_sources) {

  source_catalogue_t *cropped_sky_models = make_cropped_sky_models(RA0, MWA_LAT_RAD,
                                                    n_points, n_gauss, n_shapes,
                                                    num_sources);

  beam_settings_t *beam_settings = malloc(sizeof(beam_settings_t));
  beam_settings->beamtype = ANALY_DIPOLE;

  woden_settings_t *woden_settings = make_woden_settings(RA0, MWA_LAT_RAD);
  woden_settings->beamtype = ANALY_DIPOLE;

  visibility_set_t *visibility_set = test_calculate_visibilities(cropped_sky_models,
                                          beam_settings, woden_settings, RA0, MWA_LAT_RAD,
                                          beam_settings->beamtype);

  float gain1xx = (n_points + n_gauss + n_shapes)*num_sources;
  float gain1yy = (n_points + n_gauss + n_shapes)*num_sources;

  float gain2xx = 0.6277595 * (n_points + n_gauss + n_shapes)*num_sources;
  float gain2yy = 0.3825515 * (n_points + n_gauss + n_shapes)*num_sources;
  //
  test_comp_phase_centre_twogains(visibility_set, gain1xx, gain1yy,
                                  gain2xx, gain2yy);

  // for (size_t visi = 0; visi < NUM_VISI; visi++) {
  //   printf("%.7f %.1f %.1f %.1f %.1f %.1f %.7f %.1f\n",
  //           visibility_set->sum_visi_XX_real[visi],
  //           visibility_set->sum_visi_XX_imag[visi],
  //           visibility_set->sum_visi_XY_real[visi],
  //           visibility_set->sum_visi_XY_imag[visi],
  //           visibility_set->sum_visi_YX_real[visi],
  //           visibility_set->sum_visi_YX_imag[visi],
  //           visibility_set->sum_visi_YY_real[visi],
  //           visibility_set->sum_visi_YY_imag[visi]);
  // }

}

//Test with a single SOURCE, single COMPONENT
void test_calculate_visibilities_EDA2Beam_OneSource_SinglePoint(void) {
  int n_points = 1;
  int n_gauss = 0;
  int n_shapes = 0;
  int num_sources = 1;
  test_calculate_visibilities_EDA2Beam(n_points, n_gauss, n_shapes, num_sources);
}

void test_calculate_visibilities_EDA2Beam_OneSource_SingleEDA2(void) {
  int n_points = 0;
  int n_gauss = 1;
  int n_shapes = 0;
  int num_sources = 1;
  test_calculate_visibilities_EDA2Beam(n_points, n_gauss, n_shapes, num_sources);

}

void test_calculate_visibilities_EDA2Beam_OneSource_SingleShape(void) {
  int n_points = 0;
  int n_gauss = 0;
  int n_shapes = 1;
  int num_sources = 1;
  test_calculate_visibilities_EDA2Beam(n_points, n_gauss, n_shapes, num_sources);
}

void test_calculate_visibilities_EDA2Beam_OneSource_SingleAll(void) {
  int n_points = 1;
  int n_gauss = 1;
  int n_shapes = 1;
  int num_sources = 1;
  test_calculate_visibilities_EDA2Beam(n_points, n_gauss, n_shapes, num_sources);
}


//Test with a three SOURCEs, single COMPONENT
void test_calculate_visibilities_EDA2Beam_ThreeSource_SinglePoint(void) {
  int n_points = 1;
  int n_gauss = 0;
  int n_shapes = 0;
  int num_sources = 3;
  test_calculate_visibilities_EDA2Beam(n_points, n_gauss, n_shapes, num_sources);

}

void test_calculate_visibilities_EDA2Beam_ThreeSource_SingleEDA2(void) {
  int n_points = 0;
  int n_gauss = 1;
  int n_shapes = 0;
  int num_sources = 3;
  test_calculate_visibilities_EDA2Beam(n_points, n_gauss, n_shapes, num_sources);
}

void test_calculate_visibilities_EDA2Beam_ThreeSource_SingleShape(void) {
  int n_points = 0;
  int n_gauss = 0;
  int n_shapes = 1;
  int num_sources = 3;
  test_calculate_visibilities_EDA2Beam(n_points, n_gauss, n_shapes, num_sources);
}

void test_calculate_visibilities_EDA2Beam_ThreeSource_SingleAll(void) {
  int n_points = 1;
  int n_gauss = 1;
  int n_shapes = 1;
  int num_sources = 3;
  test_calculate_visibilities_EDA2Beam(n_points, n_gauss, n_shapes, num_sources);
}


//Test with three SOURCEs, three COPMONENTs
void test_calculate_visibilities_EDA2Beam_ThreeSource_ThreePoint(void) {
  int n_points = 3;
  int n_gauss = 0;
  int n_shapes = 0;
  int num_sources = 3;
  test_calculate_visibilities_EDA2Beam(n_points, n_gauss, n_shapes, num_sources);

}

void test_calculate_visibilities_EDA2Beam_ThreeSource_ThreeEDA2(void) {
  int n_points = 0;
  int n_gauss = 3;
  int n_shapes = 0;
  int num_sources = 3;
  test_calculate_visibilities_EDA2Beam(n_points, n_gauss, n_shapes, num_sources);
}

void test_calculate_visibilities_EDA2Beam_ThreeSource_ThreeShape(void) {
  int n_points = 0;
  int n_gauss = 0;
  int n_shapes = 3;
  int num_sources = 3;
  test_calculate_visibilities_EDA2Beam(n_points, n_gauss, n_shapes, num_sources);
}

void test_calculate_visibilities_EDA2Beam_ThreeSource_ThreeAll(void) {
  int n_points = 3;
  int n_gauss = 3;
  int n_shapes = 3;
  int num_sources = 3;
  test_calculate_visibilities_EDA2Beam(n_points, n_gauss, n_shapes, num_sources);
}



// Run the test with unity
int main(void)
{
    UNITY_BEGIN();
    //Test with a single SOURCE, single COMPONENT
    RUN_TEST(test_calculate_visibilities_EDA2Beam_OneSource_SinglePoint);
    RUN_TEST(test_calculate_visibilities_EDA2Beam_OneSource_SingleEDA2);
    RUN_TEST(test_calculate_visibilities_EDA2Beam_OneSource_SingleShape);
    RUN_TEST(test_calculate_visibilities_EDA2Beam_OneSource_SingleAll);

    //Test with three SOURCEs, single COPMONENT
    RUN_TEST(test_calculate_visibilities_EDA2Beam_ThreeSource_SinglePoint);
    RUN_TEST(test_calculate_visibilities_EDA2Beam_ThreeSource_SingleEDA2);
    RUN_TEST(test_calculate_visibilities_EDA2Beam_ThreeSource_SingleShape);
    RUN_TEST(test_calculate_visibilities_EDA2Beam_ThreeSource_SingleAll);

    //Test with three SOURCEs, three COPMONENTs
    RUN_TEST(test_calculate_visibilities_EDA2Beam_ThreeSource_ThreePoint);
    RUN_TEST(test_calculate_visibilities_EDA2Beam_ThreeSource_ThreeEDA2);
    RUN_TEST(test_calculate_visibilities_EDA2Beam_ThreeSource_ThreeShape);
    RUN_TEST(test_calculate_visibilities_EDA2Beam_ThreeSource_ThreeAll);

    return UNITY_END();
}
