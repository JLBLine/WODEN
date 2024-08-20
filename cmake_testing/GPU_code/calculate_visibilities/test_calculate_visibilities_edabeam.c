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

  source_catalogue_t *cropped_sky_models = make_cropped_sky_models(RA0, -0.46606083776035967,
                                                    n_points, n_gauss, n_shapes,
                                                    num_sources);

  beam_settings_t *beam_settings = malloc(sizeof(beam_settings_t));
  beam_settings->beamtype = ANALY_DIPOLE;

  woden_settings_t *woden_settings = make_woden_settings(RA0, -0.46606083776035967);
  woden_settings->beamtype = ANALY_DIPOLE;

  visibility_set_t *visibility_set = test_calculate_visibilities(cropped_sky_models,
                                          beam_settings, woden_settings, RA0, -0.46606083776035967,
                                          beam_settings->beamtype);

  int num_comps = (n_points + n_gauss + n_shapes)*num_sources;

  double _Complex gain1x = 1.0 + I*0.0;
  double _Complex gain1y = 1.0 + I*0.0;
  double _Complex gain2x = 0.792312729527 + I*0.0;
  double _Complex gain2y = 0.618507509916 + I*0.0;

  test_comp_phase_centre_twogains(visibility_set, num_comps, gain1x, gain1y,
                                  gain2x, gain2y, woden_settings);

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
  test_comp_phase_centre_twogains(visibility_set, num_comps, gain1x, gain1y,
                                  gain2x, gain2y, woden_settings);

  free_visi_set_inputs(visibility_set);
  free_visi_set_outputs(visibility_set);

}

//Test with a single SOURCE, single COMPONENT
void test_calculate_visibilities_EDA2Beam_OneSource_SinglePoint(void) {
  int n_points = 1;
  int n_gauss = 0;
  int n_shapes = 0;
  int num_sources = 1;
  test_calculate_visibilities_EDA2Beam(n_points, n_gauss, n_shapes, num_sources);
}

void test_calculate_visibilities_EDA2Beam_OneSource_SingleGauss(void) {
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

void test_calculate_visibilities_EDA2Beam_ThreeSource_SingleGauss(void) {
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
void test_calculate_visibilities_EDA2Beam_ThreeSource_FivePoint(void) {
  int n_points = 5;
  int n_gauss = 0;
  int n_shapes = 0;
  int num_sources = 3;
  test_calculate_visibilities_EDA2Beam(n_points, n_gauss, n_shapes, num_sources);

}

void test_calculate_visibilities_EDA2Beam_ThreeSource_FiveEDA2(void) {
  int n_points = 0;
  int n_gauss = 5;
  int n_shapes = 0;
  int num_sources = 3;
  test_calculate_visibilities_EDA2Beam(n_points, n_gauss, n_shapes, num_sources);
}

void test_calculate_visibilities_EDA2Beam_ThreeSource_FiveShape(void) {
  int n_points = 0;
  int n_gauss = 0;
  int n_shapes = 5;
  int num_sources = 3;
  test_calculate_visibilities_EDA2Beam(n_points, n_gauss, n_shapes, num_sources);
}

void test_calculate_visibilities_EDA2Beam_ThreeSource_FiveAll(void) {
  int n_points = 5;
  int n_gauss = 5;
  int n_shapes = 5;
  int num_sources = 3;
  test_calculate_visibilities_EDA2Beam(n_points, n_gauss, n_shapes, num_sources);
}



// Run the test with unity
int main(void)
{
    UNITY_BEGIN();
    //Test with a single SOURCE, single COMPONENT
    RUN_TEST(test_calculate_visibilities_EDA2Beam_OneSource_SinglePoint);
    RUN_TEST(test_calculate_visibilities_EDA2Beam_OneSource_SingleGauss);
    RUN_TEST(test_calculate_visibilities_EDA2Beam_OneSource_SingleShape);
    RUN_TEST(test_calculate_visibilities_EDA2Beam_OneSource_SingleAll);

    //Test with three SOURCEs, single COPMONENT
    RUN_TEST(test_calculate_visibilities_EDA2Beam_ThreeSource_SinglePoint);
    RUN_TEST(test_calculate_visibilities_EDA2Beam_ThreeSource_SingleGauss);
    RUN_TEST(test_calculate_visibilities_EDA2Beam_ThreeSource_SingleShape);
    RUN_TEST(test_calculate_visibilities_EDA2Beam_ThreeSource_SingleAll);

    //Test with three SOURCEs, three COPMONENTs
    RUN_TEST(test_calculate_visibilities_EDA2Beam_ThreeSource_FivePoint);
    RUN_TEST(test_calculate_visibilities_EDA2Beam_ThreeSource_FiveEDA2);
    RUN_TEST(test_calculate_visibilities_EDA2Beam_ThreeSource_FiveShape);
    RUN_TEST(test_calculate_visibilities_EDA2Beam_ThreeSource_FiveAll);

    return UNITY_END();
}
