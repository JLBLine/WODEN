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

void test_calculate_visibilities_GaussBeam(int n_points, int n_gauss, int n_shapes,
                                           int num_sources) {

  source_catalogue_t *cropped_sky_models = make_cropped_sky_models(RA0, -0.46606083776035967,
                                                    n_points, n_gauss, n_shapes,
                                                    num_sources);

  beam_settings_t *beam_settings = malloc(sizeof(beam_settings_t));
  beam_settings->beamtype = GAUSS_BEAM;

  //Make some fiducial settings - point beam slightly away
  double dec_point = -30.0*DD2R;
  beam_settings->gauss_sdec = sin(dec_point);
  beam_settings->gauss_cdec = cos(dec_point);
  beam_settings->gauss_ha = 5.0*DD2R;

  beam_settings->beam_FWHM_rad = 60*DD2R;
  beam_settings->beam_ref_freq = 180e+6;

  woden_settings_t *woden_settings = make_woden_settings(RA0, -0.46606083776035967);
  woden_settings->beamtype = GAUSS_BEAM;

  visibility_set_t *visibility_set = test_calculate_visibilities(cropped_sky_models,
                                          beam_settings, woden_settings, RA0, -0.46606083776035967,
                                          beam_settings->beamtype);

  // user_precision_t gain = (n_points + n_gauss + n_shapes)*num_sources;
  // test_comp_phase_centre(visibility_set, gain);

  int num_comps = (n_points + n_gauss + n_shapes)*num_sources;

  double _Complex gain1 = 0.985035846188 + I*0.0;
  double _Complex gain2 = 0.579610294748 + I*0.0;

  test_comp_phase_centre_twogains(visibility_set, num_comps, gain1, gain1,
                                  gain2, gain2, woden_settings);

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
  test_comp_phase_centre_twogains(visibility_set, num_comps, gain1, gain1,
                                  gain2, gain2, woden_settings);

  free_visi_set_inputs(visibility_set);
  free_visi_set_outputs(visibility_set);



}

//Test with a single SOURCE, single COMPONENT
void test_calculate_visibilities_GaussBeam_OneSource_SinglePoint(void) {
  int n_points = 1;
  int n_gauss = 0;
  int n_shapes = 0;
  int num_sources = 1;
  test_calculate_visibilities_GaussBeam(n_points, n_gauss, n_shapes, num_sources);
}

void test_calculate_visibilities_GaussBeam_OneSource_SingleGauss(void) {
  int n_points = 0;
  int n_gauss = 1;
  int n_shapes = 0;
  int num_sources = 1;
  test_calculate_visibilities_GaussBeam(n_points, n_gauss, n_shapes, num_sources);

}

void test_calculate_visibilities_GaussBeam_OneSource_SingleShape(void) {
  int n_points = 0;
  int n_gauss = 0;
  int n_shapes = 1;
  int num_sources = 1;
  test_calculate_visibilities_GaussBeam(n_points, n_gauss, n_shapes, num_sources);
}

void test_calculate_visibilities_GaussBeam_OneSource_SingleAll(void) {
  int n_points = 1;
  int n_gauss = 1;
  int n_shapes = 1;
  int num_sources = 1;
  test_calculate_visibilities_GaussBeam(n_points, n_gauss, n_shapes, num_sources);
}


//Test with a three SOURCEs, single COMPONENT
void test_calculate_visibilities_GaussBeam_ThreeSource_SinglePoint(void) {
  int n_points = 1;
  int n_gauss = 0;
  int n_shapes = 0;
  int num_sources = 3;
  test_calculate_visibilities_GaussBeam(n_points, n_gauss, n_shapes, num_sources);

}

void test_calculate_visibilities_GaussBeam_ThreeSource_SingleGauss(void) {
  int n_points = 0;
  int n_gauss = 1;
  int n_shapes = 0;
  int num_sources = 3;
  test_calculate_visibilities_GaussBeam(n_points, n_gauss, n_shapes, num_sources);
}

void test_calculate_visibilities_GaussBeam_ThreeSource_SingleShape(void) {
  int n_points = 0;
  int n_gauss = 0;
  int n_shapes = 1;
  int num_sources = 3;
  test_calculate_visibilities_GaussBeam(n_points, n_gauss, n_shapes, num_sources);
}

void test_calculate_visibilities_GaussBeam_ThreeSource_SingleAll(void) {
  int n_points = 1;
  int n_gauss = 1;
  int n_shapes = 1;
  int num_sources = 3;
  test_calculate_visibilities_GaussBeam(n_points, n_gauss, n_shapes, num_sources);
}


//Test with three SOURCEs, three COPMONENTs
void test_calculate_visibilities_GaussBeam_ThreeSource_ThreePoint(void) {
  int n_points = 3;
  int n_gauss = 0;
  int n_shapes = 0;
  int num_sources = 3;
  test_calculate_visibilities_GaussBeam(n_points, n_gauss, n_shapes, num_sources);

}

void test_calculate_visibilities_GaussBeam_ThreeSource_ThreeGauss(void) {
  int n_points = 0;
  int n_gauss = 3;
  int n_shapes = 0;
  int num_sources = 3;
  test_calculate_visibilities_GaussBeam(n_points, n_gauss, n_shapes, num_sources);
}

void test_calculate_visibilities_GaussBeam_ThreeSource_ThreeShape(void) {
  int n_points = 0;
  int n_gauss = 0;
  int n_shapes = 3;
  int num_sources = 3;
  test_calculate_visibilities_GaussBeam(n_points, n_gauss, n_shapes, num_sources);
}

void test_calculate_visibilities_GaussBeam_ThreeSource_ThreeAll(void) {
  int n_points = 3;
  int n_gauss = 3;
  int n_shapes = 3;
  int num_sources = 3;
  test_calculate_visibilities_GaussBeam(n_points, n_gauss, n_shapes, num_sources);
}



// Run the test with unity
int main(void)
{
    UNITY_BEGIN();
    //Test with a single SOURCE, single COMPONENT
    RUN_TEST(test_calculate_visibilities_GaussBeam_OneSource_SinglePoint);
    RUN_TEST(test_calculate_visibilities_GaussBeam_OneSource_SingleGauss);
    RUN_TEST(test_calculate_visibilities_GaussBeam_OneSource_SingleShape);
    RUN_TEST(test_calculate_visibilities_GaussBeam_OneSource_SingleAll);

    //Test with three SOURCEs, single COPMONENT
    RUN_TEST(test_calculate_visibilities_GaussBeam_ThreeSource_SinglePoint);
    RUN_TEST(test_calculate_visibilities_GaussBeam_ThreeSource_SingleGauss);
    RUN_TEST(test_calculate_visibilities_GaussBeam_ThreeSource_SingleShape);
    RUN_TEST(test_calculate_visibilities_GaussBeam_ThreeSource_SingleAll);

    //Test with three SOURCEs, three COPMONENTs
    RUN_TEST(test_calculate_visibilities_GaussBeam_ThreeSource_ThreePoint);
    RUN_TEST(test_calculate_visibilities_GaussBeam_ThreeSource_ThreeGauss);
    RUN_TEST(test_calculate_visibilities_GaussBeam_ThreeSource_ThreeShape);
    RUN_TEST(test_calculate_visibilities_GaussBeam_ThreeSource_ThreeAll);

    return UNITY_END();
}
