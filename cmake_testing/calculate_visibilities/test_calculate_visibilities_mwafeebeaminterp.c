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

void test_calculate_visibilities_MWAFEEBeamInterp(int n_points, int n_gauss, int n_shapes,
                                           int num_sources) {

  source_catalogue_t *cropped_sky_models = make_cropped_sky_models(RA0, -0.46606083776035967,
                                                    n_points, n_gauss, n_shapes,
                                                    num_sources);

  woden_settings_t *woden_settings = make_woden_settings(RA0, -0.46606083776035967);
  woden_settings->beamtype = FEE_BEAM_INTERP;

  for (int i = 0; i < 16; i++) {
    woden_settings->FEE_ideal_delays[i] = 0.0;
  }

  beam_settings_t *beam_settings = malloc(sizeof(beam_settings_t));
  beam_settings->beamtype = FEE_BEAM_INTERP;

  visibility_set_t *visibility_set = test_calculate_visibilities(cropped_sky_models,
                                          beam_settings, woden_settings, RA0, -0.46606083776035967,
                                          beam_settings->beamtype);

  // printf("WODEN settings %d %d %d\n",woden_settings->num_time_steps,
  //                                    woden_settings->num_freqs,
  //                                    woden_settings->num_baselines );
  // for (int visi = 0; visi < NUM_VISI; visi++) {
  //   printf("%.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f\n",
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
  double gain1xx_re = 1.000000061109 * multiplier;
  double gain1xx_im = 0.0 * multiplier;
  double gain1xy_re = -0.000479380595 * multiplier;
  double gain1xy_im = -0.000005223917 * multiplier;
  double gain1yx_re = -0.000479380595 * multiplier;
  double gain1yx_im = 0.000005223917 * multiplier;
  double gain1yy_re = 1.000000062556 * multiplier;
  double gain1yy_im = 0.0 * multiplier;

  double gain2xx_re = 0.040796295396 * multiplier;
  double gain2xx_im = 0.0 * multiplier;
  double gain2xy_re = -0.002259936299 * multiplier;
  double gain2xy_im = 0.000231286354 * multiplier;
  double gain2yx_re = -0.002259936299 * multiplier;
  double gain2yx_im = -0.000231286354 * multiplier;
  double gain2yy_re = 0.020706076999 * multiplier;
  double gain2yy_im = 0.0 * multiplier;

  #ifdef DOUBLE_PRECISION
    double TOL = 1e-7;
  #else
    double TOL = 3e-2;
  #endif

  test_comp_phase_centre_allgains(visibility_set,
                                  gain1xx_re, gain1xx_im,
                                  gain1xy_re, gain1xy_im,
                                  gain1yx_re, gain1yx_im,
                                  gain1yy_re, gain1yy_im,
                                  gain2xx_re, gain2xx_im,
                                  gain2xy_re, gain2xy_im,
                                  gain2yx_re, gain2yx_im,
                                  gain2yy_re, gain2yy_im,
                                  woden_settings, TOL);

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
  test_comp_phase_centre_allgains(visibility_set,
                                  gain1xx_re, gain1xx_im,
                                  gain1xy_re, gain1xy_im,
                                  gain1yx_re, gain1yx_im,
                                  gain1yy_re, gain1yy_im,
                                  gain2xx_re, gain2xx_im,
                                  gain2xy_re, gain2xy_im,
                                  gain2yx_re, gain2yx_im,
                                  gain2yy_re, gain2yy_im,
                                  woden_settings, TOL);

  free_fee_beam(beam_settings->fee_beam);
  free(beam_settings);
  free_visi_set_inputs(visibility_set);
  free_visi_set_outputs(visibility_set);
  free(woden_settings);
}

//Test with a single SOURCE, single COMPONENT
void test_calculate_visibilities_MWAFEEBeamInterp_OneSource_SinglePoint(void) {
  int n_points = 1;
  int n_gauss = 0;
  int n_shapes = 0;
  int num_sources = 1;
  test_calculate_visibilities_MWAFEEBeamInterp(n_points, n_gauss, n_shapes, num_sources);
}

void test_calculate_visibilities_MWAFEEBeamInterp_OneSource_SingleGauss(void) {
  int n_points = 0;
  int n_gauss = 1;
  int n_shapes = 0;
  int num_sources = 1;
  test_calculate_visibilities_MWAFEEBeamInterp(n_points, n_gauss, n_shapes, num_sources);

}

void test_calculate_visibilities_MWAFEEBeamInterp_OneSource_SingleShape(void) {
  int n_points = 0;
  int n_gauss = 0;
  int n_shapes = 1;
  int num_sources = 1;
  test_calculate_visibilities_MWAFEEBeamInterp(n_points, n_gauss, n_shapes, num_sources);
}

void test_calculate_visibilities_MWAFEEBeamInterp_OneSource_SingleAll(void) {
  int n_points = 1;
  int n_gauss = 1;
  int n_shapes = 1;
  int num_sources = 1;
  test_calculate_visibilities_MWAFEEBeamInterp(n_points, n_gauss, n_shapes, num_sources);
}


//Test with a three SOURCEs, single COMPONENT
void test_calculate_visibilities_MWAFEEBeamInterp_ThreeSource_SinglePoint(void) {
  int n_points = 1;
  int n_gauss = 0;
  int n_shapes = 0;
  int num_sources = 3;
  test_calculate_visibilities_MWAFEEBeamInterp(n_points, n_gauss, n_shapes, num_sources);

}

void test_calculate_visibilities_MWAFEEBeamInterp_ThreeSource_SingleGauss(void) {
  int n_points = 0;
  int n_gauss = 1;
  int n_shapes = 0;
  int num_sources = 3;
  test_calculate_visibilities_MWAFEEBeamInterp(n_points, n_gauss, n_shapes, num_sources);
}

void test_calculate_visibilities_MWAFEEBeamInterp_ThreeSource_SingleShape(void) {
  int n_points = 0;
  int n_gauss = 0;
  int n_shapes = 1;
  int num_sources = 3;
  test_calculate_visibilities_MWAFEEBeamInterp(n_points, n_gauss, n_shapes, num_sources);
}

void test_calculate_visibilities_MWAFEEBeamInterp_ThreeSource_SingleAll(void) {
  int n_points = 1;
  int n_gauss = 1;
  int n_shapes = 1;
  int num_sources = 3;
  test_calculate_visibilities_MWAFEEBeamInterp(n_points, n_gauss, n_shapes, num_sources);
}


//Test with three SOURCEs, three COPMONENTs
void test_calculate_visibilities_MWAFEEBeamInterp_ThreeSource_ThreePoint(void) {
  int n_points = 3;
  int n_gauss = 0;
  int n_shapes = 0;
  int num_sources = 3;
  test_calculate_visibilities_MWAFEEBeamInterp(n_points, n_gauss, n_shapes, num_sources);

}

void test_calculate_visibilities_MWAFEEBeamInterp_ThreeSource_ThreeGauss(void) {
  int n_points = 0;
  int n_gauss = 3;
  int n_shapes = 0;
  int num_sources = 3;
  test_calculate_visibilities_MWAFEEBeamInterp(n_points, n_gauss, n_shapes, num_sources);
}

void test_calculate_visibilities_MWAFEEBeamInterp_ThreeSource_ThreeShape(void) {
  int n_points = 0;
  int n_gauss = 0;
  int n_shapes = 3;
  int num_sources = 3;
  test_calculate_visibilities_MWAFEEBeamInterp(n_points, n_gauss, n_shapes, num_sources);
}

void test_calculate_visibilities_MWAFEEBeamInterp_ThreeSource_ThreeAll(void) {
  int n_points = 3;
  int n_gauss = 3;
  int n_shapes = 3;
  int num_sources = 3;
  test_calculate_visibilities_MWAFEEBeamInterp(n_points, n_gauss, n_shapes, num_sources);
}



// Run the test with unity
int main(void)
{
    UNITY_BEGIN();

    char* mwa_fee_hdf5 = getenv("MWA_FEE_HDF5_INTERP");

    if (mwa_fee_hdf5) {
      printf("MWA_FEE_HDF5_INTERP: %s\n", mwa_fee_hdf5 );

          //Test with a single SOURCE, single COMPONENT
          RUN_TEST(test_calculate_visibilities_MWAFEEBeamInterp_OneSource_SinglePoint);
          RUN_TEST(test_calculate_visibilities_MWAFEEBeamInterp_OneSource_SingleGauss);
          RUN_TEST(test_calculate_visibilities_MWAFEEBeamInterp_OneSource_SingleShape);
          RUN_TEST(test_calculate_visibilities_MWAFEEBeamInterp_OneSource_SingleAll);

          //Test with three SOURCEs, single COPMONENT
          RUN_TEST(test_calculate_visibilities_MWAFEEBeamInterp_ThreeSource_SinglePoint);
          RUN_TEST(test_calculate_visibilities_MWAFEEBeamInterp_ThreeSource_SingleGauss);
          RUN_TEST(test_calculate_visibilities_MWAFEEBeamInterp_ThreeSource_SingleShape);
          RUN_TEST(test_calculate_visibilities_MWAFEEBeamInterp_ThreeSource_SingleAll);

          //Test with three SOURCEs, three COPMONENTs
          RUN_TEST(test_calculate_visibilities_MWAFEEBeamInterp_ThreeSource_ThreePoint);
          RUN_TEST(test_calculate_visibilities_MWAFEEBeamInterp_ThreeSource_ThreeGauss);
          RUN_TEST(test_calculate_visibilities_MWAFEEBeamInterp_ThreeSource_ThreeShape);
          RUN_TEST(test_calculate_visibilities_MWAFEEBeamInterp_ThreeSource_ThreeAll);
    }
    else {
      printf("MWA_FEE_HDF5_INTERP not found - not running test_calculate_visibilities_MWAFEEBeamInterp tests");
    }

    return UNITY_END();
}
