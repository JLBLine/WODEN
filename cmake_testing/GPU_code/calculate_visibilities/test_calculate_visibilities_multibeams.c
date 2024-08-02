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

  woden_settings->use_dipamps = 1;

  // int num_delays_per_tile = 16;
  int delays[48] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

  //These are the amplitudes for the dipoles, as read in from metafits
    // I believe that they have X - east-west, Y - north-south
  double amps[96] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
    0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4,
    0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6,
    0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

  // double amps[96] = {0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
  //   0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
  //   0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
  //   0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
  //   0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
  //   0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2};

  woden_settings->mwa_dipole_amps = amps;

  // int num_beams = NUM_ANTS;

  woden_settings->FEE_ideal_delays = delays;

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

  int num_comps = (n_points + n_gauss + n_shapes)*num_sources;

  double _Complex gain1x, leak1x, leak1y, gain1y;
  double _Complex gain2x, leak2x, leak2y, gain2y;

  gain1x = 0.947431083057 + I*0.319959908285;
  leak1x = -0.000246714867 + I*-0.000016523331;
  leak1y = -0.000249894393 + I*-0.000011798344;
  gain1y = 0.947019358550 + I*0.321176485103;
  gain2x = -0.198505432118 + I*-0.037188871831;
  leak2x = -0.001189280798 + I*-0.000248382900;
  leak2y = 0.011616550617 + I*0.003414333525;
  gain2y = -0.137450706771 + I*-0.040789246826;

  #ifdef DOUBLE_PRECISION
    double TOL = 8e-7;
  #else
    double TOL = 3e-2;
  #endif

  double antx_mult[3] = {0.2, 0.6, 1.0};
  double anty_mult[3] = {0.0, 0.4, 0.8};

  test_comp_phase_centre_allgains_multiants(visibility_set,num_comps,
                                  gain1x, leak1x, leak1y, gain1y,
                                  gain2x, leak2x, leak2y, gain2y,
                                  antx_mult, anty_mult, NUM_ANTS,
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

  test_comp_phase_centre_allgains_multiants(visibility_set,num_comps,
                                  gain1x, leak1x, leak1y, gain1y,
                                  gain2x, leak2x, leak2y, gain2y,
                                  antx_mult, anty_mult, NUM_ANTS,
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

          // //Test with three SOURCEs, single COPMONENT
          RUN_TEST(test_calculate_visibilities_MWAFEEBeamInterp_ThreeSource_SinglePoint);
          RUN_TEST(test_calculate_visibilities_MWAFEEBeamInterp_ThreeSource_SingleGauss);
          RUN_TEST(test_calculate_visibilities_MWAFEEBeamInterp_ThreeSource_SingleShape);
          RUN_TEST(test_calculate_visibilities_MWAFEEBeamInterp_ThreeSource_SingleAll);

          // //Test with three SOURCEs, three COPMONENTs
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
