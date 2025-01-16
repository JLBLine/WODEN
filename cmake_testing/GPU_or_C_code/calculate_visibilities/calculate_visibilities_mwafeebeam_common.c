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

#include "calculate_visibilities_mwafeebeam_common.h"

void test_calculate_visibilities_MWAFEEBeam(int n_points, int n_gauss, int n_shapes,
                                           int num_sources, int do_gpu) {

  source_catalogue_t *cropped_sky_models = make_cropped_sky_models(RA0, -0.46606083776035967,
                                                    n_points, n_gauss, n_shapes,
                                                    num_sources);

  woden_settings_t *woden_settings = make_woden_settings(RA0, -0.46606083776035967);
  woden_settings->beamtype = FEE_BEAM;
  woden_settings->do_gpu = do_gpu;

  int delays[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  // int delays[16] = {2,4,6,8,2,4,6,8,2,4,6,8,2,4,6,8};
  woden_settings->FEE_ideal_delays = delays;
  
  beam_settings_t *beam_settings = malloc(sizeof(beam_settings_t));
  beam_settings->beamtype = FEE_BEAM;

  visibility_set_t *visibility_set = test_calculate_visibilities(cropped_sky_models,
                                          beam_settings, woden_settings, RA0, -0.46606083776035967,
                                          beam_settings->beamtype);

  //Small user_precision_t errors in the measurement equation mean that when at phase
  //centre, altough imaginary should be zero, you get fluctuations dependent
  //on u,v,w. So need a different tolerance here
  //The accuracy we get depends on whether we are float or double
  #ifdef DOUBLE_PRECISION
    // double TOL = 3e-9;
    double TOL = 4e-9;
  #else
    double TOL = 1e-5;
  #endif

  int num_comps = (n_points + n_gauss + n_shapes)*num_sources;

  double _Complex gain1x, leak1x, leak1y, gain1y;
  double _Complex gain2x, leak2x, leak2y, gain2y;

  //These are beam values taken from the double version of the beam
  gain1x = 0.877503615222 + I*-0.479570021294;
  leak1x = -0.000194011282 + I*-0.000034485292;
  leak1y = -0.000202113303 + I*-0.000027368942;
  gain1y = 0.877752206962 + I*-0.479114874769;
  gain2x = -0.071903981807 + I*0.041648775167;
  leak2x = 0.000032362088 + I*-0.000191146794;
  leak2y = 0.005000391879 + I*-0.002406349296;
  gain2y = -0.056871066138 + I*0.024991212027;
  
  test_comp_phase_centre_allgains(visibility_set, num_comps,
                                  gain1x, leak1x, leak1y, gain1y,
                                  gain2x, leak2x, leak2y, gain2y,
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

  test_comp_phase_centre_allgains(visibility_set, num_comps,
                                  gain1x, leak1x, leak1y, gain1y,
                                  gain2x, leak2x, leak2y, gain2y,
                                  woden_settings, TOL);

  free_fee_beam(beam_settings->fee_beam);
  free(beam_settings);
  free_visi_set_inputs(visibility_set);
  free_visi_set_outputs(visibility_set);
  free(woden_settings->mwa_dipole_amps);
  free(woden_settings);

}

//Test with a single SOURCE, single COMPONENT
void test_calculate_visibilities_MWAFEEBeam_OneSource_SinglePoint(int do_gpu) {
  int n_points = 1;
  int n_gauss = 0;
  int n_shapes = 0;
  int num_sources = 1;
  test_calculate_visibilities_MWAFEEBeam(n_points, n_gauss, n_shapes, num_sources,
                                          do_gpu);
}

void test_calculate_visibilities_MWAFEEBeam_OneSource_SingleGauss(int do_gpu) {
  int n_points = 0;
  int n_gauss = 1;
  int n_shapes = 0;
  int num_sources = 1;
  test_calculate_visibilities_MWAFEEBeam(n_points, n_gauss, n_shapes, num_sources,
                                          do_gpu);

}

void test_calculate_visibilities_MWAFEEBeam_OneSource_SingleShape(int do_gpu) {
  int n_points = 0;
  int n_gauss = 0;
  int n_shapes = 1;
  int num_sources = 1;
  test_calculate_visibilities_MWAFEEBeam(n_points, n_gauss, n_shapes, num_sources,
                                          do_gpu);
}

void test_calculate_visibilities_MWAFEEBeam_OneSource_SingleAll(int do_gpu) {
  int n_points = 1;
  int n_gauss = 1;
  int n_shapes = 1;
  int num_sources = 1;
  test_calculate_visibilities_MWAFEEBeam(n_points, n_gauss, n_shapes, num_sources,
                                          do_gpu);
}


//Test with a three SOURCEs, single COMPONENT
void test_calculate_visibilities_MWAFEEBeam_ThreeSource_SinglePoint(int do_gpu) {
  int n_points = 1;
  int n_gauss = 0;
  int n_shapes = 0;
  int num_sources = 3;
  test_calculate_visibilities_MWAFEEBeam(n_points, n_gauss, n_shapes, num_sources,
                                          do_gpu);

}

void test_calculate_visibilities_MWAFEEBeam_ThreeSource_SingleGauss(int do_gpu) {
  int n_points = 0;
  int n_gauss = 1;
  int n_shapes = 0;
  int num_sources = 3;
  test_calculate_visibilities_MWAFEEBeam(n_points, n_gauss, n_shapes, num_sources,
                                          do_gpu);
}

void test_calculate_visibilities_MWAFEEBeam_ThreeSource_SingleShape(int do_gpu) {
  int n_points = 0;
  int n_gauss = 0;
  int n_shapes = 1;
  int num_sources = 3;
  test_calculate_visibilities_MWAFEEBeam(n_points, n_gauss, n_shapes, num_sources,
                                          do_gpu);
}

void test_calculate_visibilities_MWAFEEBeam_ThreeSource_SingleAll(int do_gpu) {
  int n_points = 1;
  int n_gauss = 1;
  int n_shapes = 1;
  int num_sources = 3;
  test_calculate_visibilities_MWAFEEBeam(n_points, n_gauss, n_shapes, num_sources,
                                          do_gpu);
}


//Test with three SOURCEs, three COPMONENTs
void test_calculate_visibilities_MWAFEEBeam_ThreeSource_ThreePoint(int do_gpu) {
  int n_points = 5;
  int n_gauss = 0;
  int n_shapes = 0;
  int num_sources = 3;
  test_calculate_visibilities_MWAFEEBeam(n_points, n_gauss, n_shapes, num_sources,
                                          do_gpu);

}

void test_calculate_visibilities_MWAFEEBeam_ThreeSource_ThreeGauss(int do_gpu) {
  int n_points = 0;
  int n_gauss = 5;
  int n_shapes = 0;
  int num_sources = 3;
  test_calculate_visibilities_MWAFEEBeam(n_points, n_gauss, n_shapes, num_sources,
                                          do_gpu);
}

void test_calculate_visibilities_MWAFEEBeam_ThreeSource_ThreeShape(int do_gpu) {
  int n_points = 0;
  int n_gauss = 0;
  int n_shapes = 5;
  int num_sources = 3;
  test_calculate_visibilities_MWAFEEBeam(n_points, n_gauss, n_shapes, num_sources,
                                          do_gpu);
}

void test_calculate_visibilities_MWAFEEBeam_ThreeSource_ThreeAll(int do_gpu) {
  int n_points = 5;
  int n_gauss = 5;
  int n_shapes = 5;
  int num_sources = 3;
  test_calculate_visibilities_MWAFEEBeam(n_points, n_gauss, n_shapes, num_sources,
                                          do_gpu);
}