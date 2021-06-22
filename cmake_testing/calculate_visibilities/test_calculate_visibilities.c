#include <math.h>
#include <unity.h>
#include <stdlib.h>
#include <complex.h>

#include "constants.h"
#include "woden_struct_defs.h"
#include "FEE_primary_beam.h"
#include "shapelet_basis.h"
#include "woden_settings.h"
#include "visibility_set.h"
#include "shapelet_basis.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

//External CUDA code we're linking in
extern void calculate_visibilities(array_layout_t *array_layout,
  source_catalogue_t *cropped_sky_models, beam_settings_t *beam_settings,
  woden_settings_t *woden_settings, visibility_set_t *visibility_set,
  float *sbf);

#define NUM_BASELINES 3
#define NUM_FREQS 3
#define NUM_TIME_STEPS 2

#define UNITY_INCLUDE_FLOAT

//Different delays settings, which control the pointing of the MWA beam
float zenith_delays[16] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                           0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

void test_calculate_visibilities(source_catalogue_t *cropped_sky_models,
                                 beam_settings_t *beam_settings,
                                 woden_settings_t *woden_settings,
                                 int beamtype, float base_band_freq) {

  array_layout_t *array_layout = malloc(sizeof(array_layout_t));

  array_layout->X_diff_metres = malloc(NUM_BASELINES*sizeof(float));
  array_layout->Y_diff_metres = malloc(NUM_BASELINES*sizeof(float));
  array_layout->Z_diff_metres = malloc(NUM_BASELINES*sizeof(float));

  for (int baseline = 0; baseline < NUM_BASELINES; baseline++) {
    array_layout->X_diff_metres[baseline] = (baseline + 1) * 100;
    array_layout->Y_diff_metres[baseline] = (baseline + 1) * 100;
    array_layout->Z_diff_metres[baseline] = 0.0;
  }

  int status = 0;

  //The intial setup of the FEE beam is done on the CPU, so call it here
  if (woden_settings->beamtype == FEE_BEAM){

    beam_settings->FEE_beam = malloc(sizeof(RTS_MWA_FEE_beam_t));
      //We need the zenith beam to get the normalisation
    beam_settings->FEE_beam_zenith = malloc(sizeof(RTS_MWA_FEE_beam_t));

    float base_middle_freq = base_band_freq + (woden_settings->coarse_band_width/2.0);
  //
    printf("Middle freq is %.8e \n",base_middle_freq );
  //
    float float_zenith_delays[16] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  //
    printf("Setting up the zenith FEE beam...");
    status = RTS_MWAFEEInit(woden_settings->hdf5_beam_path, base_middle_freq,
                            beam_settings->FEE_beam_zenith, zenith_delays);
    printf(" done.\n");

    printf("Setting up the FEE beam...");
    status = RTS_MWAFEEInit(woden_settings->hdf5_beam_path, base_middle_freq,
          beam_settings->FEE_beam, woden_settings->FEE_ideal_delays);
    printf(" done.\n");

  }

  float *sbf = NULL;
  if (cropped_sky_models->num_shapelets > 0) {
    sbf = malloc( sbf_N * sbf_L * sizeof(float) );
    sbf = create_sbf(sbf);
  }

  visibility_set_t *visibility_set = setup_visibility_set(woden_settings->num_visis);

  float lsts[] = {0.0, M_PI / 2};

  fill_timefreq_visibility_set(visibility_set, woden_settings,
                               base_band_freq, lsts);

  calculate_visibilities(array_layout, cropped_sky_models, beam_settings,
                         woden_settings, visibility_set, sbf);

  for (size_t visi = 0; visi < woden_settings->num_visis; visi++) {
    printf("%.4f %.4f %.4f %.1f %.1f %.1f %.1f %.1f %.1f %.1f %.1f\n",
            visibility_set->us_metres[visi],
            visibility_set->vs_metres[visi],
            visibility_set->ws_metres[visi],
            visibility_set->sum_visi_XX_real[visi],
            visibility_set->sum_visi_XX_imag[visi],
            visibility_set->sum_visi_XY_real[visi],
            visibility_set->sum_visi_XY_imag[visi],
            visibility_set->sum_visi_YX_real[visi],
            visibility_set->sum_visi_YX_imag[visi],
            visibility_set->sum_visi_YY_real[visi],
            visibility_set->sum_visi_YY_imag[visi]);
  }
}

void test_calculate_visibilities_NoBeam_Point_PhaseCent() {
  source_catalogue_t *cropped_sky_models = malloc(sizeof(cropped_sky_models));

  cropped_sky_models->num_sources = 1;
  cropped_sky_models->num_shapelets = 0;
  cropped_sky_models->catsources = malloc(sizeof(catsource_t));

  float one_array[] = {1.0};
  float zero_array[] = {0.0};

  cropped_sky_models->catsources[0].n_comps = 1;
  cropped_sky_models->catsources[0].n_points = 1;
  cropped_sky_models->catsources[0].n_gauss = 0;
  cropped_sky_models->catsources[0].n_shapes = 0;
  cropped_sky_models->catsources[0].n_shape_coeffs = 0;

  cropped_sky_models->catsources[0].point_ref_stokesI = one_array;
  cropped_sky_models->catsources[0].point_ref_stokesQ = zero_array;
  cropped_sky_models->catsources[0].point_ref_stokesU = zero_array;
  cropped_sky_models->catsources[0].point_ref_stokesV = zero_array;
  cropped_sky_models->catsources[0].point_SIs = zero_array;
  // // cropped_sky_models->catsource[0].point_azs = zero_array;
  // // cropped_sky_models->catsource[0].point_zas = zero_array;
  //
  float ra0 = 0.0;
  float dec0 = MWA_LAT_RAD;

  float ra_array[] = {ra0};
  float dec_array[] = {dec0};
  float freq_array[] = {150e+6};
  //
  cropped_sky_models->catsources[0].point_ras = ra_array;
  cropped_sky_models->catsources[0].point_decs = dec_array;
  cropped_sky_models->catsources[0].point_ref_freqs = freq_array;

  beam_settings_t *beam_settings = malloc(sizeof(beam_settings));

  beam_settings->beamtype = NO_BEAM;

  woden_settings_t *woden_settings = malloc(sizeof(woden_settings_t));

  woden_settings->ra0 = ra0;
  woden_settings->dec0 = dec0;
  woden_settings->sdec0 = sinf(dec0);
  woden_settings->cdec0 = cosf(dec0);
  woden_settings->num_baselines = NUM_BASELINES;
  woden_settings->num_freqs = NUM_FREQS;
  woden_settings->num_time_steps = NUM_TIME_STEPS;
  woden_settings->beamtype = NO_BEAM;
  woden_settings->num_visis = NUM_BASELINES * NUM_FREQS * NUM_TIME_STEPS;
  woden_settings->coarse_band_width = 1.28e+6;

  float base_band_freq = 120e+6;

  test_calculate_visibilities(cropped_sky_models,
                              beam_settings, woden_settings,
                              beam_settings->beamtype, base_band_freq);

}



// /*
// Check whether the environment variable for the FEE hdf5 beam exists, don't run
// the test if it's missing
// */
// void check_for_env_and_run_test(float freq, float *delays, float *expected,
//                                 char *outname) {
//   char* mwa_fee_hdf5 = getenv("MWA_FEE_HDF5");
//
//   if (mwa_fee_hdf5) {
//     printf("MWA_FEE_HDF5: %s", mwa_fee_hdf5 );
//     test_RTS_CUDA_FEE_beam_VaryFreqVaryPointing(freq, delays, mwa_fee_hdf5, expected, outname);
//   }
//   else {
//     printf("MWA_FEE_HDF5 not found - not running test_RTS_FEE_beam test");
//   }
// }
//
// /*
// Run the test but vary the frequency and pointings. Compare to pre-calculated
// values that are stored in test_RTS_FEE_beam.h
// */
// void test_RTS_CUDA_FEE_beam_100MHz_zenith(void) {
//   check_for_env_and_run_test(100e+6, zenith_delays, zenith_100, "zenith_100.txt");
// }


//Run the test with unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_calculate_visibilities_NoBeam_Point_PhaseCent);

    return UNITY_END();
}
