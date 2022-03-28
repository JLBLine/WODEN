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
#include <mwa_hyperbeam.h>

// //External CUDA code we're linking in
// extern void calculate_visibilities(array_layout_t *array_layout,
//   source_catalogue_t *cropped_sky_models, beam_settings_t *beam_settings,
//   woden_settings_t *woden_settings, visibility_set_t *visibility_set,
//   user_precision_t *sbf);

//Just two LSTs to check things change with time
double lsts[] = {0.0, M_PI / 4};

user_precision_t azs[] = {0.0, 4.528359553989764};
user_precision_t zas[] = {0.0, 0.6978088917603547};

//Different delays settings, which control the pointing of the MWA beam
user_precision_t zenith_delays[16] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                           0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};


/*
Create a number of SOURCEs and input into a sky model catalogue `source_catalogue_t`
struct. For each SOURCE, populate with as many COMPONENTs as requested of
whatever combination of POINT, GAUSSIAN, and SHAPELET types
Keep everything to just Stokes I 1 Jy, and stick the source at phase centre

We'll stick with just one shapelet coeff per shapelet compoenent as there
are other tests to make sure SHAPELETs works elsewhere
*/
source_catalogue_t * make_cropped_sky_models(double ra0, double dec0,
                                             int n_points, int n_gauss,
                                             int n_shapes,
                                             int num_sources) {

  int n_comps = n_points + n_gauss + n_shapes;

  source_catalogue_t *cropped_sky_models = malloc(sizeof(cropped_sky_models));
  cropped_sky_models->num_sources = num_sources;
  cropped_sky_models->num_shapelets = n_shapes*n_comps;
  cropped_sky_models->sources = malloc(num_sources*sizeof(source_t));

  for (int cats_ind = 0; cats_ind < num_sources; cats_ind++) {
    cropped_sky_models->sources[cats_ind].n_comps = n_comps;
    cropped_sky_models->sources[cats_ind].n_points = n_points;
    cropped_sky_models->sources[cats_ind].n_gauss = n_gauss;
    cropped_sky_models->sources[cats_ind].n_shapes = n_shapes;
    cropped_sky_models->sources[cats_ind].n_shape_coeffs = n_shapes;

    if (n_points > 0) {
      cropped_sky_models->sources[cats_ind].point_components.num_primarybeam_values = n_points*NUM_FREQS*NUM_TIME_STEPS;
      cropped_sky_models->sources[cats_ind].point_components.ref_stokesI = malloc(n_points*sizeof(user_precision_t));
      cropped_sky_models->sources[cats_ind].point_components.ref_stokesQ = malloc(n_points*sizeof(user_precision_t));
      cropped_sky_models->sources[cats_ind].point_components.ref_stokesU = malloc(n_points*sizeof(user_precision_t));
      cropped_sky_models->sources[cats_ind].point_components.ref_stokesV = malloc(n_points*sizeof(user_precision_t));
      cropped_sky_models->sources[cats_ind].point_components.SIs = malloc(n_points*sizeof(user_precision_t));

      cropped_sky_models->sources[cats_ind].point_components.ras = malloc(n_points*sizeof(double));
      cropped_sky_models->sources[cats_ind].point_components.decs = malloc(n_points*sizeof(double));
      cropped_sky_models->sources[cats_ind].point_components.ref_freqs = malloc(n_points*sizeof(double));

      cropped_sky_models->sources[cats_ind].point_components.beam_has = malloc(n_points*NUM_TIME_STEPS*sizeof(double));
      cropped_sky_models->sources[cats_ind].point_components.beam_decs = malloc(n_points*NUM_TIME_STEPS*sizeof(double));

      cropped_sky_models->sources[cats_ind].point_components.azs = malloc(n_points*NUM_TIME_STEPS*sizeof(user_precision_t));
      cropped_sky_models->sources[cats_ind].point_components.zas = malloc(n_points*NUM_TIME_STEPS*sizeof(user_precision_t));

      for (int point = 0; point < n_points; point++) {
        cropped_sky_models->sources[cats_ind].point_components.ras[point] = ra0;
        cropped_sky_models->sources[cats_ind].point_components.decs[point] = dec0;
        cropped_sky_models->sources[cats_ind].point_components.ref_freqs[point] = 150e+6;

        cropped_sky_models->sources[cats_ind].point_components.ref_stokesI[point] = STOKESI;
        cropped_sky_models->sources[cats_ind].point_components.ref_stokesQ[point] = 0.0;
        cropped_sky_models->sources[cats_ind].point_components.ref_stokesU[point] = 0.0;
        cropped_sky_models->sources[cats_ind].point_components.ref_stokesV[point] = 0.0;
        cropped_sky_models->sources[cats_ind].point_components.SIs[point] = 0.0;

        for (int time = 0; time < NUM_TIME_STEPS; time++) {
          int step = point*NUM_TIME_STEPS + time;
          cropped_sky_models->sources[cats_ind].point_components.beam_has[step] = lsts[time] - RA0;
          cropped_sky_models->sources[cats_ind].point_components.beam_decs[step] = MWA_LAT_RAD;

          cropped_sky_models->sources[cats_ind].point_components.azs[step] = azs[time];
          cropped_sky_models->sources[cats_ind].point_components.zas[step] = zas[time];
        }


      }
    }

    if (n_gauss > 0) {
      cropped_sky_models->sources[cats_ind].gauss_components.num_primarybeam_values = n_gauss*NUM_FREQS*NUM_TIME_STEPS;
      cropped_sky_models->sources[cats_ind].gauss_components.ref_stokesI = malloc(n_gauss*sizeof(user_precision_t));
      cropped_sky_models->sources[cats_ind].gauss_components.ref_stokesQ = malloc(n_gauss*sizeof(user_precision_t));
      cropped_sky_models->sources[cats_ind].gauss_components.ref_stokesU = malloc(n_gauss*sizeof(user_precision_t));
      cropped_sky_models->sources[cats_ind].gauss_components.ref_stokesV = malloc(n_gauss*sizeof(user_precision_t));
      cropped_sky_models->sources[cats_ind].gauss_components.SIs = malloc(n_gauss*sizeof(user_precision_t));

      cropped_sky_models->sources[cats_ind].gauss_components.ras = malloc(n_gauss*sizeof(double));
      cropped_sky_models->sources[cats_ind].gauss_components.decs = malloc(n_gauss*sizeof(double));
      cropped_sky_models->sources[cats_ind].gauss_components.ref_freqs = malloc(n_gauss*sizeof(double));

      cropped_sky_models->sources[cats_ind].gauss_components.majors = malloc(n_gauss*sizeof(user_precision_t));
      cropped_sky_models->sources[cats_ind].gauss_components.minors = malloc(n_gauss*sizeof(user_precision_t));
      cropped_sky_models->sources[cats_ind].gauss_components.pas = malloc(n_gauss*sizeof(user_precision_t));

      cropped_sky_models->sources[cats_ind].gauss_components.beam_has = malloc(n_gauss*NUM_TIME_STEPS*sizeof(double));
      cropped_sky_models->sources[cats_ind].gauss_components.beam_decs = malloc(n_gauss*NUM_TIME_STEPS*sizeof(double));

      cropped_sky_models->sources[cats_ind].gauss_components.azs = malloc(n_gauss*NUM_TIME_STEPS*sizeof(user_precision_t));
      cropped_sky_models->sources[cats_ind].gauss_components.zas = malloc(n_gauss*NUM_TIME_STEPS*sizeof(user_precision_t));

      for (int gauss = 0; gauss < n_gauss; gauss++) {
        cropped_sky_models->sources[cats_ind].gauss_components.ras[gauss] = ra0;
        cropped_sky_models->sources[cats_ind].gauss_components.decs[gauss] = dec0;
        cropped_sky_models->sources[cats_ind].gauss_components.ref_freqs[gauss] = 150e+6;

        cropped_sky_models->sources[cats_ind].gauss_components.ref_stokesI[gauss] = STOKESI;
        cropped_sky_models->sources[cats_ind].gauss_components.ref_stokesQ[gauss] = 0.0;
        cropped_sky_models->sources[cats_ind].gauss_components.ref_stokesU[gauss] = 0.0;
        cropped_sky_models->sources[cats_ind].gauss_components.ref_stokesV[gauss] = 0.0;
        cropped_sky_models->sources[cats_ind].gauss_components.SIs[gauss] = 0.0;

        cropped_sky_models->sources[cats_ind].gauss_components.majors[gauss] = 1e-10;
        cropped_sky_models->sources[cats_ind].gauss_components.minors[gauss] = 1e-10;
        cropped_sky_models->sources[cats_ind].gauss_components.pas[gauss] = 0;

        for (int time = 0; time < NUM_TIME_STEPS; time++) {
          int step = gauss*NUM_TIME_STEPS + time;
          cropped_sky_models->sources[cats_ind].gauss_components.beam_has[step] = lsts[time] - RA0;
          cropped_sky_models->sources[cats_ind].gauss_components.beam_decs[step] = MWA_LAT_RAD;

          cropped_sky_models->sources[cats_ind].gauss_components.azs[step] = azs[time];
          cropped_sky_models->sources[cats_ind].gauss_components.zas[step] = zas[time];
        }
      }
    }

    if (n_shapes > 0) {
      cropped_sky_models->sources[cats_ind].shape_components.num_primarybeam_values = n_shapes*NUM_FREQS*NUM_TIME_STEPS;
      cropped_sky_models->sources[cats_ind].shape_components.ref_stokesI = malloc(n_shapes*sizeof(user_precision_t));
      cropped_sky_models->sources[cats_ind].shape_components.ref_stokesQ = malloc(n_shapes*sizeof(user_precision_t));
      cropped_sky_models->sources[cats_ind].shape_components.ref_stokesU = malloc(n_shapes*sizeof(user_precision_t));
      cropped_sky_models->sources[cats_ind].shape_components.ref_stokesV = malloc(n_shapes*sizeof(user_precision_t));
      cropped_sky_models->sources[cats_ind].shape_components.SIs = malloc(n_shapes*sizeof(user_precision_t));

      cropped_sky_models->sources[cats_ind].shape_components.ras = malloc(n_shapes*sizeof(double));
      cropped_sky_models->sources[cats_ind].shape_components.decs = malloc(n_shapes*sizeof(double));
      cropped_sky_models->sources[cats_ind].shape_components.ref_freqs = malloc(n_shapes*sizeof(double));

      cropped_sky_models->sources[cats_ind].shape_components.majors = malloc(n_shapes*sizeof(user_precision_t));
      cropped_sky_models->sources[cats_ind].shape_components.minors = malloc(n_shapes*sizeof(user_precision_t));
      cropped_sky_models->sources[cats_ind].shape_components.pas = malloc(n_shapes*sizeof(user_precision_t));

      cropped_sky_models->sources[cats_ind].shape_components.shape_coeffs = malloc(n_shapes*sizeof(user_precision_t));
      cropped_sky_models->sources[cats_ind].shape_components.n1s = malloc(n_shapes*sizeof(user_precision_t));
      cropped_sky_models->sources[cats_ind].shape_components.n2s = malloc(n_shapes*sizeof(user_precision_t));
      cropped_sky_models->sources[cats_ind].shape_components.param_indexes = malloc(n_shapes*sizeof(user_precision_t));

      cropped_sky_models->sources[cats_ind].shape_components.beam_has = malloc(n_shapes*NUM_TIME_STEPS*sizeof(double));
      cropped_sky_models->sources[cats_ind].shape_components.beam_decs = malloc(n_shapes*NUM_TIME_STEPS*sizeof(double));

      cropped_sky_models->sources[cats_ind].shape_components.azs = malloc(n_shapes*NUM_TIME_STEPS*sizeof(user_precision_t));
      cropped_sky_models->sources[cats_ind].shape_components.zas = malloc(n_shapes*NUM_TIME_STEPS*sizeof(user_precision_t));

      for (int shape = 0; shape < n_shapes; shape++) {
        cropped_sky_models->sources[cats_ind].shape_components.ras[shape] = ra0;
        cropped_sky_models->sources[cats_ind].shape_components.decs[shape] = dec0;
        cropped_sky_models->sources[cats_ind].shape_components.ref_freqs[shape] = 150e+6;

        cropped_sky_models->sources[cats_ind].shape_components.ref_stokesI[shape] = STOKESI;
        cropped_sky_models->sources[cats_ind].shape_components.ref_stokesQ[shape] = 0.0;
        cropped_sky_models->sources[cats_ind].shape_components.ref_stokesU[shape] = 0.0;
        cropped_sky_models->sources[cats_ind].shape_components.ref_stokesV[shape] = 0.0;
        cropped_sky_models->sources[cats_ind].shape_components.SIs[shape] = 0.0;

        cropped_sky_models->sources[cats_ind].shape_components.majors[shape] = 1e-10;
        cropped_sky_models->sources[cats_ind].shape_components.minors[shape] = 1e-10;
        cropped_sky_models->sources[cats_ind].shape_components.pas[shape] = 0;

        //These settings basically make this shapelet a GAUSSIAN
        cropped_sky_models->sources[cats_ind].shape_components.shape_coeffs[shape] = 1.0;
        cropped_sky_models->sources[cats_ind].shape_components.n1s[shape] = 0.0;
        cropped_sky_models->sources[cats_ind].shape_components.n2s[shape] = 0.0;
        cropped_sky_models->sources[cats_ind].shape_components.param_indexes[shape] = 0.0;

        for (int time = 0; time < NUM_TIME_STEPS; time++) {
          int step = shape*NUM_TIME_STEPS + time;
          cropped_sky_models->sources[cats_ind].shape_components.beam_has[step] = lsts[time] - RA0;
          cropped_sky_models->sources[cats_ind].shape_components.beam_decs[step] = MWA_LAT_RAD;

          cropped_sky_models->sources[cats_ind].shape_components.azs[step] = azs[time];
          cropped_sky_models->sources[cats_ind].shape_components.zas[step] = zas[time];
        }
      }
    }

  }
  return cropped_sky_models;
}

/*
Based on what is in the sky model, free stuff
*/
void free_sky_model(source_catalogue_t *cropped_sky_models) {

  for (int cats_ind = 0; cats_ind < cropped_sky_models->num_sources; cats_ind++) {

    int n_points = cropped_sky_models->sources[cats_ind].n_points;
    int n_gauss = cropped_sky_models->sources[cats_ind].n_gauss;
    int n_shapes = cropped_sky_models->sources[cats_ind].n_shapes;

    if (n_points > 0) {
      free(cropped_sky_models->sources[cats_ind].point_components.ref_stokesI);
      free(cropped_sky_models->sources[cats_ind].point_components.ref_stokesQ);
      free(cropped_sky_models->sources[cats_ind].point_components.ref_stokesU);
      free(cropped_sky_models->sources[cats_ind].point_components.ref_stokesV);
      free(cropped_sky_models->sources[cats_ind].point_components.SIs);

      free(cropped_sky_models->sources[cats_ind].point_components.ras);
      free(cropped_sky_models->sources[cats_ind].point_components.decs);
      free(cropped_sky_models->sources[cats_ind].point_components.ref_freqs);

      free(cropped_sky_models->sources[cats_ind].point_components.beam_has);
      free(cropped_sky_models->sources[cats_ind].point_components.beam_decs);

      free(cropped_sky_models->sources[cats_ind].point_components.azs);
      free(cropped_sky_models->sources[cats_ind].point_components.zas);
    }

    if (n_gauss > 0) {
      free(cropped_sky_models->sources[cats_ind].gauss_components.ref_stokesI);
      free(cropped_sky_models->sources[cats_ind].gauss_components.ref_stokesQ);
      free(cropped_sky_models->sources[cats_ind].gauss_components.ref_stokesU);
      free(cropped_sky_models->sources[cats_ind].gauss_components.ref_stokesV);
      free(cropped_sky_models->sources[cats_ind].gauss_components.SIs);

      free(cropped_sky_models->sources[cats_ind].gauss_components.ras);
      free(cropped_sky_models->sources[cats_ind].gauss_components.decs);
      free(cropped_sky_models->sources[cats_ind].gauss_components.ref_freqs);

      free(cropped_sky_models->sources[cats_ind].gauss_components.majors);
      free(cropped_sky_models->sources[cats_ind].gauss_components.minors);
      free(cropped_sky_models->sources[cats_ind].gauss_components.pas);

      free(cropped_sky_models->sources[cats_ind].gauss_components.beam_has);
      free(cropped_sky_models->sources[cats_ind].gauss_components.beam_decs);

      free(cropped_sky_models->sources[cats_ind].gauss_components.azs);
      free(cropped_sky_models->sources[cats_ind].gauss_components.zas);
    }

    if (n_shapes > 0) {
      free(cropped_sky_models->sources[cats_ind].shape_components.ref_stokesI);
      free(cropped_sky_models->sources[cats_ind].shape_components.ref_stokesQ);
      free(cropped_sky_models->sources[cats_ind].shape_components.ref_stokesU);
      free(cropped_sky_models->sources[cats_ind].shape_components.ref_stokesV);
      free(cropped_sky_models->sources[cats_ind].shape_components.SIs);

      free(cropped_sky_models->sources[cats_ind].shape_components.ras);
      free(cropped_sky_models->sources[cats_ind].shape_components.decs);
      free(cropped_sky_models->sources[cats_ind].shape_components.ref_freqs);

      free(cropped_sky_models->sources[cats_ind].shape_components.majors);
      free(cropped_sky_models->sources[cats_ind].shape_components.minors);
      free(cropped_sky_models->sources[cats_ind].shape_components.pas);

      free(cropped_sky_models->sources[cats_ind].shape_components.shape_coeffs);
      free(cropped_sky_models->sources[cats_ind].shape_components.n1s);
      free(cropped_sky_models->sources[cats_ind].shape_components.n2s);
      free(cropped_sky_models->sources[cats_ind].shape_components.param_indexes);

      free(cropped_sky_models->sources[cats_ind].shape_components.beam_has);
      free(cropped_sky_models->sources[cats_ind].shape_components.beam_decs);

      free(cropped_sky_models->sources[cats_ind].shape_components.azs);
      free(cropped_sky_models->sources[cats_ind].shape_components.zas);

    }
  }
  free(cropped_sky_models->sources);
  free(cropped_sky_models);
}

woden_settings_t * make_woden_settings(double ra0, double dec0) {

    woden_settings_t *woden_settings = malloc(sizeof(woden_settings_t));

    woden_settings->ra0 = ra0;
    woden_settings->dec0 = dec0;
    woden_settings->sdec0 = sin(dec0);
    woden_settings->cdec0 = cos(dec0);
    woden_settings->num_baselines = NUM_BASELINES;
    woden_settings->num_freqs = NUM_FREQS;
    woden_settings->num_time_steps = NUM_TIME_STEPS;
    woden_settings->num_visis = NUM_BASELINES * NUM_FREQS * NUM_TIME_STEPS;
    woden_settings->coarse_band_width = 1.28e+6;
    //Make the fine channel width insanely small so beam changes little
    //with frequency - that way we can test for just one gain value per time
    woden_settings->frequency_resolution = 1e-6;
    woden_settings->latitude = MWA_LAT_RAD;

    return woden_settings;
}

/*
Due to the way I've setup the array layout and phase centre,
the uvw should be a simple reduced product between phase centre and lsts.
Create some expected u,v,w arrays and compare to what we've found
*/
void test_uvw(visibility_set_t *visibility_set,  double *lsts,
              double ra0, double dec0) {

  int num_visis = NUM_BASELINES * NUM_FREQS * NUM_TIME_STEPS;
  double *expec_u = malloc(num_visis * sizeof(double));
  double *expec_v = malloc(num_visis * sizeof(double));
  double *expec_w = malloc(num_visis * sizeof(double));

  int index = 0;
  double xy_length;
  double ha0;
  for (int time = 0; time < NUM_TIME_STEPS; time++) {
    ha0 = lsts[time] - ra0;
    for (int freq = 0; freq < NUM_FREQS; freq++) {
      for (int baseline = 0; baseline < NUM_BASELINES; baseline++) {
        xy_length = (baseline + 1) * 100;
        expec_u[index] = (cos(ha0) + sin(ha0))*xy_length;
        expec_v[index] = xy_length * sin(MWA_LAT_RAD)*(-cos(ha0) + sin(ha0));
        expec_w[index] = xy_length * cos(MWA_LAT_RAD)*(cos(ha0) - sin(ha0));

        index ++;
      }
    }
  }

  #ifdef DOUBLE_PRECISION
    double TOL = 1e-12;
  #else
    double TOL = 1e-5;
  #endif

  for (int visi = 0; visi < num_visis; visi++) {
    // printf("%.1f %.1f %.1f %.1f %.1f %.1f\n",
    //         visibility_set->us_metres[visi], expec_u[visi],
    //         visibility_set->vs_metres[visi], expec_v[visi],
    //         visibility_set->ws_metres[visi], expec_w[visi] );
    TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_u[visi],
                              visibility_set->us_metres[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_v[visi],
                              visibility_set->vs_metres[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_w[visi],
                              visibility_set->ws_metres[visi]);

  }

  free(expec_u);
  free(expec_v);
  free(expec_w);
}

/*
Pump many many many settings into the function we are trying to test
Checkthe u,v,w are correct and return the visibility_set for further testing
*/
visibility_set_t * test_calculate_visibilities(source_catalogue_t *cropped_sky_models,
                                 beam_settings_t *beam_settings,
                                 woden_settings_t *woden_settings,
                                 double ra0, double dec0,
                                 int beamtype) {

  double base_band_freq = BASE_BAND_FREQ;
  // user_precision_t base_band_freq = 120e+6;

  array_layout_t *array_layout = malloc(sizeof(array_layout_t));

  array_layout->X_diff_metres = malloc(NUM_BASELINES*sizeof(double));
  array_layout->Y_diff_metres = malloc(NUM_BASELINES*sizeof(double));
  array_layout->Z_diff_metres = malloc(NUM_BASELINES*sizeof(double));

  for (int baseline = 0; baseline < NUM_BASELINES; baseline++) {
    array_layout->X_diff_metres[baseline] = (baseline + 1) * 100;
    array_layout->Y_diff_metres[baseline] = (baseline + 1) * 100;
    array_layout->Z_diff_metres[baseline] = 0.0;
  }

  user_precision_t *sbf = NULL;
  if (cropped_sky_models->num_shapelets > 0) {
    sbf = malloc( sbf_N * sbf_L * sizeof(user_precision_t) );
    sbf = create_sbf(sbf);
  }

  visibility_set_t *visibility_set = setup_visibility_set(woden_settings->num_visis);

  fill_timefreq_visibility_set(visibility_set, woden_settings,
                               base_band_freq, lsts);

  if (beam_settings->beamtype == FEE_BEAM) {
    char* mwa_fee_hdf5 = getenv("MWA_FEE_HDF5");

    if (mwa_fee_hdf5) {
      printf("MWA_FEE_HDF5: %s\n", mwa_fee_hdf5 );

      int status =  new_fee_beam(mwa_fee_hdf5, &beam_settings->fee_beam,
                                 beam_settings->hyper_error_str);
      if (status != 0) {
        printf("hyperbeam error %d %s\n", status, beam_settings->hyper_error_str );
      }

    } else {
      printf("MWA_FEE_HDF5 not found - not running test_hyperbeam test");
    }
  }
  else if (beam_settings->beamtype == FEE_BEAM_INTERP) {
    char* mwa_fee_hdf5 = getenv("MWA_FEE_HDF5_INTERP");

    if (mwa_fee_hdf5) {
      printf("MWA_FEE_HDF5: %s\n", mwa_fee_hdf5 );

      int status =  new_fee_beam(mwa_fee_hdf5, &beam_settings->fee_beam,
                                 beam_settings->hyper_error_str);
      if (status != 0) {
        printf("hyperbeam error %d %s\n", status, beam_settings->hyper_error_str );
      }

    } else {
      printf("MWA_FEE_HDF5 not found - not running test_hyperbeam test");
    }
  }

  printf("Calling CUDA\n");
  calculate_visibilities(array_layout, cropped_sky_models, beam_settings,
                         woden_settings, visibility_set, sbf);
  printf("CUDA has finished\n");

  //Be free my pretties!
  if (cropped_sky_models->num_shapelets > 0) {
    free(sbf);
  }

  free_sky_model(cropped_sky_models);
  free(array_layout->X_diff_metres);
  free(array_layout->Y_diff_metres);
  free(array_layout->Z_diff_metres);
  // free(array_layout);

  test_uvw(visibility_set, lsts, ra0, dec0);

  return visibility_set;

}

/*
This test works for beams that change with time, but have no cross-pol info
Can just test whether the XX/YY real values match expected
*/
void test_comp_phase_centre_twogains(visibility_set_t *visibility_set,
                                     double gain1xx, double gain1yy,
                                     double gain2xx, double gain2yy) {

  double *expec_gainsxx = malloc(NUM_VISI*sizeof(double));
  double *expec_gainsyy = malloc(NUM_VISI*sizeof(double));

  double gainsxx[] = {gain1xx, gain2xx};
  double gainsyy[] = {gain1yy, gain2yy};

  for (int time = 0; time < NUM_TIME_STEPS; time++) {
    for (int visi = 0; visi < NUM_VISI/2; visi++) {
      expec_gainsxx[time*(NUM_VISI/2) + visi] = gainsxx[time];
      expec_gainsyy[time*(NUM_VISI/2) + visi] = gainsyy[time];
    }

  }

  #ifdef DOUBLE_PRECISION
    double TOL = 1e-8;
  #else
    double TOL = 1e-5;
  #endif

  //We're testing all COMPONENT types with this function, where the values of
  //the real vary slightly for baseline (I've set the major/minor to be small
  //enough to be close to 1.0 but not quite)
  for (int visi = 0; visi < NUM_VISI; visi++) {
  //   printf("%.1f %.1f %.1f %.1f %.1f %.1f %.1f %.1f\n",
  //           visibility_set->sum_visi_XX_real[visi],
  //           visibility_set->sum_visi_XX_imag[visi],
  //           visibility_set->sum_visi_XY_real[visi],
  //           visibility_set->sum_visi_XY_imag[visi],
  //           visibility_set->sum_visi_YX_real[visi],
  //           visibility_set->sum_visi_YX_imag[visi],
  //           visibility_set->sum_visi_YY_imag[visi],
  //           visibility_set->sum_visi_YY_real[visi] );

    // printf("Expec GAIN %.16f %.16f\n", expec_gainsxx[visi]
    //                             visibility_set->sum_visi_XX_real[visi] );

    TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_gainsxx[visi],
                              visibility_set->sum_visi_XX_real[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0,
                              visibility_set->sum_visi_XX_imag[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0,
                              visibility_set->sum_visi_XY_real[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0,
                              visibility_set->sum_visi_XY_imag[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0,
                              visibility_set->sum_visi_YX_real[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0,
                              visibility_set->sum_visi_YX_imag[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_gainsyy[visi],
                              visibility_set->sum_visi_YY_real[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0,
                              visibility_set->sum_visi_YY_imag[visi]);

  }
  free(expec_gainsxx);
  free(expec_gainsyy);
}
