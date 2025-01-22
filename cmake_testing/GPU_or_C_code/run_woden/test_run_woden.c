#include <unity.h>
// #include <stdlib.h>
// #include <math.h>

// #include "constants.h"
// #include "beam_settings.h"
// #include "woden_struct_defs.h"
#include "woden.h"
#include "calculate_visibilities_common_common.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

#ifdef DOUBLE_PRECISION
  double TOL = 1e-6;
#else
  double TOL = 1e-4;
#endif



void check_visi_contents(visibility_set_t *visibility_set, int num_cross,
                         double _Complex g1x, double _Complex D1x,
                         double _Complex D1y, double _Complex g1y,
                         double _Complex g2x, double _Complex D2x,
                         double _Complex D2y, double _Complex g2y) {
  double _Complex g1x_conj, D1x_conj, D1y_conj, g1y_conj;
  double _Complex g2x_conj, D2x_conj, D2y_conj, g2y_conj;
  double _Complex xx1, xy1, yx1, yy1;
  double _Complex xx2, xy2, yx2, yy2;

  g1x_conj = conj(g1x);
  D1x_conj = conj(D1x);
  D1y_conj = conj(D1y);
  g1y_conj = conj(g1y);

  g2x_conj = conj(g2x);
  D2x_conj = conj(D2x);
  D2y_conj = conj(D2y);
  g2y_conj = conj(g2y);

  xx1 = (g1x*g1x_conj + D1x*D1x_conj);
  xy1 = (g1x*D1y_conj + D1x*g1y_conj);
  yx1 = (D1y*g1x_conj + g1y*D1x_conj);
  yy1 = (D1y*D1y_conj + g1y*g1y_conj);

  xx2 = (g2x*g2x_conj + D2x*D2x_conj);
  xy2 = (g2x*D2y_conj + D2x*g2y_conj);
  yx2 = (D2y*g2x_conj + g2y*D2x_conj);
  yy2 = (D2y*D2y_conj + g2y*g2y_conj);


  for (int visi = 0; visi < num_cross/2; visi++) {

    TEST_ASSERT_DOUBLE_WITHIN(TOL, creal(xx1), visibility_set->sum_visi_XX_real[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, cimag(xx1), visibility_set->sum_visi_XX_imag[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, creal(xy1), visibility_set->sum_visi_XY_real[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, cimag(xy1), visibility_set->sum_visi_XY_imag[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, creal(yx1), visibility_set->sum_visi_YX_real[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, cimag(yx1), visibility_set->sum_visi_YX_imag[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, creal(yy1), visibility_set->sum_visi_YY_real[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, cimag(yy1), visibility_set->sum_visi_YY_imag[visi]);

  }

  for (int visi = num_cross/2; visi < num_cross; visi++) {

    TEST_ASSERT_DOUBLE_WITHIN(TOL, creal(xx2), visibility_set->sum_visi_XX_real[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, cimag(xx2), visibility_set->sum_visi_XX_imag[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, creal(xy2), visibility_set->sum_visi_XY_real[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, cimag(xy2), visibility_set->sum_visi_XY_imag[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, creal(yx2), visibility_set->sum_visi_YX_real[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, cimag(yx2), visibility_set->sum_visi_YX_imag[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, creal(yy2), visibility_set->sum_visi_YY_real[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, cimag(yy2), visibility_set->sum_visi_YY_imag[visi]);

  }
}

/*
Test `fill_primary_beam_settings` when `beamtype = MWA_ANALY`
*/
void test_run_woden(void) {

  //Setup array layout as defined in calculate_visibilities_common_common.c
  array_layout_t *array_layout = malloc(sizeof(array_layout_t));
  fill_array_layout(array_layout);


  double ra0 = 0.0;
  double dec0 = -0.46606083776035967;
  int num_bands = 3;

  //Setup woden_settingsas defined in calculate_visibilities_common_common.c,
  //then add some arguments
  woden_settings_t *woden_settings = make_woden_settings(ra0, dec0);
  woden_settings->do_gpu = 1;
  woden_settings->verbose = 1;
  woden_settings->num_bands = num_bands;
  woden_settings->off_cardinal_dipoles = 0;

  //To do the most code coverage, we need to do the MWA FEE beam
  //If the user doesn't have the MWA_FEE_HDF5 env var set, then
  //just default to NO_BEAM


  char* mwa_fee_hdf5 = getenv("MWA_FEE_HDF5");

  if (mwa_fee_hdf5) {
    printf("MWA_FEE_HDF5: %s\n", mwa_fee_hdf5 );
    woden_settings->beamtype = FEE_BEAM;
    woden_settings->hdf5_beam_path = mwa_fee_hdf5;
    // int delays[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    woden_settings->FEE_ideal_delays = malloc(16*sizeof(int));
    for (int i = 0; i < 16; i++)
    {
      woden_settings->FEE_ideal_delays[i] = 0;
    }

    
  } 
  else{
    printf("MWA_FEE_HDF5 not found - running test_woden with NO_BEAM\n");
    woden_settings->beamtype = NO_BEAM;
  }

  woden_settings->band_nums = malloc(num_bands*sizeof(int));
  for (int i = 0; i < num_bands; i++)
  {
    woden_settings->band_nums[i] = i + 1;
  }

  woden_settings->base_low_freq = BASE_BAND_FREQ;

  //Not doing shapelets so shapelet basis function array can be a NULL
  user_precision_t *sbf = NULL;

   
  visibility_set_t *visibility_sets = malloc(num_bands*sizeof(visibility_set_t));
  visibility_set_t *visibility_set;

  for (int i = 0; i < num_bands; i++)
  {
    
    visibility_set = setup_visibility_set(woden_settings->num_cross);
    visibility_sets[i] = *visibility_set;
  }

  //Just do a single point source at the phase centre
  int num_sources = 1;
  int n_points = 1;
  int n_comps = 1;

  source_catalogue_t *cropped_sky_models = malloc(sizeof(cropped_sky_models));
  cropped_sky_models->num_sources = num_sources;
  cropped_sky_models->num_shapelets = 0;
  cropped_sky_models->sources = malloc(num_sources*sizeof(source_t));

  cropped_sky_models->sources[0].n_points = n_points;
  cropped_sky_models->sources[0].n_gauss = 0;
  cropped_sky_models->sources[0].n_shapes = 0;
  cropped_sky_models->sources[0].n_shape_coeffs = 0;

  cropped_sky_models->sources[0].n_point_powers = n_points;
  cropped_sky_models->sources[0].n_point_curves = 0;
  cropped_sky_models->sources[0].n_point_lists = 0;
  cropped_sky_models->sources[0].n_gauss_powers = 0;
  cropped_sky_models->sources[0].n_gauss_curves = 0;
  cropped_sky_models->sources[0].n_gauss_lists = 0;
  cropped_sky_models->sources[0].n_shape_powers = 0;
  cropped_sky_models->sources[0].n_shape_curves = 0;
  cropped_sky_models->sources[0].n_shape_lists = 0;

  components_t *comps = &cropped_sky_models->sources[0].point_components;

  comps->ras = malloc(n_comps*sizeof(double));
  comps->decs = malloc(n_comps*sizeof(double));
  comps->azs = malloc(n_comps*NUM_TIME_STEPS*sizeof(user_precision_t));
  comps->zas = malloc(n_comps*NUM_TIME_STEPS*sizeof(user_precision_t));
  comps->power_ref_freqs = malloc(sizeof(double));
  comps->power_ref_stokesI = malloc(sizeof(user_precision_t));
  comps->power_SIs = malloc(sizeof(user_precision_t));
  comps->power_comp_inds = malloc(sizeof(int));
  comps->n_stokesV_power = 0;
  comps->n_stokesV_curve = 0;
  comps->n_stokesV_list = 0;
  comps->n_linpol_power = 0;
  comps->n_linpol_curve = 0;
  comps->n_linpol_list = 0;
  comps->n_linpol_p_list = 0;
  comps->n_linpol_angles = 0;
  comps->n_stokesV_pol_frac = 0;
  comps->n_linpol_pol_frac = 0;

  comps->ras[0] = ra0;
  comps->decs[0] = dec0;
  comps->power_ref_freqs[0] = 200e+6;
  comps->power_ref_stokesI[0] = 1.0;
  comps->power_SIs[0] = 0.0;
  comps->power_comp_inds[0] = 0;
  
  user_precision_t azs[] = {0.0, 4.528359553989764};
  user_precision_t zas[] = {0.0, 0.6978088917603547};

  for (int time = 0; time < NUM_TIME_STEPS; time++) {

    comps->azs[time] = azs[time];
    comps->zas[time] = zas[time];
  }

  
  
  //No matter what, always have one POL_FRACTION source for both pols
  comps->do_QUV = 0;

  run_woden(woden_settings, visibility_sets,
             cropped_sky_models, array_layout,
             sbf);

  //Run uvw test defined in calculate_visibilities_common_common.c
  for (int i = 0; i < num_bands; i++)
  {
    // visibility_set_t *visibility_set = &visibility_sets[i];
    test_uvw(&visibility_sets[i], woden_settings->lsts,
              ra0, dec0, woden_settings);
  }

  double _Complex g1x, D1x, D1y, g1y;
  double _Complex g2x, D2x, D2y, g2y;
  
  if (mwa_fee_hdf5) {

    g1x = 0.87750363 -0.47957003*I;
    D1x = -0.00019401 -0.00003449*I;
    D1y = -0.00020211 -0.00002737*I;
    g1y = 0.87775218 -0.47911486*I;
    g2x = -0.07190398 +0.04164878*I;
    D2x = 0.00003236 -0.00019115*I;
    D2y = 0.00500040 -0.00240635*I;
    g2y = -0.05687106 +0.02499121*I;

    check_visi_contents(&visibility_sets[0], woden_settings->num_cross,
                         g1x, D1x, D1y, g1y, g2x, D2x, D2y, g2y);

    g1x = 0.89373189 -0.44860148*I;
    D1x = -0.00020943 -0.00003043*I;
    D1y = -0.00021026 -0.00003153*I;
    g1y = 0.89396733 -0.44813216*I;
    g2x = -0.08054522 +0.04309910*I;
    D2x = -0.00000170 -0.00014954*I;
    D2y = 0.00559733 -0.00245659*I;
    g2y = -0.06276022 +0.02624865*I;

    check_visi_contents(&visibility_sets[1], woden_settings->num_cross,
                         g1x, D1x, D1y, g1y, g2x, D2x, D2y, g2y);

    g1x = 0.90929735 -0.41614705*I;
    D1x = -0.00020053 -0.00004292*I;
    D1y = -0.00020929 -0.00003362*I;
    g1y = 0.90951550 -0.41566998*I;
    g2x = -0.08919068 +0.04381648*I;
    D2x = -0.00003944 -0.00011095*I;
    D2y = 0.00619277 -0.00245741*I;
    g2y = -0.06864354 +0.02692343*I;

    check_visi_contents(&visibility_sets[2], woden_settings->num_cross,
                         g1x, D1x, D1y, g1y, g2x, D2x, D2y, g2y);

  } else {
    g1x = 1.0 +0.0*I;
    D1x = 0.0 +0.0*I;
    D1y = 0.0 +0.0*I;
    g1y = 1.0 +0.0*I;

    g2x = 1.0 +0.0*I;
    D2x = 0.0 +0.0*I;
    D2y = 0.0 +0.0*I;
    g2y = 1.0 +0.0*I;

    for (int i = 0; i < num_bands; i++)
    {
      check_visi_contents(&visibility_sets[i], woden_settings->num_cross,
                         g1x, D1x, D1y, g1y, g2x, D2x, D2y, g2y);
    }
  }

  free(comps->ras);
  free(comps->decs);
  free(comps->azs);
  free(comps->zas);
  free(comps->power_ref_freqs);
  free(comps->power_ref_stokesI);
  free(comps->power_SIs);
  free(comps->power_comp_inds);
  free(cropped_sky_models->sources);
  free(cropped_sky_models);


  for (int i = 0; i < num_bands; i++)
  {
    free_visi_set_inputs(&visibility_sets[i]);
    free_visi_set_outputs(&visibility_sets[i]);
  }

  free(array_layout->X_diff_metres);
  free(array_layout->Y_diff_metres);
  free(array_layout->Z_diff_metres);
  // free(array_layout);

  free(woden_settings->mwa_dipole_amps);
  free(woden_settings->FEE_ideal_delays);
  free(woden_settings->band_nums);
  free(woden_settings);


}


//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_run_woden);

    return UNITY_END();
}
