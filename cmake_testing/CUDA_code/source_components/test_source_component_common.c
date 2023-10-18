#include "test_source_component_common.h"
#include "common_testing_functions.h"
#include <mwa_hyperbeam.h>

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

// void sincos(user_precision_t x, user_precision_t *sin, user_precision_t *cos);

//External CUDA code we're linking in
extern void test_source_component_common(int num_of_each_flux_type,
           components_t components,
           double *freqs, woden_settings_t *woden_settings,
           beam_settings_t *beam_settings,
           user_precision_complex_t *gxs, user_precision_complex_t *Dxs,
           user_precision_complex_t *Dys, user_precision_complex_t *gys,
           user_precision_t *extrap_flux_I, user_precision_t *extrap_flux_Q,
           user_precision_t *extrap_flux_U, user_precision_t *extrap_flux_V,
           double *ls, double *ms, double *ns,
           e_component_type comptype);

double TOL;

/*
Test that l,m,n and beam values are calculated correctly by `source_component_common`
for a constant Dec=0, latitude=0, and all beam types
There are other tests for the l,m,n coords and beam functions, so no need to
test millions of scenarios here, so stick with Dec=0
*/
void test_source_component_common_ConstantDecChooseBeams(int beamtype, char* mwa_fee_hdf5,
                                                         e_component_type comptype) {

  //Set up some test condition inputs
  int num_times = 3;
  int num_freqs = 2;

  int num_components = num_powers + num_curves + num_lists;

  int num_beam_values = num_times*num_freqs*num_components;

  user_precision_t *zeroes = calloc(num_components, sizeof(user_precision_t));

  //Keep RA between 0 and 2*pi here but enter RAs that should return
  //negative l values
  double ras[9] = {(3*M_PI)/2, (5*M_PI)/3, (7*M_PI)/4, (11*M_PI)/6,
                   0.0, M_PI/6, M_PI/4, M_PI/3, M_PI/2};

  double ra0 = 0.0*DD2R;
  double dec0 = 0.0*DD2R;

  double *decs = malloc(num_components*sizeof(double));
  for (int i = 0; i < num_components; i++) {
    decs[i] = dec0;
  }

  //Get the settings into a woden_settings_t struct
  woden_settings_t *woden_settings = malloc(sizeof(woden_settings_t));
  woden_settings->num_freqs = num_freqs;
  woden_settings->num_time_steps = num_times;
  woden_settings->ra0 = ra0;
  woden_settings->sdec0 = sin(dec0);
  woden_settings->cdec0 = cos(dec0);
  woden_settings->latitude = -0.46606083776035967;

  woden_settings->latitudes = malloc(num_times*sizeof(double));
  for (int i = 0; i < num_times; i++)
  {
    woden_settings->latitudes[i] = -0.46606083776035967;
  }
  

  woden_settings->beamtype = beamtype;

  beam_settings_t *beam_settings = malloc(sizeof(beam_settings_t));
  beam_settings->beamtype = beamtype;

  double freqs[] = {100e+6, 200e+6};

  double *beam_has = malloc(num_components*num_times*sizeof(double));
  double *beam_decs = malloc(num_components*num_times*sizeof(double));

  //Make ha/dec coords if using the Gaussian beam
  //Need ha/decs for the analytic MWA beam as well
  if (beamtype == GAUSS_BEAM || beamtype == MWA_ANALY) {
    //Stick the Gaussian beam pointed at zenith
    //We're testing at latitude=zero
    beam_settings->gauss_ha = 0.0;
    beam_settings->gauss_sdec = 0.0;
    beam_settings->gauss_cdec = 1.0;
    beam_settings->beam_FWHM_rad = 80.0*DD2R;
    beam_settings->beam_ref_freq = 150e+6;

    for (int component = 0; component < num_components; component++) {
      for (int time_ind = 0; time_ind < num_times; time_ind++) {

        beam_has[component*num_times + time_ind] = ras[component];
        beam_decs[component*num_times + time_ind] = decs[component];

      }
    }
  }

  //If FEE_BEAM, call the C code to interrogate the hdf5 file and set beam
  //things up
  if (beamtype == FEE_BEAM || beamtype == FEE_BEAM_INTERP) {

    int32_t status =  new_fee_beam(mwa_fee_hdf5, &beam_settings->fee_beam);

    TEST_ASSERT_EQUAL(status, 0);

    uint32_t num_freqs_hyper;
    uint32_t *freqs_hz;
    if (beamtype == FEE_BEAM) {
      freqs_hz = malloc(2*sizeof(uint32_t));
      freqs_hz[0] = 150e+6;
      freqs_hz[1] = 150e+6;
      num_freqs_hyper = 2;
    } else {
      freqs_hz = malloc(2*sizeof(uint32_t));
      freqs_hz[0] = 100e+6;
      freqs_hz[1] = 200e+6;
      num_freqs_hyper = 2;
    }

    beam_settings->hyper_delays = (uint32_t*)malloc(16*sizeof(uint32_t));

    for (int delay = 0; delay < 16; delay++) {
      beam_settings->hyper_delays[delay] = 0;
    }

    double amps[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    uint32_t num_tiles = 1;
    uint32_t num_amps = 16;
    uint8_t norm_to_zenith = 1;

    status = new_gpu_fee_beam(beam_settings->fee_beam,
                               freqs_hz,
                               beam_settings->hyper_delays,
                               amps,
                               num_freqs_hyper,
                               num_tiles,
                               num_amps,
                               norm_to_zenith,
                               &beam_settings->cuda_fee_beam);

    TEST_ASSERT_EQUAL(status, 0);

  }

  else if (beamtype == MWA_ANALY) {
    //Zenith pointing is your friend
    for(int i=0;i<16;i++) {
        woden_settings->FEE_ideal_delays[i] = 0.0;
    }
  }

  // //Output arrays
  user_precision_complex_t *primay_beam_J00 = calloc(num_beam_values, sizeof(user_precision_complex_t));
  user_precision_complex_t *primay_beam_J01 = calloc(num_beam_values, sizeof(user_precision_complex_t));
  user_precision_complex_t *primay_beam_J10 = calloc(num_beam_values, sizeof(user_precision_complex_t));
  user_precision_complex_t *primay_beam_J11 = calloc(num_beam_values, sizeof(user_precision_complex_t));

  double *ls = malloc(num_components*sizeof(double));
  double *ms = malloc(num_components*sizeof(double));
  double *ns = malloc(num_components*sizeof(double));

  //Space for outputs
  user_precision_t *extrap_flux_I = malloc(num_freqs*num_components*sizeof(user_precision_t));
  user_precision_t *extrap_flux_Q = malloc(num_freqs*num_components*sizeof(user_precision_t));
  user_precision_t *extrap_flux_U = malloc(num_freqs*num_components*sizeof(user_precision_t));
  user_precision_t *extrap_flux_V = malloc(num_freqs*num_components*sizeof(user_precision_t));

  //Set up the components with our values
  components_t components;

  components.ras = ras;
  components.decs = decs;
  components.azs = azs;
  components.zas = zas;
  components.beam_has = beam_has;
  components.beam_decs = beam_decs;

  int *power_comp_inds = malloc(num_powers*sizeof(int));
  int *curve_comp_inds = malloc(num_curves*sizeof(int));
  int *list_comp_inds = malloc(num_lists*sizeof(int));

  for (int i = 0; i < num_powers; i++) {
    power_comp_inds[i] = i;
  }

  for (int i = 0; i < num_curves; i++) {
    curve_comp_inds[i] = num_powers + i;
  }

  for (int i = 0; i < num_lists; i++) {
    list_comp_inds[i] = num_powers + num_curves + i;
  }


  components.power_ref_freqs = ref_freqs;
  components.power_ref_stokesI = ref_stokesI;
  components.power_ref_stokesQ = ref_stokesQ;
  components.power_ref_stokesU = ref_stokesU;
  components.power_ref_stokesV = ref_stokesV;
  components.power_SIs = ref_power_SIs;
  components.power_comp_inds = power_comp_inds;

  components.curve_ref_freqs = ref_freqs;
  components.curve_ref_stokesI = ref_stokesI;
  components.curve_ref_stokesQ = ref_stokesQ;
  components.curve_ref_stokesU = ref_stokesU;
  components.curve_ref_stokesV = ref_stokesV;
  components.curve_SIs = ref_curve_SIs;
  components.curve_qs = ref_qs;
  components.curve_comp_inds = curve_comp_inds;

  components.list_freqs = list_freqs;
  components.list_stokesI = list_stokesI;
  components.list_stokesQ = list_stokesQ;
  components.list_stokesU = list_stokesU;
  components.list_stokesV = list_stokesV;
  components.num_list_values = num_list_values;
  components.list_start_indexes = list_start_indexes;
  components.list_comp_inds = list_comp_inds;

  components.total_num_flux_entires = 0;

  for (int i = 0; i < num_lists; i++) {
    components.total_num_flux_entires += num_list_values[i];
  }

  //Don't need to test values are copied as tested elsewhere, but the function
  //that copies components from CPU to GPU needs these extra fields defined
  //for GAUSSIAN/SHAPELET components
  if (comptype == GAUSSIAN || comptype == SHAPELET) {
    components.pas = zeroes;
    components.majors = zeroes;
    components.minors = zeroes;
  }

  if (comptype == SHAPELET) {
    components.shape_coeffs = zeroes;
    components.n1s = zeroes;
    components.n2s = zeroes;
    components.param_indexes = zeroes;
  }

  components.num_primarybeam_values = num_components*woden_settings->num_freqs*woden_settings->num_time_steps;

  //Run the CUDA code
  test_source_component_common(num_powers, components, freqs,
           woden_settings, beam_settings,
           primay_beam_J00, primay_beam_J01,
           primay_beam_J10, primay_beam_J11,
           extrap_flux_I, extrap_flux_Q, extrap_flux_U, extrap_flux_V,
           ls, ms, ns, comptype);

  double l_expected[9] = {-1.0, -sqrt(3)/2.0, -sqrt(2)/2.0, -0.5,
                          0.0, 0.5, sqrt(2)/2.0, sqrt(3)/2.0, 1.0};
  double n_expected[9] = {0.0, 0.5, sqrt(2)/2.0, sqrt(3)/2.0,
                         1.0, sqrt(3)/2.0, sqrt(2)/2.0, 0.5, 0.0};

  //Check the l values match expectations
  //Both FLOAT and DOUBLE use 64 bit here so same tolerance
  TOL = 1e-15;

  for (int i = 0; i < num_components; i++) {
    TEST_ASSERT_DOUBLE_WITHIN(TOL, l_expected[i], ls[i]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, ms[i]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, n_expected[i], ns[i]);
  }

  //Depending on beamtype, check results match expectations

  //For the Gaussian beam, the way the l,m coords are set up measns we can
  //analytically predict the values, so go through results and calcluate
  //expected values, and compare to what we got
  int beam_ind = 0;
  if (beamtype == GAUSS_BEAM) {

    #ifdef DOUBLE_PRECISION
      TOL = 1e-12;
    #else
      TOL = 1e-7;
    #endif

    double fwhm_lm = sin(beam_settings->beam_FWHM_rad);
    double beam_ref_freq = 150e+6;

    for (int time_ind = 0; time_ind < num_times; time_ind++) {
      for (int freq_ind = 0; freq_ind < num_freqs; freq_ind++) {
        for (int comp_ind = 0; comp_ind < num_components; comp_ind++) {

          double std = (fwhm_lm / FWHM_FACTOR) * (beam_ref_freq / freqs[freq_ind]);
          double exp_inside = (double)ls[comp_ind] / std;
          double estimate = exp(-0.5*exp_inside*exp_inside);

          //Only real gains have value in GAUSSIAN beam model
          TEST_ASSERT_DOUBLE_WITHIN(TOL, estimate,
                                   creal(primay_beam_J00[beam_ind]));

          TEST_ASSERT_DOUBLE_WITHIN(TOL, estimate,
                                   creal(primay_beam_J11[beam_ind]));

          //Everything else should be zero
          TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(primay_beam_J10[beam_ind]));
          TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(primay_beam_J11[beam_ind]));
          TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(primay_beam_J00[beam_ind]));
          TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, creal(primay_beam_J01[beam_ind]));
          TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(primay_beam_J01[beam_ind]));
          TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, creal(primay_beam_J10[beam_ind]));

          beam_ind += 1;

        }
      }
    }
  }
  else if (beamtype == ANALY_DIPOLE) {

    #ifdef DOUBLE_PRECISION
      TOL = 1e-7;
    #else
      TOL = 1e-6;
    #endif

    for (int output = 0; output < num_beam_values; output++) {

      TEST_ASSERT_DOUBLE_WITHIN(TOL, analy_expec_J00[output], creal(primay_beam_J00[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, analy_expec_J11[output], creal(primay_beam_J11[output]));

      TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(primay_beam_J00[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, creal(primay_beam_J01[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(primay_beam_J01[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, creal(primay_beam_J10[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(primay_beam_J10[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(primay_beam_J11[output]));
    }
  }
  else if (beamtype == FEE_BEAM) {

    #ifdef DOUBLE_PRECISION
      TOL = 1e-7;
    #else
      TOL = 1e-7;
    #endif

    for (int output = 0; output < num_beam_values; output++) {
      // printf("%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n", creal(primay_beam_J00[output]), cimag(primay_beam_J00[output]),
      //         creal(primay_beam_J01[output]), cimag(primay_beam_J01[output]),
      //         creal(primay_beam_J10[output]), cimag(primay_beam_J10[output]),
      //         creal(primay_beam_J11[output]), cimag(primay_beam_J11[output]) );


      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_J00_re[output], creal(primay_beam_J00[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_J00_im[output], cimag(primay_beam_J00[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_J01_re[output], creal(primay_beam_J01[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_J01_im[output], cimag(primay_beam_J01[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_J10_re[output], creal(primay_beam_J10[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_J10_im[output], cimag(primay_beam_J10[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_J11_re[output], creal(primay_beam_J11[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_J11_im[output], cimag(primay_beam_J11[output]));
    }

    free_fee_beam(beam_settings->fee_beam);
    free_gpu_fee_beam(beam_settings->cuda_fee_beam);

  }
  else if (beamtype == FEE_BEAM_INTERP) {

    #ifdef DOUBLE_PRECISION
      TOL = 1e-7;
    #else
      TOL = 1e-7;
    #endif

    for (int output = 0; output < num_beam_values; output++) {

      // printf("%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n", creal(primay_beam_J00[output]), cimag(primay_beam_J00[output]),
      //         creal(primay_beam_J01[output]), cimag(primay_beam_J01[output]),
      //         creal(primay_beam_J10[output]), cimag(primay_beam_J10[output]),
      //         creal(primay_beam_J11[output]), cimag(primay_beam_J11[output]) );

      // printf("%.8f %.8f \n", creal(primay_beam_J00[output]),
      //                        fee_expec_interp_J00_re[output]  );


      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_interp_J00_re[output], creal(primay_beam_J00[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_interp_J00_im[output], cimag(primay_beam_J00[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_interp_J01_re[output], creal(primay_beam_J01[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_interp_J01_im[output], cimag(primay_beam_J01[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_interp_J10_re[output], creal(primay_beam_J10[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_interp_J10_im[output], cimag(primay_beam_J10[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_interp_J11_re[output], creal(primay_beam_J11[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_interp_J11_im[output], cimag(primay_beam_J11[output]));
    }

    free_fee_beam(beam_settings->fee_beam);
    free_gpu_fee_beam(beam_settings->cuda_fee_beam);
  }
  else if (beamtype == MWA_ANALY) {

    #ifdef DOUBLE_PRECISION
      TOL = 1e-12;
    #else
      TOL = 1e-7;
    #endif

    //All imaginary should be zero.Turns out for these ra,dec coords, the
    //Dy leakages are essentially zero as well, test a few values equal zero
    for (int output = 0; output < num_beam_values; output++) {

        TEST_ASSERT_DOUBLE_WITHIN(TOL, MWA_analy_expec_J00_re[output],
                                  creal(primay_beam_J00[output]));

        TEST_ASSERT_DOUBLE_WITHIN(TOL, MWA_analy_expec_J01_re[output],
                                  creal(primay_beam_J01[output]));

        TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, creal(primay_beam_J10[output]));

        TEST_ASSERT_DOUBLE_WITHIN(TOL, MWA_analy_expec_J11_re[output],
                                  creal(primay_beam_J11[output]));

        TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(primay_beam_J00[output]));
        TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(primay_beam_J01[output]));
        TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(primay_beam_J10[output]));
        TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(primay_beam_J11[output]));

        // printf("%.12f %.12f %.12f %.12f\n",creal(primay_beam_J00[output]),
        // creal(primay_beam_J01[output]),
        // creal(primay_beam_J10[output]),
        // creal(primay_beam_J11[output]) );

    }
  }


  else if (beamtype == NO_BEAM) {
    //Don't need to check beam values for NO_BEAM, these values
    //The function `source_components::get_beam_gains` which is called later
    //by the COMPONENT kernels assigns gains on one and leakage of zero
    //in the NO_BEAM case.
    ;
  }

  //Now check the fluxes were extrapolated correctly
  //Make some expected value arrays
  double *expec_flux_I = malloc(num_freqs*num_components*sizeof(double));
  double *expec_flux_Q = malloc(num_freqs*num_components*sizeof(double));
  double *expec_flux_U = malloc(num_freqs*num_components*sizeof(double));
  double *expec_flux_V = malloc(num_freqs*num_components*sizeof(double));

  CPU_extrapolate_fluxes_in_components(&components, num_powers, num_curves, num_lists,
                        freqs, num_freqs,
                        expec_flux_I, expec_flux_Q, expec_flux_U, expec_flux_V);

  #ifdef DOUBLE_PRECISION
    TOL = 1e-12;
  #else
    TOL = 1e-4;
  #endif

  for (int i = 0; i < num_freqs*(num_powers + num_curves + num_lists); i++) {
    //Check the two are within tolerace
    // printf("%d %.3f %.3f\n",i, expec_flux_I[i], extrap_flux_I[i] );
    TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_flux_I[i], extrap_flux_I[i]);
    //in the future this should be testable, for only doing Stokes I so
    //lock to zero
    // TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_flux_Q[i], extrap_flux_Q[i]);
    // TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_flux_U[i], extrap_flux_U[i]);
    // TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_flux_V[i], extrap_flux_V[i]);

    TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, extrap_flux_Q[i]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, extrap_flux_U[i]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, extrap_flux_V[i]);

  }

  //Be free my beauties
  free(zeroes);
  free(decs);
  free(primay_beam_J00);
  free(primay_beam_J01);
  free(primay_beam_J10);
  free(primay_beam_J11);
  free(ls);
  free(ms);
  free(ns);
  free(beam_has);
  free(beam_decs);
  free(power_comp_inds);
  free(curve_comp_inds);
  free(list_comp_inds);
  free(woden_settings->latitudes);

}

/*
This test checks source_component_common with beamtype=FEE_BEAM, for a given
comptype
*/
void test_source_component_common_ConstantDecFEEBeam(e_component_type comptype) {
  //Look for environment variable MWA_FEE_HDF5, which should point
  //towards mwa_full_embedded_element_pattern.h5. If we can't find it,
  //can't run the tests, so don't
  char* mwa_fee_hdf5 = getenv("MWA_FEE_HDF5");

  if (mwa_fee_hdf5) {
    printf("MWA_FEE_HDF5: %s\n", mwa_fee_hdf5 );
    test_source_component_common_ConstantDecChooseBeams(FEE_BEAM, mwa_fee_hdf5, comptype);
  }
  else {
    printf("MWA_FEE_HDF5 not found - not running MWA_FEE beam test\n");
  }
}

void test_source_component_common_ConstantDecFEEBeamPoint(void) {
  test_source_component_common_ConstantDecFEEBeam(POINT);
}
void test_source_component_common_ConstantDecFEEBeamGauss(void) {
  test_source_component_common_ConstantDecFEEBeam(GAUSSIAN);
}
void test_source_component_common_ConstantDecFEEBeamShapelet(void) {
  test_source_component_common_ConstantDecFEEBeam(SHAPELET);
}

/*
This test checks source_component_common with beamtype=ANALY_DIPOLE
*/
void test_source_component_common_ConstantDecAnalyBeam(e_component_type comptype) {
  test_source_component_common_ConstantDecChooseBeams(ANALY_DIPOLE, " ", comptype);
}

void test_source_component_common_ConstantDecAnalyBeamPoint(void){
  test_source_component_common_ConstantDecAnalyBeam(POINT);
}
void test_source_component_common_ConstantDecAnalyBeamGaussian(void){
  test_source_component_common_ConstantDecAnalyBeam(GAUSSIAN);
}
void test_source_component_common_ConstantDecAnalyBeamShapelet(void){
  test_source_component_common_ConstantDecAnalyBeam(SHAPELET);
}

/*
This test checks source_component_common with beamtype=GAUSS_BEAM
*/
void test_source_component_common_ConstantDecGaussBeam(e_component_type comptype) {
  test_source_component_common_ConstantDecChooseBeams(GAUSS_BEAM, " ", comptype);
}

void test_source_component_common_ConstantDecGaussBeamPoint(void){
  test_source_component_common_ConstantDecGaussBeam(POINT);
}
void test_source_component_common_ConstantDecGaussBeamGaussian(void){
  test_source_component_common_ConstantDecGaussBeam(GAUSSIAN);
}
void test_source_component_common_ConstantDecGaussBeamShapelet(void){
  test_source_component_common_ConstantDecGaussBeam(SHAPELET);
}

/*
This test checks source_component_common with beamtype=NO_BEAM
*/
void test_source_component_common_ConstantDecNoBeam(e_component_type comptype) {
  test_source_component_common_ConstantDecChooseBeams(NO_BEAM, " ", comptype);
}

void test_source_component_common_ConstantDecNoBeamPoint(void){
  test_source_component_common_ConstantDecNoBeam(POINT);
}
void test_source_component_common_ConstantDecNoBeamGaussian(void){
  test_source_component_common_ConstantDecNoBeam(GAUSSIAN);
}
void test_source_component_common_ConstantDecNoBeamShapelet(void){
  test_source_component_common_ConstantDecNoBeam(SHAPELET);
}

/*
This test checks source_component_common with beamtype=FEE_BEAM_INTERP
*/
void test_source_component_common_ConstantDecFEEBeamInterp(e_component_type comptype) {
  //Look for environment variable MWA_FEE_HDF5, which should point
  //towards mwa_full_embedded_element_pattern.h5. If we can't find it,
  //can't run the tests, so don't
  char* mwa_fee_hdf5 = getenv("MWA_FEE_HDF5_INTERP");

  if (mwa_fee_hdf5) {
    printf("MWA_FEE_HDF5_INTERP: %s\n", mwa_fee_hdf5 );
    test_source_component_common_ConstantDecChooseBeams(FEE_BEAM_INTERP, mwa_fee_hdf5, comptype);
  }
  else {
    printf("MWA_FEE_HDF5_INTERP not found - not running MWA_FEE beam test\n");
  }
}

void test_source_component_common_ConstantDecFEEBeamInterpPoint(void){
  test_source_component_common_ConstantDecFEEBeamInterp(POINT);
}
void test_source_component_common_ConstantDecFEEBeamInterpGaussian(void){
  test_source_component_common_ConstantDecFEEBeamInterp(GAUSSIAN);
}
void test_source_component_common_ConstantDecFEEBeamInterpShapelet(void){
  test_source_component_common_ConstantDecFEEBeamInterp(SHAPELET);
}

/*
This test checks source_component_common with beamtype=MWA_ANALY
*/
void test_source_component_common_ConstantDecMWAAnaly(e_component_type comptype) {
  test_source_component_common_ConstantDecChooseBeams(MWA_ANALY, " ", comptype);
}

void test_source_component_common_ConstantDecMWAAnalyPoint(void){
  test_source_component_common_ConstantDecMWAAnaly(POINT);
}
void test_source_component_common_ConstantDecMWAAnalyGaussian(void){
  test_source_component_common_ConstantDecMWAAnaly(GAUSSIAN);
}
void test_source_component_common_ConstantDecMWAAnalyShapelet(void){
  test_source_component_common_ConstantDecMWAAnaly(SHAPELET);
}

//Run the test with unity
int main(void)
{
    UNITY_BEGIN();

    // RUN_TEST(test_source_component_common_ConstantDecNoBeamPoint);
    // RUN_TEST(test_source_component_common_ConstantDecNoBeamGaussian);
    // RUN_TEST(test_source_component_common_ConstantDecNoBeamShapelet);

    // RUN_TEST(test_source_component_common_ConstantDecFEEBeamPoint);
    // RUN_TEST(test_source_component_common_ConstantDecFEEBeamGauss);
    // RUN_TEST(test_source_component_common_ConstantDecFEEBeamShapelet);

    // RUN_TEST(test_source_component_common_ConstantDecAnalyBeamPoint);
    // RUN_TEST(test_source_component_common_ConstantDecAnalyBeamGaussian);
    // RUN_TEST(test_source_component_common_ConstantDecAnalyBeamShapelet);

    // RUN_TEST(test_source_component_common_ConstantDecGaussBeamPoint);
    // RUN_TEST(test_source_component_common_ConstantDecGaussBeamGaussian);
    // RUN_TEST(test_source_component_common_ConstantDecGaussBeamShapelet);

    RUN_TEST(test_source_component_common_ConstantDecFEEBeamInterpPoint);
    // RUN_TEST(test_source_component_common_ConstantDecFEEBeamInterpGaussian);
    // RUN_TEST(test_source_component_common_ConstantDecFEEBeamInterpShapelet);

    // RUN_TEST(test_source_component_common_ConstantDecMWAAnalyPoint);
    // RUN_TEST(test_source_component_common_ConstantDecMWAAnalyGaussian);
    // RUN_TEST(test_source_component_common_ConstantDecMWAAnalyShapelet);

    return UNITY_END();
}
