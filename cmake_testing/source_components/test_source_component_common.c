#include "test_source_component_common.h"
#include <mwa_hyperbeam.h>

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

// void sincos(user_precision_t x, user_precision_t *sin, user_precision_t *cos);

//External CUDA code we're linking in
extern void test_source_component_common(int num_components,
           int num_shape_coeffs, components_t components,
           components_t d_components,
           double *freqs, woden_settings_t *woden_settings,
           beam_settings_t *beam_settings,
           user_precision_complex_t *gxs, user_precision_complex_t *Dxs,
           user_precision_complex_t *Dys, user_precision_complex_t *gys,
           double *ls, double *ms, double *ns);

double TOL;

/*
Test that l,m,n and beam values are calculated correctly by `source_component_common`
for a constant Dec=0, latitude=0, and all beam types
There are other tests for the l,m,n coords and beam functions, so no need to
test millions of scenarios here, so stick with Dec=0
*/
void test_source_component_common_ConstantDecChooseBeams(int beamtype, char* mwa_fee_hdf5) {

  //Set up some test condition inputs
  int num_times = 3;
  int num_freqs = 2;
  int num_components = 9;
  // int num_baselines = 5;

  int num_beam_values = num_times*num_freqs*num_components;

  user_precision_t *zeroes = calloc(num_components, sizeof(user_precision_t));

  user_precision_t *ref_stokesI = malloc(num_components*sizeof(user_precision_t));
  user_precision_t *ref_stokesQ = malloc(num_components*sizeof(user_precision_t));
  user_precision_t *ref_stokesU = malloc(num_components*sizeof(user_precision_t));
  user_precision_t *ref_stokesV = malloc(num_components*sizeof(user_precision_t));
  user_precision_t *SIs = malloc(num_components*sizeof(user_precision_t));

  ref_stokesQ = zeroes;
  ref_stokesU = zeroes;
  ref_stokesV = zeroes;
  SIs = zeroes;

  //Keep RA between 0 and 2*pi here but enter RAs that should return
  //negative l values
  double ras[9] = {(3*M_PI)/2, (5*M_PI)/3, (7*M_PI)/4, (11*M_PI)/6,
                   0.0, M_PI/6, M_PI/4, M_PI/3, M_PI/2};

  double ra0 = 0.0*DD2R;
  double dec0 = 0.0*DD2R;

  double *decs = malloc(num_components*sizeof(double));
  double *ref_freqs = malloc(num_components*sizeof(double));

  for (int i = 0; i < num_components; i++) {
    decs[i] = dec0;
    ref_freqs[i] = 150e+6;
    ref_stokesI[i] = 1.0;
  }

  //Get the settings into a woden_settings_t struct
  woden_settings_t *woden_settings = malloc(sizeof(woden_settings_t));
  woden_settings->num_freqs = num_freqs;
  woden_settings->num_time_steps = num_times;
  woden_settings->ra0 = ra0;
  woden_settings->sdec0 = sin(dec0);
  woden_settings->cdec0 = cos(dec0);
  woden_settings->latitude = MWA_LAT_RAD;

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

    int32_t status =  new_fee_beam(mwa_fee_hdf5, &beam_settings->fee_beam,
                               beam_settings->hyper_error_str);

    TEST_ASSERT_EQUAL(status, 0);

    uint32_t num_freqs_hyper;
    uint32_t *freqs_hz;
    if (beamtype == FEE_BEAM) {
      freqs_hz = malloc(sizeof(uint32_t));
      freqs_hz[0] = 150e+6;
      num_freqs_hyper = 1;
    } else {
      freqs_hz = malloc(sizeof(uint32_t));
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

    status = new_cuda_fee_beam(beam_settings->fee_beam,
                            freqs_hz,
                            beam_settings->hyper_delays,
                            amps,
                            num_freqs_hyper,
                            num_tiles,
                            num_amps,
                            norm_to_zenith,
                            &beam_settings->cuda_fee_beam,
                            beam_settings->hyper_error_str);

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

  //Set up the components with our values
  components_t components;
  components_t d_components;

  components.ras = ras;
  components.decs = decs;
  components.azs = azs;
  components.zas = zas;
  components.beam_has = beam_has;
  components.beam_decs = beam_decs;
  components.ref_freqs = ref_freqs;
  components.ref_stokesI = ref_stokesI;
  components.ref_stokesQ = ref_stokesQ;
  components.ref_stokesU = ref_stokesU;
  components.ref_stokesV = ref_stokesV;
  components.SIs = SIs;


  components.num_primarybeam_values = num_components*woden_settings->num_freqs*woden_settings->num_time_steps;

  //not testing for shapelets here
  int num_shape_coeffs = 0;

  //Run the CUDA code
  test_source_component_common(num_components, num_shape_coeffs,
           components, d_components,
           freqs, woden_settings,
           beam_settings,
           primay_beam_J00, primay_beam_J01,
           primay_beam_J10, primay_beam_J11,
           ls, ms, ns);

  double l_expected[9] = {-1.0, -sqrt(3)/2.0, -sqrt(2)/2.0, -0.5,
                          0.0, 0.5, sqrt(2)/2.0, sqrt(3)/2.0, 1.0};
  double n_expected[9] = {0.0, 0.5, sqrt(2)/2.0, sqrt(3)/2.0,
                         1.0, sqrt(3)/2.0, sqrt(2)/2.0, 0.5, 0.0};

  //Check the l values match expectations
  //Both FLOAT and DOUBLE use 64 bit here so same tolerance
  TOL = 1e-15;

  for (int i = 0; i < num_components; i++) {
    // printf("%.6e %.6e\n",l_expected[i], ls[i] );
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
    free_cuda_fee_beam(beam_settings->cuda_fee_beam);

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
    free_cuda_fee_beam(beam_settings->cuda_fee_beam);
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
    //Don't need to calculate beam values for NO_BEAM, these values
    //should still be zero, as we initialised the array with calloc
    //The function `source_components::get_beam_gains` which is called later
    //by the COMPONENT kernels assigns gains on one and leakage of zero
    //in the NO_BEAM case

    TOL = 1e-15;

    for (int output = 0; output < num_beam_values; output++) {

      TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, creal(primay_beam_J00[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(primay_beam_J00[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, creal(primay_beam_J01[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(primay_beam_J01[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, creal(primay_beam_J10[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(primay_beam_J10[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, creal(primay_beam_J11[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(primay_beam_J11[output]));

    }
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

}

/*
This test checks source_component_common with beamtype=FEE_BEAM
*/
void test_source_component_common_ConstantDecFEEBeam(void) {
  //Look for environment variable MWA_FEE_HDF5, which should point
  //towards mwa_full_embedded_element_pattern.h5. If we can't find it,
  //can't run the tests, so don't
  char* mwa_fee_hdf5 = getenv("MWA_FEE_HDF5");

  if (mwa_fee_hdf5) {
    printf("MWA_FEE_HDF5: %s", mwa_fee_hdf5 );
    test_source_component_common_ConstantDecChooseBeams(FEE_BEAM, mwa_fee_hdf5);
  }
  else {
    printf("MWA_FEE_HDF5 not found - not running MWA_FEE beam test");
  }
}

/*
This test checks source_component_common with beamtype=ANALY_DIPOLE
*/
void test_source_component_common_ConstantDecAnalyBeam(void) {
  test_source_component_common_ConstantDecChooseBeams(ANALY_DIPOLE, " ");
}

/*
This test checks source_component_common with beamtype=GAUSS_BEAM
*/
void test_source_component_common_ConstantDecGaussBeam(void) {
  test_source_component_common_ConstantDecChooseBeams(GAUSS_BEAM, " ");
}

/*
This test checks source_component_common with beamtype=NO_BEAM
*/
void test_source_component_common_ConstantDecNoBeam(void) {
  test_source_component_common_ConstantDecChooseBeams(NO_BEAM, " ");
}

/*
This test checks source_component_common with beamtype=FEE_BEAM_INTERP
*/
void test_source_component_common_ConstantDecFEEBeamInterp(void) {
  //Look for environment variable MWA_FEE_HDF5, which should point
  //towards mwa_full_embedded_element_pattern.h5. If we can't find it,
  //can't run the tests, so don't
  char* mwa_fee_hdf5 = getenv("MWA_FEE_HDF5_INTERP");

  if (mwa_fee_hdf5) {
    printf("MWA_FEE_HDF5_INTERP: %s\n", mwa_fee_hdf5 );
    test_source_component_common_ConstantDecChooseBeams(FEE_BEAM_INTERP, mwa_fee_hdf5);
  }
  else {
    printf("MWA_FEE_HDF5_INTERP not found - not running MWA_FEE beam test");
  }
}

/*
This test checks source_component_common with beamtype=MWA_ANALY
*/
void test_source_component_common_ConstantDecMWAAnaly(void) {
  test_source_component_common_ConstantDecChooseBeams(MWA_ANALY, " ");
}

//Run the test with unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_source_component_common_ConstantDecNoBeam);
    RUN_TEST(test_source_component_common_ConstantDecGaussBeam);
    RUN_TEST(test_source_component_common_ConstantDecAnalyBeam);
    RUN_TEST(test_source_component_common_ConstantDecFEEBeam);
    RUN_TEST(test_source_component_common_ConstantDecFEEBeamInterp);
    RUN_TEST(test_source_component_common_ConstantDecMWAAnaly);

    return UNITY_END();
}
