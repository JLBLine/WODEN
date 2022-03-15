#include "test_source_component_common.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

// void sincos(user_precision_t x, user_precision_t *sin, user_precision_t *cos);

//External CUDA code we're linking in
extern void test_source_component_common(int num_components,
           user_precision_complex_t *primay_beam_J00, user_precision_complex_t *primay_beam_J01,
           user_precision_complex_t *primay_beam_J10, user_precision_complex_t *primay_beam_J11,
           double *freqs, double *ls, double *ms, double *ns,
           double *ras, double *decs, user_precision_t *azs, user_precision_t *zas,
           user_precision_t *sin_para_angs, user_precision_t *cos_para_angs,
           double *beam_has, double *beam_decs,
           woden_settings_t *woden_settings,
           beam_settings_t *beam_settings);

extern void get_HDFBeam_normalisation(RTS_MWA_FEE_beam_t *FEE_beam_zenith,
                RTS_MWA_FEE_beam_t *FEE_beam);

extern void free_FEE_primary_beam_from_GPU(RTS_MWA_FEE_beam_t *primary_beam);

extern void copy_FEE_primary_beam_to_GPU(RTS_MWA_FEE_beam_t *FEE_beam);

//MWA FEE interp CUDA code
extern void multifreq_get_MWAFEE_normalisation(beam_settings_t *beam_settings);

extern void test_run_and_map_multifreq_calc_CUDA_FEE_beam(beam_settings_t *beam_settings,
    user_precision_t *azs, user_precision_t *zas, double latitude,
    user_precision_complex_t *primay_beam_J00,
    user_precision_complex_t *primay_beam_J01,
    user_precision_complex_t *primay_beam_J10,
    user_precision_complex_t *primay_beam_J11,
    int num_freqs, int NUM_COMPS, int num_times,
    int rotation, int scaling);

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

  double ra0 = 0.0*DD2R;
  double dec0 = 0.0*DD2R;

  double *decs = malloc(num_components*sizeof(double));
  double *zeroes = calloc(num_components, sizeof(double));

  //Keep RA between 0 and 2*pi here but enter RAs that should return
  //negative l values
  double ras[9] = {(3*M_PI)/2, (5*M_PI)/3, (7*M_PI)/4, (11*M_PI)/6,
                   0.0, M_PI/6, M_PI/4, M_PI/3, M_PI/2};

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

  beam_settings_t *beam_settings = malloc(sizeof(beam_settings_t));
  beam_settings->beamtype = beamtype;

  /*********************************************************************
  Code used to generate the az / za and parallactic angles is below
  Stick into permanent arrays so test doesn't rely on erfa.
  **********************************************************************

  user_precision_t *zas = malloc(num_components*num_times*sizeof(user_precision_t));
  user_precision_t *azs = malloc(num_components*num_times*sizeof(user_precision_t));
  user_precision_t *sin_para_angs = malloc(num_components*num_times*sizeof(user_precision_t));
  user_precision_t *cos_para_angs = malloc(num_components*num_times*sizeof(user_precision_t));

  //Let's say the sources are stationary on the sky to keep the maths in the
  //test to a minimum
  double erfa_az, el, ha, para_angle;
  double lst = 0.0;
  double latitude = 0.0;

  for (int component = 0; component < num_components; component++) {
    for (int time_ind = 0; time_ind < num_times; time_ind++) {

      double ha = lst - ras[component];

      eraHd2ae(ha, (double)decs[component], latitude, &erfa_az, &el );

      azs[component*num_times + time_ind] = (user_precision_t)erfa_az;
      zas[component*num_times + time_ind] = M_PI / 2. - (user_precision_t)el;

      para_angle = eraHd2pa(ha, (double)decs[component], latitude);

      sin_para_angs[component*num_times + time_ind] = (user_precision_t)sin(para_angle);
      cos_para_angs[component*num_times + time_ind] = (user_precision_t)cos(para_angle);

      printf("%d %.8f %.8f %.8f %.8f\n",component*num_times + time_ind,
                                 (user_precision_t)erfa_az, M_PI / 2. - (user_precision_t)el,
                                 (user_precision_t)sin(para_angle), (user_precision_t)cos(para_angle));
    }
  }
  */

  double freqs[] = {100e+6, 200e+6};

  double *beam_has = malloc(num_components*num_times*sizeof(double));
  double *beam_decs = malloc(num_components*num_times*sizeof(double));

  //Make ha/dec coords if using the Gaussian beam
  if (beamtype == GAUSS_BEAM) {
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

  beam_settings->FEE_beam = malloc(sizeof(RTS_MWA_FEE_beam_t));
  //If FEE_BEAM, call the C code to interrogate the hdf5 file and set beam
  //things up
  if (beamtype == FEE_BEAM) {

    //Get a zenith pointing beam for normalisation purposes
    RTS_MWA_FEE_beam_t *FEE_beam_zenith = malloc(sizeof(RTS_MWA_FEE_beam_t));

    double base_middle_freq = 150e+6;
//
    user_precision_t user_precision_t_zenith_delays[16] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    //
    printf("\n\tSetting up the zenith FEE beam...");
    RTS_MWAFEEInit(mwa_fee_hdf5, base_middle_freq, FEE_beam_zenith, user_precision_t_zenith_delays);
    printf(" done.\n");

    printf("\tSetting up the FEE beam...");
    RTS_MWAFEEInit(mwa_fee_hdf5, base_middle_freq, beam_settings->FEE_beam, user_precision_t_zenith_delays);
    printf(" done.\n");

    printf("\tGetting FEE beam normalisation...");
    get_HDFBeam_normalisation(FEE_beam_zenith, beam_settings->FEE_beam);
    //Free the zenith pointing as done with it now
    free_FEE_primary_beam_from_GPU(FEE_beam_zenith);
    printf(" done.\n");

    printf("\tCopying the FEE beam across to the GPU...");
    copy_FEE_primary_beam_to_GPU(beam_settings->FEE_beam);
    printf(" done.\n");

  }

  else if (beamtype == FEE_BEAM_INTERP) {

    woden_settings->base_low_freq = 100e+6;
    woden_settings->num_freqs = 2;
    woden_settings->hdf5_beam_path = mwa_fee_hdf5;



    //Setup up the MWA FEE beams on the CPU
    multifreq_RTS_MWAFEEInit(beam_settings,  woden_settings, freqs);

    //Send them to GPU and calculate normalisations
    multifreq_get_MWAFEE_normalisation(beam_settings);

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

  //Run the CUDA code
  test_source_component_common(num_components,
             primay_beam_J00, primay_beam_J01,
             primay_beam_J10, primay_beam_J11,
             freqs, ls, ms, ns,
             ras, decs, azs, zas,
             sin_para_angs, cos_para_angs,
             beam_has, beam_decs,
             woden_settings,
             beam_settings);

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
      TOL = 3e-2;
    #endif

    for (int output = 0; output < num_beam_values; output++) {
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_J00_re[output], creal(primay_beam_J00[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_J00_im[output], cimag(primay_beam_J00[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_J01_re[output], creal(primay_beam_J01[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_J01_im[output], cimag(primay_beam_J01[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_J10_re[output], creal(primay_beam_J10[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_J10_im[output], cimag(primay_beam_J10[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_J11_re[output], creal(primay_beam_J11[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_J11_im[output], cimag(primay_beam_J11[output]));
    }

  }
  else if (beamtype == FEE_BEAM_INTERP) {

    #ifdef DOUBLE_PRECISION
      TOL = 1e-7;
    #else
      TOL = 3e-2;
    #endif

    for (int output = 0; output < num_beam_values; output++) {
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_interp_J00_re[output], creal(primay_beam_J00[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_interp_J00_im[output], cimag(primay_beam_J00[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_interp_J01_re[output], creal(primay_beam_J01[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_interp_J01_im[output], cimag(primay_beam_J01[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_interp_J10_re[output], creal(primay_beam_J10[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_interp_J10_im[output], cimag(primay_beam_J10[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_interp_J11_re[output], creal(primay_beam_J11[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_interp_J11_im[output], cimag(primay_beam_J11[output]));
    }

    for (int freq_ind = 0; freq_ind < beam_settings->num_MWAFEE; freq_ind++) {

      RTS_MWA_FEE_beam_t *FEE_beam = &beam_settings->FEE_beams[freq_ind];
      RTS_MWA_FEE_beam_t *FEE_beam_zenith = &beam_settings->FEE_beam_zeniths[freq_ind];

      RTS_freeHDFBeam(FEE_beam);
      RTS_freeHDFBeam(FEE_beam_zenith);

    }
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
