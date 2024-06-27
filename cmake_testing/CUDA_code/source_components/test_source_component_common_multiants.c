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
  //WARNING azs and zas are stored in a header so we don't have to link/call
  //erfa or pal. Which means you CANNOT CHANGE NUMBER OF TIMES HERE
  int num_times = 3;
  int num_freqs = 2;
  int num_beams = 3;

  int num_components = num_powers + num_curves + num_lists;

  int num_beam_values = num_times*num_freqs*num_components*num_beams;

  double ra0 = 0.0*DD2R;
  double dec0 = 0.0*DD2R;

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
  woden_settings->use_dipamps = 1;
  woden_settings->num_ants = num_beams;
  woden_settings->num_autos = 0;
  woden_settings->do_autos = 0;

  user_precision_t *zeroes = calloc(num_components, sizeof(user_precision_t));

  //Keep RA between 0 and 2*pi here but enter RAs that should return
  //negative l values
  double ras[9] = {(3*M_PI)/2, (5*M_PI)/3, (7*M_PI)/4, (11*M_PI)/6,
                   0.0, M_PI/6, M_PI/4, M_PI/3, M_PI/2};

  // double ras[9] = {(3*M_PI)/2, (3*M_PI)/2, (3*M_PI)/2, (3*M_PI)/2,
  //                  (3*M_PI)/2, (3*M_PI)/2, (3*M_PI)/2, (3*M_PI)/2, (3*M_PI)/2};

  

  double *decs = malloc(num_components*sizeof(double));
  for (int i = 0; i < num_components; i++) {
    decs[i] = dec0;
  }

  

  beam_settings_t *beam_settings = malloc(sizeof(beam_settings_t));
  beam_settings->beamtype = beamtype;

  double freqs[] = {100e+6, 200e+6};

  double *beam_has = malloc(num_components*num_times*sizeof(double));
  double *beam_decs = malloc(num_components*num_times*sizeof(double));


  //If FEE_BEAM, call the C code to interrogate the hdf5 file and set beam
  //things up
  if (beamtype == FEE_BEAM || beamtype == FEE_BEAM_INTERP) {

    int32_t status =  new_fee_beam(mwa_fee_hdf5, &beam_settings->fee_beam);

    TEST_ASSERT_EQUAL(0, status);

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

    int num_delays_per_tile = 16;

    beam_settings->hyper_delays = (uint32_t*)malloc(num_beams*num_delays_per_tile*sizeof(uint32_t));

    for (int delay = 0; delay < num_delays_per_tile; delay++) {
      for (int tile = 0; tile < num_beams; tile++) {
        beam_settings->hyper_delays[tile*num_delays_per_tile + delay] = 0;
      }
    }

    //These are the amplitudes for the dipoles, as read in from metafits
    // I believe that they have X - east-west, Y - north-south
    double amps[96] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                      0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                      0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4,
                      0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6,
                      0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8,
                      1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

    // double amps[96] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    //                   1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    //                   1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    //                   1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    //                   1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    //                   1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

    uint32_t num_amps = 32;
    uint8_t norm_to_zenith = 1;

    status = new_gpu_fee_beam(beam_settings->fee_beam,
                               freqs_hz,
                               beam_settings->hyper_delays,
                               amps,
                               num_freqs_hyper,
                               num_beams,
                               num_amps,
                               norm_to_zenith,
                               &beam_settings->cuda_fee_beam);

    // printf("STATUS %d\n", status);

    if (status != 0) {
      handle_hyperbeam_error(__FILE__, __LINE__, "new_gpu_fee_beam");
    }

    TEST_ASSERT_EQUAL(0, status);

    

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


  double antx_mult[3] = {0.2, 0.6, 1.0};
  double anty_mult[3] = {0.0, 0.4, 0.8};

  //Depending on beamtype, check results match expectations

  //For the Gaussian beam, the way the l,m coords are set up measns we can
  //analytically predict the values, so go through results and calcluate
  //expected values, and compare to what we got
  int beam_ind = 0;
  int expec_ind = 0;
  
  if (beamtype == FEE_BEAM) {

    #ifdef DOUBLE_PRECISION
      TOL = 1e-7;
    #else
      TOL = 1e-7;
    #endif

    for (int ant = 0; ant < num_beams; ant ++) {
      for (int time = 0; time < num_times; time ++) {
        for (int freq = 0; freq < num_freqs; freq ++) {
          for (int comp = 0; comp < num_components; comp ++) {
  
            beam_ind = ant*num_freqs*num_times*num_components + num_freqs*time*num_components + num_components*freq + comp;
            expec_ind =  num_freqs*time*num_components + num_components*freq + comp;

            TEST_ASSERT_DOUBLE_WITHIN(TOL, antx_mult[ant]*fee_expec_J00_re[expec_ind],
                                      creal(primay_beam_J00[beam_ind]) );
            TEST_ASSERT_DOUBLE_WITHIN(TOL, antx_mult[ant]*fee_expec_J00_im[expec_ind],
                                      cimag(primay_beam_J00[beam_ind]) );
            TEST_ASSERT_DOUBLE_WITHIN(TOL, antx_mult[ant]*fee_expec_J01_re[expec_ind],
                                      creal(primay_beam_J01[beam_ind]) );
            TEST_ASSERT_DOUBLE_WITHIN(TOL, antx_mult[ant]*fee_expec_J01_im[expec_ind],
                                      cimag(primay_beam_J01[beam_ind]) );

            TEST_ASSERT_DOUBLE_WITHIN(TOL, anty_mult[ant]*fee_expec_J10_re[expec_ind],
                                      creal(primay_beam_J10[beam_ind]) );
            TEST_ASSERT_DOUBLE_WITHIN(TOL, anty_mult[ant]*fee_expec_J10_im[expec_ind],
                                      cimag(primay_beam_J10[beam_ind]) );
            TEST_ASSERT_DOUBLE_WITHIN(TOL, anty_mult[ant]*fee_expec_J11_re[expec_ind],
                                      creal(primay_beam_J11[beam_ind]) );
            TEST_ASSERT_DOUBLE_WITHIN(TOL, anty_mult[ant]*fee_expec_J11_im[expec_ind],
                                      cimag(primay_beam_J11[beam_ind]) );
  
          }
        }
      }
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

    for (int ant = 0; ant < num_beams; ant ++) {
      for (int time = 0; time < num_times; time ++) {
        for (int freq = 0; freq < num_freqs; freq ++) {
          for (int comp = 0; comp < num_components; comp ++) {
  
            int beam_ind = ant*num_freqs*num_times*num_components + num_freqs*time*num_components + num_components*freq + comp;

            int expec_ind =  num_freqs*time*num_components + num_components*freq + comp;

            TEST_ASSERT_DOUBLE_WITHIN(TOL, antx_mult[ant]*fee_expec_interp_J00_re[expec_ind],
                                      creal(primay_beam_J00[beam_ind]) );
            TEST_ASSERT_DOUBLE_WITHIN(TOL, antx_mult[ant]*fee_expec_interp_J00_im[expec_ind],
                                      cimag(primay_beam_J00[beam_ind]) );
            TEST_ASSERT_DOUBLE_WITHIN(TOL, antx_mult[ant]*fee_expec_interp_J01_re[expec_ind],
                                      creal(primay_beam_J01[beam_ind]) );
            TEST_ASSERT_DOUBLE_WITHIN(TOL, antx_mult[ant]*fee_expec_interp_J01_im[expec_ind],
                                      cimag(primay_beam_J01[beam_ind]) );

            TEST_ASSERT_DOUBLE_WITHIN(TOL, anty_mult[ant]*fee_expec_interp_J10_re[expec_ind],
                                      creal(primay_beam_J10[beam_ind]) );
            TEST_ASSERT_DOUBLE_WITHIN(TOL, anty_mult[ant]*fee_expec_interp_J10_im[expec_ind],
                                      cimag(primay_beam_J10[beam_ind]) );
            TEST_ASSERT_DOUBLE_WITHIN(TOL, anty_mult[ant]*fee_expec_interp_J11_re[expec_ind],
                                      creal(primay_beam_J11[beam_ind]) );
            TEST_ASSERT_DOUBLE_WITHIN(TOL, anty_mult[ant]*fee_expec_interp_J11_im[expec_ind],
                                      cimag(primay_beam_J11[beam_ind]) );
  
          }
        }
      }
    }

    free_fee_beam(beam_settings->fee_beam);
    free_gpu_fee_beam(beam_settings->cuda_fee_beam);
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


//Run the test with unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_source_component_common_ConstantDecFEEBeamPoint);
    RUN_TEST(test_source_component_common_ConstantDecFEEBeamGauss);
    RUN_TEST(test_source_component_common_ConstantDecFEEBeamShapelet);

    RUN_TEST(test_source_component_common_ConstantDecFEEBeamInterpPoint);
    RUN_TEST(test_source_component_common_ConstantDecFEEBeamInterpGaussian);
    RUN_TEST(test_source_component_common_ConstantDecFEEBeamInterpShapelet);

    return UNITY_END();
}
