#include <unity.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "read_and_write.h"
#include "create_sky_model.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

#define UNITY_INCLUDE_FLOAT

/*
Check that the values in the json that aren't beam specific are read in
correctly
*/
void check_observation_params(woden_settings_t *woden_settings) {

  TEST_ASSERT_EQUAL_FLOAT(-26.70331944*DD2R,woden_settings->latitude);
  TEST_ASSERT_EQUAL_FLOAT(0.44312771*DD2R, woden_settings->lst_base);
  TEST_ASSERT_EQUAL_FLOAT(0.0, woden_settings->ra0);
  TEST_ASSERT_EQUAL_FLOAT(-27*DD2R, woden_settings->dec0);
  TEST_ASSERT_EQUAL_INT(16, woden_settings->num_freqs);
  TEST_ASSERT_EQUAL_FLOAT(40000.0, woden_settings->frequency_resolution);
  TEST_ASSERT_EQUAL_FLOAT(1.67035e+08, woden_settings->base_low_freq);
  TEST_ASSERT_EQUAL_FLOAT(1.28e+06, woden_settings->coarse_band_width);
  TEST_ASSERT_EQUAL_INT(4, woden_settings->num_time_steps);
  TEST_ASSERT_EQUAL_FLOAT(2.0, woden_settings->time_res);
  TEST_ASSERT_EQUAL_STRING("srclist_singlepoint.txt", woden_settings->cat_filename);
  TEST_ASSERT_EQUAL_FLOAT(2457278.2010995, woden_settings->jd_date);
  TEST_ASSERT_EQUAL_INT(5000 ,woden_settings->chunking_size);

  TEST_ASSERT_EQUAL_INT(1, woden_settings->array_layout_file);
  TEST_ASSERT_EQUAL_STRING("example_array_layout.txt", woden_settings->array_layout_file_path);
  TEST_ASSERT_EQUAL_INT(3, woden_settings->num_bands);

  int expect_band_nums[] = {1,4,9};
  TEST_ASSERT_EQUAL_INT_ARRAY(expect_band_nums, woden_settings->band_nums, 3);

  TEST_ASSERT_EQUAL_INT(CROP_COMPONENTS, woden_settings->sky_crop_type);

}

/*
Read in a json settings file and check values are good
Check that MWA FEE beam specific parameters are read in correctly
*/
void test_read_json_settings_MWAFEE(void) {

  int status=0;
  //Read in the settings from the controlling json file
  woden_settings_t *woden_settings = malloc( sizeof(woden_settings_t) );
  status = read_json_settings("run_woden_MWAFEE.json", woden_settings);
  //
  TEST_ASSERT_EQUAL_INT(0, status);

  //Check generic observation params have been read in correctly
  check_observation_params(woden_settings);

  //Check MWA FEE beam specific values
  float expect_delays[] = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(expect_delays, woden_settings->FEE_ideal_delays, 16);
  TEST_ASSERT_EQUAL_INT(FEE_BEAM, woden_settings->beamtype);
  TEST_ASSERT_EQUAL_STRING("/home/jline/software/useful/mwa_full_embedded_element_pattern.h5",
                            woden_settings->hdf5_beam_path);

}

/*
Read in a json settings file and check values are good
Check that Gaussian beam parameters are read in correctly
*/
void test_read_json_settings_GaussBeam(char *json_path, int default_gauss) {

  int status=0;
  //Read in the settings from the controlling json file
  woden_settings_t *woden_settings = malloc( sizeof(woden_settings_t) );
  status = read_json_settings(json_path, woden_settings);

  TEST_ASSERT_EQUAL_INT(0, status);

  //Check generic observation params have been read in correctly
  check_observation_params(woden_settings);

  TEST_ASSERT_EQUAL_FLOAT(56.13129905*DD2R, woden_settings->gauss_ra_point);
  TEST_ASSERT_EQUAL_FLOAT(-39.47577940*DD2R, woden_settings->gauss_dec_point);
  TEST_ASSERT_EQUAL_INT(GAUSS_BEAM, woden_settings->beamtype);

  //If no arguments provided for FWHM and ref freq should have default values
  if (default_gauss) {

    TEST_ASSERT_EQUAL_FLOAT(20.0, woden_settings->gauss_beam_FWHM);
    TEST_ASSERT_EQUAL_FLOAT(150e+6, woden_settings->gauss_beam_ref_freq);
  }
  //If not, should have values specified
  else {
    TEST_ASSERT_EQUAL_FLOAT(60.0, woden_settings->gauss_beam_FWHM);
    TEST_ASSERT_EQUAL_FLOAT(75e+6, woden_settings->gauss_beam_ref_freq);

  }
}


/*
Runs without the Gaussian FWHM and ref frequency specified, so should have
default values
*/
void test_read_json_settings_GaussBeamDefault(void) {
  test_read_json_settings_GaussBeam("run_woden_gaussian_default.json", 1);
}

/*
Runs with the Gaussian FWHM and ref frequency specified, so should have
bespoke values
*/
void test_read_json_settings_GaussBeamBespoke(void) {
  test_read_json_settings_GaussBeam("run_woden_gaussian_bespoke.json", 0);
}

/*
Read in a json settings file and check values are good
Check that EDA2 beam specific parameters are read in correctly
*/
void test_read_json_settings_EDA2(void) {

  int status=0;
  //Read in the settings from the controlling json file
  woden_settings_t *woden_settings = malloc( sizeof(woden_settings_t) );
  status = read_json_settings("run_woden_EDA2.json", woden_settings);

  TEST_ASSERT_EQUAL_INT(0, status);

  //Check generic observation params have been read in correctly
  check_observation_params(woden_settings);

  //Check EDA2 beam specific value
  TEST_ASSERT_EQUAL_INT(ANALY_DIPOLE, woden_settings->beamtype);

}

/*
Read in a json settings file and check values are good
Check that using no beam is selected
*/
void test_read_json_settings_NoBeam(void) {

  int status=0;
  //Read in the settings from the controlling json file
  woden_settings_t *woden_settings = malloc( sizeof(woden_settings_t) );
  status = read_json_settings("run_woden_nobeam.json", woden_settings);

  TEST_ASSERT_EQUAL_INT(0, status);

  //Check generic observation params have been read in correctly
  check_observation_params(woden_settings);

  //Check EDA2 beam specific value
  TEST_ASSERT_EQUAL_INT(NO_BEAM, woden_settings->beamtype);

}

/*
Read in a json settings file that had mulitple primary beam types selected
This should error
*/
void test_read_json_settings_MultiBeam(void) {

  int status=0;
  //Read in the settings from the controlling json file
  woden_settings_t *woden_settings = malloc( sizeof(woden_settings_t) );
  status = read_json_settings("run_woden_multiple_beams.json", woden_settings);

  TEST_ASSERT_EQUAL_INT(1, status);

}

/*
Read in a json settings file that only has 15 delay settings for MWA FEE beam
This should error
*/
void test_read_json_settings_MWAFEEBadDelay(void) {

  int status=0;
  //Read in the settings from the controlling json file
  woden_settings_t *woden_settings = malloc( sizeof(woden_settings_t) );
  status = read_json_settings("run_woden_MWAFEE_baddelay.json", woden_settings);

  TEST_ASSERT_EQUAL_INT(1, status);

}

/*
Read in a json settings file that has no hdf5_beam_path set for an MWA FEE sim
This should error
*/
void test_read_json_settings_MWAFEENoPath(void) {

  int status=0;
  //Read in the settings from the controlling json file
  woden_settings_t *woden_settings = malloc( sizeof(woden_settings_t) );
  status = read_json_settings("run_woden_MWAFEE_nopath.json", woden_settings);

  TEST_ASSERT_EQUAL_INT(1, status);

}



//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_read_json_settings_MWAFEE);
    RUN_TEST(test_read_json_settings_GaussBeamDefault);
    RUN_TEST(test_read_json_settings_GaussBeamBespoke);
    RUN_TEST(test_read_json_settings_EDA2);
    RUN_TEST(test_read_json_settings_NoBeam);
    RUN_TEST(test_read_json_settings_MultiBeam);
    RUN_TEST(test_read_json_settings_MWAFEEBadDelay);
    RUN_TEST(test_read_json_settings_MWAFEENoPath);

    return UNITY_END();
}
