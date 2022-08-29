#include <unity.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "array_layout.h"
#include "woden_struct_defs.h"
#include "test_RTS_XYZ_common.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

/*
`calc_XYZ_diffs` uses woden_settings to find an array text file, read in e,n,h
coords, transform to X,Y,Z, precess them back to J2000 if requested, and then
calculates the baseline length in X,Y,Z

Test by reading in a known set of e,n,h and a given example julian date, and
make sure outputs match expectations
*/
void test_calc_XYZ_diffs_GivesCorrectValues(int do_precession, int do_autos)
{

  woden_settings_t *woden_settings = malloc(sizeof(woden_settings_t));
  //Set up where woden_settings correctly
  woden_settings->array_layout_file_path="example_array_layout.txt";
  woden_settings->latitude = MWA_LAT_RAD;
  woden_settings->latitude_obs_epoch_base = MWA_LAT_RAD;
  woden_settings->jd_date = 2457278.2010995;
  woden_settings->lst_obs_epoch_base = LST_BEFORE;

  if (do_precession) {
    woden_settings->lst_base = LST_AFTER;
    woden_settings->latitude = MWA_LAT_RAD;
  }
  else{
    woden_settings->lst_base = LST_BEFORE;
  }

  woden_settings->do_precession = do_precession;

  woden_settings->num_time_steps = 2;
  woden_settings->lsts = lsts;
  woden_settings->time_res = 8;
  woden_settings->mjds = mjds;
  woden_settings->num_freqs = 3;

  woden_settings->do_autos = do_autos;

  //Call the function we are testing
  array_layout_t * array_layout;
  array_layout = calc_XYZ_diffs(woden_settings, do_precession);

  //Test some basic quantities
  TEST_ASSERT_EQUAL_INT(28, array_layout->num_baselines);
  TEST_ASSERT_EQUAL_INT(8, array_layout->num_tiles);

  int expec_num_visi, expec_num_cross, expec_num_auto;

  expec_num_cross = woden_settings->num_time_steps*array_layout->num_baselines*woden_settings->num_freqs;

  if (do_autos) {
    expec_num_auto = woden_settings->num_time_steps*array_layout->num_tiles*woden_settings->num_freqs;
  } else {
    expec_num_auto = 0;
  }

  expec_num_visi = expec_num_auto + expec_num_cross;

  TEST_ASSERT_EQUAL_INT(expec_num_visi, woden_settings->num_visis);
  TEST_ASSERT_EQUAL_INT(expec_num_cross, woden_settings->num_cross);
  TEST_ASSERT_EQUAL_INT(expec_num_auto, woden_settings->num_autos);

  //These are the values in "example_array_layout.txt"
  double expec_east[] = {84., 12., 86., 780., 813., 899., 460., 810.};
  double expec_north[] = {112., 202., 377., 227., 561., 600., 70., 524.};
  double expec_height[] = {3., 3., 5., 5., 7., 1., 0., 8.};

  //Check read them in correctly
  TEST_ASSERT_EQUAL_DOUBLE_ARRAY(expec_east, array_layout->ant_east, 8);
  TEST_ASSERT_EQUAL_DOUBLE_ARRAY(expec_north, array_layout->ant_north, 8);
  TEST_ASSERT_EQUAL_DOUBLE_ARRAY(expec_height, array_layout->ant_height, 8);

  // for (int ant = 0; ant < array_layout->num_tiles; ant++) {
  //   printf("%.12f %.12f %.12f\n",array_layout->ant_X[ant],
  //                             array_layout->ant_Y[ant],
  //                             array_layout->ant_Z[ant] );
  // }
  //
  // for (int baseline = 0; baseline < array_layout->num_baselines; baseline++) {
  //   printf("%.12f %.12f %.12f\n",array_layout->X_diff_metres[baseline],
  //                             array_layout->Y_diff_metres[baseline],
  //                             array_layout->Z_diff_metres[baseline] );
  // }

  //Check e,n,h are coverted to X,Y,Z correctly, and precessed back to J2000 if
  //requested
  //Explicitly make things float for checking so this works for both float
  //and double precision compiled versions of the code
  if (do_precession == 0) {

    for (int i = 0; i < 16; i++) {
      TEST_ASSERT_EQUAL_DOUBLE(expec_X_noprec[i],
                               array_layout->ant_X[i]);
      TEST_ASSERT_EQUAL_DOUBLE(expec_Y_noprec[i],
                               array_layout->ant_Y[i]);
      TEST_ASSERT_EQUAL_DOUBLE(expec_Z_noprec[i],
                               array_layout->ant_Z[i]);
    }

    for (int i = 0; i < 56; i++) {
      TEST_ASSERT_EQUAL_DOUBLE(expec_X_diffs_noprec[i],
                               array_layout->X_diff_metres[i]);
      TEST_ASSERT_EQUAL_DOUBLE(expec_Y_diffs_noprec[i],
                               array_layout->Y_diff_metres[i]);
      TEST_ASSERT_EQUAL_DOUBLE(expec_Z_diffs_noprec[i],
                               array_layout->Z_diff_metres[i]);
    }

    // // printf("Precess lst %.7f\n",woden_settings->lst_base );
    // TEST_ASSERT_EQUAL_DOUBLE(LST_BEFORE, woden_settings->lst_base);

  } else if (do_precession == 1) {
    for (int i = 0; i < 16; i++) {
      TEST_ASSERT_EQUAL_DOUBLE(expec_X_prec[i],
                              array_layout->ant_X[i]);
      TEST_ASSERT_EQUAL_DOUBLE(expec_Y_prec[i],
                              array_layout->ant_Y[i]);
      TEST_ASSERT_EQUAL_DOUBLE(expec_Z_prec[i],
                              array_layout->ant_Z[i]);
    }

    for (int i = 0; i < 56; i++) {
      TEST_ASSERT_EQUAL_DOUBLE(expec_X_diffs_prec[i],
                              array_layout->X_diff_metres[i]);
      TEST_ASSERT_EQUAL_DOUBLE(expec_Y_diffs_prec[i],
                              array_layout->Y_diff_metres[i]);
      TEST_ASSERT_EQUAL_DOUBLE(expec_Z_diffs_prec[i],
                              array_layout->Z_diff_metres[i]);
    }

  }
}

void test_calc_XYZ_diffs_GivesCorrectValuesNoPrecessNoAutos(void)
{
  int do_precession = 0;
  int do_autos = 0;
  test_calc_XYZ_diffs_GivesCorrectValues(do_precession, do_autos);
}

void test_calc_XYZ_diffs_GivesCorrectValuesPrecessNoAutos(void)
{
  int do_precession = 1;
  int do_autos = 0;
  test_calc_XYZ_diffs_GivesCorrectValues(do_precession, do_autos);
}
//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_calc_XYZ_diffs_GivesCorrectValuesNoPrecessNoAutos);
    RUN_TEST(test_calc_XYZ_diffs_GivesCorrectValuesPrecessNoAutos);

    return UNITY_END();
}
