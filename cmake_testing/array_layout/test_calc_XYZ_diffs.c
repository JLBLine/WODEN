#include <unity.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "array_layout.h"
#include "woden_struct_defs.h"
#include "test_RTS_XYZ_common.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

#define UNITY_INCLUDE_FLOAT

/*
`calc_XYZ_diffs` uses woden_settings to find an array text file, read in e,n,h
coords, transform to X,Y,Z, precess them back to J2000 if requested, and then
calculates the baseline length in X,Y,Z

Test by reading in a known set of e,n,h and a given example julian date, and
make sure outputs match expectations
*/
void test_calc_XYZ_diffs_GivesCorrectValues(int do_precession)
{

  woden_settings_t *woden_settings = malloc(sizeof(woden_settings_t));
  //Set up where woden_settings correctly
  woden_settings->array_layout_file_path="example_array_layout.txt";
  woden_settings->latitude = MWA_LAT_RAD;
  woden_settings->jd_date = JD_DATE;
  woden_settings->lst_base = LST_BEFORE;

  //Call the function we are testing
  array_layout_t * array_layout;
  array_layout = calc_XYZ_diffs(woden_settings, do_precession);

  //Test some basic quantities
  TEST_ASSERT_EQUAL_INT(28, array_layout->num_baselines);
  TEST_ASSERT_EQUAL_INT(8, array_layout->num_tiles);

  //These are the values in "example_array_layout.txt"
  float expec_east[] = {84, 12, 86, 780, 813, 899, 460, 810};
  float expec_north[] = {112, 202, 377, 227, 561, 600, 70, 524};
  float expec_height[] = {3, 3, 5, 5, 7, 1, 0, 8};

  //Check read them in correctly
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_east, array_layout->ant_east, 8);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_north, array_layout->ant_north, 8);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_height, array_layout->ant_height, 8);

  // for (int ant = 0; ant < array_layout->num_tiles; ant++) {
  //   printf("%.5f %.5f %.5f\n",array_layout->ant_X[ant],
  //                             array_layout->ant_Y[ant],
  //                             array_layout->ant_Z[ant] );
  // }

  // for (int baseline = 0; baseline < array_layout->num_baselines; baseline++) {
  //   printf("%.5f %.5f %.5f\n",array_layout->X_diff_metres[baseline],
  //                             array_layout->Y_diff_metres[baseline],
  //                             array_layout->Z_diff_metres[baseline] );
  // }

  //Check e,n,h are coverted to X,Y,Z correctly, and precessed back to J2000 if
  //requested
  if (do_precession == 0) {
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_X_noprec, array_layout->ant_X, 8);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_Y_noprec, array_layout->ant_Y, 8);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_Z_noprec, array_layout->ant_Z, 8);

    TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_X_diffs_noprec,
                                  array_layout->X_diff_metres, 28);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_Y_diffs_noprec,
                                  array_layout->Y_diff_metres, 28);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_Z_diffs_noprec,
                                  array_layout->Z_diff_metres, 28);
    // printf("Precess lst %.7f\n",woden_settings->lst_base );
    TEST_ASSERT_EQUAL_FLOAT(LST_BEFORE, woden_settings->lst_base);

  } else if (do_precession == 1) {
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_X_prec, array_layout->ant_X, 8);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_Y_prec, array_layout->ant_Y, 8);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_Z_prec, array_layout->ant_Z, 8);

    TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_X_diffs_prec,
                                  array_layout->X_diff_metres, 28);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_Y_diffs_prec,
                                  array_layout->Y_diff_metres, 28);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_Z_diffs_prec,
                                  array_layout->Z_diff_metres, 28);
    // printf("Precess lst %.7f\n",woden_settings->lst_base );
    TEST_ASSERT_EQUAL_FLOAT(LST_AFTER, woden_settings->lst_base);
  }
}

void test_calc_XYZ_diffs_GivesCorrectValuesNoPrecess(void)
{
  test_calc_XYZ_diffs_GivesCorrectValues(0);
}

void test_calc_XYZ_diffs_GivesCorrectValuesPrecess(void)
{
  test_calc_XYZ_diffs_GivesCorrectValues(1);
}
//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_calc_XYZ_diffs_GivesCorrectValuesNoPrecess);
    RUN_TEST(test_calc_XYZ_diffs_GivesCorrectValuesPrecess);

    return UNITY_END();
}