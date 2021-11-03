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
`RTS_PrecessXYZtoJ2000` uses woden_settings and takes the current X,Y,Z coords
in `array_layout` and precesses them back to J2000. Also updates the LST.

Test by reading giving a known set of X,Y,Z and a given example julian date, and
make sure outputs match expectations
*/
void test_RTS_PrecessXYZtoJ2000_GivesCorrectValues(void)
{

  woden_settings_t *woden_settings = malloc(sizeof(woden_settings_t));
  //Set up where woden_settings correctly
  woden_settings->array_layout_file_path="example_array_layout.txt";
  woden_settings->latitude = MWA_LAT_RAD;
  woden_settings->jd_date = 2457278.2010995;
  woden_settings->lst_base = LST_BEFORE;

  //Set up array_layout with some precalculated values
  array_layout_t * array_layout = malloc(sizeof(array_layout_t));

  array_layout->ant_X = expec_X_noprec;
  array_layout->ant_Y = expec_Y_noprec;
  array_layout->ant_Z = expec_Z_noprec;

  array_layout->num_baselines = 28;
  array_layout->num_tiles = 8;

  //Call the function we are testing
  RTS_PrecessXYZtoJ2000(array_layout, woden_settings);

  //Check rotations give correct answers

  for (int i = 0; i < 8; i++) {
    TEST_ASSERT_EQUAL_FLOAT((float)expec_X_prec[i],
                            (float)array_layout->ant_X[i]);
    TEST_ASSERT_EQUAL_FLOAT((float)expec_Y_prec[i],
                            (float)array_layout->ant_Y[i]);
    TEST_ASSERT_EQUAL_FLOAT((float)expec_Z_prec[i],
                            (float)array_layout->ant_Z[i]);
  }

  TEST_ASSERT_EQUAL_FLOAT(LST_AFTER, woden_settings->lst_base);

}

//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_RTS_PrecessXYZtoJ2000_GivesCorrectValues);

    return UNITY_END();
}
