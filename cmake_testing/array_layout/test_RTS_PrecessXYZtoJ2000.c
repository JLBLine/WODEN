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
  woden_settings->latitude = MWA_LAT_RAD;
  woden_settings->jd_date = 2457278.2010995;
  woden_settings->lst_obs_epoch_base = LST_BEFORE;
  woden_settings->lst_base = LST_AFTER;
  woden_settings->do_precession = 1;

  woden_settings->num_time_steps = 2;
  woden_settings->lsts = lsts;
  woden_settings->time_res = 8;
  woden_settings->mjds = mjds;

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

  for (int i = 0; i < array_layout->num_tiles*woden_settings->num_time_steps; i++) {
    TEST_ASSERT_EQUAL_DOUBLE(expec_X_prec[i],
                             array_layout->ant_X[i]);
    TEST_ASSERT_EQUAL_DOUBLE(expec_Y_prec[i],
                             array_layout->ant_Y[i]);
    TEST_ASSERT_EQUAL_DOUBLE(expec_Z_prec[i],
                             array_layout->ant_Z[i]);

  }


  // for (int i = 0; i < array_layout->num_tiles*woden_settings->num_time_steps; i++) {
  //   printf("%.12f\n", array_layout->ant_X[i]);
  // }
  // printf("----\n");
  // for (int i = 0; i < array_layout->num_tiles*woden_settings->num_time_steps; i++) {
  //   printf("%.12f\n", array_layout->ant_Y[i]);
  // }
  // printf("----\n");
  // for (int i = 0; i < array_layout->num_tiles*woden_settings->num_time_steps; i++) {
  //   printf("%.12f\n", array_layout->ant_Z[i]);
  // }

  // TEST_ASSERT_EQUAL_DOUBLE(LST_AFTER, woden_settings->lst_base);

}

//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_RTS_PrecessXYZtoJ2000_GivesCorrectValues);

    return UNITY_END();
}
