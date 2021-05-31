#include <unity.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "read_and_write.h"
#include "woden_struct_defs.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

#define UNITY_INCLUDE_FLOAT

/*
RTS_ENH2XYZ_local transforms a local east, north, height coord into X,Y,Z
coords, which can be used to calculate u,v,w. Test for multiple coords and two
example latitudes
*/

void test_calc_XYZ_diffs_GivesCorrectValues(int do_precession)
{

  woden_settings_t *woden_settings = malloc(sizeof(woden_settings_t));
  //Set up where woden_settings correctly
  woden_settings->array_layout_file_path="example_array_layout.txt";
  woden_settings->latitude = MWA_LAT_RAD;
  woden_settings->jd_date = 2457278.2010995;
  woden_settings->lst_base = 0.44312771*DD2R;

  //Call the function we are testing
  array_layout_t * array_layout;
  array_layout = calc_XYZ_diffs(woden_settings, do_precession);

  printf("Number of baselines: %d\n",array_layout->num_baselines);

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

  float expec_X_noprec[] = {53.00956, 93.45293 , 173.87949, 106.47388,
                      258.35040, 270.51578, 31.45595, 242.61703};

  float expec_Y_noprec[] = {84.00000, 12.00000, 86.00000, 780.00000,
                      813.00000, 899.00000, 460.00000, 810.00000};

  float expec_Z_noprec[] = {98.70657, 179.10765, 334.54434, 200.54254,
                      498.02115, 535.55786, 62.53418, 464.51801};


  float expec_X_prec[] = {53.16233, 93.72652, 174.39215, 106.80007,
                    259.13132, 271.35620, 31.56328, 243.34676};

  float expec_Y_prec[] = {83.99356, 11.98839, 85.97833, 779.98694,
                    812.96777, 898.96552, 459.99597, 809.96991};

  float expec_Z_prec[] = {98.62984, 178.96541, 334.28296, 200.41980,
                    497.66794, 535.19043, 62.50969, 464.18869};


  for (int ant = 0; ant < array_layout->num_tiles; ant++) {
    printf("%.5f %.5f %.5f\n",array_layout->ant_X[ant],
                              array_layout->ant_Y[ant],
                              array_layout->ant_Z[ant] );

  if (do_precession == 0) {
    // Do some stuff
  }


  }

  // int num_coords = 10;
  //
  // float X,Y,Z;
  // float east, north, height;
  //
  // //Loop over 10 example coords
  // for (int coord = 0; coord < num_coords; coord++) {
  //   east = (coord + 1)*10;
  //   north = (coord + 2)*10;
  //   height = (coord + 3)*10;
  //
  //   //test a latitude of 0.0
  //   RTS_ENH2XYZ_local(east, north, height, 0.0,
  //                     &X, &Y, &Z);
  //   // printf("%.5f %.5f %.5f\n",X, Y, Z );
  //
  //   TEST_ASSERT_EQUAL_FLOAT(height, X);
  //   TEST_ASSERT_EQUAL_FLOAT(east, Y);
  //   TEST_ASSERT_EQUAL_FLOAT(north, Z);
  //
  //   //test a latitude of -30.0
  //   RTS_ENH2XYZ_local(east, north, height, -30.0*DD2R,
  //                     &X, &Y, &Z);
  //
  //   TEST_ASSERT_EQUAL_FLOAT(0.5*north + (sqrt(3)/2)*height, X);
  //   TEST_ASSERT_EQUAL_FLOAT(east, Y);
  //   TEST_ASSERT_EQUAL_FLOAT((sqrt(3)/2)*north + -0.5*height, Z);
  //
  // }
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
