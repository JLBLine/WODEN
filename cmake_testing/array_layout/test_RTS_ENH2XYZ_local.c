#include <unity.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "array_layout.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

#define UNITY_INCLUDE_FLOAT

/*
RTS_ENH2XYZ_local transforms a local east, north, height coord into X,Y,Z
coords, which can be used to calculate u,v,w. Test for multiple coords and two
example latitudes
*/

void test_RTS_ENH2XYZ_local_GivesCorrectValues(void)
{
  int num_coords = 10;

  double X,Y,Z;
  double east, north, height;

  double TOL = 1e-13;

  //Loop over 10 example coords
  for (int coord = 0; coord < num_coords; coord++) {
    east = (coord + 1.)*10.;
    north = (coord + 2.)*10.;
    height = (coord + 3.)*10.;

    //test a latitude of 0.0
    RTS_ENH2XYZ_local(east, north, height, 0.0,
                      &X, &Y, &Z);

    TEST_ASSERT_DOUBLE_WITHIN(TOL, height, X);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, east, Y);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, north, Z);

    //test a latitude of -30.0
    RTS_ENH2XYZ_local(east, north, height, -30.0*DD2R,
                      &X, &Y, &Z);

    TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.5*north + (sqrt(3)/2)*height, X);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, east, Y);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, (sqrt(3)/2)*north + -0.5*height, Z);

  }
}

//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_RTS_ENH2XYZ_local_GivesCorrectValues);

    return UNITY_END();
}
