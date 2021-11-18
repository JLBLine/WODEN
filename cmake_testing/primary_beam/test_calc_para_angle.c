#include <unity.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "primary_beam.h"
#include "woden_struct_defs.h"
#include "expected_para_angles.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

#define UNITY_INCLUDE_FLOAT

/*
Check calc does what is expected give a predefined sky model `src`
*/
void test_calc_para_angle() {
  //Make sky model
  catsource_t *src = make_sky_model();

  //Function being tested
  int num_time_steps = 3;
  calc_para_angle(src, lsts, MWA_LAT_RAD, num_time_steps);

  #ifdef DOUBLE_PRECISION
    double TOL = 1e-15;
  #else
    double TOL = 1e-7;
  #endif

  for (int ang = 0; ang < 9; ang++) {
    TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_point_sin_para[ang],
                                   src->sin_point_para_angs[ang]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_point_cos_para[ang],
                                   src->cos_point_para_angs[ang]);

    TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_gauss_sin_para[ang],
                                   src->sin_gauss_para_angs[ang]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_gauss_cos_para[ang],
                                   src->cos_gauss_para_angs[ang]);

    TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_shape_sin_para[ang],
                                   src->sin_shape_para_angs[ang]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_shape_cos_para[ang],
                                   src->cos_shape_para_angs[ang]);
  }
}

//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_calc_para_angle);

    return UNITY_END();
}
