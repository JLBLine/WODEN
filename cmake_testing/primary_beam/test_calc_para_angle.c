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
  calc_para_angle(src, lsts, src->point_decs[0], num_time_steps);

  TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_point_sin_para,
                                src->sin_point_para_angs, 9);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_point_cos_para,
                                src->cos_point_para_angs, 9);

  TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_gauss_sin_para,
                                src->sin_gauss_para_angs, 9);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_gauss_cos_para,
                                src->cos_gauss_para_angs, 9);

  TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_shape_sin_para,
                                src->sin_shape_para_angs, 9);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_shape_cos_para,
                                src->cos_shape_para_angs, 9);

}

//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_calc_para_angle);

    return UNITY_END();
}