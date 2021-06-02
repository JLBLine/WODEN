#include <unity.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "chunk_sky_model.h"
#include "woden_struct_defs.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

#define UNITY_INCLUDE_FLOAT


/*
`create_sky_model:horizon_test` is designed to be used by
`create_sky_model:crop_sky_model`.

`horizon_test` takes information on a single COMPONENT of a SOURCE in a sky
model and tests whether it is above or below the horizon. Depending on
whether we are cropping the sky model by SOURCE or by COMPONENT, it updates
various counters to effec the sky model cropping. SHAPELETs are complicating
factors as a single position can match multiple basis function parameters (of
any length) so `horizon_test` does some logic to count how many SHAPELET
parameters are being retained.
*/


/*
Test when a POINT/GAUSSIAN/SHAPELET type component is above horizon and cropping by
SOURCE
*/
void test_horizon_test_CropSourceAbove(void)
{

  //Stick a component above zenith which is not a SHAPELET
  double za = 0;
  int num_shape_coeff_retained = 0;
  int num_shape_coeff_component = 0;
  float *shape_param_indexes;
  int shape = 0;
  //We want to crop by whole SOURCE not by COMPONENT
  e_sky_crop sky_crop_type = CROP_SOURCES;
  //All components are currently above horizon
  e_horizon all_comps_above_horizon = ABOVE;
  //This number is used when sky_crop_type = CROP_COMPONENTS
  int num_comp_retained = 0;

  horizon_test(za, sky_crop_type,
       &all_comps_above_horizon, &num_comp_retained,
       &num_shape_coeff_retained, num_shape_coeff_component,
       shape_param_indexes, shape);

  //We should get that everything is still above horizon
  TEST_ASSERT_EQUAL_INT(ABOVE, all_comps_above_horizon);
  //SHAPELET stuff should still be zero
  TEST_ASSERT_EQUAL_INT(0, num_shape_coeff_retained);
  TEST_ASSERT_EQUAL_INT(0, num_shape_coeff_component);

}

/*
Test when a POINT/GAUSSIAN/SHAPELET type component is below horizon and cropping by
SOURCE
*/
void test_horizon_test_CropSourceBelow(void)
{

  //Stick a component below horizon which is not a SHAPELET
  double za = M_PI;
  int num_shape_coeff_retained = 0;
  int num_shape_coeff_component = 0;
  float *shape_param_indexes;
  int shape = 0;
  //We want to crop by whole SOURCE not by COMPONENT
  e_sky_crop sky_crop_type = CROP_SOURCES;
  //All components are currently above horizon
  e_horizon all_comps_above_horizon = ABOVE;
  //This number is used when sky_crop_type = CROP_COMPONENTS
  int num_comp_retained = 0;

  horizon_test(za, sky_crop_type,
       &all_comps_above_horizon, &num_comp_retained,
       &num_shape_coeff_retained, num_shape_coeff_component,
       shape_param_indexes, shape);

  //We should get that everything is still above horizon
  TEST_ASSERT_EQUAL_INT(BELOW, all_comps_above_horizon);
  //SHAPELET stuff should still be zero
  TEST_ASSERT_EQUAL_INT(0, num_shape_coeff_retained);
  TEST_ASSERT_EQUAL_INT(0, num_shape_coeff_component);

}

/*
Test when a POINT/GAUSSIAN type component is above horizon and cropping by
COMPONENT
*/
void test_horizon_test_PointCropComponentAbove(void)
{

  //Stick a component above zenith which is not a SHAPELET
  double za = 0;
  int num_shape_coeff_retained = 0;
  int num_shape_coeff_component = 0;
  float *shape_param_indexes;
  int shape = 0;
  //We want to crop by COMPONENT not whole SOURCE
  e_sky_crop sky_crop_type = CROP_COMPONENTS;
  //Say in cropping in this SOURCE so far has found 1 COMPONENT
  //above the horizon, we should have these values
  int num_comp_retained = 1;
  e_horizon all_comps_above_horizon = ABOVE;

  horizon_test(za, sky_crop_type,
       &all_comps_above_horizon, &num_comp_retained,
       &num_shape_coeff_retained, num_shape_coeff_component,
       shape_param_indexes, shape);

  //We should get that everything is still above horizon
  TEST_ASSERT_EQUAL_INT(ABOVE, all_comps_above_horizon);
  //We should increment the number of retained COMPONENTs
  TEST_ASSERT_EQUAL_INT(2, num_comp_retained);
  //SHAPELET stuff should still be zero
  TEST_ASSERT_EQUAL_INT(0, num_shape_coeff_retained);
  TEST_ASSERT_EQUAL_INT(0, num_shape_coeff_component);

}

/*
Test when a POINT/GAUSSIAN type component is below horizon and cropping by
COMPONENT
*/
void test_horizon_test_PointCropComponentBelow(void)
{

  //Stick a component above zenith which is not a SHAPELET
  double za = M_PI;
  int num_shape_coeff_retained = 0;
  int num_shape_coeff_component = 0;
  float *shape_param_indexes;
  int shape = 0;
  //We want to crop by COMPONENT not whole SOURCE
  e_sky_crop sky_crop_type = CROP_COMPONENTS;
  //Say in cropping in this SOURCE so far has found 1 COMPONENT
  //above the horizon, we should have these values
  int num_comp_retained = 1;
  e_horizon all_comps_above_horizon = ABOVE;

  horizon_test(za, sky_crop_type,
       &all_comps_above_horizon, &num_comp_retained,
       &num_shape_coeff_retained, num_shape_coeff_component,
       shape_param_indexes, shape);

  //Number of retained COMPONENTs should not have incremented
  TEST_ASSERT_EQUAL_INT(1, num_comp_retained);
  //SHAPELET stuff should still be zero
  TEST_ASSERT_EQUAL_INT(0, num_shape_coeff_retained);
  TEST_ASSERT_EQUAL_INT(0, num_shape_coeff_component);

}

/*
Test that when cropping by COMPONENT with multiple SHAPELET components in a
SOURCE works. We'll test with 3 COMPONENTs. This function will test each
COMPONENT index (shape) against expected outcomes
*/

void test_horizon_test_ShapeletCropComponent(double za, int shape,
                                             int num_shape_coeff_retained,
                                             int num_comp_retained,
                                             int out_num_shape_coeff_retained,
                                             int out_num_comp_retained) {

  //The param index array maps shapelet basis functions to their respective
  //COMPONENT info (like ra,dec,flux etc)
  float shape_param_indexes[] = {0, 0, 1, 1, 1, 1, 2};
  int num_shape_coeff_component = 7;
  //We want to crop by COMPONENT not whole SOURCE
  e_sky_crop sky_crop_type = CROP_COMPONENTS;
  e_horizon all_comps_above_horizon = ABOVE;

  //Run function for testing
  horizon_test(za, sky_crop_type,
       &all_comps_above_horizon, &num_comp_retained,
       &num_shape_coeff_retained, num_shape_coeff_component,
       shape_param_indexes, shape);

  //Test outcomes are as expected
  TEST_ASSERT_EQUAL_INT(out_num_shape_coeff_retained, num_shape_coeff_retained);
  TEST_ASSERT_EQUAL_INT(out_num_comp_retained, num_comp_retained);
}

/*
Test when SHAPELET COMPONENT 0 is above horizon
*/
void test_horizon_test_ShapeletCropComponent0Above(void){

  //Initial conditions
  double za = 0.0;
  int shape = 0;
  int num_shape_coeff_retained = 0;
  int num_comp_retained = 0;
  //Outcomes based on conditions and what is in shape_param_indexes
  int out_num_shape_coeff_retained = 2;
  int out_num_comp_retained = 1;

  test_horizon_test_ShapeletCropComponent(za, shape,
                                          num_shape_coeff_retained,
                                          num_comp_retained,
                                          out_num_shape_coeff_retained,
                                          out_num_comp_retained);
}

/*
Test when SHAPELET COMPONENT 0 is below horizon
*/
void test_horizon_test_ShapeletCropComponent0Below(void){

  //Initial conditions
  double za = M_PI;
  int shape = 0;
  int num_shape_coeff_retained = 0;
  int num_comp_retained = 0;
  //Outcomes based on conditions and what is in shape_param_indexes
  int out_num_shape_coeff_retained = 0;
  int out_num_comp_retained = 0;

  test_horizon_test_ShapeletCropComponent(za, shape,
                                          num_shape_coeff_retained,
                                          num_comp_retained,
                                          out_num_shape_coeff_retained,
                                          out_num_comp_retained);
}

/*
Test when SHAPELET COMPONENT 1 is above horizon
*/
void test_horizon_test_ShapeletCropComponent1Above(void){

  //Initial conditions
  double za = 0.0;
  int shape = 1;
  //Let's say COMPONENT 0 was retained
  int num_shape_coeff_retained = 2;
  int num_comp_retained = 1;
  //Outcomes based on conditions and what is in shape_param_indexes
  //Should get 4 extra coeffs retained, and one COMPONENT
  int out_num_shape_coeff_retained = 6;
  int out_num_comp_retained = 2;

  test_horizon_test_ShapeletCropComponent(za, shape,
                                          num_shape_coeff_retained,
                                          num_comp_retained,
                                          out_num_shape_coeff_retained,
                                          out_num_comp_retained);
}

/*
Test when SHAPELET COMPONENT 0 is below horizon
*/
void test_horizon_test_ShapeletCropComponent1Below(void){

  //Initial conditions
  double za = M_PI;
  int shape = 1;
  //Let's say COMPONENT 0 was retained
  int num_shape_coeff_retained = 2;
  int num_comp_retained = 1;
  //Nothing should change as below horizon
  int out_num_shape_coeff_retained = 2;
  int out_num_comp_retained = 1;

  test_horizon_test_ShapeletCropComponent(za, shape,
                                          num_shape_coeff_retained,
                                          num_comp_retained,
                                          out_num_shape_coeff_retained,
                                          out_num_comp_retained);
}

/*
Test when SHAPELET COMPONENT 2 is above horizon
*/
void test_horizon_test_ShapeletCropComponent2Above(void){

  //Initial conditions
  double za = 0.0;
  int shape = 2;
  //Let's say COMPONENT 0 was retained
  int num_shape_coeff_retained = 2;
  int num_comp_retained = 1;
  //Outcomes based on conditions and what is in shape_param_indexes
  //Should get 4 extra coeffs retained, and one COMPONENT
  int out_num_shape_coeff_retained = 3;
  int out_num_comp_retained = 2;

  test_horizon_test_ShapeletCropComponent(za, shape,
                                          num_shape_coeff_retained,
                                          num_comp_retained,
                                          out_num_shape_coeff_retained,
                                          out_num_comp_retained);
}

/*
Test when SHAPELET COMPONENT 2 is below horizon
*/
void test_horizon_test_ShapeletCropComponent2Below(void){

  //Initial conditions
  double za = M_PI;
  int shape = 2;
  //Let's say COMPONENT 0 was retained
  int num_shape_coeff_retained = 2;
  int num_comp_retained = 1;
  //Nothing should change as below horizon
  int out_num_shape_coeff_retained = 2;
  int out_num_comp_retained = 1;

  test_horizon_test_ShapeletCropComponent(za, shape,
                                          num_shape_coeff_retained,
                                          num_comp_retained,
                                          out_num_shape_coeff_retained,
                                          out_num_comp_retained);
}

//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_horizon_test_CropSourceAbove);
    RUN_TEST(test_horizon_test_CropSourceBelow);
    RUN_TEST(test_horizon_test_PointCropComponentAbove);
    RUN_TEST(test_horizon_test_PointCropComponentBelow);
    RUN_TEST(test_horizon_test_ShapeletCropComponent0Above);
    RUN_TEST(test_horizon_test_ShapeletCropComponent0Below);
    RUN_TEST(test_horizon_test_ShapeletCropComponent1Above);
    RUN_TEST(test_horizon_test_ShapeletCropComponent1Below);
    RUN_TEST(test_horizon_test_ShapeletCropComponent2Above);
    RUN_TEST(test_horizon_test_ShapeletCropComponent2Below);

    return UNITY_END();
}
