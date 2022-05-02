#include <unity.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "woden_settings.h"
#include "expected_skymodel_outcomes.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

// #ifdef DOUBLE_PRECISION
//   double TOL = 1e-15;
// #else
//   double TOL = 1e-7;
// #endif

#ifdef DOUBLE_PRECISION
  double TOL_HERE = 1e-15;
#else
  double TOL_HERE = 1e-7;
#endif

void test_read_yaml_skymodel_SinglePointEmpty(void) {
  test_read_skymodel_SinglePoint("srclist_empty_line.yaml", POWER_LAW);
}

void test_read_yaml_skymodel_SinglePointComments(void) {
  test_read_skymodel_SinglePoint("srclist_comment.yaml", POWER_LAW);
}

//POWER LAW flux behaviour tests------------------------------------------------
void test_read_yaml_skymodel_SinglePointPower(void) {
  test_read_skymodel_SinglePoint("srclist_singlepoint_power.yaml", POWER_LAW);
}

void test_read_yaml_skymodel_SingleGaussianPower(void) {
  test_read_skymodel_SingleGaussian("srclist_singlegauss_power.yaml", POWER_LAW);
}

void test_read_yaml_skymodel_SingleShapeletPower(void) {
  test_read_skymodel_SingleShapelet("srclist_singleshape_power.yaml", POWER_LAW);
}

void test_read_yaml_skymodel_ThreeSourcesPower(void) {
  test_read_skymodel_ThreeSources("srclist_threesources_power.yaml", POWER_LAW);
}

void test_read_yaml_skymodel_ThreeComponentsPower(void) {
  test_read_skymodel_ThreeComponents("srclist_threecomponents_power.yaml", POWER_LAW);
}

//LIST flux behaviour tests------------------------------------------------
void test_read_yaml_skymodel_SinglePointList(void) {
  test_read_skymodel_SinglePoint("srclist_singlepoint_list.yaml", LIST);
}

void test_read_yaml_skymodel_SingleGaussianList(void) {
  test_read_skymodel_SingleGaussian("srclist_singlegauss_list.yaml", LIST);
}

void test_read_yaml_skymodel_SingleShapeletList(void) {
  test_read_skymodel_SingleShapelet("srclist_singleshape_list.yaml", LIST);
}

void test_read_yaml_skymodel_ThreeSourcesList(void) {
  test_read_skymodel_ThreeSources("srclist_threesources_list.yaml", LIST);
}

void test_read_yaml_skymodel_ThreeComponentsList(void) {
  test_read_skymodel_ThreeComponents("srclist_threecomponents_list.yaml", LIST);
}

//Curve power law behaviour-----------------------------------------------------
void test_read_yaml_skymodel_SinglePointCurve(void) {
  test_read_skymodel_SinglePoint("srclist_singlepoint_curve.yaml", CURVED_POWER_LAW);
}

void test_read_yaml_skymodel_SingleGaussianCurve(void) {
  test_read_skymodel_SingleGaussian("srclist_singlegauss_curve.yaml", CURVED_POWER_LAW);
}

void test_read_yaml_skymodel_SingleShapeletCurve(void) {
  test_read_skymodel_SingleShapelet("srclist_singleshape_curve.yaml", CURVED_POWER_LAW);
}

void test_read_yaml_skymodel_ThreeSourcesCurve(void) {
  test_read_skymodel_ThreeSources("srclist_threesources_curve.yaml", CURVED_POWER_LAW);
}

void test_read_yaml_skymodel_ThreeComponentsCurve(void) {
  test_read_skymodel_ThreeComponents("srclist_threecomponents_curve.yaml", CURVED_POWER_LAW);
}





void test_read_text_skymodel_MultiSourceComponents(void) {
  // // Read in the source catalogue
  source_catalogue_t *raw_srccat = malloc( sizeof(source_catalogue_t) );
  int status=0;
  status = read_skymodel("srclist_mulitple_source-components.yaml", raw_srccat);

  //Zero means we're happy
  TEST_ASSERT_EQUAL_INT(0, status);

  //Expected overall numbers
  int num_sources = 10;
  int num_shapelets = 9;

  //Expected numbers for all sources. These are obviously hardcoded to match
  // what exists in srclist_mulitple_source-components.yaml
  int n_comps[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 18};
  int n_points[] = {1, 1, 1, 0, 0, 0, 0, 0, 0, 6};
  int n_gauss[] = {0, 0, 0, 1, 1, 1, 0, 0, 0, 6};
  int n_shapes[] = {0, 0, 0, 0, 0, 0, 1, 1, 1, 6};
  int n_shape_components[] = {0, 0, 0, 0, 0, 0, 4, 4, 4, 24};
  //
  int n_point_powers[] = {1, 0, 0, 0, 0, 0, 0, 0, 0, 2};
  int n_point_curves[] = {0, 0, 1, 0, 0, 0, 0, 0, 0, 2};
  int n_point_lists[] = {0, 1, 0, 0, 0, 0, 0, 0, 0, 2};

  int n_gauss_powers[] = {0, 0, 0, 1, 0, 0, 0, 0, 0, 2};
  int n_gauss_curves[] = {0, 0, 0, 0, 1, 0, 0, 0, 0, 2};
  int n_gauss_lists[] = {0, 0, 0, 0, 0, 1, 0, 0, 0, 2};

  int n_shape_powers[] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 2};
  int n_shape_curves[] = {0, 0, 0, 0, 0, 0, 1, 0, 0, 2};
  int n_shape_lists[] = {0, 0, 0, 0, 0, 0, 0, 1, 0, 2};

  for (int source_index = 0; source_index < num_sources; source_index++) {
    check_single_source_numbers(raw_srccat,
      num_sources,  num_shapelets, n_comps[source_index],
      n_points[source_index], n_point_lists[source_index],
      n_point_powers[source_index], n_point_curves[source_index],
      n_gauss[source_index], n_gauss_lists[source_index],
      n_gauss_powers[source_index], n_gauss_curves[source_index],
      n_shapes[source_index], n_shape_lists[source_index],
      n_shape_powers[source_index], n_shape_curves[source_index],
      n_shape_components[source_index], source_index);
  }
  //
  // //Check that the POINT components are located where expected
  // //Index of array if the SOURCE index. Second number is the COMPONENT in that
  // //particular SOURCE, followed by flux type, then flux index in that COMPONENT

  check_single_component_point(raw_srccat->sources[0], 0, POWER_LAW, 0);
  check_single_component_point(raw_srccat->sources[1], 0, LIST, 0);
  check_single_component_point(raw_srccat->sources[2], 0, CURVED_POWER_LAW, 0);

  check_single_component_gauss(raw_srccat->sources[3], 0, POWER_LAW, 0);
  check_single_component_gauss(raw_srccat->sources[4], 0, CURVED_POWER_LAW, 0);
  check_single_component_gauss(raw_srccat->sources[5], 0, LIST, 0);

  check_single_component_shapelet(raw_srccat->sources[6], 0, CURVED_POWER_LAW, 0);
  check_single_component_shapelet(raw_srccat->sources[7], 0, LIST, 0);
  check_single_component_shapelet(raw_srccat->sources[8], 0, POWER_LAW, 0);

  // I've put two of each component type in the final SOURCE, so check
  //the whole shebang gets read in correctly
  //I stuffed up and put the two CURVED_POWER_LAW gaussian components at
  //the end of the srclist, but actually having a changed up order is a good
  //test

  check_single_component_point(raw_srccat->sources[9], 0, POWER_LAW, 0);
  check_single_component_point(raw_srccat->sources[9], 1, LIST, 0);
  check_single_component_point(raw_srccat->sources[9], 2, CURVED_POWER_LAW, 0);
  check_single_component_point(raw_srccat->sources[9], 3, POWER_LAW, 1);
  check_single_component_point(raw_srccat->sources[9], 4, LIST, 1);
  check_single_component_point(raw_srccat->sources[9], 5, CURVED_POWER_LAW, 1);


  check_single_component_gauss(raw_srccat->sources[9], 0, POWER_LAW, 0);
  check_single_component_gauss(raw_srccat->sources[9], 1, LIST, 0);
  check_single_component_gauss(raw_srccat->sources[9], 2, POWER_LAW, 1);
  check_single_component_gauss(raw_srccat->sources[9], 3, LIST, 1);
  check_single_component_gauss(raw_srccat->sources[9], 4, CURVED_POWER_LAW, 0);
  check_single_component_gauss(raw_srccat->sources[9], 5, CURVED_POWER_LAW, 1);

  check_single_component_shapelet(raw_srccat->sources[9], 0, CURVED_POWER_LAW, 0);
  check_single_component_shapelet(raw_srccat->sources[9], 1, LIST, 0);
  check_single_component_shapelet(raw_srccat->sources[9], 2, POWER_LAW, 0);
  check_single_component_shapelet(raw_srccat->sources[9], 3, CURVED_POWER_LAW, 1);
  check_single_component_shapelet(raw_srccat->sources[9], 4, LIST, 1);
  check_single_component_shapelet(raw_srccat->sources[9], 5, POWER_LAW, 1);

}

void test_read_yaml_skymodel_MissingFile(void) {
  // Read in the source catalogue
  source_catalogue_t *raw_srccat = malloc( sizeof(source_catalogue_t) );
  int status=0;
  status = read_yaml_skymodel("not_a_file.yaml", raw_srccat);
  //status=1 means this catalogue read failed, so this fail is a pass
  TEST_ASSERT_EQUAL_INT(1, status);
};










//Run test using unity
int main(void)
{
    UNITY_BEGIN();


    RUN_TEST(test_read_yaml_skymodel_SinglePointEmpty);
    RUN_TEST(test_read_yaml_skymodel_SinglePointComments);
    RUN_TEST(test_read_yaml_skymodel_MissingFile);
    // // //RUN_TEST(test_read_source_catalogue_BadSpell);

    RUN_TEST(test_read_yaml_skymodel_SinglePointPower);
    RUN_TEST(test_read_yaml_skymodel_SingleGaussianPower);
    RUN_TEST(test_read_yaml_skymodel_SingleShapeletPower);
    RUN_TEST(test_read_yaml_skymodel_ThreeSourcesPower);
    RUN_TEST(test_read_yaml_skymodel_ThreeComponentsPower);

    RUN_TEST(test_read_yaml_skymodel_SinglePointList);
    RUN_TEST(test_read_yaml_skymodel_SingleGaussianList);
    RUN_TEST(test_read_yaml_skymodel_SingleShapeletList);
    RUN_TEST(test_read_yaml_skymodel_ThreeSourcesList);
    RUN_TEST(test_read_yaml_skymodel_ThreeComponentsList);

    RUN_TEST(test_read_yaml_skymodel_SinglePointCurve);
    RUN_TEST(test_read_yaml_skymodel_SingleGaussianCurve);
    RUN_TEST(test_read_yaml_skymodel_SingleShapeletCurve);
    RUN_TEST(test_read_yaml_skymodel_ThreeSourcesCurve);
    RUN_TEST(test_read_yaml_skymodel_ThreeComponentsCurve);

    RUN_TEST(test_read_text_skymodel_MultiSourceComponents);

    return UNITY_END();
}
