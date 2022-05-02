#include <unity.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "create_sky_model.h"
#include "expected_skymodel_outcomes.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */



void test_read_skymodel_MissingTxt(void) {
  // Read in the source catalogue
  source_catalogue_t *raw_srccat = malloc( sizeof(source_catalogue_t) );
  int status=0;
  status = read_skymodel("not_a_file.txt", raw_srccat);
  //status=1 means this catalogue read failed, so this fail is a pass
  TEST_ASSERT_EQUAL_INT(1, status);
};

void test_read_skymodel_MissingYaml(void) {
  // Read in the source catalogue
  source_catalogue_t *raw_srccat = malloc( sizeof(source_catalogue_t) );
  int status=0;
  status = read_skymodel("not_a_file.yaml", raw_srccat);
  //status=1 means this catalogue read failed, so this fail is a pass
  TEST_ASSERT_EQUAL_INT(1, status);
};

void test_read_skymodel_NotYamlTxt(void) {
  // Read in the source catalogue
  source_catalogue_t *raw_srccat = malloc( sizeof(source_catalogue_t) );
  int status=0;
  status = read_skymodel("what_are_you_trying_to_pull.dat", raw_srccat);
  //status=1 means this catalogue read failed, so this fail is a pass
  TEST_ASSERT_EQUAL_INT(1, status);
};

//Just check that a single point source YAML file reads in
//Full testing for the yaml reading function is handled in another test suite
void test_read_skymodel_PointPowerYaml(void) {
  // Read in the source catalogue
  source_catalogue_t *raw_srccat = malloc( sizeof(source_catalogue_t) );
  int status=0;
  status = read_skymodel("srclist_singlepoint.yaml", raw_srccat);

  TEST_ASSERT_EQUAL_INT(0, status);

  int num_sources = 1;
  int num_shapelets = 0;
  int n_comps = 1;
  int n_points = 1;
  int n_gauss = 0;
  int n_shapes = 0;
  int n_shape_components = 0;

  int n_point_lists = 0;
  int n_point_powers = 1;
  int n_point_curves = 0;
  int n_gauss_lists = 0;
  int n_gauss_powers = 0;
  int n_gauss_curves = 0;
  int n_shape_lists = 0;
  int n_shape_powers = 0;
  int n_shape_curves = 0;

  check_single_source_numbers(raw_srccat, num_sources,  num_shapelets,  n_comps,
          n_points, n_point_lists, n_point_powers,  n_point_curves,
          n_gauss, n_gauss_lists, n_gauss_powers,  n_gauss_curves,
          n_shapes, n_shape_lists, n_shape_powers,  n_shape_curves,
          n_shape_components, 0);

  // free_source_catalogue(raw_srccat);
};

//Just check that a single point source YAML file reads in
//Full testing for the yaml reading function is handled in another test suite
void test_read_skymodel_PointPowerTxt(void) {
  // Read in the source catalogue
  source_catalogue_t *raw_srccat = malloc( sizeof(source_catalogue_t) );
  int status=0;
  status = read_skymodel("srclist_singlepoint_power.txt", raw_srccat);

  TEST_ASSERT_EQUAL_INT(0, status);

  int num_sources = 1;
  int num_shapelets = 0;
  int n_comps = 1;
  int n_points = 1;
  int n_gauss = 0;
  int n_shapes = 0;
  int n_shape_components = 0;

  int n_point_lists = 0;
  int n_point_powers = 1;
  int n_point_curves = 0;
  int n_gauss_lists = 0;
  int n_gauss_powers = 0;
  int n_gauss_curves = 0;
  int n_shape_lists = 0;
  int n_shape_powers = 0;
  int n_shape_curves = 0;

  check_single_source_numbers(raw_srccat, num_sources,  num_shapelets,  n_comps,
          n_points, n_point_lists, n_point_powers,  n_point_curves,
          n_gauss, n_gauss_lists, n_gauss_powers,  n_gauss_curves,
          n_shapes, n_shape_lists, n_shape_powers,  n_shape_curves,
          n_shape_components, 0);

  // free_source_catalogue(raw_srccat);
};

//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_read_skymodel_MissingTxt);
    RUN_TEST(test_read_skymodel_MissingYaml);
    RUN_TEST(test_read_skymodel_NotYamlTxt);
    RUN_TEST(test_read_skymodel_PointPowerYaml);
    RUN_TEST(test_read_skymodel_PointPowerTxt);

    return UNITY_END();
}
