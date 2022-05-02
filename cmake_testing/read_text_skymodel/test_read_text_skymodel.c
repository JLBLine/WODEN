#include <unity.h>
#include <stdlib.h>
#include <math.h>
//
#include "constants.h"
#include "expected_skymodel_outcomes.h"
#include "read_text_skymodel.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

#ifdef DOUBLE_PRECISION
  double TOL_HERE = 1e-15;
#else
  double TOL_HERE = 1e-7;
#endif

void test_read_text_skymodel_SinglePointGood(void) {
  test_read_skymodel_SinglePoint("srclist_singlepoint.txt", POWER_LAW);
}

void test_read_text_skymodel_SinglePointEmpty(void) {
  test_read_skymodel_SinglePoint("srclist_empty_line.txt", POWER_LAW);
}

void test_read_text_skymodel_SinglePointComments(void) {
  test_read_skymodel_SinglePoint("srclist_comment.txt", POWER_LAW);
}

void test_read_text_skymodel_SingleGaussian(void) {
  test_read_skymodel_SingleGaussian("srclist_singlegauss.txt", POWER_LAW);
}

void test_read_text_skymodel_SingleShapelet(void) {
  test_read_skymodel_SingleShapelet("srclist_singleshape.txt", POWER_LAW);
}

void test_read_text_skymodel_ThreeSources(void) {
  test_read_skymodel_ThreeSources("srclist_threesources.txt", POWER_LAW);
}

void test_read_text_skymodel_ThreeComponents(void) {
  test_read_skymodel_ThreeComponents("srclist_threecomponents.txt", POWER_LAW);
}


void test_read_text_skymodel_MultiSourceComponents(void) {
  // // Read in the source catalogue
  source_catalogue_t *raw_srccat = malloc( sizeof(source_catalogue_t) );
  int status=0;
  status = read_text_skymodel("srclist_mulitple_source-components.txt", raw_srccat);

  //Zero means we're happy
  TEST_ASSERT_EQUAL_INT(0, status);
  //
  //Expected overall numbers
  int num_sources = 5;
  int num_shapelets = 4;

  //Expected numbers for all sources
  // int source_indexes[] = {0}
  int n_comps[] = {1, 1, 1, 3, 6};
  int n_points[] = {1, 0, 0, 1, 2};
  int n_gauss[] = {0, 1, 0, 1, 2};
  int n_shapes[] = {0, 0, 1, 1, 2};
  int n_shape_components[] = {0, 0, 4, 4, 8};

  //Everything is a POWER LAW so set lists and curves to zero
  int n_point_lists=0, n_point_curves=0;
  int n_gauss_lists=0, n_gauss_curves=0;
  int n_shape_lists=0, n_shape_curves=0;

  for (int source_index = 0; source_index < num_sources; source_index++) {
    check_single_source_numbers(raw_srccat,
      num_sources,  num_shapelets, n_comps[source_index],
      n_points[source_index], n_point_lists, n_points[source_index], n_point_curves,
      n_gauss[source_index], n_gauss_lists, n_gauss[source_index], n_gauss_curves,
      n_shapes[source_index], n_shape_lists, n_shapes[source_index], n_shape_curves,
      n_shape_components[source_index], source_index);
  }

  //Check that the POINT components are located where expected
  //Index of array if the SOURCE index. Second number is the COMPONENT in that
  //particular SOURCE
  check_single_component_point(raw_srccat->sources[0], 0, POWER_LAW, 0);
  check_single_component_point(raw_srccat->sources[3], 0, POWER_LAW, 0);
  check_single_component_point(raw_srccat->sources[4], 0, POWER_LAW, 0);
  check_single_component_point(raw_srccat->sources[4], 1, POWER_LAW, 1);
  //
  //Check that the GAUSSIAN components are located where expected
  //Index of array if the SOURCE index. Second number is the COMPONENT in that
  //particular SOURCE
  check_single_component_gauss(raw_srccat->sources[1], 0, POWER_LAW, 0);
  check_single_component_gauss(raw_srccat->sources[3], 0, POWER_LAW, 0);
  check_single_component_gauss(raw_srccat->sources[4], 0, POWER_LAW, 0);
  check_single_component_gauss(raw_srccat->sources[4], 1, POWER_LAW, 1);
  //
  //Check that the SHAPELET components are located where expected
  //Index of array if the SOURCE index. Second number is the COMPONENT in that
  //particular SOURCE
  check_single_component_shapelet(raw_srccat->sources[2], 0, POWER_LAW, 0);
  check_single_component_shapelet(raw_srccat->sources[3], 0, POWER_LAW, 0);
  check_single_component_shapelet(raw_srccat->sources[4], 0, POWER_LAW, 0);
  check_single_component_shapelet(raw_srccat->sources[4], 1, POWER_LAW, 1);

}

/*
Test whether a single point source is read in correctly when the SED is
specified by the LINEAR keyword, not the FREQ keyword
Use this to test in some allmost broken srclists
*/
void test_read_text_skymodel_Linear(void)
{
  // Read in the source catalogue
  source_catalogue_t *raw_srccat = malloc( sizeof(source_catalogue_t) );
  int status=0;
  status = read_text_skymodel("srclist_linear.txt", raw_srccat);

  //Zero means we're happy
  TEST_ASSERT_EQUAL_INT(0, status);

  //Expected values
  int num_sources = 1;
  int num_shapelets = 0;
  int n_comps = 1;
  int n_points = 1;
  int n_gauss = 0;
  int n_shapes = 0;
  int n_shape_components = 0;

  //Everything is a POWER LAW so set lists and curves to zero
  int n_point_lists=0, n_point_curves=0;
  int n_gauss_lists=0, n_gauss_curves=0;
  int n_shape_lists=0, n_shape_curves=0;

  //Check these values
  check_single_source_numbers(raw_srccat, num_sources,  num_shapelets,  n_comps,
          n_points, n_point_lists, n_points,  n_point_curves,
          n_gauss, n_gauss_lists, n_gauss,  n_gauss_curves,
          n_shapes, n_shape_lists, n_shapes,  n_shape_curves,
          n_shape_components, 0);

  // source_t source = raw_srccat->sources[0];

  TEST_ASSERT_DOUBLE_WITHIN(1e-15, 15*17.0*DD2R, raw_srccat->sources[0].point_components.ras[0]);
  TEST_ASSERT_DOUBLE_WITHIN(1e-15, -23.4545756877980*DD2R, raw_srccat->sources[0].point_components.decs[0]);
  TEST_ASSERT_DOUBLE_WITHIN(1e-15, 1748492837, raw_srccat->sources[0].point_components.power_ref_freqs[0]);
  TEST_ASSERT_DOUBLE_WITHIN(TOL_HERE, 1.0, raw_srccat->sources[0].point_components.power_ref_stokesI[0]);
  TEST_ASSERT_DOUBLE_WITHIN(TOL_HERE, 2.0, raw_srccat->sources[0].point_components.power_ref_stokesQ[0]);
  TEST_ASSERT_DOUBLE_WITHIN(TOL_HERE, 3.0, raw_srccat->sources[0].point_components.power_ref_stokesU[0]);
  TEST_ASSERT_DOUBLE_WITHIN(TOL_HERE, 4.0, raw_srccat->sources[0].point_components.power_ref_stokesV[0]);
  TEST_ASSERT_DOUBLE_WITHIN(TOL_HERE, -0.198236, raw_srccat->sources[0].point_components.power_SIs[0]);
}

void test_read_text_skymodel_MissingFile(void) {
  // Read in the source catalogue
  source_catalogue_t *raw_srccat = malloc( sizeof(source_catalogue_t) );
  int status=0;
  status = read_text_skymodel("not_a_file.txt", raw_srccat);
  //status=1 means this catalogue read failed, so this fail is a pass
  TEST_ASSERT_EQUAL_INT(1, status);
};

void test_read_text_skymodel_BadSpell(void) {
  // Read in the source catalogue
  source_catalogue_t *raw_srccat = malloc( sizeof(source_catalogue_t) );
  int status=0;
  status = read_text_skymodel("srclist_badspell.txt", raw_srccat);
  //status=1 means this catalogue read failed, so this fail is a pass
  TEST_ASSERT_EQUAL_INT(1, status);
};

void test_read_text_skymodel_BadCoeff(void) {
  // Read in the source catalogue
  source_catalogue_t *raw_srccat = malloc( sizeof(source_catalogue_t) );
  int status=0;
  status = read_text_skymodel("srclist_badcoeff.txt", raw_srccat);
  //status=1 means this catalogue read failed, so this fail is a pass
  TEST_ASSERT_EQUAL_INT(1, status);
};

//Test if the number of components for a SOURCE are missing
void test_read_text_skymodel_NoCompNumbers(void) {
  // Read in the source catalogue
  source_catalogue_t *raw_srccat = malloc( sizeof(source_catalogue_t) );
  int status=0;
  status = read_text_skymodel("srclist_no-comp_numbers.txt", raw_srccat);
  //status=1 means this catalogue read failed, so this fail is a pass
  TEST_ASSERT_EQUAL_INT(1, status);
};
//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_read_text_skymodel_SinglePointGood);
    RUN_TEST(test_read_text_skymodel_SinglePointEmpty);
    RUN_TEST(test_read_text_skymodel_SinglePointComments);
    RUN_TEST(test_read_text_skymodel_SingleGaussian);
    RUN_TEST(test_read_text_skymodel_SingleShapelet);
    RUN_TEST(test_read_text_skymodel_ThreeSources);
    RUN_TEST(test_read_text_skymodel_ThreeComponents);
    RUN_TEST(test_read_text_skymodel_MultiSourceComponents);
    RUN_TEST(test_read_text_skymodel_MissingFile);
    RUN_TEST(test_read_text_skymodel_BadSpell);
    RUN_TEST(test_read_text_skymodel_BadCoeff);
    RUN_TEST(test_read_text_skymodel_NoCompNumbers);
    RUN_TEST(test_read_text_skymodel_Linear);

    return UNITY_END();
}
