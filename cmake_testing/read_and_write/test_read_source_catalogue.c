#include <unity.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "read_and_write.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

#define UNITY_INCLUDE_FLOAT

/*
Check that the counts of number of SOURCEs and COMPONENT types in the
source_catalogue_t matches those provide through argument
*/
void check_single_source_numbers(source_catalogue_t *raw_srccat,
                             int num_sources,  int num_shapelets,  int n_comps,
                             int n_points,  int n_gauss,  int n_shapes,
                             int n_shape_coeffs, int source_index) {

  // printf("%d %d %d %d %d %d %d %d\n",source_index,
  //                                 raw_srccat->num_sources,
  //                                 raw_srccat->num_shapelets,
  //                                 raw_srccat->catsources[source_index].n_comps,
  //                                 raw_srccat->catsources[source_index].n_points,
  //                                 raw_srccat->catsources[source_index].n_gauss,
  //                                 raw_srccat->catsources[source_index].n_shapes,
  //                                 raw_srccat->catsources[source_index].n_shape_coeffs );
  //
  TEST_ASSERT_EQUAL_INT(num_sources, raw_srccat->num_sources);
  TEST_ASSERT_EQUAL_INT(num_shapelets, raw_srccat->num_shapelets);
  TEST_ASSERT_EQUAL_INT(n_comps, raw_srccat->catsources[source_index].n_comps);

  TEST_ASSERT_EQUAL_INT(n_points, raw_srccat->catsources[source_index].n_points);
  TEST_ASSERT_EQUAL_INT(n_gauss, raw_srccat->catsources[source_index].n_gauss);
  TEST_ASSERT_EQUAL_INT(n_shapes, raw_srccat->catsources[source_index].n_shapes);
  TEST_ASSERT_EQUAL_INT(n_shape_coeffs, raw_srccat->catsources[source_index].n_shape_coeffs);
}

/*
Checks that the POINT values for the COMPONENT at the given index are correct
Assumes the catalogue being read in is one of the specific sky models made
for testing in the WODEN/cmake_testing/read_and_write dir
*/
void check_single_component_point(catsource_t catsource, int comp_index) {
  TEST_ASSERT_EQUAL_FLOAT(15*1.0*DD2R, catsource.point_ras[comp_index]);
  TEST_ASSERT_EQUAL_FLOAT(-30.0*DD2R, catsource.point_decs[comp_index]);
  TEST_ASSERT_EQUAL_FLOAT(150e+6, catsource.point_ref_freqs[comp_index]);
  TEST_ASSERT_EQUAL_FLOAT(1.0, catsource.point_ref_stokesI[comp_index]);
  TEST_ASSERT_EQUAL_FLOAT(0.0, catsource.point_ref_stokesQ[comp_index]);
  TEST_ASSERT_EQUAL_FLOAT(0.0, catsource.point_ref_stokesU[comp_index]);
  TEST_ASSERT_EQUAL_FLOAT(0.0, catsource.point_ref_stokesV[comp_index]);
  TEST_ASSERT_EQUAL_FLOAT(-0.8, catsource.point_SIs[comp_index]);
}

/*
Checks that the GAUSSIAN values for the COMPONENT at the given index are correct
Assumes the catalogue being read in is one of the specific sky models made
for testing in the WODEN/cmake_testing/read_and_write dir
*/
void check_single_component_gauss(catsource_t catsource, int comp_index) {
  TEST_ASSERT_EQUAL_FLOAT(15*2.0*DD2R, catsource.gauss_ras[comp_index]);
  TEST_ASSERT_EQUAL_FLOAT(-30.0*DD2R, catsource.gauss_decs[comp_index]);
  TEST_ASSERT_EQUAL_FLOAT(150e+6, catsource.gauss_ref_freqs[comp_index]);
  TEST_ASSERT_EQUAL_FLOAT(2.0, catsource.gauss_ref_stokesI[comp_index]);
  TEST_ASSERT_EQUAL_FLOAT(0.0, catsource.gauss_ref_stokesQ[comp_index]);
  TEST_ASSERT_EQUAL_FLOAT(0.0, catsource.gauss_ref_stokesU[comp_index]);
  TEST_ASSERT_EQUAL_FLOAT(0.0, catsource.gauss_ref_stokesV[comp_index]);
  TEST_ASSERT_EQUAL_FLOAT(-0.8, catsource.gauss_SIs[comp_index]);
  TEST_ASSERT_EQUAL_FLOAT(-10.0*DD2R, catsource.gauss_pas[comp_index]);
  TEST_ASSERT_EQUAL_FLOAT((3.0/60.0)*DD2R,
                          catsource.gauss_majors[comp_index]);
  TEST_ASSERT_EQUAL_FLOAT((6.0/60.0)*DD2R,
                          catsource.gauss_minors[comp_index]);
}

/*
Checks that the SHAPELET values for the COMPONENT at the given index are correct
Assumes the catalogue being read in is one of the specific sky models made
for testing in the WODEN/cmake_testing/read_and_write dir
*/
void check_single_component_shapelet(catsource_t catsource, int comp_index) {
  TEST_ASSERT_EQUAL_FLOAT(15*3.0*DD2R, catsource.shape_ras[comp_index]);
  TEST_ASSERT_EQUAL_FLOAT(20.0*DD2R, catsource.shape_decs[comp_index]);
  TEST_ASSERT_EQUAL_FLOAT(70e+6, catsource.shape_ref_freqs[comp_index]);
  TEST_ASSERT_EQUAL_FLOAT(3.0, catsource.shape_ref_stokesI[comp_index]);
  TEST_ASSERT_EQUAL_FLOAT(0.0, catsource.shape_ref_stokesQ[comp_index]);
  TEST_ASSERT_EQUAL_FLOAT(0.0, catsource.shape_ref_stokesU[comp_index]);
  TEST_ASSERT_EQUAL_FLOAT(0.0, catsource.shape_ref_stokesV[comp_index]);
  TEST_ASSERT_EQUAL_FLOAT(-0.8, catsource.shape_SIs[comp_index]);
  TEST_ASSERT_EQUAL_FLOAT(56.0*DD2R, catsource.shape_pas[comp_index]);
  TEST_ASSERT_EQUAL_FLOAT((7.0/60.0)*DD2R, catsource.shape_majors[comp_index]);
  TEST_ASSERT_EQUAL_FLOAT((5.0/60.0)*DD2R, catsource.shape_minors[comp_index]);

  float expected_n1s[] = {0.0, 14.0, 41.0, 37.0};
  float expected_n2s[] = {0.0, 2.0, -15.0, 7.0};
  float expected_coeffs[] = {0.48255952, -0.18494293, -0.08973978, -0.22137849};
  float expected_param_indexes[] = {comp_index, comp_index, comp_index, comp_index};

  //Increment the array pointer to the start of the parameters for this
  //particular component index (each set of params is 4 long)

  TEST_ASSERT_EQUAL_FLOAT_ARRAY(expected_n1s,
                                catsource.shape_n1s + (4*comp_index), 4);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(expected_n2s,
                                catsource.shape_n2s + (4*comp_index), 4);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(expected_coeffs,
                                catsource.shape_coeffs + (4*comp_index), 4);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(expected_param_indexes,
                                catsource.shape_param_indexes + (4*comp_index), 4);
}

/*
Test whether a single point source is read in correctly
Use this to test in some allmost broken srclists
*/
void test_read_source_catalogue_SinglePoint(char *srclist)
{
  // Read in the source catalogue
  source_catalogue_t *raw_srccat = malloc( sizeof(source_catalogue_t) );
  int status=0;
  status = read_source_catalogue(srclist, raw_srccat);

  //Zero means we're happy
  TEST_ASSERT_EQUAL_INT(0, status);

  //Expected values
  int num_sources = 1;
  int num_shapelets = 0;
  int n_comps = 1;
  int n_points = 1;
  int n_gauss = 0;
  int n_shapes = 0;
  int n_shape_coeffs = 0;

  //Check these values
  check_single_source_numbers(raw_srccat,
                              num_sources,  num_shapelets,  n_comps,
                              n_points,  n_gauss,  n_shapes,
                              n_shape_coeffs, 0);
  //
  // //Check the POINT specific values in the first COMPONENT
  // //of the first SOURCE are correct
  check_single_component_point(raw_srccat->catsources[0], 0);
}

void test_read_source_catalogue_SinglePointGood(void) {
  test_read_source_catalogue_SinglePoint("srclist_singlepoint.txt");
}

void test_read_source_catalogue_SinglePointEmpty(void) {
  test_read_source_catalogue_SinglePoint("srclist_empty_line.txt");
}

void test_read_source_catalogue_SinglePointComments(void) {
  test_read_source_catalogue_SinglePoint("srclist_comment.txt");
}


/*
Test whether a single Gaussian source is read in correctly
*/
void test_read_source_catalogue_SingleGaussian(void)
{
  // Read in the source catalogue
  source_catalogue_t *raw_srccat = malloc( sizeof(source_catalogue_t) );
  int status=0;
  status = read_source_catalogue("srclist_singlegauss.txt", raw_srccat);

  //Zero means we're happy
  TEST_ASSERT_EQUAL_INT(0, status);

  //Expected values
  int num_sources = 1;
  int num_shapelets = 0;
  int n_comps = 1;
  int n_points = 0;
  int n_gauss = 1;
  int n_shapes = 0;
  int n_shape_coeffs = 0;

  //Check these values
  check_single_source_numbers(raw_srccat,
                              num_sources,  num_shapelets,  n_comps,
                              n_points,  n_gauss,  n_shapes,
                              n_shape_coeffs, 0);

  // Check that the GAUSSIAN specific values in the first COMPONENT
  //of the first SOURCE are correct
  check_single_component_gauss(raw_srccat->catsources[0], 0);

}

/*
Test whether a single shapelet source is read in correctly
*/
void test_read_source_catalogue_SingleShapelet(void)
{
  // Read in the source catalogue
  source_catalogue_t *raw_srccat = malloc( sizeof(source_catalogue_t) );
  int status=0;
  status = read_source_catalogue("srclist_singleshape.txt", raw_srccat);

  //Zero means we're happy
  TEST_ASSERT_EQUAL_INT(0, status);

  //Expected values
  int num_sources = 1;
  int num_shapelets = 1;
  int n_comps = 1;
  int n_points = 0;
  int n_gauss = 0;
  int n_shapes = 1;
  int n_shape_coeffs = 4;

  //Check these values
  check_single_source_numbers(raw_srccat,
                              num_sources,  num_shapelets,  n_comps,
                              n_points,  n_gauss,  n_shapes,
                              n_shape_coeffs, 0);

  // Check that the SHAPELET specific values in the first COMPONENT
  //of the first SOURCE are correct
  check_single_component_shapelet(raw_srccat->catsources[0], 0);


}

/*
Test whether three separate SOURCEs, each with a single COMPONENT,
ordered as POINT, GAUSSIAN, SHAPELET, read in correctly
*/
void test_read_source_catalogue_ThreeSources(void) {
  // // Read in the source catalogue
  source_catalogue_t *raw_srccat = malloc( sizeof(source_catalogue_t) );
  int status=0;
  status = read_source_catalogue("srclist_threesources.txt", raw_srccat);

  //Zero means we're happy
  TEST_ASSERT_EQUAL_INT(0, status);
  //
  //Expected overall numbers
  int num_sources = 3;
  int num_shapelets = 1;

  //Expected numbers for all sources
  // int source_indexes[] = {0}
  int n_comps[] = {1, 1, 1};
  int n_points[] = {1, 0, 0};
  int n_gauss[] = {0, 1, 0};
  int n_shapes[] = {0, 0, 1};
  int n_shape_coeffs[] = {0, 0, 4};

  for (int source_index = 0; source_index < num_sources; source_index++) {
    check_single_source_numbers(raw_srccat,
                                num_sources,  num_shapelets,
                                n_comps[source_index], n_points[source_index],
                                n_gauss[source_index],  n_shapes[source_index],
                                n_shape_coeffs[source_index], source_index);
  }

  // Check that the POINT specific values in the first COMPONENT
  //of the first SOURCE are correct
  check_single_component_point(raw_srccat->catsources[0], 0);

  // Check that the GAUSSIAN specific values in the first COMPONENT
  //of the second SOURCE are correct
  check_single_component_gauss(raw_srccat->catsources[1], 0);

  // Check that the SHAPELET specific values in the first COMPONENT
  //of the third SOURCE are correct
  check_single_component_shapelet(raw_srccat->catsources[2], 0);

}

void test_read_source_catalogue_ThreeComponents(void) {
  // // Read in the source catalogue
  source_catalogue_t *raw_srccat = malloc( sizeof(source_catalogue_t) );
  int status=0;
  status = read_source_catalogue("srclist_threecomponents.txt", raw_srccat);

  //Zero means we're happy
  TEST_ASSERT_EQUAL_INT(0, status);

  //Expected overall numbers
  int num_sources = 1;
  int num_shapelets = 1;

  int n_comps = 3;
  int n_points = 1;
  int n_gauss = 1;
  int n_shapes = 1;
  int n_shape_coeffs = 4;

  int source_index = 0;

  check_single_source_numbers(raw_srccat,
                              num_sources,  num_shapelets,
                              n_comps, n_points,
                              n_gauss,  n_shapes,
                              n_shape_coeffs, source_index);

  // Check that the POINT specific values in the first COMPONENT
  //of the first SOURCE are correct
  check_single_component_point(raw_srccat->catsources[0], 0);

  // Check that the GAUSSIAN specific values in the first COMPONENT
  //of the first SOURCE are correct
  check_single_component_gauss(raw_srccat->catsources[0], 0);

  // Check that the SHAPELET specific values in the first COMPONENT
  //of the first SOURCE are correct
  check_single_component_shapelet(raw_srccat->catsources[0], 0);

}

void test_read_source_catalogue_MultiSourceComponents(void) {
  // // Read in the source catalogue
  source_catalogue_t *raw_srccat = malloc( sizeof(source_catalogue_t) );
  int status=0;
  status = read_source_catalogue("srclist_mulitple_source-components.txt", raw_srccat);

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
  int n_shape_coeffs[] = {0, 0, 4, 4, 8};

  for (int source_index = 0; source_index < num_sources; source_index++) {
    check_single_source_numbers(raw_srccat,
                                num_sources,  num_shapelets,
                                n_comps[source_index], n_points[source_index],
                                n_gauss[source_index],  n_shapes[source_index],
                                n_shape_coeffs[source_index], source_index);
  }

  //Check that the POINT components are located where expected
  //Index of array if the SOURCE index. Second number is the COMPONENT in that
  //particular SOURCE
  check_single_component_point(raw_srccat->catsources[0], 0);
  check_single_component_point(raw_srccat->catsources[3], 0);
  check_single_component_point(raw_srccat->catsources[4], 0);
  check_single_component_point(raw_srccat->catsources[4], 1);
  //
  //Check that the GAUSSIAN components are located where expected
  //Index of array if the SOURCE index. Second number is the COMPONENT in that
  //particular SOURCE
  check_single_component_gauss(raw_srccat->catsources[1], 0);
  check_single_component_gauss(raw_srccat->catsources[3], 0);
  check_single_component_gauss(raw_srccat->catsources[4], 0);
  check_single_component_gauss(raw_srccat->catsources[4], 1);
  //
  //Check that the SHAPELET components are located where expected
  //Index of array if the SOURCE index. Second number is the COMPONENT in that
  //particular SOURCE
  check_single_component_shapelet(raw_srccat->catsources[2], 0);
  check_single_component_shapelet(raw_srccat->catsources[3], 0);
  check_single_component_shapelet(raw_srccat->catsources[4], 0);
  check_single_component_shapelet(raw_srccat->catsources[4], 1);

}

void test_read_source_catalogue_MissingFile(void) {
  // Read in the source catalogue
  source_catalogue_t *raw_srccat = malloc( sizeof(source_catalogue_t) );
  int status=0;
  status = read_source_catalogue("not_a_file.txt", raw_srccat);
  //status=1 means this catalogue read failed, so this fail is a pass
  TEST_ASSERT_EQUAL_INT(1, status);
};

void test_read_source_catalogue_BadSpell(void) {
  // Read in the source catalogue
  source_catalogue_t *raw_srccat = malloc( sizeof(source_catalogue_t) );
  int status=0;
  status = read_source_catalogue("srclist_badspell.txt", raw_srccat);
  //status=1 means this catalogue read failed, so this fail is a pass
  TEST_ASSERT_EQUAL_INT(1, status);
};

void test_read_source_catalogue_BadCoeff(void) {
  // Read in the source catalogue
  source_catalogue_t *raw_srccat = malloc( sizeof(source_catalogue_t) );
  int status=0;
  status = read_source_catalogue("srclist_badcoeff.txt", raw_srccat);
  //status=1 means this catalogue read failed, so this fail is a pass
  TEST_ASSERT_EQUAL_INT(1, status);
};

//Test if the number of components for a SOURCE are missing
void test_read_source_catalogue_NoCompNumbers(void) {
  // Read in the source catalogue
  source_catalogue_t *raw_srccat = malloc( sizeof(source_catalogue_t) );
  int status=0;
  status = read_source_catalogue("srclist_no-comp_numbers.txt", raw_srccat);
  //status=1 means this catalogue read failed, so this fail is a pass
  TEST_ASSERT_EQUAL_INT(1, status);
};
//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_read_source_catalogue_SinglePointGood);
    RUN_TEST(test_read_source_catalogue_SinglePointEmpty);
    RUN_TEST(test_read_source_catalogue_SinglePointComments);
    RUN_TEST(test_read_source_catalogue_SingleGaussian);
    RUN_TEST(test_read_source_catalogue_SingleShapelet);
    RUN_TEST(test_read_source_catalogue_ThreeSources);
    RUN_TEST(test_read_source_catalogue_ThreeComponents);
    RUN_TEST(test_read_source_catalogue_MultiSourceComponents);
    RUN_TEST(test_read_source_catalogue_MissingFile);
    RUN_TEST(test_read_source_catalogue_BadSpell);
    RUN_TEST(test_read_source_catalogue_BadCoeff);
    RUN_TEST(test_read_source_catalogue_NoCompNumbers);

    return UNITY_END();
}
