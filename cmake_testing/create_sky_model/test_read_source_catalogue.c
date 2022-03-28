#include <unity.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "create_sky_model.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

#ifdef DOUBLE_PRECISION
  double TOL = 1e-15;
#else
  double TOL = 1e-7;
#endif

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
  //                                 raw_srccat->sources[source_index].n_comps,
  //                                 raw_srccat->sources[source_index].n_points,
  //                                 raw_srccat->sources[source_index].n_gauss,
  //                                 raw_srccat->sources[source_index].n_shapes,
  //                                 raw_srccat->sources[source_index].n_shape_components );
  //
  TEST_ASSERT_EQUAL_INT(num_sources, raw_srccat->num_sources);
  TEST_ASSERT_EQUAL_INT(num_shapelets, raw_srccat->num_shapelets);
  TEST_ASSERT_EQUAL_INT(n_comps, raw_srccat->sources[source_index].n_comps);

  TEST_ASSERT_EQUAL_INT(n_points, raw_srccat->sources[source_index].n_points);
  TEST_ASSERT_EQUAL_INT(n_gauss, raw_srccat->sources[source_index].n_gauss);
  TEST_ASSERT_EQUAL_INT(n_shapes, raw_srccat->sources[source_index].n_shapes);
  TEST_ASSERT_EQUAL_INT(n_shape_coeffs, raw_srccat->sources[source_index].n_shape_coeffs);
}

/*
Checks that the POINT values for the COMPONENT at the given index are correct
Assumes the catalogue being read in is one of the specific sky models made
for testing in the WODEN/cmake_testing/read_and_write dir
*/
void check_single_component_point(source_t source, int comp_index) {
  TEST_ASSERT_DOUBLE_WITHIN(1e-15, 15*1.0*DD2R, source.point_components.ras[comp_index]);
  TEST_ASSERT_DOUBLE_WITHIN(1e-15, -30.0*DD2R, source.point_components.decs[comp_index]);
  TEST_ASSERT_DOUBLE_WITHIN(1e-15, 150e+6, source.point_components.ref_freqs[comp_index]);
  TEST_ASSERT_DOUBLE_WITHIN(TOL, 1.0, source.point_components.ref_stokesI[comp_index]);
  TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, source.point_components.ref_stokesQ[comp_index]);
  TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, source.point_components.ref_stokesU[comp_index]);
  TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, source.point_components.ref_stokesV[comp_index]);
  TEST_ASSERT_DOUBLE_WITHIN(TOL, -0.8, source.point_components.SIs[comp_index]);
}

/*
Checks that the GAUSSIAN values for the COMPONENT at the given index are correct
Assumes the catalogue being read in is one of the specific sky models made
for testing in the WODEN/cmake_testing/read_and_write dir
*/
void check_single_component_gauss(source_t source, int comp_index) {
  TEST_ASSERT_DOUBLE_WITHIN(1e-15, 15*2.0*DD2R, source.gauss_components.ras[comp_index]);
  TEST_ASSERT_DOUBLE_WITHIN(1e-15, -30.0*DD2R, source.gauss_components.decs[comp_index]);
  TEST_ASSERT_DOUBLE_WITHIN(1e-15, 150e+6, source.gauss_components.ref_freqs[comp_index]);
  TEST_ASSERT_DOUBLE_WITHIN(TOL, 2.0, source.gauss_components.ref_stokesI[comp_index]);
  TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, source.gauss_components.ref_stokesQ[comp_index]);
  TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, source.gauss_components.ref_stokesU[comp_index]);
  TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, source.gauss_components.ref_stokesV[comp_index]);
  TEST_ASSERT_DOUBLE_WITHIN(TOL, -0.8, source.gauss_components.SIs[comp_index]);
  TEST_ASSERT_DOUBLE_WITHIN(TOL, -10.0*DD2R, source.gauss_components.pas[comp_index]);
  TEST_ASSERT_DOUBLE_WITHIN(TOL, (3.0/60.0)*DD2R,
                          source.gauss_components.majors[comp_index]);
  TEST_ASSERT_DOUBLE_WITHIN(TOL, (6.0/60.0)*DD2R,
                          source.gauss_components.minors[comp_index]);
}

/*
Checks that the SHAPELET values for the COMPONENT at the given index are correct
Assumes the catalogue being read in is one of the specific sky models made
for testing in the WODEN/cmake_testing/read_and_write dir
*/
void check_single_component_shapelet(source_t source, int comp_index) {
  TEST_ASSERT_DOUBLE_WITHIN(1e-15, 15*3.0*DD2R, source.shape_components.ras[comp_index]);
  TEST_ASSERT_DOUBLE_WITHIN(1e-15, 20.0*DD2R, source.shape_components.decs[comp_index]);
  TEST_ASSERT_DOUBLE_WITHIN(1e-15, 70e+6, source.shape_components.ref_freqs[comp_index]);
  TEST_ASSERT_DOUBLE_WITHIN(TOL, 3.0, source.shape_components.ref_stokesI[comp_index]);
  TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, source.shape_components.ref_stokesQ[comp_index]);
  TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, source.shape_components.ref_stokesU[comp_index]);
  TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, source.shape_components.ref_stokesV[comp_index]);
  TEST_ASSERT_DOUBLE_WITHIN(TOL, -0.8, source.shape_components.SIs[comp_index]);
  TEST_ASSERT_DOUBLE_WITHIN(TOL, 56.0*DD2R, source.shape_components.pas[comp_index]);
  TEST_ASSERT_DOUBLE_WITHIN(TOL, (7.0/60.0)*DD2R, source.shape_components.majors[comp_index]);
  TEST_ASSERT_DOUBLE_WITHIN(TOL, (5.0/60.0)*DD2R, source.shape_components.minors[comp_index]);

  double expected_n1s[] = {0.0, 14.0, 41.0, 37.0};
  double expected_n2s[] = {0.0, 2.0, -15.0, 7.0};
  double expected_coeffs[] = {0.48255952, -0.18494293, -0.08973978, -0.22137849};
  double expected_param_indexes[] = {comp_index, comp_index, comp_index, comp_index};

  for (int param_index = 0; param_index < 4; param_index++) {
    TEST_ASSERT_DOUBLE_WITHIN(TOL, expected_n1s[param_index],
                      source.shape_components.n1s[4*comp_index + param_index]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, expected_n2s[param_index],
                      source.shape_components.n2s[4*comp_index + param_index]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, expected_coeffs[param_index],
                      source.shape_components.shape_coeffs[4*comp_index + param_index]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, expected_param_indexes[param_index],
                      source.shape_components.param_indexes[4*comp_index + param_index]);
  }
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
  int n_shape_components = 0;

  //Check these values
  check_single_source_numbers(raw_srccat,
                              num_sources,  num_shapelets,  n_comps,
                              n_points,  n_gauss,  n_shapes,
                              n_shape_components, 0);
  //
  // //Check the POINT specific values in the first COMPONENT
  // //of the first SOURCE are correct
  check_single_component_point(raw_srccat->sources[0], 0);
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
  int n_shape_components = 0;

  //Check these values
  check_single_source_numbers(raw_srccat,
                              num_sources,  num_shapelets,  n_comps,
                              n_points,  n_gauss,  n_shapes,
                              n_shape_components, 0);

  // Check that the GAUSSIAN specific values in the first COMPONENT
  //of the first SOURCE are correct
  check_single_component_gauss(raw_srccat->sources[0], 0);

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
  int n_shape_components = 4;

  //Check these values
  check_single_source_numbers(raw_srccat,
                              num_sources,  num_shapelets,  n_comps,
                              n_points,  n_gauss,  n_shapes,
                              n_shape_components, 0);

  // Check that the SHAPELET specific values in the first COMPONENT
  //of the first SOURCE are correct
  check_single_component_shapelet(raw_srccat->sources[0], 0);


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
  int n_shape_components[] = {0, 0, 4};

  for (int source_index = 0; source_index < num_sources; source_index++) {
    check_single_source_numbers(raw_srccat,
                                num_sources,  num_shapelets,
                                n_comps[source_index], n_points[source_index],
                                n_gauss[source_index],  n_shapes[source_index],
                                n_shape_components[source_index], source_index);
  }

  // Check that the POINT specific values in the first COMPONENT
  //of the first SOURCE are correct
  check_single_component_point(raw_srccat->sources[0], 0);

  // Check that the GAUSSIAN specific values in the first COMPONENT
  //of the second SOURCE are correct
  check_single_component_gauss(raw_srccat->sources[1], 0);

  // Check that the SHAPELET specific values in the first COMPONENT
  //of the third SOURCE are correct
  check_single_component_shapelet(raw_srccat->sources[2], 0);

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
  int n_shape_components = 4;

  int source_index = 0;

  check_single_source_numbers(raw_srccat,
                              num_sources,  num_shapelets,
                              n_comps, n_points,
                              n_gauss,  n_shapes,
                              n_shape_components, source_index);

  // Check that the POINT specific values in the first COMPONENT
  //of the first SOURCE are correct
  check_single_component_point(raw_srccat->sources[0], 0);

  // Check that the GAUSSIAN specific values in the first COMPONENT
  //of the first SOURCE are correct
  check_single_component_gauss(raw_srccat->sources[0], 0);

  // Check that the SHAPELET specific values in the first COMPONENT
  //of the first SOURCE are correct
  check_single_component_shapelet(raw_srccat->sources[0], 0);

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
  int n_shape_components[] = {0, 0, 4, 4, 8};

  for (int source_index = 0; source_index < num_sources; source_index++) {
    check_single_source_numbers(raw_srccat,
                                num_sources,  num_shapelets,
                                n_comps[source_index], n_points[source_index],
                                n_gauss[source_index],  n_shapes[source_index],
                                n_shape_components[source_index], source_index);
  }

  //Check that the POINT components are located where expected
  //Index of array if the SOURCE index. Second number is the COMPONENT in that
  //particular SOURCE
  check_single_component_point(raw_srccat->sources[0], 0);
  check_single_component_point(raw_srccat->sources[3], 0);
  check_single_component_point(raw_srccat->sources[4], 0);
  check_single_component_point(raw_srccat->sources[4], 1);
  //
  //Check that the GAUSSIAN components are located where expected
  //Index of array if the SOURCE index. Second number is the COMPONENT in that
  //particular SOURCE
  check_single_component_gauss(raw_srccat->sources[1], 0);
  check_single_component_gauss(raw_srccat->sources[3], 0);
  check_single_component_gauss(raw_srccat->sources[4], 0);
  check_single_component_gauss(raw_srccat->sources[4], 1);
  //
  //Check that the SHAPELET components are located where expected
  //Index of array if the SOURCE index. Second number is the COMPONENT in that
  //particular SOURCE
  check_single_component_shapelet(raw_srccat->sources[2], 0);
  check_single_component_shapelet(raw_srccat->sources[3], 0);
  check_single_component_shapelet(raw_srccat->sources[4], 0);
  check_single_component_shapelet(raw_srccat->sources[4], 1);

}

/*
Test whether a single point source is read in correctly when the SED is
specified by the LINEAR keyword, not the FREQ keyword
Use this to test in some allmost broken srclists
*/
void test_read_source_catalogue_Linear(void)
{
  // Read in the source catalogue
  source_catalogue_t *raw_srccat = malloc( sizeof(source_catalogue_t) );
  int status=0;
  status = read_source_catalogue("srclist_linear.txt", raw_srccat);

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

  //Check these values
  check_single_source_numbers(raw_srccat,
                              num_sources,  num_shapelets,  n_comps,
                              n_points,  n_gauss,  n_shapes,
                              n_shape_components, 0);

  // source_t source = raw_srccat->sources[0];

  TEST_ASSERT_DOUBLE_WITHIN(1e-15, 15*17.0*DD2R, raw_srccat->sources[0].point_components.ras[0]);
  TEST_ASSERT_DOUBLE_WITHIN(1e-15, -23.4545756877980*DD2R, raw_srccat->sources[0].point_components.decs[0]);
  TEST_ASSERT_DOUBLE_WITHIN(1e-15, 1748492837, raw_srccat->sources[0].point_components.ref_freqs[0]);
  TEST_ASSERT_DOUBLE_WITHIN(TOL, 1.0, raw_srccat->sources[0].point_components.ref_stokesI[0]);
  TEST_ASSERT_DOUBLE_WITHIN(TOL, 2.0, raw_srccat->sources[0].point_components.ref_stokesQ[0]);
  TEST_ASSERT_DOUBLE_WITHIN(TOL, 3.0, raw_srccat->sources[0].point_components.ref_stokesU[0]);
  TEST_ASSERT_DOUBLE_WITHIN(TOL, 4.0, raw_srccat->sources[0].point_components.ref_stokesV[0]);
  TEST_ASSERT_DOUBLE_WITHIN(TOL, -0.198236, raw_srccat->sources[0].point_components.SIs[0]);
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
    RUN_TEST(test_read_source_catalogue_Linear);

    return UNITY_END();
}
