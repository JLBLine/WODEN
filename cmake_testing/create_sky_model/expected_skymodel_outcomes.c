#include <unity.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "constants.h"
#include "expected_skymodel_outcomes.h"

#ifdef DOUBLE_PRECISION
  double TOL = 1e-15;
#else
  double TOL = 5e-5;
#endif


double POINT0_RA = 15*DD2R;
double POINT0_DEC = -30.0*DD2R;
//Power law / curved power law expected outcomes
double POINT0_POW_FREQ = 150e+6;
double POINT0_POW_I = 1.0;
double POINT0_POW_Q = 0.0;
double POINT0_POW_U = 0.0;
double POINT0_POW_V = 0.0;
double POINT0_POW_SI = -0.8;
double POINT0_CURVE_Q = 0.2;

//List flux type outcomes
double POINT0_LIST_FREQ[] = {120e+6, 170e+6, 180e+6, 190e+6};
double POINT0_LIST_I[] = {1., 5., 10., 4.};
double POINT0_LIST_Q[] = {-2., 1., 0., 0.};
double POINT0_LIST_U[] = {0., 2., 0., 3.};
double POINT0_LIST_V[] = {0., 3., 0, 0.};
int POINT0_NUM_LIST_ENTRIES = 4;

//First gaussian component values-----------------------------------------------
double GAUSS0_RA = 15*2.0*DD2R;
double GAUSS0_DEC = -30.0*DD2R;
double GAUSS0_MAJ = (3.0/60.0)*DD2R;
double GAUSS0_MIN = (6.0/60.0)*DD2R;
double GAUSS0_PA = -10.0*DD2R;

//Power law / curved power law expected outcomes
double GAUSS0_POW_FREQ = 150e+6;
double GAUSS0_POW_I = 2.0;
double GAUSS0_POW_Q = 0.0;
double GAUSS0_POW_U = 0.0;
double GAUSS0_POW_V = 0.0;
double GAUSS0_POW_SI = -0.8;
double GAUSS0_CURVE_Q = -0.3;

//List flux type outcomes
double GAUSS0_LIST_FREQ[] = {100e+6, 150e+6, 220e+6};
double GAUSS0_LIST_I[] = {5., 4., 10.};
double GAUSS0_LIST_Q[] = {3., 0., 3.};
double GAUSS0_LIST_U[] = {0., -5., 4.};
double GAUSS0_LIST_V[] = {0., 0., -1.0};
int GAUSS0_NUM_LIST_ENTRIES = 3;

//First shapelet component values-----------------------------------------------
double SHAPE0_RA = 15*3.0*DD2R;
double SHAPE0_DEC = 20.0*DD2R;
double SHAPE0_MAJ = (7.0/60.0)*DD2R;
double SHAPE0_MIN = (5.0/60.0)*DD2R;
double SHAPE0_PA = 56.0*DD2R;
double SHAPE0_N1S[] = {0.0, 14.0, 41.0, 37.0};
double SHAPE0_N2S[] = {0.0, 2.0, -15.0, 7.0};
double SHAPE0_COEFFS[] = {0.48255952, -0.18494293, -0.08973978, -0.22137849};

//Power law / curved power law expected outcomes
double SHAPE0_POW_FREQ = 70e+6;
double SHAPE0_POW_I = 3.0;
double SHAPE0_POW_Q = 0.0;
double SHAPE0_POW_U = 0.0;
double SHAPE0_POW_V = 0.0;
double SHAPE0_POW_SI = -0.8;
double SHAPE0_CURVE_Q = 0.3245;

//List flux type outcomes
double SHAPE0_LIST_FREQ[] = {110e+6, 120e+6, 130e+6, 150e+6, 180e+6};
double SHAPE0_LIST_I[] = {4.456, 4.75, 3.45, 532.897, 4.0};
double SHAPE0_LIST_Q[] = {902.3, 0.0, 54.786, 0.0, 0.0 };
double SHAPE0_LIST_U[] = {234.234, -432.987, 0.0, 0.0, 0.0};
double SHAPE0_LIST_V[] = {-13.234, 0.0, 0.0, -9.3824, 0.0};

int SHAPE0_NUM_LIST_ENTRIES = 5;

/*
Check that the counts of number of SOURCEs and COMPONENT types in the
source_catalogue_t matches those provide through argument
*/
void check_single_source_numbers(source_catalogue_t *raw_srccat,
       int num_sources,  int num_shapelets,  int n_comps,
       int n_points, int n_point_lists, int n_point_powers, int n_point_curves,
       int n_gauss, int n_gauss_lists, int n_gauss_powers, int n_gauss_curves,
       int n_shapes, int n_shape_lists, int n_shape_powers, int n_shape_curves,
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
  TEST_ASSERT_EQUAL_INT(n_point_lists,
                               raw_srccat->sources[source_index].n_point_lists);
  TEST_ASSERT_EQUAL_INT(n_point_powers,
                              raw_srccat->sources[source_index].n_point_powers);
  TEST_ASSERT_EQUAL_INT(n_point_curves,
                              raw_srccat->sources[source_index].n_point_curves);

  TEST_ASSERT_EQUAL_INT(n_gauss, raw_srccat->sources[source_index].n_gauss);
  TEST_ASSERT_EQUAL_INT(n_gauss_lists,
                               raw_srccat->sources[source_index].n_gauss_lists);
  TEST_ASSERT_EQUAL_INT(n_gauss_powers,
                              raw_srccat->sources[source_index].n_gauss_powers);
  TEST_ASSERT_EQUAL_INT(n_gauss_curves,
                              raw_srccat->sources[source_index].n_gauss_curves);

  TEST_ASSERT_EQUAL_INT(n_shapes, raw_srccat->sources[source_index].n_shapes);
  TEST_ASSERT_EQUAL_INT(n_shape_lists,
                               raw_srccat->sources[source_index].n_shape_lists);
  TEST_ASSERT_EQUAL_INT(n_shape_powers,
                              raw_srccat->sources[source_index].n_shape_powers);
  TEST_ASSERT_EQUAL_INT(n_shape_curves,
                              raw_srccat->sources[source_index].n_shape_curves);

  TEST_ASSERT_EQUAL_INT(n_shape_coeffs, raw_srccat->sources[source_index].n_shape_coeffs);
}

/*
Checks that the POINT values for the COMPONENT at the given index are correct
Assumes the catalogue being read in is one of the specific sky models made
for testing in the WODEN/cmake_testing/read_and_write dir
*/
void check_single_component_point(source_t source, int comp_index,
                                  e_flux_type flux_type, int flux_ind) {

  TEST_ASSERT_DOUBLE_WITHIN(1e-15, POINT0_RA, source.point_components.ras[comp_index]);
  TEST_ASSERT_DOUBLE_WITHIN(1e-15, POINT0_DEC, source.point_components.decs[comp_index]);

  if (flux_type == POWER_LAW) {
    TEST_ASSERT_DOUBLE_WITHIN(1e-15, POINT0_POW_FREQ, source.point_components.power_ref_freqs[flux_ind]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, POINT0_POW_I, source.point_components.power_ref_stokesI[flux_ind]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, POINT0_POW_Q, source.point_components.power_ref_stokesQ[flux_ind]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, POINT0_POW_U, source.point_components.power_ref_stokesU[flux_ind]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, POINT0_POW_V, source.point_components.power_ref_stokesV[flux_ind]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, POINT0_POW_SI, source.point_components.power_SIs[flux_ind]);
    TEST_ASSERT_EQUAL_INT(comp_index, source.point_components.power_comp_inds[flux_ind]);


  }
  else if (flux_type == CURVED_POWER_LAW) {

    TEST_ASSERT_DOUBLE_WITHIN(1e-15, POINT0_POW_FREQ, source.point_components.curve_ref_freqs[flux_ind]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, POINT0_POW_I, source.point_components.curve_ref_stokesI[flux_ind]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, POINT0_POW_Q, source.point_components.curve_ref_stokesQ[flux_ind]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, POINT0_POW_U, source.point_components.curve_ref_stokesU[flux_ind]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, POINT0_POW_V, source.point_components.curve_ref_stokesV[flux_ind]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, POINT0_POW_SI, source.point_components.curve_SIs[flux_ind]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, POINT0_CURVE_Q, source.point_components.curve_qs[flux_ind]);

    TEST_ASSERT_EQUAL_INT(comp_index, source.point_components.curve_comp_inds[flux_ind]);

  }
  else if (flux_type == LIST) {
    // printf("MADE IT HERE FOR TESTING LIST\n");
    for (int list_ind = 0; list_ind < source.point_components.num_list_values[flux_ind]; list_ind++) {

      // printf("Looky looky %.8f %.8f\n",POINT0_LIST_FREQ[list_ind],
      //               source.point_components.list_freqs[list_ind] );
      TEST_ASSERT_DOUBLE_WITHIN(TOL, POINT0_LIST_FREQ[list_ind],
                       source.point_components.list_freqs[list_ind]);
      TEST_ASSERT_DOUBLE_WITHIN(TOL, POINT0_LIST_I[list_ind],
                       source.point_components.list_stokesI[list_ind]);
      TEST_ASSERT_DOUBLE_WITHIN(TOL, POINT0_LIST_Q[list_ind],
                       source.point_components.list_stokesQ[list_ind]);
      TEST_ASSERT_DOUBLE_WITHIN(TOL, POINT0_LIST_U[list_ind],
                       source.point_components.list_stokesU[list_ind]);
      TEST_ASSERT_DOUBLE_WITHIN(TOL, POINT0_LIST_V[list_ind],
                       source.point_components.list_stokesV[list_ind]);
    }

    TEST_ASSERT_EQUAL_INT(POINT0_NUM_LIST_ENTRIES*flux_ind,
                            source.point_components.list_start_indexes[flux_ind]);

    TEST_ASSERT_EQUAL_INT(POINT0_NUM_LIST_ENTRIES,
                            source.point_components.num_list_values[flux_ind]);

    TEST_ASSERT_EQUAL_INT(comp_index,
                              source.point_components.list_comp_inds[flux_ind]);
  }
}

/*
Checks that the GAUSSIAN values for the COMPONENT at the given index are correct
Assumes the catalogue being read in is one of the specific sky models made
for testing in the WODEN/cmake_testing/read_and_write dir
*/
void check_single_component_gauss(source_t source, int comp_index,
                                  e_flux_type flux_type, int flux_ind) {
  TEST_ASSERT_DOUBLE_WITHIN(1e-15, GAUSS0_RA, source.gauss_components.ras[comp_index]);
  TEST_ASSERT_DOUBLE_WITHIN(1e-15, GAUSS0_DEC, source.gauss_components.decs[comp_index]);
  TEST_ASSERT_DOUBLE_WITHIN(TOL, GAUSS0_PA, source.gauss_components.pas[comp_index]);
  TEST_ASSERT_DOUBLE_WITHIN(TOL, GAUSS0_MAJ,
                          source.gauss_components.majors[comp_index]);
  TEST_ASSERT_DOUBLE_WITHIN(TOL, GAUSS0_MIN,
                          source.gauss_components.minors[comp_index]);

  if (flux_type == POWER_LAW) {
    TEST_ASSERT_DOUBLE_WITHIN(1e-15, GAUSS0_POW_FREQ, source.gauss_components.power_ref_freqs[flux_ind]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, GAUSS0_POW_I, source.gauss_components.power_ref_stokesI[flux_ind]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, GAUSS0_POW_Q, source.gauss_components.power_ref_stokesQ[flux_ind]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, GAUSS0_POW_U, source.gauss_components.power_ref_stokesU[flux_ind]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, GAUSS0_POW_V, source.gauss_components.power_ref_stokesV[flux_ind]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, GAUSS0_POW_SI, source.gauss_components.power_SIs[flux_ind]);
    TEST_ASSERT_EQUAL_INT(comp_index, source.gauss_components.power_comp_inds[flux_ind]);


  }
  else if (flux_type == CURVED_POWER_LAW) {

    TEST_ASSERT_DOUBLE_WITHIN(1e-15, GAUSS0_POW_FREQ, source.gauss_components.curve_ref_freqs[flux_ind]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, GAUSS0_POW_I, source.gauss_components.curve_ref_stokesI[flux_ind]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, GAUSS0_POW_Q, source.gauss_components.curve_ref_stokesQ[flux_ind]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, GAUSS0_POW_U, source.gauss_components.curve_ref_stokesU[flux_ind]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, GAUSS0_POW_V, source.gauss_components.curve_ref_stokesV[flux_ind]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, GAUSS0_POW_SI, source.gauss_components.curve_SIs[flux_ind]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, GAUSS0_CURVE_Q, source.gauss_components.curve_qs[flux_ind]);

    TEST_ASSERT_EQUAL_INT(comp_index, source.gauss_components.curve_comp_inds[flux_ind]);

  }
  else if (flux_type == LIST) {
    // printf("MADE IT HERE FOR TESTING LIST\n");
    for (int list_ind = 0; list_ind < source.gauss_components.num_list_values[flux_ind]; list_ind++) {

      // printf("Looky looky %.8f %.8f\n",GAUSS0_LIST_FREQ[list_ind],
      //               source.gauss_components.list_freqs[list_ind] );
      TEST_ASSERT_DOUBLE_WITHIN(TOL, GAUSS0_LIST_FREQ[list_ind],
                       source.gauss_components.list_freqs[list_ind]);
      TEST_ASSERT_DOUBLE_WITHIN(TOL, GAUSS0_LIST_I[list_ind],
                       source.gauss_components.list_stokesI[list_ind]);
      TEST_ASSERT_DOUBLE_WITHIN(TOL, GAUSS0_LIST_Q[list_ind],
                       source.gauss_components.list_stokesQ[list_ind]);
      TEST_ASSERT_DOUBLE_WITHIN(TOL, GAUSS0_LIST_U[list_ind],
                       source.gauss_components.list_stokesU[list_ind]);
      TEST_ASSERT_DOUBLE_WITHIN(TOL, GAUSS0_LIST_V[list_ind],
                       source.gauss_components.list_stokesV[list_ind]);
    }

    TEST_ASSERT_EQUAL_INT(GAUSS0_NUM_LIST_ENTRIES*flux_ind,
                            source.gauss_components.list_start_indexes[flux_ind]);

    TEST_ASSERT_EQUAL_INT(GAUSS0_NUM_LIST_ENTRIES,
                            source.gauss_components.num_list_values[flux_ind]);

    TEST_ASSERT_EQUAL_INT(comp_index,
                               source.gauss_components.list_comp_inds[flux_ind]);
  }
}

/*
Checks that the SHAPELET values for the COMPONENT at the given index are correct
Assumes the catalogue being read in is one of the specific sky models made
for testing in the WODEN/cmake_testing/read_and_write dir
*/
void check_single_component_shapelet(source_t source, int comp_index,
                                     e_flux_type flux_type, int flux_ind) {
  TEST_ASSERT_DOUBLE_WITHIN(1e-15, SHAPE0_RA, source.shape_components.ras[comp_index]);
  TEST_ASSERT_DOUBLE_WITHIN(1e-15, SHAPE0_DEC, source.shape_components.decs[comp_index]);
  TEST_ASSERT_DOUBLE_WITHIN(TOL, SHAPE0_PA, source.shape_components.pas[comp_index]);
  TEST_ASSERT_DOUBLE_WITHIN(TOL, SHAPE0_MAJ,
                          source.shape_components.majors[comp_index]);
  TEST_ASSERT_DOUBLE_WITHIN(TOL, SHAPE0_MIN,
                          source.shape_components.minors[comp_index]);

  double expected_param_indexes[] = {comp_index, comp_index, comp_index, comp_index};

  for (int param_index = 0; param_index < 4; param_index++) {
    TEST_ASSERT_DOUBLE_WITHIN(TOL, SHAPE0_N1S[param_index],
                      source.shape_components.n1s[4*comp_index + param_index]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, SHAPE0_N2S[param_index],
                      source.shape_components.n2s[4*comp_index + param_index]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, SHAPE0_COEFFS[param_index],
                      source.shape_components.shape_coeffs[4*comp_index + param_index]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, expected_param_indexes[param_index],
                      source.shape_components.param_indexes[4*comp_index + param_index]);
  }

  if (flux_type == POWER_LAW) {
    TEST_ASSERT_DOUBLE_WITHIN(1e-15, SHAPE0_POW_FREQ, source.shape_components.power_ref_freqs[flux_ind]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, SHAPE0_POW_I, source.shape_components.power_ref_stokesI[flux_ind]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, SHAPE0_POW_Q, source.shape_components.power_ref_stokesQ[flux_ind]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, SHAPE0_POW_U, source.shape_components.power_ref_stokesU[flux_ind]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, SHAPE0_POW_V, source.shape_components.power_ref_stokesV[flux_ind]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, SHAPE0_POW_SI, source.shape_components.power_SIs[flux_ind]);
    TEST_ASSERT_EQUAL_INT(comp_index, source.shape_components.power_comp_inds[flux_ind]);


  }
  else if (flux_type == CURVED_POWER_LAW) {

    TEST_ASSERT_DOUBLE_WITHIN(1e-15, SHAPE0_POW_FREQ, source.shape_components.curve_ref_freqs[flux_ind]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, SHAPE0_POW_I, source.shape_components.curve_ref_stokesI[flux_ind]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, SHAPE0_POW_Q, source.shape_components.curve_ref_stokesQ[flux_ind]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, SHAPE0_POW_U, source.shape_components.curve_ref_stokesU[flux_ind]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, SHAPE0_POW_V, source.shape_components.curve_ref_stokesV[flux_ind]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, SHAPE0_POW_SI, source.shape_components.curve_SIs[flux_ind]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, SHAPE0_CURVE_Q, source.shape_components.curve_qs[flux_ind]);

    TEST_ASSERT_EQUAL_INT(comp_index, source.shape_components.curve_comp_inds[flux_ind]);

  }
  else if (flux_type == LIST) {
    // printf("MADE IT HERE FOR TESTING LIST\n");
    for (int list_ind = 0; list_ind < source.shape_components.num_list_values[flux_ind]; list_ind++) {

      // printf("Looky looky %.8f %.8f\n",SHAPE0_LIST_FREQ[list_ind],
      //               source.shape_components.list_freqs[list_ind] );

      // #ifdef DOUBLE_PRECISION
        TEST_ASSERT_DOUBLE_WITHIN(TOL, SHAPE0_LIST_FREQ[list_ind],
                         source.shape_components.list_freqs[list_ind]);
        TEST_ASSERT_DOUBLE_WITHIN(TOL, SHAPE0_LIST_I[list_ind],
                         source.shape_components.list_stokesI[list_ind]);
        TEST_ASSERT_DOUBLE_WITHIN(TOL, SHAPE0_LIST_Q[list_ind],
                         source.shape_components.list_stokesQ[list_ind]);
        TEST_ASSERT_DOUBLE_WITHIN(TOL, SHAPE0_LIST_U[list_ind],
                         source.shape_components.list_stokesU[list_ind]);
        TEST_ASSERT_DOUBLE_WITHIN(TOL, SHAPE0_LIST_V[list_ind],
                         source.shape_components.list_stokesV[list_ind]);
    }

    TEST_ASSERT_EQUAL_INT(SHAPE0_NUM_LIST_ENTRIES*flux_ind,
                            source.shape_components.list_start_indexes[flux_ind]);

    TEST_ASSERT_EQUAL_INT(SHAPE0_NUM_LIST_ENTRIES,
                            source.shape_components.num_list_values[flux_ind]);

    TEST_ASSERT_EQUAL_INT(comp_index,
                              source.shape_components.list_comp_inds[flux_ind]);
  }
}

void set_flux_component_type_nums(e_flux_type flux_type,
      int n_points, int n_gauss, int n_shapes,
      int * n_point_lists, int * n_point_powers,  int * n_point_curves,
      int * n_gauss_lists, int * n_gauss_powers,  int * n_gauss_curves,
      int * n_shape_lists, int * n_shape_powers,  int * n_shape_curves){

  * n_point_lists = 0;
  * n_point_powers = 0;
  * n_point_curves = 0;
  * n_gauss_lists = 0;
  * n_gauss_powers = 0;
  * n_gauss_curves = 0;
  * n_shape_lists = 0;
  * n_shape_powers = 0;
  * n_shape_curves = 0;

  if (flux_type == POWER_LAW) {
    * n_point_powers = n_points;
    * n_gauss_powers = n_gauss;
    * n_shape_powers = n_shapes;

  }
  else if (flux_type == CURVED_POWER_LAW) {
    * n_point_curves = n_points;
    * n_gauss_curves = n_gauss;
    * n_shape_curves = n_shapes;
  }

  else if (flux_type == LIST) {
    * n_point_lists = n_points;
    * n_gauss_lists = n_gauss;
    * n_shape_lists = n_shapes;

  }
}

/*
Test whether a single point source is read in correctly
Use this to test in some allmost broken srclists
*/
void test_read_skymodel_SinglePoint(char *srclist, e_flux_type flux_type)
{


  // Read in the source catalogue
  source_catalogue_t *raw_srccat = malloc( sizeof(source_catalogue_t) );
  int status=0;
  status = read_skymodel(srclist, raw_srccat);

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

  int n_point_lists = 0;
  int n_point_powers = 0;
  int n_point_curves = 0;
  int n_gauss_lists = 0;
  int n_gauss_powers = 0;
  int n_gauss_curves = 0;
  int n_shape_lists = 0;
  int n_shape_powers = 0;
  int n_shape_curves = 0;

  set_flux_component_type_nums(flux_type, n_points, n_gauss, n_shapes,
        &n_point_lists, &n_point_powers,  &n_point_curves,
        &n_gauss_lists, &n_gauss_powers,  &n_gauss_curves,
        &n_shape_lists, &n_shape_powers,  &n_shape_curves);

  //Check these values
  check_single_source_numbers(raw_srccat, num_sources,  num_shapelets,  n_comps,
          n_points, n_point_lists, n_point_powers,  n_point_curves,
          n_gauss, n_gauss_lists, n_gauss_powers,  n_gauss_curves,
          n_shapes, n_shape_lists, n_shape_powers,  n_shape_curves,
          n_shape_components, 0);

  // //Check the POINT specific values in the first COMPONENT
  // //of the first SOURCE are correct
  check_single_component_point(raw_srccat->sources[0], 0, flux_type, 0);

  // free_source_catalogue(raw_srccat);

  // free(raw_srccat);
}

/*
Test whether a single Gaussian source is read in correctly
*/
void test_read_skymodel_SingleGaussian(char *srclist, e_flux_type flux_type)
{
  // Read in the source catalogue
  source_catalogue_t *raw_srccat = malloc( sizeof(source_catalogue_t) );
  int status=0;
  status = read_skymodel(srclist, raw_srccat);

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

  int n_point_lists = 0;
  int n_point_powers = 0;
  int n_point_curves = 0;
  int n_gauss_lists = 0;
  int n_gauss_powers = 0;
  int n_gauss_curves = 0;
  int n_shape_lists = 0;
  int n_shape_powers = 0;
  int n_shape_curves = 0;

  // printf("OVER HERE MATE %d %d %d\n", raw_srccat->sources[0].n_gauss_lists,
  //                               raw_srccat->sources[0].n_gauss_powers,
  //                               raw_srccat->sources[0].n_gauss_curves);

  set_flux_component_type_nums(flux_type, n_points, n_gauss, n_shapes,
        &n_point_lists, &n_point_powers,  &n_point_curves,
        &n_gauss_lists, &n_gauss_powers,  &n_gauss_curves,
        &n_shape_lists, &n_shape_powers,  &n_shape_curves);

  //Check these values
  check_single_source_numbers(raw_srccat, num_sources,  num_shapelets,  n_comps,
          n_points, n_point_lists, n_point_powers,  n_point_curves,
          n_gauss, n_gauss_lists, n_gauss_powers,  n_gauss_curves,
          n_shapes, n_shape_lists, n_shape_powers,  n_shape_curves,
          n_shape_components, 0);

  // Check that the GAUSSIAN specific values in the first COMPONENT
  //of the first SOURCE are correct
  check_single_component_gauss(raw_srccat->sources[0], 0, flux_type, 0);

  // free_source_catalogue(raw_srccat);

}




/*
Test whether a single shapelet source is read in correctly
*/
void test_read_skymodel_SingleShapelet(char *srclist, e_flux_type flux_type)
{
  // Read in the source catalogue
  source_catalogue_t *raw_srccat = malloc( sizeof(source_catalogue_t) );
  int status=0;
  status = read_skymodel(srclist, raw_srccat);

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

  int n_point_lists = 0;
  int n_point_powers = 0;
  int n_point_curves = 0;
  int n_gauss_lists = 0;
  int n_gauss_powers = 0;
  int n_gauss_curves = 0;
  int n_shape_lists = 0;
  int n_shape_powers = 0;
  int n_shape_curves = 0;

  set_flux_component_type_nums(flux_type, n_points, n_gauss, n_shapes,
        &n_point_lists, &n_point_powers,  &n_point_curves,
        &n_gauss_lists, &n_gauss_powers,  &n_gauss_curves,
        &n_shape_lists, &n_shape_powers,  &n_shape_curves);

  //Check these values
  check_single_source_numbers(raw_srccat, num_sources,  num_shapelets,  n_comps,
          n_points, n_point_lists, n_point_powers,  n_point_curves,
          n_gauss, n_gauss_lists, n_gauss_powers,  n_gauss_curves,
          n_shapes, n_shape_lists, n_shape_powers,  n_shape_curves,
          n_shape_components, 0);

  // Check that the SHAPELET specific values in the first COMPONENT
  //of the first SOURCE are correct
  check_single_component_shapelet(raw_srccat->sources[0], 0, flux_type, 0);

  // free_source_catalogue(raw_srccat);

}



/*
Test whether three separate SOURCEs, each with a single COMPONENT,
ordered as POINT, GAUSSIAN, SHAPELET, read in correctly
*/
void test_read_skymodel_ThreeSources(char *srclist, e_flux_type flux_type) {
  // Read in the source catalogue
  source_catalogue_t *raw_srccat = malloc( sizeof(source_catalogue_t) );
  int status=0;
  status = read_skymodel(srclist, raw_srccat);

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

  //These will change depending on which COMPONENT we have iterated to testing
  int n_point_lists = 0;
  int n_point_powers = 0;
  int n_point_curves = 0;
  int n_gauss_lists = 0;
  int n_gauss_powers = 0;
  int n_gauss_curves = 0;
  int n_shape_lists = 0;
  int n_shape_powers = 0;
  int n_shape_curves = 0;

  for (int source_index = 0; source_index < num_sources; source_index++) {

    set_flux_component_type_nums(flux_type, n_points[source_index],
          n_gauss[source_index], n_shapes[source_index],
          &n_point_lists, &n_point_powers,  &n_point_curves,
          &n_gauss_lists, &n_gauss_powers,  &n_gauss_curves,
          &n_shape_lists, &n_shape_powers,  &n_shape_curves);

    //Check these values
    check_single_source_numbers(raw_srccat, num_sources, num_shapelets,
            n_comps[source_index],
            n_points[source_index], n_point_lists, n_point_powers,  n_point_curves,
            n_gauss[source_index], n_gauss_lists, n_gauss_powers,  n_gauss_curves,
            n_shapes[source_index], n_shape_lists, n_shape_powers,  n_shape_curves,
            n_shape_components[source_index], source_index);
  }

  // Check that the POINT specific values in the first COMPONENT
  //of the first SOURCE are correct
  check_single_component_point(raw_srccat->sources[0], 0, flux_type, 0);

  // Check that the GAUSSIAN specific values in the first COMPONENT
  //of the second SOURCE are correct
  check_single_component_gauss(raw_srccat->sources[1], 0, flux_type, 0);

  // Check that the SHAPELET specific values in the first COMPONENT
  //of the third SOURCE are correct
  check_single_component_shapelet(raw_srccat->sources[2], 0, flux_type, 0);

  // free_source_catalogue(raw_srccat);

}



void test_read_skymodel_ThreeComponents(char *srclist, e_flux_type flux_type) {
  // Read in the source catalogue
  source_catalogue_t *raw_srccat = malloc( sizeof(source_catalogue_t) );
  int status=0;
  status = read_skymodel(srclist, raw_srccat);

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

  int n_point_lists = 0;
  int n_point_powers = 0;
  int n_point_curves = 0;
  int n_gauss_lists = 0;
  int n_gauss_powers = 0;
  int n_gauss_curves = 0;
  int n_shape_lists = 0;
  int n_shape_powers = 0;
  int n_shape_curves = 0;

  set_flux_component_type_nums(flux_type, n_points, n_gauss, n_shapes,
        &n_point_lists, &n_point_powers,  &n_point_curves,
        &n_gauss_lists, &n_gauss_powers,  &n_gauss_curves,
        &n_shape_lists, &n_shape_powers,  &n_shape_curves);

  int source_index = 0;

  check_single_source_numbers(raw_srccat, num_sources,  num_shapelets,  n_comps,
          n_points, n_point_lists, n_point_powers,  n_point_curves,
          n_gauss, n_gauss_lists, n_gauss_powers,  n_gauss_curves,
          n_shapes, n_shape_lists, n_shape_powers,  n_shape_curves,
                              n_shape_components, source_index);

  // Check that the POINT specific values in the first COMPONENT
  //of the first SOURCE are correct
  check_single_component_point(raw_srccat->sources[0], 0, flux_type, 0);

  // Check that the GAUSSIAN specific values in the first COMPONENT
  //of the first SOURCE are correct
  check_single_component_gauss(raw_srccat->sources[0], 0, flux_type, 0);

  // Check that the SHAPELET specific values in the first COMPONENT
  //of the first SOURCE are correct
  check_single_component_shapelet(raw_srccat->sources[0], 0, flux_type, 0);

  // free_source_catalogue(raw_srccat);

}
