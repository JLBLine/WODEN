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
`chunk_sky_model::null_*_comps` are functions that NULLs out specific parts of
the sky model. Test by setting up a polpulated `source_t` and seeing if
correct parts become NULL
*/

user_precision_t one_array[] = {1};
double one_array_double[] = {1};
int one_array_int[] = {1};

/*
Set things inside a components_t to be one or arrays of one
*/
void set_component_arrays_to_one(components_t * components, e_component_type comptype){
  components->ras = one_array_double;
  components->decs = one_array_double;

  //Power law stuff
  components->power_ref_freqs = one_array_double;
  components->power_ref_stokesI = one_array;
  components->power_ref_stokesQ = one_array;
  components->power_ref_stokesU = one_array;
  components->power_ref_stokesV = one_array;
  components->power_SIs = one_array;
  components->power_comp_inds = one_array_int;

  //Curved power law stuff
  components->curve_ref_stokesI = one_array;
  components->curve_ref_stokesQ = one_array;
  components->curve_ref_stokesU = one_array;
  components->curve_ref_stokesV = one_array;
  components->curve_ref_freqs = one_array_double;
  components->curve_SIs = one_array;
  components->curve_qs = one_array;
  components->curve_comp_inds = one_array_int;
  //List flux things
  components->list_freqs = one_array_double;
  components->list_stokesI = one_array;
  components->list_stokesQ = one_array;
  components->list_stokesU = one_array;
  components->list_stokesV = one_array;
  components->num_list_values = one_array_int;
  components->list_start_indexes = one_array_int;
  components->list_comp_inds = one_array_int;

  //WODEN things
  components->azs = one_array;
  components->zas = one_array;
  components->beam_has = one_array_double;
  components->beam_decs = one_array_double;


  if (comptype == GAUSSIAN || comptype == SHAPELET) {
    components->majors = one_array;
    components->minors = one_array;
    components->pas = one_array;
  }

  if (comptype == SHAPELET) {
    components->shape_coeffs = one_array;
    components->n1s = one_array;
    components->n2s = one_array;
    components->param_indexes = one_array;

  }
}

/*
Make the polpulated source_t struct. Just stick everything equal to 1.0
*/
source_t * make_sky_model(void) {

  source_t *src = malloc(sizeof(source_t));

  src->n_comps = 1;
  src->n_points = 1;
  src->n_point_lists = 1;
  src->n_point_powers = 1;
  src->n_point_curves = 1;
  src->n_gauss = 1;
  src->n_gauss_lists = 1;
  src->n_gauss_powers = 1;
  src->n_gauss_curves = 1;
  src->n_shapes = 1;
  src->n_shape_lists = 1;
  src->n_shape_powers = 1;
  src->n_shape_curves = 1;
  src->n_shape_coeffs = 1;

  set_component_arrays_to_one(&src->point_components, POINT);
  set_component_arrays_to_one(&src->gauss_components, GAUSSIAN);
  set_component_arrays_to_one(&src->shape_components, SHAPELET);

  return src;

}

/*
Set things inside a components_t to be one or arrays of one
*/
void check_component_arrays_are_one(components_t components, e_component_type comptype){
  TEST_ASSERT_EQUAL_DOUBLE_ARRAY(one_array_double, components.ras, 1);
  TEST_ASSERT_EQUAL_DOUBLE_ARRAY(one_array_double, components.decs, 1);

  //Power law stuff
  TEST_ASSERT_EQUAL_DOUBLE_ARRAY(one_array_double, components.power_ref_freqs, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(one_array, components.power_ref_stokesI, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(one_array, components.power_ref_stokesQ, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(one_array, components.power_ref_stokesU, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(one_array, components.power_ref_stokesV, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(one_array, components.power_SIs, 1);
  TEST_ASSERT_EQUAL_INT_ARRAY(one_array_int, components.power_comp_inds, 1);

  //Curved power law stuff
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(one_array, components.curve_ref_stokesI, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(one_array, components.curve_ref_stokesQ, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(one_array, components.curve_ref_stokesU, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(one_array, components.curve_ref_stokesV, 1);
  TEST_ASSERT_EQUAL_DOUBLE_ARRAY(one_array_double, components.curve_ref_freqs, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(one_array, components.curve_SIs, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(one_array, components.curve_qs, 1);
  TEST_ASSERT_EQUAL_INT_ARRAY(one_array_int, components.curve_comp_inds, 1);
  //List flux things
  TEST_ASSERT_EQUAL_DOUBLE_ARRAY(one_array_double, components.list_freqs, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(one_array, components.list_stokesI, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(one_array, components.list_stokesQ, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(one_array, components.list_stokesU, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(one_array, components.list_stokesV, 1);
  TEST_ASSERT_EQUAL_INT_ARRAY(one_array_int, components.num_list_values, 1);
  TEST_ASSERT_EQUAL_INT_ARRAY(one_array_int, components.list_start_indexes, 1);
  TEST_ASSERT_EQUAL_INT_ARRAY(one_array_int, components.list_comp_inds, 1);

  // //WODEN things
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(one_array, components.azs, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(one_array, components.zas, 1);
  TEST_ASSERT_EQUAL_DOUBLE_ARRAY(one_array_double, components.beam_has, 1);
  TEST_ASSERT_EQUAL_DOUBLE_ARRAY(one_array_double, components.beam_decs, 1);
  

  if (comptype == GAUSSIAN || comptype == SHAPELET) {
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(one_array, components.majors, 1);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(one_array, components.minors, 1);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(one_array, components.pas, 1);
  }

  if (comptype == SHAPELET) {
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(one_array, components.shape_coeffs, 1);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(one_array, components.n1s, 1);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(one_array, components.n2s, 1);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(one_array, components.param_indexes, 1);

  }
}




void assert_made_null(components_t components, e_component_type comptype){
  TEST_ASSERT_NULL(components.ras);
  TEST_ASSERT_NULL(components.decs);

  //Power law stuff
  TEST_ASSERT_NULL(components.power_ref_freqs);
  TEST_ASSERT_NULL(components.power_ref_stokesI);
  TEST_ASSERT_NULL(components.power_ref_stokesQ);
  TEST_ASSERT_NULL(components.power_ref_stokesU);
  TEST_ASSERT_NULL(components.power_ref_stokesV);
  TEST_ASSERT_NULL(components.power_SIs);
  TEST_ASSERT_NULL(components.power_comp_inds);

  //Curved power law stuff
  TEST_ASSERT_NULL(components.curve_ref_stokesI);
  TEST_ASSERT_NULL(components.curve_ref_stokesQ);
  TEST_ASSERT_NULL(components.curve_ref_stokesU);
  TEST_ASSERT_NULL(components.curve_ref_stokesV);
  TEST_ASSERT_NULL(components.curve_ref_freqs);
  TEST_ASSERT_NULL(components.curve_SIs);
  TEST_ASSERT_NULL(components.curve_qs);
  TEST_ASSERT_NULL(components.curve_comp_inds);
  //List flux things
  TEST_ASSERT_NULL(components.list_freqs);
  TEST_ASSERT_NULL(components.list_stokesI);
  TEST_ASSERT_NULL(components.list_stokesQ);
  TEST_ASSERT_NULL(components.list_stokesU);
  TEST_ASSERT_NULL(components.list_stokesV);
  TEST_ASSERT_NULL(components.num_list_values);
  TEST_ASSERT_NULL(components.list_start_indexes);
  TEST_ASSERT_NULL(components.list_comp_inds);

  // //WODEN things
  TEST_ASSERT_NULL(components.azs);
  TEST_ASSERT_NULL(components.zas);
  TEST_ASSERT_NULL(components.beam_has);
  TEST_ASSERT_NULL(components.beam_decs);


  if (comptype == GAUSSIAN || comptype == SHAPELET) {
    TEST_ASSERT_NULL(components.majors);
    TEST_ASSERT_NULL(components.minors);
    TEST_ASSERT_NULL(components.pas);
  }

  if (comptype == SHAPELET) {
    TEST_ASSERT_NULL(components.shape_coeffs);
    TEST_ASSERT_NULL(components.n1s);
    TEST_ASSERT_NULL(components.n2s);
    TEST_ASSERT_NULL(components.param_indexes);

  }
}

/*
Check that POINT values have NOT be NULL-ed
*/
void assert_point_retained(source_t *src) {

  TEST_ASSERT_EQUAL_INT(1, src->n_points);
  TEST_ASSERT_EQUAL_INT(1, src->n_point_lists);
  TEST_ASSERT_EQUAL_INT(1, src->n_point_powers);
  TEST_ASSERT_EQUAL_INT(1, src->n_point_curves);

  check_component_arrays_are_one(src->point_components, POINT);
}

/*
Check that GAUSSIAN values have NOT be NULL-ed
*/
void assert_gauss_retained(source_t *src) {

  TEST_ASSERT_EQUAL_INT(1, src->n_gauss);
  TEST_ASSERT_EQUAL_INT(1, src->n_gauss_lists);
  TEST_ASSERT_EQUAL_INT(1, src->n_gauss_powers);
  TEST_ASSERT_EQUAL_INT(1, src->n_gauss_curves);

  check_component_arrays_are_one(src->gauss_components, GAUSSIAN);
}

/*
Check that SHAPELET values have NOT be NULL-ed
*/
void assert_shape_retained(source_t *src) {

  TEST_ASSERT_EQUAL_INT(1, src->n_shapes);
  TEST_ASSERT_EQUAL_INT(1, src->n_shape_lists);
  TEST_ASSERT_EQUAL_INT(1, src->n_shape_powers);
  TEST_ASSERT_EQUAL_INT(1, src->n_shape_curves);

  check_component_arrays_are_one(src->shape_components, SHAPELET);
}

/*
Test POINT source nulling works correctly
*/
void test_null_point_comps_DoesTheNull(void) {
  //Make sky model
  source_t *src =  make_sky_model();

  //Function to be tested
  null_components(src, POINT);

  //Test things are nulled
  TEST_ASSERT_EQUAL_INT(0, src->n_points);
  TEST_ASSERT_EQUAL_INT(0, src->n_point_lists);
  TEST_ASSERT_EQUAL_INT(0, src->n_point_powers);
  TEST_ASSERT_EQUAL_INT(0, src->n_point_curves);

  assert_made_null(src->point_components, POINT);

  // Check the GAUSS and SHAPELET stuff is left alone
  assert_gauss_retained(src);
  assert_shape_retained(src);

}

/*
Test GAUSS source nulling works correctly
*/
void test_null_gauss_comps_DoesTheNull(void) {
  //Make sky model
  source_t *src =  make_sky_model();

  //Function to be tested
  null_components(src, GAUSSIAN);

  //Test things are nulled
  TEST_ASSERT_EQUAL_INT(0, src->n_gauss);
  TEST_ASSERT_EQUAL_INT(0, src->n_gauss_lists);
  TEST_ASSERT_EQUAL_INT(0, src->n_gauss_powers);
  TEST_ASSERT_EQUAL_INT(0, src->n_gauss_curves);
  assert_made_null(src->gauss_components, GAUSSIAN);

  //Check the POINT and SHAPELET stuff is left alone
  assert_point_retained(src);
  assert_shape_retained(src);

}


/*
Test SHAPELET source nulling works correctly
*/
void test_null_shape_comps_DoesTheNull(void) {
  //Make sky model
  source_t *src =  make_sky_model();

  //Function to be tested
  null_components(src, SHAPELET);

  //Test things are nulled
  TEST_ASSERT_EQUAL_INT(0, src->n_shapes);
  TEST_ASSERT_EQUAL_INT(0, src->n_shape_coeffs);
  TEST_ASSERT_EQUAL_INT(0, src->n_shape_lists);
  TEST_ASSERT_EQUAL_INT(0, src->n_shape_powers);
  TEST_ASSERT_EQUAL_INT(0, src->n_shape_curves);
  assert_made_null(src->shape_components, SHAPELET);

  //Check the POINT and GAUSSIAN stuff is left alone
  assert_point_retained(src);
  assert_gauss_retained(src);

}

//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_null_point_comps_DoesTheNull);
    RUN_TEST(test_null_gauss_comps_DoesTheNull);
    RUN_TEST(test_null_shape_comps_DoesTheNull);

    return UNITY_END();
}
