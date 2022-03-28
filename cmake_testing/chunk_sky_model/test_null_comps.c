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

/*
Make the polpulated source_t struct. Just stick everything equal to 1.0
*/
source_t * make_sky_model(void) {

  source_t *src = malloc(sizeof(source_t));

  src->n_comps = 1;
  src->n_points= 1;
  src->n_gauss = 1;
  src->n_shapes = 1;
  src->n_shape_coeffs = 1;

  src->point_components.ras = one_array_double;
  src->point_components.decs = one_array_double;
  src->point_components.ref_freqs = one_array_double;
  src->point_components.ref_stokesI = one_array;
  src->point_components.ref_stokesQ = one_array;
  src->point_components.ref_stokesU = one_array;
  src->point_components.ref_stokesV = one_array;
  src->point_components.SIs = one_array;
  src->point_components.azs = one_array;
  src->point_components.zas = one_array;
  src->point_components.beam_has = one_array_double;
  src->point_components.beam_decs = one_array_double;

  src->gauss_components.ras = one_array_double;
  src->gauss_components.decs = one_array_double;
  src->gauss_components.ref_freqs = one_array_double;
  src->gauss_components.ref_stokesI = one_array;
  src->gauss_components.ref_stokesQ = one_array;
  src->gauss_components.ref_stokesU = one_array;
  src->gauss_components.ref_stokesV = one_array;
  src->gauss_components.SIs = one_array;
  src->gauss_components.majors = one_array;
  src->gauss_components.minors = one_array;
  src->gauss_components.pas = one_array;
  src->gauss_components.azs = one_array;
  src->gauss_components.zas = one_array;
  src->gauss_components.beam_has = one_array_double;
  src->gauss_components.beam_decs = one_array_double;

  src->shape_components.ras = one_array_double;
  src->shape_components.decs = one_array_double;
  src->shape_components.ref_freqs = one_array_double;
  src->shape_components.ref_stokesI = one_array;
  src->shape_components.ref_stokesQ = one_array;
  src->shape_components.ref_stokesU = one_array;
  src->shape_components.ref_stokesV = one_array;
  src->shape_components.SIs = one_array;
  src->shape_components.majors = one_array;
  src->shape_components.minors = one_array;
  src->shape_components.pas = one_array;
  src->shape_components.azs = one_array;
  src->shape_components.zas = one_array;
  src->shape_components.beam_has = one_array_double;
  src->shape_components.beam_decs = one_array_double;
  src->shape_components.shape_coeffs = one_array;
  src->shape_components.n1s = one_array;
  src->shape_components.n2s = one_array;
  src->shape_components.param_indexes = one_array;

  return src;

}


/*
Check that POINT values have NOT be NULL-ed
*/
void assert_point_retained(source_t *src) {

  TEST_ASSERT_EQUAL_INT(1, src->n_points);

  TEST_ASSERT_EQUAL_DOUBLE_ARRAY( one_array_double, src->point_components.ras, 1);
  TEST_ASSERT_EQUAL_DOUBLE_ARRAY( one_array_double, src->point_components.decs, 1);
  TEST_ASSERT_EQUAL_DOUBLE_ARRAY( one_array_double, src->point_components.ref_freqs, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->point_components.ref_stokesI, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->point_components.ref_stokesQ, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->point_components.ref_stokesU, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->point_components.ref_stokesV, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->point_components.SIs, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->point_components.azs, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->point_components.zas, 1);
  TEST_ASSERT_EQUAL_DOUBLE_ARRAY( one_array_double, src->point_components.beam_has, 1);
  TEST_ASSERT_EQUAL_DOUBLE_ARRAY( one_array_double, src->point_components.beam_decs, 1);
}


/*
Check that GAUSSIAN values have NOT be NULL-ed
*/
void assert_gauss_retained(source_t *src) {

  TEST_ASSERT_EQUAL_INT(1, src->n_gauss);

  TEST_ASSERT_EQUAL_DOUBLE_ARRAY( one_array_double, src->gauss_components.ras, 1);
  TEST_ASSERT_EQUAL_DOUBLE_ARRAY( one_array_double, src->gauss_components.decs, 1);
  TEST_ASSERT_EQUAL_DOUBLE_ARRAY( one_array_double, src->gauss_components.ref_freqs, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->gauss_components.ref_stokesI, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->gauss_components.ref_stokesQ, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->gauss_components.ref_stokesU, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->gauss_components.ref_stokesV, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->gauss_components.SIs, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->gauss_components.majors, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->gauss_components.minors, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->gauss_components.pas, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->gauss_components.azs, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->gauss_components.zas, 1);
  TEST_ASSERT_EQUAL_DOUBLE_ARRAY( one_array_double, src->gauss_components.beam_has, 1);
  TEST_ASSERT_EQUAL_DOUBLE_ARRAY( one_array_double, src->gauss_components.beam_decs, 1);
}

/*
Check that SHAPELET values have NOT be NULL-ed
*/
void assert_shape_retained(source_t *src) {

  TEST_ASSERT_EQUAL_INT(1, src->n_shapes);
  TEST_ASSERT_EQUAL_INT(1, src->n_shape_coeffs);

  TEST_ASSERT_EQUAL_DOUBLE_ARRAY( one_array_double, src->shape_components.ras, 1);
  TEST_ASSERT_EQUAL_DOUBLE_ARRAY( one_array_double, src->shape_components.decs, 1);
  TEST_ASSERT_EQUAL_DOUBLE_ARRAY( one_array_double, src->shape_components.ref_freqs, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->shape_components.ref_stokesI, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->shape_components.ref_stokesQ, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->shape_components.ref_stokesU, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->shape_components.ref_stokesV, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->shape_components.SIs, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->shape_components.majors, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->shape_components.minors, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->shape_components.pas, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->shape_components.azs, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->shape_components.zas, 1);
  TEST_ASSERT_EQUAL_DOUBLE_ARRAY( one_array_double, src->shape_components.beam_has, 1);
  TEST_ASSERT_EQUAL_DOUBLE_ARRAY( one_array_double, src->shape_components.beam_decs, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->shape_components.shape_coeffs, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->shape_components.n1s, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->shape_components.n2s, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->shape_components.param_indexes, 1);
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

  TEST_ASSERT_NULL(src->point_components.ras);
  TEST_ASSERT_NULL(src->point_components.decs);
  TEST_ASSERT_NULL(src->point_components.ref_freqs);
  TEST_ASSERT_NULL(src->point_components.ref_stokesI);
  TEST_ASSERT_NULL(src->point_components.ref_stokesQ);
  TEST_ASSERT_NULL(src->point_components.ref_stokesU);
  TEST_ASSERT_NULL(src->point_components.ref_stokesV);
  TEST_ASSERT_NULL(src->point_components.SIs);
  TEST_ASSERT_NULL(src->point_components.azs);
  TEST_ASSERT_NULL(src->point_components.zas);
  TEST_ASSERT_NULL(src->point_components.beam_has);
  TEST_ASSERT_NULL(src->point_components.beam_decs);

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

  TEST_ASSERT_NULL(src->gauss_components.ras);
  TEST_ASSERT_NULL(src->gauss_components.decs);
  TEST_ASSERT_NULL(src->gauss_components.ref_freqs);
  TEST_ASSERT_NULL(src->gauss_components.ref_stokesI);
  TEST_ASSERT_NULL(src->gauss_components.ref_stokesQ);
  TEST_ASSERT_NULL(src->gauss_components.ref_stokesU);
  TEST_ASSERT_NULL(src->gauss_components.ref_stokesV);
  TEST_ASSERT_NULL(src->gauss_components.SIs);
  TEST_ASSERT_NULL(src->gauss_components.majors);
  TEST_ASSERT_NULL(src->gauss_components.minors);
  TEST_ASSERT_NULL(src->gauss_components.pas);
  TEST_ASSERT_NULL(src->gauss_components.azs);
  TEST_ASSERT_NULL(src->gauss_components.zas);
  TEST_ASSERT_NULL(src->gauss_components.beam_has);
  TEST_ASSERT_NULL(src->gauss_components.beam_decs);

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

  TEST_ASSERT_NULL(src->shape_components.ras);
  TEST_ASSERT_NULL(src->shape_components.decs);
  TEST_ASSERT_NULL(src->shape_components.ref_freqs);
  TEST_ASSERT_NULL(src->shape_components.ref_stokesI);
  TEST_ASSERT_NULL(src->shape_components.ref_stokesQ);
  TEST_ASSERT_NULL(src->shape_components.ref_stokesU);
  TEST_ASSERT_NULL(src->shape_components.ref_stokesV);
  TEST_ASSERT_NULL(src->shape_components.SIs);
  TEST_ASSERT_NULL(src->shape_components.majors);
  TEST_ASSERT_NULL(src->shape_components.minors);
  TEST_ASSERT_NULL(src->shape_components.pas);
  TEST_ASSERT_NULL(src->shape_components.azs);
  TEST_ASSERT_NULL(src->shape_components.zas);
  TEST_ASSERT_NULL(src->shape_components.beam_has);
  TEST_ASSERT_NULL(src->shape_components.beam_decs);
  TEST_ASSERT_NULL(src->shape_components.shape_coeffs);
  TEST_ASSERT_NULL(src->shape_components.n1s);
  TEST_ASSERT_NULL(src->shape_components.n2s);
  TEST_ASSERT_NULL(src->shape_components.param_indexes);

  //Check the POINT and SHAPELET stuff is left alone
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
