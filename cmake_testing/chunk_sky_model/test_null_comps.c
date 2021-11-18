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
the sky model. Test by setting up a polpulated `catsource_t` and seeing if
correct parts become NULL
*/

user_precision_t one_array[] = {1};
double one_array_double[] = {1};

/*
Make the polpulated catsource_t struct. Just stick everything equal to 1.0
*/
catsource_t * make_sky_model(void) {

  catsource_t *src = malloc(sizeof(catsource_t));

  src->n_comps = 1;
  src->n_points= 1;
  src->n_gauss = 1;
  src->n_shapes = 1;
  src->n_shape_coeffs = 1;

  src->point_ras = one_array_double;
  src->point_decs = one_array_double;
  src->point_ref_freqs = one_array_double;
  src->point_ref_stokesI = one_array;
  src->point_ref_stokesQ = one_array;
  src->point_ref_stokesU = one_array;
  src->point_ref_stokesV = one_array;
  src->point_SIs = one_array;
  src->point_azs = one_array;
  src->point_zas = one_array;
  src->cos_point_para_angs = one_array;
  src->sin_point_para_angs = one_array;
  src->point_gaussbeam_has = one_array_double;
  src->point_gaussbeam_decs = one_array_double;

  src->gauss_ras = one_array_double;
  src->gauss_decs = one_array_double;
  src->gauss_ref_freqs = one_array_double;
  src->gauss_ref_stokesI = one_array;
  src->gauss_ref_stokesQ = one_array;
  src->gauss_ref_stokesU = one_array;
  src->gauss_ref_stokesV = one_array;
  src->gauss_SIs = one_array;
  src->gauss_majors = one_array;
  src->gauss_minors = one_array;
  src->gauss_pas = one_array;
  src->gauss_azs = one_array;
  src->gauss_zas = one_array;
  src->cos_gauss_para_angs = one_array;
  src->sin_gauss_para_angs = one_array;
  src->gauss_gaussbeam_has = one_array_double;
  src->gauss_gaussbeam_decs = one_array_double;

  src->shape_ras = one_array_double;
  src->shape_decs = one_array_double;
  src->shape_ref_freqs = one_array_double;
  src->shape_ref_stokesI = one_array;
  src->shape_ref_stokesQ = one_array;
  src->shape_ref_stokesU = one_array;
  src->shape_ref_stokesV = one_array;
  src->shape_SIs = one_array;
  src->shape_majors = one_array;
  src->shape_minors = one_array;
  src->shape_pas = one_array;
  src->shape_azs = one_array;
  src->shape_zas = one_array;
  src->cos_shape_para_angs = one_array;
  src->sin_shape_para_angs = one_array;
  src->shape_gaussbeam_has = one_array_double;
  src->shape_gaussbeam_decs = one_array_double;
  src->shape_coeffs = one_array;
  src->shape_n1s = one_array;
  src->shape_n2s = one_array;
  src->shape_param_indexes = one_array;

  return src;

}


/*
Check that POINT values have NOT be NULL-ed
*/
void assert_point_retained(catsource_t *src) {

  TEST_ASSERT_EQUAL_INT(1, src->n_points);

  TEST_ASSERT_EQUAL_DOUBLE_ARRAY( one_array_double, src->point_ras, 1);
  TEST_ASSERT_EQUAL_DOUBLE_ARRAY( one_array_double, src->point_decs, 1);
  TEST_ASSERT_EQUAL_DOUBLE_ARRAY( one_array_double, src->point_ref_freqs, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->point_ref_stokesI, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->point_ref_stokesQ, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->point_ref_stokesU, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->point_ref_stokesV, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->point_SIs, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->point_azs, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->point_zas, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->cos_point_para_angs, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->sin_point_para_angs, 1);
  TEST_ASSERT_EQUAL_DOUBLE_ARRAY( one_array_double, src->point_gaussbeam_has, 1);
  TEST_ASSERT_EQUAL_DOUBLE_ARRAY( one_array_double, src->point_gaussbeam_decs, 1);
}


/*
Check that GAUSSIAN values have NOT be NULL-ed
*/
void assert_gauss_retained(catsource_t *src) {

  TEST_ASSERT_EQUAL_INT(1, src->n_gauss);

  TEST_ASSERT_EQUAL_DOUBLE_ARRAY( one_array_double, src->gauss_ras, 1);
  TEST_ASSERT_EQUAL_DOUBLE_ARRAY( one_array_double, src->gauss_decs, 1);
  TEST_ASSERT_EQUAL_DOUBLE_ARRAY( one_array_double, src->gauss_ref_freqs, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->gauss_ref_stokesI, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->gauss_ref_stokesQ, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->gauss_ref_stokesU, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->gauss_ref_stokesV, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->gauss_SIs, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->gauss_majors, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->gauss_minors, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->gauss_pas, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->gauss_azs, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->gauss_zas, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->cos_gauss_para_angs, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->sin_gauss_para_angs, 1);
  TEST_ASSERT_EQUAL_DOUBLE_ARRAY( one_array_double, src->gauss_gaussbeam_has, 1);
  TEST_ASSERT_EQUAL_DOUBLE_ARRAY( one_array_double, src->gauss_gaussbeam_decs, 1);
}

/*
Check that SHAPELET values have NOT be NULL-ed
*/
void assert_shape_retained(catsource_t *src) {

  TEST_ASSERT_EQUAL_INT(1, src->n_shapes);
  TEST_ASSERT_EQUAL_INT(1, src->n_shape_coeffs);

  TEST_ASSERT_EQUAL_DOUBLE_ARRAY( one_array_double, src->shape_ras, 1);
  TEST_ASSERT_EQUAL_DOUBLE_ARRAY( one_array_double, src->shape_decs, 1);
  TEST_ASSERT_EQUAL_DOUBLE_ARRAY( one_array_double, src->shape_ref_freqs, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->shape_ref_stokesI, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->shape_ref_stokesQ, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->shape_ref_stokesU, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->shape_ref_stokesV, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->shape_SIs, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->shape_majors, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->shape_minors, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->shape_pas, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->shape_azs, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->shape_zas, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->cos_shape_para_angs, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->sin_shape_para_angs, 1);
  TEST_ASSERT_EQUAL_DOUBLE_ARRAY( one_array_double, src->shape_gaussbeam_has, 1);
  TEST_ASSERT_EQUAL_DOUBLE_ARRAY( one_array_double, src->shape_gaussbeam_decs, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->shape_coeffs, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->shape_n1s, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->shape_n2s, 1);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY( one_array, src->shape_param_indexes, 1);
}

/*
Test POINT source nulling works correctly
*/
void test_null_point_comps_DoesTheNull(void) {
  //Make sky model
  catsource_t *src =  make_sky_model();

  //Function to be tested
  null_point_comps(src);

  //Test things are nulled
  TEST_ASSERT_EQUAL_INT(0, src->n_points);

  TEST_ASSERT_NULL(src->point_ras);
  TEST_ASSERT_NULL(src->point_decs);
  TEST_ASSERT_NULL(src->point_ref_freqs);
  TEST_ASSERT_NULL(src->point_ref_stokesI);
  TEST_ASSERT_NULL(src->point_ref_stokesQ);
  TEST_ASSERT_NULL(src->point_ref_stokesU);
  TEST_ASSERT_NULL(src->point_ref_stokesV);
  TEST_ASSERT_NULL(src->point_SIs);
  TEST_ASSERT_NULL(src->point_azs);
  TEST_ASSERT_NULL(src->point_zas);
  TEST_ASSERT_NULL(src->cos_point_para_angs);
  TEST_ASSERT_NULL(src->sin_point_para_angs);
  TEST_ASSERT_NULL(src->point_gaussbeam_has);
  TEST_ASSERT_NULL(src->point_gaussbeam_decs);

  //Check the GAUSS and SHAPELET stuff is left alone
  assert_gauss_retained(src);
  assert_shape_retained(src);

}

/*
Test GAUSS source nulling works correctly
*/
void test_null_gauss_comps_DoesTheNull(void) {
  //Make sky model
  catsource_t *src =  make_sky_model();

  //Function to be tested
  null_gauss_comps(src);

  //Test things are nulled
  TEST_ASSERT_EQUAL_INT(0, src->n_gauss);

  TEST_ASSERT_NULL(src->gauss_ras);
  TEST_ASSERT_NULL(src->gauss_decs);
  TEST_ASSERT_NULL(src->gauss_ref_freqs);
  TEST_ASSERT_NULL(src->gauss_ref_stokesI);
  TEST_ASSERT_NULL(src->gauss_ref_stokesQ);
  TEST_ASSERT_NULL(src->gauss_ref_stokesU);
  TEST_ASSERT_NULL(src->gauss_ref_stokesV);
  TEST_ASSERT_NULL(src->gauss_SIs);
  TEST_ASSERT_NULL(src->gauss_majors);
  TEST_ASSERT_NULL(src->gauss_minors);
  TEST_ASSERT_NULL(src->gauss_pas);
  TEST_ASSERT_NULL(src->gauss_azs);
  TEST_ASSERT_NULL(src->gauss_zas);
  TEST_ASSERT_NULL(src->cos_gauss_para_angs);
  TEST_ASSERT_NULL(src->sin_gauss_para_angs);
  TEST_ASSERT_NULL(src->gauss_gaussbeam_has);
  TEST_ASSERT_NULL(src->gauss_gaussbeam_decs);

  //Check the POINT and SHAPELET stuff is left alone
  assert_point_retained(src);
  assert_shape_retained(src);

}


/*
Test SHAPELET source nulling works correctly
*/
void test_null_shape_comps_DoesTheNull(void) {
  //Make sky model
  catsource_t *src =  make_sky_model();

  //Function to be tested
  null_shapelet_comps(src);

  //Test things are nulled
  TEST_ASSERT_EQUAL_INT(0, src->n_shapes);
  TEST_ASSERT_EQUAL_INT(0, src->n_shape_coeffs);

  TEST_ASSERT_NULL(src->shape_ras);
  TEST_ASSERT_NULL(src->shape_decs);
  TEST_ASSERT_NULL(src->shape_ref_freqs);
  TEST_ASSERT_NULL(src->shape_ref_stokesI);
  TEST_ASSERT_NULL(src->shape_ref_stokesQ);
  TEST_ASSERT_NULL(src->shape_ref_stokesU);
  TEST_ASSERT_NULL(src->shape_ref_stokesV);
  TEST_ASSERT_NULL(src->shape_SIs);
  TEST_ASSERT_NULL(src->shape_majors);
  TEST_ASSERT_NULL(src->shape_minors);
  TEST_ASSERT_NULL(src->shape_pas);
  TEST_ASSERT_NULL(src->shape_azs);
  TEST_ASSERT_NULL(src->shape_zas);
  TEST_ASSERT_NULL(src->cos_shape_para_angs);
  TEST_ASSERT_NULL(src->sin_shape_para_angs);
  TEST_ASSERT_NULL(src->shape_gaussbeam_has);
  TEST_ASSERT_NULL(src->shape_gaussbeam_decs);
  TEST_ASSERT_NULL(src->shape_coeffs);
  TEST_ASSERT_NULL(src->shape_n1s);
  TEST_ASSERT_NULL(src->shape_n2s);
  TEST_ASSERT_NULL(src->shape_param_indexes);

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
