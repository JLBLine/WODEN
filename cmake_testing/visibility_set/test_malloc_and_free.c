#include <unity.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "woden_struct_defs.h"
#include "visibility_set.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

#define UNITY_INCLUDE_FLOAT

void allocate_visi_set(visibility_set_t *visibility_set) {
  visibility_set->us_metres = malloc(sizeof(float));
  visibility_set->vs_metres = malloc(sizeof(float));
  visibility_set->ws_metres = malloc(sizeof(float));
  visibility_set->sum_visi_XX_real = malloc(sizeof(float));
  visibility_set->sum_visi_XX_imag = malloc(sizeof(float));
  visibility_set->sum_visi_XY_real = malloc(sizeof(float));
  visibility_set->sum_visi_XY_imag = malloc(sizeof(float));
  visibility_set->sum_visi_YX_real = malloc(sizeof(float));
  visibility_set->sum_visi_YX_imag = malloc(sizeof(float));
  visibility_set->sum_visi_YY_real = malloc(sizeof(float));
  visibility_set->sum_visi_YY_imag = malloc(sizeof(float));
  visibility_set->allsteps_sha0s = malloc(sizeof(float));
  visibility_set->allsteps_cha0s = malloc(sizeof(float));
  visibility_set->allsteps_lsts = malloc(sizeof(float));
  visibility_set->allsteps_wavelengths = malloc(sizeof(float));
  visibility_set->channel_frequencies = malloc(sizeof(float));
}

/*
Check some malloc and free functions
*/
void test_setup_visibility_set() {
  //Output container
  visibility_set_t *visibility_set = setup_visibility_set(1);

  TEST_ASSERT_NOT_NULL(visibility_set->us_metres);
  TEST_ASSERT_NOT_NULL(visibility_set->vs_metres);
  TEST_ASSERT_NOT_NULL(visibility_set->ws_metres);
  TEST_ASSERT_NOT_NULL(visibility_set->sum_visi_XX_real);
  TEST_ASSERT_NOT_NULL(visibility_set->sum_visi_XX_imag);
  TEST_ASSERT_NOT_NULL(visibility_set->sum_visi_XY_real);
  TEST_ASSERT_NOT_NULL(visibility_set->sum_visi_XY_imag);
  TEST_ASSERT_NOT_NULL(visibility_set->sum_visi_YX_real);
  TEST_ASSERT_NOT_NULL(visibility_set->sum_visi_YX_imag);
  TEST_ASSERT_NOT_NULL(visibility_set->sum_visi_YY_real);
  TEST_ASSERT_NOT_NULL(visibility_set->sum_visi_YY_imag);

  free(visibility_set->us_metres);
  free(visibility_set->vs_metres);
  free(visibility_set->ws_metres);
  free(visibility_set->sum_visi_XX_real);
  free(visibility_set->sum_visi_XX_imag);
  free(visibility_set->sum_visi_XY_real);
  free(visibility_set->sum_visi_XY_imag);
  free(visibility_set->sum_visi_YX_real);
  free(visibility_set->sum_visi_YX_imag);
  free(visibility_set->sum_visi_YY_real);
  free(visibility_set->sum_visi_YY_imag);

}

void test_free_visi_set_inputs() {

  //Allocate some memory in the visibility_set
  visibility_set_t *visibility_set = malloc(sizeof(visibility_set_t));
  allocate_visi_set(visibility_set);

  free_visi_set_inputs(visibility_set);

  //Should be able to free these guys still
  free(visibility_set->us_metres);
  free(visibility_set->vs_metres);
  free(visibility_set->ws_metres);
  free(visibility_set->sum_visi_XX_real);
  free(visibility_set->sum_visi_XX_imag);
  free(visibility_set->sum_visi_XY_real);
  free(visibility_set->sum_visi_XY_imag);
  free(visibility_set->sum_visi_YX_real);
  free(visibility_set->sum_visi_YX_imag);
  free(visibility_set->sum_visi_YY_real);
  free(visibility_set->sum_visi_YY_imag);

}
void test_free_visi_set_outputs() {

  //Allocate some memory in the visibility_set
  visibility_set_t *visibility_set = malloc(sizeof(visibility_set_t));
  allocate_visi_set(visibility_set);

  free_visi_set_outputs(visibility_set);

  //Should be able to free these guys still
  free(visibility_set->allsteps_sha0s);
  free(visibility_set->allsteps_cha0s);
  free(visibility_set->allsteps_lsts);
  free(visibility_set->allsteps_wavelengths);
  free(visibility_set->channel_frequencies);

}

//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_setup_visibility_set);
    RUN_TEST(test_free_visi_set_inputs);
    RUN_TEST(test_free_visi_set_outputs);

    return UNITY_END();
}
