#include <unity.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "chunk_sky_model.h"
#include "woden_struct_defs.h"
#include "fill_chunks_common.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

#define UNITY_INCLUDE_FLOAT

/*
The GPU will never have enough memory to simulate millions of COMPONENTs at
once, so we need to split the sky model up before running the simulation.
`fill_chunk_src_*` functions do just that, given the sky model above the
horizon (`cropped_src`) and a chunking index, they return a portion of the
sky model containing specific COMPONENT types.
*/

void test_create_chunked_sky_models(int chunking_size,
                                    int num_points, int num_gauss,
                                    int num_shapes, int num_coeff_per_shape,
                                    int num_list_values,
                                    int num_time_steps,
                                    int num_baselines, int num_freqs) {

  //Make the sky model given the input params
  source_t *cropped_src = make_sky_model(num_points, num_gauss, num_shapes,
                                        num_coeff_per_shape,
                                        num_list_values,
                                        num_time_steps);

  printf("====================================================\n");
  printf("cropped_src->n_comps %d\n", cropped_src->n_comps);
  printf("cropped_src->n_points %d\n", cropped_src->n_points);
  printf("cropped_src->n_point_lists %d\n", cropped_src->n_point_lists);
  printf("cropped_src->n_point_powers %d\n", cropped_src->n_point_powers);
  printf("cropped_src->n_point_curves %d\n", cropped_src->n_point_curves);
  printf("cropped_src->n_gauss %d\n", cropped_src->n_gauss);
  printf("cropped_src->n_gauss_lists %d\n", cropped_src->n_gauss_lists);
  printf("cropped_src->n_gauss_powers %d\n", cropped_src->n_gauss_powers);
  printf("cropped_src->n_gauss_curves %d\n", cropped_src->n_gauss_curves);
  printf("cropped_src->n_shapes %d\n", cropped_src->n_shapes);
  printf("cropped_src->n_shape_lists %d\n", cropped_src->n_shape_lists);
  printf("cropped_src->n_shape_powers %d\n", cropped_src->n_shape_powers);
  printf("cropped_src->n_shape_curves %d\n", cropped_src->n_shape_curves);
  printf("cropped_src->n_shape_coeffs %d\n", cropped_src->n_shape_coeffs);
  printf("====================================================\n");


  //Set some dummy baseline, frequency simulation settings
  // int num_visis = num_baselines * num_freqs * num_time_steps;

  woden_settings_t *woden_settings = malloc(sizeof(woden_settings_t));
  woden_settings->num_freqs = num_freqs;
  woden_settings->num_time_steps = num_time_steps;
  woden_settings->chunking_size = chunking_size;
  woden_settings->num_baselines = num_baselines;
  //
  // // //
  // int num_comps_to_chunk = cropped_src->n_points + cropped_src->n_gauss;
  // int comps_per_chunk = (int)floorf((user_precision_t)chunking_size / (user_precision_t)(num_baselines * num_freqs * num_time_steps));
  // int num_comp_chunks = (int)ceilf((user_precision_t)num_comps_to_chunk / (user_precision_t)comps_per_chunk);
  //
  // int num_coeffs_to_chunk = cropped_src->n_shape_coeffs;
  // int coeffs_per_chunk = (int)floorf((user_precision_t)chunking_size / (user_precision_t)(num_baselines * num_freqs * num_time_steps));
  // int num_coeff_chunks = (int)ceilf((user_precision_t)num_coeffs_to_chunk / (user_precision_t)coeffs_per_chunk);
  // //


  int comps_per_chunk = (int)floorf((float)chunking_size / (float)(num_baselines * num_freqs * num_time_steps));

  if (comps_per_chunk < 1) {
    comps_per_chunk = 1;
  }
  //Numbers of things we are splitting up
  //How many chunks we need for the POINT/GAUSS, and the SHAPELETs
  int num_point_chunks = (int)ceilf((float)(NUM_FLUX_TYPES*num_points) / (float)comps_per_chunk);
  int num_gauss_chunks = (int)ceilf((float)(NUM_FLUX_TYPES*num_gauss) / (float)comps_per_chunk);
  int num_shape_chunks = (int)ceilf((float)cropped_src->n_shape_coeffs / (float)comps_per_chunk);

  //
  // //Function to be testsd
  source_catalogue_t *cropped_sky_models = create_chunked_sky_models(cropped_src, woden_settings);
  // //
  //Counters to see how far many point/gauss have already been assigned
  //to previous chunks

  int point_comp_accum = 0;
  int point_power_accum = 0;
  int point_curve_accum = 0;
  int point_list_accum = 0;


  for (int chunk_ind = 0; chunk_ind < num_point_chunks; chunk_ind++) {

    //Check outputs are as expected
    check_pointgauss_chunking(chunk_ind, comps_per_chunk, num_point_chunks,
                              num_list_values,
                              &point_comp_accum, &point_power_accum,
                              &point_curve_accum, &point_list_accum,
                              POINT,
                              cropped_src,
                              &cropped_sky_models->sources[chunk_ind]);
  }

  int gauss_comp_accum = 0;
  int gauss_power_accum = 0;
  int gauss_curve_accum = 0;
  int gauss_list_accum = 0;


  for (int chunk_ind = 0; chunk_ind <  num_gauss_chunks; chunk_ind++) {

    //Check outputs are as expected
    check_pointgauss_chunking(chunk_ind, comps_per_chunk, num_gauss_chunks,
                              num_list_values,
                              &gauss_comp_accum, &gauss_power_accum,
                              &gauss_curve_accum, &gauss_list_accum,
                              GAUSSIAN,
                              cropped_src,
                              &cropped_sky_models->sources[num_point_chunks + chunk_ind]);
  }

  int num_pg_to_chunk = num_point_chunks + num_gauss_chunks;

  for (int chunk_ind = 0; chunk_ind < num_shape_chunks; chunk_ind++) {

    check_shapelet_chunking(chunk_ind, comps_per_chunk, num_time_steps,
                            num_coeff_per_shape, cropped_src,
                            &cropped_sky_models->sources[num_pg_to_chunk + chunk_ind]);

  } //END iteration over all chunks
  // free_sky_model(cropped_sky_models);


  // TEST_ASSERT_EQUAL_INT(cropped_src->n_comps, );
  TEST_ASSERT_EQUAL_INT(cropped_src->n_points, point_comp_accum);
  TEST_ASSERT_EQUAL_INT(cropped_src->n_point_lists, point_power_accum);
  TEST_ASSERT_EQUAL_INT(cropped_src->n_point_powers, point_curve_accum);
  TEST_ASSERT_EQUAL_INT(cropped_src->n_point_curves, point_list_accum);

  TEST_ASSERT_EQUAL_INT(cropped_src->n_gauss, gauss_comp_accum);
  TEST_ASSERT_EQUAL_INT(cropped_src->n_gauss_lists, gauss_power_accum);
  TEST_ASSERT_EQUAL_INT(cropped_src->n_gauss_powers, gauss_curve_accum);
  TEST_ASSERT_EQUAL_INT(cropped_src->n_gauss_curves, gauss_list_accum);

}


void test_create_chunked_sky_models_P3_G0000_S00_00_C1_Time001(void) {
  int num_points = 3;
  int num_gauss = 0;
  int num_shapes = 0;
  int num_coeff_per_shape = 0;
  int chunking_size = 2;

  int num_time_steps = 1;
  int num_baselines = 1;
  int num_freqs = 1;

  int num_list_values = 3;
  test_create_chunked_sky_models(chunking_size, num_points, num_gauss,
                                 num_shapes,  num_coeff_per_shape,
                                 num_list_values,
                                 num_time_steps, num_baselines, num_freqs);
}



void test_create_chunked_sky_models_P4544_G1736_S20_14_C1e8_Time014(void){

  int num_points = 4544;
  int num_gauss = 1736;
  int num_shapes = 20;
  int num_coeff_per_shape = 14;

  int chunking_size = 1e8;

  int num_time_steps = 14;
  int num_baselines = 8128;
  int num_freqs = 16;

  int num_list_values = 4;
  test_create_chunked_sky_models(chunking_size, num_points, num_gauss,
                                 num_shapes,  num_coeff_per_shape,
                                 num_list_values,
                                 num_time_steps, num_baselines, num_freqs);
}

void test_create_chunked_sky_models_P1000_G0000_S00_00_C1e8_Time014(void){

  int num_points = 1000;
  int num_gauss = 0;
  int num_shapes = 0;
  int num_coeff_per_shape = 0;

  int chunking_size = 1e9;

  int num_time_steps = 14;
  int num_baselines = 8128;
  int num_freqs = 16;

  int num_list_values = 6;
  test_create_chunked_sky_models(chunking_size, num_points, num_gauss,
                                 num_shapes,  num_coeff_per_shape,
                                 num_list_values,
                                 num_time_steps, num_baselines, num_freqs);
}

void test_create_chunked_sky_models_P0000_G1000_S00_00_C1e8_Time014(void){

  int num_points = 0;
  int num_gauss = 1000;
  int num_shapes = 0;
  int num_coeff_per_shape = 0;

  int chunking_size = 1e9;

  int num_time_steps = 14;
  int num_baselines = 8128;
  int num_freqs = 16;

  int num_list_values = 5;
  test_create_chunked_sky_models(chunking_size, num_points, num_gauss,
                                 num_shapes,  num_coeff_per_shape,
                                 num_list_values,
                                 num_time_steps, num_baselines, num_freqs);
}

void test_create_chunked_sky_models_P0000_G0000_S67_78_C1e8_Time014(void){

  int num_points = 0;
  int num_gauss = 0;
  int num_shapes = 67;
  int num_coeff_per_shape = 78;

  int chunking_size = 1e9;

  int num_time_steps = 14;
  int num_baselines = 8128;
  int num_freqs = 16;

  int num_list_values = 2;
  test_create_chunked_sky_models(chunking_size, num_points, num_gauss,
                                 num_shapes,  num_coeff_per_shape,
                                 num_list_values,
                                 num_time_steps, num_baselines, num_freqs);
}

void test_create_chunked_sky_models_P0050_G0050_S10_25_C4576_Time014(void){

  int num_points = 50;
  int num_gauss = 50;
  int num_shapes = 10;
  int num_coeff_per_shape = 25;

  int chunking_size = 4576;

  int num_time_steps = 14;
  int num_baselines = 8128;
  int num_freqs = 16;

  int num_list_values = 8;
  test_create_chunked_sky_models(chunking_size, num_points, num_gauss,
                                 num_shapes,  num_coeff_per_shape,
                                 num_list_values,
                                 num_time_steps, num_baselines, num_freqs);
}

void test_create_chunked_sky_models_P87760_G12207_S121_678_C1e10_Time056(void){

  int num_points = 87760;
  int num_gauss = 12207;
  int num_shapes = 121;
  int num_coeff_per_shape = 678;

  long int chunking_size = 1e10;

  int num_time_steps = 56;
  int num_baselines = 8128;
  int num_freqs = 32;

  int num_list_values = 16;
  test_create_chunked_sky_models(chunking_size, num_points, num_gauss,
                                 num_shapes,  num_coeff_per_shape,
                                 num_list_values,
                                 num_time_steps, num_baselines, num_freqs);
}

//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_create_chunked_sky_models_P3_G0000_S00_00_C1_Time001);
    RUN_TEST(test_create_chunked_sky_models_P1000_G0000_S00_00_C1e8_Time014);
    RUN_TEST(test_create_chunked_sky_models_P0000_G1000_S00_00_C1e8_Time014);
    RUN_TEST(test_create_chunked_sky_models_P0000_G0000_S67_78_C1e8_Time014);
    RUN_TEST(test_create_chunked_sky_models_P0050_G0050_S10_25_C4576_Time014);
    RUN_TEST(test_create_chunked_sky_models_P4544_G1736_S20_14_C1e8_Time014);
    RUN_TEST(test_create_chunked_sky_models_P87760_G12207_S121_678_C1e10_Time056);

    return UNITY_END();
}
