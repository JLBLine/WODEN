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
                                    int num_time_steps,
                                    int num_baselines, int num_freqs) {

  //Make the sky model given the input params
  catsource_t *cropped_src = make_sky_model(num_points, num_gauss, num_shapes,
                                            num_coeff_per_shape, num_time_steps);

  //Set some dummy baseline, frequency simulation settings
  int num_visis = num_baselines * num_freqs * num_time_steps;

  woden_settings_t *woden_settings = malloc(sizeof(woden_settings_t));
  woden_settings->num_freqs = num_freqs;
  woden_settings->num_time_steps = num_time_steps;
  woden_settings->chunking_size = chunking_size;
  woden_settings->num_baselines = num_baselines;

  // //
  int num_comps_to_chunk = cropped_src->n_points + cropped_src->n_gauss;
  int comps_per_chunk = (int)floorf((float)chunking_size / (float)(num_baselines * num_freqs * num_time_steps));
  int num_comp_chunks = (int)ceilf((float)num_comps_to_chunk / (float)comps_per_chunk);

  int num_coeffs_to_chunk = cropped_src->n_shape_coeffs;
  int coeffs_per_chunk = (int)floorf((float)chunking_size / (float)(num_baselines * num_freqs * num_time_steps));
  int num_coeff_chunks = (int)ceilf((float)num_coeffs_to_chunk / (float)coeffs_per_chunk);
  //

  //Function to be testsd
  source_catalogue_t *cropped_sky_models = create_chunked_sky_models(cropped_src, woden_settings);
  //
  //Counters to see how far many point/gauss have already been assigned
  //to previous chunks
  int point_accum = 0;
  int gauss_accum = 0;

  for (int chunk_ind = 0; chunk_ind < num_comp_chunks; chunk_ind++) {

    //Check outputs are as expected
    check_pointgauss_chunking(chunk_ind, comps_per_chunk,
                              num_time_steps,
                              &point_accum, &gauss_accum,
                              cropped_src,
                              &cropped_sky_models->catsources[chunk_ind]);
  }

  for (int chunk_ind = num_comps_to_chunk; chunk_ind < num_comp_chunks + num_coeff_chunks; chunk_ind++) {

    check_shapelet_chunking(chunk_ind, coeffs_per_chunk, num_time_steps,
                            num_coeff_per_shape, cropped_src,
                            &cropped_sky_models->catsources[chunk_ind]);

  } //END iteration over all chunks
  // free_sky_model(cropped_sky_models);
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

  test_create_chunked_sky_models(chunking_size, num_points, num_gauss,
                                 num_shapes,  num_coeff_per_shape,
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

  test_create_chunked_sky_models(chunking_size, num_points, num_gauss,
                                 num_shapes,  num_coeff_per_shape,
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

  test_create_chunked_sky_models(chunking_size, num_points, num_gauss,
                                 num_shapes,  num_coeff_per_shape,
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

  test_create_chunked_sky_models(chunking_size, num_points, num_gauss,
                                 num_shapes,  num_coeff_per_shape,
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

  test_create_chunked_sky_models(chunking_size, num_points, num_gauss,
                                 num_shapes,  num_coeff_per_shape,
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

  test_create_chunked_sky_models(chunking_size, num_points, num_gauss,
                                 num_shapes,  num_coeff_per_shape,
                                 num_time_steps, num_baselines, num_freqs);
}

//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_create_chunked_sky_models_P1000_G0000_S00_00_C1e8_Time014);
    RUN_TEST(test_create_chunked_sky_models_P0000_G1000_S00_00_C1e8_Time014);
    RUN_TEST(test_create_chunked_sky_models_P0000_G0000_S67_78_C1e8_Time014);
    RUN_TEST(test_create_chunked_sky_models_P0050_G0050_S10_25_C4576_Time014);
    RUN_TEST(test_create_chunked_sky_models_P4544_G1736_S20_14_C1e8_Time014);
    RUN_TEST(test_create_chunked_sky_models_P87760_G12207_S121_678_C1e10_Time056);

    return UNITY_END();
}
