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

void test_fill_chunk_src_with_pointgauss(int chunking_size,
                                         int num_points, int num_gauss,
                                         int num_list_values,
                                         int num_time_steps,
                                         int num_baselines, int num_freqs) {

  //Not testing SHAPELETs, set these to zero
  int num_shapes = 0;
  int num_coeff_per_shape = 0;
  source_t *cropped_src = make_sky_model(num_points, num_gauss, num_shapes,
                                        num_coeff_per_shape, num_list_values,
                                        num_time_steps);

  woden_settings_t *woden_settings = malloc(sizeof(woden_settings_t));
  woden_settings->num_freqs = num_freqs;
  woden_settings->num_time_steps = num_time_steps;

  int comps_per_chunk = (int)floorf((float)chunking_size / (float)(num_baselines * num_freqs * num_time_steps));

  int num_point_chunks = (int)ceilf((float)(NUM_FLUX_TYPES*num_points) / (float)comps_per_chunk);
  int num_gauss_chunks = (int)ceilf((float)(NUM_FLUX_TYPES*num_gauss) / (float)comps_per_chunk);

  source_t *temp_cropped_src = malloc(sizeof(source_t));
  //
  //Counters to see how far many point/gauss have already been assigned
  //to previous chunks
  int point_comp_accum = 0;
  int point_power_accum = 0;
  int point_curve_accum = 0;
  int point_list_accum = 0;

  for (int chunk_ind = 0; chunk_ind < num_point_chunks; chunk_ind++) {

    fill_chunk_src_with_pointgauss(temp_cropped_src, cropped_src, chunk_ind,
                                   comps_per_chunk, woden_settings, POINT);

    check_pointgauss_chunking(chunk_ind, comps_per_chunk, num_point_chunks,
                              num_list_values,
                              &point_comp_accum, &point_power_accum,
                              &point_curve_accum, &point_list_accum,
                              POINT,
                              cropped_src,
                              temp_cropped_src);
  }

  int gauss_comp_accum = 0;
  int gauss_power_accum = 0;
  int gauss_curve_accum = 0;
  int gauss_list_accum = 0;
  for (int chunk_ind = 0; chunk_ind < num_gauss_chunks; chunk_ind++) {

    fill_chunk_src_with_pointgauss(temp_cropped_src, cropped_src, chunk_ind,
                                   comps_per_chunk, woden_settings, GAUSSIAN);

    check_pointgauss_chunking(chunk_ind, comps_per_chunk, num_point_chunks,
                              num_list_values,
                              &gauss_comp_accum, &gauss_power_accum,
                              &gauss_curve_accum, &gauss_list_accum,
                              GAUSSIAN,
                              cropped_src,
                              temp_cropped_src);
  }

  // free_sky_model(cropped_src);
}


void test_fill_chunk_src_with_pointgauss_Point100_Gauss000_Chunk1000_Time004(void){
  int chunking_size = 1000;
  int num_points = 100;
  int num_gauss = 0;
  int num_time_steps = 4;
  int num_baselines = 5;
  int num_freqs = 2;

  int num_list_values = 2;
  test_fill_chunk_src_with_pointgauss(chunking_size, num_points, num_gauss,
                                      num_list_values,
                                      num_time_steps, num_baselines, num_freqs);
}

void test_fill_chunk_src_with_pointgauss_Point000_Gauss100_Chunk1000_Time004(void){
  int chunking_size = 1000;
  int num_points = 0;
  int num_gauss = 100;
  int num_time_steps = 4;
  int num_baselines = 5;
  int num_freqs = 2;

  int num_list_values = 3;
  test_fill_chunk_src_with_pointgauss(chunking_size, num_points, num_gauss,
                                      num_list_values,
                                      num_time_steps, num_baselines, num_freqs);
}

void test_fill_chunk_src_with_pointgauss_Point100_Gauss100_Chunk1000_Time004(void){
  int chunking_size = 1000;
  int num_points = 100;
  int num_gauss = 100;
  int num_time_steps = 4;
  int num_baselines = 5;
  int num_freqs = 2;

  int num_list_values = 4;
  test_fill_chunk_src_with_pointgauss(chunking_size, num_points, num_gauss,
                                      num_list_values,
                                      num_time_steps, num_baselines, num_freqs);
}

void test_fill_chunk_src_with_pointgauss_Point100_Gauss000_Chunk1000_Time003(void){
  int chunking_size = 1000;
  int num_points = 100;
  int num_gauss = 0;
  int num_time_steps = 3;
  int num_baselines = 5;
  int num_freqs = 2;

  int num_list_values = 5;
  test_fill_chunk_src_with_pointgauss(chunking_size, num_points, num_gauss,
                                      num_list_values,
                                      num_time_steps, num_baselines, num_freqs);
}

void test_fill_chunk_src_with_pointgauss_Point000_Gauss100_Chunk1000_Time003(void){
  int chunking_size = 1000;
  int num_points = 0;
  int num_gauss = 100;
  int num_time_steps = 3;
  int num_baselines = 5;
  int num_freqs = 2;

  int num_list_values = 6;
  test_fill_chunk_src_with_pointgauss(chunking_size, num_points, num_gauss,
                                      num_list_values,
                                      num_time_steps, num_baselines, num_freqs);
}

void test_fill_chunk_src_with_pointgauss_Point100_Gauss100_Chunk1000_Time003(void){
  int chunking_size = 1000;
  int num_points = 100;
  int num_gauss = 100;
  int num_time_steps = 3;
  int num_baselines = 5;
  int num_freqs = 2;

  int num_list_values = 5;
  test_fill_chunk_src_with_pointgauss(chunking_size, num_points, num_gauss,
                                      num_list_values,
                                      num_time_steps, num_baselines, num_freqs);
}

void test_fill_chunk_src_with_pointgauss_Point100_Gauss000_Chunk0173_Time005(void){
  int chunking_size = 173;
  int num_points = 100;
  int num_gauss = 0;
  int num_time_steps = 5;
  int num_baselines = 5;
  int num_freqs = 2;

  int num_list_values = 4;
  test_fill_chunk_src_with_pointgauss(chunking_size, num_points, num_gauss,
                                      num_list_values,
                                      num_time_steps, num_baselines, num_freqs);
}

void test_fill_chunk_src_with_pointgauss_Point000_Gauss100_Chunk0173_Time005(void){
  int chunking_size = 173;
  int num_points = 0;
  int num_gauss = 100;
  int num_time_steps = 5;
  int num_baselines = 5;
  int num_freqs = 2;

  int num_list_values = 3;
  test_fill_chunk_src_with_pointgauss(chunking_size, num_points, num_gauss,
                                      num_list_values,
                                      num_time_steps, num_baselines, num_freqs);
}

void test_fill_chunk_src_with_pointgauss_Point100_Gauss100_Chunk0173_Time005(void){
  int chunking_size = 173;
  int num_points = 100;
  int num_gauss = 100;
  int num_time_steps = 5;
  int num_baselines = 5;
  int num_freqs = 2;

  int num_list_values = 2;
  test_fill_chunk_src_with_pointgauss(chunking_size, num_points, num_gauss,
                                      num_list_values,
                                      num_time_steps, num_baselines, num_freqs);
}

void test_fill_chunk_src_with_pointgauss_Point10743_Gauss00000_Chunk4345_Time012(void){
  int chunking_size = 4345;
  int num_points = 10743;
  int num_gauss = 0;
  int num_time_steps = 12;
  int num_baselines = 5;
  int num_freqs = 2;

  int num_list_values = 16;
  test_fill_chunk_src_with_pointgauss(chunking_size, num_points, num_gauss,
                                      num_list_values,
                                      num_time_steps, num_baselines, num_freqs);
}

void test_fill_chunk_src_with_pointgauss_Point00000_Gauss16789_Chunk4345_Time012(void){
  int chunking_size = 4345;
  int num_points = 0;
  int num_gauss = 16789;
  int num_time_steps = 12;
  int num_baselines = 5;
  int num_freqs = 2;

  int num_list_values = 20;
  test_fill_chunk_src_with_pointgauss(chunking_size, num_points, num_gauss,
                                      num_list_values,
                                      num_time_steps, num_baselines, num_freqs);
}

void test_fill_chunk_src_with_pointgauss_Point10743_Gauss16789_Chunk4345_Time012(void){
  int chunking_size = 4345;
  int num_points = 10743;
  int num_gauss = 16789;
  int num_time_steps = 12;
  int num_baselines = 5;
  int num_freqs = 2;

  int num_list_values = 10;
  test_fill_chunk_src_with_pointgauss(chunking_size, num_points, num_gauss,
                                      num_list_values,
                                      num_time_steps, num_baselines, num_freqs);
}

//
void test_fill_chunk_src_with_pointgauss_Point4544_Gauss1736_Chunk1e8_Time014(void){
  int chunking_size = 1e8;
  int num_points = 4544;
  int num_gauss = 1736;
  int num_time_steps = 14;
  int num_baselines = 8128;
  int num_freqs = 16;

  int num_list_values = 12;
  test_fill_chunk_src_with_pointgauss(chunking_size, num_points, num_gauss,
                                      num_list_values,
                                      num_time_steps, num_baselines, num_freqs);
}

//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_fill_chunk_src_with_pointgauss_Point100_Gauss000_Chunk1000_Time004);
    RUN_TEST(test_fill_chunk_src_with_pointgauss_Point000_Gauss100_Chunk1000_Time004);
    RUN_TEST(test_fill_chunk_src_with_pointgauss_Point100_Gauss100_Chunk1000_Time004);
    RUN_TEST(test_fill_chunk_src_with_pointgauss_Point100_Gauss000_Chunk1000_Time003);
    RUN_TEST(test_fill_chunk_src_with_pointgauss_Point000_Gauss100_Chunk1000_Time003);
    RUN_TEST(test_fill_chunk_src_with_pointgauss_Point100_Gauss100_Chunk1000_Time003);
    RUN_TEST(test_fill_chunk_src_with_pointgauss_Point100_Gauss000_Chunk0173_Time005);
    RUN_TEST(test_fill_chunk_src_with_pointgauss_Point000_Gauss100_Chunk0173_Time005);
    RUN_TEST(test_fill_chunk_src_with_pointgauss_Point100_Gauss100_Chunk0173_Time005);
    RUN_TEST(test_fill_chunk_src_with_pointgauss_Point10743_Gauss00000_Chunk4345_Time012);
    RUN_TEST(test_fill_chunk_src_with_pointgauss_Point00000_Gauss16789_Chunk4345_Time012);
    RUN_TEST(test_fill_chunk_src_with_pointgauss_Point4544_Gauss1736_Chunk1e8_Time014);

    return UNITY_END();
}
