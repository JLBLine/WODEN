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

void test_fill_chunk_src_with_shapelets(int chunking_size,
                                         int num_shapes, int num_coeff_per_shape,
                                         int num_time_steps) {

  //Not testing POINT/GAUSS, set these to zero
  int num_points = 0;
  int num_gauss = 0;
  catsource_t *cropped_src = make_sky_model(num_points, num_gauss, num_shapes,
                                            num_coeff_per_shape, num_time_steps);

  //Set some dummy baseline, frequency simulation settings
  int num_baselines = 5;
  int num_freqs = 2;
  int num_visis = num_baselines * num_freqs * num_time_steps;

  int num_coeffs_to_chunk = cropped_src->n_shape_coeffs;
  int coeffs_per_chunk = (int)floorf((float)chunking_size / (float)(num_baselines * num_freqs * num_time_steps));
  int num_chunks = (int)ceilf((float)num_coeffs_to_chunk / (float)coeffs_per_chunk);

  printf("Num visis %d\n",num_visis );
  printf("Chunking size %d\n",chunking_size );
  printf("coeffs per chunk %d\n", coeffs_per_chunk );
  printf("Num chunks %d\n", num_chunks );
  // printf("Is num_visis*num_coeffs %d < num_visis*coeffs_per_chunk*num_chunks %d\n",
         // num_visis*num_coeffs_to_chunk, num_visis*coeffs_per_chunk*num_chunks);

  catsource_t *temp_cropped_src = malloc(sizeof(catsource_t));
  woden_settings_t *woden_settings = malloc(sizeof(woden_settings_t));
  woden_settings->num_freqs = num_freqs;
  woden_settings->num_time_steps = num_time_steps;

  for (int chunk_ind = 0; chunk_ind < num_chunks; chunk_ind++) {

    fill_chunk_src_with_shapelets(temp_cropped_src, cropped_src, chunk_ind,
                                   coeffs_per_chunk, woden_settings);

    // Check the outputs are correct for this chunk

    check_shapelet_chunking(chunk_ind, coeffs_per_chunk, num_time_steps,
                            num_coeff_per_shape, cropped_src,
                            temp_cropped_src);

  } //END iteration over all chunks
  free_sky_model(cropped_src);
}


void test_fill_chunk_src_with_shapelets_Shapelet010_Coeff10_Chunk1000_Time004(void){
  int chunking_size = 1000;
  int num_shapes = 10;
  int num_coeff_per_shape = 10;
  int num_time_steps = 4;

  test_fill_chunk_src_with_shapelets(chunking_size, num_shapes, num_coeff_per_shape,
                                      num_time_steps);
}

void test_fill_chunk_src_with_shapelets_Shapelet010_Coeff13_Chunk1000_Time004(void){
  int chunking_size = 1000;
  int num_shapes = 10;
  int num_coeff_per_shape = 13;
  int num_time_steps = 4;

  test_fill_chunk_src_with_shapelets(chunking_size, num_shapes, num_coeff_per_shape,
                                      num_time_steps);
}

void test_fill_chunk_src_with_shapelets_Shapelet193_Coeff1266_Chunk3913_Time007(void){
  int chunking_size = 3913;
  int num_shapes = 193;
  int num_coeff_per_shape = 1266;
  int num_time_steps = 7;

  test_fill_chunk_src_with_shapelets(chunking_size, num_shapes, num_coeff_per_shape,
                                      num_time_steps);
}

//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_fill_chunk_src_with_shapelets_Shapelet010_Coeff10_Chunk1000_Time004);
    RUN_TEST(test_fill_chunk_src_with_shapelets_Shapelet010_Coeff13_Chunk1000_Time004);
    RUN_TEST(test_fill_chunk_src_with_shapelets_Shapelet193_Coeff1266_Chunk3913_Time007);

    return UNITY_END();
}
