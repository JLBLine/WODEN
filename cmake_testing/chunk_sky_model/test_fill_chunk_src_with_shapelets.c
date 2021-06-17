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

    //How far through all coeffs are we with this chunk
    int chunk_coeff_ind = chunk_ind*coeffs_per_chunk;
    //How many POINT sources are left after this chunk index
    int coeff_remainder = cropped_src->n_shape_coeffs - chunk_coeff_ind;

    //Things to use in logic below
    int expected_n_coeffs = 0;
    int chunk_remainder = 0;

    if (coeff_remainder > 0) { //There are SHAPELET sources, how many should there be?
      if (coeff_remainder >= coeffs_per_chunk) { //SHAPELET sources fill the whole chunk
        expected_n_coeffs = coeffs_per_chunk;
      } else { //Not enough shapelet coeffs to fill the chunk
        expected_n_coeffs = coeff_remainder;
      }
    }
    else {
      printf("SHOULD NOT BE HERE\n");
    }

    if (expected_n_coeffs > 0) {
      // printf("Found POINT, chunk_ind %d expected_n_points %d\n",chunk_ind, expected_n_points );

      TEST_ASSERT_EQUAL_INT(expected_n_coeffs, temp_cropped_src->n_shape_coeffs);
      TEST_ASSERT_EQUAL_INT(cropped_src->n_shapes, temp_cropped_src->n_shapes);
      //
      //As we only split over basis function coeff info, all of these arrrays
      //should just be pointer copies
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(cropped_src->shape_ras,
                              temp_cropped_src->shape_ras, cropped_src->n_shapes);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(cropped_src->shape_decs,
                              temp_cropped_src->shape_decs, cropped_src->n_shapes);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(cropped_src->shape_ref_freqs,
                              temp_cropped_src->shape_ref_freqs, cropped_src->n_shapes);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(cropped_src->shape_ref_stokesI,
                              temp_cropped_src->shape_ref_stokesI, cropped_src->n_shapes);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(cropped_src->shape_ref_stokesQ,
                              temp_cropped_src->shape_ref_stokesQ, cropped_src->n_shapes);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(cropped_src->shape_ref_stokesU,
                              temp_cropped_src->shape_ref_stokesU, cropped_src->n_shapes);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(cropped_src->shape_ref_stokesV,
                              temp_cropped_src->shape_ref_stokesV, cropped_src->n_shapes);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(cropped_src->shape_SIs,
                              temp_cropped_src->shape_SIs, cropped_src->n_shapes);
      //
      //As we only split over basis function coeff info, all of these arrrays
      //should just be pointer copies
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(cropped_src->shape_azs,
                      temp_cropped_src->shape_azs, cropped_src->n_shapes*num_time_steps);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(cropped_src->shape_zas,
                      temp_cropped_src->shape_zas, cropped_src->n_shapes*num_time_steps);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(cropped_src->cos_shape_para_angs,
            temp_cropped_src->cos_shape_para_angs, cropped_src->n_shapes*num_time_steps);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(cropped_src->sin_shape_para_angs,
            temp_cropped_src->sin_shape_para_angs, cropped_src->n_shapes*num_time_steps);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(cropped_src->shape_gaussbeam_has,
            temp_cropped_src->shape_gaussbeam_has, cropped_src->n_shapes*num_time_steps);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(cropped_src->shape_gaussbeam_decs,
            temp_cropped_src->shape_gaussbeam_decs, cropped_src->n_shapes*num_time_steps);


      //THESE are the arrays that should actually be split up
      //With the way I've set up the sky model creation, the coeff splitting
      //should yield arrays which are the index of cropped_src integer divided
      //by number of coeffs per shapelet

      float *expec_repeat_shape_array = malloc(expected_n_coeffs*sizeof(float));

      int new_ind = 0;
      for (int orig_index = chunk_coeff_ind; orig_index < chunk_coeff_ind + expected_n_coeffs; orig_index++) {
        expec_repeat_shape_array[new_ind] = (float)(orig_index / num_coeff_per_shape);
        new_ind ++;
      }

      TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_repeat_shape_array,
                      temp_cropped_src->shape_coeffs, expected_n_coeffs);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_repeat_shape_array,
                      temp_cropped_src->shape_n1s, expected_n_coeffs);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_repeat_shape_array,
                      temp_cropped_src->shape_n2s, expected_n_coeffs);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_repeat_shape_array,
                      temp_cropped_src->shape_param_indexes, expected_n_coeffs);

      free(expec_repeat_shape_array);

    }

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
