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
                                         int num_time_steps,
                                         int num_baselines, int num_freqs) {

  //Not testing SHAPELETs, set these to zero
  int num_shapes = 0;
  int num_coeff_per_shape = 0;
  catsource_t *cropped_src = make_sky_model(num_points, num_gauss, num_shapes,
                                            num_coeff_per_shape, num_time_steps);

  //Set some dummy baseline, frequency simulation settings
  int num_visis = num_baselines * num_freqs * num_time_steps;

  woden_settings_t *woden_settings = malloc(sizeof(woden_settings_t));
  woden_settings->num_freqs = num_freqs;
  woden_settings->num_time_steps = num_time_steps;

  // //
  int num_comps_to_chunk = cropped_src->n_points + cropped_src->n_gauss;
  int comps_per_chunk = (int)floorf((float)chunking_size / (float)(num_baselines * num_freqs * num_time_steps));
  int num_chunks = (int)ceilf((float)num_comps_to_chunk / (float)comps_per_chunk);
  //
  printf("Num visis %d\n",num_visis );
  printf("Chunking size %d\n",chunking_size );
  printf("Comps per chunk %d\n", comps_per_chunk );
  printf("Num chunks %d\n", num_chunks );
  // printf("Is num_visis*num_comps %d < num_visis*comps_per_chunk*num_chunks %d\n",
  //        num_visis*num_comps_to_chunk, num_visis*comps_per_chunk*num_chunks);

  catsource_t *temp_cropped_src = malloc(sizeof(catsource_t));
  //
  //Counters to see how far many point/gauss have already been assigned
  //to previous chunks
  int point_accum = 0;
  int gauss_accum = 0;

  for (int chunk_ind = 0; chunk_ind < num_chunks; chunk_ind++) {

    fill_chunk_src_with_pointgauss(temp_cropped_src, cropped_src, chunk_ind,
                                   comps_per_chunk, woden_settings);

    // printf("Chunk %d has %d points\n",chunk_ind, temp_cropped_src->n_points );

    //How far through all components are we with this chunk
    int chunk_comp_ind = chunk_ind*comps_per_chunk;
    //How many POINT sources are left after this chunk index
    int point_remainder = cropped_src->n_points - chunk_comp_ind;

    //Things to use in logic below
    int expected_n_points = 0;
    int expected_n_gauss = 0;
    int chunk_remainder = 0;
    int gauss_remainder = 0;

    if (point_remainder > 0) { //There are POINT sources, how many should there be?

      if (point_remainder >= comps_per_chunk) { //POINT sources fill the whole chunk
        expected_n_points = comps_per_chunk;
      } else { //Not enough point sources to fill the chunk
        expected_n_points = point_remainder;

        //How many component can fit into rest of chunk?
        chunk_remainder = comps_per_chunk - expected_n_points;
        //See if GAUSS fill the remainder or are smaller
        if (cropped_src->n_gauss > chunk_remainder) {
          expected_n_gauss = chunk_remainder;
        } else if ( cropped_src->n_gauss > 0) {
          expected_n_gauss = cropped_src->n_gauss;
        } else {
          expected_n_gauss = 0;
        }
      } //END there are POINTs, not enough point sources to fill the chunk

    } else { //There are no POINT sources in chunk
      expected_n_points = 0;

      gauss_remainder = num_comps_to_chunk - chunk_comp_ind;

      if (gauss_remainder >= comps_per_chunk) { //POINT sources fill the whole chunk
        expected_n_gauss = comps_per_chunk;
      } else { //Not enough point sources to fill the chunk
        expected_n_gauss = gauss_remainder;
      }
    } //END there are no POINT sources

    // printf("expect point, gauss %d %d\n",expected_n_points,expected_n_gauss );
    if (expected_n_points > 0) {
      // printf("Found POINT, chunk_ind %d expected_n_points %d\n",chunk_ind, expected_n_points );

      TEST_ASSERT_EQUAL_INT(expected_n_points, temp_cropped_src->n_points);

      float *expec_index_point_array = malloc(expected_n_points*sizeof(float));
      make_index_array(expec_index_point_array, expected_n_points, point_accum);

      float *expec_repeat_point_array = malloc(num_time_steps*expected_n_points*sizeof(float));
      make_repeat_array(expec_repeat_point_array, expected_n_points,
                        num_time_steps, point_accum);

      //Check POINT source params were split correctly
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_index_point_array,
                              temp_cropped_src->point_ras, expected_n_points);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_index_point_array,
                              temp_cropped_src->point_decs, expected_n_points);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_index_point_array,
                              temp_cropped_src->point_ref_freqs, expected_n_points);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_index_point_array,
                              temp_cropped_src->point_ref_stokesI, expected_n_points);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_index_point_array,
                              temp_cropped_src->point_ref_stokesQ, expected_n_points);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_index_point_array,
                              temp_cropped_src->point_ref_stokesU, expected_n_points);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_index_point_array,
                              temp_cropped_src->point_ref_stokesV, expected_n_points);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_index_point_array,
                              temp_cropped_src->point_SIs, expected_n_points);
      //
      //Check POINT source prinary beam params were split correctly
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_repeat_point_array,
                      temp_cropped_src->point_azs, expected_n_points*num_time_steps);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_repeat_point_array,
                      temp_cropped_src->point_zas, expected_n_points*num_time_steps);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_repeat_point_array,
            temp_cropped_src->cos_point_para_angs, expected_n_points*num_time_steps);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_repeat_point_array,
            temp_cropped_src->sin_point_para_angs, expected_n_points*num_time_steps);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_repeat_point_array,
            temp_cropped_src->point_gaussbeam_has, expected_n_points*num_time_steps);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_repeat_point_array,
            temp_cropped_src->point_gaussbeam_decs, expected_n_points*num_time_steps);

      free(expec_index_point_array);
      free(expec_repeat_point_array);

    }

    if (expected_n_gauss > 0) {
      // printf("Found GAUSS, chunk_ind %d expected_n_gausss %d\n",chunk_ind, expected_n_gauss );
      TEST_ASSERT_EQUAL_INT(expected_n_gauss, temp_cropped_src->n_gauss);

      float *expec_index_gauss_array = malloc(expected_n_gauss*sizeof(float));
      make_index_array(expec_index_gauss_array, expected_n_gauss, gauss_accum);

      float *expec_repeat_gauss_array = malloc(num_time_steps*expected_n_gauss*sizeof(float));
      make_repeat_array(expec_repeat_gauss_array, expected_n_gauss,
                        num_time_steps, gauss_accum);

      //Check GAUSS source params were split correctly
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_index_gauss_array,
                                    temp_cropped_src->gauss_ras, expected_n_gauss);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_index_gauss_array,
                                    temp_cropped_src->gauss_decs, expected_n_gauss);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_index_gauss_array,
                                    temp_cropped_src->gauss_ref_freqs, expected_n_gauss);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_index_gauss_array,
                                    temp_cropped_src->gauss_ref_stokesI, expected_n_gauss);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_index_gauss_array,
                                    temp_cropped_src->gauss_ref_stokesQ, expected_n_gauss);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_index_gauss_array,
                                    temp_cropped_src->gauss_ref_stokesU, expected_n_gauss);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_index_gauss_array,
                                    temp_cropped_src->gauss_ref_stokesV, expected_n_gauss);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_index_gauss_array,
                                    temp_cropped_src->gauss_SIs, expected_n_gauss);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_index_gauss_array,
                                    temp_cropped_src->gauss_pas, expected_n_gauss);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_index_gauss_array,
                                    temp_cropped_src->gauss_majors, expected_n_gauss);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_index_gauss_array,
                                    temp_cropped_src->gauss_minors, expected_n_gauss);

      //Check GAUSS source prinary beam params were split correctly
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_repeat_gauss_array,
                      temp_cropped_src->gauss_azs, expected_n_gauss*num_time_steps);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_repeat_gauss_array,
                      temp_cropped_src->gauss_zas, expected_n_gauss*num_time_steps);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_repeat_gauss_array,
            temp_cropped_src->cos_gauss_para_angs, expected_n_gauss*num_time_steps);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_repeat_gauss_array,
            temp_cropped_src->sin_gauss_para_angs, expected_n_gauss*num_time_steps);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_repeat_gauss_array,
            temp_cropped_src->gauss_gaussbeam_has, expected_n_gauss*num_time_steps);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_repeat_gauss_array,
            temp_cropped_src->gauss_gaussbeam_decs, expected_n_gauss*num_time_steps);

      free(expec_index_gauss_array);
      free(expec_repeat_gauss_array);
    }
    point_accum += expected_n_points;
    gauss_accum += expected_n_gauss;


  } //END iteration over all chunks
  free_sky_model(cropped_src);
}


void test_fill_chunk_src_with_pointgauss_Point100_Gauss000_Chunk1000_Time004(void){
  int chunking_size = 1000;
  int num_points = 100;
  int num_gauss = 0;
  int num_time_steps = 4;
  int num_baselines = 5;
  int num_freqs = 2;

  test_fill_chunk_src_with_pointgauss(chunking_size, num_points, num_gauss,
                                      num_time_steps, num_baselines, num_freqs);
}

void test_fill_chunk_src_with_pointgauss_Point000_Gauss100_Chunk1000_Time004(void){
  int chunking_size = 1000;
  int num_points = 0;
  int num_gauss = 100;
  int num_time_steps = 4;
  int num_baselines = 5;
  int num_freqs = 2;

  test_fill_chunk_src_with_pointgauss(chunking_size, num_points, num_gauss,
                                        num_time_steps, num_baselines, num_freqs);
}

void test_fill_chunk_src_with_pointgauss_Point100_Gauss100_Chunk1000_Time004(void){
  int chunking_size = 1000;
  int num_points = 100;
  int num_gauss = 100;
  int num_time_steps = 4;
  int num_baselines = 5;
  int num_freqs = 2;

  test_fill_chunk_src_with_pointgauss(chunking_size, num_points, num_gauss,
                                        num_time_steps, num_baselines, num_freqs);
}

void test_fill_chunk_src_with_pointgauss_Point100_Gauss000_Chunk1000_Time003(void){
  int chunking_size = 1000;
  int num_points = 100;
  int num_gauss = 0;
  int num_time_steps = 3;
  int num_baselines = 5;
  int num_freqs = 2;

  test_fill_chunk_src_with_pointgauss(chunking_size, num_points, num_gauss,
                                        num_time_steps, num_baselines, num_freqs);
}

void test_fill_chunk_src_with_pointgauss_Point000_Gauss100_Chunk1000_Time003(void){
  int chunking_size = 1000;
  int num_points = 0;
  int num_gauss = 100;
  int num_time_steps = 3;
  int num_baselines = 5;
  int num_freqs = 2;

  test_fill_chunk_src_with_pointgauss(chunking_size, num_points, num_gauss,
                                        num_time_steps, num_baselines, num_freqs);
}

void test_fill_chunk_src_with_pointgauss_Point100_Gauss100_Chunk1000_Time003(void){
  int chunking_size = 1000;
  int num_points = 100;
  int num_gauss = 100;
  int num_time_steps = 3;
  int num_baselines = 5;
  int num_freqs = 2;

  test_fill_chunk_src_with_pointgauss(chunking_size, num_points, num_gauss,
                                        num_time_steps, num_baselines, num_freqs);
}

void test_fill_chunk_src_with_pointgauss_Point100_Gauss000_Chunk0173_Time005(void){
  int chunking_size = 173;
  int num_points = 100;
  int num_gauss = 0;
  int num_time_steps = 5;
  int num_baselines = 5;
  int num_freqs = 2;

  test_fill_chunk_src_with_pointgauss(chunking_size, num_points, num_gauss,
                                        num_time_steps, num_baselines, num_freqs);
}

void test_fill_chunk_src_with_pointgauss_Point000_Gauss100_Chunk0173_Time005(void){
  int chunking_size = 173;
  int num_points = 0;
  int num_gauss = 100;
  int num_time_steps = 5;
  int num_baselines = 5;
  int num_freqs = 2;

  test_fill_chunk_src_with_pointgauss(chunking_size, num_points, num_gauss,
                                        num_time_steps, num_baselines, num_freqs);
}

void test_fill_chunk_src_with_pointgauss_Point100_Gauss100_Chunk0173_Time005(void){
  int chunking_size = 173;
  int num_points = 100;
  int num_gauss = 100;
  int num_time_steps = 5;
  int num_baselines = 5;
  int num_freqs = 2;

  test_fill_chunk_src_with_pointgauss(chunking_size, num_points, num_gauss,
                                        num_time_steps, num_baselines, num_freqs);
}

void test_fill_chunk_src_with_pointgauss_Point10743_Gauss00000_Chunk4345_Time012(void){
  int chunking_size = 4345;
  int num_points = 10743;
  int num_gauss = 0;
  int num_time_steps = 12;
  int num_baselines = 5;
  int num_freqs = 2;

  test_fill_chunk_src_with_pointgauss(chunking_size, num_points, num_gauss,
                                        num_time_steps, num_baselines, num_freqs);
}

void test_fill_chunk_src_with_pointgauss_Point00000_Gauss16789_Chunk4345_Time012(void){
  int chunking_size = 4345;
  int num_points = 0;
  int num_gauss = 16789;
  int num_time_steps = 12;
  int num_baselines = 5;
  int num_freqs = 2;

  test_fill_chunk_src_with_pointgauss(chunking_size, num_points, num_gauss,
                                        num_time_steps, num_baselines, num_freqs);
}

void test_fill_chunk_src_with_pointgauss_Point10743_Gauss16789_Chunk4345_Time012(void){
  int chunking_size = 4345;
  int num_points = 10743;
  int num_gauss = 16789;
  int num_time_steps = 12;
  int num_baselines = 5;
  int num_freqs = 2;

  test_fill_chunk_src_with_pointgauss(chunking_size, num_points, num_gauss,
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

  test_fill_chunk_src_with_pointgauss(chunking_size, num_points, num_gauss,
                                        num_time_steps, num_baselines, num_freqs);
}

//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    // RUN_TEST(test_fill_chunk_src_with_pointgauss_Point100_Gauss000_Chunk1000_Time004);
    // RUN_TEST(test_fill_chunk_src_with_pointgauss_Point000_Gauss100_Chunk1000_Time004);
    // RUN_TEST(test_fill_chunk_src_with_pointgauss_Point100_Gauss100_Chunk1000_Time004);
    // RUN_TEST(test_fill_chunk_src_with_pointgauss_Point100_Gauss000_Chunk1000_Time003);
    // RUN_TEST(test_fill_chunk_src_with_pointgauss_Point000_Gauss100_Chunk1000_Time003);
    // RUN_TEST(test_fill_chunk_src_with_pointgauss_Point100_Gauss100_Chunk1000_Time003);
    // RUN_TEST(test_fill_chunk_src_with_pointgauss_Point100_Gauss000_Chunk0173_Time005);
    // RUN_TEST(test_fill_chunk_src_with_pointgauss_Point000_Gauss100_Chunk0173_Time005);
    // RUN_TEST(test_fill_chunk_src_with_pointgauss_Point100_Gauss100_Chunk0173_Time005);
    // RUN_TEST(test_fill_chunk_src_with_pointgauss_Point10743_Gauss00000_Chunk4345_Time012);
    // RUN_TEST(test_fill_chunk_src_with_pointgauss_Point00000_Gauss16789_Chunk4345_Time012);
    RUN_TEST(test_fill_chunk_src_with_pointgauss_Point4544_Gauss1736_Chunk1e8_Time014);

    return UNITY_END();
}
