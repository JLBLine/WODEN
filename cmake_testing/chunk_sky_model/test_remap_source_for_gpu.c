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


//Check the remapping has worked for a set of components_t
void check_components_remapping(source_t *chunked_source, source_t *remapped_source,
                                e_component_type comptype, int num_time_steps,
                                int num_list_values, int * list_accum){

  int n_powers = 0;
  int n_curves = 0;
  int n_lists = 0;

  int new_comp_index = 0;

  components_t *remapped_comps = NULL;
  components_t *chunked_comps = NULL;

  if (comptype == POINT) {

    remapped_comps = &remapped_source->point_components;
    chunked_comps = &chunked_source->point_components;

    n_powers = chunked_source->n_point_powers;
    n_curves = chunked_source->n_point_curves;
    n_lists = chunked_source->n_point_lists;
  }
  else if (comptype == GAUSSIAN) {

    remapped_comps = &remapped_source->gauss_components;
    chunked_comps = &chunked_source->gauss_components;

    n_powers = chunked_source->n_gauss_powers;
    n_curves = chunked_source->n_gauss_curves;
    n_lists = chunked_source->n_gauss_lists;
  }

  //These things should have just been pointer copied

//   remapped_components->power_ref_stokesI = chunked_components->power_ref_stokesI;
// remapped_components->power_ref_stokesQ = chunked_components->power_ref_stokesQ;
// remapped_components->power_ref_stokesU = chunked_components->power_ref_stokesU;
// remapped_components->power_ref_stokesV = chunked_components->power_ref_stokesV;
// remapped_components->power_ref_freqs = chunked_components->power_ref_freqs;
// remapped_components->power_SIs = chunked_components->power_SIs;


  //Things that were actually remapped

  double orig_value;

  //Check the POWER_LAW flux type models
  for (int pow_ind = 0; pow_ind < n_powers; pow_ind++) {

    orig_value = (double)chunked_comps->power_comp_inds[pow_ind];

    TEST_ASSERT_EQUAL_INT(new_comp_index, remapped_comps->power_comp_inds[pow_ind]);

    TEST_ASSERT_EQUAL_DOUBLE(orig_value, remapped_comps->ras[new_comp_index]);
    TEST_ASSERT_EQUAL_DOUBLE(orig_value, remapped_comps->decs[new_comp_index]);

    if (comptype == GAUSSIAN) {
      TEST_ASSERT_EQUAL_USER((user_precision_t)orig_value, remapped_comps->majors[new_comp_index]);
      TEST_ASSERT_EQUAL_USER((user_precision_t)orig_value, remapped_comps->minors[new_comp_index]);
      TEST_ASSERT_EQUAL_USER((user_precision_t)orig_value, remapped_comps->pas[new_comp_index]);
    }

    for (int time_ind = 0; time_ind < num_time_steps; time_ind++) {
      TEST_ASSERT_EQUAL_USER((user_precision_t)orig_value,
                remapped_comps->azs[new_comp_index*num_time_steps + time_ind]);
      TEST_ASSERT_EQUAL_USER((user_precision_t)orig_value,
                remapped_comps->zas[new_comp_index*num_time_steps + time_ind]);
      TEST_ASSERT_EQUAL_USER((user_precision_t)orig_value,
                remapped_comps->beam_has[new_comp_index*num_time_steps + time_ind]);
      TEST_ASSERT_EQUAL_USER((user_precision_t)orig_value,
                remapped_comps->beam_decs[new_comp_index*num_time_steps + time_ind]);
    }

    new_comp_index += 1;
  }

  //Check the POWER_LAW flux type models
  for (int cur_ind = 0; cur_ind < n_curves; cur_ind++) {

    orig_value = (double)chunked_comps->curve_comp_inds[cur_ind];

    // printf("DOUBLE check %d %d %d\n",new_comp_index, remapped_comps->curve_comp_inds[cur_ind], cur_ind );

    TEST_ASSERT_EQUAL_INT(new_comp_index, remapped_comps->curve_comp_inds[cur_ind]);

    TEST_ASSERT_EQUAL_DOUBLE(orig_value, remapped_comps->ras[new_comp_index]);
    TEST_ASSERT_EQUAL_DOUBLE(orig_value, remapped_comps->decs[new_comp_index]);

    if (comptype == GAUSSIAN) {
      TEST_ASSERT_EQUAL_USER((user_precision_t)orig_value, remapped_comps->majors[new_comp_index]);
      TEST_ASSERT_EQUAL_USER((user_precision_t)orig_value, remapped_comps->minors[new_comp_index]);
      TEST_ASSERT_EQUAL_USER((user_precision_t)orig_value, remapped_comps->pas[new_comp_index]);
    }

    for (int time_ind = 0; time_ind < num_time_steps; time_ind++) {
      TEST_ASSERT_EQUAL_USER((user_precision_t)orig_value,
                remapped_comps->azs[new_comp_index*num_time_steps + time_ind]);
      TEST_ASSERT_EQUAL_USER((user_precision_t)orig_value,
                remapped_comps->zas[new_comp_index*num_time_steps + time_ind]);
      TEST_ASSERT_EQUAL_USER((user_precision_t)orig_value,
                remapped_comps->beam_has[new_comp_index*num_time_steps + time_ind]);
      TEST_ASSERT_EQUAL_USER((user_precision_t)orig_value,
                remapped_comps->beam_decs[new_comp_index*num_time_steps + time_ind]);
    }

    new_comp_index += 1;
  }

  //Check the POWER_LAW flux type models
  for (int list_ind = 0; list_ind < n_lists; list_ind++) {


    orig_value = (double)chunked_comps->list_comp_inds[list_ind];

    TEST_ASSERT_EQUAL_INT(new_comp_index, remapped_comps->list_comp_inds[list_ind]);
    TEST_ASSERT_EQUAL_INT(list_ind*num_list_values,
                          remapped_comps->list_start_indexes[list_ind]);
    TEST_ASSERT_EQUAL_INT(num_list_values,
                          remapped_comps->num_list_values[list_ind]);

    TEST_ASSERT_EQUAL_DOUBLE(orig_value, remapped_comps->ras[new_comp_index]);
    TEST_ASSERT_EQUAL_DOUBLE(orig_value, remapped_comps->decs[new_comp_index]);

    if (comptype == GAUSSIAN) {
      TEST_ASSERT_EQUAL_USER((user_precision_t)orig_value,
                                        remapped_comps->majors[new_comp_index]);
      TEST_ASSERT_EQUAL_USER((user_precision_t)orig_value,
                                        remapped_comps->minors[new_comp_index]);
      TEST_ASSERT_EQUAL_USER((user_precision_t)orig_value,
                                        remapped_comps->pas[new_comp_index]);
    }



    for (int lis = 0; lis < num_list_values; lis++) {
      TEST_ASSERT_EQUAL_DOUBLE((double)* list_accum,
                    remapped_comps->list_freqs[list_ind*num_list_values + lis]);
      TEST_ASSERT_EQUAL_USER((user_precision_t)* list_accum,
                    remapped_comps->list_stokesI[list_ind*num_list_values + lis]);
      TEST_ASSERT_EQUAL_USER((user_precision_t)* list_accum,
                    remapped_comps->list_stokesQ[list_ind*num_list_values + lis]);
      TEST_ASSERT_EQUAL_USER((user_precision_t)* list_accum,
                    remapped_comps->list_stokesU[list_ind*num_list_values + lis]);
      TEST_ASSERT_EQUAL_USER((user_precision_t)* list_accum,
                    remapped_comps->list_stokesV[list_ind*num_list_values + lis]);

    }

    for (int time_ind = 0; time_ind < num_time_steps; time_ind++) {
      TEST_ASSERT_EQUAL_USER((user_precision_t)orig_value,
                remapped_comps->azs[new_comp_index*num_time_steps + time_ind]);
      TEST_ASSERT_EQUAL_USER((user_precision_t)orig_value,
                remapped_comps->zas[new_comp_index*num_time_steps + time_ind]);
      TEST_ASSERT_EQUAL_USER((user_precision_t)orig_value,
                remapped_comps->beam_has[new_comp_index*num_time_steps + time_ind]);
      TEST_ASSERT_EQUAL_USER((user_precision_t)orig_value,
                remapped_comps->beam_decs[new_comp_index*num_time_steps + time_ind]);
    }

    new_comp_index += 1;
    * list_accum += 1;
  }

}

//Check the remapping has worked correctly for the whole source_t
void check_source_remapping(source_t *chunked_source, source_t *remapped_source,
                            int num_time_steps, int num_list_values,
                            int * point_list_accum, int * gauss_list_accum){

  check_components_remapping(chunked_source, remapped_source, POINT, num_time_steps, num_list_values, point_list_accum);
  check_components_remapping(chunked_source, remapped_source, GAUSSIAN, num_time_steps, num_list_values, gauss_list_accum);

  TEST_ASSERT_EQUAL_INT(chunked_source->n_points,
                        remapped_source->n_points);
  TEST_ASSERT_EQUAL_INT(chunked_source->n_point_lists,
                        remapped_source->n_point_lists);
  TEST_ASSERT_EQUAL_INT(chunked_source->n_point_powers,
                        remapped_source->n_point_powers);
  TEST_ASSERT_EQUAL_INT(chunked_source->n_point_curves,
                        remapped_source->n_point_curves);
  TEST_ASSERT_EQUAL_INT(chunked_source->n_gauss,
                        remapped_source->n_gauss);
  TEST_ASSERT_EQUAL_INT(chunked_source->n_gauss_lists,
                        remapped_source->n_gauss_lists);
  TEST_ASSERT_EQUAL_INT(chunked_source->n_gauss_powers,
                        remapped_source->n_gauss_powers);
  TEST_ASSERT_EQUAL_INT(chunked_source->n_gauss_curves,
                        remapped_source->n_gauss_curves);

}



/*
Gotta remap the chunked sources in a way that allows easy indexing on the GPU,
without using massive amounts of memory
*/

void test_remap_source_for_gpu(int chunking_size,
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

  // printf("====================================================\n");
  // printf("cropped_src->n_comps %d\n", cropped_src->n_comps);
  // printf("cropped_src->n_points %d\n", cropped_src->n_points);
  // printf("cropped_src->n_point_lists %d\n", cropped_src->n_point_lists);
  // printf("cropped_src->n_point_powers %d\n", cropped_src->n_point_powers);
  // printf("cropped_src->n_point_curves %d\n", cropped_src->n_point_curves);
  // printf("cropped_src->n_gauss %d\n", cropped_src->n_gauss);
  // printf("cropped_src->n_gauss_lists %d\n", cropped_src->n_gauss_lists);
  // printf("cropped_src->n_gauss_powers %d\n", cropped_src->n_gauss_powers);
  // printf("cropped_src->n_gauss_curves %d\n", cropped_src->n_gauss_curves);
  // printf("cropped_src->n_shapes %d\n", cropped_src->n_shapes);
  // printf("cropped_src->n_shape_lists %d\n", cropped_src->n_shape_lists);
  // printf("cropped_src->n_shape_powers %d\n", cropped_src->n_shape_powers);
  // printf("cropped_src->n_shape_curves %d\n", cropped_src->n_shape_curves);
  // printf("cropped_src->n_shape_coeffs %d\n", cropped_src->n_shape_coeffs);
  // printf("====================================================\n");


  //Set some dummy baseline, frequency simulation settings
  // int num_visis = num_baselines * num_freqs * num_time_steps;

  woden_settings_t *woden_settings = malloc(sizeof(woden_settings_t));
  woden_settings->num_freqs = num_freqs;
  woden_settings->num_time_steps = num_time_steps;
  woden_settings->chunking_size = chunking_size;
  woden_settings->num_baselines = num_baselines;


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

  //Need to keep track of the number of LIST type flux components so we
  //can predict array values
  int point_list_accum = 0;
  int gauss_list_accum = 0;

  int num_pg_chunks = num_point_chunks + num_gauss_chunks;

  for (int chunk_ind = 0; chunk_ind < num_pg_chunks; chunk_ind++) {

    source_t *remapped_source = malloc(sizeof(source_t));

    //Code we are testing
    remap_source_for_gpu(remapped_source, &cropped_sky_models->sources[chunk_ind],
                         num_time_steps, MWA_ANALY);

    // printf("remapped_source->n_points %d\n", remapped_source->n_points );
    // printf("remapped_source->n_point_lists %d\n", remapped_source->n_point_lists );
    // printf("remapped_source->n_point_powers %d\n", remapped_source->n_point_powers );
    // printf("remapped_source->n_point_curves %d\n", remapped_source->n_point_curves );
    // printf("remapped_source->n_gauss %d\n", remapped_source->n_gauss );
    // printf("remapped_source->n_gauss_lists %d\n", remapped_source->n_gauss_lists );
    // printf("remapped_source->n_gauss_powers %d\n", remapped_source->n_gauss_powers );
    // printf("remapped_source->n_gauss_curves %d\n", remapped_source->n_gauss_curves );
    // printf("----------------------------------------------------\n");

    //Check the code has worked correctly
    check_source_remapping(&cropped_sky_models->sources[chunk_ind],
                          remapped_source, num_time_steps, num_list_values,
                          &point_list_accum, &gauss_list_accum);

    //free up things that are going to have a malloc in the next loop
    free_remapped_source_for_gpu(remapped_source, MWA_ANALY);
  }

  for (int chunk_ind = 0; chunk_ind < num_shape_chunks; chunk_ind++) {

    source_t *remapped_source = malloc(sizeof(source_t));

    //Code we are testing
    remap_source_for_gpu(remapped_source, &cropped_sky_models->sources[num_pg_chunks +chunk_ind],
                         num_time_steps, MWA_ANALY);

    //The SHAPELET part should work out same as `create_chunked_sky_models`
    //so we can use the same test here
    check_shapelet_chunking(chunk_ind, comps_per_chunk, num_time_steps,
                            num_coeff_per_shape, cropped_src,
                            &cropped_sky_models->sources[num_pg_chunks + chunk_ind]);

    //free up things that are going to have a malloc in the next loop
    free_remapped_source_for_gpu(remapped_source, MWA_ANALY);
  }

}


void test_remap_source_for_gpu_P3_G0000_S00_00_C1_Time001(void) {
  int num_points = 3;
  int num_gauss = 0;
  int num_shapes = 0;
  int num_coeff_per_shape = 0;
  int chunking_size = 2;

  int num_time_steps = 1;
  int num_baselines = 1;
  int num_freqs = 1;

  int num_list_values = 3;
  test_remap_source_for_gpu(chunking_size, num_points, num_gauss,
                                 num_shapes,  num_coeff_per_shape,
                                 num_list_values,
                                 num_time_steps, num_baselines, num_freqs);
}



void test_remap_source_for_gpu_P4544_G1736_S20_14_C1e8_Time014(void){

  int num_points = 4544;
  int num_gauss = 1736;
  int num_shapes = 20;
  int num_coeff_per_shape = 14;

  int chunking_size = 1e8;

  int num_time_steps = 14;
  int num_baselines = 8128;
  int num_freqs = 16;

  int num_list_values = 4;
  test_remap_source_for_gpu(chunking_size, num_points, num_gauss,
                                 num_shapes,  num_coeff_per_shape,
                                 num_list_values,
                                 num_time_steps, num_baselines, num_freqs);
}

void test_remap_source_for_gpu_P1000_G0000_S00_00_C1e8_Time014(void){

  int num_points = 1000;
  int num_gauss = 0;
  int num_shapes = 0;
  int num_coeff_per_shape = 0;

  int chunking_size = 1e9;

  int num_time_steps = 14;
  int num_baselines = 8128;
  int num_freqs = 16;

  int num_list_values = 6;
  test_remap_source_for_gpu(chunking_size, num_points, num_gauss,
                                 num_shapes,  num_coeff_per_shape,
                                 num_list_values,
                                 num_time_steps, num_baselines, num_freqs);
}

void test_remap_source_for_gpu_P0000_G1000_S00_00_C1e8_Time014(void){

  int num_points = 0;
  int num_gauss = 1000;
  int num_shapes = 0;
  int num_coeff_per_shape = 0;

  int chunking_size = 1e9;

  int num_time_steps = 14;
  int num_baselines = 8128;
  int num_freqs = 16;

  int num_list_values = 5;
  test_remap_source_for_gpu(chunking_size, num_points, num_gauss,
                                 num_shapes,  num_coeff_per_shape,
                                 num_list_values,
                                 num_time_steps, num_baselines, num_freqs);
}

void test_remap_source_for_gpu_P0000_G0000_S67_78_C1e8_Time014(void){

  int num_points = 0;
  int num_gauss = 0;
  int num_shapes = 67;
  int num_coeff_per_shape = 78;

  int chunking_size = 1e9;

  int num_time_steps = 14;
  int num_baselines = 8128;
  int num_freqs = 16;

  int num_list_values = 2;
  test_remap_source_for_gpu(chunking_size, num_points, num_gauss,
                                 num_shapes,  num_coeff_per_shape,
                                 num_list_values,
                                 num_time_steps, num_baselines, num_freqs);
}

void test_remap_source_for_gpu_P0050_G0050_S10_25_C4576_Time014(void){

  int num_points = 50;
  int num_gauss = 50;
  int num_shapes = 10;
  int num_coeff_per_shape = 25;

  int chunking_size = 4576;

  int num_time_steps = 14;
  int num_baselines = 8128;
  int num_freqs = 16;

  int num_list_values = 8;
  test_remap_source_for_gpu(chunking_size, num_points, num_gauss,
                                 num_shapes,  num_coeff_per_shape,
                                 num_list_values,
                                 num_time_steps, num_baselines, num_freqs);
}

void test_remap_source_for_gpu_P87760_G12207_S121_678_C1e10_Time056(void){

  int num_points = 87760;
  int num_gauss = 12207;
  int num_shapes = 121;
  int num_coeff_per_shape = 678;

  long int chunking_size = 1e10;

  int num_time_steps = 56;
  int num_baselines = 8128;
  int num_freqs = 32;

  int num_list_values = 16;
  test_remap_source_for_gpu(chunking_size, num_points, num_gauss,
                                 num_shapes,  num_coeff_per_shape,
                                 num_list_values,
                                 num_time_steps, num_baselines, num_freqs);
}

//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_remap_source_for_gpu_P3_G0000_S00_00_C1_Time001);
    RUN_TEST(test_remap_source_for_gpu_P1000_G0000_S00_00_C1e8_Time014);
    RUN_TEST(test_remap_source_for_gpu_P0000_G1000_S00_00_C1e8_Time014);
    RUN_TEST(test_remap_source_for_gpu_P0000_G0000_S67_78_C1e8_Time014);
    RUN_TEST(test_remap_source_for_gpu_P0050_G0050_S10_25_C4576_Time014);
    RUN_TEST(test_remap_source_for_gpu_P4544_G1736_S20_14_C1e8_Time014);
    RUN_TEST(test_remap_source_for_gpu_P87760_G12207_S121_678_C1e10_Time056);

    return UNITY_END();
}
