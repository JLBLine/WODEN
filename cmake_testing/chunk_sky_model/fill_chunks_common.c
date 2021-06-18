#include <stdlib.h>
#include <math.h>
#include <unity.h>

#include "constants.h"
#include "woden_struct_defs.h"

/*
The GPU will never have enough memory to simulate millions of COMPONENTs at
once, so we need to split the sky model up before running the simulation.
`fill_chunk_src_*` functions do just that, given the sky model above the
horizon (`cropped_src`) and a chunking index, they return a portion of the
sky model containing specific COMPONENT types.

As SHAPELETs take more memory, we separate them out from POINT and GAUSSIANs,
so have two chunking functions.
*/

//Populate an array with index values. Add an offset of `offset`
void make_index_array(float *index_array, int array_length, int offset) {
  for (int index = 0; index < array_length; index++) {
    index_array[index] = (float)(index + offset);
  }
}

/*
Given a maximum value `num_value`, populate an array from 0 to `num_value`-1,
repeating each value `num_repeat times`, also add an `offset` value
i.e. if num_value = 2, num_repeat_times = 3, offset=4 makes an array
repeat_array[] = {4, 4, 4, 5, 5, 5}

*/
void make_repeat_array(float *repeat_array, int num_value, int num_repeat,
                       int offset) {
  int index = 0;
  for (int value = 0; value < num_value; value++) {
    for (int repeat = 0; repeat < num_repeat; repeat++) {
      repeat_array[index] = offset + value;
      index ++;
    }

  }
}


/*
Make the polpulated catsource_t struct. For each COMPONENT type, assign the
index value of the COMPONENT to each array, and have a set number of coeffs,
n1s, and n2s per SHAPELET source. Should make it easy to test whether chunking
has worked sensibly

Values associated with beams (para_angs, gaussbeam, zas, azs) are for every
time step, so give them a longer index array
*/
catsource_t * make_sky_model(int num_points, int num_gauss,
                             int num_shapes, int num_coeff_per_shape,
                             int num_time_steps) {

  catsource_t *cropped_src = malloc(sizeof(catsource_t));

  //Overall COMPONENT stats
  cropped_src->n_comps = num_points*num_gauss*num_shapes;
  cropped_src->n_points= num_points;
  cropped_src->n_gauss = num_gauss;
  cropped_src->n_shapes = num_shapes;
  cropped_src->n_shape_coeffs = num_shapes*num_coeff_per_shape;

  //Make index array for POINT source
  float *index_point_array = malloc(num_points*sizeof(float));
  make_index_array(index_point_array, num_points, 0);

  //Populate POINT intrinsic properties
  cropped_src->point_ras = index_point_array;
  cropped_src->point_decs = index_point_array;
  cropped_src->point_ref_freqs = index_point_array;
  cropped_src->point_ref_stokesI = index_point_array;
  cropped_src->point_ref_stokesQ = index_point_array;
  cropped_src->point_ref_stokesU = index_point_array;
  cropped_src->point_ref_stokesV = index_point_array;
  cropped_src->point_SIs = index_point_array;

  //Make repeating array for POINT primary beam related values
  float *repeat_point_array = malloc(num_time_steps*num_points*sizeof(float));
  make_repeat_array(repeat_point_array, num_points, num_time_steps, 0);
  cropped_src->point_azs = repeat_point_array;
  cropped_src->point_zas = repeat_point_array;
  cropped_src->cos_point_para_angs = repeat_point_array;
  cropped_src->sin_point_para_angs = repeat_point_array;
  cropped_src->point_gaussbeam_has = repeat_point_array;
  cropped_src->point_gaussbeam_decs = repeat_point_array;

  //Repeat process for GAUSSIAN and SHAPELETs

  float *index_gauss_array = malloc(num_gauss*sizeof(float));
  make_index_array(index_gauss_array, num_gauss, 0);

  float *repeat_gauss_array = malloc(num_time_steps*num_gauss*sizeof(float));
  make_repeat_array(repeat_gauss_array, num_gauss, num_time_steps, 0);

  cropped_src->gauss_ras = index_gauss_array;
  cropped_src->gauss_decs = index_gauss_array;
  cropped_src->gauss_ref_freqs = index_gauss_array;
  cropped_src->gauss_ref_stokesI = index_gauss_array;
  cropped_src->gauss_ref_stokesQ = index_gauss_array;
  cropped_src->gauss_ref_stokesU = index_gauss_array;
  cropped_src->gauss_ref_stokesV = index_gauss_array;
  cropped_src->gauss_SIs = index_gauss_array;
  cropped_src->gauss_majors = index_gauss_array;
  cropped_src->gauss_minors = index_gauss_array;
  cropped_src->gauss_pas = index_gauss_array;
  cropped_src->gauss_azs = repeat_gauss_array;
  cropped_src->gauss_zas = repeat_gauss_array;
  cropped_src->cos_gauss_para_angs = repeat_gauss_array;
  cropped_src->sin_gauss_para_angs = repeat_gauss_array;
  cropped_src->gauss_gaussbeam_has = repeat_gauss_array;
  cropped_src->gauss_gaussbeam_decs = repeat_gauss_array;

  float *index_shape_array = malloc(num_shapes*sizeof(float));
  make_index_array(index_shape_array, num_shapes, 0);

  float *repeat_shape_array = malloc(num_time_steps*num_shapes*sizeof(float));
  make_repeat_array(repeat_shape_array, num_shapes, num_time_steps, 0);

  cropped_src->shape_ras = index_shape_array;
  cropped_src->shape_decs = index_shape_array;
  cropped_src->shape_ref_freqs = index_shape_array;
  cropped_src->shape_ref_stokesI = index_shape_array;
  cropped_src->shape_ref_stokesQ = index_shape_array;
  cropped_src->shape_ref_stokesU = index_shape_array;
  cropped_src->shape_ref_stokesV = index_shape_array;
  cropped_src->shape_SIs = index_shape_array;
  cropped_src->shape_majors = index_shape_array;
  cropped_src->shape_minors = index_shape_array;
  cropped_src->shape_pas = index_shape_array;
  cropped_src->shape_azs = repeat_shape_array;
  cropped_src->shape_zas = repeat_shape_array;
  cropped_src->cos_shape_para_angs = repeat_shape_array;
  cropped_src->sin_shape_para_angs = repeat_shape_array;
  cropped_src->shape_gaussbeam_has = repeat_shape_array;
  cropped_src->shape_gaussbeam_decs = repeat_shape_array;

  //These arrays also contain repeating values
  float *repeat_coeffs_array = malloc(num_coeff_per_shape*num_shapes*sizeof(float));
  make_repeat_array(repeat_coeffs_array, num_shapes, num_coeff_per_shape, 0);

  cropped_src->shape_coeffs = repeat_coeffs_array;
  cropped_src->shape_n1s = repeat_coeffs_array;
  cropped_src->shape_n2s = repeat_coeffs_array;
  cropped_src->shape_param_indexes = repeat_coeffs_array;

  return cropped_src;

}

void free_sky_model(catsource_t *cropped_src) {

  //These point to the unique arrays we made in `make_sky_model`
  //so freeing just these is enough
  free(cropped_src->point_ras);
  free(cropped_src->point_azs);
  free(cropped_src->gauss_ras);
  free(cropped_src->gauss_azs);
  free(cropped_src->shape_ras);
  free(cropped_src->shape_azs);
  free(cropped_src->shape_coeffs);

  free(cropped_src);
}


void check_pointgauss_chunking(int chunk_ind, int comps_per_chunk,
                             int num_time_steps,
                             int * point_accum, int * gauss_accum,
                             catsource_t *cropped_src,
                             catsource_t *temp_cropped_src) {

  //How many components overall to be chunked
  int num_comps_to_chunk = cropped_src->n_points + cropped_src->n_gauss;

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
    make_index_array(expec_index_point_array, expected_n_points, * point_accum);

    float *expec_repeat_point_array = malloc(num_time_steps*expected_n_points*sizeof(float));
    make_repeat_array(expec_repeat_point_array, expected_n_points,
                      num_time_steps, * point_accum);

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
    make_index_array(expec_index_gauss_array, expected_n_gauss, * gauss_accum);

    float *expec_repeat_gauss_array = malloc(num_time_steps*expected_n_gauss*sizeof(float));
    make_repeat_array(expec_repeat_gauss_array, expected_n_gauss,
                      num_time_steps, * gauss_accum);

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
  * point_accum += expected_n_points;
  * gauss_accum += expected_n_gauss;

}




void check_shapelet_chunking(int chunk_ind, int coeffs_per_chunk,
                             int num_time_steps,
                             int num_coeff_per_shape,
                             catsource_t *cropped_src,
                             catsource_t *temp_cropped_src){

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
}
