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
void make_index_array(user_precision_t *index_array, int array_length, int offset) {
  for (int index = 0; index < array_length; index++) {
    index_array[index] = (user_precision_t)(index + offset);
  }
}

//Populate an array with index values. Add an offset of `offset`
void make_index_array_double(double *index_array, int array_length, int offset) {
  for (int index = 0; index < array_length; index++) {
    index_array[index] = (double)(index + offset);
  }
}

/*
Given a maximum value `num_value`, populate an array from 0 to `num_value`-1,
repeating each value `num_repeat times`, also add an `offset` value
i.e. if num_value = 2, num_repeat_times = 3, offset=4 makes an array
repeat_array[] = {4, 4, 4, 5, 5, 5}

*/
void make_repeat_array(user_precision_t *repeat_array, int num_value, int num_repeat,
                       int offset) {
  int index = 0;
  for (int value = 0; value < num_value; value++) {
    for (int repeat = 0; repeat < num_repeat; repeat++) {
      repeat_array[index] = offset + value;
      index ++;
    }

  }
}

void make_repeat_array_double(double *repeat_array, int num_value, int num_repeat,
                       int offset) {
  int index = 0;
  for (int value = 0; value < num_value; value++) {
    for (int repeat = 0; repeat < num_repeat; repeat++) {
      repeat_array[index] = (double)(offset + value);
      index ++;
    }

  }
}


/*
Make the polpulated source_t struct. For each COMPONENT type, assign the
index value of the COMPONENT to each array, and have a set number of coeffs,
n1s, and n2s per SHAPELET source. Should make it easy to test whether chunking
has worked sensibly

Values associated with beams (para_angs, gaussbeam, zas, azs) are for every
time step, so give them a longer index array
*/
source_t * make_sky_model(int num_points, int num_gauss,
                             int num_shapes, int num_coeff_per_shape,
                             int num_time_steps) {

  source_t *cropped_src = malloc(sizeof(source_t));

  //Overall COMPONENT stats
  cropped_src->n_comps = num_points*num_gauss*num_shapes;
  cropped_src->n_points= num_points;
  cropped_src->n_gauss = num_gauss;
  cropped_src->n_shapes = num_shapes;
  cropped_src->n_shape_coeffs = num_shapes*num_coeff_per_shape;

  //Make index array for POINT source
  user_precision_t *index_point_array = malloc(num_points*sizeof(user_precision_t));
  make_index_array(index_point_array, num_points, 0);

  double *index_point_array_double = malloc(num_points*sizeof(double));
  make_index_array_double(index_point_array_double, num_points, 0);

  //Populate POINT intrinsic properties
  cropped_src->point_components.ras = index_point_array_double;
  cropped_src->point_components.decs = index_point_array_double;
  cropped_src->point_components.ref_freqs = index_point_array_double;
  cropped_src->point_components.ref_stokesI = index_point_array;
  cropped_src->point_components.ref_stokesQ = index_point_array;
  cropped_src->point_components.ref_stokesU = index_point_array;
  cropped_src->point_components.ref_stokesV = index_point_array;
  cropped_src->point_components.SIs = index_point_array;

  //Make repeating array for POINT primary beam related values
  user_precision_t *repeat_point_array = malloc(num_time_steps*num_points*sizeof(user_precision_t));
  make_repeat_array(repeat_point_array, num_points, num_time_steps, 0);

  double *repeat_point_array_double = malloc(num_time_steps*num_points*sizeof(double));
  make_repeat_array_double(repeat_point_array_double, num_points, num_time_steps, 0);

  cropped_src->point_components.azs = repeat_point_array;
  cropped_src->point_components.zas = repeat_point_array;
  cropped_src->point_components.cos_para_angs = repeat_point_array;
  cropped_src->point_components.sin_para_angs = repeat_point_array;
  cropped_src->point_components.beam_has = repeat_point_array_double;
  cropped_src->point_components.beam_decs = repeat_point_array_double;

  //Repeat process for GAUSSIAN and SHAPELETs

  user_precision_t *index_gauss_array = malloc(num_gauss*sizeof(user_precision_t));
  make_index_array(index_gauss_array, num_gauss, 0);

  double *index_gauss_array_double = malloc(num_gauss*sizeof(double));
  make_index_array_double(index_gauss_array_double, num_gauss, 0);

  user_precision_t *repeat_gauss_array = malloc(num_time_steps*num_gauss*sizeof(user_precision_t));
  make_repeat_array(repeat_gauss_array, num_gauss, num_time_steps, 0);

  double *repeat_gauss_array_double = malloc(num_time_steps*num_gauss*sizeof(double));
  make_repeat_array_double(repeat_gauss_array_double, num_gauss, num_time_steps, 0);

  cropped_src->gauss_components.ras = index_gauss_array_double;
  cropped_src->gauss_components.decs = index_gauss_array_double;
  cropped_src->gauss_components.ref_freqs = index_gauss_array_double;
  cropped_src->gauss_components.ref_stokesI = index_gauss_array;
  cropped_src->gauss_components.ref_stokesQ = index_gauss_array;
  cropped_src->gauss_components.ref_stokesU = index_gauss_array;
  cropped_src->gauss_components.ref_stokesV = index_gauss_array;
  cropped_src->gauss_components.SIs = index_gauss_array;
  cropped_src->gauss_components.majors = index_gauss_array;
  cropped_src->gauss_components.minors = index_gauss_array;
  cropped_src->gauss_components.pas = index_gauss_array;
  cropped_src->gauss_components.azs = repeat_gauss_array;
  cropped_src->gauss_components.zas = repeat_gauss_array;
  cropped_src->gauss_components.cos_para_angs = repeat_gauss_array;
  cropped_src->gauss_components.sin_para_angs = repeat_gauss_array;
  cropped_src->gauss_components.beam_has = repeat_gauss_array_double;
  cropped_src->gauss_components.beam_decs = repeat_gauss_array_double;

  user_precision_t *index_shape_array = malloc(num_shapes*sizeof(user_precision_t));
  make_index_array(index_shape_array, num_shapes, 0);

  double *index_shape_array_double = malloc(num_shapes*sizeof(double));
  make_index_array_double(index_shape_array_double, num_shapes, 0);

  user_precision_t *repeat_shape_array = malloc(num_time_steps*num_shapes*sizeof(user_precision_t));
  make_repeat_array(repeat_shape_array, num_shapes, num_time_steps, 0);

  double *repeat_shape_array_double = malloc(num_time_steps*num_shapes*sizeof(double));
  make_repeat_array_double(repeat_shape_array_double, num_shapes, num_time_steps, 0);

  cropped_src->shape_components.ras = index_shape_array_double;
  cropped_src->shape_components.decs = index_shape_array_double;
  cropped_src->shape_components.ref_freqs = index_shape_array_double;
  cropped_src->shape_components.ref_stokesI = index_shape_array;
  cropped_src->shape_components.ref_stokesQ = index_shape_array;
  cropped_src->shape_components.ref_stokesU = index_shape_array;
  cropped_src->shape_components.ref_stokesV = index_shape_array;
  cropped_src->shape_components.SIs = index_shape_array;
  cropped_src->shape_components.majors = index_shape_array;
  cropped_src->shape_components.minors = index_shape_array;
  cropped_src->shape_components.pas = index_shape_array;
  cropped_src->shape_components.azs = repeat_shape_array;
  cropped_src->shape_components.zas = repeat_shape_array;
  cropped_src->shape_components.cos_para_angs = repeat_shape_array;
  cropped_src->shape_components.sin_para_angs = repeat_shape_array;
  cropped_src->shape_components.beam_has = repeat_shape_array_double;
  cropped_src->shape_components.beam_decs = repeat_shape_array_double;

  //These arrays also contain repeating values
  user_precision_t *repeat_coeffs_array = malloc(num_coeff_per_shape*num_shapes*sizeof(user_precision_t));
  make_repeat_array(repeat_coeffs_array, num_shapes, num_coeff_per_shape, 0);

  cropped_src->shape_components.shape_coeffs = repeat_coeffs_array;
  cropped_src->shape_components.n1s = repeat_coeffs_array;
  cropped_src->shape_components.n2s = repeat_coeffs_array;
  cropped_src->shape_components.param_indexes = repeat_coeffs_array;

  return cropped_src;

}

void free_sky_model(source_t *cropped_src) {

  //These point to the unique arrays we made in `make_sky_model`
  //so freeing just these is enough
  free(cropped_src->point_components.ras);
  free(cropped_src->point_components.ref_stokesI);
  free(cropped_src->point_components.azs);
  free(cropped_src->gauss_components.ras);
  free(cropped_src->gauss_components.ref_stokesI);
  free(cropped_src->gauss_components.azs);
  free(cropped_src->shape_components.ras);
  free(cropped_src->shape_components.ref_stokesI);
  free(cropped_src->shape_components.azs);
  free(cropped_src->shape_components.shape_coeffs);

  free(cropped_src);
}


void check_pointgauss_chunking(int chunk_ind, int comps_per_chunk,
                             int num_time_steps,
                             int * point_accum, int * gauss_accum,
                             source_t *cropped_src,
                             source_t *temp_cropped_src) {

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

    //Some things are a choice between float and double at compilation, and
    //some things are always set to double, so setup two different expectation
    //arrays
    user_precision_t *expec_index_point_array = malloc(expected_n_points*sizeof(user_precision_t));
    make_index_array(expec_index_point_array, expected_n_points, * point_accum);

    user_precision_t *expec_repeat_point_array = malloc(num_time_steps*expected_n_points*sizeof(user_precision_t));
    make_repeat_array(expec_repeat_point_array, expected_n_points,
                      num_time_steps, * point_accum);

    double *expec_index_point_array_double = malloc(expected_n_points*sizeof(double));
    make_index_array_double(expec_index_point_array_double, expected_n_points, * point_accum);

    double *expec_repeat_point_array_double = malloc(num_time_steps*expected_n_points*sizeof(double));
    make_repeat_array_double(expec_repeat_point_array_double, expected_n_points,
                      num_time_steps, * point_accum);


    //Check POINT source params were split correctly
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(expec_index_point_array_double,
                            temp_cropped_src->point_components.ras, expected_n_points);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(expec_index_point_array_double,
                            temp_cropped_src->point_components.decs, expected_n_points);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(expec_index_point_array_double,
                            temp_cropped_src->point_components.ref_freqs, expected_n_points);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_index_point_array,
                            temp_cropped_src->point_components.ref_stokesI, expected_n_points);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_index_point_array,
                            temp_cropped_src->point_components.ref_stokesQ, expected_n_points);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_index_point_array,
                            temp_cropped_src->point_components.ref_stokesU, expected_n_points);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_index_point_array,
                            temp_cropped_src->point_components.ref_stokesV, expected_n_points);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_index_point_array,
                            temp_cropped_src->point_components.SIs, expected_n_points);
    //
    //Check POINT source prinary beam params were split correctly
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_repeat_point_array,
                    temp_cropped_src->point_components.azs, expected_n_points*num_time_steps);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_repeat_point_array,
                    temp_cropped_src->point_components.zas, expected_n_points*num_time_steps);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_repeat_point_array,
          temp_cropped_src->point_components.cos_para_angs, expected_n_points*num_time_steps);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_repeat_point_array,
          temp_cropped_src->point_components.sin_para_angs, expected_n_points*num_time_steps);

    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(expec_repeat_point_array_double,
         temp_cropped_src->point_components.beam_has, expected_n_points*num_time_steps);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(expec_repeat_point_array_double,
         temp_cropped_src->point_components.beam_decs, expected_n_points*num_time_steps);

    free(expec_index_point_array);
    free(expec_repeat_point_array);
    free(expec_index_point_array_double);
    free(expec_repeat_point_array_double);

  }

  if (expected_n_gauss > 0) {
    // printf("Found GAUSS, chunk_ind %d expected_n_gausss %d\n",chunk_ind, expected_n_gauss );
    TEST_ASSERT_EQUAL_INT(expected_n_gauss, temp_cropped_src->n_gauss);

    user_precision_t *expec_index_gauss_array = malloc(expected_n_gauss*sizeof(user_precision_t));
    make_index_array(expec_index_gauss_array, expected_n_gauss, * gauss_accum);

    user_precision_t *expec_repeat_gauss_array = malloc(num_time_steps*expected_n_gauss*sizeof(user_precision_t));
    make_repeat_array(expec_repeat_gauss_array, expected_n_gauss,
                      num_time_steps, * gauss_accum);

    double *expec_index_gauss_array_double = malloc(expected_n_gauss*sizeof(double));
    make_index_array_double(expec_index_gauss_array_double, expected_n_gauss, * gauss_accum);

    double *expec_repeat_gauss_array_double = malloc(num_time_steps*expected_n_gauss*sizeof(double));
    make_repeat_array_double(expec_repeat_gauss_array_double, expected_n_gauss,
                      num_time_steps, * gauss_accum);

    //Check GAUSS source params were split correctly
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(expec_index_gauss_array_double,
                            temp_cropped_src->gauss_components.ras, expected_n_gauss);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(expec_index_gauss_array_double,
                            temp_cropped_src->gauss_components.decs, expected_n_gauss);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(expec_index_gauss_array_double,
                            temp_cropped_src->gauss_components.ref_freqs, expected_n_gauss);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_index_gauss_array,
                                  temp_cropped_src->gauss_components.ref_stokesI, expected_n_gauss);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_index_gauss_array,
                                  temp_cropped_src->gauss_components.ref_stokesQ, expected_n_gauss);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_index_gauss_array,
                                  temp_cropped_src->gauss_components.ref_stokesU, expected_n_gauss);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_index_gauss_array,
                                  temp_cropped_src->gauss_components.ref_stokesV, expected_n_gauss);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_index_gauss_array,
                                  temp_cropped_src->gauss_components.SIs, expected_n_gauss);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_index_gauss_array,
                                  temp_cropped_src->gauss_components.pas, expected_n_gauss);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_index_gauss_array,
                                  temp_cropped_src->gauss_components.majors, expected_n_gauss);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_index_gauss_array,
                                  temp_cropped_src->gauss_components.minors, expected_n_gauss);

    //Check GAUSS source prinary beam params were split correctly
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_repeat_gauss_array,
                    temp_cropped_src->gauss_components.azs, expected_n_gauss*num_time_steps);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_repeat_gauss_array,
                    temp_cropped_src->gauss_components.zas, expected_n_gauss*num_time_steps);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_repeat_gauss_array,
          temp_cropped_src->gauss_components.cos_para_angs, expected_n_gauss*num_time_steps);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_repeat_gauss_array,
          temp_cropped_src->gauss_components.sin_para_angs, expected_n_gauss*num_time_steps);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(expec_repeat_gauss_array_double,
         temp_cropped_src->gauss_components.beam_has, expected_n_gauss*num_time_steps);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(expec_repeat_gauss_array_double,
         temp_cropped_src->gauss_components.beam_decs, expected_n_gauss*num_time_steps);

    free(expec_index_gauss_array);
    free(expec_repeat_gauss_array);
    free(expec_index_gauss_array_double);
    free(expec_repeat_gauss_array_double);
  }
  * point_accum += expected_n_points;
  * gauss_accum += expected_n_gauss;

}




void check_shapelet_chunking(int chunk_ind, int coeffs_per_chunk,
                             int num_time_steps,
                             int num_coeff_per_shape,
                             source_t *cropped_src,
                             source_t *temp_cropped_src){

  //How far through all coeffs are we with this chunk
  int chunk_coeff_ind = chunk_ind*coeffs_per_chunk;
  //How many POINT sources are left after this chunk index
  int coeff_remainder = cropped_src->n_shape_coeffs - chunk_coeff_ind;

  //Things to use in logic below
  int expected_n_coeffs = 0;
  // int chunk_remainder = 0;

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
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(cropped_src->shape_components.ras,
                            temp_cropped_src->shape_components.ras, cropped_src->n_shapes);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(cropped_src->shape_components.decs,
                            temp_cropped_src->shape_components.decs, cropped_src->n_shapes);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(cropped_src->shape_components.ref_freqs,
                            temp_cropped_src->shape_components.ref_freqs, cropped_src->n_shapes);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(cropped_src->shape_components.ref_stokesI,
                            temp_cropped_src->shape_components.ref_stokesI, cropped_src->n_shapes);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(cropped_src->shape_components.ref_stokesQ,
                            temp_cropped_src->shape_components.ref_stokesQ, cropped_src->n_shapes);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(cropped_src->shape_components.ref_stokesU,
                            temp_cropped_src->shape_components.ref_stokesU, cropped_src->n_shapes);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(cropped_src->shape_components.ref_stokesV,
                            temp_cropped_src->shape_components.ref_stokesV, cropped_src->n_shapes);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(cropped_src->shape_components.SIs,
                            temp_cropped_src->shape_components.SIs, cropped_src->n_shapes);
    //
    //As we only split over basis function coeff info, all of these arrrays
    //should just be pointer copies
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(cropped_src->shape_components.azs,
                    temp_cropped_src->shape_components.azs, cropped_src->n_shapes*num_time_steps);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(cropped_src->shape_components.zas,
                    temp_cropped_src->shape_components.zas, cropped_src->n_shapes*num_time_steps);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(cropped_src->shape_components.cos_para_angs,
          temp_cropped_src->shape_components.cos_para_angs, cropped_src->n_shapes*num_time_steps);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(cropped_src->shape_components.sin_para_angs,
          temp_cropped_src->shape_components.sin_para_angs, cropped_src->n_shapes*num_time_steps);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(cropped_src->shape_components.beam_has,
          temp_cropped_src->shape_components.beam_has, cropped_src->n_shapes*num_time_steps);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(cropped_src->shape_components.beam_decs,
          temp_cropped_src->shape_components.beam_decs, cropped_src->n_shapes*num_time_steps);

    //THESE are the arrays that should actually be split up
    //With the way I've set up the sky model creation, the coeff splitting
    //should yield arrays which are the index of cropped_src integer divided
    //by number of coeffs per shapelet

    user_precision_t *expec_repeat_shape_array = malloc(expected_n_coeffs*sizeof(user_precision_t));

    int new_ind = 0;
    for (int orig_index = chunk_coeff_ind; orig_index < chunk_coeff_ind + expected_n_coeffs; orig_index++) {
      expec_repeat_shape_array[new_ind] = (user_precision_t)(orig_index / num_coeff_per_shape);
      new_ind ++;
    }

    TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_repeat_shape_array,
                    temp_cropped_src->shape_components.shape_coeffs, expected_n_coeffs);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_repeat_shape_array,
                    temp_cropped_src->shape_components.n1s, expected_n_coeffs);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_repeat_shape_array,
                    temp_cropped_src->shape_components.n2s, expected_n_coeffs);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_repeat_shape_array,
                    temp_cropped_src->shape_components.param_indexes, expected_n_coeffs);

    free(expec_repeat_shape_array);

  }
}
