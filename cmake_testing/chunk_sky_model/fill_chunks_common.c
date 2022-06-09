#include <stdlib.h>
#include <math.h>
#include <unity.h>

#include "constants.h"
#include "woden_struct_defs.h"
#include "fill_chunks_common.h"

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
void make_index_array_user(user_precision_t *index_array, int array_length, int offset) {
  for (int index = 0; index < array_length; index++) {
    index_array[index] = (user_precision_t)(index + offset);
    // printf("MAKING AN INDEX ARRAY %f\n", index_array[index] );
  }
}

//Populate an array with index values. Add an offset of `offset`
void make_index_array_int(int *index_array, int array_length, int offset) {
  for (int index = 0; index < array_length; index++) {
    index_array[index] = (index + offset);
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

void make_repeat_array_int(int *repeat_array, int num_value, int num_repeat,
                       int offset) {
  int index = 0;
  for (int value = 0; value < num_value; value++) {
    for (int repeat = 0; repeat < num_repeat; repeat++) {
      repeat_array[index] = offset + value;
      // printf("MAKIG REPEAT ARRAY %d\n",repeat_array[index] );
      index ++;
    }

  }
}

void fill_components(components_t *components, int num_comps,
                     int num_coeff_per_shape,
                     int num_list_values, int num_time_steps,
                     e_component_type comptype){
    //Make index arrays
    user_precision_t *index_array_user = malloc(num_comps*sizeof(user_precision_t));
    make_index_array_user(index_array_user, num_comps, 0);

    double *index_array_double = malloc(num_comps*sizeof(double));
    make_index_array_double(index_array_double, num_comps, 0);

    user_precision_t *triple_index_array_user = malloc(NUM_FLUX_TYPES*num_comps*sizeof(user_precision_t));
    make_index_array_user(triple_index_array_user, NUM_FLUX_TYPES*num_comps, 0);

    double *triple_index_array_double = malloc(NUM_FLUX_TYPES*num_comps*sizeof(double));
    make_index_array_double(triple_index_array_double, NUM_FLUX_TYPES*num_comps, 0);

    int *index_array_int = malloc(num_comps*sizeof(int));
    make_index_array_int(index_array_int, num_comps, 0);

    //Populate POINT intrinsic properties
    //Need enough to cover all the different flux types so use triple arrays
    components->ras = triple_index_array_double;
    components->decs = triple_index_array_double;

    //Power law type stuff
    components->power_ref_freqs = index_array_double;
    components->power_ref_stokesI = index_array_user;
    components->power_ref_stokesQ = index_array_user;
    components->power_ref_stokesU = index_array_user;
    components->power_ref_stokesV = index_array_user;
    components->power_SIs = index_array_user;
    components->power_comp_inds = index_array_int;

    //Curved power law type stuff
    components->curve_ref_freqs = index_array_double;
    components->curve_ref_stokesI = index_array_user;
    components->curve_ref_stokesQ = index_array_user;
    components->curve_ref_stokesU = index_array_user;
    components->curve_ref_stokesV = index_array_user;
    components->curve_SIs = index_array_user;

    //Set things up so the curved power law components are after the power law
    components->curve_comp_inds = malloc(num_comps*sizeof(int));
    for (int i = 0; i < num_comps; i++) {
      components->curve_comp_inds[i] = num_comps + i;
    }

    //List freq things
    user_precision_t *repeat_array_list = malloc(num_list_values*num_comps*sizeof(user_precision_t));
    make_repeat_array(repeat_array_list, num_comps, num_list_values, 0);

    double *repeat_array_list_double = malloc(num_list_values*num_comps*sizeof(double));
    make_repeat_array_double(repeat_array_list_double, num_comps, num_list_values, 0);

    components->list_freqs = repeat_array_list_double;
    components->list_stokesI = repeat_array_list;
    components->list_stokesQ = repeat_array_list;
    components->list_stokesU = repeat_array_list;
    components->list_stokesV = repeat_array_list;
    components->num_list_values = malloc(num_comps*sizeof(int));
    components->list_start_indexes = malloc(num_comps*sizeof(int));
    components->list_comp_inds = malloc(num_comps*sizeof(int));

    //Set things up so the curved power law components are after the power law

    for (int i = 0; i < num_comps; i++) {
      components->num_list_values[i] = num_list_values;
      components->list_start_indexes[i] = i*num_list_values;
      components->list_comp_inds[i] = 2*num_comps + i;
    }

    //Make repeating array for primary beam related values
    user_precision_t *repeat_array = malloc(NUM_FLUX_TYPES*num_time_steps*num_comps*sizeof(user_precision_t));
    make_repeat_array(repeat_array, NUM_FLUX_TYPES*num_comps, num_time_steps, 0);

    double *repeat_array_double = malloc(NUM_FLUX_TYPES*num_time_steps*num_comps*sizeof(double));
    make_repeat_array_double(repeat_array_double, NUM_FLUX_TYPES*num_comps, num_time_steps, 0);

    components->azs = repeat_array;
    components->zas = repeat_array;
    components->beam_has = repeat_array_double;
    components->beam_decs = repeat_array_double;

    if (comptype == GAUSSIAN || comptype == SHAPELET) {
      components->majors = triple_index_array_user;
      components->minors = triple_index_array_user;
      components->pas = triple_index_array_user;
    }

    if (comptype == SHAPELET) {

      // printf("LENGTH OF SHAPE COEFFS %d %d\n",num_coeff_per_shape,num_comps );

        user_precision_t *repeat_coeffs_array = malloc(NUM_FLUX_TYPES*num_coeff_per_shape*num_comps*sizeof(user_precision_t));
        make_repeat_array(repeat_coeffs_array, NUM_FLUX_TYPES*num_comps, num_coeff_per_shape, 0);

        components->shape_coeffs = repeat_coeffs_array;
        components->n1s = repeat_coeffs_array;
        components->n2s = repeat_coeffs_array;
        components->param_indexes = repeat_coeffs_array;

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
                          int num_list_values,
                          int num_time_steps) {

  source_t *cropped_src = malloc(sizeof(source_t));

  //Overall COMPONENT stats - multiply by three as three flux model behaviours
  cropped_src->n_comps = NUM_FLUX_TYPES*(num_points + num_gauss + num_shapes);

  cropped_src->n_points = NUM_FLUX_TYPES*num_points;
  cropped_src->n_point_powers = num_points;
  cropped_src->n_point_curves = num_points;
  cropped_src->n_point_lists = num_points;

  cropped_src->n_gauss = NUM_FLUX_TYPES*num_gauss;
  cropped_src->n_gauss_powers = num_gauss;
  cropped_src->n_gauss_curves = num_gauss;
  cropped_src->n_gauss_lists = num_gauss;

  cropped_src->n_shapes = NUM_FLUX_TYPES*num_shapes;
  cropped_src->n_shape_powers = num_shapes;
  cropped_src->n_shape_curves = num_shapes;
  cropped_src->n_shape_lists = num_shapes;
  cropped_src->n_shape_coeffs = NUM_FLUX_TYPES*num_shapes*num_coeff_per_shape;

  fill_components(&cropped_src->point_components, num_points, 0, num_list_values,
                    num_time_steps, POINT);

  fill_components(&cropped_src->gauss_components, num_gauss, 0, num_list_values,
                    num_time_steps, GAUSSIAN);

  fill_components(&cropped_src->shape_components, num_shapes, num_coeff_per_shape,
                    num_list_values,
                    num_time_steps, SHAPELET);

  return cropped_src;

}

void free_sky_model(source_t *cropped_src) {

  //These point to the unique arrays we made in `make_sky_model`
  //so freeing just these is enough
  free(cropped_src->point_components.ras);
  free(cropped_src->point_components.power_ref_stokesI);
  free(cropped_src->point_components.azs);
  free(cropped_src->gauss_components.ras);
  free(cropped_src->gauss_components.power_ref_stokesI);
  free(cropped_src->gauss_components.azs);
  free(cropped_src->shape_components.ras);
  free(cropped_src->shape_components.power_ref_stokesI);
  free(cropped_src->shape_components.azs);
  free(cropped_src->shape_components.shape_coeffs);

  free(cropped_src);
}


void check_pointgauss_chunking(int chunk_ind, int comps_per_chunk, int num_chunks,
                             int num_list_values,
                             int * n_comp_accum, int * n_power_accum,
                             int * n_curve_accum, int * n_list_accum,
                             e_component_type comptype,
                             source_t *cropped_src,
                             source_t *temp_cropped_src) {



  int total_n_powers = 0;
  int total_n_curves = 0;
  int total_n_lists = 0;
  int total_n_comps = 0;

  int found_n_powers = 0;
  int found_n_curves = 0;
  int found_n_lists = 0;

  components_t *components = NULL;
  components_t *temp_components = NULL;

  if (comptype == POINT) {
    total_n_powers = cropped_src->n_point_powers;
    total_n_curves = cropped_src->n_point_curves;
    total_n_lists = cropped_src->n_point_lists;
    total_n_comps = cropped_src->n_points;

    found_n_powers = temp_cropped_src->n_point_powers;
    found_n_curves = temp_cropped_src->n_point_curves;
    found_n_lists = temp_cropped_src->n_point_lists;

    components = &cropped_src->point_components;
    temp_components = &temp_cropped_src->point_components;
  }

  else if (comptype == GAUSSIAN) {
    total_n_powers = cropped_src->n_gauss_powers;
    total_n_curves = cropped_src->n_gauss_curves;
    total_n_lists = cropped_src->n_gauss_lists;
    total_n_comps = cropped_src->n_gauss;

    found_n_powers = temp_cropped_src->n_gauss_powers;
    found_n_curves = temp_cropped_src->n_gauss_curves;
    found_n_lists = temp_cropped_src->n_gauss_lists;

    components = &cropped_src->gauss_components;
    temp_components = &temp_cropped_src->gauss_components;
  }

  int expec_n_comps = 0;
  int expec_n_powers = 0;
  int expec_n_curves = 0;
  int expec_n_lists = 0;


  // //How far through all components are we with this chunk
  int chunk_comp_ind = chunk_ind*comps_per_chunk;
  // //How many components are left after this chunk index
  // int point_remainder = total_n_comps - chunk_comp_ind;

  // printf("%d %d %d\n",chunk_ind, chunk_comp_ind, point_remainder );



  if (total_n_powers >= chunk_comp_ind + comps_per_chunk) {
    expec_n_powers = comps_per_chunk;
  }
  else {
    //if chunk_comp_ind >= total_n_powers, we should already have chunked
    //all the power-law components
    if (chunk_comp_ind >= total_n_powers ) {
      expec_n_powers = 0;
    } else {
      expec_n_powers = total_n_powers - chunk_comp_ind;
    }

    //update the index where this comp begins with however many power-law
    //components we found

    chunk_comp_ind += expec_n_powers;

    // printf("HERE %d %d %d\n", expec_n_powers, total_n_curves, comps_per_chunk + chunk_comp_ind - total_n_powers);

    //If there are enough curved power laws to fill the rest of this chunk
    //take off how ever many power-law components we already have
    if (total_n_curves >= comps_per_chunk + chunk_comp_ind - total_n_powers ) {
      expec_n_curves = comps_per_chunk - expec_n_powers;
    }

    else {
      //There are some curve components left, but not enough to fill the chunk
      if (total_n_curves > chunk_comp_ind - total_n_powers && total_n_curves != 0) {

        expec_n_curves = total_n_curves - chunk_comp_ind + total_n_powers;
      } else {
        expec_n_curves = 0;
      }

      //Update the lower index with how many curved power law sources we added

      chunk_comp_ind += expec_n_curves;

      // printf("HUH %d %d\n",total_n_lists, comps_per_chunk + chunk_comp_ind - total_n_powers - total_n_curves );

      if (total_n_lists >= comps_per_chunk + chunk_comp_ind - total_n_powers - total_n_curves ) {
        expec_n_lists = comps_per_chunk - expec_n_powers - expec_n_curves;
      }
      //There are some curve components left, but not enough to fill the chunk
      else if (total_n_lists > chunk_comp_ind - total_n_powers - total_n_curves && total_n_lists != 0) {

        expec_n_lists = total_n_lists - chunk_comp_ind + total_n_powers + total_n_curves;
      } else {
        expec_n_lists = 0;
      }
    }
  }
  expec_n_comps = expec_n_powers + expec_n_curves + expec_n_lists;

  if (comptype == POINT) {
    TEST_ASSERT_EQUAL_INT(expec_n_comps, temp_cropped_src->n_points);
    TEST_ASSERT_EQUAL_INT(expec_n_powers, temp_cropped_src->n_point_powers);
    TEST_ASSERT_EQUAL_INT(expec_n_curves, temp_cropped_src->n_point_curves);
    TEST_ASSERT_EQUAL_INT(expec_n_lists, temp_cropped_src->n_point_lists);
    // printf("found/predic n_comps %d %d\n", temp_cropped_src->n_points, expec_n_comps);
    // printf("found/predic n_powers %d %d\n",
    //         temp_cropped_src->n_point_powers, expec_n_powers);
    // printf("found/predic n_curves %d %d\n",
    //         temp_cropped_src->n_point_curves, expec_n_curves);
    // printf("found/predic n_lists %d %d\n",
    //         temp_cropped_src->n_point_lists, expec_n_lists);
  }

  else if (comptype == GAUSSIAN) {
    TEST_ASSERT_EQUAL_INT(expec_n_comps, temp_cropped_src->n_gauss);
    TEST_ASSERT_EQUAL_INT(expec_n_powers, temp_cropped_src->n_gauss_powers);
    TEST_ASSERT_EQUAL_INT(expec_n_curves, temp_cropped_src->n_gauss_curves);
    TEST_ASSERT_EQUAL_INT(expec_n_lists, temp_cropped_src->n_gauss_lists);
    // printf("found/predic n_comps %d %d\n", temp_cropped_src->n_gauss, expec_n_comps);
    // printf("found/predic n_powers %d %d\n",
    //         temp_cropped_src->n_gauss_powers, expec_n_powers);
    // printf("found/predic n_curves %d %d\n",
    //         temp_cropped_src->n_gauss_curves, expec_n_curves);
    // printf("found/predic n_lists %d %d\n",
    //         temp_cropped_src->n_gauss_lists, expec_n_lists);
  }

  //Pointers should have just been copied for these arrays
  TEST_ASSERT_EQUAL_DOUBLE_ARRAY(components->ras, temp_components->ras,
                            total_n_comps);
  TEST_ASSERT_EQUAL_DOUBLE_ARRAY(components->decs, temp_components->decs,
                            total_n_comps);

  TEST_ASSERT_EQUAL_USER_ARRAY(components->azs, temp_components->azs, total_n_comps);
  TEST_ASSERT_EQUAL_USER_ARRAY(components->zas, temp_components->zas, total_n_comps);
  TEST_ASSERT_EQUAL_DOUBLE_ARRAY(components->beam_has,
                                 temp_components->beam_has, total_n_comps);
  TEST_ASSERT_EQUAL_DOUBLE_ARRAY(components->beam_decs,
                                 temp_components->beam_decs, total_n_comps);

  TEST_ASSERT_EQUAL_DOUBLE_ARRAY(components->list_freqs,
                                 temp_components->list_freqs, total_n_lists);
  TEST_ASSERT_EQUAL_USER_ARRAY(components->list_stokesI,
                                 temp_components->list_stokesI, total_n_lists);
  TEST_ASSERT_EQUAL_USER_ARRAY(components->list_stokesQ,
                                 temp_components->list_stokesQ, total_n_lists);
  TEST_ASSERT_EQUAL_USER_ARRAY(components->list_stokesU,
                                 temp_components->list_stokesU, total_n_lists);
  TEST_ASSERT_EQUAL_USER_ARRAY(components->list_stokesV,
                                 temp_components->list_stokesV, total_n_lists);

  if (comptype == GAUSSIAN) {
    TEST_ASSERT_EQUAL_USER_ARRAY(components->majors, temp_components->majors,
                                                                total_n_comps);
    TEST_ASSERT_EQUAL_USER_ARRAY(components->minors, temp_components->minors,
                                                                total_n_comps);
    TEST_ASSERT_EQUAL_USER_ARRAY(components->pas, temp_components->pas,
                                                                total_n_comps);
  }

  //These things have had pointer arithmatic added to make copying them to
  //the GPU easier later on

  //This should test if we should have POWER_LAW components, or if we've
  //found some when we shouldn't
  if (found_n_powers > 0 || expec_n_powers > 0) {
    // printf("WTF %d %d %d\n", found_n_powers, expec_n_powers, * n_power_accum );
    user_precision_t *expec_power_values_user = malloc(expec_n_powers*sizeof(user_precision_t));
    make_index_array_user(expec_power_values_user, expec_n_powers, * n_power_accum);

    double *expec_power_values_double = malloc(expec_n_powers*sizeof(double));
    make_index_array_double(expec_power_values_double, expec_n_powers, * n_power_accum);

    // for (int i = 0; i < expec_n_powers; i++) {
      // printf("LOOP o truth %.1f %.1f\n", expec_power_values_user[i], temp_components->power_ref_stokesI[i] );
      // printf("LOOP o truth %.1f \n", expec_power_values_user[i] );
    // }

    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(expec_power_values_double,
                              temp_components->power_ref_freqs, expec_n_powers);

    TEST_ASSERT_EQUAL_USER_ARRAY(expec_power_values_user,
                              temp_components->power_ref_stokesI, expec_n_powers);
    TEST_ASSERT_EQUAL_USER_ARRAY(expec_power_values_user,
                              temp_components->power_ref_stokesQ, expec_n_powers);
    TEST_ASSERT_EQUAL_USER_ARRAY(expec_power_values_user,
                              temp_components->power_ref_stokesU, expec_n_powers);
    TEST_ASSERT_EQUAL_USER_ARRAY(expec_power_values_user,
                              temp_components->power_ref_stokesV, expec_n_powers);
    TEST_ASSERT_EQUAL_USER_ARRAY(expec_power_values_user,
                              temp_components->power_SIs, expec_n_powers);

    free(expec_power_values_user);
    free(expec_power_values_double);

  }

  //This should test if we should have CURVED_POWER_LAW components, or if we've
  //found some when we shouldn't
  if (found_n_curves > 0 || expec_n_curves > 0) {
    // printf("WTF2 %d %d %d\n", found_n_curves, expec_n_curves, * n_curve_accum );

    user_precision_t *expec_curve_values_user = malloc(expec_n_curves*sizeof(user_precision_t));
    make_index_array_user(expec_curve_values_user, expec_n_curves, * n_curve_accum);

    double *expec_curve_values_double = malloc(expec_n_curves*sizeof(double));
    make_index_array_double(expec_curve_values_double, expec_n_curves, * n_curve_accum);

    int *expec_curve_comp_inds = malloc(expec_n_curves*sizeof(int));
    make_index_array_int(expec_curve_comp_inds, expec_n_curves,
                            total_n_powers + * n_curve_accum);

    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(expec_curve_values_double,
                              temp_components->curve_ref_freqs, expec_n_curves);

    TEST_ASSERT_EQUAL_USER_ARRAY(expec_curve_values_user,
                              temp_components->curve_ref_stokesI, expec_n_curves);
    TEST_ASSERT_EQUAL_USER_ARRAY(expec_curve_values_user,
                              temp_components->curve_ref_stokesQ, expec_n_curves);
    TEST_ASSERT_EQUAL_USER_ARRAY(expec_curve_values_user,
                              temp_components->curve_ref_stokesU, expec_n_curves);
    TEST_ASSERT_EQUAL_USER_ARRAY(expec_curve_values_user,
                              temp_components->curve_ref_stokesV, expec_n_curves);
    TEST_ASSERT_EQUAL_USER_ARRAY(expec_curve_values_user,
                              temp_components->curve_SIs, expec_n_curves);

    TEST_ASSERT_EQUAL_INT_ARRAY(expec_curve_comp_inds, temp_components->curve_comp_inds,
                                expec_n_curves);

    free(expec_curve_values_user);
    free(expec_curve_values_double);
    free(expec_curve_comp_inds);

  }

    //This should test if we should have listD_POWER_LAW components, or if we've
    //found some when we shouldn't
    if (found_n_lists > 0 || expec_n_lists > 0) {
      // printf("WTF3 %d %d %d\n", found_n_lists, expec_n_lists, * n_list_accum );

      int *expec_list_comp_inds = malloc(expec_n_lists*sizeof(int));
      make_index_array_int(expec_list_comp_inds, expec_n_lists,
                              total_n_powers + total_n_curves + * n_list_accum);

      int *expec_num_list_values = malloc(expec_n_lists*sizeof(int));
      make_repeat_array_int(expec_num_list_values, 1, expec_n_lists, num_list_values);

      int *list_start_indexes = malloc(expec_n_lists*sizeof(int));

      for (int expec = 0; expec < expec_n_lists; expec++) {
        list_start_indexes[expec] = (expec + * n_list_accum)*num_list_values;
        // printf("ROIGHT %d\n",list_start_indexes[expec] );
      }

      TEST_ASSERT_EQUAL_INT_ARRAY(expec_list_comp_inds, temp_components->list_comp_inds,
                                  expec_n_lists);
      TEST_ASSERT_EQUAL_INT_ARRAY(expec_num_list_values, temp_components->num_list_values,
                                  expec_n_lists);
      TEST_ASSERT_EQUAL_INT_ARRAY(list_start_indexes, temp_components->list_start_indexes,
                                  expec_n_lists);

      free(expec_list_comp_inds);
      free(expec_num_list_values);
      free(list_start_indexes);

    }


  * n_comp_accum += expec_n_comps;
  * n_power_accum += expec_n_powers;
  * n_curve_accum += expec_n_curves;
  * n_list_accum += expec_n_lists;
  // printf("------------------------------------------------------\n");
}




void check_shapelet_chunking(int chunk_ind, int coeffs_per_chunk,
                             int num_time_steps,
                             int num_coeff_per_shape,
                             source_t *cropped_src,
                             source_t *temp_cropped_src){

  //How far through all coeffs are we with this chunk
  int chunk_coeff_ind = chunk_ind*coeffs_per_chunk;
  //How many SHAPELET coeffs are left after this chunk index
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

    //As we only split over basis function coeff info, all of these arrrays
    //should just be pointer copies

    TEST_ASSERT_EQUAL_INT(expected_n_coeffs, temp_cropped_src->n_shape_coeffs);
    TEST_ASSERT_EQUAL_INT(cropped_src->n_shapes, temp_cropped_src->n_shapes);
    //
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(cropped_src->shape_components.ras,
                            temp_cropped_src->shape_components.ras, cropped_src->n_shapes);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(cropped_src->shape_components.decs,
                            temp_cropped_src->shape_components.decs, cropped_src->n_shapes);
    //Power law flux stuff
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(cropped_src->shape_components.power_ref_freqs,
        temp_cropped_src->shape_components.power_ref_freqs, cropped_src->n_shape_powers);
    TEST_ASSERT_EQUAL_USER_ARRAY(cropped_src->shape_components.power_ref_stokesI,
        temp_cropped_src->shape_components.power_ref_stokesI, cropped_src->n_shape_powers);
    TEST_ASSERT_EQUAL_USER_ARRAY(cropped_src->shape_components.power_ref_stokesQ,
        temp_cropped_src->shape_components.power_ref_stokesQ, cropped_src->n_shape_powers);
    TEST_ASSERT_EQUAL_USER_ARRAY(cropped_src->shape_components.power_ref_stokesU,
        temp_cropped_src->shape_components.power_ref_stokesU, cropped_src->n_shape_powers);
    TEST_ASSERT_EQUAL_USER_ARRAY(cropped_src->shape_components.power_ref_stokesV,
        temp_cropped_src->shape_components.power_ref_stokesV, cropped_src->n_shape_powers);
    TEST_ASSERT_EQUAL_USER_ARRAY(cropped_src->shape_components.power_SIs,
        temp_cropped_src->shape_components.power_SIs, cropped_src->n_shape_powers);

    //Curved power law flux stuff
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(cropped_src->shape_components.curve_ref_freqs,
        temp_cropped_src->shape_components.curve_ref_freqs, cropped_src->n_shape_curves);
    TEST_ASSERT_EQUAL_USER_ARRAY(cropped_src->shape_components.curve_ref_stokesI,
        temp_cropped_src->shape_components.curve_ref_stokesI, cropped_src->n_shape_curves);
    TEST_ASSERT_EQUAL_USER_ARRAY(cropped_src->shape_components.curve_ref_stokesQ,
        temp_cropped_src->shape_components.curve_ref_stokesQ, cropped_src->n_shape_curves);
    TEST_ASSERT_EQUAL_USER_ARRAY(cropped_src->shape_components.curve_ref_stokesU,
        temp_cropped_src->shape_components.curve_ref_stokesU, cropped_src->n_shape_curves);
    TEST_ASSERT_EQUAL_USER_ARRAY(cropped_src->shape_components.curve_ref_stokesV,
        temp_cropped_src->shape_components.curve_ref_stokesV, cropped_src->n_shape_curves);
    TEST_ASSERT_EQUAL_USER_ARRAY(cropped_src->shape_components.curve_SIs,
        temp_cropped_src->shape_components.curve_SIs, cropped_src->n_shape_curves);

    //Curved power law flux stuff
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(cropped_src->shape_components.list_freqs,
        temp_cropped_src->shape_components.list_freqs, cropped_src->n_shape_lists);
    TEST_ASSERT_EQUAL_USER_ARRAY(cropped_src->shape_components.list_stokesI,
        temp_cropped_src->shape_components.list_stokesI, cropped_src->n_shape_lists);
    TEST_ASSERT_EQUAL_USER_ARRAY(cropped_src->shape_components.list_stokesQ,
        temp_cropped_src->shape_components.list_stokesQ, cropped_src->n_shape_lists);
    TEST_ASSERT_EQUAL_USER_ARRAY(cropped_src->shape_components.list_stokesU,
        temp_cropped_src->shape_components.list_stokesU, cropped_src->n_shape_lists);
    TEST_ASSERT_EQUAL_USER_ARRAY(cropped_src->shape_components.list_stokesV,
        temp_cropped_src->shape_components.list_stokesV, cropped_src->n_shape_lists);


    TEST_ASSERT_EQUAL_USER_ARRAY(cropped_src->shape_components.azs,
                    temp_cropped_src->shape_components.azs, cropped_src->n_shapes*num_time_steps);
    TEST_ASSERT_EQUAL_USER_ARRAY(cropped_src->shape_components.zas,
                    temp_cropped_src->shape_components.zas, cropped_src->n_shapes*num_time_steps);
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
      // printf("%d %.1f\n",chunk_coeff_ind, expec_repeat_shape_array[new_ind]);
      new_ind ++;
    }

    // for (int i = 0; i < expected_n_coeffs; i++) {
    //   printf("F U SHAPE %d %.1f %.1f\n",chunk_ind, expec_repeat_shape_array[i],
    //                   temp_cropped_src->shape_components.shape_coeffs[i] );
    // }

    // printf("Expec %d\n",expected_n_coeffs );

    TEST_ASSERT_EQUAL_USER_ARRAY(expec_repeat_shape_array,
                    temp_cropped_src->shape_components.shape_coeffs, expected_n_coeffs);
    TEST_ASSERT_EQUAL_USER_ARRAY(expec_repeat_shape_array,
                    temp_cropped_src->shape_components.n1s, expected_n_coeffs);
    TEST_ASSERT_EQUAL_USER_ARRAY(expec_repeat_shape_array,
                    temp_cropped_src->shape_components.n2s, expected_n_coeffs);
    TEST_ASSERT_EQUAL_USER_ARRAY(expec_repeat_shape_array,
                    temp_cropped_src->shape_components.param_indexes, expected_n_coeffs);

    free(expec_repeat_shape_array);

  }
}
