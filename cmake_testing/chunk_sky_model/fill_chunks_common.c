#include <stdlib.h>
#include <math.h>

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
