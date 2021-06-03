#include <unity.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "chunk_sky_model.h"
#include "woden_struct_defs.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

#define UNITY_INCLUDE_FLOAT

/*
The GPU will never have enough memory to simulate millions of COMPONENTs at
once, so we need to split the sky model up before running the simulation.
`fill_chunk_src_*` functions do just that, given the sky model above the
horizon (`cropped_src`) and a chunking index, they return a portion of the
sky model containing specific COMPONENT types.

As SHAPELETs take more memory, we separate them out from POINT and GAUSSIANs,
so have two chunking functions.
*/

//Populate an array with index values
void make_index_array(float *index_array, int array_length) {
  for (int index = 0; index < array_length; index++) {
    index_array[index] = index;
  }
}

/*
Given a maximum value `num_value`, populate an array from 0 to `num_value`-1,
repeating each value `num_repeat times`
i.e. if num_value = 2, num_repeat_times = 3, make array
repeat_array[] = {0, 0, 0, 1, 1, 1}

*/
void make_repeat_array(float *repeat_array, int num_value, int num_repeat) {
  int index = 0;
  for (int value = 0; value < num_value; value++) {
    for (int repeat = 0; repeat < num_repeat; repeat++) {
      repeat_array[index] = value;
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
  make_index_array(index_point_array, num_points);

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
  make_repeat_array(repeat_point_array, num_points, num_time_steps);
  cropped_src->point_azs = repeat_point_array;
  cropped_src->point_zas = repeat_point_array;
  cropped_src->cos_point_para_angs = repeat_point_array;
  cropped_src->sin_point_para_angs = repeat_point_array;
  cropped_src->point_gaussbeam_has = repeat_point_array;
  cropped_src->point_gaussbeam_decs = repeat_point_array;

  //Repeat process for GAUSSIAN and SHAPELETs

  float *index_gauss_array = malloc(num_gauss*sizeof(float));
  make_index_array(index_gauss_array, num_gauss);

  float *repeat_gauss_array = malloc(num_time_steps*num_gauss*sizeof(float));
  make_repeat_array(repeat_gauss_array, num_gauss, num_time_steps);

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
  make_index_array(index_shape_array, num_shapes);

  float *repeat_shape_array = malloc(num_time_steps*num_shapes*sizeof(float));
  make_repeat_array(repeat_shape_array, num_shapes, num_time_steps);

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
  make_repeat_array(repeat_coeffs_array, num_shapes, num_coeff_per_shape);

  cropped_src->shape_coeffs = repeat_coeffs_array;
  cropped_src->shape_n1s = repeat_coeffs_array;
  cropped_src->shape_n2s = repeat_coeffs_array;
  cropped_src->shape_param_indexes = repeat_coeffs_array;

  free(index_point_array);
  free(repeat_point_array);
  free(index_gauss_array);
  free(repeat_gauss_array);
  free(index_shape_array);
  free(repeat_shape_array);

  return cropped_src;

}

void test_fill_chunk_src_with_pointgauss(int chunking_size,
                                         int num_points, int num_gauss,
                                         int num_time_steps) {

  //Not testing SHAPELETs, set these to zero
  int num_shapes = 0;
  int num_coeff_per_shape = 0;
  catsource_t *cropped_src = make_sky_model(num_points, num_gauss, num_shapes,
                                            num_coeff_per_shape, num_time_steps);

  // printf("%.1f %.1f\n",cropped_src->cos_point_para_angs[75], cropped_src->point_ras[75] );

  //Set some dummy baseline, frequency simulation settings
  int num_baselines = 5;
  int num_freqs = 2;
  int num_visis = num_baselines * num_freqs * num_time_steps;

  //
  int num_comps_to_chunk = cropped_src->n_points + cropped_src->n_gauss;
  int comps_per_chunk = (int)floorf((float)chunking_size / (float)(num_baselines * num_freqs * num_time_steps));
  int num_chunks = (int)ceilf((float)num_comps_to_chunk / (float)comps_per_chunk);

  printf("Num visis %d\n",num_visis );
  printf("Chunking size %d\n",chunking_size );
  printf("Comps per chunk %d\n", comps_per_chunk );
  printf("Num chunks %d\n", num_chunks );
  printf("Is num_visis*num_comps %d < num_visis*comps_per_chunk*num_chunks %d\n",num_visis*num_comps_to_chunk, num_visis*comps_per_chunk*num_chunks);


  catsource_t *temp_cropped_src = malloc(sizeof(catsource_t));

  for (int chunk_ind = 0; chunk_ind < num_chunks; chunk_ind++) {

    fill_chunk_src_with_pointgauss(temp_cropped_src, cropped_src, chunk_ind,
                                   comps_per_chunk, num_time_steps);

    printf("Chunk %d has %d points\n",chunk_ind, temp_cropped_src->n_points );


  }

}


void test_fill_chunk_src_with_pointgauss_Point100_Gauss000_Chunk1000_Time004(void){
  int chunking_size = 1000;
  int num_points = 100;
  int num_gauss = 0;
  int num_time_steps = 4;

  test_fill_chunk_src_with_pointgauss(chunking_size, num_points, num_gauss,
                                      num_time_steps);
}



//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_fill_chunk_src_with_pointgauss_Point100_Gauss000_Chunk1000_Time004);

    return UNITY_END();
}
