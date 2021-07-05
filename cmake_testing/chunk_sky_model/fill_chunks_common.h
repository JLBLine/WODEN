#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "woden_struct_defs.h"


//Populate an array with index values. Add an offset of `offset`
void make_index_array(float *index_array, int array_length, int offset);

/*
Given a maximum value `num_value`, populate an array from 0 to `num_value`-1,
repeating each value `num_repeat times`, also add an `offset` value
i.e. if num_value = 2, num_repeat_times = 3, offset=4 makes an array
repeat_array[] = {4, 4, 4, 5, 5, 5}

*/
void make_repeat_array(float *repeat_array, int num_value, int num_repeat,
                       int offset);


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
                             int num_time_steps);

//Frees the created sky model
void free_sky_model(catsource_t *cropped_src);

/*
Check the point/gaussian chunking has worked for a particular chunk
*/
void check_pointgauss_chunking(int chunk_ind, int comps_per_chunk,
                             int num_time_steps,
                             int * point_accum, int * gauss_accum,
                             catsource_t *cropped_src,
                             catsource_t *temp_cropped_src);

/*
Check the shapelet chunking has worked for a particular chunk
*/
void check_shapelet_chunking(int chunk_ind, int coeffs_per_chunk,
                             int num_time_steps,
                             int num_coeff_per_shape,
                             catsource_t *cropped_src,
                             catsource_t *temp_cropped_src);
