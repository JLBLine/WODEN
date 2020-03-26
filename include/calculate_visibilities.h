#include <math.h>
#include "read_and_write.h"

extern "C" void copy_XYZ_to_GPU(float *d_X_diff, float *d_Y_diff, float *d_Z_diff,
                                float *X_diff_metres, float *Y_diff_metres, float *Z_diff_metres,
                                const int num_baselines);


extern "C" void calculate_visibilities(float *X_diff_metres, float *Y_diff_metres, float *Z_diff_metres,
                    catsource_t catsource, float *angles_array,
                    const int num_baselines, const int num_visis,
                    visibility_set_t *visibility_set, float *sbf);
