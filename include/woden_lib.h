#define _USE_MATH_DEFINES
#define M_PI_2_2_LN_2 7.11941466249375271693034 /* pi^2 / (2 log_e(2)) */

#include <math.h>
#include "read_and_write.h"

extern "C" void copy_XYZ_to_GPU(float *d_X_diff, float *d_Y_diff, float *d_Z_diff,
                                float *X_diff_metres, float *Y_diff_metres, float *Z_diff_metres,
                                const int num_baselines);

// float *d_point_ras, float *d_point_decs, float *d_point_fluxes,
extern "C" void Atomic_time_step(float *X_diff_metres, float *Y_diff_metres, float *Z_diff_metres,
                    catsource_t catsource, float *angles_array,
                    const int num_baselines, const int num_visis,
                    visibility_set_t *visibility_set, float *sbf);

// extern "C" void assign_pointsource_on_GPU(catsource_t src);
