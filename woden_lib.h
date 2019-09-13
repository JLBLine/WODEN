#define _USE_MATH_DEFINES
#include <math.h>
#include "read_and_write.h"

#define VELC 299792458.0
#define SOLAR2SIDEREAL 1.00274
#define D2R M_PI/180.0
#define DEFAULT_SI -0.8

extern "C" void copy_XYZ_to_GPU(float *d_X_diff, float *d_Y_diff, float *d_Z_diff,
                                float *X_diff_metres, float *Y_diff_metres, float *Z_diff_metres,
                                const int num_baselines);

// float *d_point_ras, float *d_point_decs, float *d_point_fluxes,
extern "C" void Atomic_time_step(float *X_diff_metres, float *Y_diff_metres, float *Z_diff_metres,
                    catsource_t catsource, float *angles_array,
                    const int num_baselines, const int num_time_steps,
                    visibility_set_t visibility_set, float *sbf2);

extern "C" void assign_pointsource_on_GPU(catsource_t src);


//const float sdec0, const float cdec0, const float sha0, const float cha0, const float ra0,
