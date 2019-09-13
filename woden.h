#include <math.h>

#define VELC 299792458.0
#define SOLAR2SIDEREAL 1.00274
#define D2R M_PI/180.0

typedef struct _d_XYZ_t {
    float *d_X_diff;
    float *d_Y_diff;
    float *d_Z_diff;
} d_XYZ_t;
