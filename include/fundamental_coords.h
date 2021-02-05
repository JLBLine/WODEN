#include <math.h>
#include <stdint.h>
#include <fitsio.h>

__device__ void calc_uvw(float *d_X_diff, float *d_Y_diff, float *d_Z_diff,
           float d_sdec0, float d_cdec0, float d_sha0, float d_cha0,
           int iBaseline, int num_baselines,
           float * u, float * v, float * w);

__global__ void kern_calc_uvw(float *d_X_diff, float *d_Y_diff, float *d_Z_diff,
           float *d_u_metres, float *d_v_metres, float *d_w_metres,
           float *d_u, float *d_v, float *d_w, float *d_wavelengths,
           float sdec0, float cdec0,
           float *d_cha0s, float *d_sha0s,
           int num_visis, int num_baselines);

__device__ void calc_lmn(float ra0, float sdec0, float cdec0,
           float d_ra, float d_dec,
           float * l, float * m, float * n);

__global__ void kern_calc_lmn(float ra0, float sdec0, float cdec0,
           float *d_ras, float *d_decs,
           float *d_l, float *d_m, float *d_n,
           int num_components);

__global__ void kern_calc_uvw_shapelet(float *d_X_diff, float *d_Y_diff, float *d_Z_diff,
      float *d_u_s_metres, float *d_v_s_metres, float *d_w_s_metres,
      float *d_lsts, float *d_ras, float *d_decs,
      const int num_baselines, const int num_visis,
      const int num_shapes);
