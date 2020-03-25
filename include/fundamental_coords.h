#include <math.h>
#include <stdint.h>
#include <fitsio.h>


__global__ void kern_calc_uvw(float *d_X_diff, float *d_Y_diff, float *d_Z_diff,
      float *d_u_metres, float *d_v_metres, float *d_w_metres,
      float *d_u, float *d_v, float *d_w, float *d_wavelengths,
      float *d_angles_array, float *d_cha0s, float *d_sha0s);

__global__ void kern_calc_lmn(float *d_angles_array, float *d_ras, float *d_decs,
                         float *d_l, float *d_m, float *d_n);

// __device__ void extrap_uvw_flux_calc_lmn(float *d_angles_array,
//                 float *d_ras, float *d_decs, float *d_fluxes, float *d_freqs,
//                 float *d_u_metres, float *d_v_metres, float *d_w_metres, float *d_wavelengths,
//                 int iComponent, int iBaseline,
//                 float * l, float * m, float * n, float * u, float * v, float * w, float * extrap_flux);
//
// __device__ void extrap_uvw_flux(float *d_angles_array,
//                 float *d_u_metres, float *d_v_metres, float *d_w_metres, float *d_wavelengths,
//                 float *d_freqs, float *d_fluxes,
//                 int iComponent, int iBaseline,
//                 float * u, float * v, float * w, float * extrap_flux);
