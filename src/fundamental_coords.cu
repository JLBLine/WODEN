#include <stdio.h>
#include <stdlib.h>
#include <cuComplex.h>
#include <complex.h>
#include <math.h>
#include "constants.h"

__global__ void kern_calc_uvw(float *d_X_diff, float *d_Y_diff, float *d_Z_diff,
      float *d_u_metres, float *d_v_metres, float *d_w_metres,
      float *d_u, float *d_v, float *d_w, float *d_wavelengths,
      float *d_angles_array, float *d_cha0s, float *d_sha0s) {
  //

  float u, v, w;

  // Start by computing which baseline we're going to do
  const int iBaseline = threadIdx.x + (blockDim.x*blockIdx.x);
  float d_sdec0 = d_angles_array[0];
  float d_cdec0 = d_angles_array[1];
  float d_sha0 = d_sha0s[iBaseline];
  float d_cha0 = d_cha0s[iBaseline];

  int mod_baseline = iBaseline - 8128*floorf(iBaseline / 8128.0);

  u = (d_sha0*d_X_diff[mod_baseline]) + (d_cha0*d_Y_diff[mod_baseline]);
  v = (d_sdec0*d_sha0*d_Y_diff[mod_baseline]) + (d_cdec0*d_Z_diff[mod_baseline]) - (d_sdec0*d_cha0*d_X_diff[mod_baseline]);
  w = (d_cdec0*d_cha0*d_X_diff[mod_baseline]) - (d_cdec0*d_sha0*d_Y_diff[mod_baseline]) + (d_sdec0*d_Z_diff[mod_baseline]);

  d_u_metres[iBaseline] = u;
  d_v_metres[iBaseline] = v;
  d_w_metres[iBaseline] = w;

  float d_wavelength = d_wavelengths[iBaseline];

  d_u[iBaseline] = u / d_wavelength;
  d_v[iBaseline] = v / d_wavelength;
  d_w[iBaseline] = w / d_wavelength;

}

__device__ void calc_lmn(float *d_angles_array, float d_ra, float d_dec,
                         float * l, float * m, float * n){
  float d_sdec0 = d_angles_array[0];
  float d_cdec0 = d_angles_array[1];
  float d_ra0 = d_angles_array[2];

  float cdec;
  float sdec;
  float cdra;
  float sdra;

  cdec = cosf(d_dec);
  sdec = sinf(d_dec);
  cdra = cosf((d_ra - d_ra0));
  sdra = sinf((d_ra - d_ra0));

  * l = cdec*sdra;
  * m = sdec*d_cdec0 - cdec*d_sdec0*cdra;
  * n = sdec*d_sdec0 + cdec*d_cdec0*cdra;

}

__global__ void kern_calc_lmn(float *d_angles_array, float *d_ras, float *d_decs,
                         float *d_l, float *d_m, float *d_n){
  const int iComponent = threadIdx.x + (blockDim.x*blockIdx.x);

  float l, m, n;

  calc_lmn(d_angles_array, d_ras[iComponent], d_decs[iComponent],
            &l, &m, &n);

  d_l[iComponent] = l;
  d_m[iComponent] = m;
  d_n[iComponent] = n;

}
