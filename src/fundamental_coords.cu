#include <stdio.h>
#include <stdlib.h>
#include <cuComplex.h>
#include <complex.h>
#include <math.h>
#include "constants.h"

__device__ void calc_uvw(float *d_X_diff, float *d_Y_diff, float *d_Z_diff,
           float sdec0, float cdec0, float sha0, float cha0,
           int iBaseline, int num_baselines,
           float * u, float * v, float * w) {

  int mod_baseline = iBaseline - num_baselines*floorf((float)iBaseline / (float)num_baselines);

  * u = (sha0*d_X_diff[mod_baseline]) + (cha0*d_Y_diff[mod_baseline]);
  * v = -(sdec0*cha0*d_X_diff[mod_baseline]) + (sdec0*sha0*d_Y_diff[mod_baseline]) + (cdec0*d_Z_diff[mod_baseline]);
  * w = (cdec0*cha0*d_X_diff[mod_baseline]) - (cdec0*sha0*d_Y_diff[mod_baseline]) + (sdec0*d_Z_diff[mod_baseline]);

}

__global__ void kern_calc_uvw(float *d_X_diff, float *d_Y_diff, float *d_Z_diff,
           float *d_u_metres, float *d_v_metres, float *d_w_metres,
           float *d_u, float *d_v, float *d_w, float *d_wavelengths,
           float sdec0, float cdec0,
           float *d_cha0s, float *d_sha0s,
           int num_visis, int num_baselines){
  // Start by computing which baseline we're going to do
  const int iBaseline = threadIdx.x + (blockDim.x*blockIdx.x);

  if (iBaseline < num_visis){
    float u, v, w;

    float d_sha0 = d_sha0s[iBaseline];
    float d_cha0 = d_cha0s[iBaseline];

    calc_uvw(d_X_diff, d_Y_diff, d_Z_diff,
               sdec0, cdec0, d_sha0, d_cha0,
               iBaseline, num_baselines,
               &u, &v, &w);

    d_u_metres[iBaseline] = u;
    d_v_metres[iBaseline] = v;
    d_w_metres[iBaseline] = w;

    float d_wavelength = d_wavelengths[iBaseline];

    d_u[iBaseline] = u / d_wavelength;
    d_v[iBaseline] = v / d_wavelength;
    d_w[iBaseline] = w / d_wavelength;
  }
}

__global__ void kern_calc_uvw_shapelet(float *d_X_diff, float *d_Y_diff, float *d_Z_diff,
      float *d_u_shapes, float *d_v_shapes, float *d_w_shapes, float *d_wavelengths,
      float *d_lsts, float *d_ras, float *d_decs,
      const int num_baselines, const int num_visis,
      const int num_shapes) {
  // Start by computing which baseline we're going to do
  const int iBaseline = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iComponent = threadIdx.y + (blockDim.y*blockIdx.y);

  if(iBaseline < num_visis && iComponent < num_shapes) {

    float u_shape, v_shape, w_shape;
    //TODO do the sin/cos outside of the GPU kernel?
    float d_sdec0 = sinf(d_decs[iComponent]);
    float d_cdec0 = cosf(d_decs[iComponent]);
    float d_sha0 = sinf(d_lsts[iBaseline] - d_ras[iComponent]);
    float d_cha0 = cosf(d_lsts[iBaseline] - d_ras[iComponent]);
    float d_wavelength = d_wavelengths[iBaseline];

    calc_uvw(d_X_diff, d_Y_diff, d_Z_diff,
               d_sdec0, d_cdec0, d_sha0, d_cha0,
               iBaseline, num_baselines,
               &u_shape, &v_shape, &w_shape);

    d_u_shapes[num_visis*iComponent + iBaseline] = u_shape / d_wavelength;
    d_v_shapes[num_visis*iComponent + iBaseline] = v_shape / d_wavelength;
    d_w_shapes[num_visis*iComponent + iBaseline] = w_shape / d_wavelength;
  }
}

__device__ void calc_lmn(float ra0, float sdec0, float cdec0,
                         float ra, float dec,
                         float * l, float * m, float * n){
  float cdec;
  float sdec;
  float cdra;
  float sdra;

  cdec = cosf(dec);
  sdec = sinf(dec);
  cdra = cosf((ra - ra0));
  sdra = sinf((ra - ra0));

  * l = cdec*sdra;
  * m = sdec*cdec0 - cdec*sdec0*cdra;
  * n = sdec*sdec0 + cdec*cdec0*cdra;

}

__global__ void kern_calc_lmn(float ra0, float sdec0, float cdec0,
                              float *d_ras, float *d_decs,
                              float *d_l, float *d_m, float *d_n,
                              int num_components){

  const int iComponent = threadIdx.x + (blockDim.x*blockIdx.x);

  if (iComponent < num_components){
    float l, m, n;

    calc_lmn(ra0, sdec0, cdec0,
             d_ras[iComponent], d_decs[iComponent],
             &l, &m, &n);

    d_l[iComponent] = l;
    d_m[iComponent] = m;
    d_n[iComponent] = n;

  }
}
