#include <stdio.h>
#include <stdlib.h>
#include <cuComplex.h>
#include <complex.h>
#include <math.h>
#include "constants.h"
#include "cudacheck.h"
#include "woden_precision_defs.h"

__device__ void calc_uvw(double *d_X_diff, double *d_Y_diff,
                         double *d_Z_diff,
                         double sdec0, double cdec0,
                         double sha0, double cha0,
                         int iBaseline, int num_baselines,
                         user_precision_t * u, user_precision_t * v,
                         user_precision_t * w) {

  int mod_baseline = iBaseline - num_baselines*floorf((float)iBaseline / (float)num_baselines);

  * u = (sha0*d_X_diff[mod_baseline]) + (cha0*d_Y_diff[mod_baseline]);
  * v = -(sdec0*cha0*d_X_diff[mod_baseline]) + (sdec0*sha0*d_Y_diff[mod_baseline]) + (cdec0*d_Z_diff[mod_baseline]);
  * w = (cdec0*cha0*d_X_diff[mod_baseline]) - (cdec0*sha0*d_Y_diff[mod_baseline]) + (sdec0*d_Z_diff[mod_baseline]);

}

__global__ void kern_calc_uvw(double *d_X_diff, double *d_Y_diff,
           double *d_Z_diff, user_precision_t *d_u_metres,
           user_precision_t *d_v_metres, user_precision_t *d_w_metres,
           user_precision_t *d_u, user_precision_t *d_v, user_precision_t *d_w, user_precision_t *d_wavelengths,
           double sdec0, double cdec0,
           double *d_cha0s, double *d_sha0s,
           int num_visis, int num_baselines){
  // Start by computing which baseline we're going to do
  const int iBaseline = threadIdx.x + (blockDim.x*blockIdx.x);

  if (iBaseline < num_visis){
    user_precision_t u, v, w;

    double d_sha0 = d_sha0s[iBaseline];
    double d_cha0 = d_cha0s[iBaseline];

    calc_uvw(d_X_diff, d_Y_diff, d_Z_diff,
               sdec0, cdec0, d_sha0, d_cha0,
               iBaseline, num_baselines,
               &u, &v, &w);

    d_u_metres[iBaseline] = u;
    d_v_metres[iBaseline] = v;
    d_w_metres[iBaseline] = w;

    user_precision_t d_wavelength = d_wavelengths[iBaseline];

    d_u[iBaseline] = u / d_wavelength;
    d_v[iBaseline] = v / d_wavelength;
    d_w[iBaseline] = w / d_wavelength;
  }
}

/*TODO: this might be faster to just loop over the inside the kernel? */
__global__ void kern_calc_uvw_shapelet(double *d_X_diff,
      double *d_Y_diff, double *d_Z_diff,
      user_precision_t *d_u_shapes, user_precision_t *d_v_shapes,
      user_precision_t *d_w_shapes, user_precision_t *d_wavelengths,
      double *d_lsts, double *d_ras, double *d_decs,
      const int num_baselines, const int num_visis,
      const int num_shapes) {
  // Start by computing which baseline we're going to do
  const int iBaseline = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iComponent = threadIdx.y + (blockDim.y*blockIdx.y);

  if(iBaseline < num_visis && iComponent < num_shapes) {

    user_precision_t u_shape, v_shape, w_shape;
    //TODO do the sin/cos outside of the GPU kernel?
    double d_sdec0 = sin(d_decs[iComponent]);
    double d_cdec0 = cos(d_decs[iComponent]);
    double d_sha0 = sin(d_lsts[iBaseline] - d_ras[iComponent]);
    double d_cha0 = cos(d_lsts[iBaseline] - d_ras[iComponent]);
    user_precision_t d_wavelength = d_wavelengths[iBaseline];

    calc_uvw(d_X_diff, d_Y_diff, d_Z_diff,
               d_sdec0, d_cdec0, d_sha0, d_cha0,
               iBaseline, num_baselines,
               &u_shape, &v_shape, &w_shape);

    d_u_shapes[num_visis*iComponent + iBaseline] = u_shape / d_wavelength;
    d_v_shapes[num_visis*iComponent + iBaseline] = v_shape / d_wavelength;
    d_w_shapes[num_visis*iComponent + iBaseline] = w_shape / d_wavelength;
  }
}

__device__ void calc_lmn(double ra0, double sdec0,
                         double cdec0,
                         double ra, double dec,
                         double * l, double * m, double * n){
  double cdec;
  double sdec;
  double cdra;
  double sdra;

  cdec = cos(dec);
  sdec = sin(dec);
  cdra = cos((ra - ra0));
  sdra = sin((ra - ra0));

  * l = cdec*sdra;
  * m = sdec*cdec0 - cdec*sdec0*cdra;
  * n = sdec*sdec0 + cdec*cdec0*cdra;

  //Note we could calculate n this way, which gives exactly zero at the horizon,
  //but anything below the horizon should have a negative n, and this makes n
  //positive everywhere
  // double temp_n = sqrt(1.0 -temp_l*temp_l - temp_m*temp_m );
}

__global__ void kern_calc_lmn(double ra0, double sdec0,
                              double cdec0,
                              double *d_ras, double *d_decs,
                              double *d_l, double *d_m, double *d_n,
                              int num_components){

  const int iComponent = threadIdx.x + (blockDim.x*blockIdx.x);

  if (iComponent < num_components){
    double l, m, n;

    calc_lmn(ra0, sdec0, cdec0,
             d_ras[iComponent], d_decs[iComponent],
             &l, &m, &n);

    d_l[iComponent] = l;
    d_m[iComponent] = m;
    d_n[iComponent] = n;

  }
}

/*******************************************************************************
                 Functions below to be used in unit tests
*******************************************************************************/

extern "C" void test_kern_calc_lmn(double ra0, double dec0,
                                   double *ras, double *decs, int num_coords,
                                   double * ls, double * ms, double * ns) {

  double *d_ls = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_ls, num_coords*sizeof(double) ) );

  double *d_ms = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_ms, num_coords*sizeof(double) ) );

  double *d_ns = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_ns, num_coords*sizeof(double) ) );


  double *d_ras = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_ras, num_coords*sizeof(double) ) );
  cudaErrorCheckCall( cudaMemcpy(d_ras, ras,
                           num_coords*sizeof(double), cudaMemcpyHostToDevice ) );

  double *d_decs = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_decs, num_coords*sizeof(double) ) );
  cudaErrorCheckCall( cudaMemcpy(d_decs, decs,
                           num_coords*sizeof(double), cudaMemcpyHostToDevice ) );

  dim3 grid, threads;

  threads.x = 128;
  grid.x = (int)ceil( (float)num_coords / (float)threads.x );

  cudaErrorCheckKernel("kern_calc_lmn",
          kern_calc_lmn, grid, threads,
          ra0, sin(dec0), cos(dec0),
          d_ras, d_decs, d_ls, d_ms, d_ns,
          num_coords);

  cudaErrorCheckCall( cudaMemcpy(ls, d_ls,
                             num_coords*sizeof(double),cudaMemcpyDeviceToHost) );
  cudaErrorCheckCall( cudaMemcpy(ms, d_ms,
                             num_coords*sizeof(double),cudaMemcpyDeviceToHost) );
  cudaErrorCheckCall( cudaMemcpy(ns, d_ns,
                             num_coords*sizeof(double),cudaMemcpyDeviceToHost) );

  cudaErrorCheckCall( cudaFree(d_ls) );
  cudaErrorCheckCall( cudaFree(d_ms) );
  cudaErrorCheckCall( cudaFree(d_ns) );

  cudaErrorCheckCall( cudaFree(d_ras) );
  cudaErrorCheckCall( cudaFree(d_decs) );

}

extern "C" void test_kern_calc_uvw(double *X_diff,
   double *Y_diff, double *Z_diff,
   user_precision_t *u_metres, user_precision_t *v_metres, user_precision_t *w_metres,
   user_precision_t *us, user_precision_t *vs, user_precision_t *ws,
   user_precision_t *wavelengths,
   double dec0, double *cha0s, double *sha0s,
   int num_visis, int num_baselines) {

  double *d_X_diff = NULL;
  double *d_Y_diff = NULL;
  double *d_Z_diff = NULL;

  cudaErrorCheckCall( cudaMalloc( (void**)&d_X_diff,
                                     num_baselines*sizeof(double) ) );
  cudaErrorCheckCall( cudaMemcpy( d_X_diff, X_diff,
             num_baselines*sizeof(double), cudaMemcpyHostToDevice ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_Y_diff,
                                     num_baselines*sizeof(double) ) );
  cudaErrorCheckCall( cudaMemcpy( d_Y_diff, Y_diff,
             num_baselines*sizeof(double), cudaMemcpyHostToDevice ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_Z_diff,
                                     num_baselines*sizeof(double) ) );
  cudaErrorCheckCall( cudaMemcpy( d_Z_diff, Z_diff,
             num_baselines*sizeof(double), cudaMemcpyHostToDevice ) );

  double *d_sha0s = NULL;
  double *d_cha0s = NULL;
  user_precision_t *d_wavelengths = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_sha0s,
                                         num_visis*sizeof(double) ) );
  cudaErrorCheckCall( cudaMemcpy( d_sha0s, sha0s,
                 num_visis*sizeof(double), cudaMemcpyHostToDevice ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_cha0s,
                                         num_visis*sizeof(double) ) );
  cudaErrorCheckCall( cudaMemcpy( d_cha0s, cha0s,
                 num_visis*sizeof(double), cudaMemcpyHostToDevice ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_wavelengths,
                                         num_visis*sizeof(user_precision_t) ) );
  cudaErrorCheckCall( cudaMemcpy( d_wavelengths, wavelengths,
                 num_visis*sizeof(user_precision_t), cudaMemcpyHostToDevice ) );

  user_precision_t *d_u_metres = NULL;
  user_precision_t *d_v_metres = NULL;
  user_precision_t *d_w_metres = NULL;
  user_precision_t *d_us = NULL;
  user_precision_t *d_vs = NULL;
  user_precision_t *d_ws = NULL;

  cudaErrorCheckCall( cudaMalloc( (void**)&d_u_metres, num_visis*sizeof(user_precision_t) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_v_metres, num_visis*sizeof(user_precision_t) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_w_metres, num_visis*sizeof(user_precision_t) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_us, num_visis*sizeof(user_precision_t) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_vs, num_visis*sizeof(user_precision_t) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_ws, num_visis*sizeof(user_precision_t) ) );

  dim3 grid, threads;

  threads.x = 128;
  grid.x = (int)ceil( (float)num_visis / (float)threads.x );

  cudaErrorCheckKernel("kern_calc_uvw",
          kern_calc_uvw, grid, threads,
          d_X_diff, d_Y_diff, d_Z_diff,
          d_u_metres, d_v_metres, d_w_metres,
          d_us, d_vs, d_ws, d_wavelengths,
          sin(dec0), cos(dec0),
          d_cha0s, d_sha0s,
          num_visis, num_baselines);

  cudaErrorCheckCall( cudaMemcpy(us, d_us,
                   num_visis*sizeof(user_precision_t),cudaMemcpyDeviceToHost) );
  cudaErrorCheckCall( cudaMemcpy(vs, d_vs,
                   num_visis*sizeof(user_precision_t),cudaMemcpyDeviceToHost) );
  cudaErrorCheckCall( cudaMemcpy(ws, d_ws,
                   num_visis*sizeof(user_precision_t),cudaMemcpyDeviceToHost) );

  cudaErrorCheckCall( cudaMemcpy(u_metres, d_u_metres,
                   num_visis*sizeof(user_precision_t),cudaMemcpyDeviceToHost) );
  cudaErrorCheckCall( cudaMemcpy(v_metres, d_v_metres,
                   num_visis*sizeof(user_precision_t),cudaMemcpyDeviceToHost) );
  cudaErrorCheckCall( cudaMemcpy(w_metres, d_w_metres,
                   num_visis*sizeof(user_precision_t),cudaMemcpyDeviceToHost) );

  cudaErrorCheckCall( cudaFree(d_us) );
  cudaErrorCheckCall( cudaFree(d_vs) );
  cudaErrorCheckCall( cudaFree(d_ws) );

  cudaErrorCheckCall( cudaFree(d_u_metres) );
  cudaErrorCheckCall( cudaFree(d_v_metres) );
  cudaErrorCheckCall( cudaFree(d_w_metres) );

  cudaErrorCheckCall( cudaFree(d_sha0s) );
  cudaErrorCheckCall( cudaFree(d_cha0s) );

  cudaErrorCheckCall( cudaFree(d_X_diff) );
  cudaErrorCheckCall( cudaFree(d_Y_diff) );
  cudaErrorCheckCall( cudaFree(d_Z_diff) );

}

extern "C" void test_kern_calc_uvw_shapelet(double *X_diff,
                     double *Y_diff, double *Z_diff,
                     user_precision_t *u_shapes, user_precision_t *v_shapes,
                     user_precision_t *w_shapes, user_precision_t *wavelengths,
                     double *lsts, double *ras, double *decs,
                     int num_baselines, int num_visis, int num_shapes) {

  double *d_X_diff = NULL;
  double *d_Y_diff = NULL;
  double *d_Z_diff = NULL;

  cudaErrorCheckCall( cudaMalloc( (void**)&d_X_diff,
                                     num_baselines*sizeof(double) ) );
  cudaErrorCheckCall( cudaMemcpy( d_X_diff, X_diff,
             num_baselines*sizeof(double), cudaMemcpyHostToDevice ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_Y_diff,
                                     num_baselines*sizeof(double) ) );
  cudaErrorCheckCall( cudaMemcpy( d_Y_diff, Y_diff,
             num_baselines*sizeof(double), cudaMemcpyHostToDevice ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_Z_diff,
                                     num_baselines*sizeof(double) ) );
  cudaErrorCheckCall( cudaMemcpy( d_Z_diff, Z_diff,
             num_baselines*sizeof(double), cudaMemcpyHostToDevice ) );

  double *d_lsts = NULL;
  double *d_ras = NULL;
  double *d_decs = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_lsts, num_visis*sizeof(double) ) );
  cudaErrorCheckCall( cudaMemcpy( d_lsts, lsts,
                      num_visis*sizeof(double), cudaMemcpyHostToDevice ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_ras, num_visis*sizeof(double) ) );
  cudaErrorCheckCall( cudaMemcpy( d_ras, ras,
                      num_visis*sizeof(double), cudaMemcpyHostToDevice ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_decs, num_visis*sizeof(double) ) );
  cudaErrorCheckCall( cudaMemcpy( d_decs, decs,
                      num_visis*sizeof(double), cudaMemcpyHostToDevice ) );

  user_precision_t *d_wavelengths = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_wavelengths, num_visis*sizeof(user_precision_t) ) );
  cudaErrorCheckCall( cudaMemcpy( d_wavelengths, wavelengths,
                      num_visis*sizeof(user_precision_t), cudaMemcpyHostToDevice ) );

  user_precision_t *d_u_shapes = NULL;
  user_precision_t *d_v_shapes = NULL;
  user_precision_t *d_w_shapes = NULL;

  cudaErrorCheckCall( cudaMalloc( (void**)&d_u_shapes,
                              num_shapes*num_visis*sizeof(user_precision_t) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_v_shapes,
                              num_shapes*num_visis*sizeof(user_precision_t) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_w_shapes,
                              num_shapes*num_visis*sizeof(user_precision_t) ) );

  dim3 grid, threads;

  threads.x = 64;
  grid.x = (int)ceil( (float)num_visis / (float)threads.x );

  threads.y = 2;
  grid.y = (int)ceil( (float)num_shapes / (float)threads.y );

  cudaErrorCheckKernel("kern_calc_uvw_shapelet",
          kern_calc_uvw_shapelet, grid, threads,
          d_X_diff, d_Y_diff, d_Z_diff,
          d_u_shapes, d_v_shapes, d_w_shapes, d_wavelengths,
          d_lsts, d_ras, d_decs,
          num_baselines, num_visis, num_shapes);

  cudaErrorCheckCall( cudaMemcpy(u_shapes, d_u_shapes,
         num_shapes*num_visis*sizeof(user_precision_t),cudaMemcpyDeviceToHost) );
  cudaErrorCheckCall( cudaMemcpy(v_shapes, d_v_shapes,
         num_shapes*num_visis*sizeof(user_precision_t),cudaMemcpyDeviceToHost) );
  cudaErrorCheckCall( cudaMemcpy(w_shapes, d_w_shapes,
         num_shapes*num_visis*sizeof(user_precision_t),cudaMemcpyDeviceToHost) );

  cudaErrorCheckCall( cudaFree(d_u_shapes) );
  cudaErrorCheckCall( cudaFree(d_v_shapes) );
  cudaErrorCheckCall( cudaFree(d_w_shapes) );

  cudaErrorCheckCall( cudaFree(d_lsts) );
  cudaErrorCheckCall( cudaFree(d_ras) );
  cudaErrorCheckCall( cudaFree(d_decs) );

  cudaErrorCheckCall( cudaFree(d_X_diff) );
  cudaErrorCheckCall( cudaFree(d_Y_diff) );
  cudaErrorCheckCall( cudaFree(d_Z_diff) );

}
