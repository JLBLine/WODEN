#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include "constants.h"
#include "cudacheck.h"
#include "woden_precision_defs.h"

#include "gpu_macros.h"

__device__ void calc_uvw(double *d_X_diff, double *d_Y_diff,
                         double *d_Z_diff,
                         double sdec0, double cdec0,
                         double sha0, double cha0,
                         int iBaseline, int num_baselines, int num_times,
                         int num_freqs,
                         user_precision_t * u, user_precision_t * v,
                         user_precision_t * w) {

  int mod_baseline = iBaseline - num_baselines*floorf((float)iBaseline / (float)num_baselines);
  int time_ind = floorf(((float)iBaseline - (float)mod_baseline) / ((float)num_freqs*(float)num_baselines));

  int xyz_ind = time_ind*num_baselines + mod_baseline;

  // printf("Inputs %d %d %d\n", num_baselines, num_times, num_freqs);
  //
  // printf("part 1 %d\n",(int)((float)iBaseline - (float)mod_baseline) );

  * u = (sha0*d_X_diff[xyz_ind]) + (cha0*d_Y_diff[xyz_ind]);
  * v = -(sdec0*cha0*d_X_diff[xyz_ind]) + (sdec0*sha0*d_Y_diff[xyz_ind]) + (cdec0*d_Z_diff[xyz_ind]);
  * w = (cdec0*cha0*d_X_diff[xyz_ind]) - (cdec0*sha0*d_Y_diff[xyz_ind]) + (sdec0*d_Z_diff[xyz_ind]);

}

__global__ void kern_calc_uvw(double *d_X_diff, double *d_Y_diff,
           double *d_Z_diff, user_precision_t *d_u_metres,
           user_precision_t *d_v_metres, user_precision_t *d_w_metres,
           user_precision_t *d_u, user_precision_t *d_v, user_precision_t *d_w, user_precision_t *d_wavelengths,
           double sdec0, double cdec0,
           double *d_cha0s, double *d_sha0s,
           int num_cross, int num_baselines, int num_times, int num_freqs){
  // Start by computing which baseline we're going to do
  const int iBaseline = threadIdx.x + (blockDim.x*blockIdx.x);

  if (iBaseline < num_cross){
    user_precision_t u, v, w;

    double d_sha0 = d_sha0s[iBaseline];
    double d_cha0 = d_cha0s[iBaseline];

    calc_uvw(d_X_diff, d_Y_diff, d_Z_diff,
               sdec0, cdec0, d_sha0, d_cha0,
               iBaseline, num_baselines, num_times, num_freqs,
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
__global__ void kern_calc_uv_shapelet(double *d_X_diff,
      double *d_Y_diff, double *d_Z_diff,
      user_precision_t *d_u_shapes, user_precision_t *d_v_shapes,
      double *d_lsts, double *d_ras, double *d_decs,
      const int num_baselines, const int num_times, const int num_shapes) {
  // Start by computing which baseline we're going to do
  const int iBaseline = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iComponent = threadIdx.y + (blockDim.y*blockIdx.y);

  if(iBaseline < num_baselines*num_times && iComponent < num_shapes) {

    int mod_baseline = iBaseline - num_baselines*floorf((float)iBaseline / (float)num_baselines);
    int time_ind = floorf(((float)iBaseline - (float)mod_baseline) / (float)num_baselines);
    int xyz_ind = time_ind*num_baselines + mod_baseline;

    // user_precision_t u_shape, v_shape, w_shape;
    //TODO do the sin/cos outside of the GPU kernel?
    double sdec0 = sin(d_decs[iComponent]);
    double cdec0 = cos(d_decs[iComponent]);
    double sha0 = sin(d_lsts[time_ind] - d_ras[iComponent]);
    double cha0 = cos(d_lsts[time_ind] - d_ras[iComponent]);

    user_precision_t u_shape = (sha0*d_X_diff[xyz_ind]) + (cha0*d_Y_diff[xyz_ind]);
    user_precision_t v_shape = -(sdec0*cha0*d_X_diff[xyz_ind]) + (sdec0*sha0*d_Y_diff[xyz_ind]) + (cdec0*d_Z_diff[xyz_ind]);

    int stripe = num_baselines*num_times*iComponent + time_ind*num_baselines + mod_baseline;

    d_u_shapes[stripe] = u_shape;
    d_v_shapes[stripe] = v_shape;
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
  ( gpuMalloc( (void**)&d_ls, num_coords*sizeof(double) ) );

  double *d_ms = NULL;
  ( gpuMalloc( (void**)&d_ms, num_coords*sizeof(double) ) );

  double *d_ns = NULL;
  ( gpuMalloc( (void**)&d_ns, num_coords*sizeof(double) ) );


  double *d_ras = NULL;
  ( gpuMalloc( (void**)&d_ras, num_coords*sizeof(double) ) );
  ( gpuMemcpy(d_ras, ras,
                           num_coords*sizeof(double), gpuMemcpyHostToDevice ) );

  double *d_decs = NULL;
  ( gpuMalloc( (void**)&d_decs, num_coords*sizeof(double) ) );
  ( gpuMemcpy(d_decs, decs,
                           num_coords*sizeof(double), gpuMemcpyHostToDevice ) );

  dim3 grid, threads;

  threads.x = 128;
  grid.x = (int)ceil( (float)num_coords / (float)threads.x );

  gpuErrorCheckKernel("kern_calc_lmn",
          kern_calc_lmn, grid, threads,
          ra0, sin(dec0), cos(dec0),
          d_ras, d_decs, d_ls, d_ms, d_ns,
          num_coords);

  ( gpuMemcpy(ls, d_ls,
                             num_coords*sizeof(double),gpuMemcpyDeviceToHost) );
  ( gpuMemcpy(ms, d_ms,
                             num_coords*sizeof(double),gpuMemcpyDeviceToHost) );
  ( gpuMemcpy(ns, d_ns,
                             num_coords*sizeof(double),gpuMemcpyDeviceToHost) );

  ( gpuFree(d_ls) );
  ( gpuFree(d_ms) );
  ( gpuFree(d_ns) );

  ( gpuFree(d_ras) );
  ( gpuFree(d_decs) );

}

extern "C" void test_kern_calc_uvw(double *X_diff,
   double *Y_diff, double *Z_diff,
   user_precision_t *u_metres, user_precision_t *v_metres, user_precision_t *w_metres,
   user_precision_t *us, user_precision_t *vs, user_precision_t *ws,
   user_precision_t *wavelengths,
   double dec0, double *cha0s, double *sha0s,
   int num_cross, int num_baselines, int num_times, int num_freqs) {

  double *d_X_diff = NULL;
  double *d_Y_diff = NULL;
  double *d_Z_diff = NULL;

  ( gpuMalloc( (void**)&d_X_diff,
                            num_times*num_baselines*sizeof(double) ) );
  ( gpuMemcpy( d_X_diff, X_diff,
    num_times*num_baselines*sizeof(double), gpuMemcpyHostToDevice ) );
  ( gpuMalloc( (void**)&d_Y_diff,
                            num_times*num_baselines*sizeof(double) ) );
  ( gpuMemcpy( d_Y_diff, Y_diff,
    num_times*num_baselines*sizeof(double), gpuMemcpyHostToDevice ) );
  ( gpuMalloc( (void**)&d_Z_diff,
                            num_times*num_baselines*sizeof(double) ) );
  ( gpuMemcpy( d_Z_diff, Z_diff,
    num_times*num_baselines*sizeof(double), gpuMemcpyHostToDevice ) );

  double *d_sha0s = NULL;
  double *d_cha0s = NULL;
  user_precision_t *d_wavelengths = NULL;
  ( gpuMalloc( (void**)&d_sha0s,
                                         num_cross*sizeof(double) ) );
  ( gpuMemcpy( d_sha0s, sha0s,
                 num_cross*sizeof(double), gpuMemcpyHostToDevice ) );
  ( gpuMalloc( (void**)&d_cha0s,
                                         num_cross*sizeof(double) ) );
  ( gpuMemcpy( d_cha0s, cha0s,
                 num_cross*sizeof(double), gpuMemcpyHostToDevice ) );
  ( gpuMalloc( (void**)&d_wavelengths,
                                         num_cross*sizeof(user_precision_t) ) );
  ( gpuMemcpy( d_wavelengths, wavelengths,
                 num_cross*sizeof(user_precision_t), gpuMemcpyHostToDevice ) );

  user_precision_t *d_u_metres = NULL;
  user_precision_t *d_v_metres = NULL;
  user_precision_t *d_w_metres = NULL;
  user_precision_t *d_us = NULL;
  user_precision_t *d_vs = NULL;
  user_precision_t *d_ws = NULL;

  ( gpuMalloc( (void**)&d_u_metres, num_cross*sizeof(user_precision_t) ) );
  ( gpuMalloc( (void**)&d_v_metres, num_cross*sizeof(user_precision_t) ) );
  ( gpuMalloc( (void**)&d_w_metres, num_cross*sizeof(user_precision_t) ) );
  ( gpuMalloc( (void**)&d_us, num_cross*sizeof(user_precision_t) ) );
  ( gpuMalloc( (void**)&d_vs, num_cross*sizeof(user_precision_t) ) );
  ( gpuMalloc( (void**)&d_ws, num_cross*sizeof(user_precision_t) ) );

  dim3 grid, threads;

  threads.x = 128;
  grid.x = (int)ceil( (float)num_cross / (float)threads.x );

  gpuErrorCheckKernel("kern_calc_uvw",
          kern_calc_uvw, grid, threads,
          d_X_diff, d_Y_diff, d_Z_diff,
          d_u_metres, d_v_metres, d_w_metres,
          d_us, d_vs, d_ws, d_wavelengths,
          sin(dec0), cos(dec0),
          d_cha0s, d_sha0s,
          num_cross, num_baselines, num_times, num_freqs);

  ( gpuMemcpy(us, d_us,
                   num_cross*sizeof(user_precision_t),gpuMemcpyDeviceToHost) );
  ( gpuMemcpy(vs, d_vs,
                   num_cross*sizeof(user_precision_t),gpuMemcpyDeviceToHost) );
  ( gpuMemcpy(ws, d_ws,
                   num_cross*sizeof(user_precision_t),gpuMemcpyDeviceToHost) );

  ( gpuMemcpy(u_metres, d_u_metres,
                   num_cross*sizeof(user_precision_t),gpuMemcpyDeviceToHost) );
  ( gpuMemcpy(v_metres, d_v_metres,
                   num_cross*sizeof(user_precision_t),gpuMemcpyDeviceToHost) );
  ( gpuMemcpy(w_metres, d_w_metres,
                   num_cross*sizeof(user_precision_t),gpuMemcpyDeviceToHost) );

  ( gpuFree(d_us) );
  ( gpuFree(d_vs) );
  ( gpuFree(d_ws) );

  ( gpuFree(d_u_metres) );
  ( gpuFree(d_v_metres) );
  ( gpuFree(d_w_metres) );

  ( gpuFree(d_sha0s) );
  ( gpuFree(d_cha0s) );

  ( gpuFree(d_X_diff) );
  ( gpuFree(d_Y_diff) );
  ( gpuFree(d_Z_diff) );

}

extern "C" void test_kern_calc_uv_shapelet(double *X_diff,
                     double *Y_diff, double *Z_diff,
                     user_precision_t *u_shapes, user_precision_t *v_shapes,
                     double *lsts,
                     double *ras, double *decs,
                     int num_baselines, int num_times, int num_shapes) {

  double *d_X_diff = NULL;
  double *d_Y_diff = NULL;
  double *d_Z_diff = NULL;

  ( gpuMalloc( (void**)&d_X_diff,
                                    num_times*num_baselines*sizeof(double) ) );
  ( gpuMemcpy( d_X_diff, X_diff,
            num_times*num_baselines*sizeof(double), gpuMemcpyHostToDevice ) );
  ( gpuMalloc( (void**)&d_Y_diff,
                                    num_times*num_baselines*sizeof(double) ) );
  ( gpuMemcpy( d_Y_diff, Y_diff,
            num_times*num_baselines*sizeof(double), gpuMemcpyHostToDevice ) );
  ( gpuMalloc( (void**)&d_Z_diff,
                                    num_times*num_baselines*sizeof(double) ) );
  ( gpuMemcpy( d_Z_diff, Z_diff,
            num_times*num_baselines*sizeof(double), gpuMemcpyHostToDevice ) );

  double *d_lsts = NULL;
  double *d_ras = NULL;
  double *d_decs = NULL;
  ( gpuMalloc( (void**)&d_lsts, num_times*sizeof(double) ) );
  ( gpuMemcpy( d_lsts, lsts,
                      num_times*sizeof(double), gpuMemcpyHostToDevice ) );
  ( gpuMalloc( (void**)&d_ras, num_shapes*sizeof(double) ) );
  ( gpuMemcpy( d_ras, ras,
                      num_shapes*sizeof(double), gpuMemcpyHostToDevice ) );
  ( gpuMalloc( (void**)&d_decs, num_shapes*sizeof(double) ) );
  ( gpuMemcpy( d_decs, decs,
                      num_shapes*sizeof(double), gpuMemcpyHostToDevice ) );

  user_precision_t *d_u_shapes = NULL;
  user_precision_t *d_v_shapes = NULL;

  ( gpuMalloc( (void**)&d_u_shapes,
                num_shapes*num_baselines*num_times*sizeof(user_precision_t) ) );
  ( gpuMalloc( (void**)&d_v_shapes,
                num_shapes*num_baselines*num_times*sizeof(user_precision_t) ) );

  dim3 grid, threads;

  threads.x = 64;
  grid.x = (int)ceilf( (float)(num_baselines*num_times) / (float)threads.x );

  threads.y = 2;
  grid.y = (int)ceilf( (float)num_shapes / (float)threads.y );

  gpuErrorCheckKernel("kern_calc_uv_shapelet",
          kern_calc_uv_shapelet, grid, threads,
          d_X_diff, d_Y_diff, d_Z_diff,
          d_u_shapes, d_v_shapes,
          d_lsts, d_ras, d_decs,
          num_baselines, num_times, num_shapes);

  ( gpuMemcpy(u_shapes, d_u_shapes,
                   num_shapes*num_baselines*num_times*sizeof(user_precision_t),
                                                    gpuMemcpyDeviceToHost) );
  ( gpuMemcpy(v_shapes, d_v_shapes,
                   num_shapes*num_baselines*num_times*sizeof(user_precision_t),
                                                    gpuMemcpyDeviceToHost) );

  ( gpuFree(d_u_shapes) );
  ( gpuFree(d_v_shapes) );

  ( gpuFree(d_lsts) );
  ( gpuFree(d_ras) );
  ( gpuFree(d_decs) );

  ( gpuFree(d_X_diff) );
  ( gpuFree(d_Y_diff) );
  ( gpuFree(d_Z_diff) );

}
