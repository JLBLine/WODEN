#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "woden_precision_defs.h"
#include "woden_struct_defs.h"
#include "fundamental_coords_gpu.h"


extern "C" void test_calc_uvw_gpu(double *X_diff, double *Y_diff, double *Z_diff,
   user_precision_t *u_metres, user_precision_t *v_metres, user_precision_t *w_metres,
   user_precision_t *us, user_precision_t *vs, user_precision_t *ws,
   user_precision_t *wavelengths, double *cha0s, double *sha0s,
   woden_settings_t *woden_settings) {

  int num_cross = woden_settings->num_cross;
  int num_baselines = woden_settings->num_baselines;
  int num_times = woden_settings->num_time_steps;

  double *d_X_diff = NULL;
  double *d_Y_diff = NULL;
  double *d_Z_diff = NULL;

  gpuMalloc( (void**)&d_X_diff, num_times*num_baselines*sizeof(double) );
  gpuMemcpy( d_X_diff, X_diff, num_times*num_baselines*sizeof(double), gpuMemcpyHostToDevice );
  gpuMalloc( (void**)&d_Y_diff, num_times*num_baselines*sizeof(double) );
  gpuMemcpy( d_Y_diff, Y_diff, num_times*num_baselines*sizeof(double), gpuMemcpyHostToDevice );
  gpuMalloc( (void**)&d_Z_diff, num_times*num_baselines*sizeof(double) );
  gpuMemcpy( d_Z_diff, Z_diff, num_times*num_baselines*sizeof(double), gpuMemcpyHostToDevice );

  double *d_sha0s = NULL;
  double *d_cha0s = NULL;
  user_precision_t *d_wavelengths = NULL;
  gpuMalloc( (void**)&d_sha0s,num_cross*sizeof(double) );
  gpuMemcpy( d_sha0s, sha0s, num_cross*sizeof(double), gpuMemcpyHostToDevice );
  gpuMalloc( (void**)&d_cha0s, num_cross*sizeof(double) );
  gpuMemcpy( d_cha0s, cha0s, num_cross*sizeof(double), gpuMemcpyHostToDevice );
  gpuMalloc( (void**)&d_wavelengths, num_cross*sizeof(user_precision_t) );
  gpuMemcpy( d_wavelengths, wavelengths, num_cross*sizeof(user_precision_t),
                                                        gpuMemcpyHostToDevice );

  user_precision_t *d_u_metres = NULL;
  user_precision_t *d_v_metres = NULL;
  user_precision_t *d_w_metres = NULL;
  user_precision_t *d_us = NULL;
  user_precision_t *d_vs = NULL;
  user_precision_t *d_ws = NULL;

  gpuMalloc( (void**)&d_u_metres, num_cross*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_v_metres, num_cross*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_w_metres, num_cross*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_us, num_cross*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_vs, num_cross*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_ws, num_cross*sizeof(user_precision_t) );

  calc_uvw_gpu(d_X_diff, d_Y_diff, d_Z_diff, d_u_metres, d_v_metres, d_w_metres,
                   d_us, d_vs, d_ws, d_wavelengths,
                   d_cha0s, d_sha0s, woden_settings);

  gpuMemcpy(us, d_us,
                   num_cross*sizeof(user_precision_t),gpuMemcpyDeviceToHost);
  gpuMemcpy(vs, d_vs,
                   num_cross*sizeof(user_precision_t),gpuMemcpyDeviceToHost);
  gpuMemcpy(ws, d_ws,
                   num_cross*sizeof(user_precision_t),gpuMemcpyDeviceToHost);

  gpuMemcpy(u_metres, d_u_metres,
                   num_cross*sizeof(user_precision_t),gpuMemcpyDeviceToHost);
  gpuMemcpy(v_metres, d_v_metres,
                   num_cross*sizeof(user_precision_t),gpuMemcpyDeviceToHost);
  gpuMemcpy(w_metres, d_w_metres,
                   num_cross*sizeof(user_precision_t),gpuMemcpyDeviceToHost);

  gpuFree(d_us);
  gpuFree(d_vs);
  gpuFree(d_ws);

  gpuFree(d_u_metres);
  gpuFree(d_v_metres);
  gpuFree(d_w_metres);

  gpuFree(d_sha0s);
  gpuFree(d_cha0s);

  gpuFree(d_X_diff);
  gpuFree(d_Y_diff);
  gpuFree(d_Z_diff);

}



extern "C" void test_calc_uv_shapelet_gpu(double *X_diff,
                     double *Y_diff, double *Z_diff,
                     user_precision_t *u_shapes, user_precision_t *v_shapes,
                     double *lsts,
                     double *ras, double *decs,
                     int num_baselines, int num_times, int num_shapes) {

  double *d_X_diff = NULL;
  double *d_Y_diff = NULL;
  double *d_Z_diff = NULL;

  gpuMalloc( (void**)&d_X_diff, num_times*num_baselines*sizeof(double) );
  gpuMemcpy( d_X_diff, X_diff, num_times*num_baselines*sizeof(double), gpuMemcpyHostToDevice );
  gpuMalloc( (void**)&d_Y_diff, num_times*num_baselines*sizeof(double) );
  gpuMemcpy( d_Y_diff, Y_diff, num_times*num_baselines*sizeof(double), gpuMemcpyHostToDevice );
  gpuMalloc( (void**)&d_Z_diff, num_times*num_baselines*sizeof(double) );
  gpuMemcpy( d_Z_diff, Z_diff, num_times*num_baselines*sizeof(double), gpuMemcpyHostToDevice );

  double *d_ras = NULL;
  double *d_decs = NULL;
  gpuMalloc( (void**)&d_ras, num_shapes*sizeof(double) );
  gpuMemcpy( d_ras, ras, num_shapes*sizeof(double), gpuMemcpyHostToDevice );
  gpuMalloc( (void**)&d_decs, num_shapes*sizeof(double) );
  gpuMemcpy( d_decs, decs, num_shapes*sizeof(double), gpuMemcpyHostToDevice );

  user_precision_t *d_u_shapes = NULL;
  user_precision_t *d_v_shapes = NULL;

  gpuMalloc( (void**)&d_u_shapes, num_shapes*num_baselines*num_times*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_v_shapes, num_shapes*num_baselines*num_times*sizeof(user_precision_t) );

//   gpuErrorCheckKernel("kern_calc_uv_shapelet",
//           kern_calc_uv_shapelet, grid, threads,
//           d_X_diff, d_Y_diff, d_Z_diff,
//           d_u_shapes, d_v_shapes,
//           d_lsts, d_ras, d_decs,
//           num_baselines, num_times, num_shapes);

  woden_settings_t *woden_settings = (woden_settings_t *)malloc(sizeof(woden_settings_t));
  woden_settings->num_baselines = num_baselines;
  woden_settings->num_time_steps = num_times;
  woden_settings->lsts = lsts;

  calc_uv_shapelet_gpu(d_u_shapes, d_v_shapes, num_shapes, 
                       d_X_diff, d_Y_diff, d_Z_diff,
                       d_ras, d_decs, woden_settings);

  gpuMemcpy(u_shapes, d_u_shapes,
            num_shapes*num_baselines*num_times*sizeof(user_precision_t),
            gpuMemcpyDeviceToHost);
  gpuMemcpy(v_shapes, d_v_shapes,
            num_shapes*num_baselines*num_times*sizeof(user_precision_t),
            gpuMemcpyDeviceToHost);

  gpuFree(d_u_shapes);
  gpuFree(d_v_shapes);
  gpuFree(d_ras);
  gpuFree(d_decs);
  gpuFree(d_X_diff);
  gpuFree(d_Y_diff);
  gpuFree(d_Z_diff);

}
