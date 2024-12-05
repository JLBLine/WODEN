#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "woden_precision_defs.h"
#include "woden_struct_defs.h"
#include "source_components_gpu.h"
#include "source_components_common.h"

extern "C" void test_kern_calc_measurement_equation(int num_components,
          int num_baselines,
          user_precision_t *us, user_precision_t *vs, user_precision_t *ws,
          double *ls, double *ms, double *ns, user_precision_complex_t *visis){

  user_precision_t *d_us = NULL;
  user_precision_t *d_vs = NULL;
  user_precision_t *d_ws = NULL;
  ( gpuMalloc( (void**)&d_us, num_baselines*sizeof(user_precision_t) ));
  ( gpuMalloc( (void**)&d_vs, num_baselines*sizeof(user_precision_t) ));
  ( gpuMalloc( (void**)&d_ws, num_baselines*sizeof(user_precision_t) ));
  ( gpuMemcpy(d_us, us, num_baselines*sizeof(user_precision_t),
                                                        gpuMemcpyHostToDevice ));
  ( gpuMemcpy(d_vs, vs, num_baselines*sizeof(user_precision_t),
                                                        gpuMemcpyHostToDevice ));
  ( gpuMemcpy(d_ws, ws, num_baselines*sizeof(user_precision_t),
                                                        gpuMemcpyHostToDevice ));

  double *d_ls = NULL;
  double *d_ms = NULL;
  double *d_ns = NULL;
  gpuMalloc( (void**)&d_ls, num_components*sizeof(double) );
  gpuMalloc( (void**)&d_ms, num_components*sizeof(double) );
  gpuMalloc( (void**)&d_ns, num_components*sizeof(double) );
  gpuMemcpy(d_ls, ls, num_components*sizeof(double), gpuMemcpyHostToDevice );
  gpuMemcpy(d_ms, ms, num_components*sizeof(double), gpuMemcpyHostToDevice );
  gpuMemcpy(d_ns, ns, num_components*sizeof(double), gpuMemcpyHostToDevice );

  user_precision_complex_t *d_visis = NULL;
  gpuMalloc( (void**)&d_visis, num_baselines*num_components*sizeof(user_precision_complex_t) );

  dim3 grid, threads;

  threads.x = 16;
  threads.y = 16;
  grid.x = (int)ceilf( (float)num_baselines / (float)threads.x );
  grid.y = (int)ceilf( (float)num_components / (float)threads.y );

  gpuErrorCheckKernel("kern_calc_measurement_equation",
                      kern_calc_measurement_equation, grid, threads,
                      num_components, num_baselines,
                      d_us, d_vs, d_ws,
                      d_ls, d_ms, d_ns,
                      (gpuUserComplex*)d_visis );

  ( gpuMemcpy(visis, (user_precision_complex_t*)d_visis, num_components*num_baselines*sizeof(user_precision_complex_t),gpuMemcpyDeviceToHost ));

  ( gpuFree( d_us ) );
  ( gpuFree( d_vs ) );
  ( gpuFree( d_ws ) );
  ( gpuFree( d_ls ) );
  ( gpuFree( d_ms ) );
  ( gpuFree( d_ns ) );
  ( gpuFree(d_visis ) );

}