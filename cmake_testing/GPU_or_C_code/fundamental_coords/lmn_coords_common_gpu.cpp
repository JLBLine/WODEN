#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "woden_precision_defs.h"
#include "woden_struct_defs.h"
#include "fundamental_coords_gpu.h"

extern "C" void test_calc_lmn_for_components_gpu(double *ls, double *ms, double *ns,
                                            double *ras, double *decs,
                                  int num_coords, woden_settings_t *woden_settings) {

  components_t *d_components = (components_t *)malloc(sizeof(components_t));

  gpuMalloc( (void**)&d_components->ras, num_coords*sizeof(double) );
  gpuMemcpy(d_components->ras, ras, num_coords*sizeof(double), gpuMemcpyHostToDevice );

  gpuMalloc( (void**)&d_components->decs, num_coords*sizeof(double) );
  gpuMemcpy(d_components->decs, decs, num_coords*sizeof(double), gpuMemcpyHostToDevice );

  calc_lmn_for_components_gpu(d_components, num_coords, woden_settings);

  gpuMemcpy(ls, d_components->ls, num_coords*sizeof(double),gpuMemcpyDeviceToHost);
  gpuMemcpy(ms, d_components->ms, num_coords*sizeof(double),gpuMemcpyDeviceToHost);
  gpuMemcpy(ns, d_components->ns, num_coords*sizeof(double),gpuMemcpyDeviceToHost);

  gpuFree(d_components->ls);
  gpuFree(d_components->ms);
  gpuFree(d_components->ns);

  gpuFree(d_components->ras);
  gpuFree(d_components->decs);

}

