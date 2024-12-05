#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "woden_precision_defs.h"
#include "woden_struct_defs.h"
#include "source_components_gpu.h"
#include "source_components_common.h"

extern "C" void test_kern_apply_beam_gains(int num_gains, user_precision_complex_t *g1xs,
          user_precision_complex_t *D1xs,
          user_precision_complex_t *D1ys, user_precision_complex_t *g1ys,
          user_precision_complex_t *g2xs, user_precision_complex_t *D2xs,
          user_precision_complex_t *D2ys, user_precision_complex_t *g2ys,
          user_precision_t *flux_Is, user_precision_t *flux_Qs,
          user_precision_t *flux_Us, user_precision_t *flux_Vs,
          user_precision_complex_t *visi_components,
          user_precision_complex_t *visi_XXs, user_precision_complex_t *visi_XYs,
          user_precision_complex_t *visi_YXs, user_precision_complex_t *visi_YYs){

  user_precision_complex_t *d_g1xs = NULL;
  user_precision_complex_t *d_D1xs = NULL;
  user_precision_complex_t *d_D1ys = NULL;
  user_precision_complex_t *d_g1ys = NULL;
  user_precision_complex_t *d_g2xs = NULL;
  user_precision_complex_t *d_D2xs = NULL;
  user_precision_complex_t *d_D2ys = NULL;
  user_precision_complex_t *d_g2ys = NULL;
  user_precision_t *d_flux_Is = NULL;
  user_precision_t *d_flux_Qs = NULL;
  user_precision_t *d_flux_Us = NULL;
  user_precision_t *d_flux_Vs = NULL;
  user_precision_complex_t *d_visi_components = NULL;
  user_precision_complex_t *d_visi_XXs = NULL;
  user_precision_complex_t *d_visi_XYs = NULL;
  user_precision_complex_t *d_visi_YXs = NULL;
  user_precision_complex_t *d_visi_YYs = NULL;

  gpuMalloc( (void**)&d_g1xs,
                                  num_gains*sizeof(user_precision_complex_t) );
  gpuMalloc( (void**)&d_D1xs,
                                  num_gains*sizeof(user_precision_complex_t) );
  gpuMalloc( (void**)&d_D1ys,
                                  num_gains*sizeof(user_precision_complex_t) );
  gpuMalloc( (void**)&d_g1ys,
                                  num_gains*sizeof(user_precision_complex_t) );
  gpuMalloc( (void**)&d_g2xs,
                                  num_gains*sizeof(user_precision_complex_t) );
  gpuMalloc( (void**)&d_D2xs,
                                  num_gains*sizeof(user_precision_complex_t) );
  gpuMalloc( (void**)&d_D2ys,
                                  num_gains*sizeof(user_precision_complex_t) );
  gpuMalloc( (void**)&d_g2ys,
                                  num_gains*sizeof(user_precision_complex_t) );
  gpuMalloc( (void**)&d_flux_Is,
                                          num_gains*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_flux_Qs,
                                          num_gains*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_flux_Us,
                                          num_gains*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_flux_Vs,
                                          num_gains*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_visi_components,
                                  num_gains*sizeof(user_precision_complex_t) );
  gpuMalloc( (void**)&d_visi_XXs,
                                  num_gains*sizeof(user_precision_complex_t) );
  gpuMalloc( (void**)&d_visi_XYs,
                                  num_gains*sizeof(user_precision_complex_t) );
  gpuMalloc( (void**)&d_visi_YXs,
                                  num_gains*sizeof(user_precision_complex_t) );
  gpuMalloc( (void**)&d_visi_YYs,
                                  num_gains*sizeof(user_precision_complex_t) );

  gpuMemcpy(d_g1xs, g1xs,
          num_gains*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_D1xs, D1xs,
          num_gains*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_D1ys, D1ys,
          num_gains*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_g1ys, g1ys,
          num_gains*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_g2xs, g2xs,
          num_gains*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_D2xs, D2xs,
          num_gains*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_D2ys, D2ys,
          num_gains*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_g2ys, g2ys,
          num_gains*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_visi_components, visi_components,
          num_gains*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_visi_XXs, visi_XXs,
          num_gains*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_visi_XYs, visi_XYs,
          num_gains*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_visi_YXs, visi_YXs,
          num_gains*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_visi_YYs, visi_YYs,
          num_gains*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice );

  gpuMemcpy(d_flux_Is, flux_Is,
                             num_gains*sizeof(user_precision_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_flux_Qs, flux_Qs,
                             num_gains*sizeof(user_precision_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_flux_Us, flux_Us,
                             num_gains*sizeof(user_precision_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_flux_Vs, flux_Vs,
                             num_gains*sizeof(user_precision_t), gpuMemcpyHostToDevice );

  dim3 grid, threads;

  threads.x = 128;
  grid.x = (int)ceil( (user_precision_t)num_gains / (user_precision_t)threads.x );

  gpuErrorCheckKernel("kern_apply_beam_gains_stokesIQUV_on_cardinal",
                      kern_apply_beam_gains_stokesIQUV_on_cardinal, grid, threads,
                      num_gains,
                      (gpuUserComplex *)d_g1xs, (gpuUserComplex *)d_D1xs,
                      (gpuUserComplex *)d_D1ys, (gpuUserComplex *)d_g1ys,
                      (gpuUserComplex *)d_g2xs, (gpuUserComplex *)d_D2xs,
                      (gpuUserComplex *)d_D2ys, (gpuUserComplex *)d_g2ys,
                      d_flux_Is, d_flux_Qs,
                      d_flux_Us, d_flux_Vs,
                      (gpuUserComplex *)d_visi_components,
                      (gpuUserComplex *)d_visi_XXs, (gpuUserComplex *)d_visi_XYs,
                      (gpuUserComplex *)d_visi_YXs, (gpuUserComplex *)d_visi_YYs );

  gpuMemcpy(visi_XXs, d_visi_XXs,
           num_gains*sizeof(user_precision_complex_t),gpuMemcpyDeviceToHost );
  gpuMemcpy(visi_XYs, d_visi_XYs,
           num_gains*sizeof(user_precision_complex_t),gpuMemcpyDeviceToHost );
  gpuMemcpy(visi_YXs, d_visi_YXs,
           num_gains*sizeof(user_precision_complex_t),gpuMemcpyDeviceToHost );
  gpuMemcpy(visi_YYs, d_visi_YYs,
           num_gains*sizeof(user_precision_complex_t),gpuMemcpyDeviceToHost );

  gpuFree( d_g1xs );
  gpuFree( d_D1xs );
  gpuFree( d_D1ys );
  gpuFree( d_g1ys );
  gpuFree( d_g2xs );
  gpuFree( d_D2xs );
  gpuFree( d_D2ys );
  gpuFree( d_g2ys );
  gpuFree( d_flux_Is );
  gpuFree( d_flux_Qs );
  gpuFree( d_flux_Us );
  gpuFree( d_flux_Vs );
  gpuFree( d_visi_components );
  gpuFree( d_visi_XXs );
  gpuFree( d_visi_XYs );
  gpuFree( d_visi_YXs );
  gpuFree( d_visi_YYs );

}