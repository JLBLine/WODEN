#include "woden_precision_defs.h"
#include "gpucomplex.h"
#include "primary_beam_gpu.h"

extern "C" void test_RTS_calculate_MWA_analytic_beam_gpu(int num_components,
     int num_time_steps, int num_freqs,
     user_precision_t *azs, user_precision_t *zas, int *delays,
     double latitude, int norm,
     double *beam_has, double *beam_decs, double *freqs,
     user_precision_complex_t *gxs, user_precision_complex_t *Dxs,
     user_precision_complex_t *Dys, user_precision_complex_t *gys) {

  //Allocate space for beam gains
  user_precision_complex_t *d_gxs = NULL;
  user_precision_complex_t *d_Dxs = NULL;
  user_precision_complex_t *d_Dys = NULL;
  user_precision_complex_t *d_gys = NULL;

  gpuMalloc( (void**)&d_gxs,
        num_freqs*num_time_steps*num_components*sizeof(user_precision_complex_t));
  gpuMalloc( (void**)&d_Dxs,
        num_freqs*num_time_steps*num_components*sizeof(user_precision_complex_t));
  gpuMalloc( (void**)&d_Dys,
        num_freqs*num_time_steps*num_components*sizeof(user_precision_complex_t));
  gpuMalloc( (void**)&d_gys,
        num_freqs*num_time_steps*num_components*sizeof(user_precision_complex_t));

  double *d_freqs = NULL;
  gpuMalloc( (void**)&d_freqs, num_freqs*sizeof(double));
  gpuMemcpy(d_freqs, freqs, num_freqs*sizeof(double), gpuMemcpyHostToDevice);

  //Run
  calculate_RTS_MWA_analytic_beam_gpu(num_components,
       num_time_steps, num_freqs,
       azs, zas, delays, latitude, norm,
       beam_has, beam_decs, d_freqs,
       (gpuUserComplex*)d_gxs, (gpuUserComplex*)d_Dxs,
       (gpuUserComplex*)d_Dys, (gpuUserComplex*)d_gys);

  gpuMemcpy(gxs, d_gxs, num_freqs*num_time_steps*num_components*sizeof(user_precision_complex_t),
                                                                          gpuMemcpyDeviceToHost);
  gpuMemcpy(Dxs, d_Dxs, num_freqs*num_time_steps*num_components*sizeof(user_precision_complex_t),
                                                                          gpuMemcpyDeviceToHost);
  gpuMemcpy(Dys, d_Dys, num_freqs*num_time_steps*num_components*sizeof(user_precision_complex_t),
                                                                          gpuMemcpyDeviceToHost);
  gpuMemcpy(gys, d_gys, num_freqs*num_time_steps*num_components*sizeof(user_precision_complex_t),
                                                                          gpuMemcpyDeviceToHost);

  gpuFree(d_gxs);
  gpuFree(d_Dxs);
  gpuFree(d_Dys);
  gpuFree(d_gys);
  gpuFree(d_freqs);


}