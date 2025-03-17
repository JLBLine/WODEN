#include "primary_beam_gpu.h"
#include "woden_struct_defs.h"

extern "C" void test_run_hyperbeam_gpu(int num_components,
           int num_time_steps, int num_freqs, int num_ants,
           uint8_t parallatic,
           struct FEEBeamGpu *gpu_fee_beam,
           double *azs, double *zas,
           double *latitudes,
           user_precision_complex_t *primay_beam_J00,
           user_precision_complex_t *primay_beam_J01,
           user_precision_complex_t *primay_beam_J10,
           user_precision_complex_t *primay_beam_J11){
  uint32_t num_azza = num_components * num_time_steps;
  int num_beam_values = num_azza * num_freqs * num_ants;

  user_precision_complex_t *d_gxs = NULL;
  user_precision_complex_t *d_Dxs = NULL;
  user_precision_complex_t *d_Dys = NULL;
  user_precision_complex_t *d_gys = NULL;

  gpuMalloc( (void**)&d_gxs, num_beam_values*sizeof(user_precision_complex_t));
  gpuMalloc( (void**)&d_Dxs, num_beam_values*sizeof(user_precision_complex_t));
  gpuMalloc( (void**)&d_Dys, num_beam_values*sizeof(user_precision_complex_t));
  gpuMalloc( (void**)&d_gys, num_beam_values*sizeof(user_precision_complex_t));

  double *reordered_azs = (double *)malloc(num_azza*sizeof(double));
  double *reordered_zas = (double *)malloc(num_azza*sizeof(double));

  int stripe_new, stripe_old;

  for (int time_ind = 0; time_ind < num_time_steps; time_ind++) {
    for (int comp_ind = 0; comp_ind < num_components; comp_ind++) {
      stripe_new = time_ind*num_components + comp_ind;
      stripe_old = comp_ind*num_time_steps + time_ind;
      reordered_azs[stripe_new] = azs[stripe_old];
      reordered_zas[stripe_new] = zas[stripe_old];
    }
  }
  run_hyperbeam_gpu(num_components,
             num_time_steps, num_freqs, num_ants,
             parallatic,
             gpu_fee_beam,
             reordered_azs, reordered_zas,
             latitudes,
             (gpuUserComplex *)d_gxs,
             (gpuUserComplex *)d_Dxs,
             (gpuUserComplex *)d_Dys,
             (gpuUserComplex *)d_gys);

  gpuMemcpy(primay_beam_J00, d_gxs, num_beam_values*sizeof(user_precision_complex_t),
                                                              gpuMemcpyDeviceToHost);

  gpuMemcpy(primay_beam_J01, d_Dxs, num_beam_values*sizeof(user_precision_complex_t),
                                                              gpuMemcpyDeviceToHost);

  gpuMemcpy(primay_beam_J10, d_Dys, num_beam_values*sizeof(user_precision_complex_t),
                                                              gpuMemcpyDeviceToHost);

  gpuMemcpy(primay_beam_J11, d_gys, num_beam_values*sizeof(user_precision_complex_t),
                                                              gpuMemcpyDeviceToHost);

  gpuFree(d_gxs);
  gpuFree(d_Dxs);
  gpuFree(d_Dys);
  gpuFree(d_gys);

  free(reordered_azs);
  free(reordered_zas);

}