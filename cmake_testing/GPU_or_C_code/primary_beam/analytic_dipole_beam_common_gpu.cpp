#include "woden_precision_defs.h"
#include "gpucomplex.h"
#include "primary_beam_gpu.h"

extern "C" void test_analytic_dipole_beam_gpu(int num_components,
     int num_time_steps, int num_freqs,
     user_precision_t *azs, user_precision_t *zas, double *freqs,
     user_precision_complex_t *analy_beam_X,
     user_precision_complex_t *analy_beam_Y) {

  user_precision_complex_t *d_analy_beam_X = NULL;
  gpuMalloc( (void**)&d_analy_beam_X,
    num_freqs*num_time_steps*num_components*sizeof(user_precision_complex_t));

  user_precision_complex_t *d_analy_beam_Y = NULL;
  gpuMalloc( (void**)&d_analy_beam_Y,
     num_freqs*num_time_steps*num_components*sizeof(user_precision_complex_t));

  double *d_freqs = NULL;
  gpuMalloc( (void**)&d_freqs, num_freqs*sizeof(double));
  gpuMemcpy(d_freqs, freqs, num_freqs*sizeof(double),
                      gpuMemcpyHostToDevice);

  calculate_analytic_dipole_beam_gpu(num_components,
      num_time_steps, num_freqs,
      azs, zas, d_freqs,
      (gpuUserComplex *)d_analy_beam_X, (gpuUserComplex *)d_analy_beam_Y);

  gpuMemcpy(analy_beam_X, d_analy_beam_X,
             num_freqs*num_time_steps*num_components*sizeof(user_precision_complex_t),
             gpuMemcpyDeviceToHost);
  gpuMemcpy(analy_beam_Y, d_analy_beam_Y,
             num_freqs*num_time_steps*num_components*sizeof(user_precision_complex_t),
             gpuMemcpyDeviceToHost);

  gpuFree(d_analy_beam_X);
  gpuFree(d_analy_beam_Y);
  gpuFree(d_freqs);

}