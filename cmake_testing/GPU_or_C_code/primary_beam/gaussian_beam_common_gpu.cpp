#include "woden_precision_defs.h"
#include "gpucomplex.h"
#include "primary_beam_gpu.h"

extern "C" void test_kern_gaussian_beam(double *beam_ls, double *beam_ms,
           double beam_ref_freq, double *freqs,
           user_precision_t fwhm_lm, user_precision_t cos_theta, user_precision_t sin_theta, user_precision_t sin_2theta,
           int num_freqs, int num_time_steps, int num_components,
           user_precision_complex_t *primay_beam_J00, user_precision_complex_t *primay_beam_J11) {

  int num_beam_hadec = num_components * num_time_steps;

  double *d_beam_ls = NULL;
  gpuMalloc( (void**)&d_beam_ls, num_beam_hadec*sizeof(double));
  gpuMemcpy(d_beam_ls, beam_ls, num_beam_hadec*sizeof(double), gpuMemcpyHostToDevice );

  double *d_beam_ms = NULL;
  gpuMalloc( (void**)&d_beam_ms, num_beam_hadec*sizeof(double));
  gpuMemcpy(d_beam_ms, beam_ms, num_beam_hadec*sizeof(double), gpuMemcpyHostToDevice );

  double *d_freqs = NULL;
  gpuMalloc( (void**)&d_freqs, num_freqs*sizeof(double) );
  gpuMemcpy(d_freqs, freqs, num_freqs*sizeof(double), gpuMemcpyHostToDevice );

  user_precision_complex_t *d_g1xs = NULL;
  gpuMalloc( (void**)&d_g1xs, num_freqs*num_beam_hadec*sizeof(user_precision_complex_t));

  user_precision_complex_t *d_g1ys = NULL;
  gpuMalloc( (void**)&d_g1ys, num_freqs*num_beam_hadec*sizeof(user_precision_complex_t));

  dim3 grid, threads;

  threads.x = 16;
  grid.x = (int)ceil( (float)num_beam_hadec / (float)threads.x );

  threads.y = 16;
  grid.y = (int)ceil( (float)num_freqs / (float)threads.y );

  gpuErrorCheckKernel("kern_gaussian_beam",
                        kern_gaussian_beam, grid, threads,
                        d_beam_ls, d_beam_ms,
                        beam_ref_freq, d_freqs,
                        fwhm_lm, cos_theta, sin_theta, sin_2theta,
                        num_freqs, num_time_steps, num_components,
                        (gpuUserComplex *)d_g1xs,
                        (gpuUserComplex *)d_g1ys);

  gpuMemcpy(primay_beam_J00, d_g1xs,
                 num_freqs*num_beam_hadec*sizeof(user_precision_complex_t), gpuMemcpyDeviceToHost);
  gpuMemcpy(primay_beam_J11, d_g1ys,
                 num_freqs*num_beam_hadec*sizeof(user_precision_complex_t), gpuMemcpyDeviceToHost);

  gpuFree(d_beam_ls);
  gpuFree(d_beam_ms);
  gpuFree(d_freqs);
  gpuFree(d_g1xs);
  gpuFree(d_g1ys);
}


extern "C" void test_calculate_gaussian_beam_gpu(int num_components, int num_time_steps,
     int num_freqs, user_precision_t ha0, user_precision_t sdec0, user_precision_t cdec0,
     user_precision_t fwhm_lm, user_precision_t cos_theta, user_precision_t sin_theta, user_precision_t sin_2theta,
     double beam_ref_freq, double *freqs,
     double *beam_has, double *beam_decs,
     user_precision_complex_t *primay_beam_J00, user_precision_complex_t *primay_beam_J11) {

  int num_beam_hadec = num_components * num_time_steps;

  double *d_freqs = NULL;
  gpuMalloc( (void**)&d_freqs, num_freqs*sizeof(double) );
  gpuMemcpy(d_freqs, freqs,
                           num_freqs*sizeof(double), gpuMemcpyHostToDevice );

  user_precision_complex_t *d_g1xs = NULL;
  gpuMalloc( (void**)&d_g1xs,
                                   num_freqs*num_beam_hadec*sizeof(user_precision_complex_t));

  user_precision_complex_t *d_g1ys = NULL;
  gpuMalloc( (void**)&d_g1ys,
                                   num_freqs*num_beam_hadec*sizeof(user_precision_complex_t));


  calculate_gaussian_beam_gpu(num_components, num_time_steps,
                         num_freqs, ha0, sdec0, cdec0,
                         fwhm_lm, cos_theta, sin_theta, sin_2theta,
                         beam_ref_freq, d_freqs,
                         beam_has, beam_decs,
                         (gpuUserComplex *)d_g1xs,
                         (gpuUserComplex *)d_g1ys);

  gpuMemcpy(primay_beam_J00, d_g1xs,
                 num_freqs*num_beam_hadec*sizeof(user_precision_complex_t), gpuMemcpyDeviceToHost);
  gpuMemcpy(primay_beam_J11, d_g1ys,
                 num_freqs*num_beam_hadec*sizeof(user_precision_complex_t), gpuMemcpyDeviceToHost);

  gpuFree(d_freqs);
  gpuFree(d_g1xs);
  gpuFree(d_g1ys);

}