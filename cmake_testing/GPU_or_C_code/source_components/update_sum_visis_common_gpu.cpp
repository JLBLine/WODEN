#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "constants.h"
#include "woden_precision_defs.h"
#include "woden_struct_defs.h"
#include "source_components_gpu.h"
#include "source_components_common.h"

extern "C" void test_kern_update_sum_visis(int num_freqs, int num_cross,
          int num_baselines, int num_components, int num_times, int beamtype,
          int use_twobeams, int num_ants, int off_cardinal_dipoles,
          user_precision_complex_t *primay_beam_J00,
          user_precision_complex_t *primay_beam_J01,
          user_precision_complex_t *primay_beam_J10,
          user_precision_complex_t *primay_beam_J11,
          user_precision_complex_t *visi_components,
          user_precision_t *flux_I, user_precision_t *flux_Q,
          user_precision_t *flux_U, user_precision_t *flux_V,
          user_precision_t *sum_visi_XX_real, user_precision_t *sum_visi_XX_imag,
          user_precision_t *sum_visi_XY_real, user_precision_t *sum_visi_XY_imag,
          user_precision_t *sum_visi_YX_real, user_precision_t *sum_visi_YX_imag,
          user_precision_t *sum_visi_YY_real, user_precision_t *sum_visi_YY_imag){

  user_precision_complex_t *d_gxs = NULL;
  user_precision_complex_t *d_Dxs = NULL;
  user_precision_complex_t *d_Dys = NULL;
  user_precision_complex_t *d_gys = NULL;
  user_precision_complex_t *d_visi_components = NULL;

  //if do_ants, need to malloc more gains

  int num_input_gains = num_freqs*num_times*num_components*num_ants;

  gpuMalloc( (void**)&d_gxs, num_input_gains*sizeof(user_precision_complex_t) );
  gpuMalloc( (void**)&d_Dxs, num_input_gains*sizeof(user_precision_complex_t) );
  gpuMalloc( (void**)&d_Dys, num_input_gains*sizeof(user_precision_complex_t) );
  gpuMalloc( (void**)&d_gys, num_input_gains*sizeof(user_precision_complex_t) );
  gpuMalloc( (void**)&d_visi_components, num_cross*sizeof(user_precision_complex_t) );

  gpuMemcpy(d_gxs, primay_beam_J00, num_input_gains*sizeof(user_precision_complex_t),
                                                        gpuMemcpyHostToDevice );
  gpuMemcpy(d_Dxs, primay_beam_J01, num_input_gains*sizeof(user_precision_complex_t),
                                                        gpuMemcpyHostToDevice );
  gpuMemcpy(d_Dys, primay_beam_J10, num_input_gains*sizeof(user_precision_complex_t),
                                                        gpuMemcpyHostToDevice );
  gpuMemcpy(d_gys, primay_beam_J11, num_input_gains*sizeof(user_precision_complex_t),
                                                        gpuMemcpyHostToDevice );
  gpuMemcpy(d_visi_components, visi_components, num_cross*sizeof(user_precision_complex_t),
                                                        gpuMemcpyHostToDevice );

  user_precision_t *d_flux_I = NULL;
  user_precision_t *d_flux_Q = NULL;
  user_precision_t *d_flux_U = NULL;
  user_precision_t *d_flux_V = NULL;

  gpuMalloc( (void**)&d_flux_I, num_components*num_times*num_freqs*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_flux_Q, num_components*num_times*num_freqs*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_flux_U, num_components*num_times*num_freqs*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_flux_V, num_components*num_times*num_freqs*sizeof(user_precision_t) );

  gpuMemcpy(d_flux_I, flux_I, num_components*num_times*num_freqs*sizeof(user_precision_t),
                                                                  gpuMemcpyHostToDevice );
  gpuMemcpy(d_flux_Q, flux_Q, num_components*num_times*num_freqs*sizeof(user_precision_t),
                                                                  gpuMemcpyHostToDevice );
  gpuMemcpy(d_flux_U, flux_U, num_components*num_times*num_freqs*sizeof(user_precision_t),
                                                                  gpuMemcpyHostToDevice );
  gpuMemcpy(d_flux_V, flux_V, num_components*num_times*num_freqs*sizeof(user_precision_t),
                                                                  gpuMemcpyHostToDevice );

  user_precision_t *d_sum_visi_XX_real = NULL;
  user_precision_t *d_sum_visi_XY_real = NULL;
  user_precision_t *d_sum_visi_YX_real = NULL;
  user_precision_t *d_sum_visi_YY_real = NULL;
  user_precision_t *d_sum_visi_XX_imag = NULL;
  user_precision_t *d_sum_visi_XY_imag = NULL;
  user_precision_t *d_sum_visi_YX_imag = NULL;
  user_precision_t *d_sum_visi_YY_imag = NULL;

  gpuMalloc( (void**)&d_sum_visi_XX_real, num_cross*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_sum_visi_XY_real, num_cross*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_sum_visi_YX_real, num_cross*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_sum_visi_YY_real, num_cross*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_sum_visi_XX_imag, num_cross*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_sum_visi_XY_imag, num_cross*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_sum_visi_YX_imag, num_cross*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_sum_visi_YY_imag, num_cross*sizeof(user_precision_t) );
  gpuMemcpy(d_sum_visi_XX_real, sum_visi_XX_real,
                    num_cross*sizeof(user_precision_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_sum_visi_XY_real, sum_visi_XY_real,
                    num_cross*sizeof(user_precision_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_sum_visi_YX_real, sum_visi_YX_real,
                    num_cross*sizeof(user_precision_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_sum_visi_YY_real, sum_visi_YY_real,
                    num_cross*sizeof(user_precision_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_sum_visi_XX_imag, sum_visi_XX_imag,
                    num_cross*sizeof(user_precision_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_sum_visi_XY_imag, sum_visi_XY_imag,
                    num_cross*sizeof(user_precision_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_sum_visi_YX_imag, sum_visi_YX_imag,
                    num_cross*sizeof(user_precision_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_sum_visi_YY_imag, sum_visi_YY_imag,
                    num_cross*sizeof(user_precision_t), gpuMemcpyHostToDevice );
  


  //These are only needed when we're actually grabbing information for two
  //antennas instead of one
  int *d_ant1_to_baseline_map = NULL;
  int *d_ant2_to_baseline_map = NULL;

  if (use_twobeams == 1) {
    gpuMalloc( (void**)&d_ant1_to_baseline_map, num_baselines*sizeof(int) );
    gpuMalloc( (void**)&d_ant2_to_baseline_map, num_baselines*sizeof(int) );
    //Fill in the indexes of antenna1 and antenna2 for all cross-correlation combos
    fill_ant_to_baseline_mapping_gpu(num_ants, d_ant1_to_baseline_map,
                                           d_ant2_to_baseline_map);
  }

  dim3 grid, threads;

  threads.x = 128;
  grid.x = (int)ceil( (float)num_cross / (float)threads.x );

  gpuErrorCheckKernel("kern_update_sum_visis_stokesIQUV",
                      kern_update_sum_visis_stokesIQUV, grid, threads,
                      num_freqs, num_baselines, num_components, num_times,
                      beamtype, off_cardinal_dipoles,
                      (gpuUserComplex *)d_gxs, (gpuUserComplex *)d_Dxs,
                      (gpuUserComplex *)d_Dys, (gpuUserComplex *)d_gys,
                      d_ant1_to_baseline_map, d_ant2_to_baseline_map, use_twobeams,
                      (gpuUserComplex *)d_visi_components,
                      d_flux_I, d_flux_Q, d_flux_U, d_flux_V,
                      d_sum_visi_XX_real, d_sum_visi_XX_imag,
                      d_sum_visi_XY_real, d_sum_visi_XY_imag,
                      d_sum_visi_YX_real, d_sum_visi_YX_imag,
                      d_sum_visi_YY_real, d_sum_visi_YY_imag );

  gpuMemcpy(sum_visi_XX_real, d_sum_visi_XX_real,
                  num_cross*sizeof(user_precision_t), gpuMemcpyDeviceToHost );
  gpuMemcpy(sum_visi_XY_real, d_sum_visi_XY_real,
                  num_cross*sizeof(user_precision_t), gpuMemcpyDeviceToHost );
  gpuMemcpy(sum_visi_YX_real, d_sum_visi_YX_real,
                  num_cross*sizeof(user_precision_t), gpuMemcpyDeviceToHost );
  gpuMemcpy(sum_visi_YY_real, d_sum_visi_YY_real,
                  num_cross*sizeof(user_precision_t), gpuMemcpyDeviceToHost );
  gpuMemcpy(sum_visi_XX_imag, d_sum_visi_XX_imag,
                  num_cross*sizeof(user_precision_t), gpuMemcpyDeviceToHost );
  gpuMemcpy(sum_visi_XY_imag, d_sum_visi_XY_imag,
                  num_cross*sizeof(user_precision_t), gpuMemcpyDeviceToHost );
  gpuMemcpy(sum_visi_YX_imag, d_sum_visi_YX_imag,
                  num_cross*sizeof(user_precision_t), gpuMemcpyDeviceToHost );
  gpuMemcpy(sum_visi_YY_imag, d_sum_visi_YY_imag,
                  num_cross*sizeof(user_precision_t), gpuMemcpyDeviceToHost );

  gpuFree( d_gxs );
  gpuFree( d_Dxs );
  gpuFree( d_Dys );
  gpuFree( d_gys );
  gpuFree( d_visi_components );
  gpuFree( d_flux_I );
  gpuFree( d_flux_Q );
  gpuFree( d_flux_U );
  gpuFree( d_flux_V );
  gpuFree( d_sum_visi_XX_real );
  gpuFree( d_sum_visi_XY_real );
  gpuFree( d_sum_visi_YX_real );
  gpuFree( d_sum_visi_YY_real );
  gpuFree( d_sum_visi_XX_imag );
  gpuFree( d_sum_visi_XY_imag );
  gpuFree( d_sum_visi_YX_imag );
  gpuFree( d_sum_visi_YY_imag );

  if (use_twobeams == 1) {
    gpuFree( d_ant1_to_baseline_map );
    gpuFree( d_ant2_to_baseline_map );
  }

}