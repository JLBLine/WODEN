#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "woden_precision_defs.h"
#include "woden_struct_defs.h"
#include "source_components_gpu.h"
#include "source_components_common.h"

extern "C" void test_kern_get_beam_gains(int num_freqs, int num_cross,
          int num_baselines, int num_components, int num_times, int beamtype,
          user_precision_complex_t *primay_beam_J00, user_precision_complex_t *primay_beam_J01,
          user_precision_complex_t *primay_beam_J10, user_precision_complex_t *primay_beam_J11,
          user_precision_complex_t *recover_g1x, user_precision_complex_t *recover_D1x,
          user_precision_complex_t *recover_D1y, user_precision_complex_t *recover_g1y,
          user_precision_complex_t *recover_g2x, user_precision_complex_t *recover_D2x,
          user_precision_complex_t *recover_D2y, user_precision_complex_t *recover_g2y,
          int use_twobeams, int num_ants){

  user_precision_complex_t *d_recover_g1x = NULL;
  user_precision_complex_t *d_recover_D1x = NULL;
  user_precision_complex_t *d_recover_D1y = NULL;
  user_precision_complex_t *d_recover_g1y = NULL;
  user_precision_complex_t *d_recover_g2x = NULL;
  user_precision_complex_t *d_recover_D2x = NULL;
  user_precision_complex_t *d_recover_D2y = NULL;
  user_precision_complex_t *d_recover_g2y = NULL;

  user_precision_complex_t *d_g1xs = NULL;
  user_precision_complex_t *d_D1xs = NULL;
  user_precision_complex_t *d_D1ys = NULL;
  user_precision_complex_t *d_g1ys = NULL;

  int num_recover_gains = num_components*num_cross;
  int num_input_gains = num_freqs*num_times*num_components*num_ants;

  gpuMalloc( (void**)&d_recover_g1x, num_recover_gains*sizeof(user_precision_complex_t) );
  gpuMalloc( (void**)&d_recover_D1x, num_recover_gains*sizeof(user_precision_complex_t) );
  gpuMalloc( (void**)&d_recover_D1y, num_recover_gains*sizeof(user_precision_complex_t) );
  gpuMalloc( (void**)&d_recover_g1y, num_recover_gains*sizeof(user_precision_complex_t) );
  gpuMalloc( (void**)&d_recover_g2x, num_recover_gains*sizeof(user_precision_complex_t) );
  gpuMalloc( (void**)&d_recover_D2x, num_recover_gains*sizeof(user_precision_complex_t) );
  gpuMalloc( (void**)&d_recover_D2y, num_recover_gains*sizeof(user_precision_complex_t) );
  gpuMalloc( (void**)&d_recover_g2y, num_recover_gains*sizeof(user_precision_complex_t) );
  gpuMalloc( (void**)&d_g1xs, num_input_gains*sizeof(user_precision_complex_t) );
  gpuMalloc( (void**)&d_D1xs, num_input_gains*sizeof(user_precision_complex_t) );
  gpuMalloc( (void**)&d_D1ys, num_input_gains*sizeof(user_precision_complex_t) );
  gpuMalloc( (void**)&d_g1ys, num_input_gains*sizeof(user_precision_complex_t) );
  gpuMemcpy(d_g1xs, primay_beam_J00, num_input_gains*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_D1xs, primay_beam_J01, num_input_gains*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_D1ys, primay_beam_J10, num_input_gains*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_g1ys, primay_beam_J11, num_input_gains*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice );

  dim3 grid, threads;

  threads.x = 128;
  grid.x = (int)ceil( (user_precision_t)num_cross / (user_precision_t)threads.x );

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

  gpuErrorCheckKernel("kern_get_beam_gains",
                      kern_get_beam_gains, grid, threads,
                      num_components, num_baselines,
                      num_freqs, num_cross, num_times, beamtype,
                      (gpuUserComplex *)d_g1xs,
                      (gpuUserComplex *)d_D1xs,
                      (gpuUserComplex *)d_D1ys,
                      (gpuUserComplex *)d_g1ys,
                      (gpuUserComplex *)d_recover_g1x, (gpuUserComplex *)d_recover_D1x,
                      (gpuUserComplex *)d_recover_D1y, (gpuUserComplex *)d_recover_g1y,
                      (gpuUserComplex *)d_recover_g2x, (gpuUserComplex *)d_recover_D2x,
                      (gpuUserComplex *)d_recover_D2y, (gpuUserComplex *)d_recover_g2y,
                      use_twobeams, num_ants,
                      d_ant1_to_baseline_map, d_ant2_to_baseline_map);

  gpuMemcpy(recover_g1x, d_recover_g1x, num_recover_gains*sizeof(user_precision_complex_t), gpuMemcpyDeviceToHost );
  gpuMemcpy(recover_D1x, d_recover_D1x, num_recover_gains*sizeof(user_precision_complex_t), gpuMemcpyDeviceToHost );
  gpuMemcpy(recover_D1y, d_recover_D1y, num_recover_gains*sizeof(user_precision_complex_t), gpuMemcpyDeviceToHost );
  gpuMemcpy(recover_g1y, d_recover_g1y, num_recover_gains*sizeof(user_precision_complex_t), gpuMemcpyDeviceToHost );
  gpuMemcpy(recover_g2x, d_recover_g2x, num_recover_gains*sizeof(user_precision_complex_t), gpuMemcpyDeviceToHost );
  gpuMemcpy(recover_D2x, d_recover_D2x, num_recover_gains*sizeof(user_precision_complex_t), gpuMemcpyDeviceToHost );
  gpuMemcpy(recover_D2y, d_recover_D2y, num_recover_gains*sizeof(user_precision_complex_t), gpuMemcpyDeviceToHost );
  gpuMemcpy(recover_g2y, d_recover_g2y, num_recover_gains*sizeof(user_precision_complex_t), gpuMemcpyDeviceToHost );

  gpuFree( d_recover_g1x );
  gpuFree( d_recover_D1x );
  gpuFree( d_recover_D1y );
  gpuFree( d_recover_g1y );
  gpuFree( d_recover_g2x );
  gpuFree( d_recover_D2x );
  gpuFree( d_recover_D2y );
  gpuFree( d_recover_g2y );
  gpuFree( d_g1xs );
  gpuFree( d_D1xs );
  gpuFree( d_D1ys );
  gpuFree( d_g1ys );

  if (use_twobeams == 1) {
    // free(ant1_to_baseline_map);
    // free(ant2_to_baseline_map);
    ( gpuFree( d_ant1_to_baseline_map ) );
    ( gpuFree( d_ant2_to_baseline_map ) );
  }

}