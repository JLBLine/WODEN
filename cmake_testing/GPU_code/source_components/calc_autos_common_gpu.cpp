#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "constants.h"
#include "woden_precision_defs.h"
#include "woden_struct_defs.h"
#include "source_components_gpu.h"
#include "source_components_common.h"

extern "C" void test_kern_calc_autos(components_t *components, int beamtype,
                                     int num_components, int num_baselines,
                                     int num_freqs, int num_times, int num_ants,
                                     int num_beams,
                                     visibility_set_t *visibility_set){

  int off_cardinal_dipoles = 0;

  int use_twobeams = 0;
  if (num_beams > 1) {
    use_twobeams = 1;
  }

  ////malloc on device and copy extrapolated fluxes
  int num_pb_values = num_beams*num_freqs*num_times*num_components;

  int num_autos = num_ants*num_freqs*num_times;
  int num_cross = num_baselines*num_freqs*num_times;
  int num_visis = num_cross + num_autos;

  components_t *d_components = (components_t* )malloc(sizeof(components_t));

  gpuMalloc( (void**)&d_components->extrap_stokesI,
                            num_components*num_freqs*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_components->extrap_stokesQ,
                            num_components*num_freqs*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_components->extrap_stokesU,
                            num_components*num_freqs*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_components->extrap_stokesV,
                            num_components*num_freqs*sizeof(user_precision_t) );

  d_components->do_QUV = components->do_QUV;


  gpuMemcpy(d_components->extrap_stokesI, components->extrap_stokesI,
     num_components*num_freqs*sizeof(user_precision_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_components->extrap_stokesQ, components->extrap_stokesQ,
     num_components*num_freqs*sizeof(user_precision_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_components->extrap_stokesU, components->extrap_stokesU,
     num_components*num_freqs*sizeof(user_precision_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_components->extrap_stokesV, components->extrap_stokesV,
     num_components*num_freqs*sizeof(user_precision_t), gpuMemcpyHostToDevice );
  //
  // //malloc on device and copy beam values
  //
  d_beam_gains_t *d_component_beam_gains = (d_beam_gains_t* )malloc(sizeof(d_beam_gains_t));
  // d_beam_gains_t d_component_beam_gains;
  gpuMalloc( (void**)&d_component_beam_gains->d_gxs,
                                num_pb_values*sizeof(user_precision_complex_t) );
  gpuMalloc( (void**)&d_component_beam_gains->d_Dxs,
                                num_pb_values*sizeof(user_precision_complex_t) );
  gpuMalloc( (void**)&d_component_beam_gains->d_Dys,
                                num_pb_values*sizeof(user_precision_complex_t) );
  gpuMalloc( (void**)&d_component_beam_gains->d_gys,
                                num_pb_values*sizeof(user_precision_complex_t) );

  gpuMemcpy(d_component_beam_gains->d_gxs, components->gxs,
                num_pb_values*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_component_beam_gains->d_Dxs, components->Dxs,
                num_pb_values*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_component_beam_gains->d_Dys, components->Dys,
                num_pb_values*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_component_beam_gains->d_gys, components->gys,
                num_pb_values*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice );


  user_precision_t *d_sum_visi_XX_real;
  user_precision_t *d_sum_visi_XX_imag;
  user_precision_t *d_sum_visi_XY_real;
  user_precision_t *d_sum_visi_XY_imag;
  user_precision_t *d_sum_visi_YX_real;
  user_precision_t *d_sum_visi_YX_imag;
  user_precision_t *d_sum_visi_YY_real;
  user_precision_t *d_sum_visi_YY_imag;

  gpuMalloc( (void**)&d_sum_visi_XX_real, num_visis*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_sum_visi_XX_imag, num_visis*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_sum_visi_XY_real, num_visis*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_sum_visi_XY_imag, num_visis*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_sum_visi_YX_real, num_visis*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_sum_visi_YX_imag, num_visis*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_sum_visi_YY_real, num_visis*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_sum_visi_YY_imag, num_visis*sizeof(user_precision_t) );


  //ensure d_sum_visi_XX_real are set entirely to zero by copying the host
  //array values, which have been set explictly to zero during chunking
  gpuMemcpy(d_sum_visi_XX_real, visibility_set->sum_visi_XX_real,
                    num_visis*sizeof(user_precision_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_sum_visi_XX_imag, visibility_set->sum_visi_XX_imag,
                    num_visis*sizeof(user_precision_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_sum_visi_XY_real, visibility_set->sum_visi_XY_real,
                    num_visis*sizeof(user_precision_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_sum_visi_XY_imag, visibility_set->sum_visi_XY_imag,
                    num_visis*sizeof(user_precision_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_sum_visi_YX_real, visibility_set->sum_visi_YX_real,
                    num_visis*sizeof(user_precision_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_sum_visi_YX_imag, visibility_set->sum_visi_YX_imag,
                    num_visis*sizeof(user_precision_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_sum_visi_YY_real, visibility_set->sum_visi_YY_real,
                    num_visis*sizeof(user_precision_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_sum_visi_YY_imag, visibility_set->sum_visi_YY_imag,
                    num_visis*sizeof(user_precision_t), gpuMemcpyHostToDevice );

  dim3 grid, threads;

  threads.x = 64;
  threads.y = 2;
  threads.z = 1;
  grid.x = (int)ceil( (float)(num_freqs*num_times) / (float)threads.x );
  grid.y = (int)ceil( (float)(num_ants) / (float)threads.y );
  grid.z = 1;

  int *d_ant_to_auto_map = NULL;

  if (use_twobeams == 1) {

    int *ant_to_auto_map = NULL;
    ant_to_auto_map = (int *)malloc(num_ants*sizeof(int));
    for (int ant = 0; ant < num_ants; ant++){
        ant_to_auto_map[ant] = ant;
    }

    gpuMalloc( (void**)&d_ant_to_auto_map, num_ants*sizeof(int) );
    gpuMemcpy(d_ant_to_auto_map, ant_to_auto_map, num_ants*sizeof(int), gpuMemcpyHostToDevice );

    free(ant_to_auto_map);
  }

  gpuErrorCheckKernel("kern_calc_autos",
                kern_calc_autos, grid, threads,
                *d_components, *d_component_beam_gains,
                beamtype, num_components, num_baselines,
                num_freqs, num_times, num_ants,
                d_sum_visi_XX_real, d_sum_visi_XX_imag,
                d_sum_visi_XY_real, d_sum_visi_XY_imag,
                d_sum_visi_YX_real, d_sum_visi_YX_imag,
                d_sum_visi_YY_real, d_sum_visi_YY_imag,
                use_twobeams, d_ant_to_auto_map,
                d_ant_to_auto_map, off_cardinal_dipoles);

  //Copy outputs onto host so we can check our answers
  gpuMemcpy(visibility_set->sum_visi_XX_real,
                         d_sum_visi_XX_real, num_visis*sizeof(user_precision_t),
                                                        gpuMemcpyDeviceToHost );
  gpuMemcpy(visibility_set->sum_visi_XY_real,
                         d_sum_visi_XY_real, num_visis*sizeof(user_precision_t),
                                                        gpuMemcpyDeviceToHost );
  gpuMemcpy(visibility_set->sum_visi_YX_real,
                         d_sum_visi_YX_real, num_visis*sizeof(user_precision_t),
                                                        gpuMemcpyDeviceToHost );
  gpuMemcpy(visibility_set->sum_visi_YY_real,
                         d_sum_visi_YY_real, num_visis*sizeof(user_precision_t),
                                                        gpuMemcpyDeviceToHost );
  gpuMemcpy(visibility_set->sum_visi_XX_imag,
                         d_sum_visi_XX_imag, num_visis*sizeof(user_precision_t),
                                                        gpuMemcpyDeviceToHost );
  gpuMemcpy(visibility_set->sum_visi_XY_imag,
                         d_sum_visi_XY_imag, num_visis*sizeof(user_precision_t),
                                                        gpuMemcpyDeviceToHost );
  gpuMemcpy(visibility_set->sum_visi_YX_imag,
                         d_sum_visi_YX_imag, num_visis*sizeof(user_precision_t),
                                                        gpuMemcpyDeviceToHost );
  gpuMemcpy(visibility_set->sum_visi_YY_imag,
                         d_sum_visi_YY_imag, num_visis*sizeof(user_precision_t),
                                                        gpuMemcpyDeviceToHost );


  gpuFree( d_components->extrap_stokesI );
  gpuFree( d_components->extrap_stokesQ );
  gpuFree( d_components->extrap_stokesU );
  gpuFree( d_components->extrap_stokesV );
  gpuFree( d_component_beam_gains->d_gxs );
  gpuFree( d_component_beam_gains->d_Dxs );
  gpuFree( d_component_beam_gains->d_Dys );
  gpuFree( d_component_beam_gains->d_gys );


  gpuFree( d_sum_visi_XX_real );
  gpuFree( d_sum_visi_XX_imag );
  gpuFree( d_sum_visi_XY_real );
  gpuFree( d_sum_visi_XY_imag );
  gpuFree( d_sum_visi_YX_real );
  gpuFree( d_sum_visi_YX_imag );
  gpuFree( d_sum_visi_YY_real );
  gpuFree( d_sum_visi_YY_imag );

  if (use_twobeams == 1) {
    gpuFree( d_ant_to_auto_map );
  }
}