#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

#include "woden_precision_defs.h"
#include "gpucomplex.h"

#include "calculate_visibilities_gpu.h"
#include "fundamental_coords_gpu.h"
#include "constants.h"
#include "source_components_gpu.h"
#include "source_components_common.h"
#include "primary_beam_gpu.h"
#include "hyperbeam_error.h"

//Helpful C code we are also using
#include "visibility_set.h"
#include "gpu_macros.h"


extern "C" calc_visi_inouts_t * create_calc_visi_inouts_gpu(array_layout_t *array_layout,
        visibility_set_t *visibility_set, visibility_set_t *d_visibility_set,
        user_precision_t *sbf, woden_settings_t *woden_settings,
        int num_shapelets, int use_twobeams) {

  calc_visi_inouts_t *d_calc_visi_inouts = (calc_visi_inouts_t *)malloc(sizeof(calc_visi_inouts_t));

  int num_cross = woden_settings->num_cross;
  int num_baselines = woden_settings->num_baselines;
  int num_time_steps = woden_settings->num_time_steps;
  int num_visis = woden_settings->num_visis;
  int num_freqs = woden_settings->num_freqs;
  int num_ants = woden_settings->num_ants;

  gpuMalloc( (void**)&d_calc_visi_inouts->X_diff, num_time_steps*num_baselines*sizeof(double) );
  gpuMemcpy( d_calc_visi_inouts->X_diff, array_layout->X_diff_metres,
      num_time_steps*num_baselines*sizeof(double), gpuMemcpyHostToDevice );
      
  gpuMalloc( (void**)&d_calc_visi_inouts->Y_diff, num_time_steps*num_baselines*sizeof(double) );
  gpuMemcpy( d_calc_visi_inouts->Y_diff, array_layout->Y_diff_metres,
      num_time_steps*num_baselines*sizeof(double), gpuMemcpyHostToDevice );

  gpuMalloc( (void**)&d_calc_visi_inouts->Z_diff, num_time_steps*num_baselines*sizeof(double) );
  gpuMemcpy( d_calc_visi_inouts->Z_diff, array_layout->Z_diff_metres,
      num_time_steps*num_baselines*sizeof(double), gpuMemcpyHostToDevice );

  gpuMalloc( (void**)&d_calc_visi_inouts->ant_X, num_time_steps*num_ants*sizeof(double) );
  gpuMemcpy( d_calc_visi_inouts->ant_X, array_layout->ant_X,
      num_time_steps*num_ants*sizeof(double), gpuMemcpyHostToDevice );
      
  gpuMalloc( (void**)&d_calc_visi_inouts->ant_Y, num_time_steps*num_ants*sizeof(double) );
  gpuMemcpy( d_calc_visi_inouts->ant_Y, array_layout->ant_Y,
      num_time_steps*num_ants*sizeof(double), gpuMemcpyHostToDevice );

  gpuMalloc( (void**)&d_calc_visi_inouts->ant_Z, num_time_steps*num_ants*sizeof(double) );
  gpuMemcpy( d_calc_visi_inouts->ant_Z, array_layout->ant_Z,
      num_time_steps*num_ants*sizeof(double), gpuMemcpyHostToDevice );

  gpuMalloc( (void**)&d_calc_visi_inouts->allsteps_sha0s, num_cross*sizeof(double) );
  gpuMemcpy( d_calc_visi_inouts->allsteps_sha0s, visibility_set->allsteps_sha0s,
                      num_cross*sizeof(double), gpuMemcpyHostToDevice );

  gpuMalloc( (void**)&d_calc_visi_inouts->allsteps_cha0s, num_cross*sizeof(double) );
  gpuMemcpy( d_calc_visi_inouts->allsteps_cha0s, visibility_set->allsteps_cha0s,
                      num_cross*sizeof(double), gpuMemcpyHostToDevice );

  gpuMalloc( (void**)&d_calc_visi_inouts->allsteps_wavelengths, num_cross*sizeof(user_precision_t) );
  gpuMemcpy( d_calc_visi_inouts->allsteps_wavelengths, visibility_set->allsteps_wavelengths,
                      num_cross*sizeof(user_precision_t), gpuMemcpyHostToDevice );

  gpuMalloc( (void**)&d_calc_visi_inouts->u_metres, num_visis*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_calc_visi_inouts->v_metres, num_visis*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_calc_visi_inouts->w_metres, num_visis*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_calc_visi_inouts->us, num_visis*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_calc_visi_inouts->vs, num_visis*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_calc_visi_inouts->ws, num_visis*sizeof(user_precision_t) );

  gpuMalloc( (void**)&d_visibility_set->sum_visi_XX_real,
                      num_visis*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_visibility_set->sum_visi_XX_imag,
                      num_visis*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_visibility_set->sum_visi_XY_real,
                      num_visis*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_visibility_set->sum_visi_XY_imag,
                      num_visis*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_visibility_set->sum_visi_YX_real,
                      num_visis*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_visibility_set->sum_visi_YX_imag,
                      num_visis*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_visibility_set->sum_visi_YY_real,
                      num_visis*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_visibility_set->sum_visi_YY_imag,
                      num_visis*sizeof(user_precision_t) );

  gpuMalloc( (void**)&d_calc_visi_inouts->freqs, num_freqs*sizeof(double) );
  gpuMemcpy( d_calc_visi_inouts->freqs, visibility_set->channel_frequencies,
                      num_freqs*sizeof(double), gpuMemcpyHostToDevice );

  //if we have shapelets in our sky model, copy the shapelet basis functions
  //into GPU memory
  if (num_shapelets > 0) {
    gpuMalloc( (void**)&(d_calc_visi_inouts->sbf), sbf_N*sbf_L*sizeof(user_precision_t) );
    gpuMemcpy( d_calc_visi_inouts->sbf, sbf, sbf_N*sbf_L*sizeof(user_precision_t),
                        gpuMemcpyHostToDevice );

    //Technically we're assigning more memory than we need here, as we're
    //assigning memory for all shapelet sources in the whole source_catalogue
    //(we could instead assign memory for each chunked source)
    //However, this is a small amount of memory, and it saves multiple
    //mallocs and frees so is probably more efficient
    gpuMalloc( (void**)&d_calc_visi_inouts->u_shapes,
        num_shapelets*num_baselines*num_time_steps*sizeof(user_precision_t));
    gpuMalloc( (void**)&d_calc_visi_inouts->v_shapes,
        num_shapelets*num_baselines*num_time_steps*sizeof(user_precision_t));
  }

  if (1) { // use_twobeams == 1

    gpuMalloc( (void**)&d_calc_visi_inouts->ant1_to_baseline_map, num_baselines*sizeof(int) );
    gpuMalloc( (void**)&d_calc_visi_inouts->ant2_to_baseline_map, num_baselines*sizeof(int) );

    fill_ant_to_baseline_mapping_gpu(num_ants, d_calc_visi_inouts->ant1_to_baseline_map,
                                     d_calc_visi_inouts->ant2_to_baseline_map);
  }
  return d_calc_visi_inouts;
}

extern "C" void set_visibilities_to_zero_gpu(visibility_set_t *d_visibility_set,
                              visibility_set_t *chunk_visibility_set, int num_visis) {
    //ensure d_visibility_set is set entirely to zero by copying the host
    //array values, which have been set explictly to zero above
    gpuMemcpy(d_visibility_set->sum_visi_XX_real,
               chunk_visibility_set->sum_visi_XX_real,
               num_visis*sizeof(user_precision_t), gpuMemcpyHostToDevice );
    gpuMemcpy(d_visibility_set->sum_visi_XX_imag,
               chunk_visibility_set->sum_visi_XX_imag,
               num_visis*sizeof(user_precision_t), gpuMemcpyHostToDevice );
    gpuMemcpy(d_visibility_set->sum_visi_XY_real,
               chunk_visibility_set->sum_visi_XY_real,
               num_visis*sizeof(user_precision_t), gpuMemcpyHostToDevice );
    gpuMemcpy(d_visibility_set->sum_visi_XY_imag,
               chunk_visibility_set->sum_visi_XY_imag,
               num_visis*sizeof(user_precision_t), gpuMemcpyHostToDevice );
    gpuMemcpy(d_visibility_set->sum_visi_YX_real,
               chunk_visibility_set->sum_visi_YX_real,
               num_visis*sizeof(user_precision_t), gpuMemcpyHostToDevice );
    gpuMemcpy(d_visibility_set->sum_visi_YX_imag,
               chunk_visibility_set->sum_visi_YX_imag,
               num_visis*sizeof(user_precision_t), gpuMemcpyHostToDevice );
    gpuMemcpy(d_visibility_set->sum_visi_YY_real,
               chunk_visibility_set->sum_visi_YY_real,
               num_visis*sizeof(user_precision_t), gpuMemcpyHostToDevice );
    gpuMemcpy(d_visibility_set->sum_visi_YY_imag,
               chunk_visibility_set->sum_visi_YY_imag,
               num_visis*sizeof(user_precision_t), gpuMemcpyHostToDevice );
}


extern "C" void copy_gpu_visi_set_to_host(visibility_set_t *d_visibility_set,
                                          calc_visi_inouts_t *d_calc_visi_inouts,
                                          visibility_set_t *chunk_visibility_set,
                                          int num_visis) {
    gpuMemcpy(chunk_visibility_set->sum_visi_XX_real,
               d_visibility_set->sum_visi_XX_real,
               num_visis*sizeof(user_precision_t), gpuMemcpyDeviceToHost );
    gpuMemcpy(chunk_visibility_set->sum_visi_XX_imag,
               d_visibility_set->sum_visi_XX_imag,
               num_visis*sizeof(user_precision_t), gpuMemcpyDeviceToHost );
    gpuMemcpy(chunk_visibility_set->sum_visi_XY_real,
               d_visibility_set->sum_visi_XY_real,
               num_visis*sizeof(user_precision_t), gpuMemcpyDeviceToHost );
    gpuMemcpy(chunk_visibility_set->sum_visi_XY_imag,
               d_visibility_set->sum_visi_XY_imag,
               num_visis*sizeof(user_precision_t), gpuMemcpyDeviceToHost );
    gpuMemcpy(chunk_visibility_set->sum_visi_YX_real,
               d_visibility_set->sum_visi_YX_real,
               num_visis*sizeof(user_precision_t), gpuMemcpyDeviceToHost );
    gpuMemcpy(chunk_visibility_set->sum_visi_YX_imag,
               d_visibility_set->sum_visi_YX_imag,
               num_visis*sizeof(user_precision_t), gpuMemcpyDeviceToHost );
    gpuMemcpy(chunk_visibility_set->sum_visi_YY_real,
               d_visibility_set->sum_visi_YY_real,
               num_visis*sizeof(user_precision_t), gpuMemcpyDeviceToHost );
    gpuMemcpy(chunk_visibility_set->sum_visi_YY_imag,
               d_visibility_set->sum_visi_YY_imag,
               num_visis*sizeof(user_precision_t), gpuMemcpyDeviceToHost );

    gpuMemcpy(chunk_visibility_set->us_metres,
                d_calc_visi_inouts->u_metres,num_visis*sizeof(user_precision_t),
                                                         gpuMemcpyDeviceToHost);
    gpuMemcpy(chunk_visibility_set->vs_metres,
                d_calc_visi_inouts->v_metres,num_visis*sizeof(user_precision_t),
                                                         gpuMemcpyDeviceToHost);
    gpuMemcpy(chunk_visibility_set->ws_metres,
                d_calc_visi_inouts->w_metres,num_visis*sizeof(user_precision_t),
                                                         gpuMemcpyDeviceToHost);
}


extern "C"  void free_calc_visi_inouts_gpu(calc_visi_inouts_t *d_calc_visi_inouts,
                                          visibility_set_t *d_visibility_set,
                                          int num_shapelets, int use_twobeams) {


  gpuFree(d_calc_visi_inouts->X_diff);
  gpuFree(d_calc_visi_inouts->Y_diff);
  gpuFree(d_calc_visi_inouts->Z_diff);
  gpuFree(d_calc_visi_inouts->ant_X);
  gpuFree(d_calc_visi_inouts->ant_Y);
  gpuFree(d_calc_visi_inouts->ant_Z);
  gpuFree(d_calc_visi_inouts->allsteps_sha0s);
  gpuFree(d_calc_visi_inouts->allsteps_cha0s);
  gpuFree(d_calc_visi_inouts->allsteps_wavelengths);
  gpuFree(d_calc_visi_inouts->u_metres);
  gpuFree(d_calc_visi_inouts->v_metres);
  gpuFree(d_calc_visi_inouts->w_metres);
  gpuFree(d_calc_visi_inouts->us);
  gpuFree(d_calc_visi_inouts->vs);
  gpuFree(d_calc_visi_inouts->ws);
  gpuFree(d_visibility_set->sum_visi_XX_real);
  gpuFree(d_visibility_set->sum_visi_XX_imag);
  gpuFree(d_visibility_set->sum_visi_XY_real);
  gpuFree(d_visibility_set->sum_visi_XY_imag);
  gpuFree(d_visibility_set->sum_visi_YX_real);
  gpuFree(d_visibility_set->sum_visi_YX_imag);
  gpuFree(d_visibility_set->sum_visi_YY_real);
  gpuFree(d_visibility_set->sum_visi_YY_imag);
  gpuFree(d_calc_visi_inouts->freqs);

  //if we have shapelets in our sky model, copy the shapelet basis functions
  //into GPU memory
  if (num_shapelets > 0) {
    gpuFree(d_calc_visi_inouts->sbf);
    gpuFree(d_calc_visi_inouts->u_shapes);
    gpuFree(d_calc_visi_inouts->v_shapes);
  }

  if (1) { // use_twobeams == 1
    gpuFree(d_calc_visi_inouts->ant1_to_baseline_map);
    gpuFree(d_calc_visi_inouts->ant2_to_baseline_map);
  }
}