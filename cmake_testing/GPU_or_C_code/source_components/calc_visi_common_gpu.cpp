#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "constants.h"
#include "woden_precision_defs.h"
#include "woden_struct_defs.h"
#include "source_components_gpu.h"
#include "source_components_common.h"

extern "C" void test_kern_calc_visi_all(int n_powers, int n_curves, int n_lists,
          int num_baselines, int num_shape_coeffs,
          int num_freqs, int num_cross, int num_times,
          e_beamtype beamtype, e_component_type comptype,
          components_t *components,
          source_t *d_chunked_source, double *d_extrap_freqs,
          user_precision_t *us, user_precision_t *vs, user_precision_t *ws,
          user_precision_t *u_shapes, user_precision_t *v_shapes,
          user_precision_t *sum_visi_XX_real, user_precision_t *sum_visi_XX_imag,
          user_precision_t *sum_visi_XY_real, user_precision_t *sum_visi_XY_imag,
          user_precision_t *sum_visi_YX_real, user_precision_t *sum_visi_YX_imag,
          user_precision_t *sum_visi_YY_real, user_precision_t *sum_visi_YY_imag,
          user_precision_t *allsteps_wavelengths, user_precision_t *sbf,
          user_precision_complex_t *gxs, user_precision_complex_t *Dxs,
          user_precision_complex_t *Dys, user_precision_complex_t *gys){

  int num_components = n_powers + n_curves + n_lists;

  user_precision_t *d_us = NULL;
  user_precision_t *d_vs = NULL;
  user_precision_t *d_ws = NULL;
  user_precision_t *d_allsteps_wavelengths = NULL;

  gpuMalloc( (void**)&d_us, num_cross*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_vs, num_cross*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_ws, num_cross*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_allsteps_wavelengths, num_cross*sizeof(user_precision_t) );

  gpuMemcpy(d_us, us, num_cross*sizeof(user_precision_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_vs, vs, num_cross*sizeof(user_precision_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_ws, ws, num_cross*sizeof(user_precision_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_allsteps_wavelengths, allsteps_wavelengths,
                      num_cross*sizeof(user_precision_t), gpuMemcpyHostToDevice );

  components_t *d_components;

  if (comptype == POINT) {
    d_components = &d_chunked_source->point_components;
    num_components = d_chunked_source->n_points;
  }
  else if (comptype == GAUSSIAN) {
    d_components = &d_chunked_source->gauss_components;
    num_components = d_chunked_source->n_gauss;
  }
  else if (comptype == SHAPELET) {
    d_components = &d_chunked_source->shape_components;
    num_components = d_chunked_source->n_shapes;
  }

  gpuMalloc( (void**)&d_components->ls, num_components*sizeof(double) );
  gpuMalloc( (void**)&d_components->ms, num_components*sizeof(double) );
  gpuMalloc( (void**)&d_components->ns, num_components*sizeof(double) );

  gpuMemcpy(d_components->ls, components->ls, num_components*sizeof(double),
            gpuMemcpyHostToDevice );
  gpuMemcpy(d_components->ms, components->ms, num_components*sizeof(double),
            gpuMemcpyHostToDevice );
  gpuMemcpy(d_components->ns, components->ns, num_components*sizeof(double),
            gpuMemcpyHostToDevice );


  //Something to store the primary beam gains (all 4 pols) in
  beam_gains_t d_beam_gains;
  int num_beam_values = num_components*num_freqs*num_times;
  d_beam_gains.use_twobeams = 0;

  gpuMalloc( (void**)&d_beam_gains.gxs, num_beam_values*sizeof(user_precision_complex_t) );
  gpuMalloc( (void**)&d_beam_gains.Dxs, num_beam_values*sizeof(user_precision_complex_t) );
  gpuMalloc( (void**)&d_beam_gains.Dys, num_beam_values*sizeof(user_precision_complex_t) );
  gpuMalloc( (void**)&d_beam_gains.gys, num_beam_values*sizeof(user_precision_complex_t) );

  gpuMemcpy(d_beam_gains.gxs, gxs, num_beam_values*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_beam_gains.Dxs, Dxs, num_beam_values*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_beam_gains.Dys, Dys, num_beam_values*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice );
  gpuMemcpy(d_beam_gains.gys, gys, num_beam_values*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice );

  // user_precision_t *d_sum_visi_XX_real = NULL;
  // user_precision_t *d_sum_visi_XY_real = NULL;
  // user_precision_t *d_sum_visi_YX_real = NULL;
  // user_precision_t *d_sum_visi_YY_real = NULL;
  // user_precision_t *d_sum_visi_XX_imag = NULL;
  // user_precision_t *d_sum_visi_XY_imag = NULL;
  // user_precision_t *d_sum_visi_YX_imag = NULL;
  // user_precision_t *d_sum_visi_YY_imag = NULL;

  visibility_set_t *d_visibility_set = (visibility_set_t *)malloc(sizeof(visibility_set_t));

  gpuMalloc( (void**)&d_visibility_set->sum_visi_XX_real, num_cross*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_visibility_set->sum_visi_XY_real, num_cross*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_visibility_set->sum_visi_YX_real, num_cross*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_visibility_set->sum_visi_YY_real, num_cross*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_visibility_set->sum_visi_XX_imag, num_cross*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_visibility_set->sum_visi_XY_imag, num_cross*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_visibility_set->sum_visi_YX_imag, num_cross*sizeof(user_precision_t) );
  gpuMalloc( (void**)&d_visibility_set->sum_visi_YY_imag, num_cross*sizeof(user_precision_t) );

  //Make sure the visis start at zero by copying across host versions, which
  //should be set to zero already
  gpuMemcpy( d_visibility_set->sum_visi_XX_real, sum_visi_XX_real,
    num_cross*sizeof(user_precision_t), gpuMemcpyHostToDevice );
  gpuMemcpy( d_visibility_set->sum_visi_XY_real, sum_visi_XY_real,
    num_cross*sizeof(user_precision_t), gpuMemcpyHostToDevice );
  gpuMemcpy( d_visibility_set->sum_visi_YX_real, sum_visi_YX_real,
    num_cross*sizeof(user_precision_t), gpuMemcpyHostToDevice );
  gpuMemcpy( d_visibility_set->sum_visi_YY_real, sum_visi_YY_real,
    num_cross*sizeof(user_precision_t), gpuMemcpyHostToDevice );
  gpuMemcpy( d_visibility_set->sum_visi_XX_imag, sum_visi_XX_imag,
    num_cross*sizeof(user_precision_t), gpuMemcpyHostToDevice );
  gpuMemcpy( d_visibility_set->sum_visi_XY_imag, sum_visi_XY_imag,
    num_cross*sizeof(user_precision_t), gpuMemcpyHostToDevice );
  gpuMemcpy( d_visibility_set->sum_visi_YX_imag, sum_visi_YX_imag,
    num_cross*sizeof(user_precision_t), gpuMemcpyHostToDevice );
  gpuMemcpy( d_visibility_set->sum_visi_YY_imag, sum_visi_YY_imag,
    num_cross*sizeof(user_precision_t), gpuMemcpyHostToDevice );

  //Stick things into the container that the main code uses
  calc_visi_inouts_t *d_calc_visi_inouts = (calc_visi_inouts_t *)malloc(sizeof(calc_visi_inouts_t));

  d_calc_visi_inouts->us = d_us;
  d_calc_visi_inouts->vs = d_vs;
  d_calc_visi_inouts->ws = d_ws;
  d_calc_visi_inouts->allsteps_wavelengths = d_allsteps_wavelengths;

  user_precision_t *d_sbf=NULL;
  user_precision_t *d_u_shapes = NULL;
  user_precision_t *d_v_shapes = NULL;

  if (comptype == SHAPELET) {

    gpuMalloc( (void**)&d_u_shapes,
             num_components*num_baselines*num_times*sizeof(user_precision_t) );
    gpuMalloc( (void**)&d_v_shapes,
             num_components*num_baselines*num_times*sizeof(user_precision_t) );

    gpuMemcpy(d_u_shapes, u_shapes,
                 num_components*num_baselines*num_times*sizeof(user_precision_t),
                                                       gpuMemcpyHostToDevice );
    gpuMemcpy(d_v_shapes, v_shapes,
                 num_components*num_baselines*num_times*sizeof(user_precision_t),
                                                       gpuMemcpyHostToDevice );

    gpuMalloc( (void**)&(d_sbf), sbf_N*sbf_L*sizeof(user_precision_t) );
    gpuMemcpy( d_sbf, sbf, sbf_N*sbf_L*sizeof(user_precision_t),
                        gpuMemcpyHostToDevice );

    d_calc_visi_inouts->u_shapes = d_u_shapes;
    d_calc_visi_inouts->v_shapes = d_v_shapes;
    d_calc_visi_inouts->sbf = d_sbf;

  }

  woden_settings_t *woden_settings = (woden_settings_t *)malloc(sizeof(woden_settings_t));
  woden_settings->num_cross = num_cross;
  woden_settings->num_baselines = num_baselines;
  woden_settings->num_time_steps = num_times;
  woden_settings->num_visis = num_cross;
  woden_settings->num_freqs = num_freqs;



  if (comptype == POINT || comptype == GAUSSIAN ) {

    calc_visi_point_or_gauss_gpu(*d_components, d_beam_gains,
                                  d_calc_visi_inouts, d_visibility_set,
                                  num_components, beamtype, comptype,
                                  woden_settings);


  }
  else if (comptype == SHAPELET) {
      calc_visi_shapelets_gpu(*d_components, d_beam_gains,
                               d_calc_visi_inouts, d_visibility_set,
                               num_components, num_shape_coeffs,
                               beamtype, woden_settings);
  }

  gpuMemcpy(sum_visi_XX_real, d_visibility_set->sum_visi_XX_real,
                             num_cross*sizeof(user_precision_t), gpuMemcpyDeviceToHost );
  gpuMemcpy(sum_visi_XY_real, d_visibility_set->sum_visi_XY_real,
                             num_cross*sizeof(user_precision_t), gpuMemcpyDeviceToHost );
  gpuMemcpy(sum_visi_YX_real, d_visibility_set->sum_visi_YX_real,
                             num_cross*sizeof(user_precision_t), gpuMemcpyDeviceToHost );
  gpuMemcpy(sum_visi_YY_real, d_visibility_set->sum_visi_YY_real,
                             num_cross*sizeof(user_precision_t), gpuMemcpyDeviceToHost );
  gpuMemcpy(sum_visi_XX_imag, d_visibility_set->sum_visi_XX_imag,
                             num_cross*sizeof(user_precision_t), gpuMemcpyDeviceToHost );
  gpuMemcpy(sum_visi_XY_imag, d_visibility_set->sum_visi_XY_imag,
                             num_cross*sizeof(user_precision_t), gpuMemcpyDeviceToHost );
  gpuMemcpy(sum_visi_YX_imag, d_visibility_set->sum_visi_YX_imag,
                             num_cross*sizeof(user_precision_t), gpuMemcpyDeviceToHost );
  gpuMemcpy(sum_visi_YY_imag, d_visibility_set->sum_visi_YY_imag,
                             num_cross*sizeof(user_precision_t), gpuMemcpyDeviceToHost );

  gpuFree( d_visibility_set->sum_visi_XX_real );
  gpuFree( d_visibility_set->sum_visi_XY_real );
  gpuFree( d_visibility_set->sum_visi_YX_real );
  gpuFree( d_visibility_set->sum_visi_YY_real );
  gpuFree( d_visibility_set->sum_visi_XX_imag );
  gpuFree( d_visibility_set->sum_visi_XY_imag );
  gpuFree( d_visibility_set->sum_visi_YX_imag );
  gpuFree( d_visibility_set->sum_visi_YY_imag );
  gpuFree( d_allsteps_wavelengths );

  gpuFree( d_beam_gains.gxs );
  gpuFree( d_beam_gains.Dxs );
  gpuFree( d_beam_gains.Dys );
  gpuFree( d_beam_gains.gys );

  gpuFree( d_us );
  gpuFree( d_vs );
  gpuFree( d_ws );

  if (comptype == SHAPELET){
    gpuFree( d_sbf);
    gpuFree( d_u_shapes);
    gpuFree( d_v_shapes);
  }
}