#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "constants.h"
#include "woden_precision_defs.h"
#include "woden_struct_defs.h"
#include "source_components_gpu.h"
#include "source_components_common.h"

extern "C" void copy_outputs_source_component_common_gpu(int num_of_each_flux_type,
           source_t *d_chunked_source, beam_gains_t *d_beam_gains,
           woden_settings_t *woden_settings,
           beam_settings_t *beam_settings,
           user_precision_complex_t *gxs, user_precision_complex_t *Dxs,
           user_precision_complex_t *Dys, user_precision_complex_t *gys,
           user_precision_t *extrap_flux_I, user_precision_t *extrap_flux_Q,
           user_precision_t *extrap_flux_U, user_precision_t *extrap_flux_V,
           double *ls, double *ms, double *ns,
           e_component_type comptype){

  int NUM_FLUX_TYPES = 3;

  int num_beam_values = NUM_FLUX_TYPES*num_of_each_flux_type*woden_settings->num_freqs*woden_settings->num_time_steps;

  if (woden_settings->use_dipamps == 1) {
    num_beam_values *= woden_settings->num_ants;
  }

  gpuMemcpy(gxs, d_beam_gains->gxs, num_beam_values*sizeof(user_precision_complex_t), gpuMemcpyDeviceToHost );

  gpuMemcpy(gys, d_beam_gains->gys, num_beam_values*sizeof(user_precision_complex_t), gpuMemcpyDeviceToHost );

  if (beam_settings->beamtype == FEE_BEAM || beam_settings->beamtype == FEE_BEAM_INTERP ||
      beam_settings->beamtype == MWA_ANALY || beam_settings->beamtype == EB_OSKAR ||
      beam_settings->beamtype == EB_LOFAR || beam_settings->beamtype == EB_MWA) {
    gpuMemcpy(Dxs, d_beam_gains->Dxs, num_beam_values*sizeof(user_precision_complex_t), gpuMemcpyDeviceToHost );
    gpuMemcpy(Dys, d_beam_gains->Dys, num_beam_values*sizeof(user_precision_complex_t), gpuMemcpyDeviceToHost );
  }

//   Just a little shorthand so don't have to keep writing out as much in the
//   memcpy below

  components_t d_components;

  if (comptype == POINT) {
    d_components = d_chunked_source->point_components;
  }
  else if (comptype == GAUSSIAN) {
    d_components = d_chunked_source->gauss_components;
  }
  else {
    d_components = d_chunked_source->shape_components;
  }


  gpuMemcpy(ls, d_components.ls, NUM_FLUX_TYPES*num_of_each_flux_type*sizeof(double),
                                                              gpuMemcpyDeviceToHost );
  gpuMemcpy(ms, d_components.ms, NUM_FLUX_TYPES*num_of_each_flux_type*sizeof(double),
                                                              gpuMemcpyDeviceToHost );
  gpuMemcpy(ns, d_components.ns, NUM_FLUX_TYPES*num_of_each_flux_type*sizeof(double),
                                                              gpuMemcpyDeviceToHost );

  //until we get RM synthesis working, do this for testing
  //do this because I don't want to cut out all the memcpying below, laaazy
  dim3 grid, threads;

  int num_things = NUM_FLUX_TYPES*num_of_each_flux_type*woden_settings->num_freqs;

  threads.x = 128;
  threads.y = 1;
  grid.x = (int)ceil( (float)num_things / (float)threads.x );
  grid.y = 1;

  if (d_components.do_QUV == 0) {

    gpuMalloc( (void**)&d_components.extrap_stokesQ, num_things*sizeof(user_precision_t) );
    gpuMalloc( (void**)&d_components.extrap_stokesU, num_things*sizeof(user_precision_t) );
    gpuMalloc( (void**)&d_components.extrap_stokesV, num_things*sizeof(user_precision_t) );

    gpuErrorCheckKernel("kern_make_zeros_user_precision",
            kern_make_zeros_user_precision, grid, threads,
            d_components.extrap_stokesQ, num_things);
    gpuErrorCheckKernel("kern_make_zeros_user_precision",
              kern_make_zeros_user_precision, grid, threads,
              d_components.extrap_stokesU, num_things);
    gpuErrorCheckKernel("kern_make_zeros_user_precision",
              kern_make_zeros_user_precision, grid, threads,
              d_components.extrap_stokesV, num_things);
  }

  gpuMemcpy(extrap_flux_I, d_components.extrap_stokesI,
  NUM_FLUX_TYPES*num_of_each_flux_type*woden_settings->num_freqs*sizeof(user_precision_t),
                                                      gpuMemcpyDeviceToHost );
  gpuMemcpy(extrap_flux_Q, d_components.extrap_stokesQ,
  NUM_FLUX_TYPES*num_of_each_flux_type*woden_settings->num_freqs*sizeof(user_precision_t),
                                                      gpuMemcpyDeviceToHost );
  gpuMemcpy(extrap_flux_U, d_components.extrap_stokesU,
  NUM_FLUX_TYPES*num_of_each_flux_type*woden_settings->num_freqs*sizeof(user_precision_t),
                                                      gpuMemcpyDeviceToHost );
  gpuMemcpy(extrap_flux_V, d_components.extrap_stokesV,
  NUM_FLUX_TYPES*num_of_each_flux_type*woden_settings->num_freqs*sizeof(user_precision_t),
                                                      gpuMemcpyDeviceToHost );

//   free_extrapolated_flux_arrays(&d_components);
//   free_beam_gains_gpu(*d_beam_gains, beam_settings->beamtype);
//   free_components_gpu(d_chunked_source, comptype);


}