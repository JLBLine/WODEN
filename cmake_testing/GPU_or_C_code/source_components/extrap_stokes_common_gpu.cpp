#include "extrap_stokes_common_gpu.h"


extern "C" double * malloc_freqs_gpu(int num_extrap_freqs, double *extrap_freqs){

  double *d_extrap_freqs = NULL;

  gpuMalloc( (void**)&d_extrap_freqs,
                                   num_extrap_freqs*sizeof(double) );
  gpuMemcpy(d_extrap_freqs, extrap_freqs,
             num_extrap_freqs*sizeof(double), gpuMemcpyHostToDevice );

  return d_extrap_freqs;

}

extern "C" void free_freqs_gpu(double *d_extrap_freqs){

  gpuFree( d_extrap_freqs );

}

extern "C" void copy_extrapolated_flux_arrays_to_host(source_t *d_chunked_source,
                                                  int num_extrap_freqs,
                                                  user_precision_t *extrap_flux_I,
                                                  user_precision_t *extrap_flux_Q,
                                                  user_precision_t *extrap_flux_U,
                                                  user_precision_t *extrap_flux_V){

  components_t d_components = d_chunked_source->point_components;

  gpuMemcpy(extrap_flux_I, d_components.extrap_stokesI,
            d_chunked_source->n_points*num_extrap_freqs*sizeof(user_precision_t),
                                                      gpuMemcpyDeviceToHost );
  gpuMemcpy(extrap_flux_Q, d_components.extrap_stokesQ,
            d_chunked_source->n_points*num_extrap_freqs*sizeof(user_precision_t),
                                                      gpuMemcpyDeviceToHost );
  gpuMemcpy(extrap_flux_U, d_components.extrap_stokesU,
            d_chunked_source->n_points*num_extrap_freqs*sizeof(user_precision_t),
                                                      gpuMemcpyDeviceToHost );
  gpuMemcpy(extrap_flux_V, d_components.extrap_stokesV,
            d_chunked_source->n_points*num_extrap_freqs*sizeof(user_precision_t),
                                                      gpuMemcpyDeviceToHost );
}






// extern "C" void test_extrap_stokes_all_models_gpu(source_t *chunked_source,
//            int num_extrap_freqs, double *extrap_freqs,
//            user_precision_t *extrap_flux_I, user_precision_t *extrap_flux_Q,
//            user_precision_t *extrap_flux_U, user_precision_t *extrap_flux_V){

//   int do_gpu = 1;

//   source_t *d_chunked_source = copy_chunked_source_to_GPU(chunked_source);

//   double *d_extrap_freqs = NULL;
//   gpuMalloc( (void**)&d_extrap_freqs,
//                                    num_extrap_freqs*sizeof(double) );
//   gpuMemcpy(d_extrap_freqs, extrap_freqs,
//              num_extrap_freqs*sizeof(double), gpuMemcpyHostToDevice );

//   malloc_extrapolated_flux_arrays_gpu(&d_chunked_source->point_components,
//                                   d_chunked_source->n_points,
//                                   num_extrap_freqs);

//   extrapolate_Stokes(d_chunked_source, d_extrap_freqs, num_extrap_freqs, POINT,
//                      do_gpu);


//   components_t d_components = d_chunked_source->point_components;

//   gpuMemcpy(extrap_flux_I, d_components.extrap_stokesI,
//             d_chunked_source->n_points*num_extrap_freqs*sizeof(user_precision_t),
//                                                       gpuMemcpyDeviceToHost );
//   gpuMemcpy(extrap_flux_Q, d_components.extrap_stokesQ,
//             d_chunked_source->n_points*num_extrap_freqs*sizeof(user_precision_t),
//                                                       gpuMemcpyDeviceToHost );
//   gpuMemcpy(extrap_flux_U, d_components.extrap_stokesU,
//             d_chunked_source->n_points*num_extrap_freqs*sizeof(user_precision_t),
//                                                       gpuMemcpyDeviceToHost );
//   gpuMemcpy(extrap_flux_V, d_components.extrap_stokesV,
//             d_chunked_source->n_points*num_extrap_freqs*sizeof(user_precision_t),
//                                                       gpuMemcpyDeviceToHost );
//   gpuFree( d_extrap_freqs );
//   free_extrapolated_flux_arrays(&d_chunked_source->point_components);
// }