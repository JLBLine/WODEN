#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "woden_precision_defs.h"
#include "woden_struct_defs.h"
#include "source_components_gpu.h"
#include "source_components_common.h"

extern "C" double * malloc_freqs_gpu(int num_extrap_freqs, double *extrap_freqs);

extern "C" void free_freqs_gpu(double *d_extrap_freqs);

extern "C" void copy_extrapolated_flux_arrays_to_host(source_t *d_chunked_source,
                                                  int num_extrap_freqs,
                                                  user_precision_t *extrap_flux_I,
                                                  user_precision_t *extrap_flux_Q,
                                                  user_precision_t *extrap_flux_U,
                                                  user_precision_t *extrap_flux_V);