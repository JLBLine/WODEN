#pragma once
#include <math.h>
// #include "read_and_write.h"
#include "woden_struct_defs.h"

extern "C" void calculate_visibilities(array_layout_t * array_layout,
  source_catalogue_t *cropped_sky_models,
  woden_settings_t *woden_settings, visibility_set_t *visibility_set,
  visibility_set_t *chunk_visibility_set, float *sbf, int num_chunks);

// typedef struct _d_sum_visi_struct_t {
//   float *d_sum_visi_XX_real;
//   float *d_sum_visi_XX_imag;
//   float *d_sum_visi_XY_real;
//   float *d_sum_visi_XY_imag;
//   float *d_sum_visi_YX_real;
//   float *d_sum_visi_YX_imag;
//   float *d_sum_visi_YY_real;
//   float *d_sum_visi_YY_imag;
//
// } d_sum_visi_struct_t;
