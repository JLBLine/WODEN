#pragma once
#include <math.h>
#include <unity.h>
#include <stdlib.h>
#include <complex.h>

#include "constants.h"
#include "woden_struct_defs.h"
#include "woden_precision_defs.h"
#include "shapelet_basis.h"
#include "common_testing_functions.h"
#include "source_components_common.h"

//external GPU code used for testng

extern void test_kern_update_sum_visis(int num_freqs, int num_cross,
          int num_baselines, int num_components, int num_times, int beamtype,
          int use_twoants, int num_ants, int off_cardinal_dipoles,
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
          user_precision_t *sum_visi_YY_real, user_precision_t *sum_visi_YY_imag);

void test_update_sum_visis_VaryGainChooseBeams(int beamtype, int do_gpu);
void test_update_sum_visis_VaryFluxesChooseBeams(int beamtype, int do_gpu);
void test_update_sum_visis_VaryVisiChooseBeams(int beamtype, int do_gpu);