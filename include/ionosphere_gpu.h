#pragma once
#include <math.h>
#include "woden_precision_defs.h"

__device__ double calc_ionospheric_phase_offset_gpu(double *d_ant_X,
           double *d_ant_Y, double *d_ant_Z,
           user_precision_t *d_azs, user_precision_t *d_zas,
           int *ant1_to_baseline_map, int *ant2_to_baseline_map,
           int num_baselines, int num_ants, int time_ind, int num_components,
           const int iBaseline, const int iComponent);

__device__ double get_phase_delay_gpu(double pp_x, double pp_y);