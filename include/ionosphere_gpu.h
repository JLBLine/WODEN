#pragma once
#include "woden_precision_defs.h"

__device__ double calc_ionospheric_phase_offset_gpu(double *d_ant_X,
           double *d_ant_Y, double *d_ant_Z,
           user_precision_t *d_azs, user_precision_t *d_zas,
           user_precision_t *d_allsteps_wavelengths,
           int *ant1_to_baseline_map, int *ant2_to_baseline_map,
           int num_baselines, int num_ants, int time_ind, int num_components,
           const int iBaseline, const int iComponent,
           user_precision_t TEC_grad_x, user_precision_t TEC_grad_y,
           user_precision_t *TEC_screen, int resolution,
           user_precision_t screen_size, user_precision_t height);

__device__ double get_phase_delay_gpu(double pp_x, double pp_y,
           double TEC_grad_x, double TEC_grad_y,
           double wavelength);

__device__ double get_phase_delay_from_TEC_gpu(double pp_x, double pp_y,
           user_precision_t *TEC_screen, int resolution,
           user_precision_t screen_size, double wavelength);