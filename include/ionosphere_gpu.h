#pragma once
#include "woden_precision_defs.h"

__device__ double calc_ionospheric_phase_offset_gpu(double *d_ant_east,
           double *d_ant_north, double *d_ant_height,
           user_precision_t *d_azs, user_precision_t *d_zas,
           user_precision_t *d_allsteps_wavelengths,
           int *ant1_to_baseline_map, int *ant2_to_baseline_map,
           int num_baselines, int num_ants, int time_ind, int num_components,
           const int iBaseline, const int iComponent,
           user_precision_t TEC_grad_east, user_precision_t TEC_grad_north,
           user_precision_t *TEC_screen, int resolution,
           user_precision_t screen_size, user_precision_t height);

__device__ double get_phase_delay_gpu(double pp_east, double pp_north,
           double TEC_grad_east, double TEC_grad_north,
           double wavelength, double pp_zer);

__device__ user_precision_t get_phase_delay_from_TEC_gpu(user_precision_t pp_east, user_precision_t pp_north,
           user_precision_t *TEC_screen, int resolution,
           user_precision_t screen_size, user_precision_t wavelength, user_precision_t pp_zer);