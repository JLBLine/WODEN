#pragma once
#include "woden_precision_defs.h"

double calc_ionospheric_phase_offset_cpu(double ant1_east, double ant1_norht, double ant1_height,
                            double ant2_east, double ant2_north, double ant2_height,
                            user_precision_t az, user_precision_t zen,
                            user_precision_t wavelength,
                            user_precision_t TEC_grad_east, user_precision_t TEC_grad_north,
                            user_precision_t *TEC_screen, int resolution,
                            user_precision_t screen_size, user_precision_t height);

double get_phase_delay_cpu(double pp_east, double pp_north,
                        double TEC_grad_east, double TEC_grad_north,
                        double wavelength);

double get_phase_delay_from_TEC_cpu(double pp_east, double pp_north,
           user_precision_t *TEC_screen, int resolution,
           user_precision_t screen_size, double wavelength);