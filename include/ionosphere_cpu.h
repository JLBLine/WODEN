#pragma once
#include "woden_precision_defs.h"

double calc_ionospheric_phase_offset_cpu(double ant1_X, double ant1_Y, double ant1_Z,
                            double ant2_X, double ant2_Y, double ant2_Z,
                            user_precision_t az, user_precision_t zen,
                            user_precision_t wavelength,
                            user_precision_t TEC_grad_x, user_precision_t TEC_grad_y);

double get_phase_delay_cpu(double pp_x, double pp_y,
                        double TEC_grad_x, double TEC_grad_y,
                        double wavelength);