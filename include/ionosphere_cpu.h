#pragma once

#include <math.h>
#include "woden_precision_defs.h"

double calc_ionospheric_phase_offset(double ant1_X, double ant1_Y, double ant1_Z,
                            double ant2_X, double ant2_Y, double ant2_Z,
                            user_precision_t az, user_precision_t zen);

double get_phase_delay(double pp_x, double pp_y);