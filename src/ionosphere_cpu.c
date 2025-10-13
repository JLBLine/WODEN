#include "ionosphere_cpu.h"

double calc_ionospheric_phase_offset_cpu(double ant1_X, double ant1_Y, double ant1_Z,
                            double ant2_X, double ant2_Y, double ant2_Z,
                            user_precision_t az, user_precision_t zen) {
    double height = 100000;

    // find pierce points
    double pp1_x = ant1_X + height * tan(zen) * sin(az);
    double pp1_y = ant1_Y + height * tan(zen) * cos(az);
    double pp2_x = ant2_X + height * tan(zen) * sin(az);
    double pp2_y = ant2_Y + height * tan(zen) * cos(az);

    return get_phase_delay_cpu(pp1_x, pp1_y) - get_phase_delay_cpu(pp2_x, pp2_y);
}

double get_phase_delay_cpu(double pp_x, double pp_y) {
    return 1 * sin(1 + pp_x * 0.003);
}