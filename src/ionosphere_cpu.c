#include "ionosphere_cpu.h"
#include <math.h>
#include "constants.h"

double calc_ionospheric_phase_offset_cpu(double ant1_X, double ant1_Y, double ant1_Z,
                            double ant2_X, double ant2_Y, double ant2_Z,
                            user_precision_t az, user_precision_t zen,
                            user_precision_t wavelength,
                            user_precision_t TEC_grad_x, user_precision_t TEC_grad_y,
                            user_precision_t *TEC_screen, int resolution,
                            user_precision_t screen_size, user_precision_t height) {

    // find pierce points
    double pp1_x = ant1_X + (height - ant1_Z) * tan(zen) * sin(az);
    double pp1_y = ant1_Y + (height - ant1_Z) * tan(zen) * cos(az);
    double pp2_x = ant2_X + (height - ant2_Z) * tan(zen) * sin(az);
    double pp2_y = ant2_Y + (height - ant2_Z) * tan(zen) * cos(az);

    // double phase1 = get_phase_delay_cpu(pp1_x, pp1_y, (double)TEC_grad_x, (double)TEC_grad_y, (double)wavelength);
    // double phase2 = get_phase_delay_cpu(pp2_x, pp2_y, (double)TEC_grad_x, (double)TEC_grad_y, (double)wavelength);

    double phase1 = get_phase_delay_from_TEC_cpu(pp1_x, pp1_y, TEC_screen, resolution, screen_size, (double)wavelength);
    double phase2 = get_phase_delay_from_TEC_cpu(pp2_x, pp2_y, TEC_screen, resolution, screen_size, (double)wavelength);

    return phase1 - phase2;
}

double get_phase_delay_cpu(double pp_x, double pp_y,
                        double TEC_grad_x, double TEC_grad_y,
                        double wavelength) {
    double TEC = TEC_grad_x * pp_x + TEC_grad_y * pp_y;

    return TEC * wavelength * TEC_TO_PHASE;
}

double get_phase_delay_from_TEC_cpu(double pp_x, double pp_y,
           user_precision_t *TEC_screen, int resolution,
           user_precision_t screen_size, double wavelength) {
    
    if (-screen_size * 0.5 > pp_x  || pp_x > screen_size * 0.5 || 
        -screen_size * 0.5 > pp_y  || pp_y > screen_size * 0.5) {
        // oops TEC too small
        return 0;
    }

    // getting point in pixel coordinates
    double x = (pp_x / screen_size + 0.5) * resolution;
    double y = (pp_y / screen_size + 0.5) * resolution;
    
    // interpolate here instead
    double TEC = TEC_screen[resolution * (int)x + (int)y];

    return TEC * wavelength * TEC_TO_PHASE;
}