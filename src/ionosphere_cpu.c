#include "ionosphere_cpu.h"
#include <math.h>
#include "constants.h"

double calc_ionospheric_phase_offset_cpu(double ant1_east, double ant1_north, double ant1_height,
                            double ant2_east, double ant2_north, double ant2_height,
                            user_precision_t az, user_precision_t zen,
                            user_precision_t wavelength,
                            user_precision_t TEC_grad_east, user_precision_t TEC_grad_north,
                            user_precision_t *TEC_screen, int resolution,
                            user_precision_t screen_size, user_precision_t height) {

    // find pierce points
    double pp1_east = ant1_east + (height - ant1_height) * tan(zen) * sin(az);
    double pp1_north = ant1_north + (height - ant1_height) * tan(zen) * cos(az);
    double pp2_east = ant2_east + (height - ant2_height) * tan(zen) * sin(az);
    double pp2_north = ant2_north + (height - ant2_height) * tan(zen) * cos(az);

    // double phase1 = get_phase_delay_cpu(pp1_east, pp1_north, (double)TEC_grad_east, (double)TEC_grad_north, (double)wavelength);
    // double phase2 = get_phase_delay_cpu(pp2_east, pp2_north, (double)TEC_grad_east, (double)TEC_grad_north, (double)wavelength);

    double phase1 = get_phase_delay_from_TEC_cpu(pp1_east, pp1_north, TEC_screen, resolution, screen_size, (double)wavelength);
    double phase2 = get_phase_delay_from_TEC_cpu(pp2_east, pp2_north, TEC_screen, resolution, screen_size, (double)wavelength);

    return phase1 - phase2;
}

double get_phase_delay_cpu(double pp_east, double pp_north,
                        double TEC_grad_east, double TEC_grad_north,
                        double wavelength) {
    double TEC = TEC_grad_east * pp_east + TEC_grad_north * pp_north;

    return TEC * wavelength * TEC_TO_PHASE;
}

double get_phase_delay_from_TEC_cpu(double pp_east, double pp_north,
           user_precision_t *TEC_screen, int resolution,
           user_precision_t screen_size, double wavelength) {
    
    if (-screen_size * 0.5 > pp_east  || pp_east > screen_size * 0.5 || 
        -screen_size * 0.5 > pp_north  || pp_north > screen_size * 0.5) {
        // oops TEC too small
        return 0;
    }

    // getting point in pixel coordinates
    double x = (pp_east / screen_size + 0.5) * resolution;
    double y = (pp_north / screen_size + 0.5) * resolution;
    
    // interpolate here instead
    double TEC = TEC_screen[resolution * (int)y + (int)x];

    return TEC * wavelength * TEC_TO_PHASE;
}