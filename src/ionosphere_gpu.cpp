#include "ionosphere_gpu.h"
#include "constants.h"
#include <math.h>

__device__ double calc_ionospheric_phase_offset_gpu(double *d_ant_X,
           double *d_ant_Y, double *d_ant_Z,
           user_precision_t *d_azs, user_precision_t *d_zas,
           user_precision_t *d_allsteps_wavelengths,
           int *ant1_to_baseline_map, int *ant2_to_baseline_map,
           int num_baselines, int num_ants, int time_ind, int num_components,
           const int iBaseline, const int iComponent,
           user_precision_t TEC_grad_x, user_precision_t TEC_grad_y,
           user_precision_t *TEC_screen, int resolution,
           user_precision_t screen_size, user_precision_t height) {

    int baseline_ind = iBaseline % num_baselines;
    int ant1 = time_ind*num_ants + ant1_to_baseline_map[baseline_ind];
    int ant2 = time_ind*num_ants + ant2_to_baseline_map[baseline_ind];

    user_precision_t wavelength = d_allsteps_wavelengths[iBaseline];

    user_precision_t ant1_X = d_ant_X[ant1];
    user_precision_t ant1_Y = d_ant_Y[ant1];
    user_precision_t ant1_Z = d_ant_Z[ant1];
    user_precision_t ant2_X = d_ant_X[ant2];
    user_precision_t ant2_Y = d_ant_Y[ant2];
    user_precision_t ant2_Z = d_ant_Z[ant2];
    user_precision_t az = d_azs[time_ind*num_components + iComponent];
    user_precision_t zen = d_zas[time_ind*num_components + iComponent];

    // find pierce points
    user_precision_t pp1_x = ant1_X + (height - ant1_Z) * tan(zen) * sin(az);
    user_precision_t pp1_y = ant1_Y + (height - ant1_Z) * tan(zen) * cos(az);
    user_precision_t pp2_x = ant2_X + (height - ant2_Z) * tan(zen) * sin(az);
    user_precision_t pp2_y = ant2_Y + (height - ant2_Z) * tan(zen) * cos(az);

    // for now just use MWA values
    user_precision_t longitude = 2.036289866851053;
    user_precision_t latitude = MWA_LAT_RAD;

    int do_spherical_TEC = 1;
    if (do_spherical_TEC) {
        pp1_x = zen * sin(az) + ant1_X / EARTH_RADIUS;
        pp1_y = zen * cos(az) + ant1_Y / EARTH_RADIUS;
        pp2_x = zen * sin(az) + ant2_X / EARTH_RADIUS;
        pp2_y = zen * cos(az) + ant2_Y / EARTH_RADIUS;
    }

    // double phase1 = get_phase_delay_gpu(pp1_x, pp1_y, (double)TEC_grad_x, (double)TEC_grad_y, wavelength);
    // double phase2 = get_phase_delay_gpu(pp2_x, pp2_y, (double)TEC_grad_x, (double)TEC_grad_y, wavelength);

    user_precision_t phase1 = get_phase_delay_from_TEC_gpu(pp1_x, pp1_y, TEC_screen, resolution, screen_size, wavelength);
    user_precision_t phase2 = get_phase_delay_from_TEC_gpu(pp2_x, pp2_y, TEC_screen, resolution, screen_size, wavelength);

    return phase1 - phase2;
}

__device__ double get_phase_delay_gpu(double pp_x, double pp_y,
           double TEC_grad_x, double TEC_grad_y,
           double wavelength) {
    double TEC = TEC_grad_x * pp_x + TEC_grad_y * pp_y;

    return TEC * wavelength * TEC_TO_PHASE;
}

__device__ user_precision_t get_phase_delay_from_TEC_gpu(user_precision_t pp_x, user_precision_t pp_y,
           user_precision_t *TEC_screen, int resolution,
           user_precision_t screen_size, user_precision_t wavelength) {
    
    user_precision_t longitude = 2.036289866851053;
    user_precision_t latitude = MWA_LAT_RAD;
    
    if (-screen_size * 0.5 > pp_x || pp_x > screen_size * 0.5 || 
        -screen_size * 0.5 > pp_y || pp_y > screen_size * 0.5) {
        // oops TEC too small
        return 0;
    }
    
    // getting point in pixel coordinates
    user_precision_t x = (pp_x / screen_size + 0.5) * resolution;
    user_precision_t y = (pp_y / screen_size + 0.5) * resolution;
    
    // interpolating
    user_precision_t TEC00 = TEC_screen[resolution * (int)x + (int)y];
    user_precision_t TEC10 = TEC_screen[resolution * (int)x + 1 + (int)y];
    user_precision_t TEC01 = TEC_screen[resolution * (int)x + (int)y + 1];
    user_precision_t TEC11 = TEC_screen[resolution * (int)x + 1 + (int)y + 1];

    // now just want the decimal part of x and y
    x = x - (int)x;
    y = y - (int)y;

    user_precision_t X0 = (1 - x) * TEC00 + x * TEC10;
    user_precision_t X1 = (1 - x) * TEC01 + x * TEC11;

    user_precision_t TEC = (1 - y) * X0 + y * X1;

    // user_precision_t TEC = TEC_screen[resolution * (int)round(x) + (int)round(y)];

    return TEC * wavelength * TEC_TO_PHASE;
}