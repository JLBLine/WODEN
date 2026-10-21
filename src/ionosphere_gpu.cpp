#include "ionosphere_gpu.h"
#include "constants.h"
#include <math.h>
#include <stdio.h>

__device__ double calc_ionospheric_phase_offset_gpu(double *d_ant_east,
           double *d_ant_norht, double *d_ant_height,
           user_precision_t *d_azs, user_precision_t *d_zas,
           user_precision_t *d_allsteps_wavelengths,
           int *ant1_to_baseline_map, int *ant2_to_baseline_map,
           int num_baselines, int num_ants, int time_ind, int num_components,
           const int iBaseline, const int iComponent,
           user_precision_t TEC_grad_east, user_precision_t TEC_grad_north,
           user_precision_t *TEC_screen, int resolution,
           user_precision_t screen_size, user_precision_t height) {

    int baseline_ind = iBaseline % num_baselines;
    int ant1 = time_ind*num_ants + ant1_to_baseline_map[baseline_ind];
    int ant2 = time_ind*num_ants + ant2_to_baseline_map[baseline_ind];

    user_precision_t wavelength = d_allsteps_wavelengths[iBaseline];

    user_precision_t ant1_east = d_ant_east[ant1];
    user_precision_t ant1_north = d_ant_norht[ant1];
    user_precision_t ant1_height = d_ant_height[ant1];
    user_precision_t ant2_east = d_ant_east[ant2];
    user_precision_t ant2_north = d_ant_norht[ant2];
    user_precision_t ant2_height = d_ant_height[ant2];
    user_precision_t az = d_azs[time_ind*num_components + iComponent];
    user_precision_t zen = d_zas[time_ind*num_components + iComponent];

    // find pierce points
    user_precision_t pp1_east = ant1_east + (height - ant1_height) * tan(zen) * sin(az);
    user_precision_t pp1_north = ant1_north + (height - ant1_height) * tan(zen) * cos(az);
    user_precision_t pp2_east = ant2_east + (height - ant2_height) * tan(zen) * sin(az);
    user_precision_t pp2_north = ant2_north + (height - ant2_height) * tan(zen) * cos(az);

    // for now just use MWA values
    user_precision_t longitude = 2.036289866851053;
    user_precision_t latitude = MWA_LAT_RAD;

    int do_spherical_TEC = 1;
    if (do_spherical_TEC) {
        // find point directly above antennas on TEC screen
        user_precision_t ant1_east_screen = (EARTH_RADIUS + height) / EARTH_RADIUS * ant1_east;
        user_precision_t ant1_north_screen = (EARTH_RADIUS + height) / EARTH_RADIUS * ant1_north;
        user_precision_t ant2_east_screen = (EARTH_RADIUS + height) / EARTH_RADIUS * ant2_east;
        user_precision_t ant2_north_screen = (EARTH_RADIUS + height) / EARTH_RADIUS * ant2_north;

        user_precision_t pp_zen = asin(EARTH_RADIUS / (EARTH_RADIUS + height) * sin(zen));
        pp1_east = ant1_east_screen + (EARTH_RADIUS + height) * (zen - pp_zen) * sin(az);
        pp1_north = ant1_north_screen + (EARTH_RADIUS + height) * (zen - pp_zen) * cos(az);
        pp2_east = ant2_east_screen + (EARTH_RADIUS + height) * (zen - pp_zen) * sin(az);
        pp2_north = ant2_north_screen + (EARTH_RADIUS + height) * (zen - pp_zen) * cos(az);
    }

    // double phase1 = get_phase_delay_gpu(pp1_east, pp1_north, (double)TEC_grad_east, (double)TEC_grad_north, wavelength);
    // double phase2 = get_phase_delay_gpu(pp2_east, pp2_north, (double)TEC_grad_east, (double)TEC_grad_north, wavelength);

    user_precision_t phase1 = get_phase_delay_from_TEC_gpu(pp1_east, pp1_north, TEC_screen, resolution, screen_size, wavelength);
    user_precision_t phase2 = get_phase_delay_from_TEC_gpu(pp2_east, pp2_north, TEC_screen, resolution, screen_size, wavelength);

    return phase1 - phase2;
}

__device__ double get_phase_delay_gpu(double pp_east, double pp_north,
           double TEC_grad_east, double TEC_grad_north,
           double wavelength) {
    double TEC = TEC_grad_east * pp_east + TEC_grad_north * pp_north;

    return TEC * wavelength * TEC_TO_PHASE;
}

__device__ user_precision_t get_phase_delay_from_TEC_gpu(user_precision_t pp_east, user_precision_t pp_north,
           user_precision_t *TEC_screen, int resolution,
           user_precision_t screen_size, user_precision_t wavelength) {
    
    user_precision_t longitude = 2.036289866851053;
    user_precision_t latitude = MWA_LAT_RAD;
    
    if (-screen_size * 0.5 > pp_east || pp_east > screen_size * 0.5 || 
        -screen_size * 0.5 > pp_north || pp_north > screen_size * 0.5) {
        // oops TEC too small
        return 0;
    }
    
    // getting point in pixel coordinates
    user_precision_t x = (pp_east / screen_size + 0.5) * resolution;
    user_precision_t y = (pp_north / screen_size + 0.5) * resolution;
    
    // interpolating
    user_precision_t TEC00 = TEC_screen[resolution * (int)x + (int)y];
    user_precision_t TEC10 = TEC_screen[resolution * ((int)x + 1) + (int)y];
    user_precision_t TEC01 = TEC_screen[resolution * (int)x + (int)y + 1];
    user_precision_t TEC11 = TEC_screen[resolution * ((int)x + 1) + (int)y + 1];

    // now just want the decimal part of x and y
    x = x - (int)x;
    y = y - (int)y;

    user_precision_t X0 = (1 - x) * TEC00 + x * TEC10;
    user_precision_t X1 = (1 - x) * TEC01 + x * TEC11;

    user_precision_t TEC = (1 - y) * X0 + y * X1;

    // user_precision_t TEC = TEC_screen[resolution * (int)round(x) + (int)round(y)];

    return TEC * wavelength * TEC_TO_PHASE;
}