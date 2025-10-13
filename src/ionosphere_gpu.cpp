#include "ionosphere_gpu.h"

__device__ double calc_ionospheric_phase_offset_gpu(double *d_ant_X,
           double *d_ant_Y, double *d_ant_Z,
           user_precision_t *d_azs, user_precision_t *d_zas,
           int *ant1_to_baseline_map, int *ant2_to_baseline_map,
           int num_baselines, int num_ants, int time_ind, int num_components,
           const int iBaseline, const int iComponent) {

    double height = 100000;
    int baseline_ind = iBaseline % num_baselines;
    int ant1 = time_ind*num_ants + ant1_to_baseline_map[baseline_ind];
    int ant2 = time_ind*num_ants + ant2_to_baseline_map[baseline_ind];

    double ant1_X = d_ant_X[ant1];
    double ant1_Y = d_ant_Y[ant1];
    double ant1_Z = d_ant_Z[ant1];
    double ant2_X = d_ant_X[ant2];
    double ant2_Y = d_ant_Y[ant2];
    double ant2_Z = d_ant_Z[ant2];
    double az = (double)d_azs[time_ind*num_components + iComponent];
    double zen = (double)d_zas[time_ind*num_components + iComponent];

    // find pierce points
    double pp1_x = ant1_X + height * tan(zen) * sin(az);
    double pp1_y = ant1_Y + height * tan(zen) * cos(az);
    double pp2_x = ant2_X + height * tan(zen) * sin(az);
    double pp2_y = ant2_Y + height * tan(zen) * cos(az);

    return get_phase_delay_gpu(pp1_x, pp1_y) - get_phase_delay_gpu(pp2_x, pp2_y);
}

__device__ double get_phase_delay_gpu(double pp_x, double pp_y) {
    return 1 * sin(1 + pp_x * 0.003);
}