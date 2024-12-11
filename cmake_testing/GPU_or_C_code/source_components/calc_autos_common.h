#include <math.h>
#include <unity.h>
#include <stdlib.h>
#include <complex.h>

#include "constants.h"
#include "woden_struct_defs.h"
#include "visibility_set.h"
#include "source_components_cpu.h"
#include "source_components_common.h"

//External CUDA code being linked in
extern void test_kern_calc_autos(components_t *components, int beamtype,
                                 int num_components, int num_baselines,
                                 int num_freqs, int num_times, int num_ants,
                                 int num_beams,
                                 visibility_set_t *visibility_set);

void test_calc_autos_cpu(components_t *components, int beamtype,
                                     int num_components, int num_baselines,
                                     int num_freqs, int num_times, int num_ants,
                                     int num_beams,
                                     visibility_set_t *visibility_set);

void test_calculate_autos(e_beamtype beamtype, int do_gpu);