#include <math.h>
#include "read_and_write.h"

extern "C" void calculate_visibilities(float *X_diff_metres, float *Y_diff_metres, float *Z_diff_metres,
                    catsource_t catsource, float *angles_array, beam_settings_t beam_settings,
                    const int num_baselines, const int num_time_steps, const int num_visis,
                    const int num_freqs,
                    visibility_set_t *visibility_set, float *sbf);
