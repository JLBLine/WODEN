#include "read_and_write.h"
#include "constants.h"
// #include "create_sky_model.h"

void calc_para_angle(catsource_t *cropped_src, float *lsts,
                     int num_time_steps);

beam_settings_t fill_primary_beam_settings(woden_settings_t *woden_settings,
                MetaFfile_t metafits, catsource_t *cropped_src,
                float *lsts, int num_time_steps);

void setup_FEE_beam(woden_settings_t *woden_settings, MetaFfile_t metafits,
                    beam_settings_t beam_settings, float base_middle_freq);
