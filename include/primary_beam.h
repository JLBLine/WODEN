#include "read_and_write.h"
#include "constants.h"
// #include "create_sky_model.h"

void calc_para_angle(catsource_t *cropped_src, float *lsts,
                     int num_time_steps);

beam_settings_t fill_primary_beam_settings(woden_settings_t *woden_settings,
                catsource_t *cropped_src,  float *lsts, int num_time_steps);
