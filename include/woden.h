#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "woden_precision_defs.h"
#include "constants.h"
#include "woden_struct_defs.h"
#include "beam_settings.h"
#include "visibility_set.h"
#include "hyperbeam_error.h"
#include "calculate_visibilities_common.h"
#include "logger.h"

int run_woden(woden_settings_t *woden_settings, visibility_set_t *visibility_sets,
             source_catalogue_t *cropped_sky_models, array_layout_t * array_layout,
             user_precision_t *sbf);