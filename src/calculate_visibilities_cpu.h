#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

#include "woden_precision_defs.h"
#include "woden_struct_defs.h"
#include "constants.h"

// #include "calculate_visibilities_cpu.h"
// #include "fundamental_coords_cpu.h"

#include "source_components_common.h"
// #include "primary_beam_cpu.h"
// #include "hyperbeam_error.h"

//Helpful C code we are also using
#include "visibility_set.h"

calc_visi_inouts_t * create_calc_visi_inouts_cpu(array_layout_t *array_layout,
        visibility_set_t *visibility_set,  user_precision_t *sbf,
        woden_settings_t *woden_settings,
        int num_shapelets, int use_twobeams);

void set_visi_set_to_zero_cpu(visibility_set_t *visibility_set, int num_visis);

void free_calc_visi_inouts_cpu(calc_visi_inouts_t *calc_visi_inouts,
                               int num_shapelets, int use_twobeams);

void free_components_cpu(components_t *components);

void free_beam_gains_cpu(beam_gains_t *beam_gains, e_beamtype beamtype);