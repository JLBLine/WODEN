#include <math.h>
#include <unity.h>
#include <stdlib.h>
#include <complex.h>

#include "constants.h"
#include "woden_struct_defs.h"
#include "visibility_set.h"
#include "source_components_cpu.h"
#include "calc_autos_common.h"

void test_calculate_autos_multiants(e_beamtype beamtype, int do_gpu);