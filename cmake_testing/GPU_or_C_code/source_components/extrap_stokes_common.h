#include <unity.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "constants.h"
#include "woden_precision_defs.h"
#include "woden_struct_defs.h"
#include "common_testing_functions.h"
#include "source_components_common.h"
#include "source_components_cpu.h"

void test_extrap_stokes_GivesCorrectValues(int do_gpu);