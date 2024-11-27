#pragma once
#include <unity.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "woden_precision_defs.h"
#include "fundamental_coords_cpu.h"

void test_calc_lmn_GivesCorrectlCoords(int do_gpu);

void test_calc_lmn_GivesCorrectmCoords(int do_gpu);