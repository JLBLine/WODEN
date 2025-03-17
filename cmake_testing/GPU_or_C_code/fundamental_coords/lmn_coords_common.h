#pragma once
#include <unity.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "woden_precision_defs.h"
#include "woden_struct_defs.h"
#include "fundamental_coords_cpu.h"

//GPU code we're linking in
extern void test_calc_lmn_for_components_gpu(double *ls, double *ms, double *ns,
                                         double *ras, double *decs,
                                         int num_coords,
                                         woden_settings_t *woden_settings);


void test_calc_lmn_GivesCorrectlCoords(int do_gpu);

void test_calc_lmn_GivesCorrectmCoords(int do_gpu);

