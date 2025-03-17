#include <math.h>
#include <unity.h>
#include <stdlib.h>
#include <stdint.h>
#include <complex.h>

#include "constants.h"
#include "woden_struct_defs.h"
#include "hyperbeam_error.h"
#include "primary_beam_cpu.h"
#include "run_hyperbeam_common.h"

void test_hyperbeam_100MHz_zenith_multiant(int do_gpu);
void test_hyperbeam_150MHz_zenith_multiant(int do_gpu);
void test_hyperbeam_200MHz_zenith_multiant(int do_gpu);
void test_hyperbeam_100MHz_off_zenith1_multiant(int do_gpu);
void test_hyperbeam_150MHz_off_zenith1_multiant(int do_gpu);
void test_hyperbeam_200MHz_off_zenith1_multiant(int do_gpu);
void test_hyperbeam_100MHz_off_zenith2_multiant(int do_gpu);
void test_hyperbeam_150MHz_off_zenith2_multiant(int do_gpu);
void test_hyperbeam_200MHz_off_zenith2_multiant(int do_gpu);