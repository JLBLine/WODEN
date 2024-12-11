#include <math.h>
#include <unity.h>
#include <stdlib.h>
#include <stdint.h>
#include <complex.h>

#include "constants.h"
#include "woden_struct_defs.h"
#include "hyperbeam_error.h"
#include "primary_beam_cpu.h"
// #include "hyperbeam_common.h"

// //External CUDA code we're linking in
extern void test_run_hyperbeam_gpu(int num_components,
           int num_times, int num_freqs, int num_beams,
           uint8_t parallatic,
           struct FEEBeamGpu *gpu_fee_beam,
           double *azs, double *zas,
           double *latitudes,
           user_precision_complex_t *primay_beam_J00,
           user_precision_complex_t *primay_beam_J01,
           user_precision_complex_t *primay_beam_J10,
           user_precision_complex_t *primay_beam_J11);


void test_hyperbeam_100MHz_zenith(int do_gpu);
void test_hyperbeam_150MHz_zenith(int do_gpu);
void test_hyperbeam_200MHz_zenith(int do_gpu);
void test_hyperbeam_100MHz_off_zenith1(int do_gpu);
void test_hyperbeam_150MHz_off_zenith1(int do_gpu);
void test_hyperbeam_200MHz_off_zenith1(int do_gpu);
void test_hyperbeam_100MHz_off_zenith2(int do_gpu);
void test_hyperbeam_150MHz_off_zenith2(int do_gpu);
void test_hyperbeam_200MHz_off_zenith2(int do_gpu);