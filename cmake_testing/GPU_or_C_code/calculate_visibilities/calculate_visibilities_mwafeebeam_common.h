#include <math.h>
#include <unity.h>
#include <stdlib.h>
#include <complex.h>

#include "calculate_visibilities_common_common.h"


void test_calculate_visibilities_MWAFEEBeam_OneSource_SinglePoint(int do_gpu);
void test_calculate_visibilities_MWAFEEBeam_OneSource_SingleGauss(int do_gpu);
void test_calculate_visibilities_MWAFEEBeam_OneSource_SingleShape(int do_gpu);
void test_calculate_visibilities_MWAFEEBeam_OneSource_SingleAll(int do_gpu);
void test_calculate_visibilities_MWAFEEBeam_ThreeSource_SinglePoint(int do_gpu);
void test_calculate_visibilities_MWAFEEBeam_ThreeSource_SingleGauss(int do_gpu);
void test_calculate_visibilities_MWAFEEBeam_ThreeSource_SingleShape(int do_gpu);
void test_calculate_visibilities_MWAFEEBeam_ThreeSource_SingleAll(int do_gpu);
void test_calculate_visibilities_MWAFEEBeam_ThreeSource_ThreePoint(int do_gpu);
void test_calculate_visibilities_MWAFEEBeam_ThreeSource_ThreeGauss(int do_gpu);
void test_calculate_visibilities_MWAFEEBeam_ThreeSource_ThreeShape(int do_gpu);
void test_calculate_visibilities_MWAFEEBeam_ThreeSource_ThreeAll(int do_gpu);

void test_calculate_visibilities_MWAFEEBeam(int n_points, int n_gauss, int n_shapes,
                                           int num_sources, int do_gpu);