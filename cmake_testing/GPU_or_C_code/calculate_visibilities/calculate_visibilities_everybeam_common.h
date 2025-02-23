#include <math.h>
#include <unity.h>
#include <stdlib.h>
#include <complex.h>

#include "calculate_visibilities_common_common.h"

void test_calculate_visibilities_EveryBeam_OneSource_SinglePoint(int do_gpu, int beamtype, 
                                                                  const char *beam_ms_path);
void test_calculate_visibilities_EveryBeam_OneSource_SingleGauss(int do_gpu, int beamtype, 
                                                                  const char *beam_ms_path);
void test_calculate_visibilities_EveryBeam_OneSource_SingleShape(int do_gpu, int beamtype, 
                                                                  const char *beam_ms_path);
void test_calculate_visibilities_EveryBeam_OneSource_SingleAll(int do_gpu, int beamtype, 
                                                                  const char *beam_ms_path);
void test_calculate_visibilities_EveryBeam_ThreeSource_SinglePoint(int do_gpu, int beamtype, 
                                                                  const char *beam_ms_path);
void test_calculate_visibilities_EveryBeam_ThreeSource_SingleGauss(int do_gpu, int beamtype, 
                                                                  const char *beam_ms_path);
void test_calculate_visibilities_EveryBeam_ThreeSource_SingleShape(int do_gpu, int beamtype, 
                                                                  const char *beam_ms_path);
void test_calculate_visibilities_EveryBeam_ThreeSource_SingleAll(int do_gpu, int beamtype, 
                                                                  const char *beam_ms_path);
void test_calculate_visibilities_EveryBeam_ThreeSource_FivePoint(int do_gpu, int beamtype, 
                                                                  const char *beam_ms_path);
void test_calculate_visibilities_EveryBeam_ThreeSource_FiveGauss(int do_gpu, int beamtype, 
                                                                  const char *beam_ms_path);
void test_calculate_visibilities_EveryBeam_ThreeSource_FiveShape(int do_gpu, int beamtype, 
                                                                  const char *beam_ms_path);
void test_calculate_visibilities_EveryBeam_ThreeSource_FiveAll(int do_gpu, int beamtype, 
                                                                  const char *beam_ms_path);