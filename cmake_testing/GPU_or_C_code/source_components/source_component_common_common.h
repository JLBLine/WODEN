#include "constants.h"
#include "woden_precision_defs.h"
#include "woden_struct_defs.h"
#include "common_testing_functions.h"
#include "source_components_common.h"
// #include "source_components_cpu.h"
#include <string.h>
#include <mwa_hyperbeam.h>
#include "hyperbeam_error.h"


void test_source_component_common_ConstantDecFEEBeam(e_component_type comptype,
                                                     int do_gpu);

void test_source_component_common_ConstantDecFEEBeamInterp(e_component_type comptype,
                                                           int do_gpu);

void test_source_component_common_ConstantDecAnalyBeam(e_component_type comptype,
                                                       int do_gpu);
                                                       
void test_source_component_common_ConstantDecGaussBeam(e_component_type comptype,
                                                       int do_gpu);
                                                       
void test_source_component_common_ConstantDecNoBeam(e_component_type comptype,
                                                       int do_gpu);
                                                       
void test_source_component_common_ConstantDecMWAAnaly(e_component_type comptype,
                                                       int do_gpu);