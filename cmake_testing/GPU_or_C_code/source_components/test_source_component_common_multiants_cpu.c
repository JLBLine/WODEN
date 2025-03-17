#include "source_component_common_multiants_common.h"
#include "hyperbeam_error.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */



void test_source_component_common_ConstantDecFEEBeamPoint_multiant_cpu(void) {
  test_source_component_common_ConstantDecFEEBeam_multiant(POINT, 0);
}
void test_source_component_common_ConstantDecFEEBeamGauss_multiant_cpu(void) {
  test_source_component_common_ConstantDecFEEBeam_multiant(GAUSSIAN, 0);
}
void test_source_component_common_ConstantDecFEEBeamShapelet_multiant_cpu(void) {
  test_source_component_common_ConstantDecFEEBeam_multiant(SHAPELET, 0);
}


void test_source_component_common_ConstantDecFEEBeamInterpPoint_multiant_cpu(void){
  test_source_component_common_ConstantDecFEEBeamInterp_multiant(POINT, 0);
}
void test_source_component_common_ConstantDecFEEBeamInterpGaussian_multiant_cpu(void){
  test_source_component_common_ConstantDecFEEBeamInterp_multiant(GAUSSIAN, 0);
}
void test_source_component_common_ConstantDecFEEBeamInterpShapelet_multiant_cpu(void){
  test_source_component_common_ConstantDecFEEBeamInterp_multiant(SHAPELET, 0);
}


//Run the test with unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_source_component_common_ConstantDecFEEBeamPoint_multiant_cpu);
    RUN_TEST(test_source_component_common_ConstantDecFEEBeamGauss_multiant_cpu);
    RUN_TEST(test_source_component_common_ConstantDecFEEBeamShapelet_multiant_cpu);
    
    RUN_TEST(test_source_component_common_ConstantDecFEEBeamInterpPoint_multiant_cpu);
    RUN_TEST(test_source_component_common_ConstantDecFEEBeamInterpGaussian_multiant_cpu);
    RUN_TEST(test_source_component_common_ConstantDecFEEBeamInterpShapelet_multiant_cpu);


    return UNITY_END();
}
