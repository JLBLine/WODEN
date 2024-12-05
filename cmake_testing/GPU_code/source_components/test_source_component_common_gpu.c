#include "source_component_common_common.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

void test_source_component_common_ConstantDecFEEBeamPoint_gpu(void) {
  test_source_component_common_ConstantDecFEEBeam(POINT, 1);
}
void test_source_component_common_ConstantDecFEEBeamGauss_gpu(void) {
  test_source_component_common_ConstantDecFEEBeam(GAUSSIAN, 1);
}
void test_source_component_common_ConstantDecFEEBeamShapelet_gpu(void) {
  test_source_component_common_ConstantDecFEEBeam(SHAPELET, 1);
}


void test_source_component_common_ConstantDecAnalyBeamPoint_gpu(void){
  test_source_component_common_ConstantDecAnalyBeam(POINT, 1);
}
void test_source_component_common_ConstantDecAnalyBeamGaussian_gpu(void){
  test_source_component_common_ConstantDecAnalyBeam(GAUSSIAN, 1);
}
void test_source_component_common_ConstantDecAnalyBeamShapelet_gpu(void){
  test_source_component_common_ConstantDecAnalyBeam(SHAPELET, 1);
}


void test_source_component_common_ConstantDecGaussBeamPoint_gpu(void){
  test_source_component_common_ConstantDecGaussBeam(POINT, 1);
}
void test_source_component_common_ConstantDecGaussBeamGaussian_gpu(void){
  test_source_component_common_ConstantDecGaussBeam(GAUSSIAN, 1);
}
void test_source_component_common_ConstantDecGaussBeamShapelet_gpu(void){
  test_source_component_common_ConstantDecGaussBeam(SHAPELET, 1);
}


void test_source_component_common_ConstantDecNoBeamPoint_gpu(void){
  test_source_component_common_ConstantDecNoBeam(POINT, 1);
}
void test_source_component_common_ConstantDecNoBeamGaussian_gpu(void){
  test_source_component_common_ConstantDecNoBeam(GAUSSIAN, 1);
}
void test_source_component_common_ConstantDecNoBeamShapelet_gpu(void){
  test_source_component_common_ConstantDecNoBeam(SHAPELET, 1);
}


void test_source_component_common_ConstantDecFEEBeamInterpPoint_gpu(void){
  test_source_component_common_ConstantDecFEEBeamInterp(POINT, 1);
}
void test_source_component_common_ConstantDecFEEBeamInterpGaussian_gpu(void){
  test_source_component_common_ConstantDecFEEBeamInterp(GAUSSIAN, 1);
}
void test_source_component_common_ConstantDecFEEBeamInterpShapelet_gpu(void){
  test_source_component_common_ConstantDecFEEBeamInterp(SHAPELET, 1);
}

void test_source_component_common_ConstantDecMWAAnalyPoint_gpu(void){
  test_source_component_common_ConstantDecMWAAnaly(POINT, 1);
}
void test_source_component_common_ConstantDecMWAAnalyGaussian_gpu(void){
  test_source_component_common_ConstantDecMWAAnaly(GAUSSIAN, 1);
}
void test_source_component_common_ConstantDecMWAAnalyShapelet_gpu(void){
  test_source_component_common_ConstantDecMWAAnaly(SHAPELET, 1);
}

//Run the test with unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_source_component_common_ConstantDecNoBeamPoint_gpu);
    RUN_TEST(test_source_component_common_ConstantDecNoBeamGaussian_gpu);
    RUN_TEST(test_source_component_common_ConstantDecNoBeamShapelet_gpu);

    RUN_TEST(test_source_component_common_ConstantDecFEEBeamPoint_gpu);
    RUN_TEST(test_source_component_common_ConstantDecFEEBeamGauss_gpu);
    RUN_TEST(test_source_component_common_ConstantDecFEEBeamShapelet_gpu);

    RUN_TEST(test_source_component_common_ConstantDecAnalyBeamPoint_gpu);
    RUN_TEST(test_source_component_common_ConstantDecAnalyBeamGaussian_gpu);
    RUN_TEST(test_source_component_common_ConstantDecAnalyBeamShapelet_gpu);

    RUN_TEST(test_source_component_common_ConstantDecGaussBeamPoint_gpu);
    RUN_TEST(test_source_component_common_ConstantDecGaussBeamGaussian_gpu);
    RUN_TEST(test_source_component_common_ConstantDecGaussBeamShapelet_gpu);

    RUN_TEST(test_source_component_common_ConstantDecFEEBeamInterpPoint_gpu);
    RUN_TEST(test_source_component_common_ConstantDecFEEBeamInterpGaussian_gpu);
    RUN_TEST(test_source_component_common_ConstantDecFEEBeamInterpShapelet_gpu);
 
    RUN_TEST(test_source_component_common_ConstantDecMWAAnalyPoint_gpu);
    RUN_TEST(test_source_component_common_ConstantDecMWAAnalyGaussian_gpu);
    RUN_TEST(test_source_component_common_ConstantDecMWAAnalyShapelet_gpu);

    return UNITY_END();
}
