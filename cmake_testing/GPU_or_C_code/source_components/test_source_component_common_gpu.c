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

void test_source_component_common_ConstantDecEveryBeamMWAPoint_gpu(void){
  test_source_component_common_ConstantDecEveryBeamMWA(POINT, 1);
}
void test_source_component_common_ConstantDecEveryBeamMWAGaussian_gpu(void){
  test_source_component_common_ConstantDecEveryBeamMWA(GAUSSIAN, 1);
}
void test_source_component_common_ConstantDecEveryBeamMWAShapelet_gpu(void){
  test_source_component_common_ConstantDecEveryBeamMWA(SHAPELET, 1);
}

void test_source_component_common_ConstantDecEveryBeamLOFARPoint_gpu(void){
  test_source_component_common_ConstantDecEveryBeamLOFAR(POINT, 1);
}
void test_source_component_common_ConstantDecEveryBeamLOFARGaussian_gpu(void){
  test_source_component_common_ConstantDecEveryBeamLOFAR(GAUSSIAN, 1);
}
void test_source_component_common_ConstantDecEveryBeamLOFARShapelet_gpu(void){
  test_source_component_common_ConstantDecEveryBeamLOFAR(SHAPELET, 1);
}

void test_source_component_common_ConstantDecEveryBeamOSKARPoint_gpu(void){
  test_source_component_common_ConstantDecEveryBeamOSKAR(POINT, 1);
}
void test_source_component_common_ConstantDecEveryBeamOSKARGaussian_gpu(void){
  test_source_component_common_ConstantDecEveryBeamOSKAR(GAUSSIAN, 1);
}
void test_source_component_common_ConstantDecEveryBeamOSKARShapelet_gpu(void){
  test_source_component_common_ConstantDecEveryBeamOSKAR(SHAPELET, 1);
}

void test_source_component_common_ConstantDecUVBeamMWAPoint_gpu(void){
  test_source_component_common_ConstantDecUVBeamMWA(POINT, 1);
}
void test_source_component_common_ConstantDecUVBeamMWAGaussian_gpu(void){
  test_source_component_common_ConstantDecUVBeamMWA(GAUSSIAN, 1);
}
void test_source_component_common_ConstantDecUVBeamMWAShapelet_gpu(void){
  test_source_component_common_ConstantDecUVBeamMWA(SHAPELET, 1);
}

void test_source_component_common_ConstantDecUVBeamHERAPoint_gpu(void){
  test_source_component_common_ConstantDecUVBeamHERA(POINT, 1);
}
void test_source_component_common_ConstantDecUVBeamHERAGaussian_gpu(void){
  test_source_component_common_ConstantDecUVBeamHERA(GAUSSIAN, 1);
}
void test_source_component_common_ConstantDecUVBeamHERAShapelet_gpu(void){
  test_source_component_common_ConstantDecUVBeamHERA(SHAPELET, 1);
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

    #if defined(USE_EVERYBEAM)
      RUN_TEST(test_source_component_common_ConstantDecEveryBeamMWAPoint_gpu);
      RUN_TEST(test_source_component_common_ConstantDecEveryBeamMWAGaussian_gpu);
      RUN_TEST(test_source_component_common_ConstantDecEveryBeamMWAShapelet_gpu);

      RUN_TEST(test_source_component_common_ConstantDecEveryBeamLOFARPoint_gpu);
      RUN_TEST(test_source_component_common_ConstantDecEveryBeamLOFARGaussian_gpu);
      RUN_TEST(test_source_component_common_ConstantDecEveryBeamLOFARShapelet_gpu);

      RUN_TEST(test_source_component_common_ConstantDecEveryBeamOSKARPoint_gpu);
      RUN_TEST(test_source_component_common_ConstantDecEveryBeamOSKARGaussian_gpu);
      RUN_TEST(test_source_component_common_ConstantDecEveryBeamOSKARShapelet_gpu);
    #endif

    RUN_TEST(test_source_component_common_ConstantDecUVBeamMWAPoint_gpu);
    RUN_TEST(test_source_component_common_ConstantDecUVBeamMWAGaussian_gpu);
    RUN_TEST(test_source_component_common_ConstantDecUVBeamMWAShapelet_gpu);

    RUN_TEST(test_source_component_common_ConstantDecUVBeamHERAPoint_gpu);
    RUN_TEST(test_source_component_common_ConstantDecUVBeamHERAGaussian_gpu);
    RUN_TEST(test_source_component_common_ConstantDecUVBeamHERAShapelet_gpu);

    return UNITY_END();
}
