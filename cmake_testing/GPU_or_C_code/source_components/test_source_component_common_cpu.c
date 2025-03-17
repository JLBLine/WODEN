#include "source_component_common_common.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

void test_source_component_common_ConstantDecFEEBeamPoint_cpu(void) {
  test_source_component_common_ConstantDecFEEBeam(POINT, 0);
}
void test_source_component_common_ConstantDecFEEBeamGauss_cpu(void) {
  test_source_component_common_ConstantDecFEEBeam(GAUSSIAN, 0);
}
void test_source_component_common_ConstantDecFEEBeamShapelet_cpu(void) {
  test_source_component_common_ConstantDecFEEBeam(SHAPELET, 0);
}


void test_source_component_common_ConstantDecAnalyBeamPoint_cpu(void){
  test_source_component_common_ConstantDecAnalyBeam(POINT, 0);
}
void test_source_component_common_ConstantDecAnalyBeamGaussian_cpu(void){
  test_source_component_common_ConstantDecAnalyBeam(GAUSSIAN, 0);
}
void test_source_component_common_ConstantDecAnalyBeamShapelet_cpu(void){
  test_source_component_common_ConstantDecAnalyBeam(SHAPELET, 0);
}


void test_source_component_common_ConstantDecGaussBeamPoint_cpu(void){
  test_source_component_common_ConstantDecGaussBeam(POINT, 0);
}
void test_source_component_common_ConstantDecGaussBeamGaussian_cpu(void){
  test_source_component_common_ConstantDecGaussBeam(GAUSSIAN, 0);
}
void test_source_component_common_ConstantDecGaussBeamShapelet_cpu(void){
  test_source_component_common_ConstantDecGaussBeam(SHAPELET, 0);
}


void test_source_component_common_ConstantDecNoBeamPoint_cpu(void){
  test_source_component_common_ConstantDecNoBeam(POINT, 0);
}
void test_source_component_common_ConstantDecNoBeamGaussian_cpu(void){
  test_source_component_common_ConstantDecNoBeam(GAUSSIAN, 0);
}
void test_source_component_common_ConstantDecNoBeamShapelet_cpu(void){
  test_source_component_common_ConstantDecNoBeam(SHAPELET, 0);
}


void test_source_component_common_ConstantDecFEEBeamInterpPoint_cpu(void){
  test_source_component_common_ConstantDecFEEBeamInterp(POINT, 0);
}
void test_source_component_common_ConstantDecFEEBeamInterpGaussian_cpu(void){
  test_source_component_common_ConstantDecFEEBeamInterp(GAUSSIAN, 0);
}
void test_source_component_common_ConstantDecFEEBeamInterpShapelet_cpu(void){
  test_source_component_common_ConstantDecFEEBeamInterp(SHAPELET, 0);
}

void test_source_component_common_ConstantDecMWAAnalyPoint_cpu(void){
  test_source_component_common_ConstantDecMWAAnaly(POINT, 0);
}
void test_source_component_common_ConstantDecMWAAnalyGaussian_cpu(void){
  test_source_component_common_ConstantDecMWAAnaly(GAUSSIAN, 0);
}
void test_source_component_common_ConstantDecMWAAnalyShapelet_cpu(void){
  test_source_component_common_ConstantDecMWAAnaly(SHAPELET, 0);
}

void test_source_component_common_ConstantDecEveryBeamMWAPoint_cpu(void){
  test_source_component_common_ConstantDecEveryBeamMWA(POINT, 0);
}
void test_source_component_common_ConstantDecEveryBeamMWAGaussian_cpu(void){
  test_source_component_common_ConstantDecEveryBeamMWA(GAUSSIAN, 0);
}
void test_source_component_common_ConstantDecEveryBeamMWAShapelet_cpu(void){
  test_source_component_common_ConstantDecEveryBeamMWA(SHAPELET, 0);
}

void test_source_component_common_ConstantDecEveryBeamLOFARPoint_cpu(void){
  test_source_component_common_ConstantDecEveryBeamLOFAR(POINT, 0);
}
void test_source_component_common_ConstantDecEveryBeamLOFARGaussian_cpu(void){
  test_source_component_common_ConstantDecEveryBeamLOFAR(GAUSSIAN, 0);
}
void test_source_component_common_ConstantDecEveryBeamLOFARShapelet_cpu(void){
  test_source_component_common_ConstantDecEveryBeamLOFAR(SHAPELET, 0);
}

void test_source_component_common_ConstantDecEveryBeamOSKARPoint_cpu(void){
  test_source_component_common_ConstantDecEveryBeamOSKAR(POINT, 0);
}
void test_source_component_common_ConstantDecEveryBeamOSKARGaussian_cpu(void){
  test_source_component_common_ConstantDecEveryBeamOSKAR(GAUSSIAN, 0);
}
void test_source_component_common_ConstantDecEveryBeamOSKARShapelet_cpu(void){
  test_source_component_common_ConstantDecEveryBeamOSKAR(SHAPELET, 0);
}

//Run the test with unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_source_component_common_ConstantDecNoBeamPoint_cpu);
    RUN_TEST(test_source_component_common_ConstantDecNoBeamGaussian_cpu);
    RUN_TEST(test_source_component_common_ConstantDecNoBeamShapelet_cpu);

    RUN_TEST(test_source_component_common_ConstantDecFEEBeamPoint_cpu);
    RUN_TEST(test_source_component_common_ConstantDecFEEBeamGauss_cpu);
    RUN_TEST(test_source_component_common_ConstantDecFEEBeamShapelet_cpu);

    RUN_TEST(test_source_component_common_ConstantDecAnalyBeamPoint_cpu);
    RUN_TEST(test_source_component_common_ConstantDecAnalyBeamGaussian_cpu);
    RUN_TEST(test_source_component_common_ConstantDecAnalyBeamShapelet_cpu);

    RUN_TEST(test_source_component_common_ConstantDecGaussBeamPoint_cpu);
    RUN_TEST(test_source_component_common_ConstantDecGaussBeamGaussian_cpu);
    RUN_TEST(test_source_component_common_ConstantDecGaussBeamShapelet_cpu);

    RUN_TEST(test_source_component_common_ConstantDecFEEBeamInterpPoint_cpu);
    RUN_TEST(test_source_component_common_ConstantDecFEEBeamInterpGaussian_cpu);
    RUN_TEST(test_source_component_common_ConstantDecFEEBeamInterpShapelet_cpu);
 
    RUN_TEST(test_source_component_common_ConstantDecMWAAnalyPoint_cpu);
    RUN_TEST(test_source_component_common_ConstantDecMWAAnalyGaussian_cpu);
    RUN_TEST(test_source_component_common_ConstantDecMWAAnalyShapelet_cpu);

    #if defined(HAVE_EVERYBEAM)

      RUN_TEST(test_source_component_common_ConstantDecEveryBeamMWAPoint_cpu);
      RUN_TEST(test_source_component_common_ConstantDecEveryBeamMWAGaussian_cpu);
      RUN_TEST(test_source_component_common_ConstantDecEveryBeamMWAShapelet_cpu);

      RUN_TEST(test_source_component_common_ConstantDecEveryBeamLOFARPoint_cpu);
      RUN_TEST(test_source_component_common_ConstantDecEveryBeamLOFARGaussian_cpu);
      RUN_TEST(test_source_component_common_ConstantDecEveryBeamLOFARShapelet_cpu);

      RUN_TEST(test_source_component_common_ConstantDecEveryBeamOSKARPoint_cpu);
      RUN_TEST(test_source_component_common_ConstantDecEveryBeamOSKARGaussian_cpu);
      RUN_TEST(test_source_component_common_ConstantDecEveryBeamOSKARShapelet_cpu);
    #endif


    return UNITY_END();
}
