#include "update_sum_visis_common.h"

/*
This test checks varying the gain with beamtype=FEE_BEAM
*/
void test_update_sum_visis_VaryGainFEEBeam_gpu(void) {
  test_update_sum_visis_VaryGainChooseBeams(FEE_BEAM, 1);
}

/*
This test checks varying the gain with beamtype=ANALY_DIPOLE
*/
void test_update_sum_visis_VaryGainAnalyBeam_gpu(void) {
  test_update_sum_visis_VaryGainChooseBeams(ANALY_DIPOLE, 1);
}

/*
This test checks varying the gain with beamtype=GAUSS_BEAM
*/
void test_update_sum_visis_VaryGainGaussBeam_gpu(void) {
  test_update_sum_visis_VaryGainChooseBeams(GAUSS_BEAM, 1);
}

/*
This test checks varying the gain with beamtype=NO_BEAM
*/
void test_update_sum_visis_VaryGainNoBeam_gpu(void) {
  test_update_sum_visis_VaryGainChooseBeams(NO_BEAM, 1);
}

/*
This test checks varying the measurement equation with beamtype=FEE_BEAM_INTERP
*/
void test_update_sum_visis_VaryGainFEEInterpBeam_gpu(void) {
  test_update_sum_visis_VaryGainChooseBeams(FEE_BEAM_INTERP, 1);
}

/*
This test checks varying the measurement equation with beamtype=MWA_ANALY
*/
void test_update_sum_visis_VaryGainMWAAnaly_gpu(void) {
  test_update_sum_visis_VaryGainChooseBeams(MWA_ANALY, 1);
}



/*
This test checks varying the gain with beamtype=FEE_BEAM
*/
void test_update_sum_visis_VaryFluxesFEEBeam_gpu(void) {
  test_update_sum_visis_VaryFluxesChooseBeams(FEE_BEAM, 1);
}

/*
This test checks varying the gain with beamtype=ANALY_DIPOLE
*/
void test_update_sum_visis_VaryFluxesAnalyBeam_gpu(void) {
  test_update_sum_visis_VaryFluxesChooseBeams(ANALY_DIPOLE, 1);
}

/*
This test checks varying the gain with beamtype=GAUSS_BEAM
*/
void test_update_sum_visis_VaryFluxesGaussBeam_gpu(void) {
  test_update_sum_visis_VaryFluxesChooseBeams(GAUSS_BEAM, 1);
}

/*
This test checks varying the gain with beamtype=NO_BEAM
*/
void test_update_sum_visis_VaryFluxesNoBeam_gpu(void) {
  test_update_sum_visis_VaryFluxesChooseBeams(NO_BEAM, 1);
}

/*
This test checks varying the measurement equation with beamtype=FEE_BEAM_INTERP
*/
void test_update_sum_visis_VaryFluxesFEEInterpBeam_gpu(void) {
  test_update_sum_visis_VaryFluxesChooseBeams(FEE_BEAM_INTERP, 1);
}

/*
This test checks varying the measurement equation with beamtype=MWA_ANALY
*/
void test_update_sum_visis_VaryFluxesMWAAnaly_gpu(void) {
  test_update_sum_visis_VaryFluxesChooseBeams(MWA_ANALY, 1);
}



/*
This test checks varying the measurement equation with beamtype=FEE_BEAM
*/
void test_update_sum_visis_VaryVisiFEEBeam_gpu(void) {
  test_update_sum_visis_VaryVisiChooseBeams(FEE_BEAM, 1);
}

/*
This test checks varying the measurement equation with beamtype=ANALY_DIPOLE
*/
void test_update_sum_visis_VaryVisiAnalyBeam_gpu(void) {
  test_update_sum_visis_VaryVisiChooseBeams(ANALY_DIPOLE, 1);
}

/*
This test checks varying the measurement equation with beamtype=GAUSS_BEAM
*/
void test_update_sum_visis_VaryVisiGaussBeam_gpu(void) {
  test_update_sum_visis_VaryVisiChooseBeams(GAUSS_BEAM, 1);
}

/*
This test checks varying the measurement equation with beamtype=NO_BEAM
*/
void test_update_sum_visis_VaryVisiNoBeam_gpu(void) {
  test_update_sum_visis_VaryVisiChooseBeams(NO_BEAM, 1);
}

/*
This test checks varying the measurement equation with beamtype=FEE_BEAM_INTERP
*/
void test_update_sum_visis_VaryVisiFEEInterpBeam_gpu(void) {
  test_update_sum_visis_VaryVisiChooseBeams(FEE_BEAM_INTERP, 1);
}

/*
This test checks varying the measurement equation with beamtype=MWA_ANALY
*/
void test_update_sum_visis_VaryVisiMWAAnaly_gpu(void) {
  test_update_sum_visis_VaryVisiChooseBeams(MWA_ANALY, 1);
}

//Run the test with unity
int main(void)
{
    UNITY_BEGIN();
    //Test while varying beam gain for all beam types
    RUN_TEST(test_update_sum_visis_VaryGainFEEBeam_gpu);
    RUN_TEST(test_update_sum_visis_VaryGainAnalyBeam_gpu);
    RUN_TEST(test_update_sum_visis_VaryGainGaussBeam_gpu);
    RUN_TEST(test_update_sum_visis_VaryGainNoBeam_gpu);
    RUN_TEST(test_update_sum_visis_VaryGainFEEInterpBeam_gpu);
    RUN_TEST(test_update_sum_visis_VaryGainMWAAnaly_gpu);

    //Test while varying component fluxes for all beam types
    RUN_TEST(test_update_sum_visis_VaryFluxesFEEBeam_gpu);
    RUN_TEST(test_update_sum_visis_VaryFluxesAnalyBeam_gpu);
    RUN_TEST(test_update_sum_visis_VaryFluxesNoBeam_gpu);
    RUN_TEST(test_update_sum_visis_VaryFluxesGaussBeam_gpu);
    RUN_TEST(test_update_sum_visis_VaryFluxesFEEInterpBeam_gpu);
    RUN_TEST(test_update_sum_visis_VaryFluxesMWAAnaly_gpu);

    //Test while varying base visibility for all beam types
    RUN_TEST(test_update_sum_visis_VaryVisiFEEBeam_gpu);
    RUN_TEST(test_update_sum_visis_VaryVisiAnalyBeam_gpu);
    RUN_TEST(test_update_sum_visis_VaryVisiNoBeam_gpu);
    RUN_TEST(test_update_sum_visis_VaryVisiGaussBeam_gpu);
    RUN_TEST(test_update_sum_visis_VaryVisiFEEInterpBeam_gpu);
    RUN_TEST(test_update_sum_visis_VaryVisiMWAAnaly_gpu);

    return UNITY_END();
}
