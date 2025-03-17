#include "update_sum_visis_common.h"

/*
This test checks varying the gain with beamtype=FEE_BEAM
*/
void test_update_sum_visis_VaryGainFEEBeam_cpu(void) {
  test_update_sum_visis_VaryGainChooseBeams(FEE_BEAM, 0);
}

/*
This test checks varying the gain with beamtype=ANALY_DIPOLE
*/
void test_update_sum_visis_VaryGainAnalyBeam_cpu(void) {
  test_update_sum_visis_VaryGainChooseBeams(ANALY_DIPOLE, 0);
}

/*
This test checks varying the gain with beamtype=GAUSS_BEAM
*/
void test_update_sum_visis_VaryGainGaussBeam_cpu(void) {
  test_update_sum_visis_VaryGainChooseBeams(GAUSS_BEAM, 0);
}

/*
This test checks varying the gain with beamtype=NO_BEAM
*/
void test_update_sum_visis_VaryGainNoBeam_cpu(void) {
  test_update_sum_visis_VaryGainChooseBeams(NO_BEAM, 0);
}

/*
This test checks varying the measurement equation with beamtype=FEE_BEAM_INTERP
*/
void test_update_sum_visis_VaryGainFEEInterpBeam_cpu(void) {
  test_update_sum_visis_VaryGainChooseBeams(FEE_BEAM_INTERP, 0);
}

/*
This test checks varying the measurement equation with beamtype=MWA_ANALY
*/
void test_update_sum_visis_VaryGainMWAAnaly_cpu(void) {
  test_update_sum_visis_VaryGainChooseBeams(MWA_ANALY, 0);
}



/*
This test checks varying the gain with beamtype=FEE_BEAM
*/
void test_update_sum_visis_VaryFluxesFEEBeam_cpu(void) {
  test_update_sum_visis_VaryFluxesChooseBeams(FEE_BEAM, 0);
}

/*
This test checks varying the gain with beamtype=ANALY_DIPOLE
*/
void test_update_sum_visis_VaryFluxesAnalyBeam_cpu(void) {
  test_update_sum_visis_VaryFluxesChooseBeams(ANALY_DIPOLE, 0);
}

/*
This test checks varying the gain with beamtype=GAUSS_BEAM
*/
void test_update_sum_visis_VaryFluxesGaussBeam_cpu(void) {
  test_update_sum_visis_VaryFluxesChooseBeams(GAUSS_BEAM, 0);
}

/*
This test checks varying the gain with beamtype=NO_BEAM
*/
void test_update_sum_visis_VaryFluxesNoBeam_cpu(void) {
  test_update_sum_visis_VaryFluxesChooseBeams(NO_BEAM, 0);
}

/*
This test checks varying the measurement equation with beamtype=FEE_BEAM_INTERP
*/
void test_update_sum_visis_VaryFluxesFEEInterpBeam_cpu(void) {
  test_update_sum_visis_VaryFluxesChooseBeams(FEE_BEAM_INTERP, 0);
}

/*
This test checks varying the measurement equation with beamtype=MWA_ANALY
*/
void test_update_sum_visis_VaryFluxesMWAAnaly_cpu(void) {
  test_update_sum_visis_VaryFluxesChooseBeams(MWA_ANALY, 0);
}



/*
This test checks varying the measurement equation with beamtype=FEE_BEAM
*/
void test_update_sum_visis_VaryVisiFEEBeam_cpu(void) {
  test_update_sum_visis_VaryVisiChooseBeams(FEE_BEAM, 0);
}

/*
This test checks varying the measurement equation with beamtype=ANALY_DIPOLE
*/
void test_update_sum_visis_VaryVisiAnalyBeam_cpu(void) {
  test_update_sum_visis_VaryVisiChooseBeams(ANALY_DIPOLE, 0);
}

/*
This test checks varying the measurement equation with beamtype=GAUSS_BEAM
*/
void test_update_sum_visis_VaryVisiGaussBeam_cpu(void) {
  test_update_sum_visis_VaryVisiChooseBeams(GAUSS_BEAM, 0);
}

/*
This test checks varying the measurement equation with beamtype=NO_BEAM
*/
void test_update_sum_visis_VaryVisiNoBeam_cpu(void) {
  test_update_sum_visis_VaryVisiChooseBeams(NO_BEAM, 0);
}

/*
This test checks varying the measurement equation with beamtype=FEE_BEAM_INTERP
*/
void test_update_sum_visis_VaryVisiFEEInterpBeam_cpu(void) {
  test_update_sum_visis_VaryVisiChooseBeams(FEE_BEAM_INTERP, 0);
}

/*
This test checks varying the measurement equation with beamtype=MWA_ANALY
*/
void test_update_sum_visis_VaryVisiMWAAnaly_cpu(void) {
  test_update_sum_visis_VaryVisiChooseBeams(MWA_ANALY, 0);
}

//Run the test with unity
int main(void)
{
    UNITY_BEGIN();
    //Test while varying beam gain for all beam types
    RUN_TEST(test_update_sum_visis_VaryGainFEEBeam_cpu);
    RUN_TEST(test_update_sum_visis_VaryGainAnalyBeam_cpu);
    RUN_TEST(test_update_sum_visis_VaryGainGaussBeam_cpu);
    RUN_TEST(test_update_sum_visis_VaryGainNoBeam_cpu);
    RUN_TEST(test_update_sum_visis_VaryGainFEEInterpBeam_cpu);
    RUN_TEST(test_update_sum_visis_VaryGainMWAAnaly_cpu);

    //Test while varying component fluxes for all beam types
    RUN_TEST(test_update_sum_visis_VaryFluxesFEEBeam_cpu);
    RUN_TEST(test_update_sum_visis_VaryFluxesAnalyBeam_cpu);
    RUN_TEST(test_update_sum_visis_VaryFluxesNoBeam_cpu);
    RUN_TEST(test_update_sum_visis_VaryFluxesGaussBeam_cpu);
    RUN_TEST(test_update_sum_visis_VaryFluxesFEEInterpBeam_cpu);
    RUN_TEST(test_update_sum_visis_VaryFluxesMWAAnaly_cpu);

    //Test while varying base visibility for all beam types
    RUN_TEST(test_update_sum_visis_VaryVisiFEEBeam_cpu);
    RUN_TEST(test_update_sum_visis_VaryVisiAnalyBeam_cpu);
    RUN_TEST(test_update_sum_visis_VaryVisiNoBeam_cpu);
    RUN_TEST(test_update_sum_visis_VaryVisiGaussBeam_cpu);
    RUN_TEST(test_update_sum_visis_VaryVisiFEEInterpBeam_cpu);
    RUN_TEST(test_update_sum_visis_VaryVisiMWAAnaly_cpu);

    return UNITY_END();
}
