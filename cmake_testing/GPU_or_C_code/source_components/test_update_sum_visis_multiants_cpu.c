#include "update_sum_visis_multiants_common.h"

/*
This test checks varying the gain with beamtype=FEE_BEAM
*/
void test_update_sum_visis_VaryGainFEEBeam_cpu(void) {
  test_update_sum_visis_multiants_VaryGainChooseBeams(FEE_BEAM, 0);
}

/*
This test checks varying the measurement equation with beamtype=FEE_BEAM_INTERP
*/
void test_update_sum_visis_VaryGainFEEInterpBeam_cpu(void) {
  test_update_sum_visis_multiants_VaryGainChooseBeams(FEE_BEAM_INTERP, 0);
}

// Useful for possible future expansion
// /*
// This test checks varying the gain with beamtype=ANALY_DIPOLE
// */
// void test_update_sum_visis_VaryGainAnalyBeam_cpu(void) {
//   test_update_sum_visis_multiants_VaryGainChooseBeams(ANALY_DIPOLE, 0);
// }

// /*
// This test checks varying the gain with beamtype=GAUSS_BEAM
// */
// void test_update_sum_visis_VaryGainGaussBeam_cpu(void) {
//   test_update_sum_visis_multiants_VaryGainChooseBeams(GAUSS_BEAM, 0);
// }

// /*
// This test checks varying the gain with beamtype=NO_BEAM
// */
// void test_update_sum_visis_VaryGainNoBeam_cpu(void) {
//   test_update_sum_visis_multiants_VaryGainChooseBeams(NO_BEAM, 0);
// }



// /*
// This test checks varying the measurement equation with beamtype=MWA_ANALY
// */
// void test_update_sum_visis_VaryGainMWAAnaly_cpu(void) {
//   test_update_sum_visis_multiants_VaryGainChooseBeams(MWA_ANALY, 0);
// }

//Run the test with unity
int main(void)
{
    UNITY_BEGIN();
    //Test while varying beam gain for all beam types
    RUN_TEST(test_update_sum_visis_VaryGainFEEBeam_cpu);
    RUN_TEST(test_update_sum_visis_VaryGainFEEInterpBeam_cpu);
    // RUN_TEST(test_update_sum_visis_VaryGainAnalyBeam_cpu);
    // RUN_TEST(test_update_sum_visis_VaryGainGaussBeam_cpu);
    // RUN_TEST(test_update_sum_visis_VaryGainNoBeam_cpu);
    // RUN_TEST(test_update_sum_visis_VaryGainMWAAnaly_cpu);

    return UNITY_END();
}
