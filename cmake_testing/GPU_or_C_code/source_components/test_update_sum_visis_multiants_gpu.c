#include "update_sum_visis_multiants_common.h"


/*
This test checks varying the gain with beamtype=FEE_BEAM
*/
void test_update_sum_visis_VaryGainFEEBeam_gpu(void) {
  test_update_sum_visis_multiants_VaryGainChooseBeams(FEE_BEAM, 1);
}

/*
This test checks varying the measurement equation with beamtype=FEE_BEAM_INTERP
*/
void test_update_sum_visis_VaryGainFEEInterpBeam_gpu(void) {
  test_update_sum_visis_multiants_VaryGainChooseBeams(FEE_BEAM_INTERP, 1);
}

// Useful for possible future expansion
// /*
// This test checks varying the gain with beamtype=ANALY_DIPOLE
// */
// void test_update_sum_visis_VaryGainAnalyBeam_gpu(void) {
//   test_update_sum_visis_multiants_VaryGainChooseBeams(ANALY_DIPOLE, 1);
// }

// /*
// This test checks varying the gain with beamtype=GAUSS_BEAM
// */
// void test_update_sum_visis_VaryGainGaussBeam_gpu(void) {
//   test_update_sum_visis_multiants_VaryGainChooseBeams(GAUSS_BEAM, 1);
// }

// /*
// This test checks varying the gain with beamtype=NO_BEAM
// */
// void test_update_sum_visis_VaryGainNoBeam_gpu(void) {
//   test_update_sum_visis_multiants_VaryGainChooseBeams(NO_BEAM, 1);
// }



// /*
// This test checks varying the measurement equation with beamtype=MWA_ANALY
// */
// void test_update_sum_visis_VaryGainMWAAnaly_gpu(void) {
//   test_update_sum_visis_multiants_VaryGainChooseBeams(MWA_ANALY, 1);
// }

//Run the test with unity
int main(void)
{
    UNITY_BEGIN();
    //Test while varying beam gain for all beam types
    RUN_TEST(test_update_sum_visis_VaryGainFEEBeam_gpu);
    RUN_TEST(test_update_sum_visis_VaryGainFEEInterpBeam_gpu);
    // RUN_TEST(test_update_sum_visis_VaryGainAnalyBeam_gpu);
    // RUN_TEST(test_update_sum_visis_VaryGainGaussBeam_gpu);
    // RUN_TEST(test_update_sum_visis_VaryGainNoBeam_gpu);
    // RUN_TEST(test_update_sum_visis_VaryGainMWAAnaly_gpu);

    return UNITY_END();
}
