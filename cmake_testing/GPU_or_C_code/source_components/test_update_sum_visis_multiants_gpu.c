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

void test_update_sum_visis_VaryGainLOFAR_gpu(void) {
  test_update_sum_visis_multiants_VaryGainChooseBeams(EB_LOFAR, 1);
}

void test_update_sum_visis_VaryGainOSKAR_gpu(void) {
  test_update_sum_visis_multiants_VaryGainChooseBeams(EB_OSKAR, 1);
}

//Run the test with unity
int main(void)
{
    UNITY_BEGIN();
    //Test while varying beam gain for all beam types
    RUN_TEST(test_update_sum_visis_VaryGainFEEBeam_gpu);
    RUN_TEST(test_update_sum_visis_VaryGainFEEInterpBeam_gpu);
    RUN_TEST(test_update_sum_visis_VaryGainLOFAR_gpu);
    RUN_TEST(test_update_sum_visis_VaryGainOSKAR_gpu);

    return UNITY_END();
}
