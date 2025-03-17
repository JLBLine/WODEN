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

void test_update_sum_visis_VaryGainEBLOFAR_cpu(void) {
  test_update_sum_visis_multiants_VaryGainChooseBeams(EB_LOFAR, 0);
}

void test_update_sum_visis_VaryGainEBOSKAR_cpu(void) {
  test_update_sum_visis_multiants_VaryGainChooseBeams(EB_OSKAR, 0);
}



//Run the test with unity
int main(void)
{
    UNITY_BEGIN();
    //Test while varying beam gain for all beam types
    RUN_TEST(test_update_sum_visis_VaryGainFEEBeam_cpu);
    RUN_TEST(test_update_sum_visis_VaryGainFEEInterpBeam_cpu);
    RUN_TEST(test_update_sum_visis_VaryGainEBLOFAR_cpu);
    RUN_TEST(test_update_sum_visis_VaryGainEBOSKAR_cpu);


    return UNITY_END();
}
