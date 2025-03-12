#include "get_beam_gains_multiants_common.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */



// /*
// This test checks the correct behaviour when simuating with beamtype=MWA_ANALY
// */
// void test_get_beam_gains_MWAAnaly_gpu(void) {
//   test_get_beam_gains_ChooseBeams(MWA_ANALY, 1);
// }

/*
This test checks the correct behaviour when simuating with beamtype=FEE_BEAM
*/
void test_get_beam_gains_FEEBeam_gpu(void) {
  test_get_beam_gains_ChooseBeams(FEE_BEAM, 1);
}


/*
This test checks the correct behaviour when simuating with beamtype=FEE_BEAM_INTERP
*/
void test_get_beam_gains_FEEBeamInterp_gpu(void) {
  test_get_beam_gains_ChooseBeams(FEE_BEAM_INTERP, 1);
}

// void test_get_beam_gains_EBMWA_gpu(void) {
//   test_get_beam_gains_ChooseBeams(EB_MWA, 1);
// }

void test_get_beam_gains_EBLofar_gpu(void) {
  test_get_beam_gains_ChooseBeams(EB_LOFAR, 1);
}

void test_get_beam_gains_EBOskar_gpu(void) {
  test_get_beam_gains_ChooseBeams(EB_OSKAR, 1);
}


//Run the test with unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_get_beam_gains_FEEBeam_gpu);
    RUN_TEST(test_get_beam_gains_FEEBeamInterp_gpu);
    // RUN_TEST(test_get_beam_gains_EBMWA_gpu);
    RUN_TEST(test_get_beam_gains_EBLofar_gpu);
    RUN_TEST(test_get_beam_gains_EBOskar_gpu);
    
    return UNITY_END();
}
