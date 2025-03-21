#include "get_beam_gains_multiants_common.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */



// /*
// This test checks the correct behaviour when simuating with beamtype=MWA_ANALY
// */
// void test_get_beam_gains_MWAAnaly_cpu(void) {
//   test_get_beam_gains_ChooseBeams(MWA_ANALY, 0);
// }

/*
This test checks the correct behaviour when simuating with beamtype=FEE_BEAM
*/
void test_get_beam_gains_FEEBeam_cpu(void) {
  test_get_beam_gains_ChooseBeams(FEE_BEAM, 0);
}


/*
This test checks the correct behaviour when simuating with beamtype=FEE_BEAM_INTERP
*/
void test_get_beam_gains_FEEBeamInterp_cpu(void) {
  test_get_beam_gains_ChooseBeams(FEE_BEAM_INTERP, 0);
}

// void test_get_beam_gains_EBMWA_cpu(void) {
//   test_get_beam_gains_ChooseBeams(EB_MWA, 0);
// }

void test_get_beam_gains_EBLOFAR_cpu(void) {
  test_get_beam_gains_ChooseBeams(EB_LOFAR, 0);
}

void test_get_beam_gains_EBOSKAR_cpu(void) {
  test_get_beam_gains_ChooseBeams(EB_OSKAR, 0);
}


//Run the test with unity
int main(void)
{
    UNITY_BEGIN();

    // RUN_TEST(test_get_beam_gains_AnalyDipoleBeam_cpu);
    RUN_TEST(test_get_beam_gains_FEEBeam_cpu);
    RUN_TEST(test_get_beam_gains_FEEBeamInterp_cpu);
    // RUN_TEST(test_get_beam_gains_EBMWA_cpu);
    RUN_TEST(test_get_beam_gains_EBLOFAR_cpu);
    RUN_TEST(test_get_beam_gains_EBOSKAR_cpu);
    
    return UNITY_END();
}
