/*
`calculate_visibilities::calculate_visibilities` is the gateway function
to all CUDA functionality in WODEN. We'll test here in one baseline, frequency,
and time configuration. We'll vary the sky model and the primary beam. By
sticking all COMPONENTs at phase centre, we can just sum the expected fluxes
in XX / YY real to check things are being lanuched.

More variations like different phase centres / array configs etc are tested
in different test suites, so really just test that the correct CUDA functions
are launched by calculate_visibilities::calculate_visibilities`
*/


#include "calculate_visibilities_everybeam_common.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */


//Test with a single SOURCE, single COMPONENT
void test_calculate_visibilities_EveryBeamLOFAR_OneSource_SinglePoint_cpu(void) {
  test_calculate_visibilities_EveryBeam_OneSource_SinglePoint(0, EB_LOFAR);
}

void test_calculate_visibilities_EveryBeamLOFAR_OneSource_SingleGauss_cpu(void) {
  test_calculate_visibilities_EveryBeam_OneSource_SingleGauss(0, EB_LOFAR);

}

void test_calculate_visibilities_EveryBeamLOFAR_OneSource_SingleShape_cpu(void) {
  test_calculate_visibilities_EveryBeam_OneSource_SingleShape(0, EB_LOFAR);
}

void test_calculate_visibilities_EveryBeamLOFAR_OneSource_SingleAll_cpu(void) {
  test_calculate_visibilities_EveryBeam_OneSource_SingleAll(0, EB_LOFAR);
}


//Test with a three SOURCEs, single COMPONENT
void test_calculate_visibilities_EveryBeamLOFAR_ThreeSource_SinglePoint_cpu(void) {
  test_calculate_visibilities_EveryBeam_ThreeSource_SinglePoint(0, EB_LOFAR);

}

void test_calculate_visibilities_EveryBeamLOFAR_ThreeSource_SingleGauss_cpu(void) {
  test_calculate_visibilities_EveryBeam_ThreeSource_SingleGauss(0, EB_LOFAR);
}

void test_calculate_visibilities_EveryBeamLOFAR_ThreeSource_SingleShape_cpu(void) {
  test_calculate_visibilities_EveryBeam_ThreeSource_SingleShape(0, EB_LOFAR);
}

void test_calculate_visibilities_EveryBeamLOFAR_ThreeSource_SingleAll_cpu(void) {
  test_calculate_visibilities_EveryBeam_ThreeSource_SingleAll(0, EB_LOFAR);
}

//Test with three SOURCEs, three COPMONENTs
void test_calculate_visibilities_EveryBeamLOFAR_ThreeSource_FivePoint_cpu(void) {
  test_calculate_visibilities_EveryBeam_ThreeSource_FivePoint(0, EB_LOFAR);

}

void test_calculate_visibilities_EveryBeamLOFAR_ThreeSource_FiveGauss_cpu(void) {
  test_calculate_visibilities_EveryBeam_ThreeSource_FiveGauss(0, EB_LOFAR);
}

void test_calculate_visibilities_EveryBeamLOFAR_ThreeSource_FiveShape_cpu(void) {
  test_calculate_visibilities_EveryBeam_ThreeSource_FiveShape(0, EB_LOFAR);
}

void test_calculate_visibilities_EveryBeamLOFAR_ThreeSource_FiveAll_cpu(void) {
  test_calculate_visibilities_EveryBeam_ThreeSource_FiveAll(0, EB_LOFAR);
}


//Run the test with unity
int main(void)
{
    UNITY_BEGIN();
    // //Test with a single SOURCE, single COMPONENT
    RUN_TEST(test_calculate_visibilities_EveryBeamLOFAR_OneSource_SinglePoint_cpu);
    RUN_TEST(test_calculate_visibilities_EveryBeamLOFAR_OneSource_SingleGauss_cpu);
    RUN_TEST(test_calculate_visibilities_EveryBeamLOFAR_OneSource_SingleShape_cpu);
    RUN_TEST(test_calculate_visibilities_EveryBeamLOFAR_OneSource_SingleAll_cpu);

    //Test with three SOURCEs, single COPMONENT
    RUN_TEST(test_calculate_visibilities_EveryBeamLOFAR_ThreeSource_SinglePoint_cpu);
    RUN_TEST(test_calculate_visibilities_EveryBeamLOFAR_ThreeSource_SingleGauss_cpu);
    RUN_TEST(test_calculate_visibilities_EveryBeamLOFAR_ThreeSource_SingleShape_cpu);
    RUN_TEST(test_calculate_visibilities_EveryBeamLOFAR_ThreeSource_SingleAll_cpu);

    //Test with three SOURCEs, five COPMONENTs
    RUN_TEST(test_calculate_visibilities_EveryBeamLOFAR_ThreeSource_FivePoint_cpu);
    RUN_TEST(test_calculate_visibilities_EveryBeamLOFAR_ThreeSource_FiveGauss_cpu);
    RUN_TEST(test_calculate_visibilities_EveryBeamLOFAR_ThreeSource_FiveShape_cpu);
    RUN_TEST(test_calculate_visibilities_EveryBeamLOFAR_ThreeSource_FiveAll_cpu);

    return UNITY_END();
}
