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


#include "calculate_visibilities_uvbeam_common.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */


//Test with a single SOURCE, single COMPONENT
void test_calculate_visibilities_UVBeamMWA_OneSource_SinglePoint_gpu(void) {
  test_calculate_visibilities_UVBeam_OneSource_SinglePoint(1, UVB_MWA);
}

void test_calculate_visibilities_UVBeamMWA_OneSource_SingleGauss_gpu(void) {
  test_calculate_visibilities_UVBeam_OneSource_SingleGauss(1, UVB_MWA);

}

void test_calculate_visibilities_UVBeamMWA_OneSource_SingleShape_gpu(void) {
  test_calculate_visibilities_UVBeam_OneSource_SingleShape(1, UVB_MWA);
}

void test_calculate_visibilities_UVBeamMWA_OneSource_SingleAll_gpu(void) {
  test_calculate_visibilities_UVBeam_OneSource_SingleAll(1, UVB_MWA);
}


//Test with a three SOURCEs, single COMPONENT
void test_calculate_visibilities_UVBeamMWA_ThreeSource_SinglePoint_gpu(void) {
  test_calculate_visibilities_UVBeam_ThreeSource_SinglePoint(1, UVB_MWA);

}

void test_calculate_visibilities_UVBeamMWA_ThreeSource_SingleGauss_gpu(void) {
  test_calculate_visibilities_UVBeam_ThreeSource_SingleGauss(1, UVB_MWA);
}

void test_calculate_visibilities_UVBeamMWA_ThreeSource_SingleShape_gpu(void) {
  test_calculate_visibilities_UVBeam_ThreeSource_SingleShape(1, UVB_MWA);
}

void test_calculate_visibilities_UVBeamMWA_ThreeSource_SingleAll_gpu(void) {
  test_calculate_visibilities_UVBeam_ThreeSource_SingleAll(1, UVB_MWA);
}

//Test with three SOURCEs, three COPMONENTs
void test_calculate_visibilities_UVBeamMWA_ThreeSource_FivePoint_gpu(void) {
  test_calculate_visibilities_UVBeam_ThreeSource_FivePoint(1, UVB_MWA);

}

void test_calculate_visibilities_UVBeamMWA_ThreeSource_FiveGauss_gpu(void) {
  test_calculate_visibilities_UVBeam_ThreeSource_FiveGauss(1, UVB_MWA);
}

void test_calculate_visibilities_UVBeamMWA_ThreeSource_FiveShape_gpu(void) {
  test_calculate_visibilities_UVBeam_ThreeSource_FiveShape(1, UVB_MWA);
}

void test_calculate_visibilities_UVBeamMWA_ThreeSource_FiveAll_gpu(void) {
  test_calculate_visibilities_UVBeam_ThreeSource_FiveAll(1, UVB_MWA);
}


//Run the test with unity
int main(void)
{
    UNITY_BEGIN();
    // //Test with a single SOURCE, single COMPONENT
    RUN_TEST(test_calculate_visibilities_UVBeamMWA_OneSource_SinglePoint_gpu);
    RUN_TEST(test_calculate_visibilities_UVBeamMWA_OneSource_SingleGauss_gpu);
    RUN_TEST(test_calculate_visibilities_UVBeamMWA_OneSource_SingleShape_gpu);
    RUN_TEST(test_calculate_visibilities_UVBeamMWA_OneSource_SingleAll_gpu);

    //Test with three SOURCEs, single COPMONENT
    RUN_TEST(test_calculate_visibilities_UVBeamMWA_ThreeSource_SinglePoint_gpu);
    RUN_TEST(test_calculate_visibilities_UVBeamMWA_ThreeSource_SingleGauss_gpu);
    RUN_TEST(test_calculate_visibilities_UVBeamMWA_ThreeSource_SingleShape_gpu);
    RUN_TEST(test_calculate_visibilities_UVBeamMWA_ThreeSource_SingleAll_gpu);

    //Test with three SOURCEs, five COPMONENTs
    RUN_TEST(test_calculate_visibilities_UVBeamMWA_ThreeSource_FivePoint_gpu);
    RUN_TEST(test_calculate_visibilities_UVBeamMWA_ThreeSource_FiveGauss_gpu);
    RUN_TEST(test_calculate_visibilities_UVBeamMWA_ThreeSource_FiveShape_gpu);
    RUN_TEST(test_calculate_visibilities_UVBeamMWA_ThreeSource_FiveAll_gpu);

    return UNITY_END();
}