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

#include "calculate_visibilities_mwafeeinterpbeam_common.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

//Test with a single SOURCE, single COMPONENT
void test_calculate_visibilities_MWAFEEBeamInterp_OneSource_SinglePoint_cpu(void) {
  test_calculate_visibilities_MWAFEEBeamInterp_OneSource_SinglePoint(0);
}

void test_calculate_visibilities_MWAFEEBeamInterp_OneSource_SingleGauss_cpu(void) {
  test_calculate_visibilities_MWAFEEBeamInterp_OneSource_SingleGauss(0);

}

void test_calculate_visibilities_MWAFEEBeamInterp_OneSource_SingleShape_cpu(void) {
  test_calculate_visibilities_MWAFEEBeamInterp_OneSource_SingleShape(0);
}

void test_calculate_visibilities_MWAFEEBeamInterp_OneSource_SingleAll_cpu(void) {
  test_calculate_visibilities_MWAFEEBeamInterp_OneSource_SingleAll(0);
}


//Test with a three SOURCEs, single COMPONENT
void test_calculate_visibilities_MWAFEEBeamInterp_ThreeSource_SinglePoint_cpu(void) {
  test_calculate_visibilities_MWAFEEBeamInterp_ThreeSource_SinglePoint(0);

}

void test_calculate_visibilities_MWAFEEBeamInterp_ThreeSource_SingleGauss_cpu(void) {
  test_calculate_visibilities_MWAFEEBeamInterp_ThreeSource_SingleGauss(0);
}

void test_calculate_visibilities_MWAFEEBeamInterp_ThreeSource_SingleShape_cpu(void) {
  test_calculate_visibilities_MWAFEEBeamInterp_ThreeSource_SingleShape(0);
}

void test_calculate_visibilities_MWAFEEBeamInterp_ThreeSource_SingleAll_cpu(void) {
  test_calculate_visibilities_MWAFEEBeamInterp_ThreeSource_SingleAll(0);
}


//Test with three SOURCEs, three COPMONENTs
void test_calculate_visibilities_MWAFEEBeamInterp_ThreeSource_ThreePoint_cpu(void) {
  test_calculate_visibilities_MWAFEEBeamInterp_ThreeSource_ThreePoint(0);

}

void test_calculate_visibilities_MWAFEEBeamInterp_ThreeSource_ThreeGauss_cpu(void) {
  test_calculate_visibilities_MWAFEEBeamInterp_ThreeSource_ThreeGauss(0);
}

void test_calculate_visibilities_MWAFEEBeamInterp_ThreeSource_ThreeShape_cpu(void) {
  test_calculate_visibilities_MWAFEEBeamInterp_ThreeSource_ThreeShape(0);
}

void test_calculate_visibilities_MWAFEEBeamInterp_ThreeSource_ThreeAll_cpu(void) {
  test_calculate_visibilities_MWAFEEBeamInterp_ThreeSource_ThreeAll(0);
}



// Run the test with unity
int main(void)
{
    UNITY_BEGIN();

    char* mwa_fee_hdf5 = getenv("MWA_FEE_HDF5");

    if (mwa_fee_hdf5) {
      printf("MWA_FEE_HDF5: %s\n", mwa_fee_hdf5 );

        //Test with a single SOURCE, single COMPONENT
        RUN_TEST(test_calculate_visibilities_MWAFEEBeamInterp_OneSource_SinglePoint_cpu);
        RUN_TEST(test_calculate_visibilities_MWAFEEBeamInterp_OneSource_SingleGauss_cpu);
        RUN_TEST(test_calculate_visibilities_MWAFEEBeamInterp_OneSource_SingleShape_cpu);
        RUN_TEST(test_calculate_visibilities_MWAFEEBeamInterp_OneSource_SingleAll_cpu);

        //Test with three SOURCEs, single COPMONENT
        RUN_TEST(test_calculate_visibilities_MWAFEEBeamInterp_ThreeSource_SinglePoint_cpu);
        RUN_TEST(test_calculate_visibilities_MWAFEEBeamInterp_ThreeSource_SingleGauss_cpu);
        RUN_TEST(test_calculate_visibilities_MWAFEEBeamInterp_ThreeSource_SingleShape_cpu);
        RUN_TEST(test_calculate_visibilities_MWAFEEBeamInterp_ThreeSource_SingleAll_cpu);

        //Test with three SOURCEs, three COPMONENTs
        RUN_TEST(test_calculate_visibilities_MWAFEEBeamInterp_ThreeSource_ThreePoint_cpu);
        RUN_TEST(test_calculate_visibilities_MWAFEEBeamInterp_ThreeSource_ThreeGauss_cpu);
        RUN_TEST(test_calculate_visibilities_MWAFEEBeamInterp_ThreeSource_ThreeShape_cpu);
        RUN_TEST(test_calculate_visibilities_MWAFEEBeamInterp_ThreeSource_ThreeAll_cpu);

    }
    else {
      printf("MWA_FEE_HDF5 not found - not running test_calculate_visibilities_MWAFEEBeamInterp tests");
    }

    return UNITY_END();
}
