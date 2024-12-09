/*
Tests that the kernel that calculates auto-correlations is doing it's job
*/
#include "calc_autos_multiants_common.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */


//Potentially useful for the future; at the moment only MWA FEE and Everybeam 
//beams use different beams for each antenna
// void test_calculate_autos_multiants_NoBeam_gpu(void) {
//   test_calculate_autos_multiants(NO_BEAM, 1);
// }

// void test_calculate_autos_multiants_GaussBeam_gpu(void) {
//   test_calculate_autos_multiants(GAUSS_BEAM, 1);
// }

// void test_calculate_autos_multiants_EDA2Beam_gpu(void) {
//   test_calculate_autos_multiants(ANALY_DIPOLE, 1);
// }

// void test_calculate_autos_multiants_MWAAnaly_gpu(void) {
//   test_calculate_autos_multiants(ANALY_DIPOLE, 1);
// }

void test_calculate_autos_multiants_MWAFEE_gpu(void) {
  test_calculate_autos_multiants(FEE_BEAM, 1);
}

void test_calculate_autos_multiants_MWAFEEInterp_gpu(void) {
  test_calculate_autos_multiants(FEE_BEAM_INTERP, 1);
}


//Run the test with unity
int main(void)
{
    UNITY_BEGIN();
    //Test with a single SOURCE, single COMPONENT

    RUN_TEST(test_calculate_autos_multiants_MWAFEE_gpu);
    RUN_TEST(test_calculate_autos_multiants_MWAFEEInterp_gpu);

    return UNITY_END();
}
