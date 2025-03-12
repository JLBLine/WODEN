/*
Tests that the kernel that calculates auto-correlations is doing it's job
*/
#include "calc_autos_common.h"

void setUp (void) {} /* Is run before eVary test, put unit init calls here. */
void tearDown (void) {} /* Is run after eVary test, put unit clean-up calls here. */

void test_calculate_autos_NoBeam_gpu(void) {
  test_calculate_autos(NO_BEAM, 1);
}

void test_calculate_autos_GaussBeam_gpu(void) {
  test_calculate_autos(GAUSS_BEAM, 1);
}

void test_calculate_autos_EDA2Beam_gpu(void) {
  test_calculate_autos(ANALY_DIPOLE, 1);
}

void test_calculate_autos_MWAAnaly_gpu(void) {
  test_calculate_autos(ANALY_DIPOLE, 1);
}

void test_calculate_autos_MWAFEE_gpu(void) {
  test_calculate_autos(FEE_BEAM, 1);
}

void test_calculate_autos_MWAFEEInterp_gpu(void) {
  test_calculate_autos(FEE_BEAM_INTERP, 1);
}

void test_calculate_autos_EBMWA_gpu(void) {
  test_calculate_autos(EB_MWA, 1);
}

void test_calculate_autos_EBLOFAR_gpu(void) {
  test_calculate_autos(EB_LOFAR, 1);
}

void test_calculate_autos_EBOSKAR_gpu(void) {
  test_calculate_autos(EB_OSKAR, 1);
}


//Run the test with unity
int main(void)
{
    UNITY_BEGIN();
    //Test with a single SOURCE, single COMPONENT

    RUN_TEST(test_calculate_autos_NoBeam_gpu);
    RUN_TEST(test_calculate_autos_GaussBeam_gpu);
    RUN_TEST(test_calculate_autos_EDA2Beam_gpu);
    RUN_TEST(test_calculate_autos_MWAAnaly_gpu);
    RUN_TEST(test_calculate_autos_MWAFEE_gpu);
    RUN_TEST(test_calculate_autos_MWAFEEInterp_gpu);
    RUN_TEST(test_calculate_autos_EBMWA_gpu);
    RUN_TEST(test_calculate_autos_EBLOFAR_gpu);
    RUN_TEST(test_calculate_autos_EBOSKAR_gpu);

    return UNITY_END();
}
