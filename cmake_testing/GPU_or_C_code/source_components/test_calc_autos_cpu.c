/*
Tests that the kernel that calculates auto-correlations is doing it's job
*/
#include "calc_autos_common.h"

void setUp (void) {} /* Is run before eVary test, put unit init calls here. */
void tearDown (void) {} /* Is run after eVary test, put unit clean-up calls here. */

void test_calculate_autos_NoBeam_cpu(void) {
  test_calculate_autos(NO_BEAM, 0);
}

void test_calculate_autos_GaussBeam_cpu(void) {
  test_calculate_autos(GAUSS_BEAM, 0);
}

void test_calculate_autos_EDA2Beam_cpu(void) {
  test_calculate_autos(ANALY_DIPOLE, 0);
}

void test_calculate_autos_MWAAnaly_cpu(void) {
  test_calculate_autos(ANALY_DIPOLE, 0);
}

void test_calculate_autos_MWAFEE_cpu(void) {
  test_calculate_autos(FEE_BEAM, 0);
}

void test_calculate_autos_MWAFEEInterp_cpu(void) {
  test_calculate_autos(FEE_BEAM_INTERP, 0);
}


//Run the test with unity
int main(void)
{
    UNITY_BEGIN();
    //Test with a single SOURCE, single COMPONENT

    RUN_TEST(test_calculate_autos_NoBeam_cpu);
    RUN_TEST(test_calculate_autos_GaussBeam_cpu);
    RUN_TEST(test_calculate_autos_EDA2Beam_cpu);
    RUN_TEST(test_calculate_autos_MWAAnaly_cpu);
    RUN_TEST(test_calculate_autos_MWAFEE_cpu);
    RUN_TEST(test_calculate_autos_MWAFEEInterp_cpu);

    return UNITY_END();
}
