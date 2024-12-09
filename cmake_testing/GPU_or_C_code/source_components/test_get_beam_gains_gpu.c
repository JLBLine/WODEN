#include <math.h>
#include <unity.h>
#include <stdlib.h>
#include <complex.h>

#include "constants.h"
#include "woden_struct_defs.h"
#include "woden_precision_defs.h"
#include "get_beam_gains_common.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

#define UNITY_INCLUDE_FLOAT

/*
This test checks the correct behaviour when simuating with beamtype=FEE_BEAM
*/
void test_get_beam_gains_FEEBeam_gpu(void) {
  test_get_beam_gains_ChooseBeams(FEE_BEAM, 1);
}

/*
This test checks the correct behaviour when simuating with beamtype=ANALY_DIPOLE
*/
void test_get_beam_gains_AnalyDipoleBeam_gpu(void) {
  test_get_beam_gains_ChooseBeams(ANALY_DIPOLE, 1);
}

/*
This test checks the correct behaviour when simuating with beamtype=GAUSS_BEAM
*/
void test_get_beam_gains_GaussBeam_gpu(void) {
  test_get_beam_gains_ChooseBeams(GAUSS_BEAM, 1);
}

/*
This test checks the correct behaviour when simuating with beamtype=NO_BEAM
*/
void test_get_beam_gains_NoBeam_gpu(void) {
  test_get_beam_gains_ChooseBeams(NO_BEAM, 1);
}

/*
This test checks the correct behaviour when simuating with beamtype=FEE_BEAM_INTERP
*/
void test_get_beam_gains_FEEBeamInterp_gpu(void) {
  test_get_beam_gains_ChooseBeams(FEE_BEAM_INTERP, 1);
}

/*
This test checks the correct behaviour when simuating with beamtype=MWA_ANALY
*/
void test_get_beam_gains_MWAAnaly_gpu(void) {
  test_get_beam_gains_ChooseBeams(MWA_ANALY, 1);
}


//Run the test with unity
int main(void)
{
    UNITY_BEGIN();
    RUN_TEST(test_get_beam_gains_FEEBeam_gpu);
    RUN_TEST(test_get_beam_gains_AnalyDipoleBeam_gpu);
    RUN_TEST(test_get_beam_gains_GaussBeam_gpu);
    RUN_TEST(test_get_beam_gains_NoBeam_gpu);
    RUN_TEST(test_get_beam_gains_FEEBeamInterp_gpu);
    RUN_TEST(test_get_beam_gains_MWAAnaly_gpu);


    return UNITY_END();
}
