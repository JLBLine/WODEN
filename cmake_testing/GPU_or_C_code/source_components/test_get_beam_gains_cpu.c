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
void test_get_beam_gains_FEEBeam_cpu(void) {
  test_get_beam_gains_ChooseBeams(FEE_BEAM, 0);
}

/*
This test checks the correct behaviour when simuating with beamtype=ANALY_DIPOLE
*/
void test_get_beam_gains_AnalyDipoleBeam_cpu(void) {
  test_get_beam_gains_ChooseBeams(ANALY_DIPOLE, 0);
}

/*
This test checks the correct behaviour when simuating with beamtype=GAUSS_BEAM
*/
void test_get_beam_gains_GaussBeam_cpu(void) {
  test_get_beam_gains_ChooseBeams(GAUSS_BEAM, 0);
}

/*
This test checks the correct behaviour when simuating with beamtype=NO_BEAM
*/
void test_get_beam_gains_NoBeam_cpu(void) {
  test_get_beam_gains_ChooseBeams(NO_BEAM, 0);
}

/*
This test checks the correct behaviour when simuating with beamtype=FEE_BEAM_INTERP
*/
void test_get_beam_gains_FEEBeamInterp_cpu(void) {
  test_get_beam_gains_ChooseBeams(FEE_BEAM_INTERP, 0);
}

/*
This test checks the correct behaviour when simuating with beamtype=MWA_ANALY
*/
void test_get_beam_gains_MWAAnaly_cpu(void) {
  test_get_beam_gains_ChooseBeams(MWA_ANALY, 0);
}


void test_get_beam_gains_EBMWA_cpu(void) {
  test_get_beam_gains_ChooseBeams(EB_MWA, 0);
}

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
    RUN_TEST(test_get_beam_gains_FEEBeam_cpu);
    RUN_TEST(test_get_beam_gains_AnalyDipoleBeam_cpu);
    RUN_TEST(test_get_beam_gains_GaussBeam_cpu);
    RUN_TEST(test_get_beam_gains_NoBeam_cpu);
    RUN_TEST(test_get_beam_gains_FEEBeamInterp_cpu);
    RUN_TEST(test_get_beam_gains_MWAAnaly_cpu);
    RUN_TEST(test_get_beam_gains_EBMWA_cpu);
    RUN_TEST(test_get_beam_gains_EBLOFAR_cpu);
    RUN_TEST(test_get_beam_gains_EBOSKAR_cpu);


    return UNITY_END();
}
