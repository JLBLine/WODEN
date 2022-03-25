#include <math.h>
#include <unity.h>
#include <stdlib.h>
#include <complex.h>

#include "constants.h"
#include "woden_struct_defs.h"
#include "test_kern_calc_visi_common.h"

// void sincos(user_precision_t x, user_precision_t *sin, user_precision_t *cos);

void setUp (void) {} /* Is run before eVary test, put unit init calls here. */
void tearDown (void) {} /* Is run after eVary test, put unit clean-up calls here. */



/*
These tests multiple directions on the sky, with constant fluxes and beam gains
Different beamtypes test how the gains are mapped
*/
void test_kern_calc_visi_gauss_VarylmnFEEBeam(void) {
  test_kern_calc_visi_Varylmn(FEE_BEAM, GAUSSIAN);
}

void test_kern_calc_visi_gauss_VarylmnAnalyBeam(void) {
  test_kern_calc_visi_Varylmn(ANALY_DIPOLE, GAUSSIAN);
}

void test_kern_calc_visi_gauss_VarylmnGaussBeam(void) {
  test_kern_calc_visi_Varylmn(GAUSS_BEAM, GAUSSIAN);
}

void test_kern_calc_visi_gauss_VarylmnNoBeam(void) {
  test_kern_calc_visi_Varylmn(NO_BEAM, GAUSSIAN);
}

void test_kern_calc_visi_gauss_VarylmnFEEInterpBeam(void) {
  test_kern_calc_visi_Varylmn(FEE_BEAM_INTERP, GAUSSIAN);
}

void test_kern_calc_visi_gauss_VarylmnMWAAnaly(void) {
  test_kern_calc_visi_Varylmn(MWA_ANALY, GAUSSIAN);
}

/*
These tests multiple directions on the sky, with varying fluxes and
constant beam gains
Different beamtypes test how the gains are mapped
*/
void test_kern_calc_visi_gauss_VarylmnVaryFluxFEEBeam(void) {
  test_kern_calc_visi_VarylmnVaryFlux(FEE_BEAM, GAUSSIAN);
}

void test_kern_calc_visi_gauss_VarylmnVaryFluxAnalyBeam(void) {
  test_kern_calc_visi_VarylmnVaryFlux(ANALY_DIPOLE, GAUSSIAN);
}

void test_kern_calc_visi_gauss_VarylmnVaryFluxGaussBeam(void) {
  test_kern_calc_visi_VarylmnVaryFlux(GAUSS_BEAM, GAUSSIAN);
}

void test_kern_calc_visi_gauss_VarylmnVaryFluxNoBeam(void) {
  test_kern_calc_visi_VarylmnVaryFlux(NO_BEAM, GAUSSIAN);
}

void test_kern_calc_visi_gauss_VarylmnVaryFluxFEEInterpBeam(void) {
  test_kern_calc_visi_VarylmnVaryFlux(FEE_BEAM_INTERP, GAUSSIAN);
}

void test_kern_calc_visi_gauss_VarylmnVaryFluxMWAAnaly(void) {
  test_kern_calc_visi_VarylmnVaryFlux(MWA_ANALY, GAUSSIAN);
}

/*
These tests multiple directions on the sky, with varying beam gains and
constant fluxe
Different beamtypes test how the gains are mapped
*/
void test_kern_calc_visi_gauss_VarylmnVaryBeamFEEBeam(void) {
  test_kern_calc_visi_VarylmnVaryBeam(FEE_BEAM, GAUSSIAN);
}

void test_kern_calc_visi_gauss_VarylmnVaryBeamAnalyBeam(void) {
  test_kern_calc_visi_VarylmnVaryBeam(ANALY_DIPOLE, GAUSSIAN);
}

void test_kern_calc_visi_gauss_VarylmnVaryBeamGaussBeam(void) {
  test_kern_calc_visi_VarylmnVaryBeam(GAUSS_BEAM, GAUSSIAN);
}

void test_kern_calc_visi_gauss_VarylmnVaryBeamNoBeam(void) {
  test_kern_calc_visi_VarylmnVaryBeam(NO_BEAM, GAUSSIAN);
}

void test_kern_calc_visi_gauss_VarylmnVaryBeamFEEInterpBeam(void) {
  test_kern_calc_visi_VarylmnVaryBeam(FEE_BEAM_INTERP, GAUSSIAN);
}

void test_kern_calc_visi_gauss_VarylmnVaryBeamMWAAnaly(void) {
  test_kern_calc_visi_VarylmnVaryBeam(MWA_ANALY, GAUSSIAN);
}




/*
These tests multiple directions on the sky, varying the position angle, major
ann minor axes
Different beamtypes test how the gains are mapped
*/
void test_kern_calc_visi_gauss_VarylmnVaryPAMajMinFEEBeam(void) {
  test_kern_calc_visi_VarylmnVaryPAMajMin(FEE_BEAM, GAUSSIAN);
}

void test_kern_calc_visi_gauss_VarylmnVaryPAMajMinAnalyBeam(void) {
  test_kern_calc_visi_VarylmnVaryPAMajMin(ANALY_DIPOLE, GAUSSIAN);
}

void test_kern_calc_visi_gauss_VarylmnVaryPAMajMinGaussBeam(void) {
  test_kern_calc_visi_VarylmnVaryPAMajMin(GAUSS_BEAM, GAUSSIAN);
}

void test_kern_calc_visi_gauss_VarylmnVaryPAMajMinNoBeam(void) {
  test_kern_calc_visi_VarylmnVaryPAMajMin(NO_BEAM, GAUSSIAN);
}

void test_kern_calc_visi_gauss_VarylmnVaryPAMajMinFEEInterpBeam(void) {
  test_kern_calc_visi_VarylmnVaryPAMajMin(FEE_BEAM_INTERP, GAUSSIAN);
}

void test_kern_calc_visi_gauss_VarylmnVaryPAMajMinMWAAnaly(void) {
  test_kern_calc_visi_VarylmnVaryPAMajMin(MWA_ANALY, GAUSSIAN);
}

//Run the test with unity
int main(void)
{
    UNITY_BEGIN();
    //Test while varying beam gain for all beam types
    RUN_TEST(test_kern_calc_visi_gauss_VarylmnFEEBeam);
    RUN_TEST(test_kern_calc_visi_gauss_VarylmnGaussBeam);
    RUN_TEST(test_kern_calc_visi_gauss_VarylmnAnalyBeam);
    RUN_TEST(test_kern_calc_visi_gauss_VarylmnNoBeam);
    RUN_TEST(test_kern_calc_visi_gauss_VarylmnFEEInterpBeam);
    RUN_TEST(test_kern_calc_visi_gauss_VarylmnMWAAnaly);

    RUN_TEST(test_kern_calc_visi_gauss_VarylmnVaryFluxFEEBeam);
    RUN_TEST(test_kern_calc_visi_gauss_VarylmnVaryFluxGaussBeam);
    RUN_TEST(test_kern_calc_visi_gauss_VarylmnVaryFluxAnalyBeam);
    RUN_TEST(test_kern_calc_visi_gauss_VarylmnVaryFluxNoBeam);
    RUN_TEST(test_kern_calc_visi_gauss_VarylmnVaryFluxFEEInterpBeam);
    RUN_TEST(test_kern_calc_visi_gauss_VarylmnVaryFluxMWAAnaly);

    RUN_TEST(test_kern_calc_visi_gauss_VarylmnVaryBeamFEEBeam);
    RUN_TEST(test_kern_calc_visi_gauss_VarylmnVaryBeamGaussBeam);
    RUN_TEST(test_kern_calc_visi_gauss_VarylmnVaryBeamAnalyBeam);
    RUN_TEST(test_kern_calc_visi_gauss_VarylmnVaryBeamNoBeam);
    RUN_TEST(test_kern_calc_visi_gauss_VarylmnVaryBeamFEEInterpBeam);
    RUN_TEST(test_kern_calc_visi_gauss_VarylmnVaryBeamMWAAnaly);

    RUN_TEST(test_kern_calc_visi_gauss_VarylmnVaryPAMajMinFEEBeam);
    RUN_TEST(test_kern_calc_visi_gauss_VarylmnVaryPAMajMinGaussBeam);
    RUN_TEST(test_kern_calc_visi_gauss_VarylmnVaryPAMajMinAnalyBeam);
    RUN_TEST(test_kern_calc_visi_gauss_VarylmnVaryPAMajMinNoBeam);
    RUN_TEST(test_kern_calc_visi_gauss_VarylmnVaryPAMajMinFEEInterpBeam);
    RUN_TEST(test_kern_calc_visi_gauss_VarylmnVaryPAMajMinMWAAnaly);

    return UNITY_END();
}
