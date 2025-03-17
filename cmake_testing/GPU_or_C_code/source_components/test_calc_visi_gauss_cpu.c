#include <math.h>
#include <unity.h>
#include <stdlib.h>
#include <complex.h>

#include "constants.h"
#include "woden_struct_defs.h"
#include "calc_visi_common.h"

// void sincos(user_precision_t x, user_precision_t *sin, user_precision_t *cos);

void setUp (void) {} /* Is run before eVary test, put unit init calls here. */
void tearDown (void) {} /* Is run after eVary test, put unit clean-up calls here. */

/*
These tests multiple directions on the sky, with constant fluxes and beam gains
Different beamtypes test how the gains are mapped
*/
void test_calc_visi_gauss_VarylmnFEEBeam_cpu(void) {
  test_calc_visi_Varylmn(FEE_BEAM, GAUSSIAN, 0);
}

void test_calc_visi_gauss_VarylmnAnalyBeam_cpu(void) {
  test_calc_visi_Varylmn(ANALY_DIPOLE, GAUSSIAN, 0);
}

void test_calc_visi_gauss_VarylmnGaussBeam_cpu(void) {
  test_calc_visi_Varylmn(GAUSS_BEAM, GAUSSIAN, 0);
}

void test_calc_visi_gauss_VarylmnNoBeam_cpu(void) {
  test_calc_visi_Varylmn(NO_BEAM, GAUSSIAN, 0);
}

void test_calc_visi_gauss_VarylmnFEEInterpBeam_cpu(void) {
  test_calc_visi_Varylmn(FEE_BEAM_INTERP, GAUSSIAN, 0);
}

void test_calc_visi_gauss_VarylmnMWAAnaly_cpu(void) {
  test_calc_visi_Varylmn(MWA_ANALY, GAUSSIAN, 0);
}

/*
These tests multiple directions on the sky, with varying fluxes and
constant beam gains
Different beamtypes test how the gains are mapped
*/
void test_calc_visi_gauss_VarylmnVaryFluxFEEBeam_cpu(void) {
  test_calc_visi_VarylmnVaryFlux(FEE_BEAM, GAUSSIAN, 0);
}

void test_calc_visi_gauss_VarylmnVaryFluxAnalyBeam_cpu(void) {
  test_calc_visi_VarylmnVaryFlux(ANALY_DIPOLE, GAUSSIAN, 0);
}

void test_calc_visi_gauss_VarylmnVaryFluxGaussBeam_cpu(void) {
  test_calc_visi_VarylmnVaryFlux(GAUSS_BEAM, GAUSSIAN, 0);
}

void test_calc_visi_gauss_VarylmnVaryFluxNoBeam_cpu(void) {
  test_calc_visi_VarylmnVaryFlux(NO_BEAM, GAUSSIAN, 0);
}

void test_calc_visi_gauss_VarylmnVaryFluxFEEInterpBeam_cpu(void) {
  test_calc_visi_VarylmnVaryFlux(FEE_BEAM_INTERP, GAUSSIAN, 0);
}

void test_calc_visi_gauss_VarylmnVaryFluxMWAAnaly_cpu(void) {
  test_calc_visi_VarylmnVaryFlux(MWA_ANALY, GAUSSIAN, 0);
}

/*
These tests multiple directions on the sky, with varying beam gains and
constant fluxe
Different beamtypes test how the gains are mapped
*/
void test_calc_visi_gauss_VarylmnVaryBeamFEEBeam_cpu(void) {
  test_calc_visi_VarylmnVaryBeam(FEE_BEAM, GAUSSIAN, 0);
}

void test_calc_visi_gauss_VarylmnVaryBeamAnalyBeam_cpu(void) {
  test_calc_visi_VarylmnVaryBeam(ANALY_DIPOLE, GAUSSIAN, 0);
}

void test_calc_visi_gauss_VarylmnVaryBeamGaussBeam_cpu(void) {
  test_calc_visi_VarylmnVaryBeam(GAUSS_BEAM, GAUSSIAN, 0);
}

void test_calc_visi_gauss_VarylmnVaryBeamNoBeam_cpu(void) {
  test_calc_visi_VarylmnVaryBeam(NO_BEAM, GAUSSIAN, 0);
}

void test_calc_visi_gauss_VarylmnVaryBeamFEEInterpBeam_cpu(void) {
  test_calc_visi_VarylmnVaryBeam(FEE_BEAM_INTERP, GAUSSIAN, 0);
}

void test_calc_visi_gauss_VarylmnVaryBeamMWAAnaly_cpu(void) {
  test_calc_visi_VarylmnVaryBeam(MWA_ANALY, GAUSSIAN, 0);
}




/*
These tests multiple directions on the sky, varying the position angle, major
ann minor axes
Different beamtypes test how the gains are mapped
*/
void test_calc_visi_gauss_VarylmnVaryPAMajMinFEEBeam_cpu(void) {
  test_calc_visi_VarylmnVaryPAMajMin(FEE_BEAM, GAUSSIAN, 0);
}

void test_calc_visi_gauss_VarylmnVaryPAMajMinAnalyBeam_cpu(void) {
  test_calc_visi_VarylmnVaryPAMajMin(ANALY_DIPOLE, GAUSSIAN, 0);
}

void test_calc_visi_gauss_VarylmnVaryPAMajMinGaussBeam_cpu(void) {
  test_calc_visi_VarylmnVaryPAMajMin(GAUSS_BEAM, GAUSSIAN, 0);
}

void test_calc_visi_gauss_VarylmnVaryPAMajMinNoBeam_cpu(void) {
  test_calc_visi_VarylmnVaryPAMajMin(NO_BEAM, GAUSSIAN, 0);
}

void test_calc_visi_gauss_VarylmnVaryPAMajMinFEEInterpBeam_cpu(void) {
  test_calc_visi_VarylmnVaryPAMajMin(FEE_BEAM_INTERP, GAUSSIAN, 0);
}

void test_calc_visi_gauss_VarylmnVaryPAMajMinMWAAnaly_cpu(void) {
  test_calc_visi_VarylmnVaryPAMajMin(MWA_ANALY, GAUSSIAN, 0);
}

//Run the test with unity
int main(void)
{
    UNITY_BEGIN();
    //Test while varying beam gain for all beam types
    RUN_TEST(test_calc_visi_gauss_VarylmnFEEBeam_cpu);
    RUN_TEST(test_calc_visi_gauss_VarylmnGaussBeam_cpu);
    RUN_TEST(test_calc_visi_gauss_VarylmnAnalyBeam_cpu);
    RUN_TEST(test_calc_visi_gauss_VarylmnNoBeam_cpu);
    RUN_TEST(test_calc_visi_gauss_VarylmnFEEInterpBeam_cpu);
    RUN_TEST(test_calc_visi_gauss_VarylmnMWAAnaly_cpu);

    RUN_TEST(test_calc_visi_gauss_VarylmnVaryFluxFEEBeam_cpu);
    RUN_TEST(test_calc_visi_gauss_VarylmnVaryFluxGaussBeam_cpu);
    RUN_TEST(test_calc_visi_gauss_VarylmnVaryFluxAnalyBeam_cpu);
    RUN_TEST(test_calc_visi_gauss_VarylmnVaryFluxNoBeam_cpu);
    RUN_TEST(test_calc_visi_gauss_VarylmnVaryFluxFEEInterpBeam_cpu);
    RUN_TEST(test_calc_visi_gauss_VarylmnVaryFluxMWAAnaly_cpu);

    RUN_TEST(test_calc_visi_gauss_VarylmnVaryBeamFEEBeam_cpu);
    RUN_TEST(test_calc_visi_gauss_VarylmnVaryBeamGaussBeam_cpu);
    RUN_TEST(test_calc_visi_gauss_VarylmnVaryBeamAnalyBeam_cpu);
    RUN_TEST(test_calc_visi_gauss_VarylmnVaryBeamNoBeam_cpu);
    RUN_TEST(test_calc_visi_gauss_VarylmnVaryBeamFEEInterpBeam_cpu);
    RUN_TEST(test_calc_visi_gauss_VarylmnVaryBeamMWAAnaly_cpu);

    RUN_TEST(test_calc_visi_gauss_VarylmnVaryPAMajMinFEEBeam_cpu);
    RUN_TEST(test_calc_visi_gauss_VarylmnVaryPAMajMinGaussBeam_cpu);
    RUN_TEST(test_calc_visi_gauss_VarylmnVaryPAMajMinAnalyBeam_cpu);
    RUN_TEST(test_calc_visi_gauss_VarylmnVaryPAMajMinNoBeam_cpu);
    RUN_TEST(test_calc_visi_gauss_VarylmnVaryPAMajMinFEEInterpBeam_cpu);
    RUN_TEST(test_calc_visi_gauss_VarylmnVaryPAMajMinMWAAnaly_cpu);

    return UNITY_END();
}
