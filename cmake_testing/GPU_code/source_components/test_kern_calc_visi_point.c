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
void test_kern_calc_visi_point_VarylmnFEEBeam(void) {
  test_kern_calc_visi_Varylmn(FEE_BEAM, POINT);
}

void test_kern_calc_visi_point_VarylmnAnalyBeam(void) {
  test_kern_calc_visi_Varylmn(ANALY_DIPOLE, POINT);
}

void test_kern_calc_visi_point_VarylmnGaussBeam(void) {
  test_kern_calc_visi_Varylmn(GAUSS_BEAM, POINT);
}

void test_kern_calc_visi_point_VarylmnNoBeam(void) {
  test_kern_calc_visi_Varylmn(NO_BEAM, POINT);
}

void test_kern_calc_visi_point_VarylmnFEEInterpBeam(void) {
  test_kern_calc_visi_Varylmn(FEE_BEAM_INTERP, POINT);
}

void test_kern_calc_visi_point_VarylmnMWAAnaly(void) {
  test_kern_calc_visi_Varylmn(MWA_ANALY, POINT);
}

/*
These tests multiple directions on the sky, with varying fluxes and
constant beam gains
Different beamtypes test how the gains are mapped
*/
void test_kern_calc_visi_point_VarylmnVaryFluxFEEBeam(void) {
  test_kern_calc_visi_VarylmnVaryFlux(FEE_BEAM, POINT);
}

void test_kern_calc_visi_point_VarylmnVaryFluxAnalyBeam(void) {
  test_kern_calc_visi_VarylmnVaryFlux(ANALY_DIPOLE, POINT);
}

void test_kern_calc_visi_point_VarylmnVaryFluxGaussBeam(void) {
  test_kern_calc_visi_VarylmnVaryFlux(GAUSS_BEAM, POINT);
}

void test_kern_calc_visi_point_VarylmnVaryFluxNoBeam(void) {
  test_kern_calc_visi_VarylmnVaryFlux(NO_BEAM, POINT);
}

void test_kern_calc_visi_point_VarylmnVaryFluxFEEInterpBeam(void) {
  test_kern_calc_visi_VarylmnVaryFlux(FEE_BEAM_INTERP, POINT);
}

void test_kern_calc_visi_point_VarylmnVaryFluxMWAAnaly(void) {
  test_kern_calc_visi_VarylmnVaryFlux(MWA_ANALY, POINT);
}

/*
These tests multiple directions on the sky, with varying beam gains and
constant fluxe
Different beamtypes test how the gains are mapped
*/
void test_kern_calc_visi_point_VarylmnVaryBeamFEEBeam(void) {
  test_kern_calc_visi_VarylmnVaryBeam(FEE_BEAM, POINT);
}

void test_kern_calc_visi_point_VarylmnVaryBeamAnalyBeam(void) {
  test_kern_calc_visi_VarylmnVaryBeam(ANALY_DIPOLE, POINT);
}

void test_kern_calc_visi_point_VarylmnVaryBeamGaussBeam(void) {
  test_kern_calc_visi_VarylmnVaryBeam(GAUSS_BEAM, POINT);
}

void test_kern_calc_visi_point_VarylmnVaryBeamNoBeam(void) {
  test_kern_calc_visi_VarylmnVaryBeam(NO_BEAM, POINT);
}

void test_kern_calc_visi_point_VarylmnVaryBeamFEEInterpBeam(void) {
  test_kern_calc_visi_VarylmnVaryBeam(FEE_BEAM_INTERP, POINT);
}

void test_kern_calc_visi_point_VarylmnVaryBeamMWAAnaly(void) {
  test_kern_calc_visi_VarylmnVaryBeam(MWA_ANALY, POINT);
}

//Run the test with unity
int main(void)
{
    UNITY_BEGIN();
    //Test while varying beam gain for all beam types
    RUN_TEST(test_kern_calc_visi_point_VarylmnFEEBeam);
    RUN_TEST(test_kern_calc_visi_point_VarylmnGaussBeam);
    RUN_TEST(test_kern_calc_visi_point_VarylmnAnalyBeam);
    RUN_TEST(test_kern_calc_visi_point_VarylmnNoBeam);
    RUN_TEST(test_kern_calc_visi_point_VarylmnFEEInterpBeam);
    RUN_TEST(test_kern_calc_visi_point_VarylmnMWAAnaly);

    RUN_TEST(test_kern_calc_visi_point_VarylmnVaryFluxFEEBeam);
    RUN_TEST(test_kern_calc_visi_point_VarylmnVaryFluxGaussBeam);
    RUN_TEST(test_kern_calc_visi_point_VarylmnVaryFluxAnalyBeam);
    RUN_TEST(test_kern_calc_visi_point_VarylmnVaryFluxNoBeam);
    RUN_TEST(test_kern_calc_visi_point_VarylmnVaryFluxFEEInterpBeam);
    RUN_TEST(test_kern_calc_visi_point_VarylmnVaryFluxMWAAnaly);

    RUN_TEST(test_kern_calc_visi_point_VarylmnVaryBeamFEEBeam);
    RUN_TEST(test_kern_calc_visi_point_VarylmnVaryBeamGaussBeam);
    RUN_TEST(test_kern_calc_visi_point_VarylmnVaryBeamAnalyBeam);
    RUN_TEST(test_kern_calc_visi_point_VarylmnVaryBeamNoBeam);
    RUN_TEST(test_kern_calc_visi_point_VarylmnVaryBeamFEEInterpBeam);
    RUN_TEST(test_kern_calc_visi_point_VarylmnVaryBeamMWAAnaly);

    return UNITY_END();
}
