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
void test_calc_visi_point_VarylmnFEEBeam_gpu(void) {
  test_calc_visi_Varylmn(FEE_BEAM, POINT, 1);
}

void test_calc_visi_point_VarylmnAnalyBeam_gpu(void) {
  test_calc_visi_Varylmn(ANALY_DIPOLE, POINT, 1);
}

void test_calc_visi_point_VarylmnGaussBeam_gpu(void) {
  test_calc_visi_Varylmn(GAUSS_BEAM, POINT, 1);
}

void test_calc_visi_point_VarylmnNoBeam_gpu(void) {
  test_calc_visi_Varylmn(NO_BEAM, POINT, 1);
}

void test_calc_visi_point_VarylmnFEEInterpBeam_gpu(void) {
  test_calc_visi_Varylmn(FEE_BEAM_INTERP, POINT, 1);
}

void test_calc_visi_point_VarylmnMWAAnaly_gpu(void) {
  test_calc_visi_Varylmn(MWA_ANALY, POINT, 1);
}

/*
These tests multiple directions on the sky, with varying fluxes and
constant beam gains
Different beamtypes test how the gains are mapped
*/
void test_calc_visi_point_VarylmnVaryFluxFEEBeam_gpu(void) {
  test_calc_visi_VarylmnVaryFlux(FEE_BEAM, POINT, 1);
}

void test_calc_visi_point_VarylmnVaryFluxAnalyBeam_gpu(void) {
  test_calc_visi_VarylmnVaryFlux(ANALY_DIPOLE, POINT, 1);
}

void test_calc_visi_point_VarylmnVaryFluxGaussBeam_gpu(void) {
  test_calc_visi_VarylmnVaryFlux(GAUSS_BEAM, POINT, 1);
}

void test_calc_visi_point_VarylmnVaryFluxNoBeam_gpu(void) {
  test_calc_visi_VarylmnVaryFlux(NO_BEAM, POINT, 1);
}

void test_calc_visi_point_VarylmnVaryFluxFEEInterpBeam_gpu(void) {
  test_calc_visi_VarylmnVaryFlux(FEE_BEAM_INTERP, POINT, 1);
}

void test_calc_visi_point_VarylmnVaryFluxMWAAnaly_gpu(void) {
  test_calc_visi_VarylmnVaryFlux(MWA_ANALY, POINT, 1);
}

/*
These tests multiple directions on the sky, with varying beam gains and
constant fluxe
Different beamtypes test how the gains are mapped
*/
void test_calc_visi_point_VarylmnVaryBeamFEEBeam_gpu(void) {
  test_calc_visi_VarylmnVaryBeam(FEE_BEAM, POINT, 1);
}

void test_calc_visi_point_VarylmnVaryBeamAnalyBeam_gpu(void) {
  test_calc_visi_VarylmnVaryBeam(ANALY_DIPOLE, POINT, 1);
}

void test_calc_visi_point_VarylmnVaryBeamGaussBeam_gpu(void) {
  test_calc_visi_VarylmnVaryBeam(GAUSS_BEAM, POINT, 1);
}

void test_calc_visi_point_VarylmnVaryBeamNoBeam_gpu(void) {
  test_calc_visi_VarylmnVaryBeam(NO_BEAM, POINT, 1);
}

void test_calc_visi_point_VarylmnVaryBeamFEEInterpBeam_gpu(void) {
  test_calc_visi_VarylmnVaryBeam(FEE_BEAM_INTERP, POINT, 1);
}

void test_calc_visi_point_VarylmnVaryBeamMWAAnaly_gpu(void) {
  test_calc_visi_VarylmnVaryBeam(MWA_ANALY, POINT, 1);
}

//Run the test with unity
int main(void)
{
    UNITY_BEGIN();
    //Test while varying beam gain for all beam types
    RUN_TEST(test_calc_visi_point_VarylmnFEEBeam_gpu);
    RUN_TEST(test_calc_visi_point_VarylmnGaussBeam_gpu);
    RUN_TEST(test_calc_visi_point_VarylmnAnalyBeam_gpu);
    RUN_TEST(test_calc_visi_point_VarylmnNoBeam_gpu);
    RUN_TEST(test_calc_visi_point_VarylmnFEEInterpBeam_gpu);
    RUN_TEST(test_calc_visi_point_VarylmnMWAAnaly_gpu);

    RUN_TEST(test_calc_visi_point_VarylmnVaryFluxFEEBeam_gpu);
    RUN_TEST(test_calc_visi_point_VarylmnVaryFluxGaussBeam_gpu);
    RUN_TEST(test_calc_visi_point_VarylmnVaryFluxAnalyBeam_gpu);
    RUN_TEST(test_calc_visi_point_VarylmnVaryFluxNoBeam_gpu);
    RUN_TEST(test_calc_visi_point_VarylmnVaryFluxFEEInterpBeam_gpu);
    RUN_TEST(test_calc_visi_point_VarylmnVaryFluxMWAAnaly_gpu);

    RUN_TEST(test_calc_visi_point_VarylmnVaryBeamFEEBeam_gpu);
    RUN_TEST(test_calc_visi_point_VarylmnVaryBeamGaussBeam_gpu);
    RUN_TEST(test_calc_visi_point_VarylmnVaryBeamAnalyBeam_gpu);
    RUN_TEST(test_calc_visi_point_VarylmnVaryBeamNoBeam_gpu);
    RUN_TEST(test_calc_visi_point_VarylmnVaryBeamFEEInterpBeam_gpu);
    RUN_TEST(test_calc_visi_point_VarylmnVaryBeamMWAAnaly_gpu);

    return UNITY_END();
}
