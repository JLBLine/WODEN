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
void test_calc_visi_shape_VarylmnFEEBeam_cpu(void) {
  test_calc_visi_Varylmn(FEE_BEAM, SHAPELET, 0);
}

void test_calc_visi_shape_VarylmnAnalyBeam_cpu(void) {
  test_calc_visi_Varylmn(ANALY_DIPOLE, SHAPELET, 0);
}

void test_calc_visi_shape_VarylmnGaussBeam_cpu(void) {
  test_calc_visi_Varylmn(GAUSS_BEAM, SHAPELET, 0);
}

void test_calc_visi_shape_VarylmnNoBeam_cpu(void) {
  test_calc_visi_Varylmn(NO_BEAM, SHAPELET, 0);
}

void test_calc_visi_shape_VarylmnFEEInterpBeam_cpu(void) {
  test_calc_visi_Varylmn(FEE_BEAM_INTERP, SHAPELET, 0);
}

void test_calc_visi_shape_VarylmnMWAAnaly_cpu(void) {
  test_calc_visi_Varylmn(MWA_ANALY, SHAPELET, 0);
}

/*
These tests multiple directions on the sky, with varying fluxes and
constant beam gains
Different beamtypes test how the gains are mapped
*/
void test_calc_visi_shape_VarylmnVaryFluxFEEBeam_cpu(void) {
  test_calc_visi_VarylmnVaryFlux(FEE_BEAM, SHAPELET, 0);
}

void test_calc_visi_shape_VarylmnVaryFluxAnalyBeam_cpu(void) {
  test_calc_visi_VarylmnVaryFlux(ANALY_DIPOLE, SHAPELET, 0);
}

void test_calc_visi_shape_VarylmnVaryFluxGaussBeam_cpu(void) {
  test_calc_visi_VarylmnVaryFlux(GAUSS_BEAM, SHAPELET, 0);
}

void test_calc_visi_shape_VarylmnVaryFluxNoBeam_cpu(void) {
  test_calc_visi_VarylmnVaryFlux(NO_BEAM, SHAPELET, 0);
}

void test_calc_visi_shape_VarylmnVaryFluxFEEInterpBeam_cpu(void) {
  test_calc_visi_VarylmnVaryFlux(FEE_BEAM_INTERP, SHAPELET, 0);
}

void test_calc_visi_shape_VarylmnVaryFluxMWAAnaly_cpu(void) {
  test_calc_visi_VarylmnVaryFlux(MWA_ANALY, SHAPELET, 0);
}

/*
These tests multiple directions on the sky, with varying beam gains and
constant fluxe
Different beamtypes test how the gains are mapped
*/
void test_calc_visi_shape_VarylmnVaryBeamFEEBeam_cpu(void) {
  test_calc_visi_VarylmnVaryBeam(FEE_BEAM, SHAPELET, 0);
}

void test_calc_visi_shape_VarylmnVaryBeamAnalyBeam_cpu(void) {
  test_calc_visi_VarylmnVaryBeam(ANALY_DIPOLE, SHAPELET, 0);
}

void test_calc_visi_shape_VarylmnVaryBeamGaussBeam_cpu(void) {
  test_calc_visi_VarylmnVaryBeam(GAUSS_BEAM, SHAPELET, 0);
}

void test_calc_visi_shape_VarylmnVaryBeamNoBeam_cpu(void) {
  test_calc_visi_VarylmnVaryBeam(NO_BEAM, SHAPELET, 0);
}

void test_calc_visi_shape_VarylmnVaryBeamFEEInterpBeam_cpu(void) {
  test_calc_visi_VarylmnVaryBeam(FEE_BEAM_INTERP, SHAPELET, 0);
}

void test_calc_visi_shape_VarylmnVaryBeamMWAAnaly_cpu(void) {
  test_calc_visi_VarylmnVaryBeam(MWA_ANALY, SHAPELET, 0);
}




/*
These tests multiple directions on the sky, varying the position angle, major
ann minor axes
Different beamtypes test how the gains are mapped
*/
void test_calc_visi_shape_VarylmnVaryPAMajMinFEEBeam_cpu(void) {
  test_calc_visi_VarylmnVaryPAMajMin(FEE_BEAM, SHAPELET, 0);
}

void test_calc_visi_shape_VarylmnVaryPAMajMinAnalyBeam_cpu(void) {
  test_calc_visi_VarylmnVaryPAMajMin(ANALY_DIPOLE, SHAPELET, 0);
}

void test_calc_visi_shape_VarylmnVaryPAMajMinGaussBeam_cpu(void) {
  test_calc_visi_VarylmnVaryPAMajMin(GAUSS_BEAM, SHAPELET, 0);
}

void test_calc_visi_shape_VarylmnVaryPAMajMinNoBeam_cpu(void) {
  test_calc_visi_VarylmnVaryPAMajMin(NO_BEAM, SHAPELET, 0);
}

void test_calc_visi_shape_VarylmnVaryPAMajMinFEEInterpBeam_cpu(void) {
  test_calc_visi_VarylmnVaryPAMajMin(FEE_BEAM_INTERP, SHAPELET, 0);
}

void test_calc_visi_shape_VarylmnVaryPAMajMinMWAAnaly_cpu(void) {
  test_calc_visi_VarylmnVaryPAMajMin(MWA_ANALY, SHAPELET, 0);
}


/*
Check varying shapelet coeff params with various beam models
*/
void test_calc_visi_shape_VarylmnMultipleCoeffFEEBeam_cpu(void) {
  test_calc_visi_shape_VarylmnMultipleCoeff(FEE_BEAM, 0);
}

void test_calc_visi_shape_VarylmnMultipleCoeffAnalyBeam_cpu(void) {
  test_calc_visi_shape_VarylmnMultipleCoeff(ANALY_DIPOLE, 0);
}

void test_calc_visi_shape_VarylmnMultipleCoeffGaussBeam_cpu(void) {
  test_calc_visi_shape_VarylmnMultipleCoeff(GAUSS_BEAM, 0);
}

void test_calc_visi_shape_VarylmnMultipleCoeffNoBeam_cpu(void) {
  test_calc_visi_shape_VarylmnMultipleCoeff(NO_BEAM, 0);
}

void test_calc_visi_shape_VarylmnMultipleCoeffFEEInterpBeam_cpu(void) {
  test_calc_visi_shape_VarylmnMultipleCoeff(FEE_BEAM_INTERP, 0);
}

void test_calc_visi_shape_VarylmnMultipleCoeffMWAAnaly_cpu(void) {
  test_calc_visi_shape_VarylmnMultipleCoeff(MWA_ANALY, 0);
}

//Run the test with unity
int main(void)
{
    UNITY_BEGIN();
    //Test while varying beam gain for all beam types
    RUN_TEST(test_calc_visi_shape_VarylmnFEEBeam_cpu);
    RUN_TEST(test_calc_visi_shape_VarylmnGaussBeam_cpu);
    RUN_TEST(test_calc_visi_shape_VarylmnAnalyBeam_cpu);
    RUN_TEST(test_calc_visi_shape_VarylmnNoBeam_cpu);
    RUN_TEST(test_calc_visi_shape_VarylmnFEEInterpBeam_cpu);
    RUN_TEST(test_calc_visi_shape_VarylmnMWAAnaly_cpu);

    RUN_TEST(test_calc_visi_shape_VarylmnVaryFluxNoBeam_cpu);
    RUN_TEST(test_calc_visi_shape_VarylmnVaryFluxFEEBeam_cpu);
    RUN_TEST(test_calc_visi_shape_VarylmnVaryFluxGaussBeam_cpu);
    RUN_TEST(test_calc_visi_shape_VarylmnVaryFluxAnalyBeam_cpu);
    RUN_TEST(test_calc_visi_shape_VarylmnVaryFluxFEEInterpBeam_cpu);
    RUN_TEST(test_calc_visi_shape_VarylmnVaryFluxMWAAnaly_cpu);

    RUN_TEST(test_calc_visi_shape_VarylmnVaryBeamFEEBeam_cpu);
    RUN_TEST(test_calc_visi_shape_VarylmnVaryBeamGaussBeam_cpu);
    RUN_TEST(test_calc_visi_shape_VarylmnVaryBeamAnalyBeam_cpu);
    RUN_TEST(test_calc_visi_shape_VarylmnVaryBeamNoBeam_cpu);
    RUN_TEST(test_calc_visi_shape_VarylmnVaryBeamFEEInterpBeam_cpu);
    RUN_TEST(test_calc_visi_shape_VarylmnVaryBeamMWAAnaly_cpu);

    RUN_TEST(test_calc_visi_shape_VarylmnVaryPAMajMinFEEBeam_cpu);
    RUN_TEST(test_calc_visi_shape_VarylmnVaryPAMajMinGaussBeam_cpu);
    RUN_TEST(test_calc_visi_shape_VarylmnVaryPAMajMinAnalyBeam_cpu);
    RUN_TEST(test_calc_visi_shape_VarylmnVaryPAMajMinNoBeam_cpu);
    RUN_TEST(test_calc_visi_shape_VarylmnVaryPAMajMinFEEInterpBeam_cpu);
    RUN_TEST(test_calc_visi_shape_VarylmnVaryPAMajMinMWAAnaly_cpu);

    RUN_TEST(test_calc_visi_shape_VarylmnMultipleCoeffFEEBeam_cpu);
    RUN_TEST(test_calc_visi_shape_VarylmnMultipleCoeffGaussBeam_cpu);
    RUN_TEST(test_calc_visi_shape_VarylmnMultipleCoeffAnalyBeam_cpu);
    RUN_TEST(test_calc_visi_shape_VarylmnMultipleCoeffNoBeam_cpu);
    RUN_TEST(test_calc_visi_shape_VarylmnMultipleCoeffFEEInterpBeam_cpu);
    RUN_TEST(test_calc_visi_shape_VarylmnMultipleCoeffMWAAnaly_cpu);

    return UNITY_END();
}
