#include <unity.h>
#include <stdlib.h>
#include <complex.h>

#include "constants.h"
#include "woden_struct_defs.h"
#include "test_kern_calc_visi_common.h"
#include "shapelet_basis.h"

void setUp (void) {} /* Is run before eVary test, put unit init calls here. */
void tearDown (void) {} /* Is run after eVary test, put unit clean-up calls here. */



/*
These tests multiple directions on the sky, with constant fluxes and beam gains
Different beamtypes test how the gains are mapped
*/
void test_kern_calc_visi_shape_VarylmnFEEBeam(void) {
  test_kern_calc_visi_Varylmn(FEE_BEAM, SHAPELET);
}

void test_kern_calc_visi_shape_VarylmnAnalyBeam(void) {
  test_kern_calc_visi_Varylmn(ANALY_DIPOLE, SHAPELET);
}

void test_kern_calc_visi_shape_VarylmnGaussBeam(void) {
  test_kern_calc_visi_Varylmn(GAUSS_BEAM, SHAPELET);
}

void test_kern_calc_visi_shape_VarylmnNoBeam(void) {
  test_kern_calc_visi_Varylmn(NO_BEAM, SHAPELET);
}

void test_kern_calc_visi_shape_VarylmnFEEInterpBeam(void) {
  test_kern_calc_visi_Varylmn(FEE_BEAM_INTERP, SHAPELET);
}

void test_kern_calc_visi_shape_VarylmnMWAAnaly(void) {
  test_kern_calc_visi_Varylmn(MWA_ANALY, SHAPELET);
}

/*
These tests multiple directions on the sky, with varying fluxes and
constant beam gains
Different beamtypes test how the gains are mapped
*/
void test_kern_calc_visi_shape_VarylmnVaryFluxFEEBeam(void) {
  test_kern_calc_visi_VarylmnVaryFlux(FEE_BEAM, SHAPELET);
}

void test_kern_calc_visi_shape_VarylmnVaryFluxAnalyBeam(void) {
  test_kern_calc_visi_VarylmnVaryFlux(ANALY_DIPOLE, SHAPELET);
}

void test_kern_calc_visi_shape_VarylmnVaryFluxGaussBeam(void) {
  test_kern_calc_visi_VarylmnVaryFlux(GAUSS_BEAM, SHAPELET);
}

void test_kern_calc_visi_shape_VarylmnVaryFluxNoBeam(void) {
  test_kern_calc_visi_VarylmnVaryFlux(NO_BEAM, SHAPELET);
}

void test_kern_calc_visi_shape_VarylmnVaryFluxFEEInterpBeam(void) {
  test_kern_calc_visi_VarylmnVaryFlux(FEE_BEAM_INTERP, SHAPELET);
}

void test_kern_calc_visi_shape_VarylmnVaryFluxMWAAnaly(void) {
  test_kern_calc_visi_VarylmnVaryFlux(MWA_ANALY, SHAPELET);
}

/*
These tests multiple directions on the sky, with varying beam gains and
constant fluxe
Different beamtypes test how the gains are mapped
*/
void test_kern_calc_visi_shape_VarylmnVaryBeamFEEBeam(void) {
  test_kern_calc_visi_VarylmnVaryBeam(FEE_BEAM, SHAPELET);
}

void test_kern_calc_visi_shape_VarylmnVaryBeamAnalyBeam(void) {
  test_kern_calc_visi_VarylmnVaryBeam(ANALY_DIPOLE, SHAPELET);
}

void test_kern_calc_visi_shape_VarylmnVaryBeamGaussBeam(void) {
  test_kern_calc_visi_VarylmnVaryBeam(GAUSS_BEAM, SHAPELET);
}

void test_kern_calc_visi_shape_VarylmnVaryBeamNoBeam(void) {
  test_kern_calc_visi_VarylmnVaryBeam(NO_BEAM, SHAPELET);
}

void test_kern_calc_visi_shape_VarylmnVaryBeamFEEInterpBeam(void) {
  test_kern_calc_visi_VarylmnVaryBeam(FEE_BEAM_INTERP, SHAPELET);
}

void test_kern_calc_visi_shape_VarylmnVaryBeamMWAAnaly(void) {
  test_kern_calc_visi_VarylmnVaryBeam(MWA_ANALY, SHAPELET);
}




/*
These tests multiple directions on the sky, varying the position angle, major
ann minor axes
Different beamtypes test how the gains are mapped
*/
void test_kern_calc_visi_shape_VarylmnVaryPAMajMinFEEBeam(void) {
  test_kern_calc_visi_VarylmnVaryPAMajMin(FEE_BEAM, SHAPELET);
}

void test_kern_calc_visi_shape_VarylmnVaryPAMajMinAnalyBeam(void) {
  test_kern_calc_visi_VarylmnVaryPAMajMin(ANALY_DIPOLE, SHAPELET);
}

void test_kern_calc_visi_shape_VarylmnVaryPAMajMinGaussBeam(void) {
  test_kern_calc_visi_VarylmnVaryPAMajMin(GAUSS_BEAM, SHAPELET);
}

void test_kern_calc_visi_shape_VarylmnVaryPAMajMinNoBeam(void) {
  test_kern_calc_visi_VarylmnVaryPAMajMin(NO_BEAM, SHAPELET);
}

void test_kern_calc_visi_shape_VarylmnVaryPAMajMinFEEInterpBeam(void) {
  test_kern_calc_visi_VarylmnVaryPAMajMin(FEE_BEAM_INTERP, SHAPELET);
}

void test_kern_calc_visi_shape_VarylmnVaryPAMajMinMWAAnaly(void) {
  test_kern_calc_visi_VarylmnVaryPAMajMin(MWA_ANALY, SHAPELET);
}

//This test varies the shapelet coeff params
void test_kern_calc_visi_shape_VarylmnMultipleCoeff(int beamtype) {

  int num_baselines = 10.0;
  int num_times = 5.0;
  int num_freqs = 3.0;

  int num_visis = num_baselines*num_times*num_freqs;

  int num_components = 25;
  int n_powers = num_components;
  int n_curves = 0;
  int n_lists = 0;
  int num_coeffs_per_component = 3;

  int num_coeffs = num_coeffs_per_component*num_components;

  //Container for many arrays to feed the GPU
  args_for_testing_t *args_ft = malloc(sizeof(args_for_testing_t));
  //Component information
  components_t components;
  //Allocate memory
  malloc_args_for_testing(args_ft, &components, num_baselines, num_times,
                          num_freqs, num_components, n_powers, n_curves, n_lists,
                          num_coeffs, SHAPELET);

  //Setup l,m args that span a decent chunk of sky
  create_lmn(components);

  int num_beam_values = num_freqs*num_times*num_components;

  //Stick the gains to one everywhere
  for (int visi = 0; visi < num_beam_values; visi++) {
    args_ft->primay_beam_J00[visi] = 1.0 + I*0.0;
    args_ft->primay_beam_J11[visi] = 1.0 + I*0.0;

    if (beamtype == FEE_BEAM || beamtype == FEE_BEAM_INTERP || beamtype == MWA_ANALY ) {
      args_ft->primay_beam_J01[visi] = 1.0 + I*0.0;
      args_ft->primay_beam_J10[visi] = 1.0 + I*0.0;
    }
    else {
      args_ft->primay_beam_J01[visi] = 0.0 + I*0.0;
      args_ft->primay_beam_J10[visi] = 0.0 + I*0.0;
    }
  }

  //Just stick Stokes I to 1.0, SI to zero, and reference freqs to 150MHz
  for (int comp = 0; comp < num_components; comp++) {
    components.power_ref_stokesI[comp] = 1.0;
    // components.power_ref_stokesQ[comp] = 0.0;
    // components.power_ref_stokesU[comp] = 0.0;
    // components.power_ref_stokesV[comp] = 0.0;
    components.power_SIs[comp] = 0.0;
    components.power_ref_freqs[comp] = 150e+6;

    //Set major,minor to 3 arcmins
    components.pas[comp] = 0.0;
    components.majors[comp] = 3.0*(DD2R / 60.0);
    components.minors[comp] = 3.0*(DD2R / 60.0);

    components.power_comp_inds[comp] = comp;
  }

  //Make up some u,v,w values and scale by wavelength in correct order
  setup_uvw_and_freqs(args_ft, num_times, num_freqs, num_baselines);

  //Set the shapelet u,v,w same as the measurement equation one (this is not
  //true in reality but works fine for testing)

  int count = 0;

  for (int comp_step = 0; comp_step < num_components; comp_step++) {
    for ( int time_step = 0; time_step < num_times; time_step++ ) {
      for (int baseline = 0; baseline < num_baselines; baseline++) {
        args_ft->u_shapes[count] = ((baseline + 1)*10);
        args_ft->v_shapes[count] = ((baseline + 1)*10);

        count ++;
      }
    }
  }

  //Stick a number of coeffs in per component
  user_precision_t sign;
  count = 0;
  for (int comp = 0; comp < num_components; comp++) {
    for (int coeff = 0; coeff < num_coeffs_per_component; coeff++) {

      if (count % 2 == 0) {
        sign = 1.0;
      } else {
        sign = -1.0;
      }

      components.n1s[count] = count;
      components.n2s[count] = count + 1;
      components.shape_coeffs[count] = sign*1e-3*(coeff + 1);
      components.param_indexes[count] = comp;

      count ++;
    }
  }
  components.do_QUV = 0;
  test_kern_calc_visi_all(n_powers, n_curves, n_lists, num_baselines, num_coeffs,
          num_freqs, num_visis, num_times, beamtype, SHAPELET,
          components, args_ft->extrap_freqs,
          args_ft->us, args_ft->vs, args_ft->ws,
          args_ft->u_shapes, args_ft->v_shapes,
          args_ft->sum_visi_XX_real, args_ft->sum_visi_XX_imag,
          args_ft->sum_visi_XY_real, args_ft->sum_visi_XY_imag,
          args_ft->sum_visi_YX_real, args_ft->sum_visi_YX_imag,
          args_ft->sum_visi_YY_real, args_ft->sum_visi_YY_imag,
          args_ft->allsteps_wavelengths, args_ft->sbf,
          args_ft->primay_beam_J00, args_ft->primay_beam_J01,
          args_ft->primay_beam_J10, args_ft->primay_beam_J11);
  //
  // //Check all results are within 0.1% of expected value
  // // double frac_tol = 1e-3;
  test_visi_outputs(num_visis, n_powers, n_curves, n_lists,
                    num_baselines, num_freqs, args_ft->extrap_freqs,
                    beamtype, args_ft, components, SHAPELET);

  free_args_for_testing( args_ft, components, SHAPELET );
}

/*
Check varying shapelet coeff params with various beam models
*/
void test_kern_calc_visi_shape_VarylmnMultipleCoeffFEEBeam(void) {
  test_kern_calc_visi_shape_VarylmnMultipleCoeff(FEE_BEAM);
}

void test_kern_calc_visi_shape_VarylmnMultipleCoeffAnalyBeam(void) {
  test_kern_calc_visi_shape_VarylmnMultipleCoeff(ANALY_DIPOLE);
}

void test_kern_calc_visi_shape_VarylmnMultipleCoeffGaussBeam(void) {
  test_kern_calc_visi_shape_VarylmnMultipleCoeff(GAUSS_BEAM);
}

void test_kern_calc_visi_shape_VarylmnMultipleCoeffNoBeam(void) {
  test_kern_calc_visi_shape_VarylmnMultipleCoeff(NO_BEAM);
}

void test_kern_calc_visi_shape_VarylmnMultipleCoeffFEEInterpBeam(void) {
  test_kern_calc_visi_shape_VarylmnMultipleCoeff(FEE_BEAM_INTERP);
}

void test_kern_calc_visi_shape_VarylmnMultipleCoeffMWAAnaly(void) {
  test_kern_calc_visi_shape_VarylmnMultipleCoeff(MWA_ANALY);
}

//Run the test with unity
int main(void)
{
    UNITY_BEGIN();
    //Test while varying beam gain for all beam types
    RUN_TEST(test_kern_calc_visi_shape_VarylmnFEEBeam);
    RUN_TEST(test_kern_calc_visi_shape_VarylmnGaussBeam);
    RUN_TEST(test_kern_calc_visi_shape_VarylmnAnalyBeam);
    RUN_TEST(test_kern_calc_visi_shape_VarylmnNoBeam);
    RUN_TEST(test_kern_calc_visi_shape_VarylmnFEEInterpBeam);
    RUN_TEST(test_kern_calc_visi_shape_VarylmnMWAAnaly);

    RUN_TEST(test_kern_calc_visi_shape_VarylmnVaryFluxNoBeam);
    RUN_TEST(test_kern_calc_visi_shape_VarylmnVaryFluxFEEBeam);
    RUN_TEST(test_kern_calc_visi_shape_VarylmnVaryFluxGaussBeam);
    RUN_TEST(test_kern_calc_visi_shape_VarylmnVaryFluxAnalyBeam);
    RUN_TEST(test_kern_calc_visi_shape_VarylmnVaryFluxFEEInterpBeam);
    RUN_TEST(test_kern_calc_visi_shape_VarylmnVaryFluxMWAAnaly);

    RUN_TEST(test_kern_calc_visi_shape_VarylmnVaryBeamFEEBeam);
    RUN_TEST(test_kern_calc_visi_shape_VarylmnVaryBeamGaussBeam);
    RUN_TEST(test_kern_calc_visi_shape_VarylmnVaryBeamAnalyBeam);
    RUN_TEST(test_kern_calc_visi_shape_VarylmnVaryBeamNoBeam);
    RUN_TEST(test_kern_calc_visi_shape_VarylmnVaryBeamFEEInterpBeam);
    RUN_TEST(test_kern_calc_visi_shape_VarylmnVaryBeamMWAAnaly);

    RUN_TEST(test_kern_calc_visi_shape_VarylmnVaryPAMajMinFEEBeam);
    RUN_TEST(test_kern_calc_visi_shape_VarylmnVaryPAMajMinGaussBeam);
    RUN_TEST(test_kern_calc_visi_shape_VarylmnVaryPAMajMinAnalyBeam);
    RUN_TEST(test_kern_calc_visi_shape_VarylmnVaryPAMajMinNoBeam);
    RUN_TEST(test_kern_calc_visi_shape_VarylmnVaryPAMajMinFEEInterpBeam);
    RUN_TEST(test_kern_calc_visi_shape_VarylmnVaryPAMajMinMWAAnaly);

    RUN_TEST(test_kern_calc_visi_shape_VarylmnMultipleCoeffFEEBeam);
    RUN_TEST(test_kern_calc_visi_shape_VarylmnMultipleCoeffGaussBeam);
    RUN_TEST(test_kern_calc_visi_shape_VarylmnMultipleCoeffAnalyBeam);
    RUN_TEST(test_kern_calc_visi_shape_VarylmnMultipleCoeffNoBeam);
    RUN_TEST(test_kern_calc_visi_shape_VarylmnMultipleCoeffFEEInterpBeam);
    RUN_TEST(test_kern_calc_visi_shape_VarylmnMultipleCoeffMWAAnaly);

    return UNITY_END();
}
