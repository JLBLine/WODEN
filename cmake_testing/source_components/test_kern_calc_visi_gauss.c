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

//External CUDA code we're linking in
extern void test_kern_calc_visi_gaussian(int num_components, int num_baselines,
          int num_freqs, int num_visis, int num_times, int beamtype,
          user_precision_t *component_freqs,
          user_precision_t *flux_I, user_precision_t *flux_Q, user_precision_t *flux_U, user_precision_t *flux_V,
          user_precision_t *SIs, user_precision_t *us, user_precision_t *vs, user_precision_t *ws,
          user_precision_t *sum_visi_XX_real, user_precision_t *sum_visi_XX_imag,
          user_precision_t *sum_visi_XY_real, user_precision_t *sum_visi_XY_imag,
          user_precision_t *sum_visi_YX_real, user_precision_t *sum_visi_YX_imag,
          user_precision_t *sum_visi_YY_real, user_precision_t *sum_visi_YY_imag,
          user_precision_t *allsteps_wavelengths,
          double *ls, double *ms, double *ns,
          user_precision_t *gauss_pas, user_precision_t *gauss_majors, user_precision_t *gauss_minors,
          user_precision_complex_t *primay_beam_J00, user_precision_complex_t *primay_beam_J01,
          user_precision_complex_t *primay_beam_J10, user_precision_complex_t *primay_beam_J11);

#define UNITY_INCLUDE_FLOAT

//Change required accuracy of outputs for different precisions
#ifdef DOUBLE_PRECISION
  //Accurate to within 0.00000000001%
  #define FRAC_TOL 1e-12
#else
  //Accurate to within 0.001%
  #define FRAC_TOL 1e-5
#endif

/*
Test the __global__ code that calculates visibilities for gauss sources
Vary the l,m,n coords but keep all other variables constant
*/
void test_kern_calc_visi_gauss_Varylmn(int beamtype) {

  int num_baselines = 10.0;
  int num_times = 5.0;
  int num_freqs = 3.0;

  int num_visis = num_baselines*num_times*num_freqs;

  int num_components = 25;

  //Container for many arrays to feed the GPU
  args_for_testing_t *args_ft = malloc(sizeof(args_for_testing_t));
  //Allocate memory
  int num_coeffs = 0; //Only applies for SHAPELET
  malloc_args_for_testing(args_ft, num_baselines, num_times,
                          num_freqs, num_components, num_coeffs, GAUSSIAN);

  //Setup l,m args that span a decent chunk of sky
  create_lmn(args_ft);

  int num_beam_values = num_freqs*num_times*num_components;

  //Stick the gains to one everywhere
  for (size_t visi = 0; visi < num_beam_values; visi++) {
    args_ft->primay_beam_J00[visi] = 1.0 + I*0.0;
    args_ft->primay_beam_J01[visi] = 1.0 + I*0.0;
    args_ft->primay_beam_J10[visi] = 1.0 + I*0.0;
    args_ft->primay_beam_J11[visi] = 1.0 + I*0.0;
  }

  //Just stick Stokes I to 1.0, SI to zero, and reference freqs to 150MHz
  for (size_t comp = 0; comp < num_components; comp++) {
    args_ft->flux_I[comp] = 1.0;
    args_ft->SIs[comp] = 0.0;
    args_ft->component_freqs[comp] = 150e+6;
    args_ft->pas[comp] = 0.0;
    //Set major,minor to 3 arcmins
    args_ft->majors[comp] = 3.0*(DD2R / 60.0);
    args_ft->minors[comp] = 3.0*(DD2R / 60.0);
  }

  //Make up some u,v,w values and scale by wavelength in correct order
  int count = 0;
  user_precision_t freq_base = 150e+6;
  user_precision_t freq_inc = 25e+6;
  user_precision_t wavelength, frequency;
  for ( int time_step = 0; time_step < num_times; time_step++ ) {
    for (int freq_step = 0; freq_step < num_freqs; freq_step++) {
      frequency = freq_base + freq_step*freq_inc;
      wavelength = VELC / frequency;
      for (int baseline = 0; baseline < num_baselines; baseline++) {
        args_ft->us[count] = ((baseline + 1)*100) / wavelength;
        args_ft->vs[count] = ((baseline + 1)*100) / wavelength;
        //ws are usually smaller than u,v
        args_ft->ws[count] = ((baseline + 1)*10) / wavelength;

        args_ft->allsteps_wavelengths[count] = wavelength;

        count ++;
      }
    }
  }

  //Run the CUDA code
  test_kern_calc_visi_gaussian(num_components, num_baselines,
          num_freqs, num_visis, num_times, beamtype,
          args_ft->component_freqs,
          args_ft->flux_I, args_ft->flux_Q, args_ft->flux_U, args_ft->flux_V,
          args_ft->SIs, args_ft->us, args_ft->vs, args_ft->ws,
          args_ft->sum_visi_XX_real, args_ft->sum_visi_XX_imag,
          args_ft->sum_visi_XY_real, args_ft->sum_visi_XY_imag,
          args_ft->sum_visi_YX_real, args_ft->sum_visi_YX_imag,
          args_ft->sum_visi_YY_real, args_ft->sum_visi_YY_imag,
          args_ft->allsteps_wavelengths,
          args_ft->ls, args_ft->ms, args_ft->ns,
          args_ft->pas, args_ft->majors, args_ft->minors,
          args_ft->primay_beam_J00, args_ft->primay_beam_J01,
          args_ft->primay_beam_J10, args_ft->primay_beam_J11);

  //Check all results are within 0.1% of expected value
  test_visi_outputs(num_visis, num_components, num_baselines, num_freqs,
                    FRAC_TOL, beamtype, args_ft, GAUSSIAN);

  free_args_for_testing( args_ft, GAUSSIAN );
}

/*
This test checks varying the measurement equation with beamtype=FEE_BEAM
*/
void test_kern_calc_visi_gauss_VarylmnFEEBeam(void) {
  test_kern_calc_visi_gauss_Varylmn(FEE_BEAM);
}

/*
This test checks varying the measurement equation with beamtype=ANALY_DIPOLE
*/
void test_kern_calc_visi_gauss_VarylmnAnalyBeam(void) {
  test_kern_calc_visi_gauss_Varylmn(ANALY_DIPOLE);
}

/*
This test checks varying the measurement equation with beamtype=GAUSS_BEAM
*/
void test_kern_calc_visi_gauss_VarylmnGaussBeam(void) {
  test_kern_calc_visi_gauss_Varylmn(GAUSS_BEAM);
}

/*
This test checks varying the measurement equation with beamtype=NO_BEAM
*/
void test_kern_calc_visi_gauss_VarylmnNoBeam(void) {
  test_kern_calc_visi_gauss_Varylmn(NO_BEAM);
}

/*
Test the __device__ code that updates the summed visibilities by grabbing the
correct beam gain and mesurement equation, multiplying and summing onto the visi
Here we keep the component visibilities and fluxes constant and vary the beam gains
Test works for all primary beam types
*/
void test_kern_calc_visi_gauss_VarylmnVaryFlux(int beamtype) {

  int num_baselines = 10.0;
  int num_times = 5.0;
  int num_freqs = 3.0;

  int num_visis = num_baselines*num_times*num_freqs;

  int num_components = 25;

  //Container for many arrays to feed the GPU
  args_for_testing_t *args_ft = malloc(sizeof(args_for_testing_t));
  //Allocate memory
  int num_coeffs = 0; //Only applies for SHAPELET
  malloc_args_for_testing(args_ft, num_baselines, num_times,
                          num_freqs, num_components, num_coeffs, GAUSSIAN);

  //Setup l,m args that span a decent chunk of sky
  create_lmn(args_ft);

  int num_beam_values = num_freqs*num_times*num_components;

  //Stick the gains to one everywhere
  for (size_t visi = 0; visi < num_beam_values; visi++) {
    args_ft->primay_beam_J00[visi] = 1.0 + I*0.0;
    args_ft->primay_beam_J01[visi] = 1.0 + I*0.0;
    args_ft->primay_beam_J10[visi] = 1.0 + I*0.0;
    args_ft->primay_beam_J11[visi] = 1.0 + I*0.0;
  }

  //Just stick Stokes I to 1.0, SI to zero, and reference freqs to 150MHz
  for (size_t comp = 0; comp < num_components; comp++) {
    args_ft->flux_I[comp] = comp;
    args_ft->SIs[comp] = -0.8;
    args_ft->component_freqs[comp] = 150e+6;
    args_ft->pas[comp] = 0.0;
    //Set major,minor to 3 arcmins
    args_ft->majors[comp] = 3.0*(DD2R / 60.0);
    args_ft->minors[comp] = 3.0*(DD2R / 60.0);
  }

  //Make up some u,v,w values and scale by wavelength in correct order
  int count = 0;
  user_precision_t freq_base = 150e+6;
  user_precision_t freq_inc = 25e+6;
  user_precision_t wavelength, frequency;

  for ( int time_step = 0; time_step < num_times; time_step++ ) {
    for (int freq_step = 0; freq_step < num_freqs; freq_step++) {
      frequency = freq_base + freq_step*freq_inc;
      wavelength = VELC / frequency;

      for (int baseline = 0; baseline < num_baselines; baseline++) {
        args_ft->us[count] = ((baseline + 1)*100) / wavelength;
        args_ft->vs[count] = ((baseline + 1)*100) / wavelength;
        //ws are usually smaller than u,v
        args_ft->ws[count] = ((baseline + 1)*10) / wavelength;

        args_ft->allsteps_wavelengths[count] = wavelength;

        count ++;
      }
    }
  }

  //Run the CUDA code
  test_kern_calc_visi_gaussian(num_components, num_baselines,
          num_freqs, num_visis, num_times, beamtype,
          args_ft->component_freqs,
          args_ft->flux_I, args_ft->flux_Q, args_ft->flux_U, args_ft->flux_V,
          args_ft->SIs, args_ft->us, args_ft->vs, args_ft->ws,
          args_ft->sum_visi_XX_real, args_ft->sum_visi_XX_imag,
          args_ft->sum_visi_XY_real, args_ft->sum_visi_XY_imag,
          args_ft->sum_visi_YX_real, args_ft->sum_visi_YX_imag,
          args_ft->sum_visi_YY_real, args_ft->sum_visi_YY_imag,
          args_ft->allsteps_wavelengths,
          args_ft->ls, args_ft->ms, args_ft->ns,
          args_ft->pas, args_ft->majors, args_ft->minors,
          args_ft->primay_beam_J00, args_ft->primay_beam_J01,
          args_ft->primay_beam_J10, args_ft->primay_beam_J11);

  test_visi_outputs(num_visis, num_components, num_baselines, num_freqs,
                    FRAC_TOL, beamtype, args_ft, GAUSSIAN);

  free_args_for_testing( args_ft, GAUSSIAN );
}

/*
This test checks varying the measurement equation with beamtype=FEE_BEAM
*/
void test_kern_calc_visi_gauss_VarylmnVaryFluxFEEBeam(void) {
  test_kern_calc_visi_gauss_VarylmnVaryFlux(FEE_BEAM);
}

/*
This test checks varying the measurement equation with beamtype=ANALY_DIPOLE
*/
void test_kern_calc_visi_gauss_VarylmnVaryFluxAnalyBeam(void) {
  test_kern_calc_visi_gauss_VarylmnVaryFlux(ANALY_DIPOLE);
}

/*
This test checks varying the measurement equation with beamtype=GAUSS_BEAM
*/
void test_kern_calc_visi_gauss_VarylmnVaryFluxGaussBeam(void) {
  test_kern_calc_visi_gauss_VarylmnVaryFlux(GAUSS_BEAM);
}

/*
This test checks varying the measurement equation with beamtype=NO_BEAM
*/
void test_kern_calc_visi_gauss_VarylmnVaryFluxNoBeam(void) {
  test_kern_calc_visi_gauss_VarylmnVaryFlux(NO_BEAM);
}


void test_kern_calc_visi_gauss_VarylmnVaryBeam(int beamtype) {

  int num_baselines = 10.0;
  int num_times = 5.0;
  int num_freqs = 3.0;

  int num_visis = num_baselines*num_times*num_freqs;

  int num_components = 25;

  //Container for many arrays to feed the GPU
  args_for_testing_t *args_ft = malloc(sizeof(args_for_testing_t));
  //Allocate memory
  int num_coeffs = 0; //Only applies for SHAPELET
  malloc_args_for_testing(args_ft, num_baselines, num_times,
                          num_freqs, num_components, num_coeffs, GAUSSIAN);

  //Setup l,m args that span a decent chunk of sky
  create_lmn(args_ft);

  int num_beam_values = num_freqs*num_times*num_components;

  //Stick the gains to one everywhere
  for (size_t beam = 0; beam < num_beam_values; beam++) {
    args_ft->primay_beam_J00[beam] = beam + I*0.0;
    args_ft->primay_beam_J01[beam] = beam + I*0.0;
    args_ft->primay_beam_J10[beam] = beam + I*0.0;
    args_ft->primay_beam_J11[beam] = beam + I*0.0;
  }

  //Just stick Stokes I to 1.0, SI to zero, and reference freqs to 150MHz
  for (size_t comp = 0; comp < num_components; comp++) {
    args_ft->flux_I[comp] = 1.0;
    args_ft->SIs[comp] = 0.0;
    args_ft->component_freqs[comp] = 150e+6;
    args_ft->pas[comp] = 0.0;
    //Set major,minor to 3 arcmins
    args_ft->majors[comp] = 3.0*(DD2R / 60.0);
    args_ft->minors[comp] = 3.0*(DD2R / 60.0);
  }

  //Make up some u,v,w values and scale by wavelength in correct order
  int count = 0;
  user_precision_t freq_base = 150e+6;
  user_precision_t freq_inc = 25e+6;
  user_precision_t wavelength, frequency;

  for ( int time_step = 0; time_step < num_times; time_step++ ) {
    for (int freq_step = 0; freq_step < num_freqs; freq_step++) {
      frequency = freq_base + freq_step*freq_inc;
      wavelength = VELC / frequency;

      for (int baseline = 0; baseline < num_baselines; baseline++) {
        args_ft->us[count] = ((baseline + 1)*100) / wavelength;
        args_ft->vs[count] = ((baseline + 1)*100) / wavelength;
        //ws are usually smaller than u,v
        args_ft->ws[count] = ((baseline + 1)*10) / wavelength;

        // args_ft->us[count] = 0.0;
        // args_ft->vs[count] = 0.0;
        // args_ft->ws[count] = 0.0;

        args_ft->allsteps_wavelengths[count] = wavelength;

        count ++;
      }
    }
  }

  //Run the CUDA code
  test_kern_calc_visi_gaussian(num_components, num_baselines,
          num_freqs, num_visis, num_times, beamtype,
          args_ft->component_freqs,
          args_ft->flux_I, args_ft->flux_Q, args_ft->flux_U, args_ft->flux_V,
          args_ft->SIs, args_ft->us, args_ft->vs, args_ft->ws,
          args_ft->sum_visi_XX_real, args_ft->sum_visi_XX_imag,
          args_ft->sum_visi_XY_real, args_ft->sum_visi_XY_imag,
          args_ft->sum_visi_YX_real, args_ft->sum_visi_YX_imag,
          args_ft->sum_visi_YY_real, args_ft->sum_visi_YY_imag,
          args_ft->allsteps_wavelengths,
          args_ft->ls, args_ft->ms, args_ft->ns,
          args_ft->pas, args_ft->majors, args_ft->minors,
          args_ft->primay_beam_J00, args_ft->primay_beam_J01,
          args_ft->primay_beam_J10, args_ft->primay_beam_J11);

  test_visi_outputs(num_visis, num_components, num_baselines, num_freqs,
                    FRAC_TOL, beamtype, args_ft, GAUSSIAN);

  free_args_for_testing( args_ft, GAUSSIAN );
}

/*
This test checks varying the measurement equation with beamtype=FEE_BEAM
*/
void test_kern_calc_visi_gauss_VarylmnVaryBeamFEEBeam(void) {
  test_kern_calc_visi_gauss_VarylmnVaryBeam(FEE_BEAM);
}

/*
This test checks varying the measurement equation with beamtype=ANALY_DIPOLE
*/
void test_kern_calc_visi_gauss_VarylmnVaryBeamAnalyBeam(void) {
  test_kern_calc_visi_gauss_VarylmnVaryBeam(ANALY_DIPOLE);
}

/*
This test checks varying the measurement equation with beamtype=GAUSS_BEAM
*/
void test_kern_calc_visi_gauss_VarylmnVaryBeamGaussBeam(void) {
  test_kern_calc_visi_gauss_VarylmnVaryBeam(GAUSS_BEAM);
}

/*
This test checks varying the measurement equation with beamtype=NO_BEAM
*/
void test_kern_calc_visi_gauss_VarylmnVaryBeamNoBeam(void) {
  test_kern_calc_visi_gauss_VarylmnVaryBeam(NO_BEAM);
}




void test_kern_calc_visi_gauss_VarylmnVaryPAMajMin(int beamtype) {

  int num_baselines = 10.0;
  int num_times = 5.0;
  int num_freqs = 3.0;

  int num_visis = num_baselines*num_times*num_freqs;

  int num_components = 25;

  //Container for many arrays to feed the GPU
  args_for_testing_t *args_ft = malloc(sizeof(args_for_testing_t));
  //Allocate memory
  int num_coeffs = 0; //Only applies for SHAPELET
  malloc_args_for_testing(args_ft, num_baselines, num_times,
                          num_freqs, num_components, num_coeffs, GAUSSIAN);

  //Setup l,m args that span a decent chunk of sky
  create_lmn(args_ft);

  int num_beam_values = num_freqs*num_times*num_components;

  //Stick the gains to one everywhere
  for (size_t visi = 0; visi < num_beam_values; visi++) {
    args_ft->primay_beam_J00[visi] = 1.0 + I*0.0;
    args_ft->primay_beam_J01[visi] = 1.0 + I*0.0;
    args_ft->primay_beam_J10[visi] = 1.0 + I*0.0;
    args_ft->primay_beam_J11[visi] = 1.0 + I*0.0;
  }

  //Just stick Stokes I to 1.0, SI to zero, and reference freqs to 150MHz
  for (size_t comp = 0; comp < num_components; comp++) {
    args_ft->flux_I[comp] = 1.0;
    args_ft->SIs[comp] = 0.0;
    args_ft->component_freqs[comp] = 150e+6;
    //Vary the PA, minor and major axis
    args_ft->pas[comp] = (comp + 1)*DD2R;
    args_ft->majors[comp] = (comp + 1)*(DD2R / 60.0);
    args_ft->minors[comp] = (comp + 2)*(DD2R / 60.0);
  }

  //Make up some u,v,w values and scale by wavelength in correct order
  int count = 0;
  user_precision_t freq_base = 150e+6;
  user_precision_t freq_inc = 25e+6;
  user_precision_t wavelength, frequency;
  for ( int time_step = 0; time_step < num_times; time_step++ ) {
    for (int freq_step = 0; freq_step < num_freqs; freq_step++) {
      frequency = freq_base + freq_step*freq_inc;
      wavelength = VELC / frequency;
      for (int baseline = 0; baseline < num_baselines; baseline++) {
        args_ft->us[count] = ((baseline + 1)*100) / wavelength;
        args_ft->vs[count] = ((baseline + 1)*100) / wavelength;
        //ws are usually smaller than u,v
        args_ft->ws[count] = ((baseline + 1)*10) / wavelength;

        args_ft->allsteps_wavelengths[count] = wavelength;

        count ++;
      }
    }
  }

  //Run the CUDA code
  test_kern_calc_visi_gaussian(num_components, num_baselines,
          num_freqs, num_visis, num_times, beamtype,
          args_ft->component_freqs,
          args_ft->flux_I, args_ft->flux_Q, args_ft->flux_U, args_ft->flux_V,
          args_ft->SIs, args_ft->us, args_ft->vs, args_ft->ws,
          args_ft->sum_visi_XX_real, args_ft->sum_visi_XX_imag,
          args_ft->sum_visi_XY_real, args_ft->sum_visi_XY_imag,
          args_ft->sum_visi_YX_real, args_ft->sum_visi_YX_imag,
          args_ft->sum_visi_YY_real, args_ft->sum_visi_YY_imag,
          args_ft->allsteps_wavelengths,
          args_ft->ls, args_ft->ms, args_ft->ns,
          args_ft->pas, args_ft->majors, args_ft->minors,
          args_ft->primay_beam_J00, args_ft->primay_beam_J01,
          args_ft->primay_beam_J10, args_ft->primay_beam_J11);

  test_visi_outputs(num_visis, num_components, num_baselines, num_freqs,
                    FRAC_TOL, beamtype, args_ft, GAUSSIAN);

  free_args_for_testing( args_ft, GAUSSIAN );
}

/*
This test checks varying the measurement equation with beamtype=FEE_BEAM
*/
void test_kern_calc_visi_gauss_VarylmnVaryPAMajMinFEEBeam(void) {
  test_kern_calc_visi_gauss_VarylmnVaryPAMajMin(FEE_BEAM);
}

/*
This test checks varying the measurement equation with beamtype=ANALY_DIPOLE
*/
void test_kern_calc_visi_gauss_VarylmnVaryPAMajMinAnalyBeam(void) {
  test_kern_calc_visi_gauss_VarylmnVaryPAMajMin(ANALY_DIPOLE);
}

/*
This test checks varying the measurement equation with beamtype=GAUSS_BEAM
*/
void test_kern_calc_visi_gauss_VarylmnVaryPAMajMinGaussBeam(void) {
  test_kern_calc_visi_gauss_VarylmnVaryPAMajMin(GAUSS_BEAM);
}

/*
This test checks varying the measurement equation with beamtype=NO_BEAM
*/
void test_kern_calc_visi_gauss_VarylmnVaryPAMajMinNoBeam(void) {
  test_kern_calc_visi_gauss_VarylmnVaryPAMajMin(NO_BEAM);
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

    RUN_TEST(test_kern_calc_visi_gauss_VarylmnVaryFluxFEEBeam);
    RUN_TEST(test_kern_calc_visi_gauss_VarylmnVaryFluxGaussBeam);
    RUN_TEST(test_kern_calc_visi_gauss_VarylmnVaryFluxAnalyBeam);
    RUN_TEST(test_kern_calc_visi_gauss_VarylmnVaryFluxNoBeam);

    RUN_TEST(test_kern_calc_visi_gauss_VarylmnVaryBeamFEEBeam);
    RUN_TEST(test_kern_calc_visi_gauss_VarylmnVaryBeamGaussBeam);
    RUN_TEST(test_kern_calc_visi_gauss_VarylmnVaryBeamAnalyBeam);
    RUN_TEST(test_kern_calc_visi_gauss_VarylmnVaryBeamNoBeam);

    RUN_TEST(test_kern_calc_visi_gauss_VarylmnVaryPAMajMinFEEBeam);
    RUN_TEST(test_kern_calc_visi_gauss_VarylmnVaryPAMajMinGaussBeam);
    RUN_TEST(test_kern_calc_visi_gauss_VarylmnVaryPAMajMinAnalyBeam);
    RUN_TEST(test_kern_calc_visi_gauss_VarylmnVaryPAMajMinNoBeam);

    return UNITY_END();
}
