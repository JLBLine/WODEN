#include <unity.h>
#include <stdlib.h>
#include <complex.h>

#include "constants.h"
#include "woden_struct_defs.h"
#include "test_kern_calc_visi_common.h"
#include "shapelet_basis.h"

void setUp (void) {} /* Is run before eVary test, put unit init calls here. */
void tearDown (void) {} /* Is run after eVary test, put unit clean-up calls here. */

//External CUDA code we're linking in
extern void test_kern_calc_visi_shapelet(int num_components, int num_baselines,
          int num_freqs, int num_visis, int num_times, int num_coeffs, int beamtype,
          double *component_freqs,
          user_precision_t *flux_I, user_precision_t *flux_Q,
          user_precision_t *flux_U, user_precision_t *flux_V,
          user_precision_t *SIs, user_precision_t *us, user_precision_t *vs, user_precision_t *ws,
          user_precision_t *sum_visi_XX_real, user_precision_t *sum_visi_XX_imag,
          user_precision_t *sum_visi_XY_real, user_precision_t *sum_visi_XY_imag,
          user_precision_t *sum_visi_YX_real, user_precision_t *sum_visi_YX_imag,
          user_precision_t *sum_visi_YY_real, user_precision_t *sum_visi_YY_imag,
          user_precision_t *allsteps_wavelengths,
          double *ls, double *ms, double *ns,
          user_precision_t *shape_pas, user_precision_t *shape_majors,
          user_precision_t *shape_minors,
          user_precision_complex_t *primay_beam_J00, user_precision_complex_t *primay_beam_J01,
          user_precision_complex_t *primay_beam_J10, user_precision_complex_t *primay_beam_J11,
          user_precision_t *u_shapes, user_precision_t *v_shapes, user_precision_t *w_shapes,
          user_precision_t *shape_n1s, user_precision_t *shape_n2s, user_precision_t *shape_coeffs,
          user_precision_t *shape_param_indexes, user_precision_t *sbf);

//Change required accuracy of outputs for different precisions
//This is a fractional tolerance, not an absolute
#ifdef DOUBLE_PRECISION
  double FRAC_TOL = 1e-12;
#else
  double FRAC_TOL = 1e-2;
#endif

/*
Test the __global__ code that calculates visibilities for shape sources
Vary the l,m,n coords but keep all other variables constant
*/
void test_kern_calc_visi_shape_Varylmn(int beamtype) {

  int num_baselines = 10.0;
  int num_times = 5.0;
  int num_freqs = 3.0;

  int num_visis = num_baselines*num_times*num_freqs;

  int num_components = 25;
  int num_coeffs = num_components;

  //Container for many arrays to feed the GPU
  args_for_testing_t *args_ft = malloc(sizeof(args_for_testing_t));
  //Allocate memory
  malloc_args_for_testing(args_ft, num_baselines, num_times,
                          num_freqs, num_components, num_coeffs, SHAPELET);

  //Setup l,m args that span a decent chunk of sky
  create_lmn(args_ft);

  int num_beam_values = num_freqs*num_times*num_components;

  //Stick the gains to one everywhere
  for (int visi = 0; visi < num_beam_values; visi++) {
    args_ft->primay_beam_J00[visi] = 1.0 + I*0.0;
    args_ft->primay_beam_J01[visi] = 1.0 + I*0.0;
    args_ft->primay_beam_J10[visi] = 1.0 + I*0.0;
    args_ft->primay_beam_J11[visi] = 1.0 + I*0.0;
  }

  //Just stick Stokes I to 1.0, SI to zero, and reference freqs to 150MHz
  for (int comp = 0; comp < num_components; comp++) {
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
  double freq_base = 150e+6;
  double freq_inc = 25e+6;
  user_precision_t wavelength;
  double frequency;
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

  //Set the shapelet u,v,w same as the measurement equation one (this is not
  //true in reality but works fine for testing)

  for (int comp = 0; comp < num_components; comp++) {
    for (int visi = 0; visi < num_visis; visi++) {
      args_ft->u_shapes[comp*num_visis + visi] = args_ft->us[visi];
      args_ft->v_shapes[comp*num_visis + visi] = args_ft->vs[visi];
      args_ft->w_shapes[comp*num_visis + visi] = args_ft->ws[visi];
    }
  }

  //This means we have a single basis function for every component
  for (int coeff = 0; coeff < num_coeffs; coeff++) {
    args_ft->shape_n1s[coeff] = 0.0;
    args_ft->shape_n2s[coeff] = 0.0;
    args_ft->shape_coeffs[coeff] = 1.0;
    args_ft->shape_param_indexes[coeff] = coeff;
  }

  //Run the CUDA code
  test_kern_calc_visi_shapelet(num_components, num_baselines,
          num_freqs, num_visis, num_times, num_coeffs, beamtype,
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
          args_ft->primay_beam_J10, args_ft->primay_beam_J11,
          args_ft->u_shapes, args_ft->v_shapes, args_ft->w_shapes,
          args_ft->shape_n1s, args_ft->shape_n2s, args_ft->shape_coeffs,
          args_ft->shape_param_indexes, args_ft->sbf);

  //SHAPELET is NOT a typo! With a single n1, n2, coeff = 0, 0, 1
  //shapelet basis function, outputs should be identical to a Gaussian
  test_visi_outputs(num_visis, num_components, num_baselines, num_freqs,
                    FRAC_TOL, beamtype, args_ft, SHAPELET);
  //
  free_args_for_testing( args_ft, SHAPELET );
}

/*
This test checks varying the measurement equation with beamtype=FEE_BEAM
*/
void test_kern_calc_visi_shape_VarylmnFEEBeam(void) {
  test_kern_calc_visi_shape_Varylmn(FEE_BEAM);
}

/*
This test checks varying the measurement equation with beamtype=ANALY_DIPOLE
*/
void test_kern_calc_visi_shape_VarylmnAnalyBeam(void) {
  test_kern_calc_visi_shape_Varylmn(ANALY_DIPOLE);
}

/*
This test checks varying the measurement equation with beamtype=GAUSS_BEAM
*/
void test_kern_calc_visi_shape_VarylmnGaussBeam(void) {
  test_kern_calc_visi_shape_Varylmn(GAUSS_BEAM);
}

/*
This test checks varying the measurement equation with beamtype=NO_BEAM
*/
void test_kern_calc_visi_shape_VarylmnNoBeam(void) {
  test_kern_calc_visi_shape_Varylmn(NO_BEAM);
}

/*
This test checks varying the measurement equation with beamtype=FEE_BEAM_INTERP
*/
void test_kern_calc_visi_shape_VarylmnFEEInterpBeam(void) {
  test_kern_calc_visi_shape_Varylmn(FEE_BEAM_INTERP);
}

/*
This test checks varying the measurement equation with beamtype=MWA_ANALY
*/
void test_kern_calc_visi_shape_VarylmnMWAAnaly(void) {
  test_kern_calc_visi_shape_Varylmn(MWA_ANALY);
}


/*
Test the __device__ code that updates the summed visibilities by grabbing the
correct beam gain and mesurement equation, multiplying and summing onto the visi
Here we keep the component visibilities and fluxes constant and vary the beam gains
Test works for all primary beam types
*/
void test_kern_calc_visi_shape_VarylmnVaryFlux(int beamtype) {

  int num_baselines = 10.0;
  int num_times = 5.0;
  int num_freqs = 3.0;

  int num_visis = num_baselines*num_times*num_freqs;

  int num_components = 25;
  int num_coeffs = num_components;

  //Container for many arrays to feed the GPU
  args_for_testing_t *args_ft = malloc(sizeof(args_for_testing_t));
  //Allocate memory
  malloc_args_for_testing(args_ft, num_baselines, num_times,
                          num_freqs, num_components, num_coeffs, SHAPELET);

  //Setup l,m args that span a decent chunk of sky
  create_lmn(args_ft);

  int num_beam_values = num_freqs*num_times*num_components;

  //Stick the gains to one everywhere
  for (int visi = 0; visi < num_beam_values; visi++) {
    args_ft->primay_beam_J00[visi] = 1.0 + I*0.0;
    args_ft->primay_beam_J01[visi] = 1.0 + I*0.0;
    args_ft->primay_beam_J10[visi] = 1.0 + I*0.0;
    args_ft->primay_beam_J11[visi] = 1.0 + I*0.0;
  }

  //Just stick Stokes I to 1.0, SI to zero, and reference freqs to 150MHz
  for (int comp = 0; comp < num_components; comp++) {
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
  double freq_base = 150e+6;
  double freq_inc = 25e+6;
  user_precision_t wavelength;
  double frequency;

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

  //Set the shapelet u,v,w same as the measurement equation one (this is not
  //true in reality but works fine for testing)
  for (int comp = 0; comp < num_components; comp++) {
    for (int visi = 0; visi < num_visis; visi++) {
      args_ft->u_shapes[comp*num_visis + visi] = args_ft->us[visi];
      args_ft->v_shapes[comp*num_visis + visi] = args_ft->vs[visi];
      args_ft->w_shapes[comp*num_visis + visi] = args_ft->ws[visi];
    }
  }

  //This means we have a single basis function for every component
  for (int coeff = 0; coeff < num_coeffs; coeff++) {
    args_ft->shape_n1s[coeff] = 0.0;
    args_ft->shape_n2s[coeff] = 0.0;
    args_ft->shape_coeffs[coeff] = 1.0;
    args_ft->shape_param_indexes[coeff] = coeff;
  }

  //Run the CUDA code
  test_kern_calc_visi_shapelet(num_components, num_baselines,
          num_freqs, num_visis, num_times, num_coeffs, beamtype,
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
          args_ft->primay_beam_J10, args_ft->primay_beam_J11,
          args_ft->u_shapes, args_ft->v_shapes, args_ft->w_shapes,
          args_ft->shape_n1s, args_ft->shape_n2s, args_ft->shape_coeffs,
          args_ft->shape_param_indexes, args_ft->sbf);

  test_visi_outputs(num_visis, num_components, num_baselines, num_freqs,
                    FRAC_TOL, beamtype, args_ft, SHAPELET);

  free_args_for_testing( args_ft, SHAPELET );
}

/*
This test checks varying the measurement equation with beamtype=FEE_BEAM
*/
void test_kern_calc_visi_shape_VarylmnVaryFluxFEEBeam(void) {
  test_kern_calc_visi_shape_VarylmnVaryFlux(FEE_BEAM);
}

/*
This test checks varying the measurement equation with beamtype=ANALY_DIPOLE
*/
void test_kern_calc_visi_shape_VarylmnVaryFluxAnalyBeam(void) {
  test_kern_calc_visi_shape_VarylmnVaryFlux(ANALY_DIPOLE);
}

/*
This test checks varying the measurement equation with beamtype=GAUSS_BEAM
*/
void test_kern_calc_visi_shape_VarylmnVaryFluxGaussBeam(void) {
  test_kern_calc_visi_shape_VarylmnVaryFlux(GAUSS_BEAM);
}

/*
This test checks varying the measurement equation with beamtype=NO_BEAM
*/
void test_kern_calc_visi_shape_VarylmnVaryFluxNoBeam(void) {
  test_kern_calc_visi_shape_VarylmnVaryFlux(NO_BEAM);
}

/*
This test checks varying the measurement equation with beamtype=FEE_BEAM_INTERP
*/
void test_kern_calc_visi_shape_VarylmnVaryFluxFEEInterpBeam(void) {
  test_kern_calc_visi_shape_VarylmnVaryFlux(FEE_BEAM_INTERP);
}

/*
This test checks varying the measurement equation with beamtype=MWA_ANALY
*/
void test_kern_calc_visi_shape_VarylmnVaryFluxMWAAnaly(void) {
  test_kern_calc_visi_shape_VarylmnVaryFlux(MWA_ANALY);
}


void test_kern_calc_visi_shape_VarylmnVaryBeam(int beamtype) {

  int num_baselines = 10.0;
  int num_times = 5.0;
  int num_freqs = 3.0;

  int num_visis = num_baselines*num_times*num_freqs;

  int num_components = 25;
  int num_coeffs = num_components;

  //Container for many arrays to feed the GPU
  args_for_testing_t *args_ft = malloc(sizeof(args_for_testing_t));
  //Allocate memory
  malloc_args_for_testing(args_ft, num_baselines, num_times,
                          num_freqs, num_components, num_coeffs, SHAPELET);

  //Setup l,m args that span a decent chunk of sky
  create_lmn(args_ft);

  int num_beam_values = num_freqs*num_times*num_components;

  //Stick the gains to one everywhere
  for (int beam = 0; beam < num_beam_values; beam++) {
    args_ft->primay_beam_J00[beam] = beam + I*0.0;
    args_ft->primay_beam_J01[beam] = beam + I*0.0;
    args_ft->primay_beam_J10[beam] = beam + I*0.0;
    args_ft->primay_beam_J11[beam] = beam + I*0.0;
  }

  //Just stick Stokes I to 1.0, SI to zero, and reference freqs to 150MHz
  for (int comp = 0; comp < num_components; comp++) {
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
  double freq_base = 150e+6;
  double freq_inc = 25e+6;
  user_precision_t wavelength;
double frequency;

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

  //Set the shapelet u,v,w same as the measurement equation one (this is not
  //true in reality but works fine for testing)
  for (int comp = 0; comp < num_components; comp++) {
    for (int visi = 0; visi < num_visis; visi++) {
      args_ft->u_shapes[comp*num_visis + visi] = args_ft->us[visi];
      args_ft->v_shapes[comp*num_visis + visi] = args_ft->vs[visi];
      args_ft->w_shapes[comp*num_visis + visi] = args_ft->ws[visi];
    }
  }

  //This means we have a single basis function for every component
  for (int coeff = 0; coeff < num_coeffs; coeff++) {
    args_ft->shape_n1s[coeff] = 0.0;
    args_ft->shape_n2s[coeff] = 0.0;
    args_ft->shape_coeffs[coeff] = 1.0;
    args_ft->shape_param_indexes[coeff] = coeff;
  }

  //Run the CUDA code
  test_kern_calc_visi_shapelet(num_components, num_baselines,
          num_freqs, num_visis, num_times, num_coeffs, beamtype,
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
          args_ft->primay_beam_J10, args_ft->primay_beam_J11,
          args_ft->u_shapes, args_ft->v_shapes, args_ft->w_shapes,
          args_ft->shape_n1s, args_ft->shape_n2s, args_ft->shape_coeffs,
          args_ft->shape_param_indexes, args_ft->sbf);

  test_visi_outputs(num_visis, num_components, num_baselines, num_freqs,
                    FRAC_TOL, beamtype, args_ft, SHAPELET);

  free_args_for_testing( args_ft, SHAPELET );
}

/*
This test checks varying the measurement equation with beamtype=FEE_BEAM
*/
void test_kern_calc_visi_shape_VarylmnVaryBeamFEEBeam(void) {
  test_kern_calc_visi_shape_VarylmnVaryBeam(FEE_BEAM);
}

/*
This test checks varying the measurement equation with beamtype=ANALY_DIPOLE
*/
void test_kern_calc_visi_shape_VarylmnVaryBeamAnalyBeam(void) {
  test_kern_calc_visi_shape_VarylmnVaryBeam(ANALY_DIPOLE);
}

/*
This test checks varying the measurement equation with beamtype=GAUSS_BEAM
*/
void test_kern_calc_visi_shape_VarylmnVaryBeamGaussBeam(void) {
  test_kern_calc_visi_shape_VarylmnVaryBeam(GAUSS_BEAM);
}

/*
This test checks varying the measurement equation with beamtype=NO_BEAM
*/
void test_kern_calc_visi_shape_VarylmnVaryBeamNoBeam(void) {
  test_kern_calc_visi_shape_VarylmnVaryBeam(NO_BEAM);
}

/*
This test checks varying the measurement equation with beamtype=FEE_BEAM_INTERP
*/
void test_kern_calc_visi_shape_VarylmnVaryBeamFEEInterpBeam(void) {
  test_kern_calc_visi_shape_VarylmnVaryBeam(FEE_BEAM_INTERP);
}

/*
This test checks varying the measurement equation with beamtype=MWA_ANALY
*/
void test_kern_calc_visi_shape_VarylmnVaryBeamMWAAnaly(void) {
  test_kern_calc_visi_shape_VarylmnVaryBeam(MWA_ANALY);
}


void test_kern_calc_visi_shape_VarylmnVaryPAMajMin(int beamtype) {

  int num_baselines = 10.0;
  int num_times = 5.0;
  int num_freqs = 3.0;

  int num_visis = num_baselines*num_times*num_freqs;

  int num_components = 25;
  int num_coeffs = num_components;

  //Container for many arrays to feed the GPU
  args_for_testing_t *args_ft = malloc(sizeof(args_for_testing_t));
  //Allocate memory
  malloc_args_for_testing(args_ft, num_baselines, num_times,
                          num_freqs, num_components, num_coeffs, SHAPELET);

  //Setup l,m args that span a decent chunk of sky
  create_lmn(args_ft);

  int num_beam_values = num_freqs*num_times*num_components;

  //Stick the gains to one everywhere
  for (int visi = 0; visi < num_beam_values; visi++) {
    args_ft->primay_beam_J00[visi] = 1.0 + I*0.0;
    args_ft->primay_beam_J01[visi] = 1.0 + I*0.0;
    args_ft->primay_beam_J10[visi] = 1.0 + I*0.0;
    args_ft->primay_beam_J11[visi] = 1.0 + I*0.0;
  }

  //Just stick Stokes I to 1.0, SI to zero, and reference freqs to 150MHz
  for (int comp = 0; comp < num_components; comp++) {
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
  double freq_base = 150e+6;
  double freq_inc = 25e+6;
  user_precision_t wavelength;
double frequency;
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

  //Set the shapelet u,v,w same as the measurement equation one (this is not
  //true in reality but works fine for testing)

  for (int comp = 0; comp < num_components; comp++) {
    for (int visi = 0; visi < num_visis; visi++) {
      args_ft->u_shapes[comp*num_visis + visi] = args_ft->us[visi];
      args_ft->v_shapes[comp*num_visis + visi] = args_ft->vs[visi];
      args_ft->w_shapes[comp*num_visis + visi] = args_ft->ws[visi];
    }
  }

  //This means we have a single basis function for every component
  for (int coeff = 0; coeff < num_coeffs; coeff++) {
    args_ft->shape_n1s[coeff] = 0.0;
    args_ft->shape_n2s[coeff] = 0.0;
    args_ft->shape_coeffs[coeff] = 1.0;
    args_ft->shape_param_indexes[coeff] = coeff;
  }

  //Run the CUDA code
  test_kern_calc_visi_shapelet(num_components, num_baselines,
          num_freqs, num_visis, num_times, num_coeffs, beamtype,
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
          args_ft->primay_beam_J10, args_ft->primay_beam_J11,
          args_ft->u_shapes, args_ft->v_shapes, args_ft->w_shapes,
          args_ft->shape_n1s, args_ft->shape_n2s, args_ft->shape_coeffs,
          args_ft->shape_param_indexes, args_ft->sbf);

  test_visi_outputs(num_visis, num_components, num_baselines, num_freqs,
                    FRAC_TOL, beamtype, args_ft, SHAPELET);

  free_args_for_testing( args_ft, SHAPELET );
}

/*
This test checks varying the measurement equation with beamtype=FEE_BEAM
*/
void test_kern_calc_visi_shape_VarylmnVaryPAMajMinFEEBeam(void) {
  test_kern_calc_visi_shape_VarylmnVaryPAMajMin(FEE_BEAM);
}

/*
This test checks varying the measurement equation with beamtype=ANALY_DIPOLE
*/
void test_kern_calc_visi_shape_VarylmnVaryPAMajMinAnalyBeam(void) {
  test_kern_calc_visi_shape_VarylmnVaryPAMajMin(ANALY_DIPOLE);
}

/*
This test checks varying the measurement equation with beamtype=GAUSS_BEAM
*/
void test_kern_calc_visi_shape_VarylmnVaryPAMajMinGaussBeam(void) {
  test_kern_calc_visi_shape_VarylmnVaryPAMajMin(GAUSS_BEAM);
}

/*
This test checks varying the measurement equation with beamtype=NO_BEAM
*/
void test_kern_calc_visi_shape_VarylmnVaryPAMajMinNoBeam(void) {
  test_kern_calc_visi_shape_VarylmnVaryPAMajMin(NO_BEAM);
}

/*
This test checks varying the measurement equation with beamtype=FEE_BEAM_INTERP
*/
void test_kern_calc_visi_shape_VarylmnVaryPAMajMinFEEInterpBeam(void) {
  test_kern_calc_visi_shape_VarylmnVaryPAMajMin(FEE_BEAM_INTERP);
}

/*
This test checks varying the measurement equation with beamtype=MWA_ANALY
*/
void test_kern_calc_visi_shape_VarylmnVaryPAMajMinMWAAnaly(void) {
  test_kern_calc_visi_shape_VarylmnVaryPAMajMin(MWA_ANALY);
}

void test_kern_calc_visi_shape_VarylmnMultipleCoeff(int beamtype) {

  int num_baselines = 10.0;
  int num_times = 5.0;
  int num_freqs = 3.0;

  int num_visis = num_baselines*num_times*num_freqs;

  int num_components = 25;
  int num_coeffs_per_component = 3;

  int num_coeffs = num_coeffs_per_component*num_components;

  //Container for many arrays to feed the GPU
  args_for_testing_t *args_ft = malloc(sizeof(args_for_testing_t));
  //Allocate memory
  malloc_args_for_testing(args_ft, num_baselines, num_times,
                          num_freqs, num_components, num_coeffs, SHAPELET);

  //Setup l,m args that span a decent chunk of sky
  create_lmn(args_ft);

  int num_beam_values = num_freqs*num_times*num_components;

  //Stick the gains to one everywhere
  for (int visi = 0; visi < num_beam_values; visi++) {
    args_ft->primay_beam_J00[visi] = 1.0 + I*0.0;
    args_ft->primay_beam_J01[visi] = 1.0 + I*0.0;
    args_ft->primay_beam_J10[visi] = 1.0 + I*0.0;
    args_ft->primay_beam_J11[visi] = 1.0 + I*0.0;
  }

  //Just stick Stokes I to 1.0, SI to zero, and reference freqs to 150MHz
  for (int comp = 0; comp < num_components; comp++) {
    args_ft->flux_I[comp] = 1.0;
    args_ft->SIs[comp] = 0.0;
    args_ft->component_freqs[comp] = 150e+6;
    args_ft->pas[comp] = 0.0;
    args_ft->majors[comp] = 3.0*(DD2R / 60.0);
    args_ft->minors[comp] = 3.0*(DD2R / 60.0);
  }

  //Make up some u,v,w values and scale by wavelength in correct order
  int count = 0;
  double freq_base = 150e+6;
  double freq_inc = 25e+6;
  user_precision_t wavelength;
  double frequency;
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

  //Set the shapelet u,v,w same as the measurement equation one (this is not
  //true in reality but works fine for testing)

  for (int comp = 0; comp < num_components; comp++) {
    for (int visi = 0; visi < num_visis; visi++) {
      args_ft->u_shapes[comp*num_visis + visi] = args_ft->us[visi];
      args_ft->v_shapes[comp*num_visis + visi] = args_ft->vs[visi];
      args_ft->w_shapes[comp*num_visis + visi] = args_ft->ws[visi];
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

      args_ft->shape_n1s[count] = count;
      args_ft->shape_n2s[count] = count + 1;
      args_ft->shape_coeffs[count] = sign*1e-3*(coeff + 1);
      args_ft->shape_param_indexes[count] = comp;

      count ++;
    }
  }

  //Run the CUDA code
  test_kern_calc_visi_shapelet(num_components, num_baselines,
          num_freqs, num_visis, num_times, num_coeffs, beamtype,
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
          args_ft->primay_beam_J10, args_ft->primay_beam_J11,
          args_ft->u_shapes, args_ft->v_shapes, args_ft->w_shapes,
          args_ft->shape_n1s, args_ft->shape_n2s, args_ft->shape_coeffs,
          args_ft->shape_param_indexes, args_ft->sbf);

  test_visi_outputs(num_visis, num_components, num_baselines, num_freqs,
                    FRAC_TOL, beamtype, args_ft, SHAPELET);

  free_args_for_testing( args_ft, SHAPELET );
}

/*
This test checks varying the measurement equation with beamtype=FEE_BEAM
*/
void test_kern_calc_visi_shape_VarylmnMultipleCoeffFEEBeam(void) {
  test_kern_calc_visi_shape_VarylmnMultipleCoeff(FEE_BEAM);
}

/*
This test checks varying the measurement equation with beamtype=ANALY_DIPOLE
*/
void test_kern_calc_visi_shape_VarylmnMultipleCoeffAnalyBeam(void) {
  test_kern_calc_visi_shape_VarylmnMultipleCoeff(ANALY_DIPOLE);
}

/*
This test checks varying the measurement equation with beamtype=GAUSS_BEAM
*/
void test_kern_calc_visi_shape_VarylmnMultipleCoeffGaussBeam(void) {
  test_kern_calc_visi_shape_VarylmnMultipleCoeff(GAUSS_BEAM);
}

/*
This test checks varying the measurement equation with beamtype=NO_BEAM
*/
void test_kern_calc_visi_shape_VarylmnMultipleCoeffNoBeam(void) {
  test_kern_calc_visi_shape_VarylmnMultipleCoeff(NO_BEAM);
}

/*
This test checks varying the measurement equation with beamtype=FEE_BEAM_INTERP
*/
void test_kern_calc_visi_shape_VarylmnMultipleCoeffFEEInterpBeam(void) {
  test_kern_calc_visi_shape_VarylmnMultipleCoeff(FEE_BEAM_INTERP);
}

/*
This test checks varying the measurement equation with beamtype=MWA_ANALY
*/
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

    RUN_TEST(test_kern_calc_visi_shape_VarylmnVaryFluxFEEBeam);
    RUN_TEST(test_kern_calc_visi_shape_VarylmnVaryFluxGaussBeam);
    RUN_TEST(test_kern_calc_visi_shape_VarylmnVaryFluxAnalyBeam);
    RUN_TEST(test_kern_calc_visi_shape_VarylmnVaryFluxNoBeam);
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
