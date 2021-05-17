#include <math.h>
#include <unity.h>
#include <stdlib.h>
#include <complex.h>

#include "constants.h"
#include "woden_struct_defs.h"

void setUp (void) {} /* Is run before eVary test, put unit init calls here. */
void tearDown (void) {} /* Is run after eVary test, put unit clean-up calls here. */

//External CUDA code we're linking in
extern void test_kern_calc_visi_point(int num_components, int num_baselines,
          int num_freqs, int num_visis, int num_times, int beamtype,
          float *component_freqs,
          float *flux_I, float *flux_Q, float *flux_U, float *flux_V,
          float *SIs, float *us, float *vs, float *ws,
          float _Complex *sum_visi_XX_real, float _Complex *sum_visi_XX_imag,
          float _Complex *sum_visi_XY_real, float _Complex *sum_visi_XY_imag,
          float _Complex *sum_visi_YX_real, float _Complex *sum_visi_YX_imag,
          float _Complex *sum_visi_YY_real, float _Complex *sum_visi_YY_imag,
          float *allsteps_wavelengths,
          float *ls, float *ms, float *ns,
          float _Complex *primay_beam_J00, float _Complex *primay_beam_J01,
          float _Complex *primay_beam_J10, float _Complex *primay_beam_J11);

#define UNITY_INCLUDE_FLOAT

/*
Test the __device__ code that updates the summed visibilities by grabbing the
correct beam gain and mesurement equation, multiplying and summing onto the visi
Here we keep the component visibilities and fluxes constant and vary the beam gains
Test works for all primary beam types
*/
void test_kern_calc_visi_point_Varylmn(int beamtype) {

  int num_baselines = 3.0;
  int num_times = 1.0;
  int num_freqs = 3.0;

  int num_visis = num_baselines*num_times*num_freqs;

  int num_components = 25;

  int num_points = 9;
  float *ls = malloc(num_components*num_times*sizeof(float));
  float *ms = malloc(num_components*num_times*sizeof(float));
  float *ns = malloc(num_components*num_times*sizeof(float));

  int count = 0;
  for (int l_ind = 0; l_ind < 5; l_ind++) {
    for (int m_ind = 0; m_ind < 5; m_ind++) {
      ls[count] = -0.5 + 0.25*l_ind;
      ms[count] = -0.5 + 0.25*m_ind;
      ns[count] = sqrt(1 - ls[count]*ls[count] - ms[count]*ms[count]);
    }
  }

  float *us = malloc(num_visis*sizeof(float));
  float *vs = malloc(num_visis*sizeof(float));
  float *ws = malloc(num_visis*sizeof(float));
  count = 0;
  for (int visi = 0; visi < num_visis; visi++) {

  }

  num_visis

  float _Complex *primay_beam_J00 = malloc(num_freqs*num_times*num_components*sizeof(float _Complex));
  float _Complex *primay_beam_J01 = malloc(num_freqs*num_times*num_components*sizeof(float _Complex));
  float _Complex *primay_beam_J10 = malloc(num_freqs*num_times*num_components*sizeof(float _Complex));
  float _Complex *primay_beam_J11 = malloc(num_freqs*num_times*num_components*sizeof(float _Complex));

  //Fill in some varying beam gains based on beam type
  for (size_t visi = 0; visi < num_components*num_times*num_freqs; visi++) {
    primay_beam_J00[visi] = 1.0 + I*0.0;
    primay_beam_J11[visi] = 1.0 + I*0.0;
    primay_beam_J01[visi] = 1.0 + I*0.0;
    primay_beam_J10[visi] = 1.0 + I*0.0;
  }

  //Set all the Stokes I to 1.0 and other stokes to zero
  float *flux_I = malloc(num_components*sizeof(float));
  float *flux_Q = calloc(num_components, sizeof(float));
  float *flux_U = calloc(num_components, sizeof(float));
  float *flux_V = calloc(num_components, sizeof(float));
  float *SIs = malloc(num_components*sizeof(float));
  float *component_freqs = malloc(num_components*sizeof(float));

  //Just stick Stokes I to 1.0, SI to zero, and reference freqs to 150MHz
  for (size_t visi = 0; visi < num_freqs*num_components; visi++) {
    flux_I[visi] = 1.0;
    SIs[visi] = 0.0;
    component_freqs = 150e+6;
  }

  //Make sure arrays to hold summed visis are initialised to zero
  float *sum_visi_XX_real = calloc(num_visis, sizeof(float));
  float *sum_visi_XX_imag = calloc(num_visis, sizeof(float));
  float *sum_visi_XY_real = calloc(num_visis, sizeof(float));
  float *sum_visi_XY_imag = calloc(num_visis, sizeof(float));
  float *sum_visi_YX_real = calloc(num_visis, sizeof(float));
  float *sum_visi_YX_imag = calloc(num_visis, sizeof(float));
  float *sum_visi_YY_real = calloc(num_visis, sizeof(float));
  float *sum_visi_YY_imag = calloc(num_visis, sizeof(float));

  //Run the CUDA code
  test_kern_update_sum_visis(num_freqs, num_visis,
          num_baselines, num_components, num_times, beamtype,
          primay_beam_J00, primay_beam_J01,
          primay_beam_J10, primay_beam_J11,
          visi_components,
          flux_I, flux_Q, flux_U, flux_V,
          sum_visi_XX_real, sum_visi_XX_imag,
          sum_visi_XY_real, sum_visi_XY_imag,
          sum_visi_YX_real, sum_visi_YX_imag,
          sum_visi_YY_real, sum_visi_YY_imag);

  test_kern_calc_visi_point(num_components, num_baselines,
          num_freqs, num_visis, num_times, beamtype,
          component_freqs,
          flux_I, flux_Q, flux_U, flux_V,
          SIs, float *us, float *vs, float *ws,
          sum_visi_XX_real, sum_visi_XX_imag,
          sum_visi_XY_real, sum_visi_XY_imag,
          sum_visi_YX_real, sum_visi_YX_imag,
          sum_visi_YY_real, sum_visi_YY_imag,
          float *allsteps_wavelengths,
          float *ls, float *ms, float *ns,
          float _Complex *primay_beam_J00, float _Complex *primay_beam_J01,
          float _Complex *primay_beam_J10, float _Complex *primay_beam_J11)

  if (beamtype == FEE_BEAM) {
    for (size_t output = 0; output < num_visis; output++) {

      // printf("%.1f %.1f %.1f %.1f %.1f %.1f %.1f %.1f\n",
      //     sum_visi_XX_real[output], sum_visi_XX_imag[output],
      //     sum_visi_XY_real[output], sum_visi_XY_imag[output],
      //     sum_visi_YX_real[output], sum_visi_YX_imag[output],
      //     sum_visi_YY_real[output], sum_visi_YY_imag[output]);

      TEST_ASSERT_EQUAL_FLOAT(expected_order[output], sum_visi_XX_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(expected_order[output], sum_visi_XY_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(expected_order[output], sum_visi_YX_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(expected_order[output], sum_visi_YY_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_XX_imag[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_XY_imag[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_YX_imag[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_YY_imag[output]);

    }
  }
  else {
    //If there is no beam, should just end up with the number of components
    //as we apply a gain of 1.0
    for (size_t output = 0; output < num_visis; output++) {

      TEST_ASSERT_EQUAL_FLOAT(num_components, sum_visi_XX_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_XY_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_YX_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(num_components, sum_visi_YY_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_XX_imag[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_XY_imag[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_YX_imag[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_YY_imag[output]);

    }
  }
  //
  // //Be free my beauties
  free( primay_beam_J00 );
  free( primay_beam_J01 );
  free( primay_beam_J10 );
  free( primay_beam_J11 );
  free( visi_components );
  free( flux_I );
  free( flux_Q );
  free( flux_U );
  free( flux_V );
  free( sum_visi_XX_real );
  free( sum_visi_XX_imag );
  free( sum_visi_XY_real );
  free( sum_visi_XY_imag );
  free( sum_visi_YX_real );
  free( sum_visi_YX_imag );
  free( sum_visi_YY_real );
  free( sum_visi_YY_imag );

}

/*
This test checks varying the gain with beamtype=FEE_BEAM
*/
void test_kern_update_sum_visis_VaryGainFEEBeam(void) {
  test_kern_update_sum_visis_VaryGainChooseBeams(FEE_BEAM);
}

/*
This test checks varying the gain with beamtype=ANALY_DIPOLE
*/
void test_kern_update_sum_visis_VaryGainAnalyBeam(void) {
  test_kern_update_sum_visis_VaryGainChooseBeams(ANALY_DIPOLE);
}

/*
This test checks varying the gain with beamtype=GAUSS_BEAM
*/
void test_kern_update_sum_visis_VaryGainGaussBeam(void) {
  test_kern_update_sum_visis_VaryGainChooseBeams(GAUSS_BEAM);
}

/*
This test checks varying the gain with beamtype=NO_BEAM
*/
void test_kern_update_sum_visis_VaryGainNoBeam(void) {
  test_kern_update_sum_visis_VaryGainChooseBeams(NO_BEAM);
}

/*
Test the __device__ code that updates the summed visibilities by grabbing the
correct beam gain and mesurement equation, multiplying and summing onto the visi
Here we keep the component visibilities and the beam gains constant and vary fluxes
Test works for all primary beam types
*/
void test_kern_update_sum_visis_VaryFluxesChooseBeams(int beamtype) {

  int num_baselines = 3.0;
  int num_times = 4.0;
  int num_freqs = 3.0;

  int num_visis = num_baselines*num_times*num_freqs;

  int num_components = 10.0;

  float _Complex *primay_beam_J00 = malloc(num_freqs*num_times*num_components*sizeof(float _Complex));
  float _Complex *primay_beam_J01 = malloc(num_freqs*num_times*num_components*sizeof(float _Complex));
  float _Complex *primay_beam_J10 = malloc(num_freqs*num_times*num_components*sizeof(float _Complex));
  float _Complex *primay_beam_J11 = malloc(num_freqs*num_times*num_components*sizeof(float _Complex));
  float _Complex *visi_components = malloc(num_visis*sizeof(float _Complex));

  //Just stick base measurement equation to 1.0
  for (size_t visi = 0; visi < num_visis; visi++) {
    visi_components[visi] = 1.0 + I*0.0;
  }

  //Set all beam gains to 1.0
  for (size_t visi = 0; visi < num_components*num_times*num_freqs; visi++) {
    primay_beam_J00[visi] = 1.0 + I*0.0;
    primay_beam_J01[visi] = 1.0 + I*0.0;
    primay_beam_J10[visi] = 1.0 + I*0.0;
    primay_beam_J11[visi] = 1.0 + I*0.0;
  }

  //Vary Stokes parameters
  float *flux_I = malloc(num_freqs*num_components*sizeof(float));
  float *flux_Q = malloc(num_freqs*num_components*sizeof(float));
  float *flux_U = malloc(num_freqs*num_components*sizeof(float));
  float *flux_V = malloc(num_freqs*num_components*sizeof(float));

  //Just stick base measurement equation to 1.0
  int count = 1;

  for (size_t visi = 0; visi < num_freqs*num_components; visi++) {
    flux_I[visi] = count;
    flux_Q[visi] = 0;
    flux_U[visi] = 0;
    flux_V[visi] = 0;
    // printf("%.3f\n",flux_I[visi] );
    count ++;
  }

  //Make sure arrays to hold summed visis are initialised to zero
  float *sum_visi_XX_real = calloc(num_visis, sizeof(float));
  float *sum_visi_XX_imag = calloc(num_visis, sizeof(float));
  float *sum_visi_XY_real = calloc(num_visis, sizeof(float));
  float *sum_visi_XY_imag = calloc(num_visis, sizeof(float));
  float *sum_visi_YX_real = calloc(num_visis, sizeof(float));
  float *sum_visi_YX_imag = calloc(num_visis, sizeof(float));
  float *sum_visi_YY_real = calloc(num_visis, sizeof(float));
  float *sum_visi_YY_imag = calloc(num_visis, sizeof(float));

  //Run the CUDA code
  test_kern_update_sum_visis(num_freqs, num_visis,
          num_baselines, num_components, num_times, beamtype,
          primay_beam_J00, primay_beam_J01,
          primay_beam_J10, primay_beam_J11,
          visi_components,
          flux_I, flux_Q, flux_U, flux_V,
          sum_visi_XX_real, sum_visi_XX_imag,
          sum_visi_XY_real, sum_visi_XY_imag,
          sum_visi_YX_real, sum_visi_YX_imag,
          sum_visi_YY_real, sum_visi_YY_imag);

  float *expected_order = malloc(num_visis*sizeof(float));

  //These expected values have no cross-pols in them
  float freq_sum[] = {55.0, 155.0, 255.0};

  int ind = 0;
  for (size_t time_ind = 0; time_ind < num_times; time_ind++) {
    for (size_t freq = 0; freq < num_freqs; freq++) {
      for (size_t baseline = 0; baseline < num_baselines; baseline++) {
        expected_order[ind] = freq_sum[freq];
        ind ++;
      }
    }
  }

  if (beamtype == ANALY_DIPOLE || beamtype == GAUSS_BEAM) {
    for (size_t output = 0; output < num_visis; output++) {

      // printf("%.1f %.1f %.1f %.1f %.1f %.1f %.1f %.1f\n",
      //     sum_visi_XX_real[output], sum_visi_XX_imag[output],
      //     sum_visi_XY_real[output], sum_visi_XY_imag[output],
      //     sum_visi_YX_real[output], sum_visi_YX_imag[output],
      //     sum_visi_YY_real[output], sum_visi_YY_imag[output]);

      TEST_ASSERT_EQUAL_FLOAT(expected_order[output], sum_visi_XX_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_XY_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_YX_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(expected_order[output], sum_visi_YY_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_XX_imag[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_XY_imag[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_YX_imag[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_YY_imag[output]);

    }
  }
  //No cross-pols so multiply expected values by 2.0
  else if (beamtype == FEE_BEAM) {
    for (size_t output = 0; output < num_visis; output++) {

      TEST_ASSERT_EQUAL_FLOAT(2*expected_order[output], sum_visi_XX_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(2*expected_order[output], sum_visi_XY_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(2*expected_order[output], sum_visi_YX_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(2*expected_order[output], sum_visi_YY_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_XX_imag[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_XY_imag[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_YX_imag[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_YY_imag[output]);

    }
  }
  else {
    //If there is no beam, should just end up with the number of components
    //as we apply a gain of 1.0
    for (size_t output = 0; output < num_visis; output++) {

      TEST_ASSERT_EQUAL_FLOAT(expected_order[output], sum_visi_XX_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_XY_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_YX_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(expected_order[output], sum_visi_YY_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_XX_imag[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_XY_imag[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_YX_imag[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_YY_imag[output]);

    }
  }
  //
  // //Be free my beauties
  free( primay_beam_J00 );
  free( primay_beam_J01 );
  free( primay_beam_J10 );
  free( primay_beam_J11 );
  free( visi_components );
  free( flux_I );
  free( flux_Q );
  free( flux_U );
  free( flux_V );
  free( sum_visi_XX_real );
  free( sum_visi_XX_imag );
  free( sum_visi_XY_real );
  free( sum_visi_XY_imag );
  free( sum_visi_YX_real );
  free( sum_visi_YX_imag );
  free( sum_visi_YY_real );
  free( sum_visi_YY_imag );

}

/*
This test checks varying the gain with beamtype=FEE_BEAM
*/
void test_kern_update_sum_visis_VaryFluxesFEEBeam(void) {
  test_kern_update_sum_visis_VaryFluxesChooseBeams(FEE_BEAM);
}

/*
This test checks varying the gain with beamtype=ANALY_DIPOLE
*/
void test_kern_update_sum_visis_VaryFluxesAnalyBeam(void) {
  test_kern_update_sum_visis_VaryFluxesChooseBeams(ANALY_DIPOLE);
}

/*
This test checks varying the gain with beamtype=GAUSS_BEAM
*/
void test_kern_update_sum_visis_VaryFluxesGaussBeam(void) {
  test_kern_update_sum_visis_VaryFluxesChooseBeams(GAUSS_BEAM);
}

/*
This test checks varying the gain with beamtype=NO_BEAM
*/
void test_kern_update_sum_visis_VaryFluxesNoBeam(void) {
  test_kern_update_sum_visis_VaryFluxesChooseBeams(NO_BEAM);
}




/*
Test the __device__ code that updates the summed visibilities by grabbing the
correct beam gain and mesurement equation, multiplying and summing onto the visi
Here we keep the fluxes and beam gains constant and vary the component visibilities
Test works for all primary beam types
*/
void test_kern_update_sum_visis_VaryVisiChooseBeams(int beamtype) {

  int num_baselines = 3.0;
  int num_times = 4.0;
  int num_freqs = 3.0;

  int num_visis = num_baselines*num_times*num_freqs;

  int num_components = 10.0;

  float _Complex *primay_beam_J00 = malloc(num_freqs*num_times*num_components*sizeof(float _Complex));
  float _Complex *primay_beam_J01 = malloc(num_freqs*num_times*num_components*sizeof(float _Complex));
  float _Complex *primay_beam_J10 = malloc(num_freqs*num_times*num_components*sizeof(float _Complex));
  float _Complex *primay_beam_J11 = malloc(num_freqs*num_times*num_components*sizeof(float _Complex));
  float _Complex *visi_components = malloc(num_visis*sizeof(float _Complex));

  //Vary the base visibilities
  int count = 1;
  for (size_t visi = 0; visi < num_visis; visi++) {
    visi_components[visi] = count + I*count;
    count ++;
  }

  //Set all beam gains to 1.0
  for (size_t visi = 0; visi < num_components*num_times*num_freqs; visi++) {
    primay_beam_J00[visi] = 1.0 + I*0.0;
    primay_beam_J01[visi] = 1.0 + I*0.0;
    primay_beam_J10[visi] = 1.0 + I*0.0;
    primay_beam_J11[visi] = 1.0 + I*0.0;
  }

  //Vary Stokes parameters
  float *flux_I = malloc(num_freqs*num_components*sizeof(float));
  float *flux_Q = calloc(num_freqs, num_components*sizeof(float));
  float *flux_U = calloc(num_freqs, num_components*sizeof(float));
  float *flux_V = calloc(num_freqs, num_components*sizeof(float));

  //Just stick base Stokes I to 1.0
  for (size_t visi = 0; visi < num_freqs*num_components; visi++) {
    flux_I[visi] = 1.0;
  }

  //Make sure arrays to hold summed visis are initialised to zero
  float *sum_visi_XX_real = calloc(num_visis, sizeof(float));
  float *sum_visi_XX_imag = calloc(num_visis, sizeof(float));
  float *sum_visi_XY_real = calloc(num_visis, sizeof(float));
  float *sum_visi_XY_imag = calloc(num_visis, sizeof(float));
  float *sum_visi_YX_real = calloc(num_visis, sizeof(float));
  float *sum_visi_YX_imag = calloc(num_visis, sizeof(float));
  float *sum_visi_YY_real = calloc(num_visis, sizeof(float));
  float *sum_visi_YY_imag = calloc(num_visis, sizeof(float));

  //Run the CUDA code
  test_kern_update_sum_visis(num_freqs, num_visis,
          num_baselines, num_components, num_times, beamtype,
          primay_beam_J00, primay_beam_J01,
          primay_beam_J10, primay_beam_J11,
          visi_components,
          flux_I, flux_Q, flux_U, flux_V,
          sum_visi_XX_real, sum_visi_XX_imag,
          sum_visi_XY_real, sum_visi_XY_imag,
          sum_visi_YX_real, sum_visi_YX_imag,
          sum_visi_YY_real, sum_visi_YY_imag);

  //These expected values have no cross-pols in them
  float *expected_order = malloc(num_visis*sizeof(float));
  for (size_t visi = 0; visi < num_visis; visi++) {
    expected_order[visi] = (visi + 1)*num_components;
  }

  if (beamtype == ANALY_DIPOLE || beamtype == GAUSS_BEAM) {
    for (size_t output = 0; output < num_visis; output++) {

      // printf("%.1f %.1f %.1f %.1f %.1f %.1f %.1f %.1f\n",
      //     sum_visi_XX_real[output], sum_visi_XX_imag[output],
      //     sum_visi_XY_real[output], sum_visi_XY_imag[output],
      //     sum_visi_YX_real[output], sum_visi_YX_imag[output],
      //     sum_visi_YY_real[output], sum_visi_YY_imag[output]);

      TEST_ASSERT_EQUAL_FLOAT(expected_order[output], sum_visi_XX_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_XY_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_YX_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(expected_order[output], sum_visi_YY_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(expected_order[output], sum_visi_XX_imag[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_XY_imag[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_YX_imag[output]);
      TEST_ASSERT_EQUAL_FLOAT(expected_order[output], sum_visi_YY_imag[output]);

    }
  }
  //No cross-pols so multiply expected values by 2.0
  else if (beamtype == FEE_BEAM) {
    for (size_t output = 0; output < num_visis; output++) {

      TEST_ASSERT_EQUAL_FLOAT(2*expected_order[output], sum_visi_XX_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(2*expected_order[output], sum_visi_XY_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(2*expected_order[output], sum_visi_YX_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(2*expected_order[output], sum_visi_YY_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(2*expected_order[output], sum_visi_XX_imag[output]);
      TEST_ASSERT_EQUAL_FLOAT(2*expected_order[output], sum_visi_XY_imag[output]);
      TEST_ASSERT_EQUAL_FLOAT(2*expected_order[output], sum_visi_YX_imag[output]);
      TEST_ASSERT_EQUAL_FLOAT(2*expected_order[output], sum_visi_YY_imag[output]);

    }
  }
  else {
    //If there is no beam, should just end up with the number of components
    //as we apply a gain of 1.0
    for (size_t output = 0; output < num_visis; output++) {

      TEST_ASSERT_EQUAL_FLOAT(expected_order[output], sum_visi_XX_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_XY_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_YX_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(expected_order[output], sum_visi_YY_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(expected_order[output], sum_visi_XX_imag[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_XY_imag[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_YX_imag[output]);
      TEST_ASSERT_EQUAL_FLOAT(expected_order[output], sum_visi_YY_imag[output]);

    }
  }
  //
  // //Be free my beauties
  free( primay_beam_J00 );
  free( primay_beam_J01 );
  free( primay_beam_J10 );
  free( primay_beam_J11 );
  free( visi_components );
  free( flux_I );
  free( flux_Q );
  free( flux_U );
  free( flux_V );
  free( sum_visi_XX_real );
  free( sum_visi_XX_imag );
  free( sum_visi_XY_real );
  free( sum_visi_XY_imag );
  free( sum_visi_YX_real );
  free( sum_visi_YX_imag );
  free( sum_visi_YY_real );
  free( sum_visi_YY_imag );

}

/*
This test checks varying the measurement equation with beamtype=FEE_BEAM
*/
void test_kern_update_sum_visis_VaryVisiFEEBeam(void) {
  test_kern_update_sum_visis_VaryVisiChooseBeams(FEE_BEAM);
}

/*
This test checks varying the measurement equation with beamtype=ANALY_DIPOLE
*/
void test_kern_update_sum_visis_VaryVisiAnalyBeam(void) {
  test_kern_update_sum_visis_VaryVisiChooseBeams(ANALY_DIPOLE);
}

/*
This test checks varying the measurement equation with beamtype=GAUSS_BEAM
*/
void test_kern_update_sum_visis_VaryVisiGaussBeam(void) {
  test_kern_update_sum_visis_VaryVisiChooseBeams(GAUSS_BEAM);
}

/*
This test checks varying the measurement equation with beamtype=NO_BEAM
*/
void test_kern_update_sum_visis_VaryVisiNoBeam(void) {
  test_kern_update_sum_visis_VaryVisiChooseBeams(NO_BEAM);
}

//Run the test with unity
int main(void)
{
    UNITY_BEGIN();
    //Test while varying beam gain for all beam types
    RUN_TEST(test_kern_update_sum_visis_VaryGainFEEBeam);
    RUN_TEST(test_kern_update_sum_visis_VaryGainAnalyBeam);
    RUN_TEST(test_kern_update_sum_visis_VaryGainGaussBeam);
    RUN_TEST(test_kern_update_sum_visis_VaryGainNoBeam);

    //Test while varying component fluxes for all beam types
    RUN_TEST(test_kern_update_sum_visis_VaryFluxesFEEBeam);
    RUN_TEST(test_kern_update_sum_visis_VaryFluxesAnalyBeam);
    RUN_TEST(test_kern_update_sum_visis_VaryFluxesNoBeam);
    RUN_TEST(test_kern_update_sum_visis_VaryFluxesGaussBeam);

    //Test while varying base visibility for all beam types
    RUN_TEST(test_kern_update_sum_visis_VaryVisiFEEBeam);
    RUN_TEST(test_kern_update_sum_visis_VaryVisiAnalyBeam);
    RUN_TEST(test_kern_update_sum_visis_VaryVisiNoBeam);
    RUN_TEST(test_kern_update_sum_visis_VaryVisiGaussBeam);

    return UNITY_END();
}
