#include <unity.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "test_uvw_coords.h"

/*Unity needs these calls, stick them here and leave them alone */
void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

/*
CUDA code we are linking in
*/
extern void test_kern_calc_uvw(float *X_diff, float *Y_diff, float *Z_diff,
           float *u_metres, float *v_metres, float *w_metres,
           float *us, float *vs, float *ws, float *wavelengths,
           float dec0,
           float *cha0s, float *sha0s,
           int num_visis, int num_baselines);

extern void test_kern_calc_uvw_shapelet(float *X_diff, float *Y_diff, float *Z_diff,
           float *u_shapes, float *v_shapes, float *w_shapes, float *wavelengths,
           float *lsts, double *ras, double *decs,
           int num_baselines, int num_visis, int num_shapes);


#define UNITY_INCLUDE_FLOAT

/*
Given the inputs, create simulation settings that woden.c would create
*/
void setup_uvw_params(int num_times, int num_baselines, int num_freqs,
                      float ra0,
                      float lst_base, float time_res,
                      float freq_res, float base_band_freq,
                      uvw_settings_t *uvw_settings) {

  for (int i = 0; i < num_baselines; i++) {
    uvw_settings->X_diff[i] = i + 1;
    uvw_settings->Y_diff[i] = i + 1;
    uvw_settings->Z_diff[i] = i + 1;
  }

  float *lsts = malloc(num_times*sizeof(float));

  for ( int time_step = 0; time_step < num_times; time_step++ ) {
    float lst = lst_base + time_step*time_res*SOLAR2SIDEREAL*DS2R;

    //Add half a time_res so we are sampling centre of each time step
    //WODEN would do this, but in a test below we want HA to be exactly zero,
    //so don't add this half time step here
    // lst += 0.5*time_res*SOLAR2SIDEREAL*DS2R;
    lsts[time_step] = lst;
  }

  float ha0, sha0, cha0, frequency, wavelength;

  for ( int time_step = 0; time_step < num_times; time_step++ ) {
    ha0 = lsts[time_step] - ra0;
    sha0 = sin(ha0);
    cha0 = cos(ha0);

    for (int freq_step = 0; freq_step < num_freqs; freq_step++) {
      frequency = base_band_freq + (freq_res*freq_step);
      wavelength = VELC / frequency;
      int step = num_baselines*(time_step*num_freqs + freq_step);

      for (int baseline = 0; baseline < num_baselines; baseline++) {
        uvw_settings->cha0s[step + baseline] = cha0;
        uvw_settings->sha0s[step + baseline] = sha0;
        uvw_settings->wavelengths[step + baseline] = wavelength;
        uvw_settings->lsts[step + baseline] = lsts[time_step];
      }//baseline loop
    }//freq loop
  }//time loop

  free(lsts);
}

/*
Do a bunch of mallocing
*/
uvw_settings_t * setup_uvw_settings(int num_baselines, int num_visis,
                                    int num_times, int num_components){

  uvw_settings_t * uvw_settings;
  uvw_settings = malloc( sizeof(uvw_settings_t) );

  uvw_settings->X_diff = malloc(num_baselines*sizeof(float));
  uvw_settings->Y_diff = malloc(num_baselines*sizeof(float));
  uvw_settings->Z_diff = malloc(num_baselines*sizeof(float));

  uvw_settings->lsts = malloc(num_times*sizeof(float));
  uvw_settings->lsts = malloc(num_visis*sizeof(float));
  uvw_settings->wavelengths = malloc(num_visis*sizeof(float));
  uvw_settings->cha0s = malloc(num_visis*sizeof(float));
  uvw_settings->sha0s = malloc(num_visis*sizeof(float));

  //kern_calc_uvw_shapelet calcuates u,v,w for each visibility and every
  //component
  uvw_settings->us = malloc(num_components*num_visis*sizeof(float));
  uvw_settings->vs = malloc(num_components*num_visis*sizeof(float));
  uvw_settings->ws = malloc(num_components*num_visis*sizeof(float));

  uvw_settings->u_metres = malloc(num_visis*sizeof(float));
  uvw_settings->v_metres = malloc(num_visis*sizeof(float));
  uvw_settings->w_metres = malloc(num_visis*sizeof(float));

  return uvw_settings;

}

/*
Do a bunch of freeing
*/
void free_uvw_settings(uvw_settings_t * uvw_settings){

  free(uvw_settings->X_diff);
  free(uvw_settings->Y_diff);
  free(uvw_settings->Z_diff);

  free(uvw_settings->lsts);
  free(uvw_settings->wavelengths);
  free(uvw_settings->cha0s);
  free(uvw_settings->sha0s);

  free(uvw_settings->us);
  free(uvw_settings->vs);
  free(uvw_settings->ws);

  free(uvw_settings->u_metres);
  free(uvw_settings->v_metres);
  free(uvw_settings->w_metres);

  free(uvw_settings);

}

/*
Checking the function fundamental_coords.cu::kern_calc_uvw
Checks that the wavelength scaling of u,v,w is happening correctly. Set HA=0
to make checking easier
*/
void test_kern_calc_uvw_ScalesByWavelength(void){

    //Setup some observation settings
    float ra0 = 0.0*DD2R;
    float dec0 = -30.0*DD2R;

    int num_baselines = 5;
    int num_freqs = 5;
    int num_times = 1;
    int num_visis = num_baselines*num_freqs*num_times;

    float lst_base = 0.0;
    float time_res = 8.0;
    float freq_res = 50e+6;
    float base_band_freq = 300e+6;

    //Something to hold all the uvw input parameter arrays
    uvw_settings_t *uvw_settings;
    //Shapelet uvw needs multiple components, here just set to one
    int num_components = 1;
    uvw_settings = setup_uvw_settings(num_baselines, num_visis,
                                      num_times, num_components);

    //Create input parameters for kern_calc_uvw
    setup_uvw_params(num_times, num_baselines, num_freqs, ra0,
                    lst_base, time_res,
                    freq_res, base_band_freq,
                    uvw_settings);

    //Run the CUDA code via fundamental_coords::test_kern_calc_uvw
    test_kern_calc_uvw(uvw_settings->X_diff, uvw_settings->Y_diff, uvw_settings->Z_diff,
           uvw_settings->u_metres, uvw_settings->v_metres, uvw_settings->w_metres,
           uvw_settings->us, uvw_settings->vs, uvw_settings->ws, uvw_settings->wavelengths,
           dec0, uvw_settings->cha0s, uvw_settings->sha0s,
           num_visis, num_baselines);

    //Create expected values
    float *u_metres_expec = malloc(num_visis*sizeof(float));
    float *v_metres_expec = malloc(num_visis*sizeof(float));
    float *w_metres_expec = malloc(num_visis*sizeof(float));
    float *us_expec = malloc(num_visis*sizeof(float));
    float *vs_expec = malloc(num_visis*sizeof(float));
    float *ws_expec = malloc(num_visis*sizeof(float));

    float cdec0 = cos(dec0);
    float sdec0 = sin(dec0);
    //Special case of expected u,v,w when phase centre is at HA = 0.0
    int index = 0;
    for ( int time_step = 0; time_step < num_times; time_step++ ) {
      for (int freq_step = 0; freq_step < num_freqs; freq_step++) {
        for (int baseline = 0; baseline < num_baselines; baseline++) {

          u_metres_expec[index] = uvw_settings->Y_diff[baseline];
          v_metres_expec[index] = -sdec0*uvw_settings->X_diff[baseline] + cdec0*uvw_settings->Y_diff[baseline];
          w_metres_expec[index] = cdec0*uvw_settings->X_diff[baseline] + sdec0*uvw_settings->Y_diff[baseline];

          us_expec[index] = u_metres_expec[index] / uvw_settings->wavelengths[index];
          vs_expec[index] = v_metres_expec[index] / uvw_settings->wavelengths[index];
          ws_expec[index] = w_metres_expec[index] / uvw_settings->wavelengths[index];

          index += 1;

      }//baseline loop
    }//freq loop
  }//time loop

  //Check they are equal
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(u_metres_expec, uvw_settings->u_metres, num_visis);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(v_metres_expec, uvw_settings->v_metres, num_visis);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(w_metres_expec, uvw_settings->w_metres, num_visis);

  TEST_ASSERT_EQUAL_FLOAT_ARRAY(us_expec, uvw_settings->us, num_visis);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(vs_expec, uvw_settings->vs, num_visis);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(ws_expec, uvw_settings->ws, num_visis);

  free_uvw_settings(uvw_settings);

}


/*
Checking the function fundamental_coords.cu::kern_calc_uvw
Checks that u,v,w coords change with time as expected
Make checking easier by setting dec phase centre dec0=0.0
*/
void test_kern_calc_uvw_RotateWithTime(void){

    //Setup some observation settings
    float ra0 = 0.0*DD2R;
    float dec0 = 0.0*DD2R;

    int num_baselines = 5;
    int num_freqs = 1;
    int num_times = 5;
    int num_visis = num_baselines*num_freqs*num_times;

    float lst_base = 0.0;
    float time_res = 8.0;
    float freq_res = 50e+6;
    float base_band_freq = 300e+6;

    //Something to hold all the uvw input parameter arrays
    uvw_settings_t *uvw_settings;
    //Shapelet uvw needs multiple components, here just set to one
    int num_components = 1;
    uvw_settings = setup_uvw_settings(num_baselines, num_visis,
                                      num_times, num_components);

    //Create input parameters for kern_calc_uvw
    //Shapelet uvw needs multiple components, here just set to one
    setup_uvw_params(num_times, num_baselines, num_freqs, ra0,
                    lst_base, time_res,
                    freq_res, base_band_freq,
                    uvw_settings);

    //Run the CUDA code via fundamental_coords::test_kern_calc_uvw
    test_kern_calc_uvw(uvw_settings->X_diff, uvw_settings->Y_diff, uvw_settings->Z_diff,
           uvw_settings->u_metres, uvw_settings->v_metres, uvw_settings->w_metres,
           uvw_settings->us, uvw_settings->vs, uvw_settings->ws, uvw_settings->wavelengths,
           dec0, uvw_settings->cha0s, uvw_settings->sha0s,
           num_visis, num_baselines);

    //Create expected values
    float *u_metres_expec = malloc(num_visis*sizeof(float));
    float *v_metres_expec = malloc(num_visis*sizeof(float));
    float *w_metres_expec = malloc(num_visis*sizeof(float));
    float *us_expec = malloc(num_visis*sizeof(float));
    float *vs_expec = malloc(num_visis*sizeof(float));
    float *ws_expec = malloc(num_visis*sizeof(float));

    float cdec0 = cos(dec0);
    float sdec0 = sin(dec0);
    //Special case of expected u,v,w when phase centre is at Dec = 0.0
    int index = 0;
    for ( int time_step = 0; time_step < num_times; time_step++ ) {
      for (int freq_step = 0; freq_step < num_freqs; freq_step++) {
        for (int baseline = 0; baseline < num_baselines; baseline++) {

          u_metres_expec[index] = uvw_settings->sha0s[index]*uvw_settings->X_diff[baseline] + uvw_settings->cha0s[index]*uvw_settings->Y_diff[baseline];
          v_metres_expec[index] = uvw_settings->Z_diff[baseline];
          w_metres_expec[index] = uvw_settings->cha0s[index]*uvw_settings->X_diff[baseline] - uvw_settings->sha0s[index]*uvw_settings->Y_diff[baseline];

          us_expec[index] = u_metres_expec[index] / uvw_settings->wavelengths[index];
          vs_expec[index] = v_metres_expec[index] / uvw_settings->wavelengths[index];
          ws_expec[index] = w_metres_expec[index] / uvw_settings->wavelengths[index];

          index += 1;

      }//baseline loop
    }//freq loop
  }//time loop

  //Check they are equal
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(u_metres_expec, uvw_settings->u_metres, num_visis);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(v_metres_expec, uvw_settings->v_metres, num_visis);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(w_metres_expec, uvw_settings->w_metres, num_visis);

  TEST_ASSERT_EQUAL_FLOAT_ARRAY(us_expec, uvw_settings->us, num_visis);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(vs_expec, uvw_settings->vs, num_visis);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(ws_expec, uvw_settings->ws, num_visis);

  free_uvw_settings(uvw_settings);

}

/*
Checking the function fundamental_coords.cu::kern_calc_uvw_shapelet
Checks that u,v,w coords change with time as expected
Make checking easier by setting dec phase centre dec0=0.0
Also checks that results are scaled by wavelength correctly
*/
void test_kern_calc_uvw_shapelet_RotateWithTimeScalesByWavelength(void){

    //Setup some observation settings
    float ra0 = 0.0*DD2R;
    float dec0 = 0.0*DD2R;

    int num_baselines = 5;
    int num_freqs = 2;
    int num_times = 5;
    int num_visis = num_baselines*num_freqs*num_times;

    float lst_base = 0.0;
    float time_res = 8.0;
    float freq_res = 50e+6;
    float base_band_freq = 300e+6;

    //Check values for 5 different shapelet components, with different RA values
    int num_components = 5;

    double *ras = malloc(num_components*sizeof(double));
    double *decs = malloc(num_components*sizeof(double));

    for (size_t i = 0; i < num_components; i++) {
      ras[i] = i*DD2R;
      decs[i] = 0.0;
    }

    //Something to hold all the uvw input parameter arrays
    uvw_settings_t *uvw_settings;
    uvw_settings = setup_uvw_settings(num_baselines, num_visis,
                                      num_times, num_components);

    //Create input parameters for kern_calc_uvw
    setup_uvw_params(num_times, num_baselines, num_freqs, ra0,
                    lst_base, time_res,
                    freq_res, base_band_freq,
                    uvw_settings);

     //Run the CUDA code via fundamental_coords::test_kern_calc_uvw_shapelet
    test_kern_calc_uvw_shapelet(uvw_settings->X_diff, uvw_settings->Y_diff, uvw_settings->Z_diff,
           uvw_settings->us, uvw_settings->vs, uvw_settings->ws, uvw_settings->wavelengths,
           uvw_settings->lsts, ras, decs,
           num_baselines, num_visis, num_components);

    //Create expected values
    float *us_expec = malloc(num_components*num_visis*sizeof(float));
    float *vs_expec = malloc(num_components*num_visis*sizeof(float));
    float *ws_expec = malloc(num_components*num_visis*sizeof(float));

    //Variables to use when making expected values
    float u_metre, v_metre, w_metre;
    float ha, sha, cha;

    int uvw_index = 0;
    for ( int comp_step = 0; comp_step < num_times; comp_step++ ) {
      int visi_index = 0;
      for ( int time_step = 0; time_step < num_times; time_step++ ) {
        //For kern_calc_uvw_shapelet, calculating u,v,w centred on ra, dec of
        //each shapelet component. Calcualte ha using the shapelet ra
        ha = uvw_settings->lsts[visi_index] - ras[comp_step];
        sha = sinf(ha);
        cha = cosf(ha);

        for (int freq_step = 0; freq_step < num_freqs; freq_step++) {
          for (int baseline = 0; baseline < num_baselines; baseline++) {
            //Special case of expected u,v,w when u,v,w coord centre is at Dec = 0.0
            u_metre = sha*uvw_settings->X_diff[baseline] + cha*uvw_settings->Y_diff[baseline];
            v_metre = uvw_settings->Z_diff[baseline];
            w_metre = cha*uvw_settings->X_diff[baseline] - sha*uvw_settings->Y_diff[baseline];

            us_expec[uvw_index] = u_metre / uvw_settings->wavelengths[visi_index];
            vs_expec[uvw_index] = v_metre / uvw_settings->wavelengths[visi_index];
            ws_expec[uvw_index] = w_metre / uvw_settings->wavelengths[visi_index];

            uvw_index += 1;
            visi_index += 1;

        }//baseline loop
      }//freq loop
    }//time loop
  }//component loop
  //
  // //Check they are equal
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(us_expec, uvw_settings->us, num_visis*num_components);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(vs_expec, uvw_settings->vs, num_visis);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(ws_expec, uvw_settings->ws, num_visis);
  //
  free_uvw_settings(uvw_settings);

}

int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_kern_calc_uvw_ScalesByWavelength);
    RUN_TEST(test_kern_calc_uvw_RotateWithTime);
    RUN_TEST(test_kern_calc_uvw_shapelet_RotateWithTimeScalesByWavelength);

    return UNITY_END();
}
