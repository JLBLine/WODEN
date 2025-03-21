#include "uvw_coords_common.h"
#include <unity.h>


/*
CUDA code we are linking in
*/
extern void test_calc_uvw_gpu(double *X_diff, double *Y_diff, double *Z_diff,
   user_precision_t *u_metres, user_precision_t *v_metres, user_precision_t *w_metres,
   user_precision_t *us, user_precision_t *vs, user_precision_t *ws,
   user_precision_t *wavelengths, double *cha0s, double *sha0s,
   woden_settings_t *woden_settings);

extern void test_calc_uv_shapelet_gpu(double *X_diff, double *Y_diff, double *Z_diff,
                     user_precision_t *u_shapes, user_precision_t *v_shapes,
                     double *lsts, double *ras, double *decs,
                     int num_baselines, int num_times, int num_shapes);

/*
Given the inputs, create simulation settings that woden.c would create
*/
void setup_uvw_params(int num_times, int num_baselines, int num_freqs,
                      user_precision_t ra0,
                      user_precision_t lst_base, user_precision_t time_res,
                      user_precision_t freq_res, user_precision_t base_band_freq,
                      uvw_settings_t *uvw_settings) {

  for (int time_ind = 0; time_ind < num_times; time_ind++) {
    for (int i = 0; i < num_baselines; i++) {
      int time_off = time_ind*num_baselines;
      uvw_settings->X_diff[time_off + i] = (i + 1)*1000.0;
      uvw_settings->Y_diff[time_off + i] = (i + 1)*1000.0;
      uvw_settings->Z_diff[time_off + i] = (i + 1)*1000.0;
      // printf("index %d value %.1f\n", time_off + i, (i + 1)*1000.0);
    }
  }


  double *lsts = malloc(num_times*sizeof(double));

  for ( int time_step = 0; time_step < num_times; time_step++ ) {
    double lst = lst_base + time_step*time_res*SOLAR2SIDEREAL*DS2R;

    //Add half a time_res so we are sampling centre of each time step
    //WODEN would do this, but in a test below we want HA to be exactly zero,
    //so don't add this half time step here
    // lst += 0.5*time_res*SOLAR2SIDEREAL*DS2R;
    lsts[time_step] = lst;
  }

  double ha0, sha0, cha0;
  double frequency;
  user_precision_t wavelength;

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

  uvw_settings->X_diff = malloc(num_times*num_baselines*sizeof(double));
  uvw_settings->Y_diff = malloc(num_times*num_baselines*sizeof(double));
  uvw_settings->Z_diff = malloc(num_times*num_baselines*sizeof(double));

  uvw_settings->lsts = malloc(num_visis*sizeof(double));
  uvw_settings->wavelengths = malloc(num_visis*sizeof(user_precision_t));
  uvw_settings->cha0s = malloc(num_visis*sizeof(double));
  uvw_settings->sha0s = malloc(num_visis*sizeof(double));

  //kern_calc_uvw_shapelet calcuates u,v,w for each visibility and every
  //component
  uvw_settings->us = malloc(num_components*num_visis*sizeof(user_precision_t));
  uvw_settings->vs = malloc(num_components*num_visis*sizeof(user_precision_t));
  uvw_settings->ws = malloc(num_components*num_visis*sizeof(user_precision_t));

  uvw_settings->u_metres = malloc(num_visis*sizeof(user_precision_t));
  uvw_settings->v_metres = malloc(num_visis*sizeof(user_precision_t));
  uvw_settings->w_metres = malloc(num_visis*sizeof(user_precision_t));

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

/* Loop over the results stored in uvw_settings and compare to the
expected results to within some tolerance TOL*/

void check_results(user_precision_t* u_metres_expec, user_precision_t* v_metres_expec,
                   user_precision_t* w_metres_expec, user_precision_t* us_expec,
                   user_precision_t* vs_expec, user_precision_t* ws_expec,
                   uvw_settings_t *uvw_settings, int num_visis){

    for (int visi = 0; visi < num_visis; visi++) {

        TEST_ASSERT_DOUBLE_WITHIN(TOL, u_metres_expec[visi],
                                 uvw_settings->u_metres[visi]);
        TEST_ASSERT_DOUBLE_WITHIN(TOL, v_metres_expec[visi],
                                 uvw_settings->v_metres[visi]);
        TEST_ASSERT_DOUBLE_WITHIN(TOL, w_metres_expec[visi],
                                 uvw_settings->w_metres[visi]);

        // printf("%.16f %.16f\n",us_expec[visi], uvw_settings->us[visi]);

        TEST_ASSERT_DOUBLE_WITHIN(TOL, us_expec[visi],
                                 uvw_settings->us[visi]);
        TEST_ASSERT_DOUBLE_WITHIN(TOL, vs_expec[visi],
                                 uvw_settings->vs[visi]);
        TEST_ASSERT_DOUBLE_WITHIN(TOL, ws_expec[visi],
                                 uvw_settings->ws[visi]);
    }
}


woden_settings_t * make_woden_settings(int num_cross, int num_baselines,
                                      int num_times, int num_freqs,
                                      double ra0, double dec0){

  woden_settings_t *woden_settings;
  woden_settings = malloc(sizeof(woden_settings_t));
  woden_settings->ra0 = ra0;
  woden_settings->dec0 = dec0;
  woden_settings->sdec0 = sin(dec0);
  woden_settings->cdec0 = cos(dec0);
  woden_settings->num_cross = num_cross;
  woden_settings->num_baselines = num_baselines;
  woden_settings->num_time_steps = num_times;
  woden_settings->num_freqs = num_freqs;

  return woden_settings;
}

/*
Checking the function fundamental_coords.cu::kern_calc_uvw
Checks that the wavelength scaling of u,v,w is happening correctly. Set HA=0
to make checking easier
*/
void test_calc_uvw_ScalesByWavelength(int do_gpu){

    //Setup some observation settings
    double ra0 = 0.0*DD2R;
    double dec0 = -30.0*DD2R;

    int num_baselines = 5;
    int num_freqs = 5;
    int num_times = 1;
    int num_visis = num_baselines*num_freqs*num_times;

    double lst_base = 0.0;
    user_precision_t time_res = 8.0;
    double freq_res = 50e+6;
    double base_band_freq = 300e+6;

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

    woden_settings_t *woden_settings = make_woden_settings(num_visis, num_baselines,
                                                          num_times, num_freqs, ra0, dec0);
    
    if (do_gpu == 1) {
      //Run the GPU code via fundamental_coords::test_calc_uvw_gpu
      test_calc_uvw_gpu(uvw_settings->X_diff, uvw_settings->Y_diff, uvw_settings->Z_diff,
           uvw_settings->u_metres, uvw_settings->v_metres, uvw_settings->w_metres,
           uvw_settings->us, uvw_settings->vs, uvw_settings->ws, uvw_settings->wavelengths,
           uvw_settings->cha0s, uvw_settings->sha0s,
           woden_settings);
    } else {
      double sdec0 = sin(dec0);
      double cdec0 = cos(dec0);
      calc_uvw_cpu(uvw_settings->X_diff, uvw_settings->Y_diff, uvw_settings->Z_diff,
           uvw_settings->u_metres, uvw_settings->v_metres, uvw_settings->w_metres,
           uvw_settings->us, uvw_settings->vs, uvw_settings->ws, uvw_settings->wavelengths,
           sdec0, cdec0, uvw_settings->cha0s, uvw_settings->sha0s,
           num_visis, num_baselines, num_times, num_freqs);
    }

    //Create expected values
    user_precision_t *u_metres_expec = malloc(num_visis*sizeof(user_precision_t));
    user_precision_t *v_metres_expec = malloc(num_visis*sizeof(user_precision_t));
    user_precision_t *w_metres_expec = malloc(num_visis*sizeof(user_precision_t));
    user_precision_t *us_expec = malloc(num_visis*sizeof(user_precision_t));
    user_precision_t *vs_expec = malloc(num_visis*sizeof(user_precision_t));
    user_precision_t *ws_expec = malloc(num_visis*sizeof(user_precision_t));

    double cdec0 = cos(dec0);
    double sdec0 = sin(dec0);
    //Special case of expected u,v,w when phase centre is at HA = 0.0
    int index = 0;
    for ( int time_step = 0; time_step < num_times; time_step++ ) {
      for (int freq_step = 0; freq_step < num_freqs; freq_step++) {
        for (int baseline = 0; baseline < num_baselines; baseline++) {

          int time_off = time_step*num_baselines;

          u_metres_expec[index] = uvw_settings->Y_diff[time_off + baseline];
          v_metres_expec[index] = -sdec0*uvw_settings->X_diff[time_off + baseline] + cdec0*uvw_settings->Y_diff[time_off + baseline];
          w_metres_expec[index] = cdec0*uvw_settings->X_diff[time_off + baseline] + sdec0*uvw_settings->Y_diff[time_off + baseline];

          us_expec[index] = u_metres_expec[index] / uvw_settings->wavelengths[index];
          vs_expec[index] = v_metres_expec[index] / uvw_settings->wavelengths[index];
          ws_expec[index] = w_metres_expec[index] / uvw_settings->wavelengths[index];

          index += 1;

      }//baseline loop
    }//freq loop
  }//time loop

  check_results(u_metres_expec, v_metres_expec, w_metres_expec,
                us_expec, vs_expec, ws_expec, uvw_settings, num_visis);

  free_uvw_settings(uvw_settings);

}


/*
Checking the function fundamental_coords.cu::kern_calc_uvw
Checks that u,v,w coords change with time as expected
Make checking easier by setting dec phase centre dec0=0.0
*/
void test_calc_uvw_RotateWithTime(int do_gpu){

    //Setup some observation settings
    double ra0 = 0.0*DD2R;
    double dec0 = 0.0*DD2R;

    int num_baselines = 5;
    int num_freqs = 1;
    int num_times = 5;
    int num_visis = num_baselines*num_freqs*num_times;

    double lst_base = 0.0;
    user_precision_t time_res = 8.0;
    double freq_res = 50e+6;
    double base_band_freq = 300e+6;

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

    woden_settings_t *woden_settings = make_woden_settings(num_visis, num_baselines,
                                                          num_times, num_freqs, ra0, dec0);
    
    if (do_gpu == 1) {
      //Run the GPU code via fundamental_coords::test_calc_uvw_gpu
      test_calc_uvw_gpu(uvw_settings->X_diff, uvw_settings->Y_diff, uvw_settings->Z_diff,
           uvw_settings->u_metres, uvw_settings->v_metres, uvw_settings->w_metres,
           uvw_settings->us, uvw_settings->vs, uvw_settings->ws, uvw_settings->wavelengths,
           uvw_settings->cha0s, uvw_settings->sha0s,
           woden_settings);
    } else {
      double sdec0 = sin(dec0);
      double cdec0 = cos(dec0);
      calc_uvw_cpu(uvw_settings->X_diff, uvw_settings->Y_diff, uvw_settings->Z_diff,
           uvw_settings->u_metres, uvw_settings->v_metres, uvw_settings->w_metres,
           uvw_settings->us, uvw_settings->vs, uvw_settings->ws, uvw_settings->wavelengths,
           sdec0, cdec0, uvw_settings->cha0s, uvw_settings->sha0s,
           num_visis, num_baselines, num_times, num_freqs);
    }

    //Create expected values
    user_precision_t *u_metres_expec = malloc(num_visis*sizeof(user_precision_t));
    user_precision_t *v_metres_expec = malloc(num_visis*sizeof(user_precision_t));
    user_precision_t *w_metres_expec = malloc(num_visis*sizeof(user_precision_t));
    user_precision_t *us_expec = malloc(num_visis*sizeof(user_precision_t));
    user_precision_t *vs_expec = malloc(num_visis*sizeof(user_precision_t));
    user_precision_t *ws_expec = malloc(num_visis*sizeof(user_precision_t));

    //Special case of expected u,v,w when phase centre is at Dec = 0.0
    int index = 0;
    for ( int time_step = 0; time_step < num_times; time_step++ ) {
      for (int freq_step = 0; freq_step < num_freqs; freq_step++) {
        for (int baseline = 0; baseline < num_baselines; baseline++) {

          int time_off = time_step*num_baselines;

          u_metres_expec[index] = uvw_settings->sha0s[index]*uvw_settings->X_diff[time_off + baseline] + uvw_settings->cha0s[index]*uvw_settings->Y_diff[time_off + baseline];
          v_metres_expec[index] = uvw_settings->Z_diff[time_off + baseline];
          w_metres_expec[index] = uvw_settings->cha0s[index]*uvw_settings->X_diff[time_off + baseline] - uvw_settings->sha0s[index]*uvw_settings->Y_diff[time_off + baseline];

          us_expec[index] = u_metres_expec[index] / uvw_settings->wavelengths[index];
          vs_expec[index] = v_metres_expec[index] / uvw_settings->wavelengths[index];
          ws_expec[index] = w_metres_expec[index] / uvw_settings->wavelengths[index];

          index += 1;

      }//baseline loop
    }//freq loop
  }//time loop

  check_results(u_metres_expec, v_metres_expec, w_metres_expec,
                us_expec, vs_expec, ws_expec, uvw_settings, num_visis);

  free_uvw_settings(uvw_settings);

}

/*
Checking the function fundamental_coords.cu::kern_calc_uvw_shapelet
Checks that u,v,w coords change with time as expected
Make checking easier by setting dec phase centre dec0=0.0
Also checks that results are scaled by wavelength correctly
*/
void test_calc_uvw_shapelet_RotateWithTime(int do_gpu){

    //Setup some observation settings
    double ra0 = 0.0*DD2R;

    int num_baselines = 5;
    int num_freqs = 2;
    int num_times = 5;
    int num_visis = num_baselines*num_freqs*num_times;

    double lst_base = 0.0;
    user_precision_t time_res = 8.0;
    double freq_res = 50e+6;
    double base_band_freq = 300e+6;

    //Check values for 5 different shapelet components, with different RA values
    int num_components = 5;

    double *ras = malloc(num_components*sizeof(double));
    double *decs = malloc(num_components*sizeof(double));

    for (int i = 0; i < num_components; i++) {
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

    if (do_gpu == 1) {

      //Run the CUDA code via fundamental_coords::test_calc_uvw_gpu_shapelet
      test_calc_uv_shapelet_gpu(uvw_settings->X_diff,
            uvw_settings->Y_diff, uvw_settings->Z_diff,
            uvw_settings->us, uvw_settings->vs,
            uvw_settings->lsts, ras, decs,
            num_baselines, num_times, num_components);
    } else {
      calc_uv_shapelet_cpu(uvw_settings->X_diff, uvw_settings->Y_diff, uvw_settings->Z_diff,
            uvw_settings->us, uvw_settings->vs,
            uvw_settings->lsts, ras, decs,
            num_baselines, num_times, num_components);
    }

    //Create expected values
    user_precision_t *us_expec = malloc(num_components*num_baselines*num_times*sizeof(user_precision_t));
    user_precision_t *vs_expec = malloc(num_components*num_baselines*num_times*sizeof(user_precision_t));
    //Variables to use when making expected values
    user_precision_t u_metre, v_metre;
    double ha, sha, cha;

    int uvw_index = 0;
    for ( int comp_step = 0; comp_step < num_components; comp_step++ ) {
      for ( int time_step = 0; time_step < num_times; time_step++ ) {
        //For kern_calc_uvw_shapelet, calculating u,v,w centred on ra, dec of
        //each shapelet component. Calcualte ha using the shapelet ra
        ha = uvw_settings->lsts[time_step] - ras[comp_step];
        sha = sin(ha);
        cha = cos(ha);

        for (int baseline = 0; baseline < num_baselines; baseline++) {
          //Special case of expected u,v,w when u,v,w coord centre is at Dec = 0.0

          int time_off = time_step*num_baselines;

          u_metre = sha*uvw_settings->X_diff[time_off + baseline] + cha*uvw_settings->Y_diff[time_off + baseline];
          v_metre = uvw_settings->Z_diff[time_off + baseline];
          int stripe = comp_step*num_baselines*num_times + time_step*num_baselines;

          us_expec[stripe + baseline] = u_metre;
          vs_expec[stripe + baseline] = v_metre;

          uvw_index += 1;

      }//baseline loop
    }//time loop
  }//component loop

  for (int ind = 0; ind < num_components*num_baselines*num_times; ind++) {

    TEST_ASSERT_DOUBLE_WITHIN(TOL, us_expec[ind], uvw_settings->us[ind]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, vs_expec[ind], uvw_settings->vs[ind]);
  }

  free_uvw_settings(uvw_settings);

}