#pragma once
#include "woden_precision_defs.h"
#include <unity.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"

#define TOL 1e-11

typedef struct _uvw_settings_t {
  double *X_diff;
  double *Y_diff;
  double *Z_diff;
  user_precision_t *wavelengths;
  double *cha0s;
  double *sha0s;
  user_precision_t *us;
  user_precision_t *vs;
  user_precision_t *ws;
  user_precision_t *u_metres;
  user_precision_t *v_metres;
  user_precision_t *w_metres;
  double *lsts;  /*!< lst value for every visibility */

} uvw_settings_t;

/*
Given the inputs, create simulation settings that woden.c would create
*/
void setup_uvw_params(int num_times, int num_baselines, int num_freqs,
                      user_precision_t ra0,
                      user_precision_t lst_base, user_precision_t time_res,
                      user_precision_t freq_res, user_precision_t base_band_freq,
                      uvw_settings_t *uvw_settings);

/*
Do a bunch of mallocing
*/
uvw_settings_t * setup_uvw_settings(int num_baselines, int num_visis,
                                    int num_times, int num_components);

/*
Do a bunch of freeing
*/
void free_uvw_settings(uvw_settings_t * uvw_settings);

/* Loop over the results stored in uvw_settings and compare to the
expected results to within some tolerance TOL*/
void check_results(user_precision_t* u_metres_expec, user_precision_t* v_metres_expec,
                   user_precision_t* w_metres_expec, user_precision_t* us_expec,
                   user_precision_t* vs_expec, user_precision_t* ws_expec,
                   uvw_settings_t *uvw_settings, int num_visis);

/*
Checking the function fundamental_coords.cu::kern_calc_uvw
Checks that the wavelength scaling of u,v,w is happening correctly. Set HA=0
to make checking easier. Either runs test through GPU or CPU code
do_GPU = 1 for GPU, 0 for CPU
*/
void test_calc_uvw_ScalesByWavelength(int do_gpu);

/*
Checking the function fundamental_coords.cu::kern_calc_uvw
Checks that u,v,w coords change with time as expected
Make checking easier by setting dec phase centre dec0=0.0
*/
void test_calc_uvw_RotateWithTime(int do_gpu);

/*
Checking the function fundamental_coords.cu::kern_calc_uvw_shapelet
Checks that u,v,w coords change with time as expected
Make checking easier by setting dec phase centre dec0=0.0
Also checks that results are scaled by wavelength correctly
*/
void test_calc_uvw_shapelet_RotateWithTime(int do_gpu);