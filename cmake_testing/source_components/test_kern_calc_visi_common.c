// #include <math.h>
// #include <unity.h>
// #include <stdlib.h>
// #include <complex.h>
//
// #include "constants.h"
// #include "woden_struct_defs.h"
#include "test_kern_calc_visi_common.h"

void sincosf(float x, float *sin, float *cos);

#define UNITY_INCLUDE_FLOAT

void malloc_args_for_testing(args_for_testing_t *args_ft,
                            int num_baselines,  int num_times,
                            int num_freqs, int num_components,
                            int num_coeffs,
                            int component_type) {

  int num_visis = num_baselines*num_times*num_freqs;
  int num_beam_values = num_freqs*num_times*num_components;

  args_ft->num_baselines = num_baselines;
  args_ft->num_times = num_times;
  args_ft->num_freqs = num_freqs;
  args_ft->num_visis = num_visis;
  args_ft->num_beam_values = num_beam_values;

  args_ft->ls = malloc(num_components*sizeof(float));
  args_ft->ms = malloc(num_components*sizeof(float));
  args_ft->ns = malloc(num_components*sizeof(float));

  //Assume we just want Stokes I
  args_ft->flux_I = malloc(num_components*sizeof(float));
  args_ft->flux_Q = calloc(num_components, sizeof(float));
  args_ft->flux_U = calloc(num_components, sizeof(float));
  args_ft->flux_V = calloc(num_components, sizeof(float));
  args_ft->SIs = malloc(num_components*sizeof(float));
  args_ft->component_freqs = malloc(num_components*sizeof(float));

  args_ft->us = malloc(num_visis*sizeof(float));
  args_ft->vs = malloc(num_visis*sizeof(float));
  args_ft->ws = malloc(num_visis*sizeof(float));
  args_ft->allsteps_wavelengths = malloc(num_visis*sizeof(float));

  args_ft->primay_beam_J00 = malloc(num_beam_values*sizeof(float _Complex));
  args_ft->primay_beam_J01 = malloc(num_beam_values*sizeof(float _Complex));
  args_ft->primay_beam_J10 = malloc(num_beam_values*sizeof(float _Complex));
  args_ft->primay_beam_J11 = malloc(num_beam_values*sizeof(float _Complex));

  //GAUSS/SHAPELET STUFF
  args_ft->pas = malloc(num_components*sizeof(float));
  args_ft->majors = malloc(num_components*sizeof(float));
  args_ft->minors = malloc(num_components*sizeof(float));

  //Make sure arrays to hold summed visis are initialised to zero
  args_ft->sum_visi_XX_real = calloc(num_visis, sizeof(float));
  args_ft->sum_visi_XX_imag = calloc(num_visis, sizeof(float));
  args_ft->sum_visi_XY_real = calloc(num_visis, sizeof(float));
  args_ft->sum_visi_XY_imag = calloc(num_visis, sizeof(float));
  args_ft->sum_visi_YX_real = calloc(num_visis, sizeof(float));
  args_ft->sum_visi_YX_imag = calloc(num_visis, sizeof(float));
  args_ft->sum_visi_YY_real = calloc(num_visis, sizeof(float));
  args_ft->sum_visi_YY_imag = calloc(num_visis, sizeof(float));

  if (component_type == SHAPELET) {
    args_ft->num_coeffs = num_coeffs;
    args_ft->sbf = malloc( sbf_N * sbf_L * sizeof(float) );
    args_ft->sbf = create_sbf(args_ft->sbf);
    args_ft->u_shapes = malloc(num_components*num_visis*sizeof(float));
    args_ft->v_shapes = malloc(num_components*num_visis*sizeof(float));
    args_ft->w_shapes = malloc(num_components*num_visis*sizeof(float));
    args_ft->shape_n1s = malloc(num_coeffs*sizeof(float));
    args_ft->shape_n2s = malloc(num_coeffs*sizeof(float));
    args_ft->shape_coeffs = malloc(num_coeffs*sizeof(float));
    args_ft->shape_param_indexes = malloc(num_coeffs*sizeof(float));
  }
}

void free_args_for_testing(args_for_testing_t *args_ft, int component_type) {

  free( args_ft->ls );
  free( args_ft->ms );
  free( args_ft->ns );

  free( args_ft->flux_I );
  free( args_ft->flux_Q );
  free( args_ft->flux_U );
  free( args_ft->flux_V );
  free( args_ft->SIs );
  free( args_ft->component_freqs );

  free( args_ft->us );
  free( args_ft->vs );
  free( args_ft->ws );
  free( args_ft->allsteps_wavelengths );

  free( args_ft->sum_visi_XX_real );
  free( args_ft->sum_visi_XX_imag );
  free( args_ft->sum_visi_XY_real );
  free( args_ft->sum_visi_XY_imag );
  free( args_ft->sum_visi_YX_real );
  free( args_ft->sum_visi_YX_imag );
  free( args_ft->sum_visi_YY_real );
  free( args_ft->sum_visi_YY_imag );

  free( args_ft->pas );
  free( args_ft->majors );
  free( args_ft->minors );

  free( args_ft->primay_beam_J00 );
  free( args_ft->primay_beam_J01 );
  free( args_ft->primay_beam_J10 );
  free( args_ft->primay_beam_J11 );

  if (component_type == SHAPELET) {
    free( args_ft->sbf );
    free( args_ft->u_shapes );
    free( args_ft->v_shapes );
    free( args_ft->w_shapes );
    free( args_ft->shape_n1s );
    free( args_ft->shape_n2s );
    free( args_ft->shape_coeffs );
    free( args_ft->shape_param_indexes );
  }

  free( args_ft );

}

/*
Setup some l,m,n coords. HARD CODED TO BE 5 by 5 GRID spanning -0.5 to 0.5
*/
void create_lmn(args_for_testing_t *args_ft) {
  int count = 0;
  for (int l_ind = 0; l_ind < 5; l_ind++) {
    for (int m_ind = 0; m_ind < 5; m_ind++) {
      args_ft->ls[count] = -0.5 + 0.25*l_ind;
      args_ft->ms[count] = -0.5 + 0.25*m_ind;
      args_ft->ns[count] = sqrt(1 - args_ft->ls[count]*args_ft->ls[count] - args_ft->ms[count]*args_ft->ms[count]);

      count ++;
    }
  }
}

float _Complex visi_env_shape(int comp, int visi, int coeff,
                              args_for_testing_t *args_ft ) {
  float pa = args_ft->pas[comp];
  float sinpa = sin(pa);
  float cospa = cos(pa);

  float u_shape = args_ft->u_shapes[comp*args_ft->num_visis + visi];
  float v_shape = args_ft->v_shapes[comp*args_ft->num_visis + visi];

  float x = (cospa*v_shape + sinpa*u_shape); // major axis
  float y = (-sinpa*v_shape + cospa*u_shape); // minor axis

  //Scales the FWHM to std to match basis functions, and account for the
  //basis functions being stored with beta = 1.0
  //Basis functions have been stored in such a way that x is in the same
  //direction as on sky, but y is opposite, so include negative here
  float const_x = (args_ft->majors[comp]*SQRT_M_PI_2_2_LN_2)/sbf_dx;
  float const_y = -(args_ft->minors[comp]*SQRT_M_PI_2_2_LN_2)/sbf_dx;

  float _Complex Ipow_lookup[] = { 1.0 + I*0.0,
                                   0.0 + I*1.0,
                                  -1.0 + I*0.0,
                                   0.0 + I*-1.0 };

  float xlow, xhigh, ylow, yhigh, u_value, v_value, f_hat, *sbf_n;

  // find the indices in the basis functions for u*beta_u and v*beta_v

  float xpos = x*const_x + sbf_c;
  float ypos = y*const_y + sbf_c;

  int xindex = (int)floor(xpos);
  int yindex = (int)floor(ypos);
  //
  int n1 = (int)args_ft->shape_n1s[coeff];
  int n2 = (int)args_ft->shape_n2s[coeff];

  // if ( n1 < 0 || n2 < 0 || n1 >= sbf_N || n2 >= sbf_N ) continue;

  f_hat = args_ft->shape_coeffs[coeff];
  //
  sbf_n = &args_ft->sbf[n1*sbf_L];
  xlow  = sbf_n[xindex];
  xhigh = sbf_n[xindex+1];
  u_value = xlow + (xhigh-xlow)*(xpos-xindex);

  sbf_n = &args_ft->sbf[n2*sbf_L];
  ylow  = sbf_n[yindex];
  yhigh = sbf_n[yindex+1];
  v_value = ylow + (yhigh-ylow)*(ypos-yindex);

  // accumulate the intensity model for baseline pair (u,v)
  float _Complex visi_env = 0.0 + I*0.0;
  visi_env = visi_env + Ipow_lookup[(n1+n2) % 4] * f_hat * u_value*v_value;

  return visi_env;
}

float _Complex visi_env_gauss(int comp, int visi,
                              args_for_testing_t *args_ft ) {
  float pa = args_ft->pas[comp];
  float u = args_ft->us[visi];
  float v = args_ft->vs[visi];
  float sinpa = sin(pa);
  float cospa = cos(pa);

  float x =  cospa*v + sinpa*u; // major axis
  float y = -sinpa*v + cospa*u; // minor axis
  float invsig_x = args_ft->majors[comp];
  float invsig_y = args_ft->minors[comp];

  float _Complex visi_env = exp( -0.5 * ( x*x*invsig_x*invsig_x*M_PI_2_2_LN_2 + y*y*invsig_y*invsig_y*M_PI_2_2_LN_2 ) ) + I*0.0;

  return visi_env;
}

/*
Basic implementation of the measurement equation to get expected visibilities
Loops over components, gets expected flux and beam gain and sum
*/
void get_expected(int visi, int num_components, int num_baselines,
                  int num_freqs, int beamtype,
                  args_for_testing_t *args_ft,
                  int component_type,
                  float * expec_re, float * expec_im) {
  float expec_re_inc, expec_im_inc, temp;
  float flux_ratio, visi_freq, flux_extrap, xx_gain;
  int time_ind, freq_ind, beam_ind;
  * expec_re = 0.0;
  * expec_im = 0.0;
  for (size_t comp = 0; comp < num_components; comp++) {
    // printf("%.5f %.5f %.5f\n", ls[comp], ms[comp], ns[comp] );
    temp = 2*M_PI*( args_ft->us[visi]*args_ft->ls[comp] + args_ft->vs[visi]*args_ft->ms[comp] + args_ft->ws[visi]*(args_ft->ns[comp]-1) );
    sincosf(temp, &(expec_im_inc), &(expec_re_inc));

    visi_freq = VELC / args_ft->allsteps_wavelengths[visi];
    flux_ratio = powf(visi_freq / args_ft->component_freqs[comp], args_ft->SIs[comp]);
    flux_extrap = flux_ratio*args_ft->flux_I[comp];

    //Do some indexing so we can call up correct instrumental gain
    time_ind = (int)floorf( (float)visi / ((float)num_baselines * (float)num_freqs));
    freq_ind = (int)floorf( ((float)visi - ((float)time_ind*(float)num_baselines * (float)num_freqs)) / (float)num_baselines);
    beam_ind = num_freqs*time_ind*num_components + (num_components*freq_ind) + comp;

    if (beamtype == NO_BEAM) {
      xx_gain = 1.0;
    } else {
      //In tests below have set gains purely real, and all pols equal, so
      //only need to call one pol and square it
      xx_gain = creal(args_ft->primay_beam_J00[beam_ind]);
      xx_gain = xx_gain*xx_gain;
    }

    expec_re_inc = expec_re_inc*flux_extrap*xx_gain;
    expec_im_inc = expec_im_inc*flux_extrap*xx_gain;

    float _Complex visi_env = 0.0 + I*0.0;

    if (component_type == POINT) {
      visi_env = 1.0 + I*0.0;
    }
    else if (component_type == GAUSSIAN) {
      visi_env = visi_env_gauss(comp, visi, args_ft);
    }
    else if (component_type == SHAPELET) {
      for (int coeff = 0; coeff < args_ft->num_coeffs; coeff++) {
        if (args_ft->shape_param_indexes[coeff] == comp) {
          // visi_env = visi_env_gauss(comp, visi, args_ft);
          visi_env += visi_env_shape(comp, visi, coeff, args_ft);

        }
      }
    }
    * expec_re += expec_re_inc*creal(visi_env) - expec_im_inc*cimag(visi_env);
    * expec_im += expec_re_inc*cimag(visi_env) + expec_im_inc*creal(visi_env);
  }
}

//Take input parameters and test whether GPU outputs match expectations
void test_visi_outputs(int num_visis, int num_components,
                       int num_baselines, int num_freqs,
                       float frac_tol,
                       int beamtype,  args_for_testing_t *args_ft,
                       int component_type) {
  float expec_re, expec_im;
  for (int visi = 0; visi < num_visis; visi++) {

      get_expected(visi, num_components, num_baselines, num_freqs,
                  beamtype, args_ft, component_type,
                  &expec_re, &expec_im);

    // if (visi == 0) {
    //   printf("%d %.1f %.1f %.1f\n",
    //         visi, expec_re, args_ft->sum_visi_XX_real[visi],
    //         args_ft->sum_visi_YY_real[visi]);
    // }

    // printf("%d %.1f %.1f %.1f\n",
    //       visi, expec_re, args_ft->sum_visi_XX_real[visi],
    //       args_ft->sum_visi_YY_real[visi]);

    // printf("%d %.5f %.5f %.5f %.5f %.5f %.5f\n",
    //       visi, expec_re, expec_im,
    //       args_ft->sum_visi_XX_real[visi], args_ft->sum_visi_XX_imag[visi],
    //       args_ft->sum_visi_YY_real[visi], args_ft->sum_visi_YY_imag[visi]);

    // printf("%.1f %.1f %.1f %.1f %.1f %.1f %.1f %.1f\n",
    //       args_ft->sum_visi_XX_real[visi], args_ft->sum_visi_XX_imag[visi],
    //       args_ft->sum_visi_XY_real[visi], args_ft->sum_visi_XY_imag[visi],
    //       args_ft->sum_visi_YX_real[visi], args_ft->sum_visi_YX_imag[visi],
    //       args_ft->sum_visi_YY_real[visi], args_ft->sum_visi_YY_imag[visi]);

    if (beamtype == FEE_BEAM) {
      //FEE beam has cross pols which double everything when using just Stokes I
      //Also means cross-pols are non-zero
      TEST_ASSERT_FLOAT_WITHIN(frac_tol*2, 2.0, args_ft->sum_visi_XX_real[visi] / expec_re);
      TEST_ASSERT_FLOAT_WITHIN(frac_tol*2, 2.0, args_ft->sum_visi_XX_imag[visi] / expec_im);
      TEST_ASSERT_FLOAT_WITHIN(frac_tol*2, 2.0, args_ft->sum_visi_XY_real[visi] / expec_re);
      TEST_ASSERT_FLOAT_WITHIN(frac_tol*2, 2.0, args_ft->sum_visi_XY_imag[visi] / expec_im);
      TEST_ASSERT_FLOAT_WITHIN(frac_tol*2, 2.0, args_ft->sum_visi_YX_real[visi] / expec_re);
      TEST_ASSERT_FLOAT_WITHIN(frac_tol*2, 2.0, args_ft->sum_visi_YX_imag[visi] / expec_im);
      TEST_ASSERT_FLOAT_WITHIN(frac_tol*2, 2.0, args_ft->sum_visi_YY_real[visi] / expec_re);
      TEST_ASSERT_FLOAT_WITHIN(frac_tol*2, 2.0, args_ft->sum_visi_YY_imag[visi] / expec_im);
    }
    else {
      //
      TEST_ASSERT_FLOAT_WITHIN(frac_tol, 1.0, args_ft->sum_visi_XX_real[visi] / expec_re);
      TEST_ASSERT_FLOAT_WITHIN(frac_tol, 1.0, args_ft->sum_visi_XX_imag[visi] / expec_im);
      TEST_ASSERT_EQUAL_FLOAT(0.0, args_ft->sum_visi_XY_real[visi]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, args_ft->sum_visi_XY_imag[visi]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, args_ft->sum_visi_YX_real[visi]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, args_ft->sum_visi_YX_imag[visi]);
      TEST_ASSERT_FLOAT_WITHIN(frac_tol, 1.0, args_ft->sum_visi_YY_real[visi] / expec_re);
      TEST_ASSERT_FLOAT_WITHIN(frac_tol, 1.0, args_ft->sum_visi_YY_imag[visi] / expec_im);

    }
  }
}
