#pragma once
#include <math.h>
#include <unity.h>
#include <stdlib.h>
#include <complex.h>

#include "constants.h"
#include "woden_struct_defs.h"
#include "woden_precision_defs.h"
#include "shapelet_basis.h"
#include "common_testing_functions.h"

//external CUDA code used for testng

extern void test_kern_calc_visi_all(int n_powers, int n_curves, int n_lists,
          int num_baselines, int num_shape_coeffs,
          int num_freqs, int num_visis, int num_times,
          e_beamtype beamtype, e_component_type comptype,
          components_t components, double *extrap_freqs,
          user_precision_t *us, user_precision_t *vs, user_precision_t *ws,
          user_precision_t *u_shapes, user_precision_t *v_shapes,
          user_precision_t *sum_visi_XX_real, user_precision_t *sum_visi_XX_imag,
          user_precision_t *sum_visi_XY_real, user_precision_t *sum_visi_XY_imag,
          user_precision_t *sum_visi_YX_real, user_precision_t *sum_visi_YX_imag,
          user_precision_t *sum_visi_YY_real, user_precision_t *sum_visi_YY_imag,
          user_precision_t *allsteps_wavelengths, user_precision_t *sbf,
          user_precision_complex_t *gxs, user_precision_complex_t *Dxs,
          user_precision_complex_t *Dys, user_precision_complex_t *gys);

/*
Something to hold all the input argument arrays for testing
*/
typedef struct _args_for_testing_t {

  int num_baselines;
  int num_times;
  int num_freqs;
  int num_visis;
  int num_beam_values;

  double *ls;
  double *ms;
  double *ns;
  user_precision_complex_t *primay_beam_J00;
  user_precision_complex_t *primay_beam_J01;
  user_precision_complex_t *primay_beam_J10;
  user_precision_complex_t *primay_beam_J11;

  user_precision_t *flux_I;
  user_precision_t *flux_Q;
  user_precision_t *flux_U;
  user_precision_t *flux_V;
  user_precision_t *SIs;
  double *component_freqs;

  user_precision_t *us;
  user_precision_t *vs;
  user_precision_t *ws;
  user_precision_t *allsteps_wavelengths;
  double *extrap_freqs;

  //GAUSS/SHAPELET STUFF
  user_precision_t *pas;
  user_precision_t *majors;
  user_precision_t *minors;

  user_precision_t *sum_visi_XX_real;
  user_precision_t *sum_visi_XX_imag;
  user_precision_t *sum_visi_XY_real;
  user_precision_t *sum_visi_XY_imag;
  user_precision_t *sum_visi_YX_real;
  user_precision_t *sum_visi_YX_imag;
  user_precision_t *sum_visi_YY_real;
  user_precision_t *sum_visi_YY_imag;

  //SHAPELET stuff
  user_precision_t *sbf;
  user_precision_t *u_shapes;
  user_precision_t *v_shapes;
  user_precision_t *w_shapes;
  user_precision_t *shape_n1s;
  user_precision_t *shape_n2s;
  user_precision_t *shape_coeffs;
  user_precision_t *shape_param_indexes;
  int num_coeffs;

} args_for_testing_t;

void malloc_args_for_testing(args_for_testing_t *args_ft,
                             components_t *components,
                             int num_baselines,  int num_times,
                             int num_freqs, int num_components,
                             int n_powers, int n_curves, int n_lists,
                             int num_coeffs,
                             e_component_type component_type);

void free_args_for_testing(args_for_testing_t *args_ft,
                           components_t components,
                           e_component_type component_type);

/*
Setup some l,m,n coords. HARD CODED TO BE 5 by 5 GRID spanning -0.5 to 0.5
*/
void create_lmn(components_t components);

/*
Basic implementation of the measurement equation to get expected visibilities
Loops over components, gets expected flux and beam gain and sum
*/
void get_expected(int visi, int num_powers, int num_curves, int num_lists,
                  int num_baselines, int num_freqs, double *extrap_freqs,
                  int beamtype,
                  args_for_testing_t *args_ft,
                  components_t components,
                  e_component_type component_type,
                  double * expec_re, double * expec_im);

//Take input parameters and test whether GPU outputs match expectations
void test_visi_outputs(int num_visis, int num_powers, int num_curves, int num_lists,
                       int num_baselines, int num_freqs, double *extrap_freqs,
                       e_beamtype beamtype,  args_for_testing_t *args_ft,
                       components_t components,
                       e_component_type component_type);

void test_kern_calc_visi_Varylmn(e_beamtype beamtype, e_component_type comptype);

void test_kern_calc_visi_VarylmnVaryBeam(e_beamtype beamtype,
                                         e_component_type comptype);

void test_kern_calc_visi_VarylmnVaryFlux(e_beamtype beamtype,
                                         e_component_type comptype);

void test_kern_calc_visi_VarylmnVaryPAMajMin(e_beamtype beamtype,
                                             e_component_type comptype);

void setup_uvw_and_freqs(args_for_testing_t *args_ft, int num_times,
                         int num_freqs, int num_baselines);
