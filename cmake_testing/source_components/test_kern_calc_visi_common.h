#include <math.h>
#include <unity.h>
#include <stdlib.h>
#include <complex.h>

#include "constants.h"
#include "woden_struct_defs.h"
#include "shapelet_basis.h"

void sincosf(float x, float *sin, float *cos);

#define UNITY_INCLUDE_FLOAT

/*
Something to hold all the input argument arrays for testing
*/
typedef struct _args_for_testing_t {

  int num_baselines;
  int num_times;
  int num_freqs;
  int num_visis;
  int num_beam_values;

  float *ls;
  float *ms;
  float *ns;
  float _Complex *primay_beam_J00;
  float _Complex *primay_beam_J01;
  float _Complex *primay_beam_J10;
  float _Complex *primay_beam_J11;

  float *flux_I;
  float *flux_Q;
  float *flux_U;
  float *flux_V;
  float *SIs;
  float *component_freqs;

  float *us;
  float *vs;
  float *ws;
  float *allsteps_wavelengths;

  //GAUSS/SHAPELET STUFF
  float *pas;
  float *majors;
  float *minors;

  float *sum_visi_XX_real;
  float *sum_visi_XX_imag;
  float *sum_visi_XY_real;
  float *sum_visi_XY_imag;
  float *sum_visi_YX_real;
  float *sum_visi_YX_imag;
  float *sum_visi_YY_real;
  float *sum_visi_YY_imag;

  //SHAPELET stuff
  float *sbf;
  float *u_shapes;
  float *v_shapes;
  float *w_shapes;
  float *shape_n1s;
  float *shape_n2s;
  float *shape_coeffs;
  float *shape_param_indexes;
  int num_coeffs;

} args_for_testing_t;

void malloc_args_for_testing(args_for_testing_t *args_ft,
                            int num_baselines,  int num_times,
                            int num_freqs, int num_components,
                            int num_coeffs,
                            int component_type);

void free_args_for_testing(args_for_testing_t *args_ft, int component_type);

/*
Setup some l,m,n coords. HARD CODED TO BE 5 by 5 GRID spanning -0.5 to 0.5
*/
void create_lmn(args_for_testing_t *args_ft);

/*
Basic implementation of the measurement equation to get expected visibilities
Loops over components, gets expected flux and beam gain and sum
*/
void get_expected(int visi, int num_components, int num_baselines,
                  int num_freqs, int beamtype,
                  args_for_testing_t *args_ft,
                  int component_type,
                  float * expec_re, float * expec_im);

//Take input parameters and test whether GPU outputs match expectations
void test_visi_outputs(int num_visis, int num_components,
                       int num_baselines, int num_freqs,
                       float frac_tol,
                       int beamtype,  args_for_testing_t *args_ft,
                       int component_type);
