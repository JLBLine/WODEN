#include <math.h>
#include <unity.h>
#include <stdlib.h>
#include <complex.h>

#include "constants.h"
#include "woden_struct_defs.h"
#include "woden_precision_defs.h"
#include "shapelet_basis.h"

// void sincos(float x, float *sin, float *cos);
// void sincos(double x, double *sin, double *cos);

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
  user_precision_t *component_freqs;

  user_precision_t *us;
  user_precision_t *vs;
  user_precision_t *ws;
  user_precision_t *allsteps_wavelengths;

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
                  user_precision_t * expec_re, user_precision_t * expec_im);

//Take input parameters and test whether GPU outputs match expectations
void test_visi_outputs(int num_visis, int num_components,
                       int num_baselines, int num_freqs,
                       user_precision_t frac_tol,
                       int beamtype,  args_for_testing_t *args_ft,
                       int component_type);
