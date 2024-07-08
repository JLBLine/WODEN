#include <unity.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "constants.h"
#include "woden_struct_defs.h"
#include "MWA_analy_expected_nside101.h"
#include "azza_radec_nside101.h"

#define FREQ 150000000

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

//External CUDA code we're linking in

extern void test_RTS_calculate_MWA_analytic_beam(int num_components,
     int num_time_steps, int num_freqs,
     user_precision_t *azs, user_precision_t *zas, int *delays,
     double latitude, int norm,
     double *beam_has, double *beam_decs, double *freqs,
     user_precision_complex_t *gxs, user_precision_complex_t *Dxs,
     user_precision_complex_t *Dys, user_precision_complex_t *gys);

#ifdef DOUBLE_PRECISION
  double TOL = 1e-8;
#else
  double TOL = 1e-6;
#endif

void test_MWA_analytic_beam_nside101() {


  int nside = 101;

  int num_times = 2;
  int num_freqs = 2;

  double freqs[] = {150e+6, 200e+6};

  int num_components = nside*nside;
  int num_azza = num_components*num_times;
  int num_beam_values = num_azza*num_freqs;

  double *beam_has = malloc(num_azza*sizeof(double));
  double *beam_decs = malloc(num_azza*sizeof(double));
  user_precision_t *all_azs = malloc(num_azza*sizeof(user_precision_t));
  user_precision_t *all_zas = malloc(num_azza*sizeof(user_precision_t));

  for (int coord = 0; coord < num_components; coord++) {
    for (int time = 0; time < num_times; time++) {

      beam_has[coord*num_times + time] = nside101_has[coord];
      beam_decs[coord*num_times + time] = nside101_decs[coord];

      all_azs[coord*num_times + time] = nside101_azs[coord];
      all_zas[coord*num_times + time] = nside101_zas[coord];
    }
  }

  user_precision_complex_t *gxs = malloc(num_beam_values*sizeof(user_precision_complex_t));
  user_precision_complex_t *Dxs = malloc(num_beam_values*sizeof(user_precision_complex_t));
  user_precision_complex_t *Dys = malloc(num_beam_values*sizeof(user_precision_complex_t));
  user_precision_complex_t *gys = malloc(num_beam_values*sizeof(user_precision_complex_t));

  int norm = 1;
  // int delays[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    // 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  int delays[] = {0, 2, 4, 8, 2, 4, 8, 12, 4, 8, 12, 16, 8, 12, 16, 20};
  //Run the CUDA code
  test_RTS_calculate_MWA_analytic_beam(num_components,
       num_times, num_freqs,
       all_azs, all_zas, delays,
       MWA_LAT_RAD, norm,
       beam_has, beam_decs, freqs,
       gxs, Dxs, Dys, gys);

  FILE *output_text;

  char filename[100];
  sprintf(filename, "MWA_analy_gains_azza%02d.txt", num_components);

  output_text = fopen(filename,"w");

  int count = 0;

  for (int time = 0; time < num_times; time ++) {
    for (int freq = 0; freq < num_freqs; freq ++) {
      for (int comp = 0; comp < num_components; comp ++) {

        int beam_ind = num_freqs*time*num_components + num_components*freq + comp;
        int coord_ind = comp*num_times + time;

        fprintf(output_text,"%.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.1f\n",
         all_azs[coord_ind], all_zas[coord_ind],
         creal(gxs[beam_ind]), cimag(gxs[beam_ind]),
         creal(Dxs[beam_ind]), cimag(Dxs[beam_ind]),
         creal(Dys[beam_ind]), cimag(Dys[beam_ind]),
         creal(gys[beam_ind]), cimag(gys[beam_ind]),
         freqs[freq] );

         count ++;
      }
    }
  }

  for (int ind = 0; ind < num_beam_values; ind++) {

    TEST_ASSERT_DOUBLE_WITHIN(TOL, creal(gxs[ind]), gx_re_expec[ind]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, creal(Dxs[ind]), Dx_re_expec[ind]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, creal(Dys[ind]), Dy_re_expec[ind]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, creal(gys[ind]), gy_re_expec[ind]);

    TEST_ASSERT_DOUBLE_WITHIN(TOL, cimag(gxs[ind]), 0.0);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, cimag(Dxs[ind]), 0.0);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, cimag(Dys[ind]), 0.0);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, cimag(gys[ind]), 0.0);
  }

  fflush(output_text);
  fclose(output_text);

  free(all_azs);
  free(all_zas);
  free(beam_has);
  free(beam_decs);
  free(gxs);
  free(Dxs);
  free(Dys);
  free(gys);
}

//Same test but for a single az,za coord - good for diagnosing bugs
void test_single_azza() {


  int num_times = 1;
  int num_freqs = 1;

  double freqs[] = {170e+6};

  int num_components = 1;
  int num_azza = num_components*num_times;
  int num_beam_values = num_azza*num_freqs;

  double *beam_has = malloc(num_azza*sizeof(double));
  double *beam_decs = malloc(num_azza*sizeof(double));
  user_precision_t *all_azs = malloc(num_azza*sizeof(user_precision_t));
  user_precision_t *all_zas = malloc(num_azza*sizeof(user_precision_t));

  beam_has[0] = 0.24457606;
  beam_decs[0] = -0.69813170;
  all_azs[0] = 3.7997899;
  all_zas[0] = 0.3080986;

  printf("azza %.12f %.12f \n",all_azs[0],all_zas[0] );

  user_precision_complex_t *gxs = malloc(num_beam_values*sizeof(user_precision_complex_t));
  user_precision_complex_t *Dxs = malloc(num_beam_values*sizeof(user_precision_complex_t));
  user_precision_complex_t *Dys = malloc(num_beam_values*sizeof(user_precision_complex_t));
  user_precision_complex_t *gys = malloc(num_beam_values*sizeof(user_precision_complex_t));

  int norm = 1;
  // int delays[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    // 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  int delays[] = {0, 2, 4, 8, 2, 4, 8, 12, 4, 8, 12, 16, 8, 12, 16, 20};
  //Run the CUDA code
  test_RTS_calculate_MWA_analytic_beam(num_components,
       num_times, num_freqs,
       all_azs, all_zas, delays,
       MWA_LAT_RAD, norm,
       beam_has, beam_decs, freqs,
       gxs, Dxs, Dys, gys);

  printf("gx gy %.8f %.8f\n", creal(gxs[0]), creal(gys[0]));

  free(all_azs);
  free(all_zas);
  free(beam_has);
  free(beam_decs);
  free(gxs);
  free(Dxs);
  free(Dys);
  free(gys);
}

//Run the test with unity
int main(void)
{
    UNITY_BEGIN();
    RUN_TEST(test_MWA_analytic_beam_nside101);
    // RUN_TEST(test_single_azza);

    return UNITY_END();
}
