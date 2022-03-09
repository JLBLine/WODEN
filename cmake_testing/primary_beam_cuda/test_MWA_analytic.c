#include <unity.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <erfa.h>

#include "constants.h"
#include "FEE_primary_beam.h"
#include "MWA_analy_expected_nside201.h"
#include "azza_radec_nside201.h"

#define FREQ 150000000

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

//External CUDA code we're linking in

extern void test_RTS_calculate_MWA_analytic_beam(int num_components,
     int num_time_steps, int num_freqs,
     user_precision_t *azs, user_precision_t *zas, user_precision_t *delays,
     double latitude, int norm,
     user_precision_t *beam_has, user_precision_t *beam_decs, double *freqs,
     user_precision_complex_t *gxs, user_precision_complex_t *Dxs,
     user_precision_complex_t *Dys, user_precision_complex_t *gys);

#ifdef DOUBLE_PRECISION
  double TOL = 1e-8;
#else
  double TOL = 1e-6;
#endif

void test_MWA_analytic_beam_nside201() {


  int nside = 201;

  int num_times = 2;
  int num_freqs = 2;

  double freqs[] = {150e+6, 200e+6};

  int num_components = nside*nside;
  int num_azza = num_components*num_times;
  int num_beam_values = num_azza*num_freqs;

  user_precision_t *beam_has = malloc(num_azza*sizeof(user_precision_t));
  user_precision_t *beam_decs = malloc(num_azza*sizeof(user_precision_t));
  user_precision_t *all_azs = malloc(num_azza*sizeof(user_precision_t));
  user_precision_t *all_zas = malloc(num_azza*sizeof(user_precision_t));

  for (int coord = 0; coord < num_components; coord++) {
    for (int time = 0; time < num_times; time++) {

      beam_has[coord*num_times + time] = nside201_has[coord];
      beam_decs[coord*num_times + time] = nside201_decs[coord];

      all_azs[coord*num_times + time] = nside201_azs[coord];
      all_zas[coord*num_times + time] = nside201_zas[coord];
    }
  }

  user_precision_complex_t *gxs = malloc(num_beam_values*sizeof(user_precision_complex_t));
  user_precision_complex_t *Dxs = malloc(num_beam_values*sizeof(user_precision_complex_t));
  user_precision_complex_t *Dys = malloc(num_beam_values*sizeof(user_precision_complex_t));
  user_precision_complex_t *gys = malloc(num_beam_values*sizeof(user_precision_complex_t));

  int norm = 1;
  // float delays[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  //                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  user_precision_t delays[] = {6, 4, 2, 0, 8, 6, 4, 2,
                               10, 8, 6, 4, 12, 10, 8, 6};
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

//Run the test with unity
int main(void)
{
    UNITY_BEGIN();
    RUN_TEST(test_MWA_analytic_beam_nside201);

    return UNITY_END();
}
