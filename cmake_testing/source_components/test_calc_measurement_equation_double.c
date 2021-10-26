#include <math.h>
#include <unity.h>
#include <stdlib.h>
#include <complex.h>

#include "constants.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

void sincosf(float x, float *sin, float *cos);
void sincos(double x, double *sin, double *cos);

//External CUDA code we're linking in
extern void test_kern_calc_measurement_equation_double(int num_components,
                  int num_baselines,
                  double *us, double *vs, double *ws,
                  double *ls, double *ms, double *ns,
                  double _Complex *visis);

#define UNITY_INCLUDE_FLOAT

/*
Use a few exact coordinates that should give exact answers, and test
whether they are within some tolerace
*/
void test_kern_calc_measurement_equation_GiveCorrectValues(void) {

  //These combinations of l,m,n, when paired with u,v,w, result in an
  //a product of 2*pi*(u*l + v*m + w*(n-1)) equal to:
  //These angles have known outcomes for sine / cosine so can test
  //the accuracy of our measurement equation.

  int all_components = 11;

  // double all_ls[] = {0.0000000000,0.0425737516,0.0645903245,0.0871449864,0.1340695840,0.1838657911,0.2100755148,0.2373397983,0.2958758548,0.3622725654,0.4003681254};
  //
  // double all_ms[] = {0.0000000000,0.0425737516,0.0645903245,0.0871449864,0.1340695840,0.1838657911,0.2100755148,0.2373397983,0.2958758548,0.3622725654,0.4003681254};
  //
  // double all_ns[] = {1.0000000000,0.9981858301,0.9958193511,0.9923766940,0.9818608319,0.9656017511,0.9548489703,0.9419870701,0.9082482905,0.8587882024,0.8242637493};

  double target_angles[] = {0.0, M_PI/6.0, M_PI/4.0, M_PI/3.0, M_PI/2.0,
                           2*M_PI/3, 3*M_PI/4, 5*M_PI/6, M_PI,
                           7*M_PI/6, 5*M_PI/4};

  double expec_ims[] = {0.0, 0.5, sqrt(2)/2.0, sqrt(3.0)/2.0, 1.0,
                       sqrt(3.0)/2.0, sqrt(2)/2.0, 0.5, 0.0,
                       -0.5, -sqrt(2)/2.0};
  double expec_res[] = {1.0, sqrt(3.0)/2.0, sqrt(2)/2.0, 0.5, 0.0,
                       -0.5, -sqrt(2)/2.0, -sqrt(3.0)/2.0, -1.0,
                       -sqrt(3.0)/2.0, -sqrt(2)/2.0};

  double all_ls[] = {0.0000000000000000,0.0425737516338956,0.0645903244635131,0.0871449863555500,
                    0.1340695840364469,0.1838657911209207,0.2100755148372292,0.2373397982598921,
                    0.2958758547680685,0.3622725654470420,0.4003681253515569};

  double all_ms[] = {0.0000000000000000,0.0425737516338956,0.0645903244635131,0.0871449863555500,
                    0.1340695840364469,0.1838657911209207,0.2100755148372292,0.2373397982598921,
                    0.2958758547680685,0.3622725654470420,0.4003681253515569};

  double all_ns[] = {1.0000000000000000,0.9981858300655398,0.9958193510729726,0.9923766939555675,
                    0.9818608319271057,0.9656017510914922,0.9548489703255412,0.9419870701468823,
                    0.9082482904638630,0.8587882024392495,0.8242637492968862};

  //We'll set up a number of u,v,w that that should result in adding multiples
  //of 2*pi to the desired target angle, for each target angle, so we have
  //to loop over each over the l coords

  int num_baselines = 5;

  double *us = malloc(num_baselines*sizeof(double));
  double *vs = malloc(num_baselines*sizeof(double));
  double *ws = malloc(num_baselines*sizeof(double));

  // double multipliers[] = {10000.0};
  double multipliers[] = {1.0, 10.0, 100.0, 1000.0, 10000.0};
  double uvw;

  int num_components = 1;

  double *ls = malloc(sizeof(double));
  double *ms = malloc(sizeof(double));
  double *ns = malloc(sizeof(double));

  //Space for outputs
  float _Complex *visis = malloc(num_baselines*num_components*sizeof(double _Complex));

  double re_diff_sum = 0.0;
  double im_diff_sum = 0.0;
  int num_in_sum = 0;

  double re_diff;
  double im_diff;

  FILE *output_text;
  // char buff[0x100];
  // snprintf(buff, sizeof(buff), "measurement_eq_outcomes.txt");
  output_text = fopen("measurement_eq_outcomes_double.txt","w");

  for (int comp = 0; comp < all_components; comp++) {
    //Setup a number of u,v,w coords that should result in adding multiples
    //of 2*pi - this means we should keep retrieving the same output from
    //the visibilities

    ls[0] = all_ls[comp];
    ms[0] = all_ms[comp];
    ns[0] = all_ns[comp];

    for (int baseline = 0; baseline < num_baselines; baseline++) {

      if (target_angles[comp] == 0) {
        uvw = 2*M_PI*multipliers[baseline];
        // uvw = 1.0;
      }
      else{
        uvw = (target_angles[comp] + 2*M_PI*multipliers[baseline]) / target_angles[comp];
      }

      // printf("%.1f\n",uvw );

      us[baseline] = 0.0;
      vs[baseline] = 0.0;
      ws[baseline] = 0.0;

      us[baseline] = uvw;
      vs[baseline] = uvw;
      ws[baseline] = uvw;

      //Ensure visis are set to zero before going into CUDA code
      visis[baseline] = 0.0 + I*0.0;

    }

    //Run the CUDA code
    test_kern_calc_measurement_equation_double(num_components, num_baselines,
                                        us, vs, ws, ls, ms, ns, visis);

    for (size_t baseline = 0; baseline < num_baselines; baseline++) {

      // printf("uvw %.1f \t Real C %.9f CUDA %.9f Diff C %.9f CUDA %.9f\n",us[baseline], expec_re, creal(visis[baseline]), expec_re - expec_res[comp], creal(visis[baseline]) - expec_res[comp]); //,
      // printf("uvw %.1f \t Imag C %.9f CUDA %.9f Diff C %.9f CUDA %.9f\n",us[baseline], expec_im, cimag(visis[baseline]), expec_im - expec_ims[comp], cimag(visis[baseline]) - expec_ims[comp]); //,

      printf("uvw %.1f \t Real expec %.9f CUDA %.9f Diff %.9f\n",us[baseline], expec_res[comp], creal(visis[baseline]),
                                                                         creal(visis[baseline]) - expec_res[comp]); //,
      printf("uvw %.1f \t Imag expec %.9f CUDA %.9f Diff %.9f\n",us[baseline], expec_ims[comp], cimag(visis[baseline]),
                                                                         cimag(visis[baseline]) - expec_ims[comp]); //,

      float tol = 1e-4;

      re_diff_sum += abs(creal(visis[baseline]) - expec_res[comp]);
      im_diff_sum += abs(cimag(visis[baseline]) - expec_ims[comp]);
      num_in_sum += 1;

      fprintf(output_text,"%.1f %.16f %.16f %.16f %.16f\n",us[baseline], expec_res[comp],
                                                           creal(visis[baseline]),
                                                           expec_ims[comp],
                                                           cimag(visis[baseline]));


      if (expec_res[comp] == 0.0 || expec_ims[comp] == 0.0) {

        TEST_ASSERT_FLOAT_WITHIN(tol, (float)expec_res[comp], creal(visis[baseline]));
        TEST_ASSERT_FLOAT_WITHIN(tol, (float)expec_ims[comp], cimag(visis[baseline]));
      } else {
        // //Check within some tolerance
        TEST_ASSERT_FLOAT_WITHIN(tol*expec_res[comp], (float)expec_res[comp], creal(visis[baseline]));
        TEST_ASSERT_FLOAT_WITHIN(tol*expec_ims[comp], (float)expec_ims[comp], cimag(visis[baseline]));

      }
    }
  }

  printf("Mean real offset %.6f imag offset %.6f\n",re_diff_sum / num_in_sum, im_diff_sum / num_in_sum );
  printf("Num samples %d\n",num_in_sum );

  fflush(output_text);
  fclose(output_text);

  //Be free my beauties
  free(us);
  free(vs);
  free(ws);
  free(ls);
  free(ms);
  free(ns);
  free(visis);
  // free(expec_res);
  // free(expec_ims);
  // free(visi_res);
  // free(visi_ims);

}

//Run the test with unity
int main(void)
{
    UNITY_BEGIN();
    // RUN_TEST(test_kern_calc_measurement_equation_ComparedToC);

    RUN_TEST(test_kern_calc_measurement_equation_GiveCorrectValues);

    return UNITY_END();
}
