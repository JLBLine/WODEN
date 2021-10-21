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
extern void test_kern_calc_measurement_equation(int num_components,
                  int num_baselines,
                  double *us, double *vs, double *ws,
                  double *ls, double *ms, double *ns,
                  double _Complex *visis);

#define UNITY_INCLUDE_FLOAT

// /*
// Setup 1000 u,v,w coords and 100 l,m,n coords, and send them off to
// kern_calc_measurement_equation to calculate the measurement equation.
// Check against a calculation made with C
// */
// void test_kern_calc_measurement_equation_ComparedToC(void) {
//
//   //Set up some test condition inputs
//   int num_baselines = 1000;
//   int num_components = 100;
//
//   //This is a reasonable range based on MWA
//   double max_u = 1000.0;
//   double max_v = 1000.0;
//   double max_w = 20.0;
//
//   // double u_inc = (2*max_u) / num_baselines;
//   // double v_inc = (2*max_v) / num_baselines;
//   // double w_inc = (2*max_w) / num_baselines;
//
//   double u_inc = 2;
//   double v_inc = 2;
//   double w_inc = 0.04;
//
//   double *us = malloc(num_baselines*sizeof(double));
//   double *vs = malloc(num_baselines*sizeof(double));
//   double *ws = malloc(num_baselines*sizeof(double));
//
//   for (size_t baseline = 0; baseline < num_baselines; baseline++) {
//     us[baseline] = -max_u + u_inc*baseline;
//     vs[baseline] = -max_v + v_inc*baseline;
//     ws[baseline] = -max_w + w_inc*baseline;
//   }
//
//   double *ls = malloc(num_components*sizeof(double));
//   double *ms = malloc(num_components*sizeof(double));
//   double *ns = malloc(num_components*sizeof(double));
//
//   // double l_inc = 1.4 / num_components;
//   // double m_inc = 1.4 / num_components;
//   double l_inc = 0.014;
//   double m_inc = 0.014;
//   double l, m;
//
//   for (size_t component = 0; component < num_components; component++) {
//     l =  -0.7 + l_inc*component;
//     m =  -0.7 + m_inc*component;
//     ls[component] = l;
//     ms[component] = m;
//     ns[component] = sqrt(1 - l*l - m*m);
//   }
//
//
//   //Space for outputs
//   double _Complex *visis = malloc(num_baselines*num_components*sizeof(double _Complex));
//
//   //Run the CUDA code
//   test_kern_calc_measurement_equation(num_components, num_baselines,
//                                       us, vs, ws, ls, ms, ns, visis);
//
//   //Make some expected value arrays
//   double *expec_res = malloc(num_baselines*num_components*sizeof(double));
//   double *expec_ims = malloc(num_baselines*num_components*sizeof(double));
//
//   double *visi_res = malloc(num_baselines*num_components*sizeof(double));
//   double *visi_ims = malloc(num_baselines*num_components*sizeof(double));
//
//   //Fill values with what should have been found
//   int ind = 0;
//   double temp;
//   double expec_re, expec_im;
//
//   for (size_t baseline = 0; baseline < num_baselines; baseline++) {
//     for (size_t comp = 0; comp < num_components; comp++) {
//
//       temp = 2.0*M_PI*( us[baseline]*ls[comp] + vs[baseline]*ms[comp] + ws[baseline]*(ns[comp]-1) );
//       sincos((double)temp, &(expec_im), &(expec_re));
//
//       expec_res[ind] = expec_re;
//       expec_ims[ind] = expec_im;
//
//       visi_res[ind] = creal(visis[ind]);
//       visi_ims[ind] = cimag(visis[ind]);
//
//       double tol = 1e-3; //5e-3
//
//       // printf("%.10f %.10f %.10f %.10f %.10f %.10f\n",us[baseline],ls[comp],vs[baseline],ms[comp],ws[baseline],ns[comp] );
//       // printf("%.7f %.7f %.7f %.7f\n",creal(visis[ind]), expec_re, tol*fabs(expec_re), creal(visis[ind]) - expec_re  );
//       // printf("%.7f %.7f %.7f %.7f\n",cimag(visis[ind]), expec_im, tol*fabs(expec_im), cimag(visis[ind]) - expec_im  );
//
//       printf("%.9f %.9f\n", creal(visis[ind]), expec_re); //, tol*fabs(expec_re), creal(visis[ind]) - expec_re  );
//       printf("%.9f %.9f\n", cimag(visis[ind]), expec_im); //, tol*fabs(expec_im), cimag(visis[ind]) - expec_im  );
//
//       double re_tol = tol*abs(expec_re);
//
//       //Check within some tolerance
//       // TEST_ASSERT_FLOAT_WITHIN(re_tol, expec_re, creal(visis[ind]));
//       // TEST_ASSERT_FLOAT_WITHIN(tol*fabs(expec_im), expec_im, cimag(visis[ind]));
//
//       ind ++;
//     }
//   }
//
//   //Be free my beauties
//   free(us);
//   free(vs);
//   free(ws);
//   free(ls);
//   free(ms);
//   free(ns);
//   free(visis);
//   free(expec_res);
//   free(expec_ims);
//   free(visi_res);
//   free(visi_ims);
//
// }

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

  // double all_ls[] = {0.0, 0.042573751633895596, 0.06459032446351305, 0.08714498635555, 0.13406958403644692,
  //                   0.18386579112092066, 0.21007551483722917, 0.2373397982598921, 0.2958758547680685,
  //                   0.362272565447042, 0.4003681253515569};
  //
  // double all_ms[] = {0.0, 0.042573751633895596, 0.06459032446351305, 0.08714498635555, 0.13406958403644692,
  //                   0.18386579112092066, 0.21007551483722917, 0.2373397982598921, 0.2958758547680685,
  //                   0.362272565447042, 0.4003681253515569};
  //
  // double all_ns[] = {1.0, 0.9981858300655398, 0.9958193510729726, 0.9923766939555675, 0.9818608319271057,
  //                   0.9656017510914922, 0.9548489703255412, 0.9419870701468823, 0.908248290463863,
  //                   0.8587882024392495, 0.8242637492968862};

  double all_ls[] = {0.0000000000,0.0425737516,0.0645903245,0.0871449864,0.1340695840,0.1838657911,0.2100755148,0.2373397983,0.2958758548,0.3622725654,0.4003681254};

  double all_ms[] = {0.0000000000,0.0425737516,0.0645903245,0.0871449864,0.1340695840,0.1838657911,0.2100755148,0.2373397983,0.2958758548,0.3622725654,0.4003681254};

  double all_ns[] = {1.0000000000,0.9981858301,0.9958193511,0.9923766940,0.9818608319,0.9656017511,0.9548489703,0.9419870701,0.9082482905,0.8587882024,0.8242637493};


  double target_angles[] = {0.0, M_PI/6.0, M_PI/4.0, M_PI/3.0, M_PI/2.0,
                           2*M_PI/3, 3*M_PI/4, 5*M_PI/6, M_PI,
                           7*M_PI/6, 5*M_PI/4};

  double expec_ims[] = {0.0, 0.5, sqrt(2)/2.0, sqrt(3.0)/2.0, 1.0,
                       sqrt(3.0)/2.0, sqrt(2)/2.0, 0.5, 0.0,
                       -0.5, -sqrt(2)/2.0};
  double expec_res[] = {1.0, sqrt(3.0)/2.0, sqrt(2)/2.0, 0.5, 0.0,
                       -0.5, -sqrt(2)/2.0, -sqrt(3.0)/2.0, -1.0,
                       -sqrt(3.0)/2.0, -sqrt(2)/2.0};

  //We'll set up a number of u,v,w that that should result in adding multiples
  //of 2*pi to the desired target angle, for each target angle, so we have
  //to loop over each over the l coords

  int num_baselines = 6;

  double *us = malloc(num_baselines*sizeof(double));
  double *vs = malloc(num_baselines*sizeof(double));
  double *ws = malloc(num_baselines*sizeof(double));

  // double multipliers[] = {1000.0, 10000.0};
  double multipliers[] = {0.0, 1.0, 10.0, 100.0, 1000.0, 10000.0};
  double uvw;

  int num_components = 1;

  double *ls = malloc(sizeof(double));
  double *ms = malloc(sizeof(double));
  double *ns = malloc(sizeof(double));

  //Space for outputs
  double _Complex *visis = malloc(num_baselines*num_components*sizeof(double _Complex));

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
    test_kern_calc_measurement_equation(num_components, num_baselines,
                                        us, vs, ws, ls, ms, ns, visis);

    //Fill values with what should have been found
    int ind = 0;
    double temp;
    double expec_re, expec_im;

    // printf("=======Expected Real %.7f, Imag %.7f=======================================\n", expec_res[comp], expec_ims[comp]);

    for (size_t baseline = 0; baseline < num_baselines; baseline++) {

      temp = 2.0*M_PI*( us[baseline]*ls[0] + vs[baseline]*ms[0] + ws[baseline]*(ns[0]-1) );

      expec_im = sin(temp);
      expec_re = cos(temp);

      //For those interested

      // printf("uvw %.1f \t Real C %.9f CUDA %.9f Diff C %.9f CUDA %.9f\n",us[baseline], expec_re, creal(visis[baseline]), expec_re - expec_res[comp], creal(visis[baseline]) - expec_res[comp]); //,
      // printf("uvw %.1f \t Imag C %.9f CUDA %.9f Diff C %.9f CUDA %.9f\n",us[baseline], expec_im, cimag(visis[baseline]), expec_im - expec_ims[comp], cimag(visis[baseline]) - expec_ims[comp]); //,

      double tol = 5e-2;

      re_diff_sum += abs(creal(visis[baseline]) - expec_res[comp]);
      im_diff_sum += abs(cimag(visis[baseline]) - expec_ims[comp]);
      num_in_sum += 1;

      fprintf(output_text,"%.1f %.10f %.10f %.10f %.10f\n",us[baseline], expec_res[comp],
                                                           creal(visis[baseline]),
                                                           expec_ims[comp],
                                                           cimag(visis[baseline]));


      // if (expec_res[comp] == 0.0 || expec_ims[comp] == 0.0) {
      //
      //   TEST_ASSERT_FLOAT_WITHIN(tol, expec_res[comp], creal(visis[baseline]));
      //   TEST_ASSERT_FLOAT_WITHIN(tol, expec_ims[comp], cimag(visis[baseline]));
      // } else {
      //   // //Check within some tolerance
      //   TEST_ASSERT_FLOAT_WITHIN(tol*expec_res[comp], expec_res[comp], creal(visis[baseline]));
      //   TEST_ASSERT_FLOAT_WITHIN(tol*expec_ims[comp], expec_ims[comp], cimag(visis[baseline]));
      //
      //   // re_diff_sum += fabs(creal(visis[baseline]) - expec_res[comp]) / expec_res[comp];
      //   // im_diff_sum += fabs(cimag(visis[baseline]) - expec_ims[comp]) / expec_ims[comp];
      //   // num_in_sum += 1;
      // }
      // ind ++;
    }
    // }
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
