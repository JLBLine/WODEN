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
                  float *us, float *vs, float *ws,
                  double *ls, double *ms, double *ns,
                  float _Complex *visis);

#define UNITY_INCLUDE_FLOAT

/*
Setup 10000 u,v,w coords and 3600 l,m,n coords, and send them off to
kern_calc_measurement_equation to calculate the measurement equation.
Check against a calculation made with C
*/
void test_kern_calc_measurement_equation_ComparedToC(void) {

  // //Set up some test condition inputs

  //Just going off reasonable values for the MWA here
  float max_u = 1000;
  float max_v = 1000;
  float max_w = 100;

  float u_inc = 20;
  float v_inc = 20;
  float w_inc = 0.2;

  int num_us = 2*max_u / u_inc;
  int num_vs = 2*max_v / v_inc;

  int num_baselines = num_us*num_vs;

  float *us = malloc(num_baselines*sizeof(float));
  float *vs = malloc(num_baselines*sizeof(float));
  float *ws = malloc(num_baselines*sizeof(float));

  for (int u_ind = 0; u_ind < num_us; u_ind++) {
    for (int v_ind = 0; v_ind < num_vs; v_ind++) {
      us[u_ind*num_vs + v_ind] = -max_u + u_inc*u_ind;
      vs[u_ind*num_vs + v_ind] = -max_v + v_inc*v_ind;
      ws[u_ind*num_vs + v_ind] = -max_w + w_inc*u_ind;

      // printf("%.3f %.3f %.3f\n",us[u_ind*num_vs + v_ind],
      //                           vs[u_ind*num_vs + v_ind],
      //                           ws[u_ind*num_vs + v_ind] );
    }
  }

  //Setup polar coords to sample to whole l,m sky
  double r_inc = 0.05;
  double phi_inc = 2.0;

  int num_phis = 360.0 / phi_inc;
  int num_rs = 1.0 / r_inc;

  int num_components = num_phis*num_rs;

  printf("Number baselines %d Number Components %d \n",num_baselines,num_components );

  double *ls = malloc(num_components*sizeof(double));
  double *ms = malloc(num_components*sizeof(double));
  double *ns = malloc(num_components*sizeof(double));

  double r, phi;
  double l, m;

  for (int r_ind = 0; r_ind < num_rs; r_ind++) {
    for (int phi_ind = 0; phi_ind < num_phis; phi_ind++) {

      r = r_inc + r_inc*r_ind;
      phi = phi_inc*phi_ind*(M_PI / 180.0);

      l = cos(phi)*r;
      m = sin(phi)*r;

      ls[r_ind*num_phis + phi_ind] = l;
      ms[r_ind*num_phis + phi_ind] = m;
      ns[r_ind*num_phis + phi_ind] = sqrt(1 - l*l - m*m);

      // if (l == 1.0) {
      //   printf("%.6f %.6f %.6f %.1f\n",ls[r_ind*num_phis + phi_ind],
      //                             ms[r_ind*num_phis + phi_ind],
      //                             ns[r_ind*num_phis + phi_ind],
      //                             sqrt(l*l + m*m + ns[r_ind*num_phis + phi_ind]*ns[r_ind*num_phis + phi_ind]));
      // }


    }
  }

  //Space for outputs
  float _Complex *visis = malloc(num_baselines*num_components*sizeof(float _Complex));

  //Run the CUDA code
  test_kern_calc_measurement_equation(num_components, num_baselines,
                                      us, vs, ws, ls, ms, ns, visis);

  //Fill values with what should have been found
  int ind = 0;
  double temp;
  double expec_re, expec_im;

  for (size_t baseline = 0; baseline < num_baselines; baseline++) {
    for (size_t comp = 0; comp < num_components; comp++) {

      //Do the measurement equation using C
      temp = 2.0*M_PI*( (double)us[baseline]*ls[comp] + (double)vs[baseline]*ms[comp] + (double)ws[baseline]*(ns[comp]-1) );
      sincos(temp, &(expec_im), &(expec_re));

      // printf("%.10f %.10f %.10f %.10f %.10f %.10f\n",us[baseline],ls[comp],vs[baseline],ms[comp],ws[baseline],ns[comp] );
      // printf("%.7f %.7f %.7f %.7f\n",creal(visis[ind]), expec_re, tol*abs(expec_re), creal(visis[ind]) - expec_re  );
      // printf("%.7f %.7f %.7f %.7f\n",cimag(visis[ind]), expec_im, tol*abs(expec_im), cimag(visis[ind]) - expec_im  );

      float tol = 1e-6;

      //Check within some tolerance
      TEST_ASSERT_FLOAT_WITHIN(tol, (float)expec_re, creal(visis[ind]));
      TEST_ASSERT_FLOAT_WITHIN(tol, (float)expec_im, cimag(visis[ind]));

      ind ++;
    }
  }

  //Be free my beauties
  free(us);
  free(vs);
  free(ws);
  free(ls);
  free(ms);
  free(ns);
  free(visis);

}

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

  float target_angles[] = {0.0, M_PI/6.0, M_PI/4.0, M_PI/3.0, M_PI/2.0,
                           2*M_PI/3, 3*M_PI/4, 5*M_PI/6, M_PI,
                           7*M_PI/6, 5*M_PI/4};

  float expec_ims[] = {0.0, 0.5, sqrt(2)/2.0, sqrt(3.0)/2.0, 1.0,
                       sqrt(3.0)/2.0, sqrt(2)/2.0, 0.5, 0.0,
                       -0.5, -sqrt(2)/2.0};
  float expec_res[] = {1.0, sqrt(3.0)/2.0, sqrt(2)/2.0, 0.5, 0.0,
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

  float *us = malloc(num_baselines*sizeof(float));
  float *vs = malloc(num_baselines*sizeof(float));
  float *ws = malloc(num_baselines*sizeof(float));

  // float multipliers[] = {1000.0, 10000.0};
  float multipliers[] = {1.0, 10.0, 100.0, 1000.0, 10000.0};
  float uvw;

  int num_components = 1;

  double *ls = malloc(sizeof(double));
  double *ms = malloc(sizeof(double));
  double *ns = malloc(sizeof(double));

  //Space for outputs
  float _Complex *visis = malloc(num_baselines*num_components*sizeof(float _Complex));

  float re_diff_sum = 0.0;
  float im_diff_sum = 0.0;
  int num_in_sum = 0;

  float re_diff;
  float im_diff;

  FILE *output_text;
  // char buff[0x100];
  // snprintf(buff, sizeof(buff), "measurement_eq_outcomes.txt");
  output_text = fopen("measurement_eq_outcomes.txt","w");

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
    // int ind = 0;
    // float temp;
    // float expec_re, expec_im;

    // printf("=======Expected Real %.7f, Imag %.7f=======================================\n", expec_res[comp], expec_ims[comp]);

    for (size_t baseline = 0; baseline < num_baselines; baseline++) {

      // temp = 2.0*M_PI*( us[baseline]*ls[0] + vs[baseline]*ms[0] + ws[baseline]*(ns[0]-1) );
      //
      // expec_im = sinf(temp);
      // expec_re = cosf(temp);

      //For those interested
      // printf("uvw %.1f \t Real C %.9f CUDA %.9f Diff C %.9f CUDA %.9f\n",us[baseline], expec_res[comp], creal(visis[baseline]), expec_res[comp] - expec_res[comp], creal(visis[baseline]) - expec_res[comp]); //,
      // printf("uvw %.1f \t Imag C %.9f CUDA %.9f Diff C %.9f CUDA %.9f\n",us[baseline], expec_ims[comp], cimag(visis[baseline]), expec_ims[comp] - expec_ims[comp], cimag(visis[baseline]) - expec_ims[comp]); //,

      float tol = 2e-3;

      re_diff_sum += fabs(creal(visis[baseline]) - expec_res[comp]);
      im_diff_sum += fabs(cimag(visis[baseline]) - expec_ims[comp]);
      num_in_sum += 1;

      fprintf(output_text,"%.1f %.16f %.16f %.16f %.16f\n",us[baseline], expec_res[comp],
                                                           creal(visis[baseline]),
                                                           expec_ims[comp],
                                                           cimag(visis[baseline]));


      if (expec_res[comp] == 0.0 || expec_ims[comp] == 0.0) {

        TEST_ASSERT_FLOAT_WITHIN(tol, expec_res[comp], creal(visis[baseline]));
        TEST_ASSERT_FLOAT_WITHIN(tol, expec_ims[comp], cimag(visis[baseline]));
      } else {
        // //Check within some tolerance
        TEST_ASSERT_FLOAT_WITHIN(tol*expec_res[comp], expec_res[comp], creal(visis[baseline]));
        TEST_ASSERT_FLOAT_WITHIN(tol*expec_ims[comp], expec_ims[comp], cimag(visis[baseline]));

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

}

//Run the test with unity
int main(void)
{
    UNITY_BEGIN();
    RUN_TEST(test_kern_calc_measurement_equation_ComparedToC);
    RUN_TEST(test_kern_calc_measurement_equation_GiveCorrectValues);

    return UNITY_END();
}
