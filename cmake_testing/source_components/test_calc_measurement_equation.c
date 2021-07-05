#include <math.h>
#include <unity.h>
#include <stdlib.h>
#include <complex.h>

#include "constants.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

void sincosf(float x, float *sin, float *cos);

//External CUDA code we're linking in
extern void test_kern_calc_measurement_equation(int num_components,
                  int num_baselines,
                  float *us, float *vs, float *ws,
                  float *ls, float *ms, float *ns,
                  float _Complex *visis);

#define UNITY_INCLUDE_FLOAT

/*
Test that the linear SI flux extrapolation code works correctly
*/
void test_kern_calc_measurement_equation_GiveCorrectValues(void) {

  //Set up some test condition inputs
  int num_baselines = 1000;
  int num_components = 100;

  //This is a reasonable range based on MWA
  float max_u = 1000.0;
  float max_v = 1000.0;
  float max_w = 20.0;

  float u_inc = (2*max_u) / num_baselines;
  float v_inc = (2*max_v) / num_baselines;
  float w_inc = (2*max_w) / num_baselines;

  float *us = malloc(num_baselines*sizeof(float));
  float *vs = malloc(num_baselines*sizeof(float));
  float *ws = malloc(num_baselines*sizeof(float));

  for (size_t baseline = 0; baseline < num_baselines; baseline++) {
    us[baseline] = -max_u + u_inc*baseline;
    vs[baseline] = -max_v + v_inc*baseline;
    ws[baseline] = -max_w + w_inc*baseline;
  }

  float *ls = malloc(num_components*sizeof(float));
  float *ms = malloc(num_components*sizeof(float));
  float *ns = malloc(num_components*sizeof(float));

  float l_inc = 1.4 / num_components;
  float m_inc = 1.4 / num_components;
  float l, m;

  for (size_t component = 0; component < num_components; component++) {
    l =  -0.7 + l_inc*component;
    m =  -0.7 + m_inc*component;
    ls[component] = l;
    ms[component] = m;
    ns[component] = sqrt(1 - l*l - m*m);
  }


  //Space for outputs
  float _Complex *visis = malloc(num_baselines*num_components*sizeof(float _Complex));

  //Run the CUDA code
  test_kern_calc_measurement_equation(num_components, num_baselines,
                                      us, vs, ws, ls, ms, ns, visis);

  //Make some expected value arrays
  float *expec_res = malloc(num_baselines*num_components*sizeof(float));
  float *expec_ims = malloc(num_baselines*num_components*sizeof(float));

  float *visi_res = malloc(num_baselines*num_components*sizeof(float));
  float *visi_ims = malloc(num_baselines*num_components*sizeof(float));

  //Fill values with what should have been found
  int ind = 0;
  float temp, expec_re, expec_im;

  for (size_t baseline = 0; baseline < num_baselines; baseline++) {
    for (size_t comp = 0; comp < num_components; comp++) {

      temp = 2*M_PI*( us[baseline]*ls[comp] + vs[baseline]*ms[comp] + ws[baseline]*(ns[comp]-1) );
      sincosf(temp, &(expec_im), &(expec_re));

      expec_res[ind] = expec_re;
      expec_ims[ind] = expec_im;

      visi_res[ind] = creal(visis[ind]);
      visi_ims[ind] = cimag(visis[ind]);

      //Check within some tolerance
      TEST_ASSERT_FLOAT_WITHIN(5e-3, expec_re, creal(visis[ind]));
      TEST_ASSERT_FLOAT_WITHIN(5e-3, expec_im, cimag(visis[ind]));

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
  free(expec_res);
  free(expec_ims);
  free(visi_res);
  free(visi_ims);

}

//Run the test with unity
int main(void)
{
    UNITY_BEGIN();
    RUN_TEST(test_kern_calc_measurement_equation_GiveCorrectValues);
    return UNITY_END();
}
