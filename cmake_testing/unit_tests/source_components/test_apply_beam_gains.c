#include <math.h>
#include <unity.h>
#include <stdlib.h>
#include <complex.h>

#include "constants.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

//External CUDA code we're linking in
extern void test_kern_apply_beam_gains(int num_gains, float _Complex *g1xs,
          float _Complex *D1xs,
          float _Complex *D1ys, float _Complex *g1ys,
          float _Complex *g2xs, float _Complex *D2xs,
          float _Complex *D2ys, float _Complex *g2ys,
          float *flux_Is, float *flux_Qs,
          float *flux_Us, float *flux_Vs,
          float _Complex *visi_components,
          float _Complex *visi_XXs, float _Complex *visi_XYs,
          float _Complex *visi_YXs, float _Complex *visi_YYs);

#define UNITY_INCLUDE_FLOAT

/*
Test that the code applying beam gains to Stokes parameter fluxes works
*/
void test_kern_apply_beam_gains_GiveCorrectValues(void) {

  //Set up some test condition inputs
  int num_gains = 16;

  //Set some fun gain params
  float _Complex g1xs[] = {1.0 + I*0.0, 1.0 + I*0.0, 1.0 + I*0.0, 1.0 + I*0.0,
                           0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0,
                           2.0 + I*0.0, 2.0 + I*0.0, 2.0 + I*0.0, 2.0 + I*0.0,
                           1.0 + I*2.0, 1.0 + I*2.0, 1.0 + I*2.0, 1.0 + I*2.0};

  float _Complex D1xs[] = {0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0,
                           1.0 + I*0.0, 1.0 + I*0.0, 1.0 + I*0.0, 1.0 + I*0.0,
                           2.0 + I*0.0, 2.0 + I*0.0, 2.0 + I*0.0, 2.0 + I*0.0,
                           3.0 + I*4.0, 3.0 + I*4.0, 3.0 + I*4.0, 3.0 + I*4.0};

  float _Complex D1ys[] = {0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0,
                           1.0 + I*0.0, 1.0 + I*0.0, 1.0 + I*0.0, 1.0 + I*0.0,
                           2.0 + I*0.0, 2.0 + I*0.0, 2.0 + I*0.0, 2.0 + I*0.0,
                           5.0 + I*6.0, 5.0 + I*6.0, 5.0 + I*6.0, 5.0 + I*6.0};

  float _Complex g1ys[] = {1.0 + I*0.0, 1.0 + I*0.0, 1.0 + I*0.0, 1.0 + I*0.0,
                           0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0,
                           2.0 + I*0.0, 2.0 + I*0.0, 2.0 + I*0.0, 2.0 + I*0.0,
                           7.0 + I*8.0, 7.0 + I*8.0, 7.0 + I*8.0, 7.0 + I*8.0};

  //Cycle through the Stokes parameters
  float *flux_Is = calloc(num_gains, sizeof(float));
  float *flux_Qs = calloc(num_gains, sizeof(float));
  float *flux_Us = calloc(num_gains, sizeof(float));
  float *flux_Vs = calloc(num_gains, sizeof(float));

  flux_Is[0] = flux_Is[4] = flux_Is[8] = flux_Is[12] = 1.0;
  flux_Qs[1] = flux_Qs[5] = flux_Qs[9] = flux_Qs[13] = 1.0;
  flux_Us[2] = flux_Us[6] = flux_Us[10] = flux_Us[14] = 1.0;
  flux_Vs[3] = flux_Vs[7] = flux_Vs[11] = flux_Vs[15] = 1.0;

  //Just keep the visibility component equal to 1.0 + 0.0j for all entries
  float _Complex *visi_components = malloc(num_gains*sizeof(float _Complex));

  for (size_t gain = 0; gain < num_gains; gain++) {
    visi_components[gain] = 1.0 + I*0.0;
  }

  //Space for outputs
  float _Complex *visi_XXs = malloc(num_gains*sizeof(float _Complex));
  float _Complex *visi_XYs = malloc(num_gains*sizeof(float _Complex));
  float _Complex *visi_YXs = malloc(num_gains*sizeof(float _Complex));
  float _Complex *visi_YYs = malloc(num_gains*sizeof(float _Complex));

  //Run the CUDA code
  test_kern_apply_beam_gains(num_gains, g1xs, D1xs, D1ys, g1ys,
                                        g1xs, D1xs, D1ys, g1ys,
                                        flux_Is, flux_Qs,
                                        flux_Us, flux_Vs,
                                        visi_components,
                                        visi_XXs, visi_XYs,
                                        visi_YXs, visi_YYs);

  //Expected outputs
  float expec_XX_re[] =  {1.0, 1.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0,
                          8.0, 0.0, 8.0, 0.0, 30.0, -20.0, 22.0, -4.0};
  float expec_XY_re[] =  {0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
                          8.0, 0.0, 8.0, 0.0, 70.0, -36.0, 62.0, -4.0};
  float expec_YX_re[] =  {0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
                          8.0, 0.0, 8.0, 0.0, 70.0, -36.0, 62.0, -4.0};
  float expec_YY_re[] =  {1.0, -1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0,
                          8.0, 0.0, 8.0, 0.0, 174.0, -52.0, 166.0, -4.0};
  float expec_XX_im[] =  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  float expec_XY_im[] =  {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0,
                          0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 8.0, -16.0};
  float expec_YX_im[] =  {0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0,
                          0.0, 0.0, 0.0, 0.0, -8.0, 0.0, -8.0, 16.0};
  float expec_YY_im[] =  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  // for (size_t gain = 0; gain < num_gains; gain++) {
  //   printf("%.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n",
  //           creal(visi_XXs[gain]),cimag(visi_XXs[gain]),
  //           creal(visi_XYs[gain]),cimag(visi_XYs[gain]),
  //           creal(visi_YXs[gain]),cimag(visi_YXs[gain]),
  //           creal(visi_YYs[gain]),cimag(visi_YYs[gain]) );
  // }

  //Check outputs match expectations
  for (size_t gain = 0; gain < num_gains; gain++) {
    TEST_ASSERT_EQUAL_FLOAT(expec_XX_re[gain], creal(visi_XXs[gain]) );
    TEST_ASSERT_EQUAL_FLOAT(expec_XX_im[gain], cimag(visi_XXs[gain]) );
    TEST_ASSERT_EQUAL_FLOAT(expec_XY_re[gain], creal(visi_XYs[gain]) );
    TEST_ASSERT_EQUAL_FLOAT(expec_XY_im[gain], cimag(visi_XYs[gain]) );
    TEST_ASSERT_EQUAL_FLOAT(expec_YX_re[gain], creal(visi_YXs[gain]) );
    TEST_ASSERT_EQUAL_FLOAT(expec_YX_im[gain], cimag(visi_YXs[gain]) );
    TEST_ASSERT_EQUAL_FLOAT(expec_YY_re[gain], creal(visi_YYs[gain]) );
    TEST_ASSERT_EQUAL_FLOAT(expec_YY_im[gain], cimag(visi_YYs[gain]) );
  }

  // //Be free my beauties
  free(flux_Is);
  free(flux_Qs);
  free(flux_Us);
  free(flux_Vs);
  free(visi_XXs);
  free(visi_XYs);
  free(visi_YXs);
  free(visi_YYs);

}

//Run the test with unity
int main(void)
{
    UNITY_BEGIN();
    RUN_TEST(test_kern_apply_beam_gains_GiveCorrectValues);
    return UNITY_END();
}
