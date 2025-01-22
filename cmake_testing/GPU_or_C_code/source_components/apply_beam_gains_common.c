#include "apply_beam_gains_common.h"

//External CUDA code we're linking in
extern void test_kern_apply_beam_gains(int num_gains, user_precision_complex_t *g1xs,
          user_precision_complex_t *D1xs,
          user_precision_complex_t *D1ys, user_precision_complex_t *g1ys,
          user_precision_complex_t *g2xs, user_precision_complex_t *D2xs,
          user_precision_complex_t *D2ys, user_precision_complex_t *g2ys,
          user_precision_t *flux_Is, user_precision_t *flux_Qs,
          user_precision_t *flux_Us, user_precision_t *flux_Vs,
          user_precision_complex_t *visi_components,
          user_precision_complex_t *visi_XXs, user_precision_complex_t *visi_XYs,
          user_precision_complex_t *visi_YXs, user_precision_complex_t *visi_YYs,
          int off_cardinal, int full_pol);

#define UNITY_INCLUDE_FLOAT

/*
Test that the code applying beam gains to Stokes parameter fluxes works
*/
void test_apply_beam_gains_GiveCorrectValues(int do_gpu, int off_cardinal, int full_pol) {

  // #ifdef DOUBLE_PRECISION
  // printf("WODEN is using DOUBLE precision\n");
  // #else
  // printf("WODEN is using FLOAT precision\n");
  // #endif

  //Set up some test condition inputs
  int num_gains = 16;

  //Set some fun gain params
  user_precision_complex_t g1xs[] = {1.0 + I*0.0, 1.0 + I*0.0, 1.0 + I*0.0, 1.0 + I*0.0,
                           0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0,
                           2.0 + I*0.0, 2.0 + I*0.0, 2.0 + I*0.0, 2.0 + I*0.0,
                           1.0 + I*2.0, 1.0 + I*2.0, 1.0 + I*2.0, 1.0 + I*2.0};

  user_precision_complex_t D1xs[] = {0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0,
                           1.0 + I*0.0, 1.0 + I*0.0, 1.0 + I*0.0, 1.0 + I*0.0,
                           2.0 + I*0.0, 2.0 + I*0.0, 2.0 + I*0.0, 2.0 + I*0.0,
                           3.0 + I*4.0, 3.0 + I*4.0, 3.0 + I*4.0, 3.0 + I*4.0};

  user_precision_complex_t D1ys[] = {0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0,
                           1.0 + I*0.0, 1.0 + I*0.0, 1.0 + I*0.0, 1.0 + I*0.0,
                           2.0 + I*0.0, 2.0 + I*0.0, 2.0 + I*0.0, 2.0 + I*0.0,
                           5.0 + I*6.0, 5.0 + I*6.0, 5.0 + I*6.0, 5.0 + I*6.0};

  user_precision_complex_t g1ys[] = {1.0 + I*0.0, 1.0 + I*0.0, 1.0 + I*0.0, 1.0 + I*0.0,
                           0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0, 0.0 + I*0.0,
                           2.0 + I*0.0, 2.0 + I*0.0, 2.0 + I*0.0, 2.0 + I*0.0,
                           7.0 + I*8.0, 7.0 + I*8.0, 7.0 + I*8.0, 7.0 + I*8.0};

  //Cycle through the Stokes parameters
  user_precision_t *flux_Is = calloc(num_gains, sizeof(user_precision_t));
  user_precision_t *flux_Qs = calloc(num_gains, sizeof(user_precision_t));
  user_precision_t *flux_Us = calloc(num_gains, sizeof(user_precision_t));
  user_precision_t *flux_Vs = calloc(num_gains, sizeof(user_precision_t));

  flux_Is[0] = flux_Is[4] = flux_Is[8] = flux_Is[12] = 1.0;
  flux_Qs[1] = flux_Qs[5] = flux_Qs[9] = flux_Qs[13] = 1.0;
  flux_Us[2] = flux_Us[6] = flux_Us[10] = flux_Us[14] = 1.0;
  flux_Vs[3] = flux_Vs[7] = flux_Vs[11] = flux_Vs[15] = 1.0;

  //Just keep the visibility component equal to 1.0 + 0.0j for all entries
  user_precision_complex_t *visi_components = malloc(num_gains*sizeof(user_precision_complex_t));

  for (size_t gain = 0; gain < num_gains; gain++) {
    visi_components[gain] = 1.0 + I*0.0;
  }

  //Space for outputs
  user_precision_complex_t *visi_XXs = malloc(num_gains*sizeof(user_precision_complex_t));
  user_precision_complex_t *visi_XYs = malloc(num_gains*sizeof(user_precision_complex_t));
  user_precision_complex_t *visi_YXs = malloc(num_gains*sizeof(user_precision_complex_t));
  user_precision_complex_t *visi_YYs = malloc(num_gains*sizeof(user_precision_complex_t));

  if (do_gpu == 1) {
    test_kern_apply_beam_gains(num_gains, g1xs, D1xs, D1ys, g1ys,
                                        g1xs, D1xs, D1ys, g1ys,
                                        flux_Is, flux_Qs,
                                        flux_Us, flux_Vs,
                                        visi_components,
                                        visi_XXs, visi_XYs,
                                        visi_YXs, visi_YYs,
                                        off_cardinal, full_pol);
  } else {
    apply_beam_gains_arrays_cpu(num_gains, g1xs, D1xs, D1ys, g1ys,
                                          g1xs, D1xs, D1ys, g1ys,
                                          flux_Is, flux_Qs,
                                          flux_Us, flux_Vs,
                                          visi_components,
                                          visi_XXs, visi_XYs,
                                          visi_YXs, visi_YYs,
                                          off_cardinal, full_pol);
  }

  //Expected outputs
  user_precision_t expec_XX_re_on[] =  {1.0, 1.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 8.0, 0.0, 8.0, 0.0, 30.0, -20.0, 22.0, -4.0};
  user_precision_t expec_XX_im_on[] =  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  user_precision_t expec_XY_re_on[] =  {0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 8.0, 0.0, 8.0, 0.0, 70.0, -36.0, 62.0, -4.0};
  user_precision_t expec_XY_im_on[] =  {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 8.0, -16.0};
  user_precision_t expec_YX_re_on[] =  {0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 8.0, 0.0, 8.0, 0.0, 70.0, -36.0, 62.0, -4.0};
  user_precision_t expec_YX_im_on[] =  {0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, -8.0, 0.0, -8.0, 16.0};
  user_precision_t expec_YY_re_on[] =  {1.0, -1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 8.0, 0.0, 8.0, 0.0, 174.0, -52.0, 166.0, -4.0};
  user_precision_t expec_YY_im_on[] =  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  user_precision_t expec_XX_re_off[] = {1.000, 0.000, 1.000, 0.000, 1.000, 0.000, 1.000, 0.000, 8.000, -8.000, 8.000, 0.000, 30.000, -22.000, 30.000, -4.000};
  user_precision_t expec_XX_im_off[] = {0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000};
  user_precision_t expec_XY_re_off[] = {0.000, -1.000, 0.000, 0.000, 0.000, -1.000, 0.000, 0.000, 8.000, -8.000, 0.000, 0.000, 70.000, -62.000, -36.000, -4.000};
  user_precision_t expec_XY_im_off[] = {0.000, 0.000, 0.000, 1.000, 0.000, 0.000, 0.000, -1.000, 0.000, 0.000, 0.000, 0.000, 8.000, -8.000, 0.000, -16.000};
  user_precision_t expec_YX_re_off[] = {0.000, -1.000, 0.000, 0.000, 0.000, -1.000, 0.000, 0.000, 8.000, -8.000, 0.000, 0.000, 70.000, -62.000, -36.000, -4.000};
  user_precision_t expec_YX_im_off[] = {0.000, 0.000, 0.000, -1.000, 0.000, 0.000, 0.000, 1.000, 0.000, 0.000, 0.000, 0.000, -8.000, 8.000, 0.000, 16.000};
  user_precision_t expec_YY_re_off[] = {1.000, 0.000, -1.000, 0.000, 1.000, 0.000, 1.000, 0.000, 8.000, -8.000, 0.000, 0.000, 174.000, -166.000, -52.000, -4.000};
  user_precision_t expec_YY_im_off[] = {0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000};

  // if (full_pol == 0 && off_cardinal == 0) {
  //   for (size_t gain = 0; gain < num_gains; gain++) {
  //   printf("%d %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n",
  //           gain,
  //           creal(visi_XXs[gain]),cimag(visi_XXs[gain]),
  //           creal(visi_XYs[gain]),cimag(visi_XYs[gain]),
  //           creal(visi_YXs[gain]),cimag(visi_YXs[gain]),
  //           creal(visi_YYs[gain]),cimag(visi_YYs[gain]) );
  //   }
  // }

  //OK, so only gains 0, 4, 8, 12 have Stokes I values. So if we're not doing full pol,
  //we only check those, and check everything else is set to zero
  //As we've only put through Stokes I at 0, 4, 8, 12, the values are the same
  //if in full pol mode or not

  int do_gain_test = 0;

  //Check outputs match expectations
  for (int gain = 0; gain < num_gains; gain++) {

    if (gain == 0 || gain == 4 || gain == 8 || gain == 12) {
      do_gain_test = 1;
    } else {
      if (full_pol == 1) {
        do_gain_test = 1;
      } else {
        do_gain_test = 0;
      }
    }

    if (do_gain_test == 1) {
      if (off_cardinal == 1) {
        // printf("off %.3f %.3f\n", expec_YY_re_off[gain], creal(visi_YYs[gain]));
        TEST_ASSERT_EQUAL_FLOAT(expec_XX_re_off[gain], creal(visi_XXs[gain]) );
        TEST_ASSERT_EQUAL_FLOAT(expec_XX_im_off[gain], cimag(visi_XXs[gain]) );
        TEST_ASSERT_EQUAL_FLOAT(expec_XY_re_off[gain], creal(visi_XYs[gain]) );
        TEST_ASSERT_EQUAL_FLOAT(expec_XY_im_off[gain], cimag(visi_XYs[gain]) );
        TEST_ASSERT_EQUAL_FLOAT(expec_YX_re_off[gain], creal(visi_YXs[gain]) );
        TEST_ASSERT_EQUAL_FLOAT(expec_YX_im_off[gain], cimag(visi_YXs[gain]) );
        TEST_ASSERT_EQUAL_FLOAT(expec_YY_re_off[gain], creal(visi_YYs[gain]) );
        TEST_ASSERT_EQUAL_FLOAT(expec_YY_im_off[gain], cimag(visi_YYs[gain]) );
      } else{
        TEST_ASSERT_EQUAL_FLOAT(expec_XX_re_on[gain], creal(visi_XXs[gain]) );
        TEST_ASSERT_EQUAL_FLOAT(expec_XX_im_on[gain], cimag(visi_XXs[gain]) );
        TEST_ASSERT_EQUAL_FLOAT(expec_XY_re_on[gain], creal(visi_XYs[gain]) );
        TEST_ASSERT_EQUAL_FLOAT(expec_XY_im_on[gain], cimag(visi_XYs[gain]) );
        TEST_ASSERT_EQUAL_FLOAT(expec_YX_re_on[gain], creal(visi_YXs[gain]) );
        TEST_ASSERT_EQUAL_FLOAT(expec_YX_im_on[gain], cimag(visi_YXs[gain]) );
        TEST_ASSERT_EQUAL_FLOAT(expec_YY_re_on[gain], creal(visi_YYs[gain]) );
        TEST_ASSERT_EQUAL_FLOAT(expec_YY_im_on[gain], cimag(visi_YYs[gain]) );
      }
    } else {
      TEST_ASSERT_EQUAL_FLOAT(0.0, creal(visi_XXs[gain]) );
      TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(visi_XXs[gain]) );
      TEST_ASSERT_EQUAL_FLOAT(0.0, creal(visi_XYs[gain]) );
      TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(visi_XYs[gain]) );
      TEST_ASSERT_EQUAL_FLOAT(0.0, creal(visi_YXs[gain]) );
      TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(visi_YXs[gain]) );
      TEST_ASSERT_EQUAL_FLOAT(0.0, creal(visi_YYs[gain]) );
      TEST_ASSERT_EQUAL_FLOAT(0.0, cimag(visi_YYs[gain]) );
    }
  }

  //Be free my beauties
  free(flux_Is);
  free(flux_Qs);
  free(flux_Us);
  free(flux_Vs);
  free(visi_XXs);
  free(visi_XYs);
  free(visi_YXs);
  free(visi_YYs);

}