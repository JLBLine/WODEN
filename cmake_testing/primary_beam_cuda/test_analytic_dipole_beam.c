#include <unity.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "constants.h"
#include "woden_precision_defs.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

//External CUDA code we're linking in
extern void test_analytic_dipole_beam(int num_components,
             int num_time_steps, int num_freqs,
             user_precision_t *azs, user_precision_t *zas, double *freqs,
             user_precision_complex_t *analy_beam_X,
             user_precision_complex_t *analy_beam_Y);


#ifdef DOUBLE_PRECISION
  double TOL = 1e-12;
#else
  double TOL = 1e-6;
#endif

/*
Test that the analytic dipole beam code returns the correct values, in the
correct order, for two time and frequency steps, with 25 directions on the sky
*/
void test_analytic_dipole_beam_GivesCorrectValues(void) {
  int num_freqs = 2;
  int num_times = 2;
  int num_az = 5;
  int num_za = 5;

  int num_components = num_az*num_za;
  int num_azza = num_az*num_za*num_times;
  int num_beam_values = num_azza*num_freqs;

  user_precision_t az_inc = 2*M_PI / num_az;
  user_precision_t za_inc = (M_PI/2) / num_za;

  user_precision_t *azs = malloc(num_azza*sizeof(user_precision_t));
  user_precision_t *zas = malloc(num_azza*sizeof(user_precision_t));

  //Function expects all az/za to include all time steps, and all
  //components. For each component, put az/za as they change with time
  //Here, we'll just stick each component to stay stationary on the sky,
  //so input the same az/za twice, one for each time step
  int azza_ind = 0;
  for (int az_ind = 0; az_ind < num_az; az_ind++) {
    for (int za_ind = 0; za_ind < num_za; za_ind++) {
      azs[azza_ind] = 0 + az_ind*az_inc;
      zas[azza_ind] = 0 + za_ind*za_inc;
      azza_ind += 1;

      azs[azza_ind] = 0 + az_ind*az_inc;
      zas[azza_ind] = 0 + za_ind*za_inc;
      azza_ind += 1;
    }
  }

  //Just test two frequencies
  double freqs[2] = {100e+6, 200e+6};

  user_precision_complex_t *analy_beam_X = malloc(num_beam_values*sizeof(user_precision_complex_t));
  user_precision_complex_t *analy_beam_Y = malloc(num_beam_values*sizeof(user_precision_complex_t));

  //Run the CUDA code
  test_analytic_dipole_beam(num_components,
               num_times, num_freqs,
               azs, zas, freqs,
               analy_beam_X, analy_beam_Y);

  //Expected real values for the X (north-south) dipole
  double expected_X[100] = { 1.0000000000000, 0.9103665168195,
            0.6699218244816, 0.3610003044827,
            0.1014449975159, 1.0000000000000, 0.9528417601434, 0.8142947351804,
            0.5946679343282, 0.3137854227599, 1.0000000000000, 0.9268204092439,
            0.7284525940295, 0.4643463853456, 0.2096881223626, 1.0000000000000,
            0.9268204092439, 0.7284525940295, 0.4643463853456, 0.2096881223626,
            1.0000000000000, 0.9528417601434, 0.8142947351804, 0.5946679343282,
            0.3137854227599, 1.0000000000000, 0.9303050339136, 0.7234573882036,
            0.4162245584976, 0.1230724172312, 1.0000000000000, 0.9737105546032,
            0.8793675930732, 0.6856376444698, 0.3806824527246, 1.0000000000000,
            0.9471192935190, 0.7866655359594, 0.5353800726217, 0.2543922787302,
            1.0000000000000, 0.9471192935190, 0.7866655359594, 0.5353800726217,
            0.2543922787302, 1.0000000000000, 0.9737105546032, 0.8793675930732,
            0.6856376444698, 0.3806824527246, 1.0000000000000, 0.9103665168195,
            0.6699218244816, 0.3610003044827, 0.1014449975159, 1.0000000000000,
            0.9528417601434, 0.8142947351804, 0.5946679343282, 0.3137854227599,
            1.0000000000000, 0.9268204092439, 0.7284525940295, 0.4643463853456,
            0.2096881223626, 1.0000000000000, 0.9268204092439, 0.7284525940295,
            0.4643463853456, 0.2096881223626, 1.0000000000000, 0.9528417601434,
            0.8142947351804, 0.5946679343282, 0.3137854227599, 1.0000000000000,
            0.9303050339136, 0.7234573882036, 0.4162245584976, 0.1230724172312,
            1.0000000000000, 0.9737105546032, 0.8793675930732, 0.6856376444698,
            0.3806824527246, 1.0000000000000, 0.9471192935190, 0.7866655359594,
            0.5353800726217, 0.2543922787302, 1.0000000000000, 0.9471192935190,
            0.7866655359594, 0.5353800726217, 0.2543922787302, 1.0000000000000,
            0.9737105546032, 0.8793675930732, 0.6856376444698, 0.3806824527246
 };

  //Expected real values for the Y (east-west) dipole
  double expected_Y[100] = { 1.0000000000000, 0.9572160026471,
            0.8280689146699, 0.6141704016472,
            0.3282829079387, 1.0000000000000, 0.9149438511715, 0.6865984726292,
            0.3922965077309, 0.1399980843763, 1.0000000000000, 0.9412935778379,
            0.7770779861836, 0.5402859766018, 0.2721978076798, 1.0000000000000,
            0.9412935778379, 0.7770779861836, 0.5402859766018, 0.2721978076798,
            1.0000000000000, 0.9149438511715, 0.6865984726292, 0.3922965077309,
            0.1399980843763, 1.0000000000000, 0.9781806001789, 0.8942425106441,
            0.7081235142838, 0.3982707083155, 1.0000000000000, 0.9349826193816,
            0.7414667198480, 0.4523083185882, 0.1698447737577, 1.0000000000000,
            0.9619094481995, 0.8391767364049, 0.6229365717453, 0.3302286261179,
            1.0000000000000, 0.9619094481995, 0.8391767364049, 0.6229365717453,
            0.3302286261179, 1.0000000000000, 0.9349826193816, 0.7414667198480,
            0.4523083185882, 0.1698447737577, 1.0000000000000, 0.9572160026471,
            0.8280689146699, 0.6141704016472, 0.3282829079387, 1.0000000000000,
            0.9149438511715, 0.6865984726292, 0.3922965077309, 0.1399980843763,
            1.0000000000000, 0.9412935778379, 0.7770779861836, 0.5402859766018,
            0.2721978076798, 1.0000000000000, 0.9412935778379, 0.7770779861836,
            0.5402859766018, 0.2721978076798, 1.0000000000000, 0.9149438511715,
            0.6865984726292, 0.3922965077309, 0.1399980843763, 1.0000000000000,
            0.9781806001789, 0.8942425106441, 0.7081235142838, 0.3982707083155,
            1.0000000000000, 0.9349826193816, 0.7414667198480, 0.4523083185882,
            0.1698447737577, 1.0000000000000, 0.9619094481995, 0.8391767364049,
            0.6229365717453, 0.3302286261179, 1.0000000000000, 0.9619094481995,
            0.8391767364049, 0.6229365717453, 0.3302286261179, 1.0000000000000,
            0.9349826193816, 0.7414667198480, 0.4523083185882, 0.1698447737577 };

  //Separate real and imaginary components
  user_precision_t *for_testing_X_re = malloc(num_beam_values*sizeof(user_precision_t));
  user_precision_t *for_testing_Y_re = malloc(num_beam_values*sizeof(user_precision_t));
  user_precision_t *for_testing_X_im = malloc(num_beam_values*sizeof(user_precision_t));
  user_precision_t *for_testing_Y_im = malloc(num_beam_values*sizeof(user_precision_t));

  for (int i = 0; i < num_beam_values; i++) {
    for_testing_X_re[i] = creal(analy_beam_X[i]);
    for_testing_X_im[i] = cimag(analy_beam_X[i]);
    for_testing_Y_re[i] = creal(analy_beam_Y[i]);
    for_testing_Y_im[i] = cimag(analy_beam_Y[i]);
  }

  //Test the arrays make sense. Analytic dipole beam is purely real, so
  //imaginary should be zero.
  for (int ind = 0; ind < num_beam_values; ind++) {
    TEST_ASSERT_DOUBLE_WITHIN(TOL, expected_X[ind], for_testing_X_re[ind]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, for_testing_X_im[ind]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, expected_Y[ind], for_testing_Y_re[ind]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, for_testing_Y_im[ind]);
  }

  //****************************************************************************
  //UNCOMMENT TO PRINT STUFF OUT AND PLOT USING PYTHON

  // user_precision_t *for_printing_X = malloc(num_beam_values*sizeof(user_precision_t));
  // user_precision_t *for_printing_Y = malloc(num_beam_values*sizeof(user_precision_t));
  // user_precision_t *for_printing_az = malloc(num_beam_values*sizeof(user_precision_t));
  // user_precision_t *for_printing_za = malloc(num_beam_values*sizeof(user_precision_t));
  //
  // for (int iCoord = 0; iCoord < num_azza; iCoord++) {
  //   for (int freq = 0; freq < num_freqs; freq++) {
  //
  //     int component = (int)floorf((float)iCoord / (float)num_times);
  //     int time_ind = iCoord - component*num_times;
  //     int beam_ind = num_freqs*time_ind*num_components + (num_components*freq) + component;
  //
  //     for_printing_X[beam_ind] = creal(analy_beam_X[beam_ind]);
  //     for_printing_Y[beam_ind] = creal(analy_beam_Y[beam_ind]);
  //     for_printing_az[beam_ind] = azs[iCoord];
  //     for_printing_za[beam_ind] = zas[iCoord];
  //   }
  // }

  // for (int i = 0; i < num_beam_values; i++) {
  //   printf("%.13f, \n", for_printing_X[i]);
  // }
  //
  // printf("\n");
  //
  // for (int i = 0; i < num_beam_values; i++) {
  //   printf("%.13f, \n", for_printing_Y[i]);
  // }

  // for (int i = 0; i < num_beam_values; i++) {
  //   printf("%.7f %.7f %.7f 0 %.7f 0\n",for_printing_az[i], for_printing_za[i],
  //                                      for_printing_X[i], for_printing_Y[i] );
  // }

  // for (int i = 0; i < num_beam_values; i++) {
  //   printf("%.7f %.7f %.7f 0 %.7f 0\n",for_printing_az[i], for_printing_za[i],
  //                                      for_printing_X[i], for_printing_Y[i] );
  // }
  //****************************************************************************

  free(for_testing_X_re);
  free(for_testing_Y_re);
  free(for_testing_X_im);
  free(for_testing_Y_im);
}



//Run the test with unity
int main(void)
{
    UNITY_BEGIN();
    RUN_TEST(test_analytic_dipole_beam_GivesCorrectValues);
    return UNITY_END();
}
