#include <unity.h>
#include <stdlib.h>
#include <math.h>

// #include <module_a.h>
#include "constants.h"
// #include "fundamental_coords.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

extern void test_kern_calc_uvw(float *X_diff, float *Y_diff, float *Z_diff,
           float *u_metres, float *v_metres, float *w_metres,
           float *us, float *vs, float *ws, float *wavelengths,
           float dec0,
           float *cha0s, float *sha0s,
           int num_visis, int num_baselines);

#define UNITY_INCLUDE_FLOAT

// void setup_lmn_params()

void test_kern_calc_uvw_GivesCorrectCoords(void){
    float ra0 = 0.0*DD2R;
    float dec0 = 0.0*DD2R;

    int num_baselines = 5;
    int num_freqs = 2;
    int num_times = 2;
    int num_visis = num_baselines*num_freqs*num_times;

    float *X_diff = malloc(num_baselines*sizeof(float));
    float *Y_diff = malloc(num_baselines*sizeof(float));
    float *Z_diff = malloc(num_baselines*sizeof(float));

    for (int i = 0; i < num_baselines; i++) {
      X_diff[i] = i + 1;
      Y_diff[i] = i + 1;
      Z_diff[i] = i + 1;
    }

    float lst_base = 0.0;
    float time_res = 8.0;
    float freq_res = 50e+6;
    float base_band_freq = 150e+6;

    float lsts[num_times];

    for ( int time_step = 0; time_step < num_times; time_step++ ) {
      float lst = lst_base + time_step*time_res*SOLAR2SIDEREAL*DS2R;

      //Add half a time_res so we are sampling centre of each time step
      // lst += 0.5*time_res*SOLAR2SIDEREAL*DS2R;
      lsts[time_step] = lst;
    }

    float *wavelengths = malloc(num_visis*sizeof(float));
    float *cha0s = malloc(num_visis*sizeof(float));
    float *sha0s = malloc(num_visis*sizeof(float));

    float ha0, sha0, cha0, frequency, wavelength;

    for ( int time_step = 0; time_step < num_times; time_step++ ) {
      ha0 = lsts[time_step] - ra0;
      sha0 = sin(ha0);
      cha0 = cos(ha0);

      for (int freq_step = 0; freq_step < num_freqs; freq_step++) {
        frequency = base_band_freq + (freq_res*freq_step);
        wavelength = VELC / frequency;
        int step = num_baselines*(time_step*num_freqs + freq_step);

        for (int baseline = 0; baseline < num_baselines; baseline++) {
          cha0s[step + baseline] = cha0;
          sha0s[step + baseline] = sha0;
          wavelengths[step + baseline] = wavelength;
        }//baseline loop
      }//freq loop
    }//time loop

    //
    float *us = malloc(num_visis*sizeof(float));
    float *vs = malloc(num_visis*sizeof(float));
    float *ws = malloc(num_visis*sizeof(float));

    float *u_metres = malloc(num_visis*sizeof(float));
    float *v_metres = malloc(num_visis*sizeof(float));
    float *w_metres = malloc(num_visis*sizeof(float));


    //
    test_kern_calc_uvw(X_diff, Y_diff, Z_diff,
           u_metres, v_metres, w_metres,
           us, vs, ws, wavelengths,
           dec0, cha0s, sha0s,
           num_visis, num_baselines);

    for (int i = 0; i < num_visis; i++) {
      printf("%.5f %.5f\n", cha0s[i], sha0s[i]);
      // printf("%.5f %.5f %.5f %5f\n",X_diff[i], Y_diff[i], Z_diff[i], wavelengths[i] );
      printf("%.5f\n", wavelengths[i] );
      printf("%.5f %.5f %.5f\n",us[i], vs[i], ws[i] );
      // printf("%.5f %.5f %.5f\n",u_metres[i], v_metres[i], w_metres[i] );
      printf("-------------------------------------------------------------\n");

    }


    //
    // float l_expected[5] = {0.0, 0.5, sqrt(2)/2.0, sqrt(3)/2.0, 1.0};
    // float n_expected[5] = {1.0, sqrt(3)/2.0, sqrt(2)/2.0, 0.5, 0.0};
    //
    // TEST_ASSERT_EQUAL_FLOAT_ARRAY(l_expected, ls, num_points);
    // TEST_ASSERT_EQUAL_FLOAT_ARRAY(zeroes, ms, num_points);
    // // TEST_ASSERT_EQUAL_FLOAT_ARRAY(n_expected, ns, num_points);
    // //Here, when n should be 0.0, we get some floating point error from the
    // //sin/cos functions
    // for (int i = 0; i < num_points; i++) {
    //   TEST_ASSERT_FLOAT_WITHIN(1e-7, n_expected[i], ns[i]);
    // }


}

int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_kern_calc_uvw_GivesCorrectCoords);

    return UNITY_END();
}
