/*
Tests that the kernel that calculates auto-correlations is doing it's job
*/
#include <math.h>
#include <unity.h>
#include <stdlib.h>
#include <complex.h>

#include "constants.h"
#include "woden_struct_defs.h"
// #include "shapelet_basis.h"
// #include "woden_settings.h"
#include "visibility_set.h"


#ifdef DOUBLE_PRECISION
  #ifdef __HIPCC__
    double TOL = 1e-11;
  #else
    double TOL = 1e-12;
  #endif
#else
  double TOL = 1e-2;
#endif

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */


//External CUDA code being linked in
extern void test_kern_calc_autos(components_t *components, int beamtype,
                                 int num_components, int num_baselines,
                                 int num_freqs, int num_times, int num_ants,
                                 int num_beams,
                                 visibility_set_t *visibility_set);


/*
For a given beam model, calculate the auto-correlations
*/
void test_calculate_autos(e_beamtype beamtype) {

  int num_comps = 4;
  int num_times = 2;
  int num_freqs = 3;

  int num_ants = 3;
  int num_baselines = (num_ants*(num_ants - 1)) / 2;

  int num_cross = num_baselines*num_times*num_freqs;
  int num_autos = num_ants*num_times*num_freqs;
  int num_visis = num_cross + num_autos;

  int num_pb_values = num_freqs*num_times*num_comps;

  components_t *components = malloc(sizeof(components_t));

  //need to fill in the following:
  components->extrap_stokesI =  malloc(num_comps*num_freqs*sizeof(user_precision_t));
  components->extrap_stokesQ =  malloc(num_comps*num_freqs*sizeof(user_precision_t));
  components->extrap_stokesU =  malloc(num_comps*num_freqs*sizeof(user_precision_t));
  components->extrap_stokesV =  malloc(num_comps*num_freqs*sizeof(user_precision_t));
  components->gxs =  malloc(num_pb_values*sizeof(user_precision_complex_t));
  components->Dxs =  malloc(num_pb_values*sizeof(user_precision_complex_t));
  components->Dys =  malloc(num_pb_values*sizeof(user_precision_complex_t));
  components->gys =  malloc(num_pb_values*sizeof(user_precision_complex_t));

  user_precision_complex_t *expec_XXs = malloc(num_autos*sizeof(user_precision_complex_t));
  user_precision_complex_t *expec_XYs = malloc(num_autos*sizeof(user_precision_complex_t));
  user_precision_complex_t *expec_YXs = malloc(num_autos*sizeof(user_precision_complex_t));
  user_precision_complex_t *expec_YYs = malloc(num_autos*sizeof(user_precision_complex_t));

  user_precision_complex_t expec_XX, expec_XY, expec_YX, expec_YY;

  int flux_ind = 0;
  int beam_ind = 0;
  double flux_value = 0.1;
  double beam_value = 0;
  int auto_stripe;

  //Flux density of componets don't change with time, just frequency
  for (int comp = 0; comp < num_comps; comp++) {
    for (int freq = 0; freq < num_freqs; freq++) {

      flux_ind = num_freqs*comp + freq;

      // expec_XX = 0.0 + 0.0*I;

      components->extrap_stokesI[flux_ind] = flux_value;
      components->extrap_stokesQ[flux_ind] = flux_value;
      components->extrap_stokesU[flux_ind] = flux_value;
      components->extrap_stokesV[flux_ind] = flux_value;

      flux_value += 0.1;

    }
  }

  // for (int i = 0; i < num_comps*num_freqs; i++) {
  //   printf("%.1f\n",components->extrap_stokesI[i] );
  // }

  user_precision_complex_t expec_flux;

  user_precision_complex_t gx, Dx, Dy, gy;

  //Beam values also change with time as well as freq and direction
  for (int time = 0; time < num_times; time++) {
    for (int freq = 0; freq < num_freqs; freq++) {

      expec_XX = 0.0 + 0.0*I;
      expec_XY = 0.0 + 0.0*I;
      expec_YX = 0.0 + 0.0*I;
      expec_YY = 0.0 + 0.0*I;

      for (int comp = 0; comp < num_comps; comp++) {

        gx = beam_value + 0.2 + 0.2*I;
        Dx = beam_value + 0.4 + 0.4*I;
        Dy = beam_value + 0.6 + 0.6*I;
        gy = beam_value + 0.8 + 0.8*I;

        components->gxs[beam_ind] = gx;
        components->Dxs[beam_ind] = Dx;
        components->Dys[beam_ind] = Dy;
        components->gys[beam_ind] = gy;

        flux_ind = num_freqs*comp + freq;

        expec_flux = components->extrap_stokesI[flux_ind] + 0*I;

        if (beamtype == NO_BEAM) {
          gx = 1 + 0.0*I;
          Dx = 0 + 0.0*I;
          Dy = 0 + 0.0*I;
          gy = 1 + 0.0*I;
        }

        //Only MWA models have leakge terms at the moment
        else if (beamtype == ANALY_DIPOLE || beamtype == GAUSS_BEAM) {
            Dx = 0 + 0.0*I;
            Dy = 0 + 0.0*I;
        }

        expec_XX += (gx*conj(gx) + Dx*conj(Dx))*expec_flux;
        expec_XX += (gx*conj(gx) - Dx*conj(Dx))*expec_flux;
        expec_XX += (gx*conj(Dx) + Dx*conj(gx))*expec_flux;
        expec_XX += ((0.0 + I*1.0)*expec_flux)*(gx*conj(Dx) - Dx*conj(gx));

        expec_XY += (gx*conj(Dy) + Dx*conj(gy))*expec_flux;
        expec_XY += (gx*conj(Dy) - Dx*conj(gy))*expec_flux;
        expec_XY += (gx*conj(gy) + Dx*conj(Dy))*expec_flux;
        expec_XY += ((0.0 + I*1.0)*expec_flux)* (gx*conj(gy) - Dx*conj(Dy));

        expec_YX += (Dy*conj(gx) + gy*conj(Dx))*expec_flux;
        expec_YX += (Dy*conj(gx) - gy*conj(Dx))*expec_flux;
        expec_YX += (Dy*conj(Dx) + gy*conj(gx))*expec_flux;
        expec_YX += ((0.0 + I*1.0)*expec_flux)* (Dy*conj(Dx) - gy*conj(gx));

        expec_YY += (Dy*conj(Dy) + gy*conj(gy))*expec_flux;
        expec_YY += (Dy*conj(Dy) - gy*conj(gy))*expec_flux;
        expec_YY += (Dy*conj(gy) + gy*conj(Dy))*expec_flux;
        expec_YY += ((0.0 + I*1.0)*expec_flux)* (Dy*conj(gy) - gy*conj(Dy));

        // printf("expec_XX %.1f %.1f\n",creal(expec_XX), cimag(expec_XX) );

        beam_value ++;
        beam_ind ++;

      }

      for (int ant = 0; ant < num_ants; ant++) {
        auto_stripe = num_ants*num_freqs*time + num_ants*freq + ant;
        expec_XXs[auto_stripe] = expec_XX;
        expec_XYs[auto_stripe] = expec_XY;
        expec_YXs[auto_stripe] = expec_YX;
        expec_YYs[auto_stripe] = expec_YY;
      }

    }
  }

  //Setup chunk_visibility_set to hold the visibility outputs of each
  visibility_set_t *visibility_set = setup_visibility_set(num_visis);

  //Set it all to zero to start
  for (int visi = 0; visi < num_visis; visi++) {
    visibility_set->sum_visi_XX_real[visi] = 0;
    visibility_set->sum_visi_XX_imag[visi] = 0;
    visibility_set->sum_visi_XY_real[visi] = 0;
    visibility_set->sum_visi_XY_imag[visi] = 0;
    visibility_set->sum_visi_YX_real[visi] = 0;
    visibility_set->sum_visi_YX_imag[visi] = 0;
    visibility_set->sum_visi_YY_imag[visi] = 0;
    visibility_set->sum_visi_YY_real[visi] = 0;
  }

  //Doing all testing with identical primary beams so set num_beams to 1
  int num_beams = 1;
  components->do_QUV = 1;
  test_kern_calc_autos(components, beamtype,
                       num_comps, num_baselines,
                       num_freqs, num_times, num_ants, num_beams,
                       visibility_set);

  //Function shouldn't touch the cross-correlation part of the output visis
  //so check that everything is zero
  for (int cross = 0; cross < num_cross; cross++) {


    TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, visibility_set->sum_visi_XX_real[cross]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, visibility_set->sum_visi_XX_imag[cross]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, visibility_set->sum_visi_XY_real[cross]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, visibility_set->sum_visi_XY_imag[cross]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, visibility_set->sum_visi_YX_real[cross]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, visibility_set->sum_visi_YX_imag[cross]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, visibility_set->sum_visi_YY_imag[cross]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, visibility_set->sum_visi_YY_real[cross]);
  }

  for (int autos = 0; autos < num_autos; autos++) {

    // printf("XX expec, calc %.1f %.1f\n",creal(expec_XXs[autos]), visibility_set->sum_visi_XX_real[num_cross + autos] );
    // printf("XY expec, calc %.1f %.1f\n",creal(expec_XYs[autos]), visibility_set->sum_visi_XY_real[num_cross + autos] );
    // printf("YX expec, calc %.1f %.1f\n",creal(expec_YXs[autos]), visibility_set->sum_visi_YX_real[num_cross + autos] );
    // printf("YY expec, calc %.1f %.1f\n",creal(expec_YYs[autos]), visibility_set->sum_visi_YY_real[num_cross + autos] );
  //
    TEST_ASSERT_DOUBLE_WITHIN(TOL, creal(expec_XXs[autos]),
                      visibility_set->sum_visi_XX_real[num_cross + autos]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, cimag(expec_XXs[autos]),
                      visibility_set->sum_visi_XX_imag[num_cross + autos]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, creal(expec_XYs[autos]),
                      visibility_set->sum_visi_XY_real[num_cross + autos]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, cimag(expec_XYs[autos]),
                      visibility_set->sum_visi_XY_imag[num_cross + autos]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, creal(expec_YXs[autos]),
                      visibility_set->sum_visi_YX_real[num_cross + autos]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, cimag(expec_YXs[autos]),
                      visibility_set->sum_visi_YX_imag[num_cross + autos]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, creal(expec_YYs[autos]),
                      visibility_set->sum_visi_YY_real[num_cross + autos]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, cimag(expec_YYs[autos]),
                      visibility_set->sum_visi_YY_imag[num_cross + autos]);
 }

  free(visibility_set->sum_visi_XX_real);
  free(visibility_set->sum_visi_XX_imag);
  free(visibility_set->sum_visi_XY_real);
  free(visibility_set->sum_visi_XY_imag);
  free(visibility_set->sum_visi_YX_real);
  free(visibility_set->sum_visi_YX_imag);
  free(visibility_set->sum_visi_YY_real);
  free(visibility_set->sum_visi_YY_imag);

  free(components->extrap_stokesI);
  free(components->extrap_stokesQ);
  free(components->extrap_stokesU);
  free(components->extrap_stokesV);
  free(components->gxs);
  free(components->Dxs);
  free(components->Dys);
  free(components->gys);

  free(expec_XXs);

}

void test_calculate_autos_NoBeam(void) {
  test_calculate_autos(NO_BEAM);
}

void test_calculate_autos_GaussBeam(void) {
  test_calculate_autos(GAUSS_BEAM);
}

void test_calculate_autos_EDA2Beam(void) {
  test_calculate_autos(ANALY_DIPOLE);
}

void test_calculate_autos_MWAAnaly(void) {
  test_calculate_autos(ANALY_DIPOLE);
}

void test_calculate_autos_MWAFEE(void) {
  test_calculate_autos(FEE_BEAM);
}

void test_calculate_autos_MWAFEEInterp(void) {
  test_calculate_autos(FEE_BEAM_INTERP);
}


//Run the test with unity
int main(void)
{
    UNITY_BEGIN();
    //Test with a single SOURCE, single COMPONENT

    RUN_TEST(test_calculate_autos_NoBeam);
    RUN_TEST(test_calculate_autos_GaussBeam);
    RUN_TEST(test_calculate_autos_EDA2Beam);
    RUN_TEST(test_calculate_autos_MWAAnaly);
    RUN_TEST(test_calculate_autos_MWAFEE);
    RUN_TEST(test_calculate_autos_MWAFEEInterp);

    return UNITY_END();
}
