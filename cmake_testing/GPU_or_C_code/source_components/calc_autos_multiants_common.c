/*
Tests that the kernel that calculates auto-correlations is doing it's job
*/
#include "calc_autos_multiants_common.h"


/*
For a given beam model, calculate the auto-correlations
*/
void test_calculate_autos_multiants(e_beamtype beamtype, int do_gpu) {

  #ifdef DOUBLE_PRECISION
    #ifdef __HIPCC__
      double TOL = 1e-10;
    #else
      double TOL = 1e-12;
    #endif
  #else
    double TOL = 1e-1;
  #endif

  int num_comps = 4;
  int num_times = 2;
  int num_freqs = 3;

  int num_ants = 3;
  int num_baselines = (num_ants*(num_ants - 1)) / 2;

  int num_cross = num_baselines*num_times*num_freqs;
  int num_autos = num_ants*num_times*num_freqs;
  int num_visis = num_cross + num_autos;

  int num_pb_values = num_ants*num_freqs*num_times*num_comps;

  double tile_x_gains[3] = {0.0, 0.4, 0.6};
  double tile_y_gains[3] = {0.2, 0.8, 1.0};

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
  double flux_value = 1;
  double beam_value = 0;
  int auto_stripe;

  //Flux density of componets don't change with time, just frequency
  for (int comp = 0; comp < num_comps; comp++) {
    for (int freq = 0; freq < num_freqs; freq++) {

      flux_ind = num_freqs*comp + freq;

      // expec_XX = 0.0 + 0.0*I;

      components->extrap_stokesI[flux_ind] = flux_value;
      components->extrap_stokesQ[flux_ind] = 0;
      components->extrap_stokesU[flux_ind] = 0;
      components->extrap_stokesV[flux_ind] = 0;

      flux_value ++;


    }
  }

  // for (int i = 0; i < num_comps*num_freqs; i++) {
  //   printf("%.1f\n",components->extrap_stokesI[i] );
  // }

  user_precision_complex_t expec_I;

  user_precision_complex_t gx, Dx, Dy, gy;

  //Beam values change with time as well as antenna, freq, and direction
  for (int ant = 0; ant < num_ants; ant++) {
    for (int time = 0; time < num_times; time++) {
      for (int freq = 0; freq < num_freqs; freq++) {

        expec_XX = 0.0 + 0.0*I;
        expec_XY = 0.0 + 0.0*I;
        expec_YX = 0.0 + 0.0*I;
        expec_YY = 0.0 + 0.0*I;

        for (int comp = 0; comp < num_comps; comp++) {

          gx = tile_x_gains[ant]*(beam_value + 0.2 + 0.2*I);
          Dx = tile_x_gains[ant]*(beam_value + 0.4 + 0.4*I);
          Dy = tile_y_gains[ant]*(beam_value + 0.6 + 0.6*I);
          gy = tile_y_gains[ant]*(beam_value + 0.8 + 0.8*I);

          components->gxs[beam_ind] = gx;
          components->Dxs[beam_ind] = Dx;
          components->Dys[beam_ind] = Dy;
          components->gys[beam_ind] = gy;

          flux_ind = num_freqs*comp + freq;

          expec_I = components->extrap_stokesI[flux_ind] + 0*I;

          if (beamtype == NO_BEAM) {
            expec_XX += expec_I;
            expec_YY += expec_I;
          }

          //Only MWA models have leakge terms at the moment
          else if (beamtype == FEE_BEAM || beamtype == FEE_BEAM_INTERP || beamtype == MWA_ANALY) {
            expec_XX += (gx*conj(gx) + Dx*conj(Dx))*expec_I;
            expec_XY += (gx*conj(Dy) + conj(gy)*Dx)*expec_I;
            expec_YX += (Dy*conj(gx) + gy*conj(Dx))*expec_I;
            expec_YY += (gy*conj(gy) + Dy*conj(Dy))*expec_I;
          }

          else{
              expec_XX += (gx*conj(gx))*expec_I;
              expec_YY += (gy*conj(gy))*expec_I;
          }
          beam_value ++;
          beam_ind ++;
        }

        //The extra gain that dipole amplitudes is added here
        auto_stripe = num_ants*num_freqs*time + num_ants*freq + ant;
        // printf("auto_stripe %d %.1f\n",auto_stripe, creal(expec_XX) );
        expec_XXs[auto_stripe] = expec_XX;
        expec_XYs[auto_stripe] = expec_XY;
        expec_YXs[auto_stripe] = expec_YX;
        expec_YYs[auto_stripe] = expec_YY;
      }
    }
  }

  // printf("beam_ind %d %d\n",beam_ind, num_pb_values );

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

  components->do_QUV = 0;
  if (do_gpu == 1) {
    test_kern_calc_autos(components, beamtype,
                        num_comps, num_baselines,
                        num_freqs, num_times, num_ants,
                        num_ants,
                        visibility_set);
  } else {
    test_calc_autos_cpu(components, beamtype,
                        num_comps, num_baselines,
                        num_freqs, num_times, num_ants, num_ants,
                        visibility_set);
  }

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

    // printf("ind XX expec, calc %d %.1f %.1f\n",autos, creal(expec_XXs[autos]), visibility_set->sum_visi_XX_real[num_cross + autos] );
    // printf("XY expec, calc %.1f %.1f\n",creal(expec_XYs[autos]), visibility_set->sum_visi_XY_real[num_cross + autos] );
    // printf("YX expec, calc %.1f %.1f\n",creal(expec_YXs[autos]), visibility_set->sum_visi_YX_real[num_cross + autos] );
    // printf("ind YY expec, calc %d %.1f %.1f\n",autos, creal(expec_YYs[autos]), visibility_set->sum_visi_YY_real[num_cross + autos] );
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
