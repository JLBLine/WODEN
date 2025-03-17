#include "update_sum_visis_common.h"


void setUp (void) {} /* Is run before eVary test, put unit init calls here. */
void tearDown (void) {} /* Is run after eVary test, put unit clean-up calls here. */


void test_update_sum_visis_cpu(int num_freqs, int num_cross,
          int num_baselines, int num_components, int num_times, int beamtype,
          int use_twobeams, int num_ants, int off_cardinal_dipoles,
          user_precision_complex_t *g1xs, user_precision_complex_t *D1xs,
          user_precision_complex_t *D1ys, user_precision_complex_t *g1ys,
          user_precision_complex_t *visi_components,
          user_precision_t *flux_I, user_precision_t *flux_Q,
          user_precision_t *flux_U, user_precision_t *flux_V,
          user_precision_t *sum_visi_XX_real, user_precision_t *sum_visi_XX_imag,
          user_precision_t *sum_visi_XY_real, user_precision_t *sum_visi_XY_imag,
          user_precision_t *sum_visi_YX_real, user_precision_t *sum_visi_YX_imag,
          user_precision_t *sum_visi_YY_real, user_precision_t *sum_visi_YY_imag){

  int *ant1_to_baseline_map = malloc(num_baselines*sizeof(int));
  int *ant2_to_baseline_map = malloc(num_baselines*sizeof(int));

  fill_ant_to_baseline_mapping_cpu(num_ants, ant1_to_baseline_map, ant2_to_baseline_map);

  for (int iBaseline = 0; iBaseline < num_cross; iBaseline++) {
    int time_ind = (int)floorf( (user_precision_t)iBaseline / ((user_precision_t)num_baselines * (user_precision_t)num_freqs));
    int freq_ind = (int)floorf( ((user_precision_t)iBaseline - ((user_precision_t)time_ind*(user_precision_t)num_baselines * (user_precision_t)num_freqs)) / (user_precision_t)num_baselines);

    for (int iComponent = 0; iComponent < num_components; iComponent++) {

      //There is a flux for every frequnecy and component
      int flux_ind = num_components*freq_ind + iComponent;

      update_sum_visis_stokesIQUV_cpu(iBaseline, iComponent, num_freqs,
             num_baselines, num_components, num_times, beamtype, off_cardinal_dipoles,
             g1xs, D1xs, D1ys, g1ys,
             ant1_to_baseline_map, ant2_to_baseline_map, use_twobeams,
             visi_components[iBaseline],
             flux_I[flux_ind], flux_Q[flux_ind],
             flux_U[flux_ind], flux_V[flux_ind],
             sum_visi_XX_real, sum_visi_XX_imag,
             sum_visi_XY_real, sum_visi_XY_imag,
             sum_visi_YX_real, sum_visi_YX_imag,
             sum_visi_YY_real, sum_visi_YY_imag);

    }
  }
  free(ant1_to_baseline_map);
  free(ant2_to_baseline_map);
}

/*
Test the __device__ code that updates the summed visibilities by grabbing the
correct beam gain and mesurement equation, multiplying and summing onto the visi
Here we keep the component visibilities and fluxes constant and vary the beam gains
Test works for all primary beam types
*/
void test_update_sum_visis_VaryGainChooseBeams(int beamtype, int do_gpu) {

  int num_baselines = 3.0;
  int num_times = 4.0;
  int num_freqs = 3.0;

  int num_visis = num_baselines*num_times*num_freqs;

  int num_components = 10.0;

  user_precision_complex_t *primay_beam_J00 = malloc(num_freqs*num_times*num_components*sizeof(user_precision_complex_t));
  user_precision_complex_t *primay_beam_J01 = malloc(num_freqs*num_times*num_components*sizeof(user_precision_complex_t));
  user_precision_complex_t *primay_beam_J10 = malloc(num_freqs*num_times*num_components*sizeof(user_precision_complex_t));
  user_precision_complex_t *primay_beam_J11 = malloc(num_freqs*num_times*num_components*sizeof(user_precision_complex_t));
  user_precision_complex_t *visi_components = malloc(num_visis*sizeof(user_precision_complex_t));

  //Just stick base measurement equation to 1.0
  for (int visi = 0; visi < num_visis; visi++) {
    visi_components[visi] = 1.0 + I*0.0;
  }

  //Fill in some varying beam gains based on beam type
  int count = 1;
  if (beamtype != NO_BEAM) {
    for (int visi = 0; visi < num_components*num_times*num_freqs; visi++) {
      primay_beam_J00[visi] = count + I*0.0;
      primay_beam_J11[visi] = count + I*0.0;
      count ++ ;
    }
  }

  //Only FEE_BEAM has cross-pols
  count = 1;
  if (beamtype == FEE_BEAM || beamtype == FEE_BEAM_INTERP || beamtype == MWA_ANALY) {
    for (int visi = 0; visi < num_components*num_times*num_freqs; visi++) {
      primay_beam_J01[visi] = count + I*0.0;
      primay_beam_J10[visi] = count + I*0.0;
      count ++ ;
    }
  }

  //Set all the Stokes I to 1.0 and other stokes to zero
  user_precision_t *flux_I = malloc(num_freqs*num_components*sizeof(user_precision_t));
  user_precision_t *flux_Q = calloc(num_freqs*num_components, sizeof(user_precision_t));
  user_precision_t *flux_U = calloc(num_freqs*num_components, sizeof(user_precision_t));
  user_precision_t *flux_V = calloc(num_freqs*num_components, sizeof(user_precision_t));

  //Just stick base measurement equation to 1.0
  for (int visi = 0; visi < num_freqs*num_components; visi++) {
    flux_I[visi] = 1.0;
  }

  //Make sure arrays to hold summed visis are initialised to zero
  user_precision_t *sum_visi_XX_real = calloc(num_visis, sizeof(user_precision_t));
  user_precision_t *sum_visi_XX_imag = calloc(num_visis, sizeof(user_precision_t));
  user_precision_t *sum_visi_XY_real = calloc(num_visis, sizeof(user_precision_t));
  user_precision_t *sum_visi_XY_imag = calloc(num_visis, sizeof(user_precision_t));
  user_precision_t *sum_visi_YX_real = calloc(num_visis, sizeof(user_precision_t));
  user_precision_t *sum_visi_YX_imag = calloc(num_visis, sizeof(user_precision_t));
  user_precision_t *sum_visi_YY_real = calloc(num_visis, sizeof(user_precision_t));
  user_precision_t *sum_visi_YY_imag = calloc(num_visis, sizeof(user_precision_t));

  int num_ants = 1;
  int use_twoants = 0;
  int off_cardinal_dipoles = 0;

  if (do_gpu == 1) {
    test_kern_update_sum_visis(num_freqs, num_visis,
          num_baselines, num_components, num_times, beamtype,
          use_twoants, num_ants, off_cardinal_dipoles,
          primay_beam_J00, primay_beam_J01,
          primay_beam_J10, primay_beam_J11,
          visi_components,
          flux_I, flux_Q, flux_U, flux_V,
          sum_visi_XX_real, sum_visi_XX_imag,
          sum_visi_XY_real, sum_visi_XY_imag,
          sum_visi_YX_real, sum_visi_YX_imag,
          sum_visi_YY_real, sum_visi_YY_imag);
  } else {
    test_update_sum_visis_cpu(num_freqs, num_visis,
          num_baselines, num_components, num_times, beamtype,
          use_twoants, num_ants, off_cardinal_dipoles,
          primay_beam_J00, primay_beam_J01,
          primay_beam_J10, primay_beam_J11,
          visi_components,
          flux_I, flux_Q, flux_U, flux_V,
          sum_visi_XX_real, sum_visi_XX_imag,
          sum_visi_XY_real, sum_visi_XY_imag,
          sum_visi_YX_real, sum_visi_YX_imag,
          sum_visi_YY_real, sum_visi_YY_imag);
  }

  //Expected values here include cross-pol gains
  double expected_order[] = { 770.0, 770.0, 770.0, 4970.0, 4970.0, 4970.0,
                             13170.0, 13170.0, 13170.0, 25370.0, 25370.0,
                             25370.0, 41570.0, 41570.0, 41570.0, 61770.0,
                             61770.0, 61770.0, 85970.0, 85970.0, 85970.0,
                             114170.0, 114170.0, 114170.0, 146370.0, 146370.0,
                             146370.0, 182570.0, 182570.0, 182570.0, 222770.0,
                             222770.0, 222770.0, 266970.0, 266970.0, 266970.0};

  //No cross-pols so divide expected values by 2.0
  if (beamtype == ANALY_DIPOLE || beamtype == GAUSS_BEAM) {
    for (int output = 0; output < num_visis; output++) {

      // TEST_ASSERT_DOUBLE_WITHIN(TOL, expected_order[output]/2, sum_visi_XX_real[output]);
      // TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, sum_visi_XY_real[output]);
      // TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, sum_visi_YX_real[output]);
      // TEST_ASSERT_DOUBLE_WITHIN(TOL, expected_order[output]/2, sum_visi_YY_real[output]);
      // TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, sum_visi_XX_imag[output]);
      // TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, sum_visi_XY_imag[output]);
      // TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, sum_visi_YX_imag[output]);
      // TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, sum_visi_YY_imag[output]);

      TEST_ASSERT_EQUAL_DOUBLE(expected_order[output]/2, sum_visi_XX_real[output]);
      TEST_ASSERT_EQUAL_DOUBLE(0.0, sum_visi_XY_real[output]);
      TEST_ASSERT_EQUAL_DOUBLE(0.0, sum_visi_YX_real[output]);
      TEST_ASSERT_EQUAL_DOUBLE(expected_order[output]/2, sum_visi_YY_real[output]);
      TEST_ASSERT_EQUAL_DOUBLE(0.0, sum_visi_XX_imag[output]);
      TEST_ASSERT_EQUAL_DOUBLE(0.0, sum_visi_XY_imag[output]);
      TEST_ASSERT_EQUAL_DOUBLE(0.0, sum_visi_YX_imag[output]);
      TEST_ASSERT_EQUAL_DOUBLE(0.0, sum_visi_YY_imag[output]);

    }
  }
  else if (beamtype == FEE_BEAM || beamtype == FEE_BEAM_INTERP || beamtype == MWA_ANALY) {
    for (int output = 0; output < num_visis; output++) {

      // printf("%.1f %.1f %.1f %.1f %.1f %.1f %.1f %.1f\n",
      //     sum_visi_XX_real[output], sum_visi_XX_imag[output],
      //     sum_visi_XY_real[output], sum_visi_XY_imag[output],
      //     sum_visi_YX_real[output], sum_visi_YX_imag[output],
      //     sum_visi_YY_real[output], sum_visi_YY_imag[output]);

      TEST_ASSERT_EQUAL_DOUBLE(expected_order[output], sum_visi_XX_real[output]);
      TEST_ASSERT_EQUAL_DOUBLE(expected_order[output], sum_visi_XY_real[output]);
      TEST_ASSERT_EQUAL_DOUBLE(expected_order[output], sum_visi_YX_real[output]);
      TEST_ASSERT_EQUAL_DOUBLE(expected_order[output], sum_visi_YY_real[output]);
      TEST_ASSERT_EQUAL_DOUBLE(0.0, sum_visi_XX_imag[output]);
      TEST_ASSERT_EQUAL_DOUBLE(0.0, sum_visi_XY_imag[output]);
      TEST_ASSERT_EQUAL_DOUBLE(0.0, sum_visi_YX_imag[output]);
      TEST_ASSERT_EQUAL_DOUBLE(0.0, sum_visi_YY_imag[output]);

    }
  }
  else {
    //If there is no beam, should just end up with the number of components
    //as we apply a gain of 1.0
    for (int output = 0; output < num_visis; output++) {

      TEST_ASSERT_EQUAL_DOUBLE(num_components, sum_visi_XX_real[output]);
      TEST_ASSERT_EQUAL_DOUBLE(0.0, sum_visi_XY_real[output]);
      TEST_ASSERT_EQUAL_DOUBLE(0.0, sum_visi_YX_real[output]);
      TEST_ASSERT_EQUAL_DOUBLE(num_components, sum_visi_YY_real[output]);
      TEST_ASSERT_EQUAL_DOUBLE(0.0, sum_visi_XX_imag[output]);
      TEST_ASSERT_EQUAL_DOUBLE(0.0, sum_visi_XY_imag[output]);
      TEST_ASSERT_EQUAL_DOUBLE(0.0, sum_visi_YX_imag[output]);
      TEST_ASSERT_EQUAL_DOUBLE(0.0, sum_visi_YY_imag[output]);

    }
  }
  //
  // //Be free my beauties
  free( primay_beam_J00 );
  free( primay_beam_J01 );
  free( primay_beam_J10 );
  free( primay_beam_J11 );
  free( visi_components );
  free( flux_I );
  free( flux_Q );
  free( flux_U );
  free( flux_V );
  free( sum_visi_XX_real );
  free( sum_visi_XX_imag );
  free( sum_visi_XY_real );
  free( sum_visi_XY_imag );
  free( sum_visi_YX_real );
  free( sum_visi_YX_imag );
  free( sum_visi_YY_real );
  free( sum_visi_YY_imag );

}

/*
Test the __device__ code that updates the summed visibilities by grabbing the
correct beam gain and mesurement equation, multiplying and summing onto the visi
Here we keep the component visibilities and the beam gains constant and vary fluxes
Test works for all primary beam types
*/
void test_update_sum_visis_VaryFluxesChooseBeams(int beamtype, int do_gpu) {

  int num_baselines = 3.0;
  int num_times = 4.0;
  int num_freqs = 3.0;

  int num_visis = num_baselines*num_times*num_freqs;

  int num_components = 10.0;

  user_precision_complex_t *primay_beam_J00 = malloc(num_freqs*num_times*num_components*sizeof(user_precision_complex_t));
  user_precision_complex_t *primay_beam_J01 = malloc(num_freqs*num_times*num_components*sizeof(user_precision_complex_t));
  user_precision_complex_t *primay_beam_J10 = malloc(num_freqs*num_times*num_components*sizeof(user_precision_complex_t));
  user_precision_complex_t *primay_beam_J11 = malloc(num_freqs*num_times*num_components*sizeof(user_precision_complex_t));
  user_precision_complex_t *visi_components = malloc(num_visis*sizeof(user_precision_complex_t));

  //Just stick base measurement equation to 1.0
  for (int visi = 0; visi < num_visis; visi++) {
    visi_components[visi] = 1.0 + I*0.0;
  }

  //Set all beam gains to 1.0
  for (int visi = 0; visi < num_components*num_times*num_freqs; visi++) {
    primay_beam_J00[visi] = 1.0 + I*0.0;
    primay_beam_J01[visi] = 1.0 + I*0.0;
    primay_beam_J10[visi] = 1.0 + I*0.0;
    primay_beam_J11[visi] = 1.0 + I*0.0;
  }

  //Vary Stokes parameters
  user_precision_t *flux_I = malloc(num_freqs*num_components*sizeof(user_precision_t));
  user_precision_t *flux_Q = malloc(num_freqs*num_components*sizeof(user_precision_t));
  user_precision_t *flux_U = malloc(num_freqs*num_components*sizeof(user_precision_t));
  user_precision_t *flux_V = malloc(num_freqs*num_components*sizeof(user_precision_t));

  //Just stick base measurement equation to 1.0
  int count = 1;

  for (int visi = 0; visi < num_freqs*num_components; visi++) {
    flux_I[visi] = count;
    flux_Q[visi] = 0;
    flux_U[visi] = 0;
    flux_V[visi] = 0;
    // printf("%.3f\n",flux_I[visi] );
    count ++;
  }

  //Make sure arrays to hold summed visis are initialised to zero
  user_precision_t *sum_visi_XX_real = calloc(num_visis, sizeof(user_precision_t));
  user_precision_t *sum_visi_XX_imag = calloc(num_visis, sizeof(user_precision_t));
  user_precision_t *sum_visi_XY_real = calloc(num_visis, sizeof(user_precision_t));
  user_precision_t *sum_visi_XY_imag = calloc(num_visis, sizeof(user_precision_t));
  user_precision_t *sum_visi_YX_real = calloc(num_visis, sizeof(user_precision_t));
  user_precision_t *sum_visi_YX_imag = calloc(num_visis, sizeof(user_precision_t));
  user_precision_t *sum_visi_YY_real = calloc(num_visis, sizeof(user_precision_t));
  user_precision_t *sum_visi_YY_imag = calloc(num_visis, sizeof(user_precision_t));

  int num_ants = 1;
  int use_twoants = 0;
  int off_cardinal_dipoles = 0;

  if (do_gpu == 1) {
    test_kern_update_sum_visis(num_freqs, num_visis,
          num_baselines, num_components, num_times, beamtype,
          use_twoants, num_ants, off_cardinal_dipoles,
          primay_beam_J00, primay_beam_J01,
          primay_beam_J10, primay_beam_J11,
          visi_components,
          flux_I, flux_Q, flux_U, flux_V,
          sum_visi_XX_real, sum_visi_XX_imag,
          sum_visi_XY_real, sum_visi_XY_imag,
          sum_visi_YX_real, sum_visi_YX_imag,
          sum_visi_YY_real, sum_visi_YY_imag);
  } else {
    test_update_sum_visis_cpu(num_freqs, num_visis,
          num_baselines, num_components, num_times, beamtype,
          use_twoants, num_ants, off_cardinal_dipoles,
          primay_beam_J00, primay_beam_J01,
          primay_beam_J10, primay_beam_J11,
          visi_components,
          flux_I, flux_Q, flux_U, flux_V,
          sum_visi_XX_real, sum_visi_XX_imag,
          sum_visi_XY_real, sum_visi_XY_imag,
          sum_visi_YX_real, sum_visi_YX_imag,
          sum_visi_YY_real, sum_visi_YY_imag);
  }

  user_precision_t *expected_order = malloc(num_visis*sizeof(user_precision_t));

  //These expected values have no cross-pols in them
  user_precision_t freq_sum[] = {55.0, 155.0, 255.0};

  int ind = 0;
  for (int time_ind = 0; time_ind < num_times; time_ind++) {
    for (int freq = 0; freq < num_freqs; freq++) {
      for (int baseline = 0; baseline < num_baselines; baseline++) {
        expected_order[ind] = freq_sum[freq];
        ind ++;
      }
    }
  }

  //These beams are fully real with no leakage terms
  if (beamtype == ANALY_DIPOLE || beamtype == GAUSS_BEAM) {
    for (int output = 0; output < num_visis; output++) {

      // printf("%.1f %.1f %.1f %.1f %.1f %.1f %.1f %.1f\n",
      //     sum_visi_XX_real[output], sum_visi_XX_imag[output],
      //     sum_visi_XY_real[output], sum_visi_XY_imag[output],
      //     sum_visi_YX_real[output], sum_visi_YX_imag[output],
      //     sum_visi_YY_real[output], sum_visi_YY_imag[output]);

      TEST_ASSERT_EQUAL_FLOAT(expected_order[output], sum_visi_XX_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_XY_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_YX_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(expected_order[output], sum_visi_YY_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_XX_imag[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_XY_imag[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_YX_imag[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_YY_imag[output]);

    }
  }
  //No cross-pols in the expected values so multiply expected values by 2.0
  else if (beamtype == FEE_BEAM || beamtype == FEE_BEAM_INTERP || beamtype == MWA_ANALY) {
    for (int output = 0; output < num_visis; output++) {

      TEST_ASSERT_EQUAL_FLOAT(2*expected_order[output], sum_visi_XX_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(2*expected_order[output], sum_visi_XY_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(2*expected_order[output], sum_visi_YX_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(2*expected_order[output], sum_visi_YY_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_XX_imag[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_XY_imag[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_YX_imag[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_YY_imag[output]);

    }
  }
  else {
    //If there is no beam, should just end up with the number of components
    //as we apply a gain of 1.0
    for (int output = 0; output < num_visis; output++) {

      TEST_ASSERT_EQUAL_FLOAT(expected_order[output], sum_visi_XX_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_XY_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_YX_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(expected_order[output], sum_visi_YY_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_XX_imag[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_XY_imag[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_YX_imag[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_YY_imag[output]);

    }
  }
  //
  // //Be free my beauties
  free( primay_beam_J00 );
  free( primay_beam_J01 );
  free( primay_beam_J10 );
  free( primay_beam_J11 );
  free( visi_components );
  free( flux_I );
  free( flux_Q );
  free( flux_U );
  free( flux_V );
  free( sum_visi_XX_real );
  free( sum_visi_XX_imag );
  free( sum_visi_XY_real );
  free( sum_visi_XY_imag );
  free( sum_visi_YX_real );
  free( sum_visi_YX_imag );
  free( sum_visi_YY_real );
  free( sum_visi_YY_imag );

}

/*
Test the __device__ code that updates the summed visibilities by grabbing the
correct beam gain and mesurement equation, multiplying and summing onto the visi
Here we keep the fluxes and beam gains constant and vary the component visibilities
Test works for all primary beam types
*/
void test_update_sum_visis_VaryVisiChooseBeams(int beamtype, int do_gpu) {

  int num_baselines = 3.0;
  int num_times = 4.0;
  int num_freqs = 3.0;

  int num_visis = num_baselines*num_times*num_freqs;

  int num_components = 10.0;

  user_precision_complex_t *primay_beam_J00 = malloc(num_freqs*num_times*num_components*sizeof(user_precision_complex_t));
  user_precision_complex_t *primay_beam_J01 = malloc(num_freqs*num_times*num_components*sizeof(user_precision_complex_t));
  user_precision_complex_t *primay_beam_J10 = malloc(num_freqs*num_times*num_components*sizeof(user_precision_complex_t));
  user_precision_complex_t *primay_beam_J11 = malloc(num_freqs*num_times*num_components*sizeof(user_precision_complex_t));
  user_precision_complex_t *visi_components = malloc(num_visis*sizeof(user_precision_complex_t));

  //Vary the base visibilities
  int count = 1;
  for (int visi = 0; visi < num_visis; visi++) {
    visi_components[visi] = count + I*count;
    count ++;
  }

  //Set all beam gains to 1.0
  for (int visi = 0; visi < num_components*num_times*num_freqs; visi++) {
    primay_beam_J00[visi] = 1.0 + I*0.0;
    primay_beam_J01[visi] = 1.0 + I*0.0;
    primay_beam_J10[visi] = 1.0 + I*0.0;
    primay_beam_J11[visi] = 1.0 + I*0.0;
  }

  //Vary Stokes parameters
  user_precision_t *flux_I = malloc(num_freqs*num_components*sizeof(user_precision_t));
  user_precision_t *flux_Q = calloc(num_freqs, num_components*sizeof(user_precision_t));
  user_precision_t *flux_U = calloc(num_freqs, num_components*sizeof(user_precision_t));
  user_precision_t *flux_V = calloc(num_freqs, num_components*sizeof(user_precision_t));

  //Just stick base Stokes I to 1.0
  for (int visi = 0; visi < num_freqs*num_components; visi++) {
    flux_I[visi] = 1.0;
  }

  //Make sure arrays to hold summed visis are initialised to zero
  user_precision_t *sum_visi_XX_real = calloc(num_visis, sizeof(user_precision_t));
  user_precision_t *sum_visi_XX_imag = calloc(num_visis, sizeof(user_precision_t));
  user_precision_t *sum_visi_XY_real = calloc(num_visis, sizeof(user_precision_t));
  user_precision_t *sum_visi_XY_imag = calloc(num_visis, sizeof(user_precision_t));
  user_precision_t *sum_visi_YX_real = calloc(num_visis, sizeof(user_precision_t));
  user_precision_t *sum_visi_YX_imag = calloc(num_visis, sizeof(user_precision_t));
  user_precision_t *sum_visi_YY_real = calloc(num_visis, sizeof(user_precision_t));
  user_precision_t *sum_visi_YY_imag = calloc(num_visis, sizeof(user_precision_t));


  int num_ants = 1;
  int use_twoants = 0;
  int off_cardinal_dipoles = 0;

  if (do_gpu == 1) {
    test_kern_update_sum_visis(num_freqs, num_visis,
          num_baselines, num_components, num_times, beamtype,
          use_twoants, num_ants, off_cardinal_dipoles,
          primay_beam_J00, primay_beam_J01,
          primay_beam_J10, primay_beam_J11,
          visi_components,
          flux_I, flux_Q, flux_U, flux_V,
          sum_visi_XX_real, sum_visi_XX_imag,
          sum_visi_XY_real, sum_visi_XY_imag,
          sum_visi_YX_real, sum_visi_YX_imag,
          sum_visi_YY_real, sum_visi_YY_imag);
  }else {
    test_update_sum_visis_cpu(num_freqs, num_visis,
          num_baselines, num_components, num_times, beamtype,
          use_twoants, num_ants, off_cardinal_dipoles,
          primay_beam_J00, primay_beam_J01,
          primay_beam_J10, primay_beam_J11,
          visi_components,
          flux_I, flux_Q, flux_U, flux_V,
          sum_visi_XX_real, sum_visi_XX_imag,
          sum_visi_XY_real, sum_visi_XY_imag,
          sum_visi_YX_real, sum_visi_YX_imag,
          sum_visi_YY_real, sum_visi_YY_imag);
  }

  //These expected values have no cross-pols in them
  user_precision_t *expected_order = malloc(num_visis*sizeof(user_precision_t));
  for (int visi = 0; visi < num_visis; visi++) {
    expected_order[visi] = (visi + 1)*num_components;
  }

  if (beamtype == ANALY_DIPOLE || beamtype == GAUSS_BEAM) {
    for (int output = 0; output < num_visis; output++) {

      // printf("%.1f %.1f %.1f %.1f %.1f %.1f %.1f %.1f\n",
      //     sum_visi_XX_real[output], sum_visi_XX_imag[output],
      //     sum_visi_XY_real[output], sum_visi_XY_imag[output],
      //     sum_visi_YX_real[output], sum_visi_YX_imag[output],
      //     sum_visi_YY_real[output], sum_visi_YY_imag[output]);

      TEST_ASSERT_EQUAL_FLOAT(expected_order[output], sum_visi_XX_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_XY_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_YX_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(expected_order[output], sum_visi_YY_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(expected_order[output], sum_visi_XX_imag[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_XY_imag[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_YX_imag[output]);
      TEST_ASSERT_EQUAL_FLOAT(expected_order[output], sum_visi_YY_imag[output]);

    }
  }
  //No cross-pols so multiply expected values by 2.0
  else if (beamtype == FEE_BEAM || beamtype == FEE_BEAM_INTERP || beamtype == MWA_ANALY) {
    for (int output = 0; output < num_visis; output++) {

      TEST_ASSERT_EQUAL_FLOAT(2*expected_order[output], sum_visi_XX_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(2*expected_order[output], sum_visi_XY_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(2*expected_order[output], sum_visi_YX_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(2*expected_order[output], sum_visi_YY_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(2*expected_order[output], sum_visi_XX_imag[output]);
      TEST_ASSERT_EQUAL_FLOAT(2*expected_order[output], sum_visi_XY_imag[output]);
      TEST_ASSERT_EQUAL_FLOAT(2*expected_order[output], sum_visi_YX_imag[output]);
      TEST_ASSERT_EQUAL_FLOAT(2*expected_order[output], sum_visi_YY_imag[output]);

    }
  }
  else {
    //If there is no beam, should just end up with the number of components
    //as we apply a gain of 1.0
    for (int output = 0; output < num_visis; output++) {

      TEST_ASSERT_EQUAL_FLOAT(expected_order[output], sum_visi_XX_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_XY_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_YX_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(expected_order[output], sum_visi_YY_real[output]);
      TEST_ASSERT_EQUAL_FLOAT(expected_order[output], sum_visi_XX_imag[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_XY_imag[output]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, sum_visi_YX_imag[output]);
      TEST_ASSERT_EQUAL_FLOAT(expected_order[output], sum_visi_YY_imag[output]);

    }
  }
  //
  // //Be free my beauties
  free( primay_beam_J00 );
  free( primay_beam_J01 );
  free( primay_beam_J10 );
  free( primay_beam_J11 );
  free( visi_components );
  free( flux_I );
  free( flux_Q );
  free( flux_U );
  free( flux_V );
  free( sum_visi_XX_real );
  free( sum_visi_XX_imag );
  free( sum_visi_XY_real );
  free( sum_visi_XY_imag );
  free( sum_visi_YX_real );
  free( sum_visi_YX_imag );
  free( sum_visi_YY_real );
  free( sum_visi_YY_imag );

}