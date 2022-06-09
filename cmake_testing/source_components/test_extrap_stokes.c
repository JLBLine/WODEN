#include <unity.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "constants.h"
#include "woden_precision_defs.h"
#include "woden_struct_defs.h"
#include "test_extrap_stokes.h"

#include "common_testing_functions.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

//External CUDA code we're linking in
extern void test_extrap_stokes_all_models(source_t *chunked_source,
           int num_extrap_freqs, double *extrap_freqs,
           user_precision_t *extrap_flux_I, user_precision_t *extrap_flux_Q,
           user_precision_t *extrap_flux_U, user_precision_t *extrap_flux_V);


// void fill_components_info(components_t *comps)

#ifdef DOUBLE_PRECISION
  double TOL = 1e-12;
#else
  double TOL = 1e-4;
#endif

/*
Test that the linear SI flux extrapolation code works correctly
*/
void test_kern_extrap_stokes_GivesCorrectValues(void) {

  source_t *chunked_source = malloc(sizeof(source_t));

  chunked_source->n_point_lists = num_lists;
  chunked_source->n_point_powers = num_powers;
  chunked_source->n_point_curves = num_curves;

  chunked_source->n_points = num_lists + num_powers + num_curves;

  components_t *comps = &chunked_source->point_components;

  //Set up some test condition inputs
  // int num_extrap_freqs = 25;
  int num_components = chunked_source->n_points;

  comps->power_ref_freqs = ref_freqs;
  comps->power_ref_stokesI = ref_stokesI;
  comps->power_ref_stokesQ = ref_stokesQ;
  comps->power_ref_stokesU = ref_stokesU;
  comps->power_ref_stokesV = ref_stokesV;
  comps->power_SIs = ref_power_SIs;

  comps->power_comp_inds = malloc(chunked_source->n_point_powers*sizeof(int));
  for (int pow_ind = 0; pow_ind < chunked_source->n_point_powers; pow_ind++) {
    comps->power_comp_inds[pow_ind] = pow_ind;
  }

  comps->curve_ref_freqs = ref_freqs;
  comps->curve_ref_stokesI = ref_stokesI;
  comps->curve_ref_stokesQ = ref_stokesQ;
  comps->curve_ref_stokesU = ref_stokesU;
  comps->curve_ref_stokesV = ref_stokesV;
  comps->curve_SIs = ref_curve_SIs;
  comps->curve_qs = ref_qs;

  comps->curve_comp_inds = malloc(chunked_source->n_point_curves*sizeof(int));
  for (int cur_ind = 0; cur_ind < chunked_source->n_point_curves; cur_ind++) {
    comps->curve_comp_inds[cur_ind] = chunked_source->n_point_powers + cur_ind;
  }


  comps->list_freqs = list_freqs;
  comps->list_stokesI = list_stokesI;
  comps->list_stokesQ = list_stokesQ;
  comps->list_stokesU = list_stokesU;
  comps->list_stokesV = list_stokesV;
  comps->num_list_values = num_list_values;
  comps->list_start_indexes = list_start_indexes;

  comps->list_comp_inds = malloc(chunked_source->n_point_lists*sizeof(int));
  for (int list_ind = 0; list_ind < chunked_source->n_point_lists; list_ind++) {
    comps->list_comp_inds[list_ind] = num_powers + num_curves + list_ind;
  }

  //This bit is normally done when remapping chunked sources before
  //copying over to the GPU

  //Sum up everything in num_list_values gives us how much we need to malloc
    int total_num_flux_entires = 0;
    for (int list_comp = 0; list_comp < num_lists; list_comp++) {
      total_num_flux_entires += comps->num_list_values[list_comp];
    }

    comps->total_num_flux_entires = total_num_flux_entires;


  //The generic function that copies sky models from CPU to GPU needs
  //things in the ra, dec etc arrays to be malloced

  comps->ras = malloc(chunked_source->n_points*sizeof(double));
  comps->decs = malloc(chunked_source->n_points*sizeof(double));

  //
  //Space for outputs
  user_precision_t *extrap_flux_I = malloc(num_extrap_freqs*num_components*sizeof(user_precision_t));
  user_precision_t *extrap_flux_Q = malloc(num_extrap_freqs*num_components*sizeof(user_precision_t));
  user_precision_t *extrap_flux_U = malloc(num_extrap_freqs*num_components*sizeof(user_precision_t));
  user_precision_t *extrap_flux_V = malloc(num_extrap_freqs*num_components*sizeof(user_precision_t));
  //
  // //Run the CUDA code
  test_extrap_stokes_all_models(chunked_source,
             num_extrap_freqs, extrap_freqs,
             extrap_flux_I, extrap_flux_Q,
             extrap_flux_U, extrap_flux_V);
  //
  //Make some expected value arrays
  double *expec_flux_I = malloc(num_extrap_freqs*num_components*sizeof(double));
  double *expec_flux_Q = malloc(num_extrap_freqs*num_components*sizeof(double));
  double *expec_flux_U = malloc(num_extrap_freqs*num_components*sizeof(double));
  double *expec_flux_V = malloc(num_extrap_freqs*num_components*sizeof(double));

  CPU_extrapolate_fluxes_in_components(comps, num_powers, num_curves, num_lists,
                        extrap_freqs, num_extrap_freqs,
                        expec_flux_I, expec_flux_Q, expec_flux_U, expec_flux_V);



  for (int i = 0; i < num_extrap_freqs*(num_powers + num_curves + num_lists); i++) {
    //Check the two are within tolerace
    // printf("%d %.3f %.3f\n",i, expec_flux_I[i], extrap_flux_I[i] );
    TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_flux_I[i], extrap_flux_I[i]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_flux_Q[i], extrap_flux_Q[i]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_flux_U[i], extrap_flux_U[i]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_flux_V[i], extrap_flux_V[i]);

  }

  FILE *output_text;

  output_text = fopen("test_extrap_stokes.txt","w");

  for (int i = 0; i < num_extrap_freqs*(num_powers + num_curves + num_lists); i++) {

    fprintf(output_text,"%.12f %.12f %.12f %.12f\n", extrap_flux_I[i],
                          extrap_flux_Q[i], extrap_flux_U[i], extrap_flux_V[i]);

  }

  //Be free my beauties
  free(extrap_flux_I);
  free(extrap_flux_Q);
  free(extrap_flux_U);
  free(extrap_flux_V);
  free(expec_flux_I);
  free(expec_flux_Q);
  free(expec_flux_U);
  free(expec_flux_V);

}

//Run the test with unity
int main(void)
{
    UNITY_BEGIN();
    RUN_TEST(test_kern_extrap_stokes_GivesCorrectValues);
    return UNITY_END();
}
