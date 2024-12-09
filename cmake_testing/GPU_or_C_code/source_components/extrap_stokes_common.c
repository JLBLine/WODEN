#include "extrap_stokes_common.h"
#include "test_extrap_stokes.h"

//External GPU code we're linking in
extern source_t * copy_chunked_source_to_GPU(source_t *chunked_source);

extern void malloc_extrapolated_flux_arrays_gpu(components_t *d_components, int num_comps,
                                     int num_freqs);

extern void free_d_components(source_t *d_chunked_source,
                                  e_component_type comptype);

extern void free_extrapolated_flux_arrays(components_t *d_components);

extern double * malloc_freqs_gpu(int num_extrap_freqs, double *extrap_freqs);

extern void free_freqs_gpu(double *d_extrap_freqs);

extern void copy_extrapolated_flux_arrays_to_host(source_t *d_chunked_source,
                                                  int num_extrap_freqs,
                                                  user_precision_t *extrap_flux_I,
                                                  user_precision_t *extrap_flux_Q,
                                                  user_precision_t *extrap_flux_U,
                                                  user_precision_t *extrap_flux_V);

#ifdef DOUBLE_PRECISION
  double TOL = 1e-11;
#else
  double TOL = 2e-3;
#endif

/*
Test that the linear SI flux extrapolation code works correctly
Many input arrays and values are stored in test_extrap_stokes.h

*/
void test_extrap_stokes_GivesCorrectValues(int do_gpu) {

  source_t *chunked_source = malloc(sizeof(source_t));

  chunked_source->n_point_lists = num_lists;
  chunked_source->n_point_powers = num_powers;
  chunked_source->n_point_curves = num_curves;

  chunked_source->n_points = num_lists + num_powers + num_curves;

  components_t *comps = &chunked_source->point_components;

  //Set up some test condition inputs
  int num_components = chunked_source->n_points;

  //Stokes I power laws---------------------------------------------------------
  // comps->power_ref_freqs = ref_freqs;
  comps->power_ref_stokesI = ref_stokesI;
  comps->power_SIs = stokesI_power_SIs;

  comps->power_comp_inds = malloc(chunked_source->n_point_powers*sizeof(int));
  for (int pow_ind = 0; pow_ind < chunked_source->n_point_powers; pow_ind++) {
    comps->power_comp_inds[pow_ind] = pow_ind;
  }


  //Stokes I curved power laws--------------------------------------------------
  // comps->curve_ref_freqs = ref_freqs;
  comps->curve_ref_stokesI = ref_stokesI;
  comps->curve_SIs = stokesI_curve_SIs;
  comps->curve_qs = stokesI_qs;

  comps->curve_comp_inds = malloc(chunked_source->n_point_curves*sizeof(int));
  for (int cur_ind = 0; cur_ind < chunked_source->n_point_curves; cur_ind++) {
    comps->curve_comp_inds[cur_ind] = chunked_source->n_point_powers + cur_ind;
  }

  //Stokes I list fluxes--------------------------------------------------------
  comps->list_freqs = list_freqs;
  comps->list_stokesI = list_stokesI;
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

  //Stokes V shizzle------------------------------------------------------------
  //Indexes of the component pols can match any of the Stokes I flux models,
  //so offset them to be sure we're getting all this indexing hell correct

  int v_ind = 0;

  comps->n_stokesV_list = num_pol_list;
  comps->stokesV_list_ref_freqs = stokesV_list_ref_freqs;
  comps->stokesV_list_ref_flux = stokesV_list_ref_flux;
  comps->stokesV_num_list_values = stokesV_num_list_values;
  comps->stokesV_list_start_indexes = stokesV_list_start_indexes;
  comps->stokesV_list_comp_inds = malloc(num_pol_list*sizeof(int));
  comps->n_stokesV_list_flux_entries = 0;
  for (int list_ind = 0; list_ind < num_pol_list; list_ind++) {
    comps->stokesV_list_comp_inds[list_ind] = v_ind;
    comps->n_stokesV_list_flux_entries += stokesV_num_list_values[list_ind];
    v_ind += 1;
  }

  comps->n_stokesV_pol_frac = num_pol_frac;
  comps->stokesV_pol_fracs = stokesV_pol_fracs;
  comps->stokesV_pol_frac_comp_inds = malloc(num_pol_frac*sizeof(int));
  for (int frac_ind = 0; frac_ind < num_pol_frac; frac_ind++) {
    comps->stokesV_pol_frac_comp_inds[frac_ind] = v_ind;
    v_ind += 1;
  }

  comps->n_stokesV_power = num_pol_power;
  comps->stokesV_power_ref_flux = ref_stokesV;
  comps->stokesV_power_SIs = stokesV_power_SIs;
  comps->stokesV_power_comp_inds = malloc(num_pol_power*sizeof(int));
  for (int pow_ind = 0; pow_ind < num_pol_power; pow_ind++) {
    comps->stokesV_power_comp_inds[pow_ind] = v_ind;
    v_ind += 1;
  }
  

  comps->n_stokesV_curve = num_pol_curv;
  comps->stokesV_curve_ref_flux = ref_stokesV;
  comps->stokesV_curve_SIs = stokesV_curve_SIs;
  comps->stokesV_curve_qs = stokesV_qs;
  comps->stokesV_curve_comp_inds = malloc(num_pol_curv*sizeof(int));
  for (int cur_ind = 0; cur_ind < num_pol_curv; cur_ind++) {
    comps->stokesV_curve_comp_inds[cur_ind] = v_ind;
    v_ind += 1;
  }

  int l_ind = 0;

  comps->n_linpol_pol_frac = num_pol_frac;
  comps->linpol_pol_fracs = linpol_pol_fracs;
  comps->linpol_pol_frac_comp_inds = malloc(num_pol_frac*sizeof(int));
  for (int frac_ind = 0; frac_ind < num_pol_frac; frac_ind++) {
    comps->linpol_pol_frac_comp_inds[frac_ind] = l_ind;
    l_ind += 1;
  }

  comps->n_linpol_power = num_pol_power;
  comps->linpol_power_ref_flux = ref_linpol;
  comps->linpol_power_SIs = linpol_power_SIs;
  comps->linpol_power_comp_inds = malloc(num_pol_power*sizeof(int));
  for (int pow_ind = 0; pow_ind < num_pol_power; pow_ind++) {
    comps->linpol_power_comp_inds[pow_ind] = l_ind;
    l_ind += 1;
  }
  

  comps->n_linpol_curve = num_pol_curv;
  comps->linpol_curve_ref_flux = ref_linpol;
  comps->linpol_curve_SIs = linpol_curve_SIs;
  comps->linpol_curve_qs = linpol_qs;
  comps->linpol_curve_comp_inds = malloc(num_pol_curv*sizeof(int));
  for (int cur_ind = 0; cur_ind < num_pol_curv; cur_ind++) {
    comps->linpol_curve_comp_inds[cur_ind] = l_ind;
    l_ind += 1;
  }

  comps->n_linpol_p_list = num_pol_list;
  comps->linpol_p_list_ref_freqs = linpol_p_list_ref_freqs;
  comps->linpol_p_list_ref_flux = linpol_p_list_ref_flux;
  comps->linpol_p_num_list_values = linpol_p_num_list_values;
  comps->linpol_p_list_start_indexes = linpol_p_list_start_indexes;
  comps->linpol_p_list_comp_inds = malloc(num_pol_list*sizeof(int));
  comps->n_linpol_p_list_flux_entries = 0;
  for (int list_ind = 0; list_ind < num_pol_list; list_ind++) {
    comps->linpol_p_list_comp_inds[list_ind] = l_ind;
    comps->n_linpol_p_list_flux_entries += linpol_p_num_list_values[list_ind];
    l_ind += 1;
  }

  comps->n_linpol_angles = num_pol_frac + num_pol_power + num_pol_curv + num_pol_list;
  comps->intr_pol_angle = intr_pol_angle;
  comps->rm_values = rms;
  comps->linpol_angle_inds = malloc(comps->n_linpol_angles*sizeof(int));
  for (int ind = 0; ind < comps->n_linpol_angles; ind++) {
    comps->linpol_angle_inds[ind] = ind;
  }

  //When you have both Q and U as lists, they don't need RM values, hence we've
  //aleady set all the RMs above

  comps->n_linpol_list = num_pol_list;
  comps->stokesQ_list_ref_freqs = stokesQ_list_ref_freqs;
  comps->stokesQ_list_ref_flux = stokesQ_list_ref_flux;
  comps->stokesQ_num_list_values = stokesQ_num_list_values;
  comps->stokesQ_list_start_indexes = stokesQ_list_start_indexes;
  comps->stokesU_list_ref_freqs = stokesU_list_ref_freqs;
  comps->stokesU_list_ref_flux = stokesU_list_ref_flux;
  comps->stokesU_num_list_values = stokesU_num_list_values;
  comps->stokesU_list_start_indexes = stokesU_list_start_indexes;

  comps->stokesQ_list_comp_inds = malloc(num_pol_list*sizeof(int));
  comps->n_stokesQ_list_flux_entries = 0;
  comps->stokesU_list_comp_inds = malloc(num_pol_list*sizeof(int));
  comps->n_stokesU_list_flux_entries = 0;
  for (int list_ind = 0; list_ind < num_pol_list; list_ind++) {
    comps->stokesQ_list_comp_inds[list_ind] = l_ind;
    comps->stokesU_list_comp_inds[list_ind] = l_ind;
    comps->n_stokesQ_list_flux_entries += stokesQ_num_list_values[list_ind];
    comps->n_stokesU_list_flux_entries += stokesU_num_list_values[list_ind];
    l_ind += 1;
  }



  comps->do_QUV = 1;

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
  if (do_gpu == 1){
    source_t *d_chunked_source = copy_chunked_source_to_GPU(chunked_source);

    malloc_extrapolated_flux_arrays_gpu(&d_chunked_source->point_components,
                                      d_chunked_source->n_points,
                                      num_extrap_freqs);

    double *d_extrap_freqs = malloc_freqs_gpu(num_extrap_freqs, extrap_freqs);

    extrapolate_Stokes(d_chunked_source, d_extrap_freqs, num_extrap_freqs, POINT,
                       do_gpu);

    copy_extrapolated_flux_arrays_to_host(d_chunked_source, num_extrap_freqs,
                                          extrap_flux_I, extrap_flux_Q,
                                          extrap_flux_U, extrap_flux_V);

    free_d_components(d_chunked_source, POINT);
    free_extrapolated_flux_arrays(&d_chunked_source->point_components);
    free_freqs_gpu(d_extrap_freqs);
  }
  // //
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
    // printf("I %d %.3f %.3f\n", i/num_extrap_freqs, expec_flux_I[i], extrap_flux_I[i] );
    // printf("Q %d %.3f %.3f\n", i/num_extrap_freqs, expec_flux_Q[i], extrap_flux_Q[i] );
    // printf("U %d %.3f %.3f\n", i/num_extrap_freqs, expec_flux_U[i], extrap_flux_U[i] );
    // printf("V %d %.3f %.3f\n", i/num_extrap_freqs, expec_flux_V[i], extrap_flux_V[i] );
    // printf("%d %.3f %.3f\n",i, expec_flux_V[i], extrap_flux_V[i] );
    TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_flux_I[i], extrap_flux_I[i]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_flux_Q[i], extrap_flux_Q[i]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_flux_U[i], extrap_flux_U[i]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_flux_V[i], extrap_flux_V[i]);
  }

  // FILE *output_text;

  // output_text = fopen("test_extrap_stokes.txt","w");

  // for (int i = 0; i < num_extrap_freqs*(num_powers + num_curves + num_lists); i++) {

  //   fprintf(output_text,"%.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f\n", extrap_flux_I[i],
  //                     extrap_flux_Q[i], extrap_flux_U[i], extrap_flux_V[i],
  //         expec_flux_I[i], expec_flux_Q[i], expec_flux_U[i], expec_flux_V[i]);

  // }

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