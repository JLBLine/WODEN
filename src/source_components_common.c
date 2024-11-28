#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
// #include "gpucomplex.h"
// #include "fundamental_coords_gpu.h"
// #include "source_components_gpu.h"
// #include "primary_beam_gpu.h"

// #include "gpu_macros.h"

#include "source_components_common.h"



void extrapolate_Stokes(source_t *d_chunked_source,
                        double *d_extrap_freqs, int num_extrap_freqs,
                        e_component_type comptype){

  components_t d_components;
  // int n_comps = 0;
  int n_powers = 0;
  int n_curves = 0;
  int n_lists = 0;

  //Choose the right components to extrapolate for
  if (comptype == POINT) {
    d_components = d_chunked_source->point_components;
    // n_comps = d_chunked_source->n_points;
    n_powers = d_chunked_source->n_point_powers;
    n_curves = d_chunked_source->n_point_curves;
    n_lists = d_chunked_source->n_point_lists;
  }
  else if (comptype == GAUSSIAN) {
    d_components = d_chunked_source->gauss_components;
    // n_comps = d_chunked_source->n_gauss;
    n_powers = d_chunked_source->n_gauss_powers;
    n_curves = d_chunked_source->n_gauss_curves;
    n_lists = d_chunked_source->n_gauss_lists;
  // } else if (comptype == SHAPELET) {
  } else {
    d_components = d_chunked_source->shape_components;
    // n_comps = d_chunked_source->n_shapes;
    n_powers = d_chunked_source->n_shape_powers;
    n_curves = d_chunked_source->n_shape_curves;
    n_lists = d_chunked_source->n_shape_lists;
  }

  if (n_powers > 0) {
    extrap_power_laws_stokesI_gpu(d_components, n_powers, d_extrap_freqs, num_extrap_freqs);
  }
  //Next up, do the CURVED_POWER_LAW types
  if (n_curves > 0) {
    extrap_curved_power_laws_stokesI_gpu(d_components, n_curves, d_extrap_freqs, num_extrap_freqs);
  }

  //Finally, do any list flux peeps
  if (n_lists > 0) {
    extrap_list_fluxes_stokesI_gpu(d_components, n_lists, d_extrap_freqs, num_extrap_freqs);
  }

  if (d_components.n_stokesV_power > 0) {
    extrap_power_laws_stokesV_gpu(d_components, d_extrap_freqs, num_extrap_freqs);
  }

  if (d_components.n_stokesV_curve > 0) {
    extrap_curved_power_laws_stokesV_gpu(d_components, d_extrap_freqs, num_extrap_freqs);
  }

  if (d_components.n_stokesV_pol_frac > 0) {
    polarisation_fraction_stokesV_gpu(d_components, d_extrap_freqs, num_extrap_freqs);
  }

  if (d_components.n_stokesV_list > 0) {
    extrap_list_fluxes_stokesV_gpu(d_components, d_extrap_freqs, num_extrap_freqs);
  }

  if (d_components.n_linpol_power > 0) {
    extrap_power_laws_linpol_gpu(d_components, d_extrap_freqs, num_extrap_freqs);
  }

  if (d_components.n_linpol_curve > 0) {
    extrap_curved_power_laws_linpol_gpu(d_components, d_extrap_freqs, num_extrap_freqs);
  }

  if (d_components.n_linpol_pol_frac > 0) {
    polarisation_fraction_linpol_gpu(d_components, d_extrap_freqs, num_extrap_freqs);
  }

  if (d_components.n_linpol_list > 0) {
    extrap_list_fluxes_linpol_gpu(d_components, d_extrap_freqs, num_extrap_freqs);
  }

  if (d_components.n_linpol_p_list > 0) {
    extrap_p_list_fluxes_linpol_gpu(d_components, d_extrap_freqs, num_extrap_freqs);
  }

  if (d_components.n_linpol_angles > 0) {
    apply_rotation_measure_gpu(d_components, d_extrap_freqs, num_extrap_freqs);
  }
}

extern "C" void source_component_common(woden_settings_t *woden_settings,
           beam_settings_t *beam_settings, double *d_freqs,
           source_t *chunked_source, source_t *d_chunked_source,
           d_beam_gains_t *d_component_beam_gains,
           e_component_type comptype,
           visibility_set_t *d_visibility_set){

  int verbose = 0;

  //Here we see if we a single primary beam for all (num_beams = 1) or
  //a primary beam per antenna (num_beams = num_ants)
  //This can be expanded in the future to have a primary beam per tile
  //for different options
  int num_beams;
  int use_twobeams = 0;
  if (woden_settings->use_dipamps == 1) {
    num_beams = woden_settings->num_ants;
    use_twobeams = 1;
  } else {
      num_beams = 1;
  }

  int num_components = 0;
  components_t *components = NULL;
  components_t *d_components = NULL;

  if (comptype == POINT) {
    num_components = d_chunked_source->n_points;
    components = &chunked_source->point_components;
    d_components = &d_chunked_source->point_components;
  } else if (comptype == GAUSSIAN) {
    num_components = d_chunked_source->n_gauss;
    components = &chunked_source->gauss_components;
    d_components = &d_chunked_source->gauss_components;
  } else if (comptype == SHAPELET) {
    num_components = d_chunked_source->n_shapes;
    components = &chunked_source->shape_components;
    d_components = &d_chunked_source->shape_components;
  }

  //Will need this later
  malloc_extrapolated_flux_arrays_gpu(d_components, num_components,
                                     woden_settings->num_freqs);

  extrapolate_Stokes(d_chunked_source, d_freqs,
                     woden_settings->num_freqs, comptype);

  int num_gains = d_components->num_primarybeam_values*num_beams;

  //TODO keep translating down from here
  
  //If we're using an everybeam model, all memory and values have already
  //been copied to GPU, so no need to allocate here
  //
  if (beam_settings->beamtype == FEE_BEAM || beam_settings->beamtype == MWA_ANALY || beam_settings->beamtype == FEE_BEAM_INTERP
      || beam_settings->beamtype == GAUSS_BEAM || beam_settings->beamtype == ANALY_DIPOLE || beam_settings->beamtype == NO_BEAM) {

    //Only some models would have had leakage terms malloced
    if (beam_settings->beamtype == FEE_BEAM || beam_settings->beamtype == MWA_ANALY || beam_settings->beamtype == FEE_BEAM_INTERP) {
      gpuMalloc( (void**)&d_component_beam_gains->d_Dxs,
                      num_gains*sizeof(gpuUserComplex) );
      gpuMalloc( (void**)&d_component_beam_gains->d_Dys,
                      num_gains*sizeof(gpuUserComplex) );
    }
    gpuMalloc( (void**)&d_component_beam_gains->d_gxs,
                      num_gains*sizeof(gpuUserComplex) );
    gpuMalloc( (void**)&d_component_beam_gains->d_gys,
                      num_gains*sizeof(gpuUserComplex) );

  }
  
  gpuMalloc( (void**)&d_components->ls, num_components*sizeof(double));
  gpuMalloc( (void**)&d_components->ms, num_components*sizeof(double));
  gpuMalloc( (void**)&d_components->ns, num_components*sizeof(double));


  dim3 grid, threads;

  threads.x = 128;
  threads.y = 1;
  threads.z = 1;
  grid.x = (int)ceil( (float)num_components / (float)threads.x );
  grid.y = 1;
  grid.z = 1;

  gpuErrorCheckKernel("kern_calc_lmn",
                        kern_calc_lmn, grid, threads,
                        woden_settings->ra0,
                        woden_settings->sdec0, woden_settings->cdec0,
                        d_components->ras, d_components->decs,
                        d_components->ls, d_components->ms, d_components->ns, num_components);

  //If using a gaussian primary beam, calculate beam values for all freqs,
  //lsts and point component locations
  if (beam_settings->beamtype == GAUSS_BEAM) {

    //TODO currently hardcoded to have beam position angle = 0.
    //Should this change with az/za?
    user_precision_t cos_theta = 1.0;
    user_precision_t sin_theta = 0.0;
    user_precision_t sin_2theta = 0.0;
    user_precision_t fwhm_lm = sin(beam_settings->beam_FWHM_rad);

    if (verbose == 1){
      printf("\tDoing Gaussian Beam\n");
    }

    calculate_gaussian_beam_gpu(num_components,
         woden_settings->num_time_steps, woden_settings->num_freqs,
         beam_settings->gauss_ha, beam_settings->gauss_sdec,
         beam_settings->gauss_cdec,
         fwhm_lm, cos_theta, sin_theta, sin_2theta,
         beam_settings->beam_ref_freq, d_freqs,
         components->beam_has,
         components->beam_decs,
         d_component_beam_gains->d_gxs, d_component_beam_gains->d_gys);

  }// end if beam == GAUSS

  else if (beam_settings->beamtype == FEE_BEAM || beam_settings->beamtype == FEE_BEAM_INTERP) {

    if (verbose == 1){
      if (beam_settings->beamtype == FEE_BEAM_INTERP) {
        printf("\tDoing the hyperbeam (interpolated)\n");
      } else {
        printf("\tDoing the hyperbeam\n");
      }
    }

    //Have to reorder the az/za from comp ind, time ind to time ind, comp ind
    //before feeding into mwa_hyperbeam
    int num_azza = woden_settings->num_time_steps*num_components;
    double *reordered_azs = (double *)malloc(num_azza*sizeof(double));
    double *reordered_zas = (double *)malloc(num_azza*sizeof(double));

    int stripe_new, stripe_old;

    for (int time_ind = 0; time_ind < woden_settings->num_time_steps; time_ind++) {
      for (int comp_ind = 0; comp_ind < num_components; comp_ind++) {
        stripe_new = time_ind*num_components + comp_ind;
        stripe_old = comp_ind*woden_settings->num_time_steps + time_ind;
        reordered_azs[stripe_new] = (double)components->azs[stripe_old];
        reordered_zas[stripe_new] = (double)components->zas[stripe_old];
      }
    }

    //Always be doing parallatic angle rotation
    uint8_t parallactic = 1;

    run_hyperbeam_gpu(num_components,
           woden_settings->num_time_steps, woden_settings->num_freqs,
           num_beams, parallactic,
           beam_settings->gpu_fee_beam,
           reordered_azs, reordered_zas,
           woden_settings->latitudes,
           d_component_beam_gains->d_gxs, d_component_beam_gains->d_Dxs,
           d_component_beam_gains->d_Dys, d_component_beam_gains->d_gys);

    free(reordered_azs);
    free(reordered_zas);
  }

  else if (beam_settings->beamtype == ANALY_DIPOLE) {
    if (verbose == 1){
      printf("\tDoing analytic_dipole (EDA2 beam)\n");
    }

    calculate_analytic_dipole_beam_gpu(num_components,
         woden_settings->num_time_steps, woden_settings->num_freqs,
         components->azs, components->zas, d_freqs,
         d_component_beam_gains->d_gxs, d_component_beam_gains->d_gys);
  }

  else if (beam_settings->beamtype == MWA_ANALY) {

    //Always normalise to zenith
    int norm = 1;
    if (verbose == 1){
      printf("\tDoing analytic MWA beam\n");
    }

    calculate_RTS_MWA_analytic_beam_gpu(num_components,
         woden_settings->num_time_steps, woden_settings->num_freqs,
         components->azs, components->zas,
         woden_settings->FEE_ideal_delays, woden_settings->latitude,
         norm, components->beam_has, components->beam_decs,
         d_freqs, d_component_beam_gains->d_gxs, d_component_beam_gains->d_Dxs,
         d_component_beam_gains->d_Dys, d_component_beam_gains->d_gys);
  }

  //Now we've calculated the beams, we can calculate the auto-correlations,
  //if so required

  if (woden_settings->do_autos){

    int num_freqs = woden_settings->num_freqs;
    int num_times = woden_settings->num_time_steps;
    int num_ants = woden_settings->num_ants;
    int num_baselines = woden_settings->num_baselines;

    threads.x = 64;
    threads.y = 2;
    threads.z = 1;
    grid.x = (int)ceil( (float)(num_freqs*num_times) / (float)threads.x );
    grid.y = (int)ceil( (float)(num_ants) / (float)threads.y );
    grid.z = 1;

    int *d_ant_to_auto_map = NULL;

    if (use_twobeams == 1) {
      int *ant_to_auto_map = NULL;
      ant_to_auto_map = (int *)malloc(num_ants*sizeof(int));
      for (int ant = 0; ant < num_ants; ant++){
          ant_to_auto_map[ant] = ant;
      }
      ( gpuMalloc( (void**)&d_ant_to_auto_map,
                                    num_ants*sizeof(int) ));
      ( gpuMemcpy(d_ant_to_auto_map, ant_to_auto_map,
                                      num_ants*sizeof(int), gpuMemcpyHostToDevice ));
      free(ant_to_auto_map);
    }

    gpuErrorCheckKernel("kern_calc_autos",
                  kern_calc_autos, grid, threads,
                  *d_components, *d_component_beam_gains,
                  beam_settings->beamtype,
                  num_components, num_baselines,
                  num_freqs, num_times, num_ants,
                  d_visibility_set->sum_visi_XX_real,
                  d_visibility_set->sum_visi_XX_imag,
                  d_visibility_set->sum_visi_XY_real,
                  d_visibility_set->sum_visi_XY_imag,
                  d_visibility_set->sum_visi_YX_real,
                  d_visibility_set->sum_visi_YX_imag,
                  d_visibility_set->sum_visi_YY_real,
                  d_visibility_set->sum_visi_YY_imag,
                  use_twobeams, d_ant_to_auto_map,
                  d_ant_to_auto_map,
                  woden_settings->off_cardinal_dipoles);

    if (use_twobeams == 1) {
    (  gpuFree( d_ant_to_auto_map ) );
    }
  }

} //END source_component_common

// extern "C" void fill_ant_to_baseline_mapping(int num_ants, int *d_ant1_to_baseline_map,
//                                                int *d_ant2_to_baseline_map){

//   int num_baselines = ((num_ants - 1)*num_ants) / 2;

//   int *ant1_to_baseline_map = NULL;
//   int *ant2_to_baseline_map = NULL;

//   ant1_to_baseline_map = (int *)malloc(num_baselines*sizeof(int));
//   ant2_to_baseline_map = (int *)malloc(num_baselines*sizeof(int));

//   //These functions only do cross correlations, so create all combos of antennas
//   //that make up all the crosses
//   int cross_index = 0;
//   for (int ant1 = 0; ant1 < num_ants-1; ant1++)
//   {
//     for (int ant2 = ant1 + 1; ant2 < num_ants; ant2++)
//     {
//       ant1_to_baseline_map[cross_index] = ant1;
//       ant2_to_baseline_map[cross_index] = ant2;

//       cross_index += 1;
//     }
//   }

//   ( gpuMemcpy(d_ant1_to_baseline_map, ant1_to_baseline_map,
//                                   num_baselines*sizeof(int), gpuMemcpyHostToDevice ));
//   ( gpuMemcpy(d_ant2_to_baseline_map, ant2_to_baseline_map,
//                                   num_baselines*sizeof(int), gpuMemcpyHostToDevice ));

//   free(ant1_to_baseline_map);
//   free(ant2_to_baseline_map);

// }