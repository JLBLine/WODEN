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
#include "source_components_cpu.h"

void extrapolate_Stokes(source_t *d_chunked_source, double *d_extrap_freqs,
                        int num_extrap_freqs, e_component_type comptype,
                        int do_gpu){

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
    if (do_gpu == 1) {
      extrap_power_laws_stokesI_gpu(d_components, n_powers, d_extrap_freqs, num_extrap_freqs);
    }
  }
  //Next up, do the CURVED_POWER_LAW types
  if (n_curves > 0) {
    if (do_gpu == 1) {
      extrap_curved_power_laws_stokesI_gpu(d_components, n_curves, d_extrap_freqs, num_extrap_freqs);
    }
  }

  //Finally, do any list flux peeps
  if (n_lists > 0) {
    if (do_gpu == 1) {
      extrap_list_fluxes_stokesI_gpu(d_components, n_lists, d_extrap_freqs, num_extrap_freqs);
    }
  }

  if (d_components.n_stokesV_power > 0) {
    if (do_gpu == 1) {
      extrap_power_laws_stokesV_gpu(d_components, d_extrap_freqs, num_extrap_freqs);
    }
  }

  if (d_components.n_stokesV_curve > 0) {
    if (do_gpu == 1) {
      extrap_curved_power_laws_stokesV_gpu(d_components, d_extrap_freqs, num_extrap_freqs);
    }
  }

  if (d_components.n_stokesV_pol_frac > 0) {
    if (do_gpu == 1) {
      polarisation_fraction_stokesV_gpu(d_components, d_extrap_freqs, num_extrap_freqs);
    }
  }

  if (d_components.n_stokesV_list > 0) {
    if (do_gpu == 1) {
      extrap_list_fluxes_stokesV_gpu(d_components, d_extrap_freqs, num_extrap_freqs);
    }
  }

  if (d_components.n_linpol_power > 0) {
    if (do_gpu == 1) {
      extrap_power_laws_linpol_gpu(d_components, d_extrap_freqs, num_extrap_freqs);
    }
  }

  if (d_components.n_linpol_curve > 0) {
    if (do_gpu == 1) {
      extrap_curved_power_laws_linpol_gpu(d_components, d_extrap_freqs, num_extrap_freqs);
    }
  }

  if (d_components.n_linpol_pol_frac > 0) {
    if (do_gpu == 1) {
      polarisation_fraction_linpol_gpu(d_components, d_extrap_freqs, num_extrap_freqs);
    }
  }

  if (d_components.n_linpol_list > 0) {
    if (do_gpu == 1) {
      extrap_list_fluxes_linpol_gpu(d_components, d_extrap_freqs, num_extrap_freqs);
    }
  }

  if (d_components.n_linpol_p_list > 0) {
    if (do_gpu == 1) {
      extrap_p_list_fluxes_linpol_gpu(d_components, d_extrap_freqs, num_extrap_freqs);
    }
  }

  if (d_components.n_linpol_angles > 0) {
    if (do_gpu == 1) {
      apply_rotation_measure_gpu(d_components, d_extrap_freqs, num_extrap_freqs);
    }
  }
}

void source_component_common(woden_settings_t *woden_settings,
           beam_settings_t *beam_settings, double *d_freqs,
           source_t *chunked_source, source_t *d_chunked_source,
           beam_gains_t *d_component_beam_gains,
           e_component_type comptype,
           visibility_set_t *d_visibility_set){

  int do_gpu = 1;
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

  if (do_gpu == 1) {
    malloc_extrapolated_flux_arrays_gpu(d_components, num_components,
                                     woden_settings->num_freqs);
    extrapolate_Stokes(d_chunked_source, d_freqs,
                     woden_settings->num_freqs, comptype, do_gpu);

    int num_gains = d_components->num_primarybeam_values*num_beams;
    malloc_beam_gains_gpu(d_component_beam_gains, beam_settings->beamtype, num_gains);
    calc_lmn_for_components_gpu(d_components, num_components, woden_settings);
  }

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

    if (do_gpu == 1){
      wrapper_calculate_gaussian_beam_gpu(num_components, cos_theta, sin_theta,
                                          sin_2theta, fwhm_lm, woden_settings,
                                          beam_settings, components,
                                          d_component_beam_gains, d_freqs);
    }

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

    if (do_gpu == 1){
      wrapper_run_hyperbeam_gpu(num_components, components, beam_settings,
                          num_beams, parallactic,
                          reordered_azs, reordered_zas,
                          d_component_beam_gains,
                          d_freqs, woden_settings);
    }
    free(reordered_azs);
    free(reordered_zas);
  }

  else if (beam_settings->beamtype == ANALY_DIPOLE) {
    if (verbose == 1){
      printf("\tDoing analytic_dipole (EDA2 beam)\n");
    }
    if (do_gpu == 1){
      wrapper_calculate_analytic_dipole_beam_gpu(num_components, components,
                          d_component_beam_gains, d_freqs, woden_settings);
    }
  }

  else if (beam_settings->beamtype == MWA_ANALY) {
    //Always normalise to zenith
    int norm = 1;
    if (verbose == 1){
      printf("\tDoing analytic MWA beam\n");
    }

    if (do_gpu == 1) {
      wrapper_calculate_RTS_MWA_analytic_beam_gpu(num_components,
                          components, norm, d_component_beam_gains,
                          d_freqs, woden_settings);
    }
  }

  //Now we've calculated the beams, we can calculate the auto-correlations,
  //if so required

  if (woden_settings->do_autos){
    if (verbose == 1){
      printf("\tCalculating auto-correlations\n");
    }

    if (do_gpu == 1){
      calc_autos_gpu(d_components, beam_settings, d_component_beam_gains,
                     d_visibility_set, woden_settings, num_components,
                     use_twobeams);
    }
  }
} //END source_component_common