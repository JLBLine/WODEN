#include "source_components_common.h"
#include "source_components_cpu.h"
#include "calculate_visibilities_cpu.h"
#include "call_everybeam_c.h"
#include <string.h>
#include "logger.h"

//NOTE this works as both a CPU and GPU function. Anything with `mem` at the
//front means you either need host or device memory for that variable. E.g.
//if `do_gpu == 1` then `mem_chunked_source` must be in device memory.
//chunked_source and extrap_freqs MUST be in GPU memory
//The thinking of calling the GPU/CPU code from the same function, is that
//for future developers, you can tell where you're going to have to develop
//both CPU/GPU code.
void extrapolate_Stokes(source_t *mem_chunked_source, double *mem_extrap_freqs,
                        int num_extrap_freqs, e_component_type comptype,
                        int do_gpu){

  components_t mem_components;
  // int n_comps = 0;
  int n_powers = 0;
  int n_curves = 0;
  int n_lists = 0;

  //Choose the right mem_components to extrapolate for
  if (comptype == POINT) {
    mem_components = mem_chunked_source->point_components;
    n_powers = mem_chunked_source->n_point_powers;
    n_curves = mem_chunked_source->n_point_curves;
    n_lists = mem_chunked_source->n_point_lists;
  }
  else if (comptype == GAUSSIAN) {
    mem_components = mem_chunked_source->gauss_components;
    n_powers = mem_chunked_source->n_gauss_powers;
    n_curves = mem_chunked_source->n_gauss_curves;
    n_lists = mem_chunked_source->n_gauss_lists;
  } else {
    mem_components = mem_chunked_source->shape_components;
    n_powers = mem_chunked_source->n_shape_powers;
    n_curves = mem_chunked_source->n_shape_curves;
    n_lists = mem_chunked_source->n_shape_lists;
  }

  if (n_powers > 0) {
    if (do_gpu == 1) {
      extrap_power_laws_stokesI_gpu(mem_components, n_powers, mem_extrap_freqs, num_extrap_freqs);
    } else {
      extrap_power_laws_stokesI_cpu(mem_components, n_powers, mem_extrap_freqs, num_extrap_freqs);
    }
  }
  //Next up, do the CURVED_POWER_LAW types
  if (n_curves > 0) {
    if (do_gpu == 1) {
      extrap_curved_power_laws_stokesI_gpu(mem_components, n_curves, mem_extrap_freqs, num_extrap_freqs);
    } else {
      extrap_curved_power_laws_stokesI_cpu(mem_components, n_curves, mem_extrap_freqs, num_extrap_freqs);
    }
  }

  //Finally, do any list flux peeps
  if (n_lists > 0) {
    if (do_gpu == 1) {
      extrap_list_fluxes_stokesI_gpu(mem_components, n_lists, mem_extrap_freqs, num_extrap_freqs);
    } else {
      extrap_list_fluxes_stokesI_cpu(mem_components, n_lists, mem_extrap_freqs, num_extrap_freqs);
    }
  }

  if (mem_components.n_stokesV_power > 0) {
    if (do_gpu == 1) {
      extrap_power_laws_stokesV_gpu(mem_components, mem_extrap_freqs, num_extrap_freqs);
    } else {
      extrap_power_laws_stokesV_cpu(mem_components, mem_extrap_freqs, num_extrap_freqs);
    }
  }

  if (mem_components.n_stokesV_curve > 0) {
    if (do_gpu == 1) {
      extrap_curved_power_laws_stokesV_gpu(mem_components, mem_extrap_freqs, num_extrap_freqs);
    } else {
      extrap_curved_power_laws_stokesV_cpu(mem_components, mem_extrap_freqs, num_extrap_freqs);
    }
  }

  if (mem_components.n_stokesV_pol_frac > 0) {
    if (do_gpu == 1) {
      polarisation_fraction_stokesV_gpu(mem_components, mem_extrap_freqs, num_extrap_freqs);
    } else {
      polarisation_fraction_stokesV_cpu(mem_components, mem_extrap_freqs, num_extrap_freqs);
    }
  }

  if (mem_components.n_stokesV_list > 0) {
    if (do_gpu == 1) {
      extrap_list_fluxes_stokesV_gpu(mem_components, mem_extrap_freqs, num_extrap_freqs);
    } else {
      extrap_list_fluxes_stokesV_cpu(mem_components, mem_extrap_freqs, num_extrap_freqs);
    }
  }

  if (mem_components.n_linpol_power > 0) {
    if (do_gpu == 1) {
      extrap_power_laws_linpol_gpu(mem_components, mem_extrap_freqs, num_extrap_freqs);
    } else {
      extrap_power_laws_linpol_cpu(mem_components, mem_extrap_freqs, num_extrap_freqs);
    }
  }

  if (mem_components.n_linpol_curve > 0) {
    if (do_gpu == 1) {
      extrap_curved_power_laws_linpol_gpu(mem_components, mem_extrap_freqs, num_extrap_freqs);
    } else {
      extrap_curved_power_laws_linpol_cpu(mem_components, mem_extrap_freqs, num_extrap_freqs);
    }
  }

  if (mem_components.n_linpol_pol_frac > 0) {
    if (do_gpu == 1) {
      polarisation_fraction_linpol_gpu(mem_components, mem_extrap_freqs, num_extrap_freqs);
    } else {
      polarisation_fraction_linpol_cpu(mem_components, mem_extrap_freqs, num_extrap_freqs);
    }
  }

  if (mem_components.n_linpol_list > 0) {
    if (do_gpu == 1) {
      extrap_list_fluxes_linpol_gpu(mem_components, mem_extrap_freqs, num_extrap_freqs);
    } else {
      extrap_list_fluxes_linpol_cpu(mem_components, mem_extrap_freqs, num_extrap_freqs);
    }
  }

  if (mem_components.n_linpol_p_list > 0) {
    if (do_gpu == 1) {
      extrap_p_list_fluxes_linpol_gpu(mem_components, mem_extrap_freqs, num_extrap_freqs);
    } else {
      extrap_p_list_fluxes_linpol_cpu(mem_components, mem_extrap_freqs, num_extrap_freqs);
    }
  }
  if (mem_components.n_linpol_angles > 0) {
    if (do_gpu == 1) {
      apply_rotation_measure_gpu(mem_components, mem_extrap_freqs, num_extrap_freqs);
    } else {
      apply_rotation_measure_cpu(mem_components, mem_extrap_freqs, num_extrap_freqs);
    }
  }
}

//TODO get cpu_freqs in here as well so can always pass to everybeam
void source_component_common(woden_settings_t *woden_settings,
           beam_settings_t *beam_settings,
           double *cpu_freqs, double *mem_freqs,
           source_t *chunked_source, source_t *mem_chunked_source,
           beam_gains_t *mem_component_beam_gains,
           e_component_type comptype,
           visibility_set_t *mem_visibility_set){

  int do_gpu = woden_settings->do_gpu;

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

  //Default behaviour with everybeam is to use a different beam for each station
  if (beam_settings->beamtype == EB_LOFAR || beam_settings->beamtype == EB_OSKAR) {
    use_twobeams = 1;
    num_beams = woden_settings->num_ants;

    //However if single_everybeam_station is set, we're only using one beam
    //for all stations
    if (woden_settings->single_everybeam_station == 1) {
      use_twobeams = 0;
      num_beams = 1;
    }
  }

  int num_components = 0;
  components_t *components = NULL;
  components_t *mem_components = NULL;

  if (comptype == POINT) {
    num_components = mem_chunked_source->n_points;
    components = &chunked_source->point_components;
    mem_components = &mem_chunked_source->point_components;
  } else if (comptype == GAUSSIAN) {
    num_components = mem_chunked_source->n_gauss;
    components = &chunked_source->gauss_components;
    mem_components = &mem_chunked_source->gauss_components;
  } else if (comptype == SHAPELET) {
    num_components = mem_chunked_source->n_shapes;
    components = &chunked_source->shape_components;
    mem_components = &mem_chunked_source->shape_components;
  }

  // int num_gains = mem_components->num_primarybeam_values*num_beams;
  int num_gains = num_components*woden_settings->num_time_steps*woden_settings->num_freqs*num_beams;
  if (do_gpu == 1) {
    malloc_extrapolated_flux_arrays_gpu(mem_components, num_components,
                                        woden_settings->num_freqs);
    malloc_beam_gains_gpu(mem_component_beam_gains, beam_settings->beamtype, num_gains);
    calc_lmn_for_components_gpu(mem_components, num_components, woden_settings);
  } else {
    malloc_extrapolated_flux_arrays_cpu(mem_components, num_components,
                                        woden_settings->num_freqs);
    malloc_beam_gains_cpu(mem_component_beam_gains, beam_settings->beamtype, num_gains);
    calc_lmn_for_components_cpu(mem_components, num_components, woden_settings);
    
  }

  extrapolate_Stokes(mem_chunked_source, mem_freqs,
                     woden_settings->num_freqs, comptype, do_gpu);

  //If using a gaussian primary beam, calculate beam values for all freqs,
  //lsts and point component locations
  if (beam_settings->beamtype == GAUSS_BEAM) {

    //TODO currently hardcoded to have beam position angle = 0.
    //Should this change with az/za?
    user_precision_t cos_theta = 1.0;
    user_precision_t sin_theta = 0.0;
    user_precision_t sin_2theta = 0.0;
    user_precision_t fwhm_lm = sin(beam_settings->beam_FWHM_rad);

    if (woden_settings->verbose == 1){
      log_message("\tDoing Gaussian Beam");
    }

    if (do_gpu == 1){
      wrapper_calculate_gaussian_beam_gpu(num_components, cos_theta, sin_theta,
                                          sin_2theta, fwhm_lm, woden_settings,
                                          beam_settings, components,
                                          mem_component_beam_gains, mem_freqs);
    } else {
      calculate_gaussian_beam_cpu(num_components, woden_settings->num_time_steps,
                                  woden_settings->num_freqs,
                                  beam_settings->gauss_ha,
                                  beam_settings->gauss_sdec,
                                  beam_settings->gauss_cdec,
                                  fwhm_lm, cos_theta, sin_theta, sin_2theta,
                                  beam_settings->beam_ref_freq, mem_freqs,
                                  components->beam_has,
                                  components->beam_decs,
                                  mem_component_beam_gains->gxs,
                                  mem_component_beam_gains->gys);
    }

  }// end if beam == GAUSS

  else if (beam_settings->beamtype == FEE_BEAM || beam_settings->beamtype == FEE_BEAM_INTERP) {

    if (woden_settings->verbose == 1){
      if (beam_settings->beamtype == FEE_BEAM_INTERP) {
        log_message("\tDoing the hyperbeam (interpolated)");
      } else {
        log_message("\tDoing the hyperbeam");
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
      wrapper_run_hyperbeam_gpu(num_components, beam_settings,
                                num_beams, parallactic,
                                reordered_azs, reordered_zas,
                                mem_component_beam_gains,
                                mem_freqs, woden_settings);
    } else {
      //Righto, so the CPU and GPU versions of hyperbeam have different APIs.
      //When you initiate the GPU, you pass in the delays and amplitudes
      //once, and then call the beam for as many directions as you want.
      //With the CPU version, you pass in the delays and amplitudes each time
      //you call for all directions. So we need to work out how many amplitudes
      //we have here; 16 for a single beam, because for a single beam we're
      //doing the perfect amp == 1 case, so we can pass the same amplitudes
      //for both polarisations. For two beams, a.k.a every tile has it's own
      //set of amplitude, we need 32 amplitudes, 16 for each polarisation.
      int num_amps;
      if (use_twobeams == 0) {
        num_amps = 16;
      } else {
        num_amps = 32;
      }

      run_hyperbeam_cpu(num_components, woden_settings->num_time_steps,
                        woden_settings->num_freqs, num_beams,
                        parallactic, mem_freqs, beam_settings->fee_beam,
                        beam_settings->hyper_delays,
                        num_amps, woden_settings->mwa_dipole_amps,
                        reordered_azs, reordered_zas,
                        woden_settings->latitudes,
                        mem_component_beam_gains->gxs, mem_component_beam_gains->Dxs,
                        mem_component_beam_gains->Dys, mem_component_beam_gains->gys);
                      

    }
    free(reordered_azs);
    free(reordered_zas);
  }

  else if (beam_settings->beamtype == ANALY_DIPOLE) {
    if (woden_settings->verbose == 1){
      log_message("\tDoing analytic_dipole (EDA2 beam)");
    }
    if (do_gpu == 1){
      wrapper_calculate_analytic_dipole_beam_gpu(num_components, components,
                          mem_component_beam_gains, mem_freqs, woden_settings);
    } else {
      calculate_analytic_dipole_beam_cpu(num_components,
                woden_settings->num_time_steps, woden_settings->num_freqs,
                mem_components->azs, mem_components->zas, mem_freqs,
                mem_component_beam_gains->gxs, mem_component_beam_gains->gys);
    }
  }

  else if (beam_settings->beamtype == MWA_ANALY) {
    //Always normalise to zenith
    int norm = 1;
    if (woden_settings->verbose == 1){
      log_message("\tDoing analytic MWA beam");
    }

    if (do_gpu == 1) {
      wrapper_calculate_RTS_MWA_analytic_beam_gpu(num_components,
                          components, norm, mem_component_beam_gains,
                          mem_freqs, woden_settings);
    } else {
      calculate_RTS_MWA_analytic_beam_cpu(num_components,
            woden_settings->num_time_steps, woden_settings->num_freqs,
            mem_components->azs, mem_components->zas, 
            woden_settings->FEE_ideal_delays,
            woden_settings->latitude, norm,
            components->beam_has, components->beam_decs, mem_freqs,
            mem_component_beam_gains->gxs, mem_component_beam_gains->Dxs,
            mem_component_beam_gains->Dys, mem_component_beam_gains->gys);
    }
  }

  //Only include this code if we've passed the -DHAVE_EVERYBEAM compilation flag
  #if defined(HAVE_EVERYBEAM)
  else if (beam_settings->beamtype == EB_MWA || beam_settings->beamtype == EB_LOFAR || beam_settings->beamtype == EB_OSKAR) {
    if (woden_settings->verbose == 1){
      if (beam_settings->beamtype == EB_MWA) {
        log_message("\tDoing EveryBeam MWA beam");
      } else if (beam_settings->beamtype == EB_LOFAR) {
        log_message("\tDoing EveryBeam LOFAR beam");
      } else {
        log_message("\tDoing EveryBeam OSKAR beam");
      }
    }

    beam_gains_t *eb_beam_gains = malloc(sizeof(beam_gains_t));
    malloc_beam_gains_cpu(eb_beam_gains, beam_settings->beamtype, num_gains);

    int eb_status = 0;
    
    int station_idxs[num_beams];
    for (int beam = 0; beam < num_beams; beam++){
      station_idxs[beam] = beam;
    }

    double mjd_sec_times[woden_settings->num_time_steps];

    for (int timei = 0; timei < woden_settings->num_time_steps; timei++) {
      mjd_sec_times[timei] = woden_settings->mjds[timei]*86400.0;
    }

    bool apply_beam_norms = woden_settings->normalise_primary_beam;
    bool rotate = true;
    bool element_only = false;
    bool iau_order = true;

    int num_times = woden_settings->num_time_steps;
    int num_freqs = woden_settings->num_freqs;

    double _Complex *jones = malloc(MAX_POLS*num_components*num_times*num_freqs*num_beams*sizeof(double _Complex));

    if (beam_settings->beamtype == EB_MWA) {
      

      double para_angles[num_components*num_times];
      double azs[num_components*num_times];
      double zas[num_components*num_times];

      for (int comp = 0; comp < num_components*num_times; comp++) {
        para_angles[comp] = (double)components->para_angles[comp];
        azs[comp] = (double)components->azs[comp];
        zas[comp] = (double)components->zas[comp];
      }

      //MWA beam is already normalised to zenith
      apply_beam_norms = false;

      run_mwa_beam(beam_settings->everybeam_telescope,
                    num_beams, station_idxs, num_components,
                    azs, zas, para_angles,
                    num_times, mjd_sec_times,
                    num_freqs, cpu_freqs,
                    apply_beam_norms, rotate, element_only, iau_order,
                    jones);
    }

    if (beam_settings->beamtype == EB_LOFAR || beam_settings->beamtype == EB_OSKAR) {

      run_phased_array_beam(beam_settings->everybeam_telescope,
                            num_beams, station_idxs, num_components,
                            woden_settings->eb_beam_ra0,
                            woden_settings->eb_beam_dec0,
                            components->ras, components->decs,
                            num_times, mjd_sec_times,
                            num_freqs, cpu_freqs,
                            apply_beam_norms,
                            rotate, element_only, iau_order,
                            jones);

    }

    if (eb_status != 0) {
      log_message("WARNING - Something went wrong runnng the EveryBeam telescope");
    }

    for (int station = 0; station < num_beams; station ++) {
      for (int time = 0; time < num_times; time ++) {
        for (int freq = 0; freq < num_freqs; freq ++) {
          for (int comp = 0; comp < num_components; comp ++) {

          int beam_ind = (num_freqs*num_times*num_components*station + num_freqs*num_components*time + num_components*freq + comp);

          int jones_index = 4*beam_ind;

          eb_beam_gains->gxs[beam_ind] = (user_precision_complex_t)jones[jones_index + 0];
          eb_beam_gains->Dxs[beam_ind] = (user_precision_complex_t)jones[jones_index + 1];
          eb_beam_gains->Dys[beam_ind] = (user_precision_complex_t)jones[jones_index + 2];
          eb_beam_gains->gys[beam_ind] = (user_precision_complex_t)jones[jones_index + 3];
          

          }
        }
      }
    }

    if (do_gpu == 1){
      copy_CPU_beam_gains_to_GPU_beam_gains(eb_beam_gains, mem_component_beam_gains, num_gains);  
    } else {
      //It seems wasteful to copy across the beam gains, but `mem_component_beam_gains`
      //memory allocation is handled elsewhere, so this is cleaner in terms of
      //memory management.
      memcpy(mem_component_beam_gains->gxs, eb_beam_gains->gxs, num_gains*sizeof(user_precision_complex_t));
      memcpy(mem_component_beam_gains->Dxs, eb_beam_gains->Dxs, num_gains*sizeof(user_precision_complex_t));
      memcpy(mem_component_beam_gains->Dys, eb_beam_gains->Dys, num_gains*sizeof(user_precision_complex_t));
      memcpy(mem_component_beam_gains->gys, eb_beam_gains->gys, num_gains*sizeof(user_precision_complex_t));

    }
    
    free_beam_gains_cpu(eb_beam_gains, beam_settings->beamtype);
    free(jones);
  }
  #endif //end if defined(HAVE_EVERYBEAM)

  //If a pyuvbeam model, already calculated beam gains on the CPU
  //So just copy them across
  //NOTE add this back in when doing pyuvbeam
  // else if (some pyuvbeam thing) {
  //   if (woden_settings->verbose == 1){
  //     log_message("\tDoing a pyuvbeam");
  //   }
  //   int num_gains = num_components*woden_settings->num_freqs*woden_settings->num_time_steps*num_beams;
  //   if (do_gpu == 1){
  //     copy_CPU_beam_gains_to_GPU(components, mem_component_beam_gains, num_gains);  
  //   } else {
  //     //It seems wasteful to copy across the beam gains, but `components` was created
  //     //on the python side, where the memory freeing is handled (hopefully) by
  //     //the garbage collector. By copying here, we can free the beam gains on
  //     //the C side regardless of whether we generate the beam gains in C or Python.
  //     memcpy(mem_component_beam_gains->gxs, components->gxs, num_gains*sizeof(user_precision_complex_t));
  //     memcpy(mem_component_beam_gains->Dxs, components->Dxs, num_gains*sizeof(user_precision_complex_t));
  //     memcpy(mem_component_beam_gains->Dys, components->Dys, num_gains*sizeof(user_precision_complex_t));
  //     memcpy(mem_component_beam_gains->gys, components->gys, num_gains*sizeof(user_precision_complex_t));

  //   }
  // }

  //Now we've calculated the beams, we can calculate the auto-correlations,
  //if so required

  if (woden_settings->do_autos){
    if (woden_settings->verbose == 1){
      log_message("\tCalculating auto-correlations");
    }

    if (do_gpu == 1){
      calc_autos_gpu(mem_components, beam_settings, mem_component_beam_gains,
                     mem_visibility_set, woden_settings, num_components,
                     use_twobeams);
    } else {
      int *ant_to_auto_map = malloc(woden_settings->num_ants*sizeof(int));
      for (int ant = 0; ant < woden_settings->num_ants; ant++){
        ant_to_auto_map[ant] = ant;
      }

      calc_autos_cpu(*mem_components, *mem_component_beam_gains,
                  beam_settings->beamtype, num_components,  woden_settings->num_baselines,
                  woden_settings->num_freqs, woden_settings->num_time_steps,
                  woden_settings->num_ants,
                  mem_visibility_set->sum_visi_XX_real, mem_visibility_set->sum_visi_XX_imag,
                  mem_visibility_set->sum_visi_XY_real, mem_visibility_set->sum_visi_XY_imag,
                  mem_visibility_set->sum_visi_YX_real, mem_visibility_set->sum_visi_YX_imag,
                  mem_visibility_set->sum_visi_YY_real, mem_visibility_set->sum_visi_YY_imag,
                  use_twobeams, ant_to_auto_map, ant_to_auto_map,
                  woden_settings->off_cardinal_dipoles);
      free(ant_to_auto_map);
    }
  }
} //END source_component_common

void fill_ant_to_baseline_mapping_cpu(int num_ants, int *ant1_to_baseline_map,
                                      int *ant2_to_baseline_map){

  // int num_baselines = ((num_ants - 1)*num_ants) / 2;
  int cross_index = 0;
  for (int ant1 = 0; ant1 < num_ants; ant1++) {
    for (int ant2 = ant1+1; ant2 < num_ants; ant2++) {
      ant1_to_baseline_map[cross_index] = ant1;
      ant2_to_baseline_map[cross_index] = ant2;
      cross_index++;
    }
  }
}