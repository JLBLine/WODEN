#include "calculate_visibilities_cpu.h"

calc_visi_inouts_t * create_calc_visi_inouts_cpu(array_layout_t *array_layout,
        visibility_set_t *visibility_set,  user_precision_t *sbf,
        woden_settings_t *woden_settings,
        int num_shapelets, int use_twobeams) {

  calc_visi_inouts_t *calc_visi_inouts = (calc_visi_inouts_t *)malloc(sizeof(calc_visi_inouts_t));

  // int num_cross = woden_settings->num_cross;
  int num_baselines = woden_settings->num_baselines;
  int num_time_steps = woden_settings->num_time_steps;
  int num_visis = woden_settings->num_visis;
  // int num_freqs = woden_settings->num_freqs;
  int num_ants = woden_settings->num_ants;

  calc_visi_inouts->allsteps_sha0s = visibility_set->allsteps_sha0s;
  calc_visi_inouts->allsteps_cha0s = visibility_set->allsteps_cha0s;
  calc_visi_inouts->allsteps_wavelengths = visibility_set->allsteps_wavelengths;

  calc_visi_inouts->u_metres = malloc(num_visis*sizeof(user_precision_t));
  calc_visi_inouts->v_metres = malloc(num_visis*sizeof(user_precision_t));
  calc_visi_inouts->w_metres = malloc(num_visis*sizeof(user_precision_t));
  calc_visi_inouts->us = malloc(num_visis*sizeof(user_precision_t));
  calc_visi_inouts->vs = malloc(num_visis*sizeof(user_precision_t));
  calc_visi_inouts->ws = malloc(num_visis*sizeof(user_precision_t));

  calc_visi_inouts->freqs = visibility_set->channel_frequencies;


  //if we have shapelets in our sky model, copy the shapelet basis functions
  //into GPU memory
  if (num_shapelets > 0) {
    calc_visi_inouts->sbf = sbf;
    //Technically we're assigning more memory than we need here, as we're
    //assigning memory for all shapelet sources in the whole source_catalogue
    //(we could instead assign memory for each chunked source)
    //However, this is a small amount of memory, and it saves multiple
    //mallocs and frees so is probably more efficient
    calc_visi_inouts->u_shapes = malloc(num_shapelets*num_baselines*num_time_steps*sizeof(user_precision_t));
    calc_visi_inouts->v_shapes = malloc(num_shapelets*num_baselines*num_time_steps*sizeof(user_precision_t));

  }

  if (use_twobeams == 1) {

    calc_visi_inouts->ant1_to_baseline_map = malloc(num_baselines*sizeof(int) );
    calc_visi_inouts->ant2_to_baseline_map = malloc(num_baselines*sizeof(int) );

    fill_ant_to_baseline_mapping_cpu(num_ants, calc_visi_inouts->ant1_to_baseline_map,
                                     calc_visi_inouts->ant2_to_baseline_map);
  }
  return calc_visi_inouts;
}

void set_visi_set_to_zero_cpu(visibility_set_t *visibility_set, int num_visis) {
  for (int visi = 0; visi < num_visis; visi++) {
      visibility_set->sum_visi_XX_real[visi] = 0;
      visibility_set->sum_visi_XX_imag[visi] = 0;
      visibility_set->sum_visi_XY_real[visi] = 0;
      visibility_set->sum_visi_XY_imag[visi] = 0;
      visibility_set->sum_visi_YX_real[visi] = 0;
      visibility_set->sum_visi_YX_imag[visi] = 0;
      visibility_set->sum_visi_YY_real[visi] = 0;
      visibility_set->sum_visi_YY_imag[visi] = 0;
  }
}

void free_calc_visi_inouts_cpu(calc_visi_inouts_t *calc_visi_inouts,
                               int num_shapelets, int use_twobeams){

  free(calc_visi_inouts->u_metres);
  free(calc_visi_inouts->v_metres);
  free(calc_visi_inouts->w_metres);
  free(calc_visi_inouts->us);
  free(calc_visi_inouts->vs);
  free(calc_visi_inouts->ws);

  //if we have shapelets in our sky model, copy the shapelet basis functions
  //into GPU memory
  if (num_shapelets > 0) {
    free(calc_visi_inouts->u_shapes);
    free(calc_visi_inouts->v_shapes);

  }

  if (use_twobeams == 1) {
    free(calc_visi_inouts->ant1_to_baseline_map);
    free(calc_visi_inouts->ant2_to_baseline_map);
  }
}

void free_components_cpu(components_t *components) {
  free(components->ls);
  free(components->ms);
  free(components->ns);

  free( components->extrap_stokesI );

  if (components->do_QUV) {
    free(components->extrap_stokesQ );
    free(components->extrap_stokesU );
    free(components->extrap_stokesV );
  }
}

void free_beam_gains_cpu(beam_gains_t *beam_gains, e_beamtype beamtype){

  free(beam_gains->gxs );
  free(beam_gains->gys );

  if (beamtype == FEE_BEAM || beamtype == FEE_BEAM_INTERP || beamtype == MWA_ANALY || beamtype == EB_OSKAR || beamtype == EB_LOFAR || beamtype == EB_MWA){
    free(beam_gains->Dxs );
    free(beam_gains->Dys );
  }
}