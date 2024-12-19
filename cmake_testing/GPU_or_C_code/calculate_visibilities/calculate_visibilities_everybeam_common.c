/*
`calculate_visibilities::calculate_visibilities` is the gateway function
to all CUDA functionality in WODEN. We'll test here in one baseline, frequency,
and time configuration. We'll vary the sky model and the primary beam. By
sticking all COMPONENTs at phase centre, we can just sum the expected fluxes
in XX / YY real to check things are being lanuched.

More variations like different phase centres / array configs etc are tested
in different test suites, so really just test that the correct CUDA functions
are launched by calculate_visibilities::calculate_visibilities`
*/



#include "calculate_visibilities_everybeam_common.h"

void set_beam_gains(source_catalogue_t *cropped_sky_models,
                    int n_points, int n_gauss, int n_shapes, int num_sources,
                    int num_times, int num_freqs,
                    double _Complex gain1x, double _Complex leak1x,
                    double _Complex leak1y, double _Complex gain1y,
                    double _Complex gain2x, double _Complex leak2x,
                    double _Complex leak2y, double _Complex gain2y,
                    double *antx_mult, double *anty_mult){
  int num_gains_per_comp = NUM_ANTS*num_times*num_freqs;

  // components_t *components = malloc(sizeof(components_t));
  // components_t *components;

  if (n_points > 0) {
    for (int source_ind = 0; source_ind < num_sources; source_ind++) {
      // components = &cropped_sky_models->sources[source_ind].point_components;
      cropped_sky_models->sources[source_ind].point_components.gxs = malloc(n_points*num_gains_per_comp*sizeof(user_precision_complex_t));
      cropped_sky_models->sources[source_ind].point_components.Dxs = malloc(n_points*num_gains_per_comp*sizeof(user_precision_complex_t));
      cropped_sky_models->sources[source_ind].point_components.Dys = malloc(n_points*num_gains_per_comp*sizeof(user_precision_complex_t));
      cropped_sky_models->sources[source_ind].point_components.gys = malloc(n_points*num_gains_per_comp*sizeof(user_precision_complex_t));
      // components->gxs = malloc(n_points*num_gains_per_comp*sizeof(user_precision_complex_t));
      // components->Dxs = malloc(n_points*num_gains_per_comp*sizeof(user_precision_complex_t));
      // components->Dys = malloc(n_points*num_gains_per_comp*sizeof(user_precision_complex_t));
      // components->gys = malloc(n_points*num_gains_per_comp*sizeof(user_precision_complex_t));
    }
  }

  if (n_gauss > 0) {
    for (int source_ind = 0; source_ind < num_sources; source_ind++) {
      cropped_sky_models->sources[source_ind].gauss_components.gxs = malloc(n_gauss*num_gains_per_comp*sizeof(user_precision_complex_t));
      cropped_sky_models->sources[source_ind].gauss_components.Dxs = malloc(n_gauss*num_gains_per_comp*sizeof(user_precision_complex_t));
      cropped_sky_models->sources[source_ind].gauss_components.Dys = malloc(n_gauss*num_gains_per_comp*sizeof(user_precision_complex_t));
      cropped_sky_models->sources[source_ind].gauss_components.gys = malloc(n_gauss*num_gains_per_comp*sizeof(user_precision_complex_t));
    }
  }

  if (n_shapes > 0) {
    for (int source_ind = 0; source_ind < num_sources; source_ind++) {
      cropped_sky_models->sources[source_ind].shape_components.gxs = malloc(n_shapes*num_gains_per_comp*sizeof(user_precision_complex_t));
      cropped_sky_models->sources[source_ind].shape_components.Dxs = malloc(n_shapes*num_gains_per_comp*sizeof(user_precision_complex_t));
      cropped_sky_models->sources[source_ind].shape_components.Dys = malloc(n_shapes*num_gains_per_comp*sizeof(user_precision_complex_t));
      cropped_sky_models->sources[source_ind].shape_components.gys = malloc(n_shapes*num_gains_per_comp*sizeof(user_precision_complex_t));
    }
  }

  int gain_ind;

  // int num_times = woden_settings->num_time_steps;
  // int num_freqs = woden_settings->num_freqs;

  user_precision_complex_t gx, Dx, Dy, gy;

  for (int source_ind = 0; source_ind < num_sources; source_ind++) {
    for (int beam_ind = 0; beam_ind < NUM_ANTS; beam_ind++) {
      for (int time_ind = 0; time_ind < num_times; time_ind++) {
        if (time_ind == 0){
          gx = gain1x;
          Dx = leak1x;
          Dy = leak1y;
          gy = gain1y;
        } else {
          gx = gain2x;
          Dx = leak2x;
          Dy = leak2y;
          gy = gain2y;
        }
        
        for (int freq_ind = 0; freq_ind < num_freqs; freq_ind++) {
          for (int comp_ind = 0; comp_ind < n_points; comp_ind++) {
            gain_ind = beam_ind*num_times*num_freqs*n_points + time_ind*num_freqs*n_points + freq_ind*n_points + comp_ind;
            cropped_sky_models->sources[source_ind].point_components.gxs[gain_ind] = antx_mult[beam_ind]*gx;
            cropped_sky_models->sources[source_ind].point_components.Dxs[gain_ind] = antx_mult[beam_ind]*Dx;
            cropped_sky_models->sources[source_ind].point_components.Dys[gain_ind] = anty_mult[beam_ind]*Dy;
            cropped_sky_models->sources[source_ind].point_components.gys[gain_ind] = anty_mult[beam_ind]*gy;
            
          }
          for (int comp_ind = 0; comp_ind < n_gauss; comp_ind++) {
            gain_ind = beam_ind*num_times*num_freqs*n_gauss + time_ind*num_freqs*n_gauss + freq_ind*n_gauss + comp_ind;
            cropped_sky_models->sources[source_ind].gauss_components.gxs[gain_ind] = antx_mult[beam_ind]*gx;
            cropped_sky_models->sources[source_ind].gauss_components.Dxs[gain_ind] = antx_mult[beam_ind]*Dx;
            cropped_sky_models->sources[source_ind].gauss_components.Dys[gain_ind] = anty_mult[beam_ind]*Dy;
            cropped_sky_models->sources[source_ind].gauss_components.gys[gain_ind] = anty_mult[beam_ind]*gy;
            
          }
          for (int comp_ind = 0; comp_ind < n_shapes; comp_ind++) {
            gain_ind = beam_ind*num_times*num_freqs*n_shapes + time_ind*num_freqs*n_shapes + freq_ind*n_shapes + comp_ind;
            cropped_sky_models->sources[source_ind].shape_components.gxs[gain_ind] = antx_mult[beam_ind]*gx;
            cropped_sky_models->sources[source_ind].shape_components.Dxs[gain_ind] = antx_mult[beam_ind]*Dx;
            cropped_sky_models->sources[source_ind].shape_components.Dys[gain_ind] = anty_mult[beam_ind]*Dy;
            cropped_sky_models->sources[source_ind].shape_components.gys[gain_ind] = anty_mult[beam_ind]*gy;
            
          }
        }
      }
    }
  }
}



void test_calculate_visibilities_EveryBeam(int n_points, int n_gauss, int n_shapes,
                                           int num_sources, int do_gpu, int beamtype) {

  source_catalogue_t *cropped_sky_models = make_cropped_sky_models(RA0, -0.46606083776035967,
                                                    n_points, n_gauss, n_shapes,
                                                    num_sources);

  woden_settings_t *woden_settings = make_woden_settings(RA0, -0.46606083776035967);
  woden_settings->beamtype = beamtype;
  woden_settings->do_gpu = do_gpu;
  // woden_settings->use_dipamps = 0;
  // woden_settings
  woden_settings->num_ants = NUM_ANTS;

  //When doing everybeam, the gains are calculated via Python and stored in
  //the sky model. So make some beam gains here in a way that we can predict
  //the output visibilities

  double antx_mult[3] = {0.1, 0.5, 0.9};
  double anty_mult[3] = {0.2, 0.6, 1.0};

  double _Complex gain1x, leak1x, leak1y, gain1y;
  double _Complex gain2x, leak2x, leak2y, gain2y;

  gain1x = 1.0 + I*0.1;
  leak1x = -0.01 + I*-0.001;
  leak1y = 0.02 + I*0.002;
  gain1y = 0.9 + I*0.2;
  gain2x = 0.5 + I*0.1;
  leak2x = -0.001 + I*-0.001;
  leak2y = 0.002 + I*0.002;
  gain2y = 0.3 + I*0.4;

  set_beam_gains(cropped_sky_models, n_points, n_gauss, n_shapes, num_sources,
                    woden_settings->num_time_steps, woden_settings->num_freqs,
                    gain1x, leak1x, leak1y, gain1y,
                    gain2x, leak2x, leak2y, gain2y,
                    antx_mult, anty_mult);
  


  beam_settings_t *beam_settings = malloc(sizeof(beam_settings_t));
  beam_settings->beamtype = beamtype;

  visibility_set_t *visibility_set;

  visibility_set = test_calculate_visibilities(cropped_sky_models,
                                          beam_settings, woden_settings, RA0, -0.46606083776035967,
                                          beam_settings->beamtype);

  // printf("WODEN settings %d %d %d\n",woden_settings->num_time_steps,
  //                                    woden_settings->num_freqs,
  //                                    woden_settings->num_baselines );
  // for (int visi = 0; visi < NUM_VISI; visi++) {
  //   printf("%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n",
  //           visibility_set->sum_visi_XX_real[visi],
  //           visibility_set->sum_visi_XX_imag[visi],
  //           visibility_set->sum_visi_XY_real[visi],
  //           visibility_set->sum_visi_XY_imag[visi],
  //           visibility_set->sum_visi_YX_real[visi],
  //           visibility_set->sum_visi_YX_imag[visi],
  //           visibility_set->sum_visi_YY_real[visi],
  //           visibility_set->sum_visi_YY_imag[visi]);
  // }

  int num_comps = (n_points + n_gauss + n_shapes)*num_sources;

  #ifdef DOUBLE_PRECISION
    double TOL = 2e-6;
  #else
    double TOL = 3e-2;
  #endif

  test_comp_phase_centre_allgains_multiants(visibility_set, num_comps,
                                  gain1x, leak1x, leak1y, gain1y,
                                  gain2x, leak2x, leak2y, gain2y,
                                  antx_mult, anty_mult, NUM_ANTS,
                                  woden_settings, TOL);

  free_visi_set_inputs(visibility_set);
  free_visi_set_outputs(visibility_set);
  // // // free_beam_gains(cropped_sky_models, n_points, n_gauss, n_shapes, num_sources);

  woden_settings->do_autos = 1;
  woden_settings->num_autos = NUM_CROSS;
  woden_settings->num_visis = woden_settings->num_cross + woden_settings->num_autos;

  // test_calculate_visibilities frees up the sky model, so we need to remake it
  cropped_sky_models = make_cropped_sky_models(RA0, -0.46606083776035967,
                                                    n_points, n_gauss, n_shapes,
                                                    num_sources);

  set_beam_gains(cropped_sky_models, n_points, n_gauss, n_shapes, num_sources,
                    woden_settings->num_time_steps, woden_settings->num_freqs,
                    gain1x, leak1x, leak1y, gain1y,
                    gain2x, leak2x, leak2y, gain2y,
                    antx_mult, anty_mult);

  // printf("We have this many visis %d %d %d\n",woden_settings->num_visis,woden_settings->num_autos,woden_settings->num_cross );
  visibility_set = test_calculate_visibilities(cropped_sky_models,
                                          beam_settings, woden_settings, RA0, -0.46606083776035967,
                                          beam_settings->beamtype);

  test_comp_phase_centre_allgains_multiants(visibility_set,num_comps,
                                  gain1x, leak1x, leak1y, gain1y,
                                  gain2x, leak2x, leak2y, gain2y,
                                  antx_mult, anty_mult, NUM_ANTS,
                                  woden_settings, TOL);

  // free_beam_gains(cropped_sky_models, n_points, n_gauss, n_shapes, num_sources);
  free(beam_settings);
  free_visi_set_inputs(visibility_set);
  free_visi_set_outputs(visibility_set);
  free(woden_settings);
}

//Test with a single SOURCE, single COMPONENT
void test_calculate_visibilities_EveryBeam_OneSource_SinglePoint(int do_gpu,int beamtype) {
  int n_points = 1;
  int n_gauss = 0;
  int n_shapes = 0;
  int num_sources = 1;
  test_calculate_visibilities_EveryBeam(n_points, n_gauss, n_shapes,
                                               num_sources, do_gpu, beamtype);
}

void test_calculate_visibilities_EveryBeam_OneSource_SingleGauss(int do_gpu, int beamtype) {
  int n_points = 0;
  int n_gauss = 1;
  int n_shapes = 0;
  int num_sources = 1;
  test_calculate_visibilities_EveryBeam(n_points, n_gauss, n_shapes,
                                               num_sources, do_gpu, beamtype);

}

void test_calculate_visibilities_EveryBeam_OneSource_SingleShape(int do_gpu, int beamtype) {
  int n_points = 0;
  int n_gauss = 0;
  int n_shapes = 1;
  int num_sources = 1;
  test_calculate_visibilities_EveryBeam(n_points, n_gauss, n_shapes,
                                               num_sources, do_gpu, beamtype);
}

void test_calculate_visibilities_EveryBeam_OneSource_SingleAll(int do_gpu, int beamtype) {
  int n_points = 1;
  int n_gauss = 1;
  int n_shapes = 1;
  int num_sources = 1;
  test_calculate_visibilities_EveryBeam(n_points, n_gauss, n_shapes,
                                               num_sources, do_gpu, beamtype);
}


//Test with a three SOURCEs, single COMPONENT
void test_calculate_visibilities_EveryBeam_ThreeSource_SinglePoint(int do_gpu, int beamtype) {
  int n_points = 1;
  int n_gauss = 0;
  int n_shapes = 0;
  int num_sources = 3;
  test_calculate_visibilities_EveryBeam(n_points, n_gauss, n_shapes,
                                               num_sources, do_gpu, beamtype);

}

void test_calculate_visibilities_EveryBeam_ThreeSource_SingleGauss(int do_gpu, int beamtype) {
  int n_points = 0;
  int n_gauss = 1;
  int n_shapes = 0;
  int num_sources = 3;
  test_calculate_visibilities_EveryBeam(n_points, n_gauss, n_shapes,
                                               num_sources, do_gpu, beamtype);
}

void test_calculate_visibilities_EveryBeam_ThreeSource_SingleShape(int do_gpu, int beamtype) {
  int n_points = 0;
  int n_gauss = 0;
  int n_shapes = 1;
  int num_sources = 3;
  test_calculate_visibilities_EveryBeam(n_points, n_gauss, n_shapes,
                                               num_sources, do_gpu, beamtype);
}

void test_calculate_visibilities_EveryBeam_ThreeSource_SingleAll(int do_gpu, int beamtype) {
  int n_points = 1;
  int n_gauss = 1;
  int n_shapes = 1;
  int num_sources = 3;
  test_calculate_visibilities_EveryBeam(n_points, n_gauss, n_shapes,
                                               num_sources, do_gpu, beamtype);
}


//Test with three SOURCEs, three COPMONENTs
void test_calculate_visibilities_EveryBeam_ThreeSource_FivePoint(int do_gpu, int beamtype) {
  int n_points = 5;
  int n_gauss = 0;
  int n_shapes = 0;
  int num_sources = 3;
  test_calculate_visibilities_EveryBeam(n_points, n_gauss, n_shapes,
                                               num_sources, do_gpu, beamtype);

}

void test_calculate_visibilities_EveryBeam_ThreeSource_FiveGauss(int do_gpu, int beamtype) {
  int n_points = 0;
  int n_gauss = 5;
  int n_shapes = 0;
  int num_sources = 3;
  test_calculate_visibilities_EveryBeam(n_points, n_gauss, n_shapes,
                                               num_sources, do_gpu, beamtype);
}

void test_calculate_visibilities_EveryBeam_ThreeSource_FiveShape(int do_gpu, int beamtype) {
  int n_points = 0;
  int n_gauss = 0;
  int n_shapes = 5;
  int num_sources = 3;
  test_calculate_visibilities_EveryBeam(n_points, n_gauss, n_shapes,
                                               num_sources, do_gpu, beamtype);
}

void test_calculate_visibilities_EveryBeam_ThreeSource_FiveAll(int do_gpu, int beamtype) {
  int n_points = 5;
  int n_gauss = 5;
  int n_shapes = 5;
  int num_sources = 3;
  test_calculate_visibilities_EveryBeam(n_points, n_gauss, n_shapes,
                                               num_sources, do_gpu, beamtype);
}