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

#define MWA_RA0 0.75723171
#define MWA_DEC0 -0.46606083776035967
double MWA_mjds[] ={56428.45056713, 56428.57556713};
double MWA_lsts[] ={0.75723171, 1.54478021};
double MWA_azs[] ={0.00000000, 4.52781227};
double MWA_zas[] ={0.00000000, 0.69969735};
double MWA_para_angles[] ={0.00000000, 1.75537303};

#define LOFAR_RA0 2.15374123
#define LOFAR_DEC0 0.84155210
double LOFAR_mjds[] ={57050.84349547, 57050.96849547};
double LOFAR_lsts[] ={1.36655619, 2.15410469};

void set_azza_para(source_catalogue_t *cropped_sky_models,
                    int n_points, int n_gauss, int n_shapes, int num_sources,
                    int beamtype){

  // components_t *components = malloc(sizeof(components_t));
  // components_t *components;

  double azs[NUM_TIME_STEPS];
  double zas[NUM_TIME_STEPS];
  double para_angles[NUM_TIME_STEPS];
  

  if (beamtype == EB_MWA) {
    for (int time_ind = 0; time_ind < NUM_TIME_STEPS; time_ind++) {
      azs[time_ind] = MWA_azs[time_ind];
      zas[time_ind] = MWA_zas[time_ind];
      para_angles[time_ind] = MWA_para_angles[time_ind];
    }
  }

  int ind;

  if (n_points > 0) {
    for (int source_ind = 0; source_ind < num_sources; source_ind++) {
      for (int comp_ind = 0; comp_ind < n_points; comp_ind++) {
        for (int time_ind = 0; time_ind < NUM_TIME_STEPS; time_ind++) {
          
          ind = comp_ind*NUM_TIME_STEPS + time_ind;
          cropped_sky_models->sources[source_ind].point_components.azs[ind] = azs[time_ind];
          cropped_sky_models->sources[source_ind].point_components.zas[ind] = zas[time_ind];
          cropped_sky_models->sources[source_ind].point_components.para_angles[ind] = para_angles[time_ind];

        }
      }
    }
  }

  if (n_gauss > 0) {
    for (int source_ind = 0; source_ind < num_sources; source_ind++) {
      for (int comp_ind = 0; comp_ind < n_gauss; comp_ind++) {
        for (int time_ind = 0; time_ind < NUM_TIME_STEPS; time_ind++) {
          
          ind = comp_ind*NUM_TIME_STEPS + time_ind;
          cropped_sky_models->sources[source_ind].gauss_components.azs[ind] = azs[time_ind];
          cropped_sky_models->sources[source_ind].gauss_components.zas[ind] = zas[time_ind];
          cropped_sky_models->sources[source_ind].gauss_components.para_angles[ind] = para_angles[time_ind];
        }
      }
    }
  }

  if (n_shapes > 0) {
    for (int source_ind = 0; source_ind < num_sources; source_ind++) {
      for (int comp_ind = 0; comp_ind < n_shapes; comp_ind++) {
        for (int time_ind = 0; time_ind < NUM_TIME_STEPS; time_ind++) {
          
          ind = comp_ind*NUM_TIME_STEPS + time_ind;
          cropped_sky_models->sources[source_ind].shape_components.azs[ind] = azs[time_ind];
          cropped_sky_models->sources[source_ind].shape_components.zas[ind] = zas[time_ind];
          cropped_sky_models->sources[source_ind].shape_components.para_angles[ind] = para_angles[time_ind];
        }
      }
    }
  }

}

void set_mjds(woden_settings_t *woden_settings, int beamtype){
  
    woden_settings->mjds = malloc(NUM_TIME_STEPS*sizeof(double));

    if (beamtype == EB_MWA) {
      for (int time_ind = 0; time_ind < NUM_TIME_STEPS; time_ind++) {
        woden_settings->mjds[time_ind] = MWA_mjds[time_ind];
      }
    }

    else if (beamtype == EB_LOFAR) {
      for (int time_ind = 0; time_ind < NUM_TIME_STEPS; time_ind++) {
        woden_settings->mjds[time_ind] = LOFAR_mjds[time_ind];
      }
    }
}



void test_calculate_visibilities_EveryBeam(int n_points, int n_gauss, int n_shapes,
                                           int num_sources, int do_gpu, int beamtype,
                                           const char *beam_ms_path) {


  double ra0, dec0;

  if (beamtype == EB_MWA) {
    ra0 = MWA_RA0;
    dec0 = MWA_DEC0;
  } else if (beamtype == EB_LOFAR) {
    ra0 = LOFAR_RA0;
    dec0 = LOFAR_DEC0;
  }


  // source_catalogue_t *cropped_sky_models = make_cropped_sky_models(RA0, -0.46606083776035967,
  //                                                   n_points, n_gauss, n_shapes,
  //                                                   num_sources);

  // woden_settings_t *woden_settings = make_woden_settings(RA0, -0.46606083776035967);

  source_catalogue_t *cropped_sky_models = make_cropped_sky_models(RA0, dec0,
                                                    n_points, n_gauss, n_shapes,
                                                    num_sources);

  woden_settings_t *woden_settings = make_woden_settings(RA0, dec0);

  woden_settings->beamtype = beamtype;
  woden_settings->do_gpu = do_gpu;
  woden_settings->use_dipamps = 0;
  woden_settings->do_autos = 0;
  woden_settings->num_ants = NUM_ANTS;
  woden_settings->beam_ms_path = beam_ms_path;
  woden_settings->hdf5_beam_path = getenv("MWA_FEE_HDF5");
  woden_settings->off_cardinal_dipoles = 0;

  if (beamtype == EB_MWA) {
    woden_settings->single_everybeam_station = 1;
  } else {
    woden_settings->single_everybeam_station = 0;
  }

  set_mjds(woden_settings, beamtype);

  set_azza_para(cropped_sky_models, n_points, n_gauss, n_shapes, num_sources,
                beamtype);

  // printf("OUTIDE BEFORE CALL %.5f %.5f %.5f\n", cropped_sky_models->sources[0].point_components.azs[0],
  //                                     cropped_sky_models->sources[1].point_components.azs[0],
  //                                     cropped_sky_models->sources[2].point_components.azs[0] );

  beam_settings_t *beam_settings = malloc(sizeof(beam_settings_t));
  beam_settings->beamtype = beamtype;

  visibility_set_t *visibility_set;

  visibility_set = test_calculate_visibilities(cropped_sky_models,
                                               beam_settings, woden_settings, RA0,
                                               dec0,
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
    double TOL = 1e-5;
  #else
    double TOL = 3e-2;
  #endif

  double _Complex *gain1x = malloc(NUM_ANTS*NUM_TIME_STEPS*sizeof(double _Complex));
  double _Complex *leak1x = malloc(NUM_ANTS*NUM_TIME_STEPS*sizeof(double _Complex));
  double _Complex *leak1y = malloc(NUM_ANTS*NUM_TIME_STEPS*sizeof(double _Complex));
  double _Complex *gain1y = malloc(NUM_ANTS*NUM_TIME_STEPS*sizeof(double _Complex));
  double _Complex *gain2x = malloc(NUM_ANTS*NUM_TIME_STEPS*sizeof(double _Complex));
  double _Complex *leak2x = malloc(NUM_ANTS*NUM_TIME_STEPS*sizeof(double _Complex));
  double _Complex *leak2y = malloc(NUM_ANTS*NUM_TIME_STEPS*sizeof(double _Complex));
  double _Complex *gain2y = malloc(NUM_ANTS*NUM_TIME_STEPS*sizeof(double _Complex));

  // -0.8775036335 0.4795700312, 0.0001940064 0.0000344838, 0.0002021146 0.0000273663, -0.8777521849 0.4791148901
  // -0.0692513958 0.0402538478, -0.0234183073 0.0134399701, 0.0233655237 -0.0105750393, -0.0530412868 0.0235209484
  
  //MWA EveryBeam is locked to zenith, so beam values vary depending on where things are pointing
  if (beamtype == EB_MWA) {
    double _Complex gain1x_mwa[NUM_ANTS*NUM_TIME_STEPS] = {-0.8775036335 + 0.4795700312*I, -0.8775036335 + 0.4795700312*I, -0.8775036335 + 0.4795700312*I,
                                                          -0.0692513958 + 0.0402538478*I, -0.0692513958 + 0.0402538478*I, -0.0692513958 + 0.0402538478*I};

    double _Complex leak1x_mwa[NUM_ANTS*NUM_TIME_STEPS] = {0.0001940064 + 0.0000344838*I, 0.0001940064 + 0.0000344838*I, 0.0001940064 + 0.0000344838*I,
                                                          -0.0234183073 + 0.0134399701*I, -0.0234183073 + 0.0134399701*I, -0.0234183073 + 0.0134399701*I};

    double _Complex leak1y_mwa[NUM_ANTS*NUM_TIME_STEPS] = {0.0002021146 + 0.0000273663*I, 0.0002021146 + 0.0000273663*I, 0.0002021146 + 0.0000273663*I,
                                                           0.0233655237 - 0.0105750393*I, 0.0233655237 - 0.0105750393*I, 0.0233655237 - 0.0105750393*I};

    double _Complex gain1y_mwa[NUM_ANTS*NUM_TIME_STEPS] = {-0.8777521849 + 0.4791148901*I, -0.8777521849 + 0.4791148901*I, -0.8777521849 + 0.4791148901*I,
                                                          -0.0530412868 + 0.0235209484*I, -0.0530412868 + 0.0235209484*I, -0.0530412868 + 0.0235209484*I};


    double _Complex gain2x_mwa[NUM_ANTS*NUM_TIME_STEPS] = {-0.8775036335 + 0.4795700312*I, -0.8775036335 + 0.4795700312*I, -0.8775036335 + 0.4795700312*I,
                                                          -0.0692513958 + 0.0402538478*I, -0.0692513958 + 0.0402538478*I, -0.0692513958 + 0.0402538478*I};

    double _Complex leak2x_mwa[NUM_ANTS*NUM_TIME_STEPS] = {0.0001940064 + 0.0000344838*I, 0.0001940064 + 0.0000344838*I, 0.0001940064 + 0.0000344838*I,
                                                          -0.0234183073 + 0.0134399701*I, -0.0234183073 + 0.0134399701*I, -0.0234183073 + 0.0134399701*I};

    double _Complex leak2y_mwa[NUM_ANTS*NUM_TIME_STEPS] = {0.0002021146 + 0.0000273663*I, 0.0002021146 + 0.0000273663*I, 0.0002021146 + 0.0000273663*I,
                                                           0.0233655237 - 0.0105750393*I, 0.0233655237 - 0.0105750393*I, 0.0233655237 - 0.0105750393*I};

    double _Complex gain2y_mwa[NUM_ANTS*NUM_TIME_STEPS] = {-0.8777521849 + 0.4791148901*I, -0.8777521849 + 0.4791148901*I, -0.8777521849 + 0.4791148901*I,
                                                          -0.0530412868 + 0.0235209484*I, -0.0530412868 + 0.0235209484*I, -0.0530412868 + 0.0235209484*I};



    for (int gain_ind = 0; gain_ind < NUM_ANTS*NUM_TIME_STEPS; gain_ind++) {
      gain1x[gain_ind] = gain1x_mwa[gain_ind];
      leak1x[gain_ind] = leak1x_mwa[gain_ind];
      leak1y[gain_ind] = leak1y_mwa[gain_ind];
      gain1y[gain_ind] = gain1y_mwa[gain_ind];

      gain2x[gain_ind] = gain2x_mwa[gain_ind];
      leak2x[gain_ind] = leak2x_mwa[gain_ind];
      leak2y[gain_ind] = leak2y_mwa[gain_ind];
      gain2y[gain_ind] = gain2y_mwa[gain_ind];
    }

  //We normalise the other EveryBeams to the phase centre, so they're all one here
  } else {

    for (int gain_ind = 0; gain_ind < NUM_ANTS*NUM_TIME_STEPS; gain_ind++) {

      gain1x[gain_ind] = 1.0 + 0.0*I;
      leak1x[gain_ind] = 0.0 + 0.0*I;
      leak1y[gain_ind] = 0.0 + 0.0*I;
      gain1y[gain_ind] = 1.0 + 0.0*I;
      gain2x[gain_ind] = 1.0 + 0.0*I;
      leak2x[gain_ind] = 0.0 + 0.0*I;
      leak2y[gain_ind] = 0.0 + 0.0*I;
      gain2y[gain_ind] = 1.0 + 0.0*I;

    }
  }

  test_comp_phase_centre_allgains_diffants(visibility_set, num_comps, 
                                gain1x, leak1x, leak1y, gain1y,
                                gain2x, leak2x, leak2y, gain2y,
                                NUM_ANTS, woden_settings, TOL);

  free_visi_set_inputs(visibility_set);
  free_visi_set_outputs(visibility_set);

  // // ----Now switch on the autos--------------------------

  woden_settings->do_autos = 1;
  woden_settings->num_autos = NUM_CROSS;
  woden_settings->num_visis = woden_settings->num_cross + woden_settings->num_autos;

  // test_calculate_visibilities frees up the sky model, so we need to remake it
  cropped_sky_models = make_cropped_sky_models(RA0, dec0,
                                               n_points, n_gauss, n_shapes,
                                               num_sources);
  set_azza_para(cropped_sky_models, n_points, n_gauss, n_shapes, num_sources,
                beamtype);

  // printf("We have this many visis %d %d %d\n",woden_settings->num_visis,woden_settings->num_autos,woden_settings->num_cross );
  visibility_set = test_calculate_visibilities(cropped_sky_models,
                                          beam_settings, woden_settings, RA0, dec0,
                                          beam_settings->beamtype);

  test_comp_phase_centre_allgains_diffants(visibility_set, num_comps, 
                                gain1x, leak1x, leak1y, gain1y,
                                gain2x, leak2x, leak2y, gain2y,
                                NUM_ANTS, woden_settings, TOL);

  free_visi_set_inputs(visibility_set);
  free_visi_set_outputs(visibility_set);

  free(gain1x);
  free(leak1x);
  free(leak1y);
  free(gain1y);
  free(gain2x);
  free(leak2x);
  free(leak2y);
  free(gain2y);
  free(beam_settings);
  free(woden_settings->mjds);
  free(woden_settings);
}

//Test with a single SOURCE, single COMPONENT
void test_calculate_visibilities_EveryBeam_OneSource_SinglePoint(int do_gpu,int beamtype,
                                                                 const char *beam_ms_path) {
  int n_points = 1;
  int n_gauss = 0;
  int n_shapes = 0;
  int num_sources = 1;
  test_calculate_visibilities_EveryBeam(n_points, n_gauss, n_shapes,
                                               num_sources, do_gpu, beamtype,
                                               beam_ms_path);
}

void test_calculate_visibilities_EveryBeam_OneSource_SingleGauss(int do_gpu, int beamtype,
                                                                 const char *beam_ms_path) {
  int n_points = 0;
  int n_gauss = 1;
  int n_shapes = 0;
  int num_sources = 1;
  test_calculate_visibilities_EveryBeam(n_points, n_gauss, n_shapes,
                                               num_sources, do_gpu, beamtype,
                                               beam_ms_path);

}

void test_calculate_visibilities_EveryBeam_OneSource_SingleShape(int do_gpu, int beamtype,
                                                                 const char *beam_ms_path) {
  int n_points = 0;
  int n_gauss = 0;
  int n_shapes = 1;
  int num_sources = 1;
  test_calculate_visibilities_EveryBeam(n_points, n_gauss, n_shapes,
                                               num_sources, do_gpu, beamtype,
                                               beam_ms_path);
}

void test_calculate_visibilities_EveryBeam_OneSource_SingleAll(int do_gpu, int beamtype,
                                                                 const char *beam_ms_path) {
  int n_points = 1;
  int n_gauss = 1;
  int n_shapes = 1;
  int num_sources = 1;
  test_calculate_visibilities_EveryBeam(n_points, n_gauss, n_shapes,
                                               num_sources, do_gpu, beamtype,
                                               beam_ms_path);
}


//Test with a three SOURCEs, single COMPONENT
void test_calculate_visibilities_EveryBeam_ThreeSource_SinglePoint(int do_gpu, int beamtype,
                                                                 const char *beam_ms_path) {
  int n_points = 1;
  int n_gauss = 0;
  int n_shapes = 0;
  int num_sources = 3;
  test_calculate_visibilities_EveryBeam(n_points, n_gauss, n_shapes,
                                               num_sources, do_gpu, beamtype,
                                               beam_ms_path);

}

void test_calculate_visibilities_EveryBeam_ThreeSource_SingleGauss(int do_gpu, int beamtype,
                                                                 const char *beam_ms_path) {
  int n_points = 0;
  int n_gauss = 1;
  int n_shapes = 0;
  int num_sources = 3;
  test_calculate_visibilities_EveryBeam(n_points, n_gauss, n_shapes,
                                               num_sources, do_gpu, beamtype,
                                               beam_ms_path);
}

void test_calculate_visibilities_EveryBeam_ThreeSource_SingleShape(int do_gpu, int beamtype,
                                                                 const char *beam_ms_path) {
  int n_points = 0;
  int n_gauss = 0;
  int n_shapes = 1;
  int num_sources = 3;
  test_calculate_visibilities_EveryBeam(n_points, n_gauss, n_shapes,
                                               num_sources, do_gpu, beamtype,
                                               beam_ms_path);
}

void test_calculate_visibilities_EveryBeam_ThreeSource_SingleAll(int do_gpu, int beamtype,
                                                                 const char *beam_ms_path) {
  int n_points = 1;
  int n_gauss = 1;
  int n_shapes = 1;
  int num_sources = 3;
  test_calculate_visibilities_EveryBeam(n_points, n_gauss, n_shapes,
                                               num_sources, do_gpu, beamtype,
                                               beam_ms_path);
}


//Test with three SOURCEs, three COPMONENTs
void test_calculate_visibilities_EveryBeam_ThreeSource_FivePoint(int do_gpu, int beamtype,
                                                                 const char *beam_ms_path) {
  int n_points = 5;
  int n_gauss = 0;
  int n_shapes = 0;
  int num_sources = 3;
  test_calculate_visibilities_EveryBeam(n_points, n_gauss, n_shapes,
                                               num_sources, do_gpu, beamtype,
                                               beam_ms_path);

}

void test_calculate_visibilities_EveryBeam_ThreeSource_FiveGauss(int do_gpu, int beamtype,
                                                                 const char *beam_ms_path) {
  int n_points = 0;
  int n_gauss = 5;
  int n_shapes = 0;
  int num_sources = 3;
  test_calculate_visibilities_EveryBeam(n_points, n_gauss, n_shapes,
                                               num_sources, do_gpu, beamtype,
                                               beam_ms_path);
}

void test_calculate_visibilities_EveryBeam_ThreeSource_FiveShape(int do_gpu, int beamtype,
                                                                 const char *beam_ms_path) {
  int n_points = 0;
  int n_gauss = 0;
  int n_shapes = 5;
  int num_sources = 3;
  test_calculate_visibilities_EveryBeam(n_points, n_gauss, n_shapes,
                                               num_sources, do_gpu, beamtype,
                                               beam_ms_path);
}

void test_calculate_visibilities_EveryBeam_ThreeSource_FiveAll(int do_gpu, int beamtype,
                                                                 const char *beam_ms_path) {
  int n_points = 5;
  int n_gauss = 5;
  int n_shapes = 5;
  int num_sources = 3;
  test_calculate_visibilities_EveryBeam(n_points, n_gauss, n_shapes,
                                               num_sources, do_gpu, beamtype,
                                               beam_ms_path);
}

void profile_lofar_everybeam(int do_gpu, int beamtype, const char *beam_ms_path) {
  int n_points = 200;
  int n_gauss = 200;
  int n_shapes = 10;
  int num_sources = 3;
  test_calculate_visibilities_EveryBeam(n_points, n_gauss, n_shapes,
                                        num_sources, do_gpu, beamtype,
                                        beam_ms_path);
}