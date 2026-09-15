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

#define SKA_RA0 0.00000000
#define SKA_DEC0 -0.52359878
double SKA_mjds[] ={60512.84237269, 60512.96737269};
double SKA_lsts[] ={4.36691696, 5.15446548};

void set_azza_para(source_catalogue_t *cropped_sky_models, int num_time_steps,
                    int n_points, int n_gauss, int n_shapes, int num_sources,
                    int beamtype){

  // components_t *components = malloc(sizeof(components_t));
  // components_t *components;

  double azs[num_time_steps];
  double zas[num_time_steps];
  double para_angles[num_time_steps];
  

  if (beamtype == EB_MWA) {
    for (int time_ind = 0; time_ind < num_time_steps; time_ind++) {
      azs[time_ind] = MWA_azs[time_ind];
      zas[time_ind] = MWA_zas[time_ind];
      para_angles[time_ind] = MWA_para_angles[time_ind];
    }
  }

  int ind;

  if (n_points > 0) {
    for (int source_ind = 0; source_ind < num_sources; source_ind++) {
      for (int comp_ind = 0; comp_ind < n_points; comp_ind++) {
        for (int time_ind = 0; time_ind < num_time_steps; time_ind++) {
          
          ind = comp_ind*num_time_steps + time_ind;
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
        for (int time_ind = 0; time_ind < num_time_steps; time_ind++) {
          
          ind = comp_ind*num_time_steps + time_ind;
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
        for (int time_ind = 0; time_ind < num_time_steps; time_ind++) {
          
          ind = comp_ind*num_time_steps + time_ind;
          cropped_sky_models->sources[source_ind].shape_components.azs[ind] = azs[time_ind];
          cropped_sky_models->sources[source_ind].shape_components.zas[ind] = zas[time_ind];
          cropped_sky_models->sources[source_ind].shape_components.para_angles[ind] = para_angles[time_ind];
        }
      }
    }
  }

}

void set_mjds(woden_settings_t *woden_settings, int beamtype, int num_time_steps){
  
    woden_settings->mjds = malloc(num_time_steps*sizeof(double));

    if (beamtype == EB_MWA) {
      for (int time_ind = 0; time_ind < num_time_steps; time_ind++) {
        woden_settings->mjds[time_ind] = MWA_mjds[time_ind];
      }
    }

    else if (beamtype == EB_LOFAR) {
      for (int time_ind = 0; time_ind < num_time_steps; time_ind++) {
        woden_settings->mjds[time_ind] = LOFAR_mjds[time_ind];
      }
    }

    else {
      for (int time_ind = 0; time_ind < num_time_steps; time_ind++) {
        woden_settings->mjds[time_ind] = SKA_mjds[time_ind];
      }
    }
}

void test_calculate_visibilities_EveryBeam(int n_points, int n_gauss, int n_shapes,
                                           int num_sources, int do_gpu, int beamtype,
                                           const char *beam_ms_path) {


  double dec0;

  if (beamtype == EB_MWA) {
    dec0 = MWA_DEC0;
  } else if (beamtype == EB_LOFAR) {
    dec0 = LOFAR_DEC0;
  } else {
    dec0 = SKA_DEC0;
  }


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
  woden_settings->eb_beam_ra0 = RA0 - DD2R;
  woden_settings->eb_beam_dec0 = dec0;
  woden_settings->normalise_primary_beam = 1;

  if (beamtype == EB_MWA) {
    woden_settings->single_everybeam_station = 1;
  } else {
    woden_settings->single_everybeam_station = 0;
  }


  int delays[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  woden_settings->FEE_ideal_delays = delays;

  set_mjds(woden_settings, beamtype, NUM_TIME_STEPS);

  set_azza_para(cropped_sky_models, NUM_TIME_STEPS,
                n_points, n_gauss, n_shapes, num_sources,
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

  double _Complex *expec_gainx = malloc(NUM_ANTS*NUM_TIME_STEPS*sizeof(double _Complex));
  double _Complex *expec_leakx = malloc(NUM_ANTS*NUM_TIME_STEPS*sizeof(double _Complex));
  double _Complex *expec_leaky = malloc(NUM_ANTS*NUM_TIME_STEPS*sizeof(double _Complex));
  double _Complex *expec_gainy = malloc(NUM_ANTS*NUM_TIME_STEPS*sizeof(double _Complex));

  // -0.8775036335 0.4795700312, 0.0001940064 0.0000344838, 0.0002021146 0.0000273663, -0.8777521849 0.4791148901
  // -0.0692513958 0.0402538478, -0.0234183073 0.0134399701, 0.0233655237 -0.0105750393, -0.0530412868 0.0235209484
  
  //MWA EveryBeam is locked to zenith, so beam values vary depending on where things are pointing
  if (beamtype == EB_MWA) {
    double _Complex gxs_mwa[NUM_ANTS*NUM_TIME_STEPS] = {-0.8775036335 + 0.4795700312*I, -0.8775036335 + 0.4795700312*I, -0.8775036335 + 0.4795700312*I,
                                                          -0.0692513958 + 0.0402538478*I, -0.0692513958 + 0.0402538478*I, -0.0692513958 + 0.0402538478*I};

    double _Complex Dxs_mwa[NUM_ANTS*NUM_TIME_STEPS] = {0.0001940064 + 0.0000344838*I, 0.0001940064 + 0.0000344838*I, 0.0001940064 + 0.0000344838*I,
                                                          -0.0234183073 + 0.0134399701*I, -0.0234183073 + 0.0134399701*I, -0.0234183073 + 0.0134399701*I};

    double _Complex Dys_mwa[NUM_ANTS*NUM_TIME_STEPS] = {0.0002021146 + 0.0000273663*I, 0.0002021146 + 0.0000273663*I, 0.0002021146 + 0.0000273663*I,
                                                           0.0233655237 - 0.0105750393*I, 0.0233655237 - 0.0105750393*I, 0.0233655237 - 0.0105750393*I};

    double _Complex gys_mwa[NUM_ANTS*NUM_TIME_STEPS] = {-0.8777521849 + 0.4791148901*I, -0.8777521849 + 0.4791148901*I, -0.8777521849 + 0.4791148901*I,
                                                          -0.0530412868 + 0.0235209484*I, -0.0530412868 + 0.0235209484*I, -0.0530412868 + 0.0235209484*I};
    for (int gain_ind = 0; gain_ind < NUM_ANTS*NUM_TIME_STEPS; gain_ind++) {
      expec_gainx[gain_ind] = gxs_mwa[gain_ind];
      expec_leakx[gain_ind] = Dxs_mwa[gain_ind];
      expec_leaky[gain_ind] = Dys_mwa[gain_ind];
      expec_gainy[gain_ind] = gys_mwa[gain_ind];

    }

  //We normalise the other EveryBeams to the phase centre, so they're all one here
  } else if (beamtype == EB_LOFAR) {

    double _Complex gxs_lofar[NUM_ANTS*NUM_TIME_STEPS] = {0.7890886973 + -0.0007822585*I, 0.7890645712 + -0.0006363356*I, 0.5740544114 + -0.0004572933*I,
                                                    0.5302808837 + -0.0005594725*I, 0.5304586333 + -0.0009421803*I, 2.2908607380 + -0.0040391421*I};

    double _Complex Dxs_lofar[NUM_ANTS*NUM_TIME_STEPS] = {0.0142538948 + -0.0000307587*I, 0.0142534621 + -0.0000281222*I, 0.0103695733 + -0.0000203572*I,
                                                    0.0104804987 + -0.0000276144*I, 0.0104839998 + -0.0000351838*I, 0.0452766390 + -0.0001513574*I};

    double _Complex Dys_lofar[NUM_ANTS*NUM_TIME_STEPS] = {-0.0116992357 + 0.0000037154*I, -0.0116988765 + 0.0000015522*I, -0.0085110799 + 0.0000010455*I,
                                                    -0.0107561224 + -0.0000095649*I, -0.0107597429 + -0.0000018091*I, -0.0464674726 + -0.0000084174*I};

    double _Complex gys_lofar[NUM_ANTS*NUM_TIME_STEPS] = {0.7938761810 + -0.0008156168*I, 0.7938519139 + -0.0006688077*I, 0.5775372637 + -0.0004808828*I,
                                                    0.5338289758 + -0.0005919146*I, 0.5340078939 + -0.0009771928*I, 2.3061887244 + -0.0041901488*I};


    for (int gain_ind = 0; gain_ind < NUM_ANTS*NUM_TIME_STEPS; gain_ind++) {
      expec_gainx[gain_ind] = gxs_lofar[gain_ind];
      expec_leakx[gain_ind] = Dxs_lofar[gain_ind];
      expec_leaky[gain_ind] = Dys_lofar[gain_ind];
      expec_gainy[gain_ind] = gys_lofar[gain_ind];
    }

  } else {
    double _Complex gxs_oskar[NUM_ANTS*NUM_TIME_STEPS] = {1.0592501418 + -0.0124013387*I, 1.0583441460 + -0.0094384205*I, 1.0609604439 + -0.0006131539*I,
                                                          1.0401339976 + 0.0118783756*I, 1.0391701685 + 0.0135437693*I, 1.0398477908 + 0.0146351146*I};

    double _Complex Dxs_oskar[NUM_ANTS*NUM_TIME_STEPS] = {-0.0083413516 + 0.0000465002*I, -0.0083334832 + 0.0000231529*I, -0.0083540150 + -0.0000464244*I,
                                                          -0.0130903130 + -0.0004361606*I, -0.0130764745 + -0.0004568921*I, -0.0130855069 + -0.0004708017*I};

    double _Complex Dys_oskar[NUM_ANTS*NUM_TIME_STEPS] = {0.0092953280 + -0.0001676934*I, 0.0092875564 + -0.0001416831*I, 0.0093110013 + -0.0000643514*I,
                                                          0.0076253644 + -0.0012037931*I, 0.0076206350 + -0.0011903459*I, 0.0076267774 + -0.0011832643*I};

    double _Complex gys_oskar[NUM_ANTS*NUM_TIME_STEPS] = {1.0593393599 + -0.0124588684*I, 1.0584328219 + -0.0094957071*I, 1.0610502723 + -0.0006697961*I,
                                                          1.0604574958 + 0.0104804413*I, 1.0594770754 + 0.0121800614*I, 1.0601700299 + 0.0132915015*I};

    for (int gain_ind = 0; gain_ind < NUM_ANTS*NUM_TIME_STEPS; gain_ind++) {
      expec_gainx[gain_ind] = gxs_oskar[gain_ind];
      expec_leakx[gain_ind] = Dxs_oskar[gain_ind];
      expec_leaky[gain_ind] = Dys_oskar[gain_ind];
      expec_gainy[gain_ind] = gys_oskar[gain_ind];
    }

  }

  test_comp_phase_centre_allgains_diffants(visibility_set, num_comps, 
                                expec_gainx, expec_leakx,
                                expec_leaky, expec_gainy,
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
  set_azza_para(cropped_sky_models, NUM_TIME_STEPS,
                n_points, n_gauss, n_shapes, num_sources,
                beamtype);

  // printf("We have this many visis %d %d %d\n",woden_settings->num_visis,woden_settings->num_autos,woden_settings->num_cross );
  visibility_set = test_calculate_visibilities(cropped_sky_models,
                                          beam_settings, woden_settings, RA0, dec0,
                                          beam_settings->beamtype);

  test_comp_phase_centre_allgains_diffants(visibility_set, num_comps, 
                                expec_gainx, expec_leakx,
                                expec_leaky, expec_gainy,
                                NUM_ANTS, woden_settings, TOL);

  free_visi_set_inputs(visibility_set);
  free_visi_set_outputs(visibility_set);

  free(expec_gainx);
  free(expec_leakx);
  free(expec_leaky);
  free(expec_gainy);

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