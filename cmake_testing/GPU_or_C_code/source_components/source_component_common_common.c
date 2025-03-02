#include "test_source_component_common.h"
#include "source_component_common_common.h"

//External CUDA code we're linking in
extern void copy_outputs_source_component_common_gpu(int num_of_each_flux_type,
           source_t *mem_chunked_source,
           beam_gains_t *mem_beam_gains,
           woden_settings_t *woden_settings,
           beam_settings_t *beam_settings,
           user_precision_complex_t *gxs, user_precision_complex_t *Dxs,
           user_precision_complex_t *Dys, user_precision_complex_t *gys,
           user_precision_t *extrap_flux_I, user_precision_t *extrap_flux_Q,
           user_precision_t *extrap_flux_U, user_precision_t *extrap_flux_V,
           double *ls, double *ms, double *ns,
           e_component_type comptype);

extern double * malloc_freqs_gpu(int num_extrap_freqs, double *extrap_freqs);
extern void free_freqs_gpu(double *d_extrap_freqs);

extern void free_components_gpu(source_t *mem_chunked_source,
                                  e_component_type comptype);

extern void free_extrapolated_flux_arrays_gpu(components_t *mem_components);

extern void free_beam_gains_gpu(beam_gains_t *mem_beam_gains, e_beamtype beamtype);

extern source_t * copy_chunked_source_to_GPU(source_t *chunked_source);

extern void malloc_extrapolated_flux_arrays_gpu(components_t *mem_components, int num_comps,
                                        int num_freqs);

extern void malloc_beam_gains_gpu(beam_gains_t *d_component_beam_gains,
                                     int beamtype, int num_gains);

extern void calc_lmn_for_components_gpu(components_t *mem_components,
                                        int num_components,
                                        woden_settings_t *woden_settings);

extern void copy_CPU_component_gains_to_GPU_beam_gains(components_t *components,
                                       beam_gains_t *d_beam_gains, int num_gains);

/*
Test that l,m,n and beam values are calculated correctly by `source_component_common`
for a constant Dec=0, latitude=0, and all beam types
There are other tests for the l,m,n coords and beam functions, so no need to
test millions of scenarios here, so stick with Dec=0
*/
void test_source_component_common_ConstantDecChooseBeams(int beamtype, char* mwa_fee_hdf5,
                                                         e_component_type comptype,
                                                         int do_gpu) {

  double TOL;

  //Set up some test condition inputs
  int num_times = 3;
  int num_freqs = 2;

  int num_components = num_powers + num_curves + num_lists;
  int num_beam_values = num_times*num_freqs*num_components;

  user_precision_t *zeroes = calloc(num_components, sizeof(user_precision_t));

  //Keep RA between 0 and 2*pi here but enter RAs that should return
  //negative l values
  double ras[9] = {(3*M_PI)/2, (5*M_PI)/3, (7*M_PI)/4, (11*M_PI)/6,
                   0.0, M_PI/6, M_PI/4, M_PI/3, M_PI/2};

  double ra0 = 0.0*DD2R;
  double dec0 = 0.0*DD2R;

  double *decs = malloc(num_components*sizeof(double));
  for (int i = 0; i < num_components; i++) {
    decs[i] = dec0;
  }

  //Get the settings into a woden_settings_t struct
  woden_settings_t *woden_settings = malloc(sizeof(woden_settings_t));
  woden_settings->num_freqs = num_freqs;
  woden_settings->num_time_steps = num_times;
  woden_settings->ra0 = ra0;
  woden_settings->dec0 = dec0;
  woden_settings->sdec0 = sin(dec0);
  woden_settings->cdec0 = cos(dec0);
  woden_settings->latitude = -0.46606083776035967;
  woden_settings->use_dipamps = 0;
  woden_settings->num_ants = 1;
  woden_settings->do_gpu = do_gpu;
  woden_settings->verbose = 1;

  woden_settings->latitudes = malloc(num_times*sizeof(double));
  woden_settings->mjds = malloc(num_times*sizeof(double));
  for (int i = 0; i < num_times; i++)
  {
    woden_settings->latitudes[i] = -0.46606083776035967;
    woden_settings->mjds[i] = 56428.45056713;
  }

  woden_settings->single_everybeam_station = 1;
  if (beamtype == EB_MWA) {
    
    woden_settings->normalise_primary_beam = 0;
    woden_settings->beam_ms_path = "../../../../test_installation/everybeam/MWA-single-timeslot.ms";
    woden_settings->hdf5_beam_path = mwa_fee_hdf5;
  } 
  
  if (beamtype == EB_LOFAR) {
    woden_settings->beam_ms_path = "../../../../test_installation/everybeam/LOFAR_HBA_MOCK.ms";
    woden_settings->normalise_primary_beam = 1;
  }

  if (beamtype == EB_OSKAR) {
    woden_settings->beam_ms_path = "../../../../test_installation/everybeam/create_OSKAR-SKA_ms/OSKAR-SKA-layout.ms";
    woden_settings->normalise_primary_beam = 1;
  }


  

  woden_settings->beamtype = beamtype;

  beam_settings_t *beam_settings = malloc(sizeof(beam_settings_t));
  beam_settings->beamtype = beamtype;

  double freqs[] = {100e+6, 200e+6};

  double *beam_has = malloc(num_components*num_times*sizeof(double));
  double *beam_decs = malloc(num_components*num_times*sizeof(double));

  //Make ha/dec coords if using the Gaussian beam
  //Need ha/decs for the analytic MWA beam as well
  if (beamtype == GAUSS_BEAM || beamtype == MWA_ANALY) {
    //Stick the Gaussian beam pointed at zenith
    //We're testing at latitude=zero
    beam_settings->gauss_ha = 0.0;
    beam_settings->gauss_sdec = 0.0;
    beam_settings->gauss_cdec = 1.0;
    beam_settings->beam_FWHM_rad = 80.0*DD2R;
    beam_settings->beam_ref_freq = 150e+6;

    for (int component = 0; component < num_components; component++) {
      for (int time_ind = 0; time_ind < num_times; time_ind++) {

        beam_has[component*num_times + time_ind] = ras[component];
        beam_decs[component*num_times + time_ind] = decs[component];

      }
    }
  }

  //If FEE_BEAM, call the C code to interrogate the hdf5 file and set beam
  //things up
  if (beamtype == FEE_BEAM || beamtype == FEE_BEAM_INTERP) {

    int32_t status =  new_fee_beam(mwa_fee_hdf5, &beam_settings->fee_beam);

    if (status != 0) {
      handle_hyperbeam_error(__FILE__, __LINE__, "new_fee_beam");
      // printf("There was an error calling new_fee_beam\n");
    }

    TEST_ASSERT_EQUAL(status, 0);

    

    uint32_t num_freqs_hyper;
    uint32_t *freqs_hz;
    if (beamtype == FEE_BEAM) {
      freqs_hz = malloc(2*sizeof(uint32_t));
      freqs_hz[0] = 150e+6;
      freqs_hz[1] = 150e+6;
      // freqs_hz[0] = 100e+6;
      // freqs_hz[1] = 200e+6;
      num_freqs_hyper = 2;
    } else {
      freqs_hz = malloc(2*sizeof(uint32_t));
      freqs_hz[0] = 100e+6;
      freqs_hz[1] = 200e+6;
      num_freqs_hyper = 2;
    }

    beam_settings->hyper_delays = (uint32_t*)malloc(16*sizeof(uint32_t));

    for (int delay = 0; delay < 16; delay++) {
      beam_settings->hyper_delays[delay] = 0;
    }

    double amps[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    uint32_t num_tiles = 1;
    uint32_t num_amps = 16;
    uint8_t norm_to_zenith = 1;

    if (do_gpu == 1) {
      status = new_gpu_fee_beam(beam_settings->fee_beam,
                               freqs_hz,
                               beam_settings->hyper_delays,
                               amps,
                               num_freqs_hyper,
                               num_tiles,
                               num_amps,
                               norm_to_zenith,
                               &beam_settings->gpu_fee_beam);

    if (status != 0) {
      handle_hyperbeam_error(__FILE__, __LINE__, "new_gpu_fee_beam");
      // printf("There was an error calling new_fee_beam\n");
    }

    TEST_ASSERT_EQUAL(status, 0);
    } else {
      //Because the hyperbeam CPU/GPU codes have different APIs, we need to
      //define the amplitudes here - we passed them into `new_gpu_fee_beam`
      //above
      woden_settings->mwa_dipole_amps = amps;
    }
  }

  if (beam_settings->beamtype == EB_LOFAR || beam_settings->beamtype == EB_MWA || beam_settings->beamtype == EB_OSKAR) {
    int eb_status = 0;
  
    const char *element_response_model;
    bool use_differential_beam = false;
    bool use_channel_frequency = true;
    bool use_local_mwa = true;
  
    if (beam_settings->beamtype == EB_LOFAR) {
      element_response_model = "hamaker";
    } else if (beam_settings->beamtype == EB_MWA) {
      element_response_model = "MWA";
    } else {
      element_response_model = "skala40_wave";
    }
  
    beam_settings->everybeam_telescope =  load_everybeam_telescope(&eb_status, 
                                                woden_settings->beam_ms_path,
                                                element_response_model,
                                                use_differential_beam,
                                                use_channel_frequency,
                                                woden_settings->hdf5_beam_path,
                                                use_local_mwa);
  
    if (eb_status != 0) {
      log_message("WARNING - Something went wrong loading the EveryBeam telescope");
    }

  }


  else if (beamtype == MWA_ANALY) {
    //Zenith pointing is your friend
    int delays[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    woden_settings->FEE_ideal_delays = delays;
  }

  // //Output arrays
  // user_precision_complex_t *gxs = calloc(num_beam_values, sizeof(user_precision_complex_t));
  // user_precision_complex_t *Dxs = calloc(num_beam_values, sizeof(user_precision_complex_t));
  // user_precision_complex_t *Dys = calloc(num_beam_values, sizeof(user_precision_complex_t));
  // user_precision_complex_t *gys = calloc(num_beam_values, sizeof(user_precision_complex_t));

  user_precision_complex_t *gxs = malloc(num_beam_values*sizeof(user_precision_complex_t));
  user_precision_complex_t *Dxs = malloc(num_beam_values*sizeof(user_precision_complex_t));
  user_precision_complex_t *Dys = malloc(num_beam_values*sizeof(user_precision_complex_t));
  user_precision_complex_t *gys = malloc(num_beam_values*sizeof(user_precision_complex_t));

  double *ls = malloc(num_components*sizeof(double));
  double *ms = malloc(num_components*sizeof(double));
  double *ns = malloc(num_components*sizeof(double));

  //Space for outputs
  user_precision_t *extrap_flux_I = malloc(num_freqs*num_components*sizeof(user_precision_t));
  user_precision_t *extrap_flux_Q = malloc(num_freqs*num_components*sizeof(user_precision_t));
  user_precision_t *extrap_flux_U = malloc(num_freqs*num_components*sizeof(user_precision_t));
  user_precision_t *extrap_flux_V = malloc(num_freqs*num_components*sizeof(user_precision_t));

  //Set up the components with our values
  components_t components;

  components.ras = ras;
  components.decs = decs;
  components.azs = azs;
  components.zas = zas;
  components.para_angles = para_angles;
  components.beam_has = beam_has;
  components.beam_decs = beam_decs;

  int *power_comp_inds = malloc(num_powers*sizeof(int));
  int *curve_comp_inds = malloc(num_curves*sizeof(int));
  int *list_comp_inds = malloc(num_lists*sizeof(int));

  for (int i = 0; i < num_powers; i++) {
    power_comp_inds[i] = i;
  }

  for (int i = 0; i < num_curves; i++) {
    curve_comp_inds[i] = num_powers + i;
  }

  for (int i = 0; i < num_lists; i++) {
    list_comp_inds[i] = num_powers + num_curves + i;
  }


  components.power_ref_freqs = ref_freqs;
  components.power_ref_stokesI = ref_stokesI;
  components.power_SIs = ref_power_SIs;
  components.power_comp_inds = power_comp_inds;

  components.curve_ref_freqs = ref_freqs;
  components.curve_ref_stokesI = ref_stokesI;
  components.curve_SIs = ref_curve_SIs;
  components.curve_qs = ref_qs;
  components.curve_comp_inds = curve_comp_inds;

  components.list_freqs = list_freqs;
  components.list_stokesI = list_stokesI;
  components.num_list_values = num_list_values;
  components.list_start_indexes = list_start_indexes;
  components.list_comp_inds = list_comp_inds;

  components.total_num_flux_entires = 0;

  for (int i = 0; i < num_lists; i++) {
    components.total_num_flux_entires += num_list_values[i];
  }

  //Don't need to test values are copied as tested elsewhere, but the function
  //that copies components from CPU to GPU needs these extra fields defined
  //for GAUSSIAN/SHAPELET components
  if (comptype == GAUSSIAN || comptype == SHAPELET) {
    components.pas = zeroes;
    components.majors = zeroes;
    components.minors = zeroes;
  }

  if (comptype == SHAPELET) {
    components.shape_coeffs = zeroes;
    components.n1s = zeroes;
    components.n2s = zeroes;
    components.param_indexes = zeroes;
  }

  components.num_primarybeam_values = num_components*woden_settings->num_freqs*woden_settings->num_time_steps;

    //Polarisation models behbeh
  components.n_stokesV_pol_frac = n_stokesV_pol_frac;
  components.stokesV_pol_fracs = stokesV_pol_fracs;
  components.stokesV_pol_frac_comp_inds = stokesV_pol_frac_comp_inds;

  components.n_stokesV_power = n_stokesV_power;
  components.stokesV_power_ref_flux = ref_stokesV;
  components.stokesV_power_SIs = stokesV_power_SIs;
  components.stokesV_power_comp_inds =stokesV_power_comp_inds;
  
  components.n_stokesV_curve = n_stokesV_curve;
  components.stokesV_curve_ref_flux = ref_stokesV;
  components.stokesV_curve_SIs = stokesV_curve_SIs;
  components.stokesV_curve_qs = stokesV_qs;
  components.stokesV_curve_comp_inds = stokesV_curve_comp_inds;

  components.n_linpol_pol_frac = n_linpol_pol_frac;
  components.linpol_pol_fracs = linpol_pol_fracs;
  components.linpol_pol_frac_comp_inds = linpol_pol_frac_comp_inds;

  components.n_linpol_power = n_linpol_power;
  components.linpol_power_ref_flux = ref_linpol;
  components.linpol_power_SIs = linpol_power_SIs;
  components.linpol_power_comp_inds = linpol_power_comp_inds;
  
  components.n_linpol_curve = n_linpol_curve;
  components.linpol_curve_ref_flux = ref_linpol;
  components.linpol_curve_SIs = linpol_curve_SIs;
  components.linpol_curve_qs = linpol_qs;
  components.linpol_curve_comp_inds = linpol_curve_comp_inds;

  components.n_stokesV_list = n_stokesV_list;
  components.stokesV_list_ref_flux = stokesV_list_ref_flux;
  components.stokesV_list_ref_freqs = stokesV_list_ref_freqs;
  components.stokesV_num_list_values = stokesV_num_list_values;
  components.stokesV_list_start_indexes = stokesV_list_start_indexes;
  components.stokesV_list_comp_inds = stokesV_list_comp_inds;

  components.n_linpol_list = n_linpol_list;
  components.stokesQ_list_ref_flux = stokesQ_list_ref_flux;
  components.stokesQ_list_ref_freqs = stokesQ_list_ref_freqs;
  components.stokesQ_num_list_values = stokesQ_num_list_values;
  components.stokesQ_list_start_indexes = stokesQ_list_start_indexes;
  components.stokesQ_list_comp_inds = stokesQ_list_comp_inds;
  components.stokesU_list_ref_flux = stokesU_list_ref_flux;
  components.stokesU_list_ref_freqs = stokesU_list_ref_freqs;
  components.stokesU_num_list_values = stokesU_num_list_values;
  components.stokesU_list_start_indexes = stokesU_list_start_indexes;
  components.stokesU_list_comp_inds = stokesU_list_comp_inds;

  components.n_linpol_p_list = n_linpol_p_list;
  components.linpol_p_list_ref_flux = linpol_p_list_ref_flux;
  components.linpol_p_list_ref_freqs = linpol_p_list_ref_freqs;
  components.linpol_p_num_list_values = linpol_p_num_list_values;
  components.linpol_p_list_start_indexes = linpol_p_list_start_indexes;
  components.linpol_p_list_comp_inds = linpol_p_list_comp_inds;

  components.n_linpol_angles = components.n_linpol_pol_frac + components.n_linpol_power + components.n_linpol_curve + components.n_linpol_p_list;
  components.intr_pol_angle = intr_pol_angle;
  components.rm_values = rms;
  components.linpol_angle_inds = linpol_angle_inds;


  components.n_stokesV_list_flux_entries = n_stokesV_list_flux_entries;
  components.n_stokesQ_list_flux_entries = n_stokesQ_list_flux_entries;
  components.n_stokesU_list_flux_entries = n_stokesU_list_flux_entries;
  components.n_linpol_p_list_flux_entries = n_linpol_p_list_flux_entries;

  components.do_QUV = 1;

  source_t *chunked_source = (source_t *)malloc(sizeof(source_t));

  int NUM_FLUX_TYPES = 3;
  int num_of_each_flux_type = num_powers;

  if (comptype == POINT) {
    chunked_source->point_components = components;
    chunked_source->n_points = NUM_FLUX_TYPES*num_of_each_flux_type;
    chunked_source->n_point_powers = num_of_each_flux_type;
    chunked_source->n_point_curves = num_of_each_flux_type;
    chunked_source->n_point_lists = num_of_each_flux_type;

    chunked_source->n_gauss = 0;
    chunked_source->n_gauss_lists = 0;
    chunked_source->n_gauss_powers = 0;
    chunked_source->n_gauss_curves = 0;
    chunked_source->n_shapes = 0;
    chunked_source->n_shape_lists = 0;
    chunked_source->n_shape_powers = 0;
    chunked_source->n_shape_curves = 0;
    chunked_source->n_shape_coeffs = 0;
  }
  else if (comptype == GAUSSIAN) {
    chunked_source->gauss_components = components;
    chunked_source->n_gauss = NUM_FLUX_TYPES*num_of_each_flux_type;
    chunked_source->n_gauss_powers = num_of_each_flux_type;
    chunked_source->n_gauss_curves = num_of_each_flux_type;
    chunked_source->n_gauss_lists = num_of_each_flux_type;

    chunked_source->n_points = 0;
    chunked_source->n_point_lists = 0;
    chunked_source->n_point_powers = 0;
    chunked_source->n_point_curves = 0;
    chunked_source->n_shapes = 0;
    chunked_source->n_shape_lists = 0;
    chunked_source->n_shape_powers = 0;
    chunked_source->n_shape_curves = 0;
    chunked_source->n_shape_coeffs = 0;
  }
  else if (comptype == SHAPELET) {
    chunked_source->shape_components = components;
    chunked_source->n_shapes = NUM_FLUX_TYPES*num_of_each_flux_type;
    chunked_source->n_shape_powers = num_of_each_flux_type;
    chunked_source->n_shape_curves = num_of_each_flux_type;
    chunked_source->n_shape_lists = num_of_each_flux_type;
    chunked_source->n_shape_coeffs = num_of_each_flux_type;

    chunked_source->n_points = 0;
    chunked_source->n_point_lists = 0;
    chunked_source->n_point_powers = 0;
    chunked_source->n_point_curves = 0;
    chunked_source->n_gauss = 0;
    chunked_source->n_gauss_lists = 0;
    chunked_source->n_gauss_powers = 0;
    chunked_source->n_gauss_curves = 0;

  }

  woden_settings->do_autos = 0;

  beam_gains_t *mem_beam_gains = malloc(sizeof(beam_gains_t));
  source_t *mem_chunked_source = malloc(sizeof(source_t));
  double *mem_freqs = NULL;
  
  //Visibility set not used in this test, but is an argument to the function
  //It's needed when calculating auto-correlations. Test needs to be 
  //expanded to test auto-correlations
  visibility_set_t *mem_visibility_set = NULL;


  //TODO modify/use this when we do pyuvbeam
  // //the everybeam gains are calculated in python, and transferred across
  // //to the C code. So define some gains here; just make them the index of
  // //the gain.
 
  // components_t *chunk_components;
  // if (comptype == POINT) {
  //   chunk_components = &chunked_source->point_components;
  // }
  // else if (comptype == GAUSSIAN) {
  //   chunk_components = &chunked_source->gauss_components;
  // }
  // else {
  //   chunk_components = &chunked_source->shape_components;
  // }

  // if (beamtype == EB_MWA || beamtype == EB_LOFAR || beamtype == EB_OSKAR) {
  //   chunk_components->gxs = malloc(num_beam_values*sizeof(user_precision_complex_t));
  //   chunk_components->Dxs = malloc(num_beam_values*sizeof(user_precision_complex_t));
  //   chunk_components->Dys = malloc(num_beam_values*sizeof(user_precision_complex_t));
  //   chunk_components->gys = malloc(num_beam_values*sizeof(user_precision_complex_t));

  //   for (int i = 0; i < num_beam_values; i++) {
  //     chunk_components->gxs[i] = i + I*i;
  //     chunk_components->Dxs[i] = i + I*i;
  //     chunk_components->Dys[i] = i + I*i;
  //     chunk_components->gys[i] = i + I*i;
  //   }
  // }

  if (do_gpu == 1) {
    //Run the CUDA code
    components_t *mem_components;
    mem_chunked_source = copy_chunked_source_to_GPU(chunked_source);
    mem_freqs = malloc_freqs_gpu(num_freqs, freqs);

    if (comptype == POINT) {
      mem_components = &mem_chunked_source->point_components;
    }
    else if (comptype == GAUSSIAN) {
      mem_components = &mem_chunked_source->gauss_components;
    }
    else {
      mem_components = &mem_chunked_source->shape_components;
    }
    

    //TODO modify/use this when we do pyuvbeam
    // if (beam_settings->beamtype == EB_LOFAR || beam_settings->beamtype == EB_OSKAR  || beam_settings->beamtype == EB_MWA) {
    //   copy_CPU_component_gains_to_GPU_beam_gains(chunk_components, mem_beam_gains,
    //                              num_beam_values);
    // }

    source_component_common(woden_settings, beam_settings, freqs, mem_freqs,
       chunked_source, mem_chunked_source, mem_beam_gains, comptype,
       mem_visibility_set);

    copy_outputs_source_component_common_gpu(num_of_each_flux_type,
           mem_chunked_source, mem_beam_gains,
           woden_settings, beam_settings,
           gxs, Dxs, Dys, gys,
           extrap_flux_I, extrap_flux_Q, extrap_flux_U, extrap_flux_V,
           ls, ms, ns, comptype);

    
    free_beam_gains_gpu(mem_beam_gains, beam_settings->beamtype);
    free_extrapolated_flux_arrays_gpu(mem_components);
    free_components_gpu(mem_chunked_source, comptype);
    free_freqs_gpu(mem_freqs);

  } else{
    mem_chunked_source = chunked_source;
    mem_freqs = freqs;

    //TODO modify/use this when we do pyuvbeam
    // if (beam_settings->beamtype == EB_LOFAR || beam_settings->beamtype == EB_OSKAR  || beam_settings->beamtype == EB_MWA) {
    //   mem_beam_gains->gxs = chunk_components->gxs;
    //   mem_beam_gains->Dxs = chunk_components->Dxs;
    //   mem_beam_gains->Dys = chunk_components->Dys;
    //   mem_beam_gains->gys = chunk_components->gys;
    // }

    //Aight so this is somewhat hacky. But for the GPU version of hyperbeam,
    //you set the beam frequencies in the call to `new_gpu_fee_beam`. This
    //doesn't happen in the CPU version, and we pass `mem_freqs` to the CPU
    //version. In the test above, and the saved expected values, we set the
    //beam freqs to {150e+6, 150e+6} when beamtype is FEE_BEAM. So here, 
    //we just change mem_freqs to match. Fortunately, we generate the expected
    //extrapolated fluxes based on `mem_freqs`, so just changing this here
    //makes things work.
    if (beamtype == FEE_BEAM){
      mem_freqs[0] = 150e+6;
      mem_freqs[1] = 150e+6;
    }

    source_component_common(woden_settings, beam_settings, freqs, mem_freqs,
       chunked_source, mem_chunked_source, mem_beam_gains, comptype,
       mem_visibility_set);

    copy_outputs_source_component_common_cpu(num_of_each_flux_type,
           mem_chunked_source, mem_beam_gains,
           woden_settings, beam_settings,
           gxs, Dxs, Dys, gys,
           extrap_flux_I, extrap_flux_Q, extrap_flux_U, extrap_flux_V,
           ls, ms, ns, comptype);
  }

  double l_expected[9] = {-1.0, -sqrt(3)/2.0, -sqrt(2)/2.0, -0.5,
                          0.0, 0.5, sqrt(2)/2.0, sqrt(3)/2.0, 1.0};
  double n_expected[9] = {0.0, 0.5, sqrt(2)/2.0, sqrt(3)/2.0,
                         1.0, sqrt(3)/2.0, sqrt(2)/2.0, 0.5, 0.0};

  //Check the l values match expectations
  //Both FLOAT and DOUBLE use 64 bit here so same tolerance
  TOL = 1e-15;

  for (int i = 0; i < num_components; i++) {
    TEST_ASSERT_DOUBLE_WITHIN(TOL, l_expected[i], ls[i]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, ms[i]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, n_expected[i], ns[i]);
  }

  //Depending on beamtype, check results match expectations

  //For the Gaussian beam, the way the l,m coords are set up measns we can
  //analytically predict the values, so go through results and calcluate
  //expected values, and compare to what we got
  int beam_ind = 0;
  if (beamtype == GAUSS_BEAM) {

    #ifdef DOUBLE_PRECISION
      TOL = 1e-12;
    #else
      TOL = 1e-7;
    #endif

    double fwhm_lm = sin(beam_settings->beam_FWHM_rad);
    double beam_ref_freq = 150e+6;

    for (int time_ind = 0; time_ind < num_times; time_ind++) {
      for (int freq_ind = 0; freq_ind < num_freqs; freq_ind++) {
        for (int comp_ind = 0; comp_ind < num_components; comp_ind++) {

          double std = (fwhm_lm / FWHM_FACTOR) * (beam_ref_freq / freqs[freq_ind]);
          double exp_inside = (double)ls[comp_ind] / std;
          double estimate = exp(-0.5*exp_inside*exp_inside);

          //Only real gains have value in GAUSSIAN beam model
          TEST_ASSERT_DOUBLE_WITHIN(TOL, estimate,
                                   creal(gxs[beam_ind]));

          TEST_ASSERT_DOUBLE_WITHIN(TOL, estimate,
                                   creal(gys[beam_ind]));

          //Everything else should be zero
          TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(gys[beam_ind]));
          TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(gxs[beam_ind]));

          //the function get_beam_gains sets the leakage terms to zero
          //further down the line, so no need to test them here
          // TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(Dys[beam_ind]));
          // TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, creal(Dxs[beam_ind]));
          // TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(Dxs[beam_ind]));
          // TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, creal(Dys[beam_ind]));

          beam_ind += 1;

        }
      }
    }
  }
  else if (beamtype == ANALY_DIPOLE) {

    #ifdef DOUBLE_PRECISION
      TOL = 1e-7;
    #else
      TOL = 1e-6;
    #endif

    for (int output = 0; output < num_beam_values; output++) {
      // printf("%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n", creal(gxs[output]), cimag(gxs[output]),
      //         creal(Dxs[output]), cimag(Dxs[output]),
      //         creal(Dys[output]), cimag(Dys[output]),
      //         creal(gys[output]), cimag(gys[output]) );

      TEST_ASSERT_DOUBLE_WITHIN(TOL, analy_expec_J00[output], creal(gxs[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, analy_expec_J11[output], creal(gys[output]));

      TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(gxs[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(gys[output]));

      //the function get_beam_gains sets the leakage terms to zero
      //further down the line, so no need to test them here
      // TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(Dys[beam_ind]));
      // TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, creal(Dxs[beam_ind]));
      // TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(Dxs[beam_ind]));
      // TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, creal(Dys[beam_ind]));
    }
  }
  else if (beamtype == FEE_BEAM) {

    #ifdef DOUBLE_PRECISION
      TOL = 1e-7;
    #else
      TOL = 1e-7;
    #endif

    // printf("num_beam_values %d\n", num_beam_values);
    for (int output = 0; output < num_beam_values; output++) {
      // printf("%d %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n", num_beam_values, creal(gxs[output]), cimag(gxs[output]),
      //         creal(Dxs[output]), cimag(Dxs[output]),
      //         creal(Dys[output]), cimag(Dys[output]),
      //         creal(gys[output]), cimag(gys[output]) );


      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_J00_re[output], creal(gxs[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_J00_im[output], cimag(gxs[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_J01_re[output], creal(Dxs[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_J01_im[output], cimag(Dxs[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_J10_re[output], creal(Dys[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_J10_im[output], cimag(Dys[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_J11_re[output], creal(gys[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_J11_im[output], cimag(gys[output]));
    }

    free_fee_beam(beam_settings->fee_beam);
    if (do_gpu == 1) {
      free_gpu_fee_beam(beam_settings->gpu_fee_beam);
    }

  }
  else if (beamtype == FEE_BEAM_INTERP) {

    #ifdef DOUBLE_PRECISION
      TOL = 1e-7;
    #else
      TOL = 1e-7;
    #endif

    for (int output = 0; output < num_beam_values; output++) {

      // printf("%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n", creal(gxs[output]), cimag(gxs[output]),
      //         creal(Dxs[output]), cimag(Dxs[output]),
      //         creal(Dys[output]), cimag(Dys[output]),
      //         creal(gys[output]), cimag(gys[output]) );

      // printf("reals %d %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n",output, creal(gxs[output]), fee_expec_interp_J00_re[output],
      //                                                     creal(Dxs[output]), fee_expec_interp_J01_re[output],
      //                                                     creal(Dys[output]), fee_expec_interp_J10_re[output],
      //                                                     creal(gxs[output]), fee_expec_interp_J00_re[output]);

      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_interp_J00_re[output], creal(gxs[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_interp_J00_im[output], cimag(gxs[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_interp_J01_re[output], creal(Dxs[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_interp_J01_im[output], cimag(Dxs[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_interp_J10_re[output], creal(Dys[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_interp_J10_im[output], cimag(Dys[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_interp_J11_re[output], creal(gys[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, fee_expec_interp_J11_im[output], cimag(gys[output]));
    }

    free_fee_beam(beam_settings->fee_beam);

    if (do_gpu == 1) {
      free_gpu_fee_beam(beam_settings->gpu_fee_beam);
    }
  }
  else if (beamtype == MWA_ANALY) {

    #ifdef DOUBLE_PRECISION
      TOL = 1e-12;
    #else
      TOL = 1e-7;
    #endif

    //All imaginary should be zero.Turns out for these ra,dec coords, the
    //Dy leakages are essentially zero as well, test a few values equal zero
    for (int output = 0; output < num_beam_values; output++) {

        TEST_ASSERT_DOUBLE_WITHIN(TOL, MWA_analy_expec_J00_re[output],
                                  creal(gxs[output]));

        TEST_ASSERT_DOUBLE_WITHIN(TOL, MWA_analy_expec_J01_re[output],
                                  creal(Dxs[output]));

        TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, creal(Dys[output]));

        TEST_ASSERT_DOUBLE_WITHIN(TOL, MWA_analy_expec_J11_re[output],
                                  creal(gys[output]));

        TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(gxs[output]));
        TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(Dxs[output]));
        TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(Dys[output]));
        TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, cimag(gys[output]));

        // printf("%.12f %.12f %.12f %.12f\n",creal(gxs[output]),
        // creal(Dxs[output]),
        // creal(Dys[output]),
        // creal(gys[output]) );

    }
  }
  else if (beamtype == EB_MWA) {

    #ifdef DOUBLE_PRECISION
      TOL = 1e-7;
    #else
      TOL = 1e-7;
    #endif

    // printf("eb_mwa_gx_real[] = {");
    // for (int output = 0; output < num_beam_values; output++) {
    //   printf("%.8f, ", creal(gxs[output]));
    // }
    // printf("}\n\n");

    // printf("eb_mwa_gx_imag[] = {");
    // for (int output = 0; output < num_beam_values; output++) {
    //   printf("%.8f, ", cimag(gxs[output]));
    // }
    // printf("}\n\n");

    // printf("eb_mwa_Dx_real[] = {");
    // for (int output = 0; output < num_beam_values; output++) {
    //   printf("%.8f, ", creal(Dxs[output]));
    // }
    // printf("}\n\n");

    // printf("eb_mwa_Dx_imag[] = {");
    // for (int output = 0; output < num_beam_values; output++) {
    //   printf("%.8f, ", cimag(Dxs[output]));
    // }
    // printf("}\n\n");

    // printf("eb_mwa_Dy_real[] = {");
    // for (int output = 0; output < num_beam_values; output++) {
    //   printf("%.8f, ", creal(Dys[output]));
    // }
    // printf("}\n\n");

    // printf("eb_mwa_Dy_imag[] = {");
    // for (int output = 0; output < num_beam_values; output++) {
    //   printf("%.8f, ", cimag(Dys[output]));
    // }
    // printf("}\n\n");

    // printf("eb_mwa_gy_real[] = {");
    // for (int output = 0; output < num_beam_values; output++) {
    //   printf("%.8f, ", creal(gys[output]));
    // }
    // printf("}\n\n");

    // printf("eb_mwa_gy_imag[] = {");
    // for (int output = 0; output < num_beam_values; output++) {
    //   printf("%.8f, ", cimag(gys[output]));
    // }
    // printf("}\n\n");


    for (int output = 0; output < num_beam_values; output++) {
      // printf("%d %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n", output, creal(gxs[output]), cimag(gxs[output]),
      //         creal(Dxs[output]), cimag(Dxs[output]),
      //         creal(Dys[output]), cimag(Dys[output]),
      //         creal(gys[output]), cimag(gys[output]) );
      //All values should just equal the index in the array


      TEST_ASSERT_DOUBLE_WITHIN(TOL, eb_mwa_gx_real[output], creal(gxs[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, eb_mwa_gx_imag[output], cimag(gxs[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, eb_mwa_Dx_real[output], creal(Dxs[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, eb_mwa_Dx_imag[output], cimag(Dxs[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, eb_mwa_Dy_real[output], creal(Dys[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, eb_mwa_Dy_imag[output], cimag(Dys[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, eb_mwa_gy_real[output], creal(gys[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, eb_mwa_gy_imag[output], cimag(gys[output]));


      //TODO modify/use this when we do pyuvbeam
      // TEST_ASSERT_DOUBLE_WITHIN(TOL, output, creal(gxs[output]));
      // TEST_ASSERT_DOUBLE_WITHIN(TOL, output, cimag(gxs[output]));
      // TEST_ASSERT_DOUBLE_WITHIN(TOL, output, creal(Dxs[output]));
      // TEST_ASSERT_DOUBLE_WITHIN(TOL, output, cimag(Dxs[output]));
      // TEST_ASSERT_DOUBLE_WITHIN(TOL, output, creal(Dys[output]));
      // TEST_ASSERT_DOUBLE_WITHIN(TOL, output, cimag(Dys[output]));
      // TEST_ASSERT_DOUBLE_WITHIN(TOL, output, creal(gys[output]));
      // TEST_ASSERT_DOUBLE_WITHIN(TOL, output, cimag(gys[output]));
    }
  }

  else if (beamtype == EB_LOFAR) {

    #ifdef DOUBLE_PRECISION
      TOL = 1e-7;
    #else
      TOL = 1e-6;
    #endif

    // printf("double eb_lofar_gx_real[] = {");
    // for (int output = 0; output < num_beam_values; output++) {
    //   printf("%.8f, ", creal(gxs[output]));
    // }
    // printf("}\n\n");

    // printf("double eb_lofar_gx_imag[] = {");
    // for (int output = 0; output < num_beam_values; output++) {
    //   printf("%.8f, ", cimag(gxs[output]));
    // }
    // printf("}\n\n");

    // printf("double eb_lofar_Dx_real[] = {");
    // for (int output = 0; output < num_beam_values; output++) {
    //   printf("%.8f, ", creal(Dxs[output]));
    // }
    // printf("}\n\n");

    // printf("double eb_lofar_Dx_imag[] = {");
    // for (int output = 0; output < num_beam_values; output++) {
    //   printf("%.8f, ", cimag(Dxs[output]));
    // }
    // printf("}\n\n");

    // printf("double eb_lofar_Dy_real[] = {");
    // for (int output = 0; output < num_beam_values; output++) {
    //   printf("%.8f, ", creal(Dys[output]));
    // }
    // printf("}\n\n");

    // printf("double eb_lofar_Dy_imag[] = {");
    // for (int output = 0; output < num_beam_values; output++) {
    //   printf("%.8f, ", cimag(Dys[output]));
    // }
    // printf("}\n\n");

    // printf("double eb_lofar_gy_real[] = {");
    // for (int output = 0; output < num_beam_values; output++) {
    //   printf("%.8f, ", creal(gys[output]));
    // }
    // printf("}\n\n");

    // printf("double eb_lofar_gy_imag[] = {");
    // for (int output = 0; output < num_beam_values; output++) {
    //   printf("%.8f, ", cimag(gys[output]));
    // }
    // printf("}\n\n");


    for (int output = 0; output < num_beam_values; output++) {
      // printf("%d %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n", output, creal(gxs[output]), cimag(gxs[output]),
      //         creal(Dxs[output]), cimag(Dxs[output]),
      //         creal(Dys[output]), cimag(Dys[output]),
      //         creal(gys[output]), cimag(gys[output]) );
      //All values should just equal the index in the array


      TEST_ASSERT_DOUBLE_WITHIN(TOL, eb_lofar_gx_real[output], creal(gxs[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, eb_lofar_gx_imag[output], cimag(gxs[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, eb_lofar_Dx_real[output], creal(Dxs[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, eb_lofar_Dx_imag[output], cimag(Dxs[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, eb_lofar_Dy_real[output], creal(Dys[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, eb_lofar_Dy_imag[output], cimag(Dys[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, eb_lofar_gy_real[output], creal(gys[output]));
      TEST_ASSERT_DOUBLE_WITHIN(TOL, eb_lofar_gy_imag[output], cimag(gys[output]));

    }
  }


  else if (beamtype == NO_BEAM) {
    //Don't need to check beam values for NO_BEAM, these values are set by
    //the function `source_components::get_beam_gains` which is called later
    //by the COMPONENT kernels assigns gains on one and leakage of zero
    //in the NO_BEAM case.
    ;
  }

  //Now check the fluxes were extrapolated correctly
  //Make some expected value arrays
  double *expec_flux_I = malloc(num_freqs*num_components*sizeof(double));
  double *expec_flux_Q = malloc(num_freqs*num_components*sizeof(double));
  double *expec_flux_U = malloc(num_freqs*num_components*sizeof(double));
  double *expec_flux_V = malloc(num_freqs*num_components*sizeof(double));

  CPU_extrapolate_fluxes_in_components(&components, num_powers, num_curves, num_lists,
                        freqs, num_freqs,
                        expec_flux_I, expec_flux_Q, expec_flux_U, expec_flux_V);

  #ifdef DOUBLE_PRECISION
    TOL = 1e-11;
  #else
    TOL = 4e-4;
  #endif

  for (int i = 0; i < num_freqs*(num_powers + num_curves + num_lists); i++) {
    //Check the two are within tolerace
    // printf("%d %.3f %.3f\n",i, expec_flux_I[i], extrap_flux_I[i] );
    TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_flux_I[i], extrap_flux_I[i]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_flux_Q[i], extrap_flux_Q[i]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_flux_U[i], extrap_flux_U[i]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_flux_V[i], extrap_flux_V[i]);

  }

  if (beam_settings->beamtype == EB_LOFAR || beam_settings->beamtype == EB_MWA || beam_settings->beamtype == EB_OSKAR) {
    destroy_everybeam_telescope(beam_settings->everybeam_telescope);
  }

  //TODO modify/use this when we do pyuvbeam
  // //Be free my beauties
  // if (beamtype == EB_MWA || beamtype == EB_LOFAR || beamtype == EB_LOFAR) {
  //   free(chunk_components->gxs);
  //   free(chunk_components->Dxs);
  //   free(chunk_components->Dys);
  //   free(chunk_components->gys);
  // }
  free(zeroes);
  free(decs);
  free(gxs);
  free(Dxs);
  free(Dys);
  free(gys);
  free(ls);
  free(ms);
  free(ns);
  free(beam_has);
  free(beam_decs);
  free(power_comp_inds);
  free(curve_comp_inds);
  free(list_comp_inds);
  free(woden_settings->latitudes);
  free(woden_settings->mjds);

}

/*
This test checks source_component_common with beamtype=FEE_BEAM, for a given
comptype
*/
void test_source_component_common_ConstantDecFEEBeam(e_component_type comptype,
                                                     int do_gpu) {
  //Look for environment variable MWA_FEE_HDF5, which should point
  //towards mwa_full_embedded_element_pattern.h5. If we can't find it,
  //can't run the tests, so don't
  char* mwa_fee_hdf5 = getenv("MWA_FEE_HDF5");

  if (mwa_fee_hdf5) {
    printf("MWA_FEE_HDF5: %s\n", mwa_fee_hdf5 );
    test_source_component_common_ConstantDecChooseBeams(FEE_BEAM, mwa_fee_hdf5,
                                                        comptype, do_gpu);
  }
  else {
    printf("MWA_FEE_HDF5 not found - not running MWA_FEE beam test\n");
  }
}

/*
This test checks source_component_common with beamtype=FEE_BEAM_INTERP
*/
void test_source_component_common_ConstantDecFEEBeamInterp(e_component_type comptype,
                                                           int do_gpu) {
  //Look for environment variable MWA_FEE_HDF5, which should point
  //towards mwa_full_embedded_element_pattern.h5. If we can't find it,
  //can't run the tests, so don't
  char* mwa_fee_hdf5 = getenv("MWA_FEE_HDF5_INTERP");

  if (mwa_fee_hdf5) {
    printf("MWA_FEE_HDF5_INTERP: %s\n", mwa_fee_hdf5 );
    test_source_component_common_ConstantDecChooseBeams(FEE_BEAM_INTERP, mwa_fee_hdf5,
                                                        comptype, do_gpu);
  }
  else {
    printf("MWA_FEE_HDF5_INTERP not found - not running MWA_FEE beam test\n");
  }
}

/*
This test checks source_component_common with beamtype=ANALY_DIPOLE
*/
void test_source_component_common_ConstantDecAnalyBeam(e_component_type comptype,
                                                       int do_gpu) {
  test_source_component_common_ConstantDecChooseBeams(ANALY_DIPOLE, " ",
                                                      comptype, do_gpu);
}

/*
This test checks source_component_common with beamtype=GAUSS_BEAM
*/
void test_source_component_common_ConstantDecGaussBeam(e_component_type comptype,
                                                       int do_gpu) {
  test_source_component_common_ConstantDecChooseBeams(GAUSS_BEAM, " ",
                                                      comptype, do_gpu);
}

/*
This test checks source_component_common with beamtype=NO_BEAM
*/
void test_source_component_common_ConstantDecNoBeam(e_component_type comptype,
                                                       int do_gpu) {
  test_source_component_common_ConstantDecChooseBeams(NO_BEAM, " ",
                                                      comptype, do_gpu);
}

/*
This test checks source_component_common with beamtype=MWA_ANALY
*/
void test_source_component_common_ConstantDecMWAAnaly(e_component_type comptype,
                                                       int do_gpu) {
  test_source_component_common_ConstantDecChooseBeams(MWA_ANALY, " ",
                                                      comptype, do_gpu);
}

/*
This test checks source_component_common with beamtype=EB_MWA
*/
void test_source_component_common_ConstantDecEveryBeamMWA(e_component_type comptype,
                                                       int do_gpu) {

  char* mwa_fee_hdf5 = getenv("MWA_FEE_HDF5");

  if (mwa_fee_hdf5) {
    printf("MWA_FEE_HDF5: %s\n", mwa_fee_hdf5 );
    test_source_component_common_ConstantDecChooseBeams(EB_MWA, mwa_fee_hdf5,
                                                        comptype, do_gpu);
  }
  else {
    printf("MWA_FEE_HDF5 not found - not running EB MWA beam test\n");
  }
}


/*
This test checks source_component_common with beamtype=EB_LOFAR
*/
void test_source_component_common_ConstantDecEveryBeamLOFAR(e_component_type comptype,
                                                       int do_gpu) {
  test_source_component_common_ConstantDecChooseBeams(EB_LOFAR, " ",
                                                      comptype, do_gpu);
}

/*
This test checks source_component_common with beamtype=EB_MWA
*/
void test_source_component_common_ConstantDecEveryBeamOSKAR(e_component_type comptype,
                                                       int do_gpu) {
  test_source_component_common_ConstantDecChooseBeams(EB_OSKAR, " ",
                                                      comptype, do_gpu);
}