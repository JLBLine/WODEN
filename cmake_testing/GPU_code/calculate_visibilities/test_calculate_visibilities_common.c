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

#include <math.h>
#include <unity.h>
#include <stdlib.h>
#include <complex.h>
#include "test_calculate_visibilities_common.h"
#include "hyperbeam_error.h"

// //External CUDA code we're linking in
// extern void calculate_visibilities(array_layout_t *array_layout,
//   source_catalogue_t *cropped_sky_models, beam_settings_t *beam_settings,
//   woden_settings_t *woden_settings, visibility_set_t *visibility_set,
//   user_precision_t *sbf);

//Just two LSTs to check things change with time
double lsts[] = {0.0, M_PI / 4};

user_precision_t azs[] = {0.0, 4.528359553989764};
user_precision_t zas[] = {0.0, 0.6978088917603547};

//Different delays settings, which control the comping of the MWA beam
user_precision_t zenith_delays[16] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};


void populate_components(components_t *comps, int n_comps,
                         double ra0, double dec0){

  comps->ras = malloc(n_comps*sizeof(double));
  comps->decs = malloc(n_comps*sizeof(double));

  comps->num_primarybeam_values = n_comps*NUM_FREQS*NUM_TIME_STEPS;
  comps->beam_has = malloc(n_comps*NUM_TIME_STEPS*sizeof(double));
  comps->beam_decs = malloc(n_comps*NUM_TIME_STEPS*sizeof(double));
  comps->azs = malloc(n_comps*NUM_TIME_STEPS*sizeof(user_precision_t));
  comps->zas = malloc(n_comps*NUM_TIME_STEPS*sizeof(user_precision_t));

  comps->majors = malloc(n_comps*sizeof(user_precision_t));
  comps->minors = malloc(n_comps*sizeof(user_precision_t));
  comps->pas = malloc(n_comps*sizeof(user_precision_t));

  comps->shape_coeffs = malloc(n_comps*sizeof(user_precision_t));
  comps->n1s = malloc(n_comps*sizeof(user_precision_t));
  comps->n2s = malloc(n_comps*sizeof(user_precision_t));
  comps->param_indexes = malloc(n_comps*sizeof(user_precision_t));


  //No matter what, always have one POWER_LAW source
  comps->power_ref_freqs = malloc(sizeof(double));
  comps->power_ref_stokesI = malloc(sizeof(user_precision_t));
  // comps->power_ref_stokesQ = malloc(sizeof(user_precision_t));
  // comps->power_ref_stokesU = malloc(sizeof(user_precision_t));
  // comps->power_ref_stokesV = malloc(sizeof(user_precision_t));
  comps->power_SIs = malloc(sizeof(user_precision_t));
  comps->power_comp_inds = malloc(sizeof(int));

  comps->power_ref_freqs[0] = REF_FREQ;
  comps->power_ref_stokesI[0] = STOKESI;
  // comps->power_ref_stokesQ[0] = 0.0;
  // comps->power_ref_stokesU[0] = 0.0;
  // comps->power_ref_stokesV[0] = 0.0;
  comps->power_SIs[0] = 0.0;
  comps->power_comp_inds[0] = 0;

  //set everything to zero first, and update if necessary
  comps->n_stokesV_power = 0;
  comps->n_stokesV_curve = 0;
  comps->n_stokesV_list = 0;
  comps->n_linpol_power = 0;
  comps->n_linpol_curve = 0;
  comps->n_linpol_list = 0;
  comps->n_linpol_p_list = 0;
  comps->n_linpol_angles = 0;
  
  //No matter what, always have one POL_FRACTION source for both pols
  comps->do_QUV = 1;
  comps->stokesV_pol_fracs = malloc(sizeof(user_precision_t));
  comps->stokesV_pol_frac_comp_inds = malloc(sizeof(int));
  comps->linpol_pol_fracs = malloc(sizeof(user_precision_t));
  comps->linpol_pol_frac_comp_inds = malloc(sizeof(int));

  comps->n_stokesV_pol_frac = 1;
  comps->stokesV_pol_fracs[0] = POL_FRAC;
  comps->stokesV_pol_frac_comp_inds[0] = 0;
  comps->n_linpol_pol_frac = 1;
  comps->linpol_pol_fracs[0] = POL_FRAC;
  comps->linpol_pol_frac_comp_inds[0] = 0;


  //This is techincally illegal, but I wrote the code so I'll allow it;
  //here we give everything a linpol angle. This means it'll grab whatever is
  //in Stokes Q, and use that to calculate Stokes Q and U via the RM. As we've
  //set RM to zero, should always leave Stokes Q as, and set U to zero. This
  //makes predicting the visibilities easier, and we have tests elsewhere to
  //make sure the RM stuff works
  comps->n_linpol_angles = n_comps;
  comps->intr_pol_angle = malloc(n_comps*sizeof(user_precision_t));
  comps->rm_values = malloc(n_comps*sizeof(user_precision_t));
  comps->linpol_angle_inds = malloc(n_comps*sizeof(int));

  for (int i = 0; i < n_comps; i++)
  {
    comps->intr_pol_angle[i] = 0;
    comps->rm_values[i] = 0;
    comps->linpol_angle_inds[i] = i;
  }
  

  if (n_comps > 1) {
    comps->curve_ref_freqs = malloc(2*sizeof(double));
    comps->curve_ref_stokesI = malloc(2*sizeof(user_precision_t));
    comps->curve_SIs = malloc(2*sizeof(user_precision_t));
    comps->curve_qs = malloc(2*sizeof(user_precision_t));
    comps->curve_comp_inds = malloc(2*sizeof(int));

    comps->curve_ref_freqs[0] = REF_FREQ;
    comps->curve_ref_freqs[1] = REF_FREQ;
    comps->curve_ref_stokesI[0] = STOKESI;
    comps->curve_ref_stokesI[1] = STOKESI;
    comps->curve_SIs[0] = 0.0;
    comps->curve_SIs[1] = 0.0;
    comps->curve_qs[0] = 0.0;
    comps->curve_qs[1] = 0.0;
    comps->curve_comp_inds[0] = 1;
    comps->curve_comp_inds[1] = 2;

    comps->list_freqs = malloc(4*sizeof(double));
    comps->list_stokesI = malloc(4*sizeof(user_precision_t));
    comps->list_comp_inds = malloc(2*sizeof(int));
    comps->num_list_values = malloc(2*sizeof(int));
    comps->list_start_indexes = malloc(2*sizeof(int));
    comps->total_num_flux_entires = 4;

    comps->list_freqs[0] = 150e+6;
    comps->list_stokesI[0] = STOKESI;
    comps->list_freqs[1] = 170e+6;
    comps->list_stokesI[1] = STOKESI;
    comps->list_freqs[2] = 150e+6;
    comps->list_stokesI[2] = STOKESI;
    comps->list_freqs[3] = 170e+6;
    comps->list_stokesI[3] = STOKESI;
    comps->list_comp_inds[0] = 3;
    comps->list_comp_inds[1] = 4;
    comps->num_list_values[0] = 2;
    comps->num_list_values[1] = 2;
    comps->list_start_indexes[0] = 0;
    comps->list_start_indexes[1] = 2;

    comps->n_stokesV_power = 1;
    comps->stokesV_power_ref_flux = malloc(sizeof(user_precision_t));
    comps->stokesV_power_SIs = malloc(sizeof(user_precision_t));
    comps->stokesV_power_comp_inds = malloc(sizeof(int));
    comps->stokesV_power_ref_flux[0] = STOKESI;
    comps->stokesV_power_SIs[0] = 0.0;
    comps->stokesV_power_comp_inds[0] = 1;

    comps->n_stokesV_curve = 1;
    comps->stokesV_curve_ref_flux = malloc(sizeof(user_precision_t));
    comps->stokesV_curve_SIs = malloc(sizeof(user_precision_t));
    comps->stokesV_curve_qs = malloc(sizeof(user_precision_t));
    comps->stokesV_curve_comp_inds = malloc(sizeof(int));
    comps->stokesV_curve_ref_flux[0] = STOKESI;
    comps->stokesV_curve_SIs[0] = 0.0;
    comps->stokesV_curve_qs[0] = 0.0;
    comps->stokesV_curve_comp_inds[0] = 2;
    
    comps->n_linpol_power = 1;
    comps->linpol_power_ref_flux = malloc(sizeof(user_precision_t));
    comps->linpol_power_SIs = malloc(sizeof(user_precision_t));
    comps->linpol_power_comp_inds = malloc(sizeof(int));
    comps->linpol_power_ref_flux[0] = STOKESI;
    comps->linpol_power_SIs[0] = 0.0;
    comps->linpol_power_comp_inds[0] = 1;

    comps->n_linpol_curve = 1;
    comps->linpol_curve_ref_flux = malloc(sizeof(user_precision_t));
    comps->linpol_curve_SIs = malloc(sizeof(user_precision_t));
    comps->linpol_curve_qs = malloc(sizeof(user_precision_t));
    comps->linpol_curve_comp_inds = malloc(sizeof(int));
    comps->linpol_curve_ref_flux[0] = STOKESI;
    comps->linpol_curve_SIs[0] = 0.0;
    comps->linpol_curve_qs[0] = 0.0;
    comps->linpol_curve_comp_inds[0] = 2;

    comps->n_stokesV_list = 1;
    comps->n_stokesV_list_flux_entries = 2;
    comps->stokesV_list_ref_freqs = malloc(2*sizeof(double));
    comps->stokesV_list_ref_flux = malloc(2*sizeof(user_precision_t));
    comps->stokesV_list_comp_inds = malloc(sizeof(int));
    comps->stokesV_num_list_values = malloc(sizeof(int));
    comps->stokesV_list_start_indexes = malloc(sizeof(int));
    comps->stokesV_list_ref_freqs[0] = 150e+6;
    comps->stokesV_list_ref_flux[0] = STOKESI;
    comps->stokesV_list_ref_freqs[1] = 170e+6;
    comps->stokesV_list_ref_flux[1] = STOKESI;
    comps->stokesV_list_comp_inds[0] = 3;
    comps->stokesV_num_list_values[0] = 2;
    comps->stokesV_list_start_indexes[0] = 0;

    comps->n_stokesV_pol_frac = 2;
    comps->stokesV_pol_fracs[1] = POL_FRAC;
    comps->stokesV_pol_frac_comp_inds[1] = 4;

    comps->n_linpol_list = 1;
    comps->n_stokesQ_list_flux_entries = 2;
    comps->stokesQ_list_ref_freqs = malloc(2*sizeof(double));
    comps->stokesQ_list_ref_flux = malloc(2*sizeof(user_precision_t));
    comps->stokesQ_list_comp_inds = malloc(sizeof(int));
    comps->stokesQ_num_list_values = malloc(sizeof(int));
    comps->stokesQ_list_start_indexes = malloc(sizeof(int));
    comps->stokesQ_list_ref_freqs[0] = 150e+6;
    comps->stokesQ_list_ref_flux[0] = STOKESI;
    comps->stokesQ_list_ref_freqs[1] = 170e+6;
    comps->stokesQ_list_ref_flux[1] = STOKESI;
    comps->stokesQ_list_comp_inds[0] = 3;
    comps->stokesQ_num_list_values[0] = 2;
    comps->stokesQ_list_start_indexes[0] = 0;
    comps->n_stokesU_list_flux_entries = 2;
    comps->stokesU_list_ref_freqs = malloc(2*sizeof(double));
    comps->stokesU_list_ref_flux = malloc(2*sizeof(user_precision_t));
    comps->stokesU_list_comp_inds = malloc(sizeof(int));
    comps->stokesU_num_list_values = malloc(sizeof(int));
    comps->stokesU_list_start_indexes = malloc(sizeof(int));
    comps->stokesU_list_ref_freqs[0] = 150e+6;
    comps->stokesU_list_ref_flux[0] = STOKESI;
    comps->stokesU_list_ref_freqs[1] = 170e+6;
    comps->stokesU_list_ref_flux[1] = STOKESI;
    comps->stokesU_list_comp_inds[0] = 3;
    comps->stokesU_num_list_values[0] = 2;
    comps->stokesU_list_start_indexes[0] = 0;

    comps->n_linpol_p_list = 1;
    comps->n_linpol_p_list_flux_entries = 2;
    comps->linpol_p_list_ref_freqs = malloc(2*sizeof(double));
    comps->linpol_p_list_ref_flux = malloc(2*sizeof(user_precision_t));
    comps->linpol_p_list_comp_inds = malloc(sizeof(int));
    comps->linpol_p_num_list_values = malloc(sizeof(int));
    comps->linpol_p_list_start_indexes = malloc(sizeof(int));
    comps->linpol_p_list_ref_freqs[0] = 150e+6;
    comps->linpol_p_list_ref_flux[0] = STOKESI;
    comps->linpol_p_list_ref_freqs[1] = 170e+6;
    comps->linpol_p_list_ref_flux[1] = STOKESI;
    comps->linpol_p_list_comp_inds[0] = 4;
    comps->linpol_p_num_list_values[0] = 2;
    comps->linpol_p_list_start_indexes[0] = 0;


  } 

  for (int comp = 0; comp < n_comps; comp++) {
    comps->ras[comp] = ra0;
    comps->decs[comp] = dec0;

    comps->majors[comp] = 1e-10;
    comps->minors[comp] = 1e-10;
    comps->pas[comp] = 0;

    comps->shape_coeffs[comp] = 1.0;
    comps->n1s[comp] = 0.0;
    comps->n2s[comp] = 0.0;
    comps->param_indexes[comp] = 0.0;

    for (int time = 0; time < NUM_TIME_STEPS; time++) {
      int step = comp*NUM_TIME_STEPS + time;
      comps->beam_has[step] = lsts[time] - RA0;
      comps->beam_decs[step] = -0.46606083776035967;

      comps->azs[step] = azs[time];
      comps->zas[step] = zas[time];
    }
  }

}


/*
Create a number of SOURCEs and input into a sky model catalogue `source_catalogue_t`
struct. For each SOURCE, populate with as many COMPONENTs as requested of
whatever combination of comp, GAUSSIAN, and SHAPELET types
Keep everything to just Stokes I 1 Jy, and stick the source at phase centre

We'll stick with just one shapelet coeff per shapelet compoenent as there
are other tests to make sure SHAPELETs works elsewhere
*/
source_catalogue_t * make_cropped_sky_models(double ra0, double dec0,
                                             int n_points, int n_gauss,
                                             int n_shapes,
                                             int num_sources) {

  int n_comps = n_points + n_gauss + n_shapes;

  source_catalogue_t *cropped_sky_models = malloc(sizeof(cropped_sky_models));
  cropped_sky_models->num_sources = num_sources;
  cropped_sky_models->num_shapelets = n_shapes*n_comps;
  cropped_sky_models->sources = malloc(num_sources*sizeof(source_t));

  for (int cats_ind = 0; cats_ind < num_sources; cats_ind++) {
    cropped_sky_models->sources[cats_ind].n_points = n_points;
    cropped_sky_models->sources[cats_ind].n_gauss = n_gauss;
    cropped_sky_models->sources[cats_ind].n_shapes = n_shapes;
    cropped_sky_models->sources[cats_ind].n_shape_coeffs = n_shapes;

    cropped_sky_models->sources[cats_ind].n_point_powers = 0;
    cropped_sky_models->sources[cats_ind].n_point_curves = 0;
    cropped_sky_models->sources[cats_ind].n_point_lists = 0;
    cropped_sky_models->sources[cats_ind].n_gauss_powers = 0;
    cropped_sky_models->sources[cats_ind].n_gauss_curves = 0;
    cropped_sky_models->sources[cats_ind].n_gauss_lists = 0;
    cropped_sky_models->sources[cats_ind].n_shape_powers = 0;
    cropped_sky_models->sources[cats_ind].n_shape_curves = 0;
    cropped_sky_models->sources[cats_ind].n_shape_lists = 0;

    if (n_points > 0) {
      populate_components(&cropped_sky_models->sources[cats_ind].point_components,
                          n_points, ra0, dec0);
      cropped_sky_models->sources[cats_ind].n_point_powers = 1;
      if (n_points > 1) {
        cropped_sky_models->sources[cats_ind].n_point_curves = 2;
        cropped_sky_models->sources[cats_ind].n_point_lists = 2;
      }
    }

    if (n_gauss > 0) {
      populate_components(&cropped_sky_models->sources[cats_ind].gauss_components,
                          n_gauss, ra0, dec0);
      cropped_sky_models->sources[cats_ind].n_gauss_powers = 1;
      if (n_gauss > 1) {
        cropped_sky_models->sources[cats_ind].n_gauss_curves = 2;
        cropped_sky_models->sources[cats_ind].n_gauss_lists = 2;
      }
    }

    if (n_shapes > 0) {
      populate_components(&cropped_sky_models->sources[cats_ind].shape_components,
                          n_shapes, ra0, dec0);
      cropped_sky_models->sources[cats_ind].n_shape_powers = 1;
      if (n_shapes > 1) {
        cropped_sky_models->sources[cats_ind].n_shape_curves = 2;
        cropped_sky_models->sources[cats_ind].n_shape_lists = 2;
      }
    }
  }
  return cropped_sky_models;
}

void free_components(components_t comps, int num_comps) {
  free(comps.ras);
  free(comps.decs);

  free(comps.beam_has);
  free(comps.beam_decs);
  free(comps.azs);
  free(comps.zas);

  free(comps.majors);
  free(comps.minors);
  free(comps.pas);

  free(comps.shape_coeffs);
  free(comps.n1s);
  free(comps.n2s);
  free(comps.param_indexes);

  free(comps.power_ref_freqs);
  free(comps.power_ref_stokesI);
  free(comps.power_SIs);
  free(comps.power_comp_inds);

  free(comps.stokesV_pol_fracs);
  free(comps.stokesV_pol_frac_comp_inds);
  free(comps.linpol_pol_fracs);
  free(comps.linpol_pol_frac_comp_inds);

  if (num_comps > 1) {
    free(comps.curve_ref_freqs);
    free(comps.curve_ref_stokesI);
    free(comps.curve_SIs);
    free(comps.curve_qs);
    free(comps.curve_comp_inds);

    free(comps.list_freqs);
    free(comps.list_stokesI);
    free(comps.list_comp_inds);
    free(comps.num_list_values);
    free(comps.list_start_indexes);

    free(comps.stokesV_power_ref_flux);
    free(comps.stokesV_power_SIs);
    free(comps.stokesV_power_comp_inds);
    free(comps.stokesV_curve_ref_flux);
    free(comps.stokesV_curve_SIs);
    free(comps.stokesV_curve_qs);
    free(comps.stokesV_curve_comp_inds);
    free(comps.linpol_power_ref_flux);
    free(comps.linpol_power_SIs);
    free(comps.linpol_power_comp_inds);
    free(comps.linpol_curve_ref_flux);
    free(comps.linpol_curve_SIs);
    free(comps.linpol_curve_qs);
    free(comps.linpol_curve_comp_inds);

    free(comps.stokesV_num_list_values);
    free(comps.stokesV_list_start_indexes);
    free(comps.stokesV_list_comp_inds);
    free(comps.stokesV_list_ref_freqs);
    free(comps.stokesV_list_ref_flux);
    free(comps.stokesQ_num_list_values);
    free(comps.stokesQ_list_start_indexes);
    free(comps.stokesQ_list_comp_inds);
    free(comps.stokesQ_list_ref_freqs);
    free(comps.stokesQ_list_ref_flux);
    free(comps.stokesU_num_list_values);
    free(comps.stokesU_list_start_indexes);
    free(comps.stokesU_list_comp_inds);
    free(comps.stokesU_list_ref_freqs);
    free(comps.stokesU_list_ref_flux);
    free(comps.linpol_p_num_list_values);
    free(comps.linpol_p_list_start_indexes);
    free(comps.linpol_p_list_comp_inds);
    free(comps.linpol_p_list_ref_freqs);
    free(comps.linpol_p_list_ref_flux);
  }
}

/*
Based on what is in the sky model, free stuff
*/
void free_sky_model(source_catalogue_t *cropped_sky_models) {

  for (int cats_ind = 0; cats_ind < cropped_sky_models->num_sources; cats_ind++) {

    int n_points = cropped_sky_models->sources[cats_ind].n_points;
    int n_gauss = cropped_sky_models->sources[cats_ind].n_gauss;
    int n_shapes = cropped_sky_models->sources[cats_ind].n_shapes;

    if (n_points > 0) {
      free_components(cropped_sky_models->sources[cats_ind].point_components,
                      n_points);
    }
    if (n_gauss > 0) {
      free_components(cropped_sky_models->sources[cats_ind].gauss_components,
                      n_gauss);
    }
    if (n_shapes > 0) {
      free_components(cropped_sky_models->sources[cats_ind].shape_components,
                      n_shapes);
    }

  }
  free(cropped_sky_models->sources);
  free(cropped_sky_models);
}

woden_settings_t * make_woden_settings(double ra0, double dec0) {

    woden_settings_t *woden_settings = malloc(sizeof(woden_settings_t));

    woden_settings->ra0 = ra0;
    woden_settings->dec0 = dec0;
    woden_settings->sdec0 = sin(dec0);
    woden_settings->cdec0 = cos(dec0);
    woden_settings->num_baselines = NUM_BASELINES;
    woden_settings->num_freqs = NUM_FREQS;
    woden_settings->num_time_steps = NUM_TIME_STEPS;
    woden_settings->num_ants = NUM_ANTS;
    woden_settings->num_autos = 0;
    woden_settings->num_cross = NUM_CROSS;
    woden_settings->num_visis = woden_settings->num_cross;
    woden_settings->coarse_band_width = 1.28e+6;
    //Make the fine channel width insanely small so beam changes little
    //with frequency - that way we can test for just one gain value per time
    woden_settings->frequency_resolution = 1e-6;
    woden_settings->latitude = -0.46606083776035967;
    double latitudes[] = {-0.46606083776035967, -0.46606083776035967};
    woden_settings->latitudes = latitudes;

    woden_settings->lsts = lsts;

    woden_settings->do_autos = 0;
    woden_settings->use_dipamps = 0;

    return woden_settings;
}

/*
Due to the way I've setup the array layout and phase centre,
the uvw should be a simple reduced product between phase centre and lsts.
Create some expected u,v,w arrays and compare to what we've found
*/
void test_uvw(visibility_set_t *visibility_set,  double *lsts,
              double ra0, double dec0,
              woden_settings_t *woden_settings) {

  // int num_visis = NUM_VISI;
  double *expec_u = malloc(NUM_CROSS * sizeof(double));
  double *expec_v = malloc(NUM_CROSS * sizeof(double));
  double *expec_w = malloc(NUM_CROSS * sizeof(double));

  int index = 0;
  double xy_length;
  double ha0;
  for (int time = 0; time < NUM_TIME_STEPS; time++) {
    ha0 = lsts[time] - ra0;
    for (int freq = 0; freq < NUM_FREQS; freq++) {
      for (int baseline = 0; baseline < NUM_BASELINES; baseline++) {
        xy_length = (baseline + 1) * 100;
        expec_u[index] = (cos(ha0) + sin(ha0))*xy_length;
        expec_v[index] = xy_length * sin(-0.46606083776035967)*(-cos(ha0) + sin(ha0));
        expec_w[index] = xy_length * cos(-0.46606083776035967)*(cos(ha0) - sin(ha0));

        index ++;
      }
    }
  }

  #ifdef DOUBLE_PRECISION
    double TOL = 1e-12;
  #else
    double TOL = 1e-5;
  #endif

  for (int visi = 0; visi < woden_settings->num_cross; visi++) {
    // printf("%.1f %.1f %.1f %.1f %.1f %.1f\n",
    //         visibility_set->us_metres[visi], expec_u[visi],
    //         visibility_set->vs_metres[visi], expec_v[visi],
    //         visibility_set->ws_metres[visi], expec_w[visi] );
    TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_u[visi],
                              visibility_set->us_metres[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_v[visi],
                              visibility_set->vs_metres[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_w[visi],
                              visibility_set->ws_metres[visi]);

  }

  int num_autos = woden_settings->num_autos;

  //This should only be hit if we have asked for autos
  //All of these should just be zero, as there is no baseline length for an
  //auto correlation
  for (int visi = 0; visi < num_autos; visi++) {
    TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0,
                              visibility_set->us_metres[num_autos + visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0,
                              visibility_set->vs_metres[num_autos + visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0,
                              visibility_set->ws_metres[num_autos + visi]);

  }


  free(expec_u);
  free(expec_v);
  free(expec_w);
}

/*
Pump many many many settings into the function we are trying to test
Checkthe u,v,w are correct and return the visibility_set for further testing
*/
visibility_set_t * test_calculate_visibilities(source_catalogue_t *cropped_sky_models,
                                 beam_settings_t *beam_settings,
                                 woden_settings_t *woden_settings,
                                 double ra0, double dec0,
                                 int beamtype) {

  double base_band_freq = BASE_BAND_FREQ;
  // user_precision_t base_band_freq = 120e+6;

  array_layout_t *array_layout = malloc(sizeof(array_layout_t));

  array_layout->X_diff_metres = malloc(NUM_TIME_STEPS*NUM_BASELINES*sizeof(double));
  array_layout->Y_diff_metres = malloc(NUM_TIME_STEPS*NUM_BASELINES*sizeof(double));
  array_layout->Z_diff_metres = malloc(NUM_TIME_STEPS*NUM_BASELINES*sizeof(double));

  for (int time_ind = 0; time_ind < NUM_TIME_STEPS; time_ind++) {
      for (int baseline = 0; baseline < NUM_BASELINES; baseline++) {
        int time_off = time_ind*NUM_BASELINES;
        array_layout->X_diff_metres[time_off + baseline] = (baseline + 1) * 100;
        array_layout->Y_diff_metres[time_off + baseline] = (baseline + 1) * 100;
        array_layout->Z_diff_metres[time_off + baseline] = 0.0;
      }
  }



  user_precision_t *sbf = NULL;
  if (cropped_sky_models->num_shapelets > 0) {
    sbf = malloc( sbf_N * sbf_L * sizeof(user_precision_t) );
    sbf = create_sbf(sbf);
  }

  // printf("SIZE OF THING %d\n",woden_settings->num_visis );

  visibility_set_t *visibility_set = setup_visibility_set(woden_settings->num_visis);

  fill_timefreq_visibility_set(visibility_set, woden_settings,
                               base_band_freq, lsts);

  if (beam_settings->beamtype == FEE_BEAM) {
    char* mwa_fee_hdf5 = getenv("MWA_FEE_HDF5");

    if (mwa_fee_hdf5) {
      printf("MWA_FEE_HDF5: %s\n", mwa_fee_hdf5 );

      int status =  new_fee_beam(mwa_fee_hdf5, &beam_settings->fee_beam);
      if (status != 0) {
        handle_hyperbeam_error(__FILE__, __LINE__, "new_fee_beam");
        // printf("Something went wrong calling new_fee_beam\n");
      } 
    } else{
      printf("MWA_FEE_HDF5 not found - not running test_hyperbeam test\n");
    }
  }
  else if (beam_settings->beamtype == FEE_BEAM_INTERP) {
    char* mwa_fee_hdf5 = getenv("MWA_FEE_HDF5_INTERP");

    if (mwa_fee_hdf5) {
      printf("MWA_FEE_HDF5_INTERP: %s\n", mwa_fee_hdf5 );

      int status =  new_fee_beam(mwa_fee_hdf5, &beam_settings->fee_beam);
      if (status != 0) {
        handle_hyperbeam_error(__FILE__, __LINE__, "new_fee_beam");
        // printf("Something went wrong calling new_fee_beam\n");
      }
    } else {
        printf("MWA_FEE_HDF5_INTERP not found - not running test_hyperbeam test\n");
    }
  }

  printf("Calling GPU\n");
  calculate_visibilities(array_layout, cropped_sky_models, beam_settings,
                         woden_settings, visibility_set, sbf);

  // for (int visi = 0; visi < woden_settings->num_visis; visi++) {
  //   printf("\t\tVisi %d: %f %f %f %f %f %f %f %f\n", visi,
  //               visibility_set->sum_visi_XX_real[visi],
  //               visibility_set->sum_visi_XX_imag[visi],
  //               visibility_set->sum_visi_XY_real[visi],
  //               visibility_set->sum_visi_XY_imag[visi],
  //               visibility_set->sum_visi_YX_real[visi],
  //               visibility_set->sum_visi_YX_imag[visi],
  //               visibility_set->sum_visi_YY_real[visi],
  //               visibility_set->sum_visi_YY_imag[visi]);
  // }

  printf("GPU has finished\n");

  //Be free my pretties!
  if (cropped_sky_models->num_shapelets > 0) {
    free(sbf);
  }

  free_sky_model(cropped_sky_models);
  free(array_layout->X_diff_metres);
  free(array_layout->Y_diff_metres);
  free(array_layout->Z_diff_metres);
  // free(array_layout);

  test_uvw(visibility_set, lsts, ra0, dec0, woden_settings);

  return visibility_set;

}


//Due to the way I've set up the sky model, the I,Q,U,V fluxes should all be
//some multiple of the number of components. Everything is also stuck at phase
//centre, so visis should be the same for all baselines. Give given the complex
//beam values, we can predict the visibilities
void predict_inst_stokes(int num_comps, double _Complex g1x, double _Complex D1x,
                                        double _Complex D1y, double _Complex g1y,
                                        double _Complex g2x, double _Complex D2x,
                                        double _Complex D2y, double _Complex g2y,
                                        double _Complex * xx, double _Complex * xy,
                                        double _Complex * yx, double _Complex * yy) {
  
  double _Complex stokesI = num_comps*STOKESI + I*0;
  double _Complex stokesQ = num_comps*STOKESI + I*0;
  // double _Complex stokesU = num_comps*STOKESI + I*0;
  double _Complex stokesU = 0 + I*0;
  double _Complex stokesV = num_comps*STOKESI + I*0;

  double _Complex g2x_conj = conj(g2x);
  double _Complex D2x_conj = conj(D2x);
  double _Complex D2y_conj = conj(D2y);
  double _Complex g2y_conj = conj(g2y);

  * xx = (g1x*g2x_conj + D1x*D2x_conj)*stokesI;
  * xx += (g1x*g2x_conj - D1x*D2x_conj)*stokesQ;
  * xx += (g1x*D2x_conj + D1x*g2x_conj)*stokesU;
  * xx += ((0.0 + I*1.0)*stokesV)*(g1x*D2x_conj - D1x*g2x_conj);

  * xy = (g1x*D2y_conj + D1x*g2y_conj)*stokesI;
  * xy += (g1x*D2y_conj - D1x*g2y_conj)*stokesQ;
  * xy += (g1x*g2y_conj + D1x*D2y_conj)*stokesU;
  * xy += ((0.0 + I*1.0)*stokesV)* (g1x*g2y_conj - D1x*D2y_conj);

  * yx = (D1y*g2x_conj + g1y*D2x_conj)*stokesI;
  * yx += (D1y*g2x_conj - g1y*D2x_conj)*stokesQ;
  * yx += (D1y*D2x_conj + g1y*g2x_conj)*stokesU;
  * yx += ((0.0 + I*1.0)*stokesV)* (D1y*D2x_conj - g1y*g2x_conj);

  * yy = (D1y*D2y_conj + g1y*g2y_conj)*stokesI;
  * yy += (D1y*D2y_conj - g1y*g2y_conj)*stokesQ;
  * yy += (D1y*g2y_conj + g1y*D2y_conj)*stokesU;
  * yy += ((0.0 + I*1.0)*stokesV)* (D1y*g2y_conj - g1y*D2y_conj);

}

/*
This test works for beams that change with time, but have no cross-pol info
Can just test whether the XX/YY real values match expected
*/
void test_comp_phase_centre_twogains(visibility_set_t *visibility_set,
                                     int num_comps,
                                     double _Complex gain1x, double _Complex gain1y,
                                     double _Complex gain2x, double _Complex gain2y,
                                     woden_settings_t *woden_settings) {

  int num_cross = woden_settings->num_cross;

  double _Complex *expec_xx = malloc(num_cross*sizeof(double _Complex));
  double _Complex *expec_xy = malloc(num_cross*sizeof(double _Complex));
  double _Complex *expec_yx = malloc(num_cross*sizeof(double _Complex));
  double _Complex *expec_yy = malloc(num_cross*sizeof(double _Complex));

  double _Complex leak = 0.0 + I*0.0;
  double _Complex xx, xy, yx, yy;

  double _Complex gainsx[] = {gain1x, gain2x};
  double _Complex gainsy[] = {gain1y, gain2y};
  
  for (int time = 0; time < NUM_TIME_STEPS; time++) {
    predict_inst_stokes(num_comps, gainsx[time], leak, leak, gainsy[time],
                                   gainsx[time], leak, leak, gainsy[time],
                        &xx, &xy, &yx, &yy);

    for (int visi = 0; visi < num_cross/2; visi++) {
      expec_xx[time*(num_cross/2) + visi] = xx;
      expec_xy[time*(num_cross/2) + visi] = xy;
      expec_yx[time*(num_cross/2) + visi] = yx;
      expec_yy[time*(num_cross/2) + visi] = yy;
    }

  }

  #ifdef DOUBLE_PRECISION
    double TOL = 1e-8;
  #else
    double TOL = 1e-5;
  #endif

  //We're testing all COMPONENT types with this function, where the values of
  //the real vary slightly for baseline (I've set the major/minor to be small
  //enough to be close to 1.0 but not quite)
  for (int visi = 0; visi < num_cross; visi++) {

    // printf("%.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n",
    //        creal(expec_xx[visi]), cimag(expec_xx[visi]),
    //        creal(expec_xy[visi]), cimag(expec_xy[visi]),
    //        creal(expec_yx[visi]), cimag(expec_yx[visi]),
    //        creal(expec_yy[visi]), cimag(expec_yy[visi]));


    // printf("%.1f %.1f %.1f %.1f %.1f %.1f %.1f %.1f\n",
    //         visibility_set->sum_visi_XX_real[visi],
    //         visibility_set->sum_visi_XX_imag[visi],
    //         visibility_set->sum_visi_XY_real[visi],
    //         visibility_set->sum_visi_XY_imag[visi],
    //         visibility_set->sum_visi_YX_real[visi],
    //         visibility_set->sum_visi_YX_imag[visi],
    //         visibility_set->sum_visi_YY_imag[visi],
    //         visibility_set->sum_visi_YY_real[visi] );

    // printf("Expec GAIN %.16f %.16f\n", expec_gainsxx[visi],
    //                             visibility_set->sum_visi_XX_real[visi] );



    TEST_ASSERT_DOUBLE_WITHIN(TOL, creal(expec_xx[visi]),
                              visibility_set->sum_visi_XX_real[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, cimag(expec_xx[visi]),
                              visibility_set->sum_visi_XX_imag[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, creal(expec_xy[visi]),
                              visibility_set->sum_visi_XY_real[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, cimag(expec_xy[visi]),
                              visibility_set->sum_visi_XY_imag[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, creal(expec_yx[visi]),
                              visibility_set->sum_visi_YX_real[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, cimag(expec_yx[visi]),
                              visibility_set->sum_visi_YX_imag[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, creal(expec_yy[visi]),
                              visibility_set->sum_visi_YY_real[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, cimag(expec_yy[visi]),
                              visibility_set->sum_visi_YY_imag[visi]);

  }

  for (int visi = 0; visi < woden_settings->num_autos; visi++) {
    // printf("%.1f %.1f %.1f %.1f %.1f %.1f %.1f %.1f\n",
            // visibility_set->sum_visi_XX_real[num_cross + visi],
            // visibility_set->sum_visi_XX_imag[num_cross + visi],
            // visibility_set->sum_visi_XY_real[num_cross + visi],
            // visibility_set->sum_visi_XY_imag[num_cross + visi],
            // visibility_set->sum_visi_YX_real[num_cross + visi],
            // visibility_set->sum_visi_YX_imag[num_cross + visi],
            // visibility_set->sum_visi_YY_imag[num_cross + visi],
            // visibility_set->sum_visi_YY_real[num_cross + visi] );

    // printf("Expec GAIN %.16f %.16f\n", expec_gainsxx[visi],
    //                             visibility_set->sum_visi_XX_real[visi] );

    TEST_ASSERT_DOUBLE_WITHIN(TOL, creal(expec_xx[visi]),
                              visibility_set->sum_visi_XX_real[num_cross + visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, cimag(expec_xx[visi]),
                              visibility_set->sum_visi_XX_imag[num_cross + visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, creal(expec_xy[visi]),
                              visibility_set->sum_visi_XY_real[num_cross + visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, cimag(expec_xy[visi]),
                              visibility_set->sum_visi_XY_imag[num_cross + visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, creal(expec_yx[visi]),
                              visibility_set->sum_visi_YX_real[num_cross + visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, cimag(expec_yx[visi]),
                              visibility_set->sum_visi_YX_imag[num_cross + visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, creal(expec_yy[visi]),
                              visibility_set->sum_visi_YY_real[num_cross + visi]);
    TEST_ASSERT_DOUBLE_WITHIN(TOL, cimag(expec_yy[visi]),
                              visibility_set->sum_visi_YY_imag[num_cross + visi]);

  }

  free(expec_xx);
  free(expec_xy);
  free(expec_yx);
  free(expec_yy);
  
}


void test_comp_phase_centre_allgains(visibility_set_t *visibility_set,
      int num_comps,
      double _Complex gain1x, double _Complex leak1x,
      double _Complex leak1y, double _Complex gain1y,
      double _Complex gain2x, double _Complex leak2x,
      double _Complex leak2y, double _Complex gain2y,
      woden_settings_t *woden_settings, double tol) {

  int num_cross = woden_settings->num_cross;

  double _Complex xx, xy, yx, yy;

  for (int visi = 0; visi < num_cross; visi++) {
    if (visi < num_cross / 2) {
      predict_inst_stokes(num_comps, gain1x, leak1x, leak1y, gain1y,
                                     gain1x, leak1x, leak1y, gain1y,
                                     &xx, &xy, &yx, &yy);
    } else {
      predict_inst_stokes(num_comps, gain2x, leak2x, leak2y, gain2y,
                                     gain2x, leak2x, leak2y, gain2y,
                                     &xx, &xy, &yx, &yy);
    }
    // printf("%.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n",
    //        creal(xx), cimag(xx), creal(xy), cimag(xy),
    //        creal(yx), cimag(yx), creal(yy), cimag(yy));

    TEST_ASSERT_DOUBLE_WITHIN(tol, creal(xx),
                                      visibility_set->sum_visi_XX_real[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(tol, cimag(xx),
                                      visibility_set->sum_visi_XX_imag[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(tol, creal(xy),
                                      visibility_set->sum_visi_XY_real[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(tol, cimag(xy),
                                      visibility_set->sum_visi_XY_imag[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(tol, creal(yx),
                                      visibility_set->sum_visi_YX_real[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(tol, cimag(yx),
                                      visibility_set->sum_visi_YX_imag[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(tol, creal(yy),
                                      visibility_set->sum_visi_YY_real[visi]);
    TEST_ASSERT_DOUBLE_WITHIN(tol, cimag(yy),
                                      visibility_set->sum_visi_YY_imag[visi]);
  }

  for (int visi = 0; visi < woden_settings->num_autos; visi++) {
    if (visi < num_cross / 2) {
      predict_inst_stokes(num_comps, gain1x, leak1x, leak1y, gain1y,
                                     gain1x, leak1x, leak1y, gain1y,
                                     &xx, &xy, &yx, &yy);
    } else {
      predict_inst_stokes(num_comps, gain2x, leak2x, leak2y, gain2y,
                                     gain2x, leak2x, leak2y, gain2y,
                                     &xx, &xy, &yx, &yy);
    }
    // printf("%.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n",
    //        creal(xx), cimag(xx), creal(xy), cimag(xy),
    //        creal(yx), cimag(yx), creal(yy), cimag(yy));
    
    TEST_ASSERT_DOUBLE_WITHIN(tol, creal(xx),
                                      visibility_set->sum_visi_XX_real[num_cross + visi]);
    TEST_ASSERT_DOUBLE_WITHIN(tol, cimag(xx),
                                      visibility_set->sum_visi_XX_imag[num_cross + visi]);
    TEST_ASSERT_DOUBLE_WITHIN(tol, creal(xy),
                                      visibility_set->sum_visi_XY_real[num_cross + visi]);
    TEST_ASSERT_DOUBLE_WITHIN(tol, cimag(xy),
                                      visibility_set->sum_visi_XY_imag[num_cross + visi]);
    TEST_ASSERT_DOUBLE_WITHIN(tol, creal(yx),
                                      visibility_set->sum_visi_YX_real[num_cross + visi]);
    TEST_ASSERT_DOUBLE_WITHIN(tol, cimag(yx),
                                      visibility_set->sum_visi_YX_imag[num_cross + visi]);
    TEST_ASSERT_DOUBLE_WITHIN(tol, creal(yy),
                                      visibility_set->sum_visi_YY_real[num_cross + visi]);
    TEST_ASSERT_DOUBLE_WITHIN(tol, cimag(yy),
                                      visibility_set->sum_visi_YY_imag[num_cross + visi]);
  }
}


void test_comp_phase_centre_allgains_multiants(visibility_set_t *visibility_set,
                                int num_comps, 
                                double _Complex gain1x, double _Complex leak1x,
                                double _Complex leak1y, double _Complex gain1y,
                                double _Complex gain2x, double _Complex leak2x,
                                double _Complex leak2y, double _Complex gain2y,
                                double *antx_mult, double *anty_mult, int num_ants,
                                woden_settings_t *woden_settings,
                                double tol) {

  int num_cross = woden_settings->num_cross;

  int num_baselines = ((num_ants - 1)*num_ants) / 2;

  int *ant1_to_baseline_map = NULL;
  int *ant2_to_baseline_map = NULL;

  ant1_to_baseline_map = (int *)malloc(num_baselines*sizeof(int));
  ant2_to_baseline_map = (int *)malloc(num_baselines*sizeof(int));

  //These functions only do cross correlations, so create all combos of antennas
  //that make up all the crosses
  int cross_index = 0;
  for (int ant1 = 0; ant1 < num_ants-1; ant1++)
  {
    for (int ant2 = ant1 + 1; ant2 < num_ants; ant2++)
    {
      ant1_to_baseline_map[cross_index] = ant1;
      ant2_to_baseline_map[cross_index] = ant2;

      cross_index += 1;
    }
  }

  //In all of the below, we've set all 16 dipole amplitudes to a constant value
  //Can multiply expected values by combinations of constants to predict effect
  //on stored values. This means leakages are multiplied by same factor as gains
  //so can say g1x_mult = D1x_mult etc
  float g1x_mult, g1y_mult, g2x_mult, g2y_mult;
  // float xx_gains, xy_gains, yx_gains, yy_gains;

  double _Complex xx, xy, yx, yy;

  // for (int time = 0; time < woden_settings->num_time_steps; time ++) {

  int out_ind = 0;

  for (int freq = 0; freq < woden_settings->num_freqs; freq ++) {
    for (int baseline = 0; baseline < num_baselines; baseline ++) {

      g1x_mult = antx_mult[ant1_to_baseline_map[baseline]];
      g1y_mult = anty_mult[ant1_to_baseline_map[baseline]];
      g2x_mult = antx_mult[ant2_to_baseline_map[baseline]];
      g2y_mult = anty_mult[ant2_to_baseline_map[baseline]];

      for (int time = 0; time < 2; time++)
      {
        if (time == 0){
          predict_inst_stokes(num_comps, g1x_mult*gain1x, g1x_mult*leak1x, g1y_mult*leak1y, g1y_mult*gain1y,
                                     g2x_mult*gain1x, g2x_mult*leak1x, g2y_mult*leak1y, g2y_mult*gain1y,
                                     &xx, &xy, &yx, &yy);
        } else {
          predict_inst_stokes(num_comps, g1x_mult*gain2x, g1x_mult*leak2x, g1y_mult*leak2y, g1y_mult*gain2y,
                                     g2x_mult*gain2x, g2x_mult*leak2x, g2y_mult*leak2y, g2y_mult*gain2y,
                                     &xx, &xy, &yx, &yy);
        }

        out_ind = time*woden_settings->num_freqs*num_baselines + freq*num_baselines + baseline;

        TEST_ASSERT_DOUBLE_WITHIN(tol, creal(xx), visibility_set->sum_visi_XX_real[out_ind]);
        TEST_ASSERT_DOUBLE_WITHIN(tol, cimag(xx), visibility_set->sum_visi_XX_imag[out_ind]);
        TEST_ASSERT_DOUBLE_WITHIN(tol, creal(xy), visibility_set->sum_visi_XY_real[out_ind]);
        TEST_ASSERT_DOUBLE_WITHIN(tol, cimag(xy), visibility_set->sum_visi_XY_imag[out_ind]);
        TEST_ASSERT_DOUBLE_WITHIN(tol, creal(yx), visibility_set->sum_visi_YX_real[out_ind]);
        TEST_ASSERT_DOUBLE_WITHIN(tol, cimag(yx), visibility_set->sum_visi_YX_imag[out_ind]);
        TEST_ASSERT_DOUBLE_WITHIN(tol, creal(yy), visibility_set->sum_visi_YY_real[out_ind]);
        TEST_ASSERT_DOUBLE_WITHIN(tol, cimag(yy), visibility_set->sum_visi_YY_imag[out_ind]);
      }
    }
  }

  if (woden_settings->do_autos == 1){


    for (int freq = 0; freq < woden_settings->num_freqs; freq ++) {
      for (int ant = 0; ant < num_ants; ant ++) {

        g1x_mult = antx_mult[ant];
        g1y_mult = anty_mult[ant];

        for (int time = 0; time < 2; time++) {
          if (time == 0){
            predict_inst_stokes(num_comps, g1x_mult*gain1x, g1x_mult*leak1x, g1y_mult*leak1y, g1y_mult*gain1y,
                                           g1x_mult*gain1x, g1x_mult*leak1x, g1y_mult*leak1y, g1y_mult*gain1y,
                                      &xx, &xy, &yx, &yy);
          } else {
            predict_inst_stokes(num_comps, g1x_mult*gain2x, g1x_mult*leak2x, g1y_mult*leak2y, g1y_mult*gain2y,
                                           g1x_mult*gain2x, g1x_mult*leak2x, g1y_mult*leak2y, g1y_mult*gain2y,
                                      &xx, &xy, &yx, &yy);
          }

  
          out_ind = num_cross + time*woden_settings->num_freqs*num_ants + freq*num_ants + ant;

          TEST_ASSERT_DOUBLE_WITHIN(tol, creal(xx), visibility_set->sum_visi_XX_real[out_ind]);
          TEST_ASSERT_DOUBLE_WITHIN(tol, cimag(xx), visibility_set->sum_visi_XX_imag[out_ind]);
          TEST_ASSERT_DOUBLE_WITHIN(tol, creal(xy), visibility_set->sum_visi_XY_real[out_ind]);
          TEST_ASSERT_DOUBLE_WITHIN(tol, cimag(xy), visibility_set->sum_visi_XY_imag[out_ind]);
          TEST_ASSERT_DOUBLE_WITHIN(tol, creal(yx), visibility_set->sum_visi_YX_real[out_ind]);
          TEST_ASSERT_DOUBLE_WITHIN(tol, cimag(yx), visibility_set->sum_visi_YX_imag[out_ind]);
          TEST_ASSERT_DOUBLE_WITHIN(tol, creal(yy), visibility_set->sum_visi_YY_real[out_ind]);
          TEST_ASSERT_DOUBLE_WITHIN(tol, cimag(yy), visibility_set->sum_visi_YY_imag[out_ind]);
        }
      }
    }
  }


  free(ant1_to_baseline_map);
  free(ant2_to_baseline_map);
}