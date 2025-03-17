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

// #include "calculate_visibilities_mwafeebeam_common.h"

// void setUp (void) {} /* Is run before every test, put unit init calls here. */
// void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

#include <math.h>
#include <unity.h>
#include <stdlib.h>
#include <complex.h>

#include "constants.h"
#include "woden_struct_defs.h"
#include "shapelet_basis.h"
#include "visibility_set.h"
#include "shapelet_basis.h"
#include "hyperbeam_error.h"
#include "calculate_visibilities_common.h"

#define NUM_ANTS 128
#define NUM_BASELINES ((NUM_ANTS - 1) * NUM_ANTS / 2)
#define NUM_FREQS 10
#define NUM_TIME_STEPS 2
#define NUM_CROSS NUM_BASELINES*NUM_FREQS*NUM_TIME_STEPS
#define NUM_VISI NUM_ANTS*NUM_FREQS*NUM_TIME_STEPS
#define RA0 0.0
#define DEC0 -0.46606083776035967
#define BASE_BAND_FREQ 120000000.0
#define STOKESI 0.3333333333333333
// #define STOKESI 1.0
#define POL_FRAC 1.0

double lsts[] = {0.0, M_PI / 4};

user_precision_t azs[] = {0.0, 4.528359553989764};
user_precision_t zas[] = {0.0, 0.6978088917603547};

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
  comps->power_SIs = malloc(sizeof(user_precision_t));
  comps->power_comp_inds = malloc(sizeof(int));

  comps->power_ref_freqs[0] = REF_FREQ;
  comps->power_ref_stokesI[0] = STOKESI;
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

  // printf("Calling calculate_visibilities\n");
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

  // printf("calculate_visibilities has finished\n");

  //Be free my pretties!
  if (cropped_sky_models->num_shapelets > 0) {
    free(sbf);
  }

  free_sky_model(cropped_sky_models);
  free(array_layout->X_diff_metres);
  free(array_layout->Y_diff_metres);
  free(array_layout->Z_diff_metres);
  // free(array_layout);

  return visibility_set;

}

void test_calculate_visibilities_MWAFEEBeam(int n_points, int n_gauss, int n_shapes,
                                           int num_sources, int do_gpu) {

  source_catalogue_t *cropped_sky_models = make_cropped_sky_models(RA0, DEC0,
                                                    n_points, n_gauss, n_shapes,
                                                    num_sources);

  woden_settings_t *woden_settings = malloc(sizeof(woden_settings_t));

  float dec0 = DEC0;
  woden_settings->ra0 = RA0;
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
  woden_settings->latitude = DEC0;
  double latitudes[] = {DEC0, DEC0};
  woden_settings->latitudes = latitudes;

  woden_settings->lsts = lsts;

  woden_settings->do_autos = 0;
  woden_settings->use_dipamps = 0;
  woden_settings->beamtype = FEE_BEAM;
  woden_settings->do_gpu = do_gpu;

  int delays[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  // int delays[16] = {2,4,6,8,2,4,6,8,2,4,6,8,2,4,6,8};
  woden_settings->FEE_ideal_delays = delays;
  
  beam_settings_t *beam_settings = malloc(sizeof(beam_settings_t));
  // beam_settings->beamtype = FEE_BEAM;
  beam_settings->beamtype = ANALY_DIPOLE;


  woden_settings->do_autos = 1;
  woden_settings->num_autos = NUM_CROSS;
  woden_settings->num_visis = woden_settings->num_cross + woden_settings->num_autos;

  cropped_sky_models = make_cropped_sky_models(RA0, DEC0,
                                                    n_points, n_gauss, n_shapes,
                                                    num_sources);

  printf("We have this many visis %d %d %d\n",woden_settings->num_visis,woden_settings->num_autos,woden_settings->num_cross );
  visibility_set_t *visibility_set = test_calculate_visibilities(cropped_sky_models,
                                          beam_settings, woden_settings, RA0, DEC0,
                                          beam_settings->beamtype);

  // test_comp_phase_centre_allgains(visibility_set, num_comps,
  //                                 gain1x, leak1x, leak1y, gain1y,
  //                                 gain2x, leak2x, leak2y, gain2y,
  //                                 woden_settings, TOL);

  // free_fee_beam(beam_settings->fee_beam);
  free(beam_settings);
  free_visi_set_inputs(visibility_set);
  free_visi_set_outputs(visibility_set);
  free(woden_settings);

}

// Run the test with unity
int main(void)
{
    // UNITY_BEGIN();

    char* mwa_fee_hdf5 = getenv("MWA_FEE_HDF5");

    if (mwa_fee_hdf5) {
      printf("MWA_FEE_HDF5: %s\n", mwa_fee_hdf5 );

        int do_gpu = 0;
        int n_points = 10;
        int n_gauss = 10;
        int n_shapes = 10;
        int num_sources = 10;
        test_calculate_visibilities_MWAFEEBeam(n_points, n_gauss, n_shapes, num_sources,
                                                do_gpu);

    }
    else {
      printf("MWA_FEE_HDF5 not found - not running test_calculate_visibilities_MWAFEEBeam tests");
    }

    // return UNITY_END();
}
