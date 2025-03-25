/*

*/

#include "calculate_visibilities_common.h"
#include "calculate_visibilities_common_common.h"
#include "calculate_visibilities_everybeam_common.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

// void change_radecs(components_t *components, int n_comps, double ra0, double dec0) {
//   for (int i = 0; i < n_comps; i++) {
//     components->ras[i] = ra0 + 0.1*DD2R;
//     components->decs[i] = dec0 + 0.1*DD2R;
//   }
// }





void populate_these_components(components_t *comps, int n_comps,
     double ra0, double dec0, int num_freqs, int num_time_steps) {

  comps->ras = malloc(n_comps*sizeof(double));
  comps->decs = malloc(n_comps*sizeof(double));

  comps->num_primarybeam_values = n_comps*num_freqs*num_time_steps;
  comps->azs = malloc(n_comps*num_time_steps*sizeof(user_precision_t));
  comps->zas = malloc(n_comps*num_time_steps*sizeof(user_precision_t));
  comps->para_angles = malloc(n_comps*num_time_steps*sizeof(user_precision_t));

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

  
  //Make sure each direction is different just incase things are cached
  for (int comp = 0; comp < n_comps; comp++) {
    comps->ras[comp] = ra0 + comp*0.01*DD2R;
    comps->decs[comp] = dec0 + comp*0.01*DD2R;

    // comps->majors[comp] = 1e-10;
    // comps->minors[comp] = 1e-10;
    // comps->pas[comp] = 0;

    // comps->shape_coeffs[comp] = 1.0;
    // comps->n1s[comp] = 0.0;
    // comps->n2s[comp] = 0.0;
    // comps->param_indexes[comp] = 0.0;

  }

}


source_catalogue_t * make_these_cropped_sky_models(double ra0, double dec0,
                      int n_points, int num_sources, int num_freqs, int num_time_steps) {

source_catalogue_t *cropped_sky_models = malloc(sizeof(cropped_sky_models));
cropped_sky_models->num_sources = num_sources;
cropped_sky_models->num_shapelets = 0.0;
cropped_sky_models->sources = malloc(num_sources*sizeof(source_t));

  for (int cats_ind = 0; cats_ind < num_sources; cats_ind++) {
    cropped_sky_models->sources[cats_ind].n_points = n_points;
    cropped_sky_models->sources[cats_ind].n_gauss = 0;
    cropped_sky_models->sources[cats_ind].n_shapes = 0;
    cropped_sky_models->sources[cats_ind].n_shape_coeffs = 0;
    cropped_sky_models->sources[cats_ind].n_point_powers = 0;
    cropped_sky_models->sources[cats_ind].n_point_curves = 0;
    cropped_sky_models->sources[cats_ind].n_point_lists = 0;
    cropped_sky_models->sources[cats_ind].n_gauss_powers = 0;
    cropped_sky_models->sources[cats_ind].n_gauss_curves = 0;
    cropped_sky_models->sources[cats_ind].n_gauss_lists = 0;
    cropped_sky_models->sources[cats_ind].n_shape_powers = 0;
    cropped_sky_models->sources[cats_ind].n_shape_curves = 0;
    cropped_sky_models->sources[cats_ind].n_shape_lists = 0;

    populate_these_components(&cropped_sky_models->sources[cats_ind].point_components,
      n_points, ra0, dec0, num_freqs, num_time_steps);

  }

  return cropped_sky_models;
}


void do_test(){
  // profile_lofar_everybeam(0, EB_LOFAR,
  //                  "../../../../test_installation/everybeam/LOFAR_HBA_MOCK.ms");

  int n_points = 50;
  int n_gauss = 0;
  int n_shapes = 0;
  int num_sources = 3;

  int num_time_steps = 10;
  int num_ants = 10;
  int num_baselines = num_ants * (num_ants - 1) / 2;
  int num_freqs = 20;

  char* beam_ms_path = "../../../../test_installation/everybeam/LOFAR_HBA_MOCK.ms";

  double LOFAR_RA0 = 2.15374123;
  double LOFAR_DEC0 = 0.84155210;

  source_catalogue_t *cropped_sky_models = make_these_cropped_sky_models(LOFAR_RA0, LOFAR_DEC0,
                                                                   n_points, num_sources,
                                                                   num_freqs, num_time_steps);

  // for (int s_ind = 0; s_ind < num_sources; s_ind++)
  // {
  //   change_radecs(&cropped_sky_models->sources[s_ind].point_components, n_points, LOFAR_RA0, LOFAR_DEC0);
  //   change_radecs(&cropped_sky_models->sources[s_ind].gauss_components, n_gauss, LOFAR_RA0, LOFAR_DEC0);
  //   change_radecs(&cropped_sky_models->sources[s_ind].shape_components, n_shapes, LOFAR_RA0, LOFAR_DEC0);
  // }
  

  

  // woden_settings_t *woden_settings = make_woden_settings(LOFAR_RA0, LOFAR_DEC0);

  woden_settings_t *woden_settings = malloc(sizeof(woden_settings_t));

  int beamtype = EB_LOFAR;

  woden_settings->beamtype = beamtype;
  woden_settings->do_gpu = 0;
  woden_settings->use_dipamps = 0;
  woden_settings->do_autos = 0;
  woden_settings->beam_ms_path = beam_ms_path;
  // woden_settings->hdf5_beam_path = getenv("MWA_FEE_HDF5");
  woden_settings->off_cardinal_dipoles = 0;
  woden_settings->eb_beam_ra0 = LOFAR_RA0 - DD2R;
  woden_settings->eb_beam_dec0 = LOFAR_DEC0;
  woden_settings->normalise_primary_beam = 1;

  woden_settings->ra0 = LOFAR_RA0;
  woden_settings->dec0 = LOFAR_DEC0;
  woden_settings->sdec0 = sin(LOFAR_DEC0);
  woden_settings->cdec0 = cos(LOFAR_DEC0);
  woden_settings->num_baselines = num_baselines;
  woden_settings->num_freqs = num_freqs;
  woden_settings->num_time_steps = num_time_steps;
  woden_settings->num_ants = num_ants;
  woden_settings->num_autos = 0;
  woden_settings->num_cross = num_baselines*num_time_steps*num_freqs;
  woden_settings->num_visis = woden_settings->num_cross;
  woden_settings->coarse_band_width = 1.28e+6;
  //Make the fine channel width insanely small so beam changes little
  //with frequency - that way we can test for just one gain value per time
  woden_settings->frequency_resolution = 10e+3;
  woden_settings->latitude = LOFAR_DEC0;
  double latitudes[num_time_steps];
  for (int i = 0; i < num_time_steps; i++) {
    latitudes[i] = LOFAR_DEC0;
  }

  woden_settings->latitudes = latitudes;

  double *lsts = malloc(num_time_steps*sizeof(double));
  for (int time_ind = 0; time_ind < num_time_steps; time_ind++) {
    lsts[time_ind] = LOFAR_RA0 + 1.0*DD2R;
  }

  woden_settings->lsts = lsts;
  woden_settings->do_autos = 0;
  woden_settings->use_dipamps = 0;

  // woden_settings->mwa_dipole_amps = malloc(16*sizeof(double));
  // for (int i = 0; i < 16; i++) {
  //   woden_settings->mwa_dipole_amps[i] = 1.0;
  // }

  woden_settings->verbose = 0;
  woden_settings->off_cardinal_dipoles = 0;
  woden_settings->normalise_primary_beam = 1;
  woden_settings->single_everybeam_station = 0;

  set_mjds(woden_settings, beamtype, num_time_steps);

  set_azza_para(cropped_sky_models, num_time_steps,
                n_points, n_gauss, n_shapes, num_sources,
                beamtype);


  double base_band_freq = 120e+6;

  array_layout_t *array_layout = malloc(sizeof(array_layout_t));


  array_layout->X_diff_metres = malloc(num_time_steps*num_baselines*sizeof(double));
  array_layout->Y_diff_metres = malloc(num_time_steps*num_baselines*sizeof(double));
  array_layout->Z_diff_metres = malloc(num_time_steps*num_baselines*sizeof(double));

  for (int time_ind = 0; time_ind < num_time_steps; time_ind++) {
      for (int baseline = 0; baseline < num_baselines; baseline++) {
        int time_off = time_ind*num_baselines;
        array_layout->X_diff_metres[time_off + baseline] = (baseline + 1) * 10;
        array_layout->Y_diff_metres[time_off + baseline] = (baseline + 1) * 10;
        array_layout->Z_diff_metres[time_off + baseline] = 0.0;
      }
  }

  user_precision_t *sbf = NULL;
  if (cropped_sky_models->num_shapelets > 0) {
    sbf = malloc( sbf_N * sbf_L * sizeof(user_precision_t) );
    sbf = create_sbf(sbf);
  }

  printf("SIZE OF THING %d\n",woden_settings->num_visis );

  visibility_set_t *visibility_set = setup_visibility_set(woden_settings->num_visis);

  fill_timefreq_visibility_set(visibility_set, woden_settings,
                               base_band_freq, lsts);

  beam_settings_t *beam_settings = malloc(sizeof(beam_settings_t));
  beam_settings->beamtype = beamtype;

  // printf("Calling calculate_visibilities\n");
  calculate_visibilities(array_layout, cropped_sky_models, beam_settings,
                         woden_settings, visibility_set, sbf);

  free_visi_set_inputs(visibility_set);
  free_visi_set_outputs(visibility_set);

  free(beam_settings);
  free(woden_settings->mjds);
  free(woden_settings);

  // // for (int visi = 0; visi < woden_settings->num_visis; visi++) {
  // //   printf("\t\tVisi %d: %f %f %f %f %f %f %f %f\n", visi,
  // //               visibility_set->sum_visi_XX_real[visi],
  // //               visibility_set->sum_visi_XX_imag[visi],
  // //               visibility_set->sum_visi_XY_real[visi],
  // //               visibility_set->sum_visi_XY_imag[visi],
  // //               visibility_set->sum_visi_YX_real[visi],
  // //               visibility_set->sum_visi_YX_imag[visi],
  // //               visibility_set->sum_visi_YY_real[visi],
  // //               visibility_set->sum_visi_YY_imag[visi]);
  // // }
  //Be free my pretties!
  if (cropped_sky_models->num_shapelets > 0) {
    free(sbf);
  }

  // free_sky_model(cropped_sky_models);

  for (int cats_ind = 0; cats_ind < num_sources; cats_ind++) {
    free(cropped_sky_models->sources[cats_ind].point_components.ras);
    free(cropped_sky_models->sources[cats_ind].point_components.decs);
    free(cropped_sky_models->sources[cats_ind].point_components.para_angles);
    free(cropped_sky_models->sources[cats_ind].point_components.param_indexes);
    free(cropped_sky_models->sources[cats_ind].point_components.power_ref_freqs);
    free(cropped_sky_models->sources[cats_ind].point_components.power_ref_stokesI);
    free(cropped_sky_models->sources[cats_ind].point_components.power_SIs);
    free(cropped_sky_models->sources[cats_ind].point_components.power_comp_inds);
    free(cropped_sky_models->sources[cats_ind].point_components.azs);
    free(cropped_sky_models->sources[cats_ind].point_components.zas);
  }

  free(array_layout->X_diff_metres);
  free(array_layout->Y_diff_metres);
  free(array_layout->Z_diff_metres);
  free(lsts);

}


//Run the test with unity
int main(void)
{
  

  UNITY_BEGIN();

  RUN_TEST(do_test);

  return UNITY_END();

}
