#include "primary_beam.h"
#include "FEE_primary_beam.h"
#include <erfa.h>
#include <complex.h>

void calc_para_angle(catsource_t *cropped_src, user_precision_t *lsts,
                     double latitude, int num_time_steps){

  //At the moment, we calculate one beam value in frequency for MWA FEE as the resolution
  //is 1.28MHz. We should interpolate this in the future (might be memory hungry)
  //So only need to malloc enough para angles to cover time, and not freq
  cropped_src->sin_point_para_angs = malloc( cropped_src->n_points * num_time_steps * sizeof(user_precision_t) );
  cropped_src->cos_point_para_angs = malloc( cropped_src->n_points * num_time_steps * sizeof(user_precision_t) );
  cropped_src->sin_gauss_para_angs = malloc( cropped_src->n_gauss * num_time_steps * sizeof(user_precision_t) );
  cropped_src->cos_gauss_para_angs = malloc( cropped_src->n_gauss * num_time_steps * sizeof(user_precision_t) );
  cropped_src->sin_shape_para_angs = malloc( cropped_src->n_shapes * num_time_steps * sizeof(user_precision_t) );
  cropped_src->cos_shape_para_angs = malloc( cropped_src->n_shapes * num_time_steps * sizeof(user_precision_t) );

  for (int point = 0; point < cropped_src->n_points; point++) {

    double para_angle;
    for ( int time_step = 0; time_step < num_time_steps; time_step++ ) {

      double ha = (double)lsts[time_step] - cropped_src->point_ras[point];
      para_angle = eraHd2pa(ha, cropped_src->point_decs[point], latitude);

      cropped_src->sin_point_para_angs[point*num_time_steps + time_step] = sin((user_precision_t)para_angle + M_PI/2.0);
      cropped_src->cos_point_para_angs[point*num_time_steps + time_step] = cos((user_precision_t)para_angle + M_PI/2.0);

    }
  }//END point component loop

  for (int gauss = 0; gauss < cropped_src->n_gauss; gauss++) {

    double para_angle;
    for ( int time_step = 0; time_step < num_time_steps; time_step++ ) {

      double ha = (double)lsts[time_step] - cropped_src->gauss_ras[gauss];
      para_angle = eraHd2pa(ha, cropped_src->gauss_decs[gauss], latitude);

      cropped_src->sin_gauss_para_angs[gauss*num_time_steps + time_step] = sin((user_precision_t)para_angle + M_PI/2.0);
      cropped_src->cos_gauss_para_angs[gauss*num_time_steps + time_step] = cos((user_precision_t)para_angle + M_PI/2.0);
    }
  }//END gauss component loop

  for (int shape = 0; shape < cropped_src->n_shapes; shape++) {
    double para_angle;
    for ( int time_step = 0; time_step < num_time_steps; time_step++ ) {

      double ha = (double)lsts[time_step] - cropped_src->shape_ras[shape];
      para_angle = eraHd2pa((double)ha, cropped_src->shape_decs[shape], latitude);

      cropped_src->sin_shape_para_angs[shape*num_time_steps + time_step] = sin((user_precision_t)para_angle + M_PI/2.0);
      cropped_src->cos_shape_para_angs[shape*num_time_steps + time_step] = cos((user_precision_t)para_angle + M_PI/2.0);
    }
  }//END shape component loop
}

beam_settings_t * fill_primary_beam_settings(woden_settings_t *woden_settings,
                            catsource_t *cropped_src, user_precision_t *lsts) {


  //Setup primary beam settings for observation
  beam_settings_t *beam_settings = malloc(sizeof(beam_settings_t));

  //Number of beam calculations needed for point components

  // beam_settings->num_point_primarybeam_values = cropped_src->n_points * woden_settings->num_time_steps * woden_settings->num_freqs;
  // beam_settings->num_gauss_primarybeam_values = cropped_src->n_gauss * woden_settings->num_time_steps * woden_settings->num_freqs;
  // beam_settings->num_shape_primarybeam_values = cropped_src->n_shapes * woden_settings->num_time_steps * woden_settings->num_freqs;


  if (woden_settings->beamtype == GAUSS_BEAM) {
    beam_settings->beamtype = GAUSS_BEAM;

    //Angles used in calculating beam centred l,m,ns
    beam_settings->gauss_sdec = sin(woden_settings->gauss_dec_point);
    beam_settings->gauss_cdec = cos(woden_settings->gauss_dec_point);
    beam_settings->gauss_ha = woden_settings->lst_base - woden_settings->gauss_ra_point;

    printf("Setting up Gaussian primary beam settings\n");
    printf("   pointing at HA, Dec = %.5fdeg, %.5fdeg\n",
               beam_settings->gauss_ha/DD2R, woden_settings->gauss_dec_point/DD2R );
    printf("   setting beam FWHM to %.5fdeg and ref freq to %.3fMHz\n",
            woden_settings->gauss_beam_FWHM,woden_settings->gauss_beam_ref_freq / 1e+6  );

    //Set constants used in beam calculation
    beam_settings->beam_FWHM_rad = woden_settings->gauss_beam_FWHM * DD2R;
    beam_settings->beam_ref_freq = woden_settings->gauss_beam_ref_freq;

    //Store all ha (which change with lst) that the beam needs to be calculated at.
    cropped_src->point_gaussbeam_has = malloc(woden_settings->num_time_steps * cropped_src->n_points * sizeof(double));
    cropped_src->point_gaussbeam_decs = malloc(woden_settings->num_time_steps * cropped_src->n_points * sizeof(double));

    cropped_src->gauss_gaussbeam_has = malloc(woden_settings->num_time_steps * cropped_src->n_gauss * sizeof(double));
    cropped_src->gauss_gaussbeam_decs = malloc(woden_settings->num_time_steps * cropped_src->n_gauss * sizeof(double));

    cropped_src->shape_gaussbeam_has = malloc(woden_settings->num_time_steps * cropped_src->n_shapes * sizeof(double));
    cropped_src->shape_gaussbeam_decs = malloc(woden_settings->num_time_steps * cropped_src->n_shapes * sizeof(double));

    //Loop over all time and point components and calculate ha
    for (int component = 0; component < cropped_src->n_points; component++) {
      for ( int time_step = 0; time_step < woden_settings->num_time_steps; time_step++ ) {
        int step = component*woden_settings->num_time_steps + time_step;
        cropped_src->point_gaussbeam_has[step] = (double)lsts[time_step] - cropped_src->point_ras[component];
        cropped_src->point_gaussbeam_decs[step] = cropped_src->point_decs[component];
        // printf("THIS THING %.6f %.6f %.8f\n",cropped_src->point_ras[component],
        //                                      cropped_src->point_decs[component], MWA_LAT_RAD );
      }
    }//point loop

    //Loop over all time and gauss components and calculate ha
    for (int component = 0; component < cropped_src->n_gauss; component++) {
      for ( int time_step = 0; time_step < woden_settings->num_time_steps; time_step++ ) {
        int step = component*woden_settings->num_time_steps + time_step;
        cropped_src->gauss_gaussbeam_has[step] = (double)lsts[time_step] - cropped_src->gauss_ras[component];
        cropped_src->gauss_gaussbeam_decs[step] = cropped_src->gauss_decs[component];
      }
    }//gausscomp loop

    //Loop over all time and shape components and calculate ha
    for (int component = 0; component < cropped_src->n_shapes; component++) {
      for ( int time_step = 0; time_step < woden_settings->num_time_steps; time_step++ ) {
        int step = component*woden_settings->num_time_steps + time_step;
        cropped_src->shape_gaussbeam_has[step] = (double)lsts[time_step] - cropped_src->shape_ras[component];
        cropped_src->shape_gaussbeam_decs[step] = cropped_src->shape_decs[component];
      }
    }//shape loop

  } // End if (woden_settings->gaussian_beam)

  else if (woden_settings->beamtype == FEE_BEAM) {
    beam_settings->beamtype = FEE_BEAM;

    //Need to rotate the FEE model which is stored in theta/phi pols by the
    //parallactic angle to obtain XX/YY
    calc_para_angle(cropped_src, lsts, (double)woden_settings->latitude, woden_settings->num_time_steps);
  }

  else if (woden_settings->beamtype == ANALY_DIPOLE) {
    beam_settings->beamtype = ANALY_DIPOLE;
  }

  else {
    printf("NO PRIMARY BEAM HAS BEEN SELECTED\n\tWill run without a primary beam\n");
    beam_settings->beamtype = NO_BEAM;
  }

  return beam_settings;

} // end of fill_primary_beam_settings()
