#include "primary_beam.h"
#include <complex.h>
#include <math.h>

beam_settings_t * fill_primary_beam_settings(woden_settings_t *woden_settings,
                            source_t *cropped_src, double *lsts) {


  //Setup primary beam settings for observation
  beam_settings_t *beam_settings = malloc(sizeof(beam_settings_t));

  //Number of beam calculations needed for point components

  // beam_settings->num_point_primarybeam_values = cropped_src->n_points * woden_settings->num_time_steps * woden_settings->num_freqs;
  // beam_settings->num_gauss_primarybeam_values = cropped_src->n_gauss * woden_settings->num_time_steps * woden_settings->num_freqs;
  // beam_settings->num_shape_primarybeam_values = cropped_src->n_shapes * woden_settings->num_time_steps * woden_settings->num_freqs;

  //Both GAUSS_BEAM and MWA_ANALY need hour angles and declinations for
  //beam calculations so set them up here
  if (woden_settings->beamtype == GAUSS_BEAM || woden_settings->beamtype == MWA_ANALY) {
    if (woden_settings->beamtype == GAUSS_BEAM){
      //Extra settings that just GAUSS_BEAM needs
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

    } else {
      beam_settings->beamtype = MWA_ANALY;
      printf("Setting up analytic MWA primary beam settings\n");
    }

    //Store all ha (which change with lst) that the beam needs to be calculated at.
    cropped_src->point_components.beam_has = malloc(woden_settings->num_time_steps * cropped_src->n_points * sizeof(double));
    cropped_src->point_components.beam_decs = malloc(woden_settings->num_time_steps * cropped_src->n_points * sizeof(double));

    cropped_src->gauss_components.beam_has = malloc(woden_settings->num_time_steps * cropped_src->n_gauss * sizeof(double));
    cropped_src->gauss_components.beam_decs = malloc(woden_settings->num_time_steps * cropped_src->n_gauss * sizeof(double));

    cropped_src->shape_components.beam_has = malloc(woden_settings->num_time_steps * cropped_src->n_shapes * sizeof(double));
    cropped_src->shape_components.beam_decs = malloc(woden_settings->num_time_steps * cropped_src->n_shapes * sizeof(double));

    //Loop over all time and point components and calculate ha
    for (int component = 0; component < cropped_src->n_points; component++) {
      for ( int time_step = 0; time_step < woden_settings->num_time_steps; time_step++ ) {
        int step = component*woden_settings->num_time_steps + time_step;
        cropped_src->point_components.beam_has[step] = lsts[time_step] - cropped_src->point_components.ras[component];
        cropped_src->point_components.beam_decs[step] = cropped_src->point_components.decs[component];
      }
    }//point loop

    //Loop over all time and gauss components and calculate ha
    for (int component = 0; component < cropped_src->n_gauss; component++) {
      for ( int time_step = 0; time_step < woden_settings->num_time_steps; time_step++ ) {
        int step = component*woden_settings->num_time_steps + time_step;
        cropped_src->gauss_components.beam_has[step] = lsts[time_step] - cropped_src->gauss_components.ras[component];
        cropped_src->gauss_components.beam_decs[step] = cropped_src->gauss_components.decs[component];
      }
    }//gausscomp loop

    //Loop over all time and shape components and calculate ha
    for (int component = 0; component < cropped_src->n_shapes; component++) {
      for ( int time_step = 0; time_step < woden_settings->num_time_steps; time_step++ ) {
        int step = component*woden_settings->num_time_steps + time_step;
        cropped_src->shape_components.beam_has[step] = lsts[time_step] - cropped_src->shape_components.ras[component];
        cropped_src->shape_components.beam_decs[step] = cropped_src->shape_components.decs[component];
      }
    }//shape loop

  } // End if (woden_settings->gaussian_beam)

  else if (woden_settings->beamtype == FEE_BEAM) {
    beam_settings->beamtype = FEE_BEAM;
  }

  else if (woden_settings->beamtype == FEE_BEAM_INTERP) {
    beam_settings->beamtype = FEE_BEAM_INTERP;
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
