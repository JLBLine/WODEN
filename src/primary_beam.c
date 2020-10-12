#include "primary_beam.h"
#include "FEE_primary_beam.h"
#include <erfa.h>

void calc_para_angle(catsource_t *cropped_src, float *lsts,
                     int num_time_steps){

  cropped_src->sin_point_para_angs = malloc( cropped_src->n_points * num_time_steps * sizeof(float) );
  cropped_src->cos_point_para_angs = malloc( cropped_src->n_points * num_time_steps * sizeof(float) );
  cropped_src->sin_gauss_para_angs = malloc( cropped_src->n_gauss * num_time_steps * sizeof(float) );
  cropped_src->cos_gauss_para_angs = malloc( cropped_src->n_gauss * num_time_steps * sizeof(float) );
  cropped_src->sin_shape_para_angs = malloc( cropped_src->n_shapes * num_time_steps * sizeof(float) );
  cropped_src->cos_shape_para_angs = malloc( cropped_src->n_shapes * num_time_steps * sizeof(float) );

  for (int point = 0; point < cropped_src->n_points; point++) {

    double para_angle;
    for ( int time_step = 0; time_step < num_time_steps; time_step++ ) {

      float ha = lsts[time_step] - cropped_src->point_ras[point];
      // para_angle = eraHd2pa((double)ha, (double)cropped_src->point_decs[point], (double)MWA_LAT_RAD);

      // printf("PARA FUN %d %.8f %.8f %.8f\n",point*num_time_steps + time_step,ha,cropped_src->point_decs[point],para_angle );

      double this_ha, this_dec;

      float alt = M_PI/2 - cropped_src->point_zas[point*num_time_steps + time_step];
      eraAe2hd((double)cropped_src->point_azs[point*num_time_steps + time_step], (double)alt, (double)MWA_LAT_RAD,
               &this_ha, &this_dec);

      para_angle = eraHd2pa(this_ha, this_dec, (double)MWA_LAT_RAD);

      cropped_src->sin_point_para_angs[point*num_time_steps + time_step] = sinf((float)para_angle + M_PI/2.0);
      cropped_src->cos_point_para_angs[point*num_time_steps + time_step] = cosf((float)para_angle + M_PI/2.0);

      // printf("PARA FUN %d %.8f %.8f %.8f\n",point*num_time_steps + time_step,
      //       para_angle,sinf((float)para_angle + M_PI/2.0),
      //       cosf((float)para_angle + M_PI/2.0) );



    }
  }//END point component loop

  // for (size_t i = 0; i < cropped_src->n_points; i++) {
  //   printf("RIGHT HERE 1 %.5f %.5f\n",cropped_src->sin_point_para_angs[i],cropped_src->cos_point_para_angs[i] );
  // }

  for (int gauss = 0; gauss < cropped_src->n_gauss; gauss++) {

    double para_angle;
    for ( int time_step = 0; time_step < num_time_steps; time_step++ ) {

      float ha = lsts[time_step] - cropped_src->gauss_ras[gauss];
      para_angle = eraHd2pa((double)ha, (double)cropped_src->gauss_decs[gauss], (double)MWA_LAT_RAD);

      cropped_src->sin_gauss_para_angs[gauss*num_time_steps + time_step] = sinf((float)para_angle + M_PI/2.0);
      cropped_src->cos_gauss_para_angs[gauss*num_time_steps + time_step] = cosf((float)para_angle + M_PI/2.0);
    }
  }//END gauss component loop

  for (int shape = 0; shape < cropped_src->n_shapes; shape++) {
    double para_angle;
    for ( int time_step = 0; time_step < num_time_steps; time_step++ ) {

      float ha = lsts[time_step] - cropped_src->shape_ras[shape];
      para_angle = eraHd2pa((double)ha, (double)cropped_src->shape_decs[shape], (double)MWA_LAT_RAD);

      cropped_src->sin_shape_para_angs[shape*num_time_steps + time_step] = sinf((float)para_angle + M_PI/2.0);
      cropped_src->cos_shape_para_angs[shape*num_time_steps + time_step] = cosf((float)para_angle + M_PI/2.0);
    }
  }//END shape component loop

}

beam_settings_t fill_primary_beam_settings(woden_settings_t *woden_settings,
                MetaFfile_t metafits, catsource_t *cropped_src,
                float *lsts, int num_time_steps) {


  //Setup primary beam settings for observation
  beam_settings_t beam_settings; //= malloc(sizeof(beam_settings_t));
  //Angles used in calculating beam style l,m,ns
  beam_settings.beam_angles_array = malloc(3*sizeof(float));
  beam_settings.beam_angles_array[0] = sinf(metafits.dec_point);
  beam_settings.beam_angles_array[1] = cosf(metafits.dec_point);
  beam_settings.beam_angles_array[2] = woden_settings->lst_base - metafits.ra_point;

  //Number of beam calculations needed for point components
  beam_settings.num_point_beam_values = cropped_src->n_points * woden_settings->num_time_steps * woden_settings->num_freqs;
  beam_settings.num_gausscomp_beam_values = cropped_src->n_gauss * woden_settings->num_time_steps * woden_settings->num_freqs;
  beam_settings.num_shape_beam_values = cropped_src->n_shapes * woden_settings->num_time_steps * woden_settings->num_freqs;


  if (woden_settings->beamtype == GAUSS_BEAM) {
    beam_settings.beamtype = GAUSS_BEAM;

    printf("Setting up Gaussian primary beam settings\n");
    printf("   setting beam FWHM to %.5fdeg and ref freq to %.3fMHz\n",
            woden_settings->gauss_beam_FWHM,woden_settings->gauss_beam_ref_freq / 1e+6  );

    //Set constants used in beam calculation
    beam_settings.beam_FWHM_rad = woden_settings->gauss_beam_FWHM * D2R;
    //TODO I cannot for the life of me work out how to cudaMalloc and Memcpy
    //a single float (argh) so put the ref freq in an array (embarrassment)
    float beam_ref_freq_array[1] = {woden_settings->gauss_beam_ref_freq};
    beam_settings.beam_ref_freq_array = malloc(sizeof(float));
    beam_settings.beam_ref_freq_array = beam_ref_freq_array;

    //Store all ha (which change with lst) that the beam needs to be calculated at.
    beam_settings.beam_point_has = malloc(woden_settings->num_time_steps * cropped_src->n_points * sizeof(float));
    beam_settings.beam_point_decs = malloc(woden_settings->num_time_steps * cropped_src->n_points * sizeof(float));

    beam_settings.beam_gausscomp_has = malloc(woden_settings->num_time_steps * cropped_src->n_gauss * sizeof(float));
    beam_settings.beam_gausscomp_decs = malloc(woden_settings->num_time_steps * cropped_src->n_gauss * sizeof(float));

    beam_settings.beam_shape_has = malloc(woden_settings->num_time_steps * cropped_src->n_shapes * sizeof(float));
    beam_settings.beam_shape_decs = malloc(woden_settings->num_time_steps * cropped_src->n_shapes * sizeof(float));

    //Loop over all time and point components and calculate ha
    for ( int time_step = 0; time_step < woden_settings->num_time_steps; time_step++ ) {
      for (int component = 0; component < cropped_src->n_points; component++) {
        int step = cropped_src->n_points*time_step + component;

        beam_settings.beam_point_has[step] = lsts[time_step] - cropped_src->point_ras[component];
        beam_settings.beam_point_decs[step] = cropped_src->point_decs[component];

        // printf("Original mothers step, ha %d %f \n",step,beam_settings.beam_point_has[step] );

      }//point loop

    //Loop over all time and gausscomp components and calculate ha
      for (int component = 0; component < cropped_src->n_gauss; component++) {
        int step = cropped_src->n_gauss*time_step + component;

        beam_settings.beam_gausscomp_has[step] = lsts[time_step] - cropped_src->gauss_ras[component];
        beam_settings.beam_gausscomp_decs[step] = cropped_src->gauss_decs[component];
      }//gausscomp loop

    //Loop over all time and shape components and calculate ha
      for (int component = 0; component < cropped_src->n_shapes; component++) {
        int step = cropped_src->n_shapes*time_step + component;

        beam_settings.beam_shape_has[step] = lsts[time_step] - cropped_src->shape_ras[component];
        beam_settings.beam_shape_decs[step] = cropped_src->shape_decs[component];
      }//shape loop

    }//gaussian beam time loop
  } // End if (woden_settings->gaussian_beam)

  else if (woden_settings->beamtype == FEE_BEAM) {
    beam_settings.beamtype = FEE_BEAM;

    //Get the parallactic angle of the zenith for every time step
    //Need to rotate the FEE model which is stored in theta/phi pols by the
    //parallactic angle to obtain XX/YY
    //There are 4 normalisations to calculate, so need 4 times num time steps
    beam_settings.para_cosrot = malloc(woden_settings->num_time_steps*MAX_POLS*sizeof(float));
    beam_settings.para_sinrot = malloc(woden_settings->num_time_steps*MAX_POLS*sizeof(float));

    double para_angle;
    for ( int time_step = 0; time_step < woden_settings->num_time_steps; time_step++ ) {

      float zenith_HA = 0.0;
      para_angle = eraHd2pa((double)zenith_HA, (double)MWA_LAT_RAD, (double)MWA_LAT_RAD);

      for (size_t pol_direction = 0; pol_direction < MAX_POLS; pol_direction++) {
        beam_settings.para_cosrot[time_step*MAX_POLS + pol_direction] = cosf((float)para_angle + M_PI/2.0);
        beam_settings.para_sinrot[time_step*MAX_POLS + pol_direction] = sinf((float)para_angle + M_PI/2.0);
      }
      // printf("PARA ANGLEEEEE %.10f %.10f %.10f\n",para_angle,cosf((float)para_angle + M_PI/2.0),sinf((float)para_angle + M_PI/2.0) );
    }

    calc_para_angle(cropped_src, lsts, num_time_steps);
  }

  else if (woden_settings->beamtype == ANALY_DIPOLE) {
    beam_settings.beamtype = ANALY_DIPOLE;
  }

  else {
    printf("BEAM TYPE %d\n",(int)woden_settings->beamtype );
  }

  return beam_settings;

} // end of fill_primary_beam_settings()

void setup_FEE_beam(woden_settings_t *woden_settings, MetaFfile_t metafits,
                    beam_settings_t beam_settings, float base_middle_freq){

  //Just use one single tile beam for all for now - will need a certain
  //number in the future to include dipole flagging
  int st = 0;
  beam_settings.FEE_beam = malloc(sizeof(copy_primary_beam_t));
  //We need the zenith beam to get the normalisation
  beam_settings.FEE_beam_zenith = malloc(sizeof(copy_primary_beam_t));

  printf("Middle freq is %f\n",base_middle_freq );

  float float_zenith_delays[16] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  printf("Setting up the zenith FEE beam...");
  RTS_HDFBeamInit(woden_settings->hdf5_beam_path, base_middle_freq, beam_settings.FEE_beam_zenith, float_zenith_delays, st);
  printf(" done.\n");

  printf("Getting FEE beam normalisation...");
  get_HDFBeam_normalisation(beam_settings, woden_settings->num_time_steps);
  printf(" done.\n");

  float *float_delays = NULL;
  float_delays = malloc(16*sizeof(float));

  for (size_t i = 0; i < 16; i++) {
   float_delays[i] = metafits.FEE_ideal_delays[i];
  }

  printf("Setting up the FEE beam...");
  RTS_HDFBeamInit(woden_settings->hdf5_beam_path, base_middle_freq, beam_settings.FEE_beam, float_delays, st);
  printf(" done.\n");

  printf("Copying the FEE beam across to the GPU...");
  copy_FEE_primary_beam_to_GPU(beam_settings, woden_settings->num_time_steps);
  printf(" done.\n");

}
