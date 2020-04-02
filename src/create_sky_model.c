#include <stdio.h>
#include <erfa.h>
#include <math.h>
#include "constants.h"
// #include "read_and_write.h"
#include "create_sky_model.h"

void convert_radec2azza(double ra, double dec, double lst,
     double * az, double * za){

  double erfa_az, el;
  double ha = ra - lst;

  eraHd2ae( ha, dec, (double)MWA_LAT_RAD, &erfa_az, &el );

  * az = erfa_az;
  * za = M_PI / 2. - el;

}

/*********************************
// Crop a sky model contained in a source_catalogue_t struct (raw_srccat) and
// crop out all sources below the horizon (at the beginning of the observation).
// Return a new single catsource_t (cropped_srccat) that contains the full
// cropped sky model

// TODO update the sky model for every time step to account for sources that
// have risen and set during the observation?
**********************************/


catsource_t * crop_sky_model(source_catalogue_t *raw_srccat, float *lsts,
              int num_time_steps, e_sky_crop sky_crop_type){

  // //A shapelet model has no set number of basis coeffs, and so can have
  // //any number of basis function calculations - try and keep track of that
  // //using this number so we can split compute across the GPU accordingly
  // int total_comp_calculations = 0;

  for (size_t src = 0; src < raw_srccat->num_sources; src++){

    //Initialise containers for azimuth and zenith angles, used in cropping
    //and beam calculations later on

    raw_srccat->catsources[src].point_azs = malloc( raw_srccat->catsources[src].n_points * sizeof(double) );
    raw_srccat->catsources[src].point_zas = malloc( raw_srccat->catsources[src].n_points * sizeof(double) );

  }

  double az, za;

  int num_point_comp_retained = 0;
  int num_gauss_comp_retained = 0;
  int num_shape_comp_retained = 0;
  int num_shape_coeff_retained = 0;
  int all_comps_above_horizon = 0;

  //used if CROP_SOURCES
  //Just include all the sources that are totally above the horizon
  // num_sources_retained remains at 0 when all components (POINT, GAUSSIAN,
  // or SHAPELET) of a source are above horizon - goes to 1 if not
  int num_sources_retained=0;
  int *cropped_src_indexes=NULL;

  //Begin checking az/za loop here
  for (size_t src = 0; src < raw_srccat->num_sources; src++){
    all_comps_above_horizon = 0;

    //Begin point source component loop
    for (size_t point = 0; point < raw_srccat->catsources[src].n_points; point++) {
      convert_radec2azza((double)raw_srccat->catsources[src].point_ras[point],
                         (double)raw_srccat->catsources[src].point_decs[point],
                         (double)lsts[0], &az, &za);
      raw_srccat->catsources[src].point_azs[point] = az;
      raw_srccat->catsources[src].point_zas[point] = za;

      //Check if component is above horizon, and flag the source if not
      if (sky_crop_type == CROP_SOURCES) {
        if (za >= M_PI / 2.0){
          //When cropping sources, if one component is below the horizon,
          //we're dumping the whole source
          all_comps_above_horizon = 1;
        }
      }

      //Count number of components above horizon so we can allocate correct
      //amount of memory later
      else if (sky_crop_type == CROP_COMPONENTS) {
        if (za < M_PI / 2.0){
          num_point_comp_retained ++;
        }
      }

      else {
          printf("ERROR: create_sky_model.c:crop_sky_model: needs a correct sky_crop_type, allowed is 0 (CROP_SOURCES), 1 (CROP_COMPONENTS), given was %d \n", sky_crop_type );
          return NULL;
      }

    }//End point source component loop

    if (sky_crop_type == CROP_SOURCES) {
      if (all_comps_above_horizon == 0) {
        num_sources_retained ++;
        num_point_comp_retained += raw_srccat->catsources[src].n_points;
        cropped_src_indexes = realloc(cropped_src_indexes,sizeof(int)*num_sources_retained);
        cropped_src_indexes[num_sources_retained - 1] = src;
      }
    }

  }//Finish checking az/za loop here

  //TODO add in GAUSSIANs here
  //TODO add in SHAPELETs here

  //Make an empty catsource_t and malloc using the numbers gather above
  catsource_t *cropped_src=NULL;
  cropped_src = malloc(sizeof(catsource_t));

  cropped_src->n_points = num_point_comp_retained;
  cropped_src->n_gauss = num_gauss_comp_retained;
  cropped_src->n_shapes = num_shape_comp_retained;
  cropped_src->n_shape_coeffs = num_shape_coeff_retained;

  cropped_src->point_ras = malloc( num_point_comp_retained * sizeof(float) );
  cropped_src->point_decs = malloc( num_point_comp_retained * sizeof(float) );
  cropped_src->point_fluxes = malloc( num_point_comp_retained * sizeof(float) );
  cropped_src->point_freqs = malloc( num_point_comp_retained * sizeof(float) );
  cropped_src->point_azs = malloc( num_point_comp_retained * num_time_steps * sizeof(double) );
  cropped_src->point_zas = malloc( num_point_comp_retained * num_time_steps * sizeof(double) );

  // cropped_src->gauss_ras = malloc( num_gauss_comp_retained * sizeof(float) );
  // cropped_src->gauss_decs = malloc( num_gauss_comp_retained * sizeof(float) );
  // cropped_src->gauss_fluxes = malloc( num_gauss_comp_retained * sizeof(float) );
  // cropped_src->gauss_freqs = malloc( num_gauss_comp_retained * sizeof(float) );
  // cropped_src->gauss_majors = malloc( num_gauss_comp_retained * sizeof(float) );
  // cropped_src->gauss_minors = malloc( num_gauss_comp_retained * sizeof(float) );
  // cropped_src->gauss_pas = malloc( num_gauss_comp_retained * sizeof(float) );

  // cropped_src->shape_ras = malloc( num_shape_comp_retained * sizeof(float) );
  // cropped_src->shape_decs = malloc( num_shape_comp_retained * sizeof(float) );
  // cropped_src->shape_fluxes = malloc( num_shape_comp_retained * sizeof(float) );
  // cropped_src->shape_freqs = malloc( num_shape_comp_retained * sizeof(float) );
  // cropped_src->shape_majors = malloc( num_shape_comp_retained * sizeof(float) );
  // cropped_src->shape_minors = malloc( num_shape_comp_retained * sizeof(float) );
  // cropped_src->shape_pas = malloc( num_shape_comp_retained * sizeof(float) );
  // cropped_src->shape_n1s = malloc( num_shape_coeff_retained * sizeof(float) );
  // cropped_src->shape_n2s = malloc( num_shape_coeff_retained * sizeof(float) );
  // cropped_src->shape_coeffs = malloc( num_shape_coeff_retained * sizeof(float) );
  // cropped_src->shape_param_indexes = malloc( num_shape_coeff_retained * sizeof(float) );

  if (sky_crop_type == CROP_SOURCES) {
    printf("Sources retained after cropping %d\n",num_sources_retained );
    printf("Point components after cropping %d\n",num_point_comp_retained );

    int crop_component_index = 0;

    //Loop over all the retained source indexes, and add all component
    //information into cropped_src
    for (size_t retained = 0; retained < num_sources_retained; retained++) {
      int src = cropped_src_indexes[retained];

      //Loop over point components
      for (size_t point = 0; point < raw_srccat->catsources[src].n_points; point++){
        cropped_src->point_ras[crop_component_index] = raw_srccat->catsources[src].point_ras[point];
        cropped_src->point_decs[crop_component_index] = raw_srccat->catsources[src].point_decs[point];
        cropped_src->point_fluxes[crop_component_index] = raw_srccat->catsources[src].point_fluxes[point];
        cropped_src->point_freqs[crop_component_index] = raw_srccat->catsources[src].point_freqs[point];

        //Calculate az/za values for each point for all time steps
        for (size_t time_step = 0; time_step < num_time_steps; time_step++) {
          convert_radec2azza((double)raw_srccat->catsources[src].point_ras[point],
                             (double)raw_srccat->catsources[src].point_decs[point],
                             (double)lsts[time_step], &az, &za);
          cropped_src->point_azs[crop_component_index*num_time_steps + time_step] = az;
          cropped_src->point_zas[crop_component_index*num_time_steps + time_step] = za;
        }

        crop_component_index ++;

      }//End point component loop
    }//End retained sources loop
  }//End if sky_crop_type == CROP_SOURCES

  else if (sky_crop_type == CROP_COMPONENTS) {
    printf("Point components after cropping %d\n",num_point_comp_retained );
    //Keep track of sources added to cropped_src for indexing
    int crop_component_index = 0;
    //Loop over all sources in uncropped source catalogue and add all
    //components above the horizon to
    for (size_t src = 0; src < raw_srccat->num_sources; src++){
      //Begin point component loop
      for (size_t point = 0; point < raw_srccat->catsources[src].n_points; point++){

        //Check if point component above horizon
        if (raw_srccat->catsources[src].point_zas[point] < M_PI / 2.0){
          cropped_src->point_ras[crop_component_index] = raw_srccat->catsources[src].point_ras[point];
          cropped_src->point_decs[crop_component_index] = raw_srccat->catsources[src].point_decs[point];
          cropped_src->point_fluxes[crop_component_index] = raw_srccat->catsources[src].point_fluxes[point];
          cropped_src->point_freqs[crop_component_index] = raw_srccat->catsources[src].point_freqs[point];

          //Calculate az/za values for each point for all time steps
          for (size_t time_step = 0; time_step < num_time_steps; time_step++) {
            convert_radec2azza((double)raw_srccat->catsources[src].point_ras[point],
                               (double)raw_srccat->catsources[src].point_decs[point],
                               (double)lsts[time_step], &az, &za);
            cropped_src->point_azs[crop_component_index*num_time_steps + time_step] = az;
            cropped_src->point_zas[crop_component_index*num_time_steps + time_step] = za;
          }//End az/za calculation loop

          crop_component_index ++;
        }//End point above horizon loop
      }//End point component loop
    }//End raw_srccat source loop
  }//End if sky_crop_type == CROP_COMPONENTS

  return cropped_src;

}
