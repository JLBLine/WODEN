#include <stdio.h>
#include <erfa.h>
#include <math.h>
#include "constants.h"
// #include "read_and_write.h"
#include "create_sky_model.h"


/*********************************
// Takes an ra, dec, lst (all rad) and returns the azimuth and zenith angle
// assuming the latitude of the MWA. Uses the ERFA library to do the
// transformation using the ha.
**********************************/
void convert_radec2azza(double ra, double dec, double lst,
     double * az, double * za){

  double erfa_az, el;
  double ha = ra - lst;

  eraHd2ae( ha, dec, (double)MWA_LAT_RAD, &erfa_az, &el );

  * az = erfa_az;
  * za = M_PI / 2. - el;

}

/*********************************
// Takes a zenith angle and checks if it's below the horizon. Depending on
// whether we are cropping all SOURCEs partially below the horizon, or
// retaining all COMPONENTs, update either all_comps_above_horizon or
// sky_crop_type. If the component is a shapelet, we also need to cycle
// through the catsource.shape_param_indexes, checking which ones match this
// shape index (int shape), and update num_shape_coeff_retained accordingly
**********************************/
void horizon_test(double za, e_sky_crop sky_crop_type,
     e_horizon * all_comps_above_horizon, int * num_comp_retained,
     int * num_shape_coeff_retained, int num_shape_coeff_component,
     float *shape_param_indexes, int shape){
  //Check if component is above horizon, and flag the source if not
  if (sky_crop_type == CROP_SOURCES) {
    if (za >= M_PI / 2.0){
      //When cropping sources, if one component is below the horizon,
      //we're dumping the whole source
       * all_comps_above_horizon = BELOW;
    }
  }

  //Count number of components above horizon so we can allocate correct
  //amount of memory later
  else if (sky_crop_type == CROP_COMPONENTS) {
    if (za < M_PI / 2.0){
      * num_comp_retained += 1;
      //If this is a shapelet component, num_shape_coeff_component will be
      //non-zero. If so, we need to add the number of coeffs to the total
      //to go into the cropped source
      if (num_shape_coeff_component > 0) {
        for (size_t param_index = 0; param_index < num_shape_coeff_component; param_index++) {
          if ( (int)shape_param_indexes[param_index] == shape ){
            * num_shape_coeff_retained += 1;
          }
        }
      }
    }
  }
  //Something has gone terribly wrong if we get to the else statement
  else {
    printf("ERROR: create_sky_model.c:crop_sky_model: needs a correct sky_crop_type, allowed is 0 (CROP_SOURCES), 1 (CROP_COMPONENTS), given was %d \n", sky_crop_type );
  }

}

/*********************************
// Crop a sky model contained in a source_catalogue_t struct (raw_srccat) and
// crop out all sources below the horizon (at the beginning of the observation).
// Return a new single catsource_t (cropped_srccat) that contains the full
// cropped sky model
// First part of the function calculates az/za for the initial time step, and
// counts how many components / sources are to be saved
// Second part mallocs a big enough catsource_t struct to contain all cropped
// components, and copies all relevant data across from raw_srccat

// Possible TODO update the sky model for every time step to account for sources that
// have risen and set during the observation?

// Possible TODO remove the az/za calculation for all time steps, and do in a
// separate function to make it more modular and clear
**********************************/


catsource_t * crop_sky_model(source_catalogue_t *raw_srccat, float *lsts,
              int num_time_steps, e_sky_crop sky_crop_type){

  for (size_t src = 0; src < raw_srccat->num_sources; src++){

    //Initialise containers for azimuth and zenith angles, used in cropping
    //later on

    raw_srccat->catsources[src].point_azs = malloc( raw_srccat->catsources[src].n_points * sizeof(double) );
    raw_srccat->catsources[src].point_zas = malloc( raw_srccat->catsources[src].n_points * sizeof(double) );

    raw_srccat->catsources[src].gauss_azs = malloc( raw_srccat->catsources[src].n_gauss * sizeof(double) );
    raw_srccat->catsources[src].gauss_zas = malloc( raw_srccat->catsources[src].n_gauss * sizeof(double) );

    raw_srccat->catsources[src].shape_azs = malloc( raw_srccat->catsources[src].n_shapes * sizeof(double) );
    raw_srccat->catsources[src].shape_zas = malloc( raw_srccat->catsources[src].n_shapes * sizeof(double) );

  }

  double az, za;
  // Counters for different component types that survive cropping
  int num_point_comp_retained = 0;
  int num_gauss_comp_retained = 0;
  int num_shape_comp_retained = 0;
  int num_shape_coeff_retained = 0;

  //used if CROP_SOURCES
  // Just include all the sources that are totally above the horizon
  // num_sources_retained counts how many sources survive cropping
  // all_comps_above_horizon is used to record
  int num_sources_retained=0;
  int *cropped_src_indexes=NULL;
  e_horizon all_comps_above_horizon;

  //Begin checking az/za loop here
  for (size_t src = 0; src < raw_srccat->num_sources; src++){
    all_comps_above_horizon = ABOVE;

    //Begin point source component loop
    for (size_t point = 0; point < raw_srccat->catsources[src].n_points; point++) {
      //Calculate az/za for all point components
      convert_radec2azza((double)raw_srccat->catsources[src].point_ras[point],
                         (double)raw_srccat->catsources[src].point_decs[point],
                         (double)lsts[0], &az, &za);
      raw_srccat->catsources[src].point_azs[point] = az;
      raw_srccat->catsources[src].point_zas[point] = za;
      //Check if components are above the horizon, and count how many
      //components survive / flag a source if a component is below the horizon
      //Last three arguments only used for shapelets so pass 0
      horizon_test(za, sky_crop_type, &all_comps_above_horizon,
                  &num_point_comp_retained, &num_shape_coeff_retained,
                  0, raw_srccat->catsources[src].shape_param_indexes, 0);

    }//End point source component loop

    //Begin gauss source component loop
    for (size_t gauss = 0; gauss < raw_srccat->catsources[src].n_gauss; gauss++) {
        //Calculate az/za for all gauss components
      convert_radec2azza((double)raw_srccat->catsources[src].gauss_ras[gauss],
                         (double)raw_srccat->catsources[src].gauss_decs[gauss],
                         (double)lsts[0], &az, &za);
      raw_srccat->catsources[src].gauss_azs[gauss] = az;
      raw_srccat->catsources[src].gauss_zas[gauss] = za;
      //Check if components are above the horizon, and count how many
      //components survive / flag a source if a component is below the horizon
      //Last three arguments only used for shapelets so pass 0
      horizon_test(za, sky_crop_type, &all_comps_above_horizon,
                  &num_gauss_comp_retained, &num_shape_coeff_retained,
                  0, raw_srccat->catsources[src].shape_param_indexes, 0);

    }//End gauss source component loop

    //Begin shape source component loop
    for (size_t shape = 0; shape < raw_srccat->catsources[src].n_shapes; shape++) {
      //Calculate az/za for all gauss components
      convert_radec2azza((double)raw_srccat->catsources[src].shape_ras[shape],
                         (double)raw_srccat->catsources[src].shape_decs[shape],
                         (double)lsts[0], &az, &za);
      raw_srccat->catsources[src].shape_azs[shape] = az;
      raw_srccat->catsources[src].shape_zas[shape] = za;
      //Check if components are above the horizon, and count how many
      //components survive / flag a source if a component is below the horizon
      //Use last three arguments to correctly identify which shapelet coeffs
      //belong to which shapelet component so we can malloc the correctly later
      horizon_test(za, sky_crop_type, &all_comps_above_horizon,
                  &num_shape_comp_retained, &num_shape_coeff_retained,
                  raw_srccat->catsources[src].n_shape_coeffs,
                  raw_srccat->catsources[src].shape_param_indexes, shape);

    }//End shape source component loop

    //After checking all components in a source, if cropping out sources,
    //check all components were above horizon. If so, update component type
    //counters, and mark down the index of the sources that survied in
    //cropped_src_indexes
    if (sky_crop_type == CROP_SOURCES) {
      if (all_comps_above_horizon == ABOVE) {
        num_sources_retained ++;
        num_point_comp_retained += raw_srccat->catsources[src].n_points;
        num_gauss_comp_retained += raw_srccat->catsources[src].n_gauss;
        num_shape_comp_retained += raw_srccat->catsources[src].n_shapes;
        num_shape_coeff_retained += raw_srccat->catsources[src].n_shape_coeffs;
        cropped_src_indexes = realloc(cropped_src_indexes,sizeof(int)*num_sources_retained);
        cropped_src_indexes[num_sources_retained - 1] = src;
      }//end if all_comps_above_horizon == ABOVE
    }//end if sky_crop_type == CROP_SOURCES

  }//Finish checking az/za loop here

  //Make an empty catsource_t and malloc using the numbers gather above
  catsource_t *cropped_src=NULL;
  cropped_src = malloc(sizeof(catsource_t));

  cropped_src->n_points = num_point_comp_retained;
  cropped_src->n_gauss = num_gauss_comp_retained;
  cropped_src->n_shapes = num_shape_comp_retained;
  cropped_src->n_shape_coeffs = num_shape_coeff_retained;

  cropped_src->point_ras = malloc( num_point_comp_retained * sizeof(float) );
  cropped_src->point_decs = malloc( num_point_comp_retained * sizeof(float) );
  cropped_src->point_ref_stokesI = malloc( num_point_comp_retained * sizeof(float) );
  cropped_src->point_ref_stokesQ = malloc( num_point_comp_retained * sizeof(float) );
  cropped_src->point_ref_stokesU = malloc( num_point_comp_retained * sizeof(float) );
  cropped_src->point_ref_stokesV = malloc( num_point_comp_retained * sizeof(float) );
  cropped_src->point_ref_freqs = malloc( num_point_comp_retained * sizeof(float) );
  cropped_src->point_SIs = malloc( num_point_comp_retained * sizeof(float) );

  cropped_src->point_azs = malloc( num_point_comp_retained * num_time_steps * sizeof(float) );
  cropped_src->point_zas = malloc( num_point_comp_retained * num_time_steps * sizeof(float) );

  cropped_src->gauss_ras = malloc( num_gauss_comp_retained * sizeof(float) );
  cropped_src->gauss_decs = malloc( num_gauss_comp_retained * sizeof(float) );
  cropped_src->gauss_ref_stokesI = malloc( num_gauss_comp_retained * sizeof(float) );
  cropped_src->gauss_ref_stokesQ = malloc( num_gauss_comp_retained * sizeof(float) );
  cropped_src->gauss_ref_stokesU = malloc( num_gauss_comp_retained * sizeof(float) );
  cropped_src->gauss_ref_stokesV = malloc( num_gauss_comp_retained * sizeof(float) );
  cropped_src->gauss_ref_freqs = malloc( num_gauss_comp_retained * sizeof(float) );
  cropped_src->gauss_SIs = malloc( num_gauss_comp_retained * sizeof(float) );
  cropped_src->gauss_majors = malloc( num_gauss_comp_retained * sizeof(float) );
  cropped_src->gauss_minors = malloc( num_gauss_comp_retained * sizeof(float) );
  cropped_src->gauss_pas = malloc( num_gauss_comp_retained * sizeof(float) );
  cropped_src->gauss_azs = malloc( num_gauss_comp_retained * num_time_steps * sizeof(float) );
  cropped_src->gauss_zas = malloc( num_gauss_comp_retained * num_time_steps * sizeof(float) );

  cropped_src->shape_ras = malloc( num_shape_comp_retained * sizeof(float) );
  cropped_src->shape_decs = malloc( num_shape_comp_retained * sizeof(float) );
  cropped_src->shape_fluxes = malloc( num_shape_comp_retained * sizeof(float) );
  cropped_src->shape_freqs = malloc( num_shape_comp_retained * sizeof(float) );
  cropped_src->shape_majors = malloc( num_shape_comp_retained * sizeof(float) );
  cropped_src->shape_minors = malloc( num_shape_comp_retained * sizeof(float) );
  cropped_src->shape_pas = malloc( num_shape_comp_retained * sizeof(float) );

  cropped_src->shape_n1s = malloc( num_shape_coeff_retained * sizeof(float) );
  cropped_src->shape_n2s = malloc( num_shape_coeff_retained * sizeof(float) );
  cropped_src->shape_coeffs = malloc( num_shape_coeff_retained * sizeof(float) );
  cropped_src->shape_param_indexes = malloc( num_shape_coeff_retained * sizeof(float) );
  cropped_src->shape_azs = malloc( num_shape_comp_retained * num_time_steps * sizeof(float) );
  cropped_src->shape_zas = malloc( num_shape_comp_retained * num_time_steps * sizeof(float) );

  printf("Num shapelets, shape coeffs %d, %d\n",num_shape_comp_retained, num_shape_coeff_retained );

  //Now add information into cropped_src
  if (sky_crop_type == CROP_SOURCES) {
    printf("Sources retained after cropping %d\n",num_sources_retained );
    printf("Point components after cropping %d\n",num_point_comp_retained );
    printf("Gaussian components after cropping %d\n",num_gauss_comp_retained );
    printf("Shapelet components after cropping %d\n",num_shape_comp_retained );
    printf("Shapelet coefficients after cropping %d\n",num_shape_coeff_retained );

    int point_crop_component_index = 0;
    int gauss_crop_component_index = 0;
    int shape_crop_component_index = 0;
    int shape_coeff_component_index = 0;

    //Loop over all the retained source indexes, and add all component
    //information into cropped_src
    for (size_t retained = 0; retained < num_sources_retained; retained++) {
      int src = cropped_src_indexes[retained];

      //Loop over point components
      for (size_t point = 0; point < raw_srccat->catsources[src].n_points; point++){
        cropped_src->point_ras[point_crop_component_index] = raw_srccat->catsources[src].point_ras[point];
        cropped_src->point_decs[point_crop_component_index] = raw_srccat->catsources[src].point_decs[point];
        cropped_src->point_ref_stokesI[point_crop_component_index] = raw_srccat->catsources[src].point_ref_stokesI[point];
        cropped_src->point_ref_stokesQ[point_crop_component_index] = raw_srccat->catsources[src].point_ref_stokesQ[point];
        cropped_src->point_ref_stokesU[point_crop_component_index] = raw_srccat->catsources[src].point_ref_stokesU[point];
        cropped_src->point_ref_stokesV[point_crop_component_index] = raw_srccat->catsources[src].point_ref_stokesV[point];
        cropped_src->point_ref_freqs[point_crop_component_index] = raw_srccat->catsources[src].point_ref_freqs[point];
        cropped_src->point_SIs[point_crop_component_index] = raw_srccat->catsources[src].point_SIs[point];

        //Calculate az/za values for each point for all time steps
        for (size_t time_step = 0; time_step < num_time_steps; time_step++) {
          convert_radec2azza((double)raw_srccat->catsources[src].point_ras[point],
                             (double)raw_srccat->catsources[src].point_decs[point],
                             (double)lsts[time_step], &az, &za);
          cropped_src->point_azs[point_crop_component_index*num_time_steps + time_step] = (float)az;
          cropped_src->point_zas[point_crop_component_index*num_time_steps + time_step] = (float)za;
        }

        point_crop_component_index ++;

      }//End point component loop

      //Loop over gauss components
      for (size_t gauss = 0; gauss < raw_srccat->catsources[src].n_gauss; gauss++){
        cropped_src->gauss_ras[gauss_crop_component_index] = raw_srccat->catsources[src].gauss_ras[gauss];
        cropped_src->gauss_decs[gauss_crop_component_index] = raw_srccat->catsources[src].gauss_decs[gauss];
        cropped_src->gauss_ref_stokesI[gauss_crop_component_index] = raw_srccat->catsources[src].gauss_ref_stokesI[gauss];
        cropped_src->gauss_ref_stokesQ[gauss_crop_component_index] = raw_srccat->catsources[src].gauss_ref_stokesQ[gauss];
        cropped_src->gauss_ref_stokesU[gauss_crop_component_index] = raw_srccat->catsources[src].gauss_ref_stokesU[gauss];
        cropped_src->gauss_ref_stokesV[gauss_crop_component_index] = raw_srccat->catsources[src].gauss_ref_stokesV[gauss];
        cropped_src->gauss_ref_freqs[gauss_crop_component_index] = raw_srccat->catsources[src].gauss_ref_freqs[gauss];
        cropped_src->gauss_SIs[gauss_crop_component_index] = raw_srccat->catsources[src].gauss_SIs[gauss];

        cropped_src->gauss_majors[gauss_crop_component_index] = raw_srccat->catsources[src].gauss_majors[gauss];
        cropped_src->gauss_minors[gauss_crop_component_index] = raw_srccat->catsources[src].gauss_minors[gauss];
        cropped_src->gauss_pas[gauss_crop_component_index] = raw_srccat->catsources[src].gauss_pas[gauss];

        //Calculate az/za values for each gauss for all time steps
        for (size_t time_step = 0; time_step < num_time_steps; time_step++) {
          convert_radec2azza((double)raw_srccat->catsources[src].gauss_ras[gauss],
                             (double)raw_srccat->catsources[src].gauss_decs[gauss],
                             (double)lsts[time_step], &az, &za);
          cropped_src->gauss_azs[gauss_crop_component_index*num_time_steps + time_step] = (float)az;
          cropped_src->gauss_zas[gauss_crop_component_index*num_time_steps + time_step] = (float)za;
        }

        gauss_crop_component_index ++;

      }//End gauss component loop
      //
      //Loop over shapelet components
      for (size_t shape = 0; shape < raw_srccat->catsources[src].n_shapes; shape++){
        cropped_src->shape_ras[shape_crop_component_index] = raw_srccat->catsources[src].shape_ras[shape];
        cropped_src->shape_decs[shape_crop_component_index] = raw_srccat->catsources[src].shape_decs[shape];
        cropped_src->shape_fluxes[shape_crop_component_index] = raw_srccat->catsources[src].shape_fluxes[shape];
        cropped_src->shape_freqs[shape_crop_component_index] = raw_srccat->catsources[src].shape_freqs[shape];

        cropped_src->shape_majors[shape_crop_component_index] = raw_srccat->catsources[src].shape_majors[shape];
        cropped_src->shape_minors[shape_crop_component_index] = raw_srccat->catsources[src].shape_minors[shape];
        cropped_src->shape_pas[shape_crop_component_index] = raw_srccat->catsources[src].shape_pas[shape];

        //Calculate az/za values for each shapelet component for all time steps
        for (size_t time_step = 0; time_step < num_time_steps; time_step++) {
          convert_radec2azza((double)raw_srccat->catsources[src].shape_ras[shape],
                             (double)raw_srccat->catsources[src].shape_decs[shape],
                             (double)lsts[time_step], &az, &za);
          cropped_src->shape_azs[shape_crop_component_index*num_time_steps + time_step] = (float)az;
          cropped_src->shape_zas[shape_crop_component_index*num_time_steps + time_step] = (float)za;
        }

        if ((int)shape == 0) {
          // Loop over the coefficients for this shapelet source
          // Only do it once as all shapelet component coeffs, n1s, n2s for one
          // source are in 1D arrays in raw_srccat->catsources[src]. As each
          // shapelet component can have any number of coeffs, n1s, n2s, we
          // relate the coeffs, n1s, n2s, to shapelet ra, dec, etc via the
          // cropped_src->shape_param_indexes array. So need to pull that
          // information out from raw_srccat->catsources[src].shape_param_indexes
          // and keep track of how many shapelet coeff components are in
          // the new cropped_src using shape_coeff_component_index
          // Do all this work now as the 1D array goes nicely into a GPU kernel
          for (size_t coeff_ind = 0; coeff_ind < raw_srccat->catsources[src].n_shape_coeffs; coeff_ind++) {
            cropped_src->shape_coeffs[shape_coeff_component_index] = raw_srccat->catsources[src].shape_coeffs[coeff_ind];
            cropped_src->shape_n1s[shape_coeff_component_index] = raw_srccat->catsources[src].shape_n1s[coeff_ind];
            cropped_src->shape_n2s[shape_coeff_component_index] = raw_srccat->catsources[src].shape_n2s[coeff_ind];

            //We do the loop of shapelet coeffs, n1s, n2s only once per source, so we may well get ahead of the
            //shape_crop_component_index. So use to shape_param_index to account for that
            int shape_param_index = raw_srccat->catsources[src].shape_param_indexes[coeff_ind] + shape_crop_component_index;

            cropped_src->shape_param_indexes[shape_coeff_component_index] = shape_param_index;
            shape_coeff_component_index += 1;
          }//end coeff,n1,n2 loop for all shapelet components in this source
        }//end if shape == 0

        //Update shapelet component index for each shapelet component added
        shape_crop_component_index += 1;
      }//End shapelet component loop
    }//End retained sources loop
  }//End if sky_crop_type == CROP_SOURCES

  else if (sky_crop_type == CROP_COMPONENTS) {
    printf("Point components after cropping %d\n",num_point_comp_retained );
    printf("Gaussian components after cropping %d\n",num_gauss_comp_retained );
    printf("Shapelet components after cropping %d\n",num_shape_comp_retained );
    printf("Shapelet coefficients after cropping %d\n",num_shape_coeff_retained );
    //Keep track of components added to cropped_src for indexing
    int point_crop_component_index = 0;
    int gauss_crop_component_index = 0;
    int shape_crop_component_index = 0;
    int shape_coeff_component_index = 0;
    //Loop over all sources in uncropped source catalogue and add all
    //components above the horizon to
    for (size_t src = 0; src < raw_srccat->num_sources; src++){

      //Begin point component loop
      for (size_t point = 0; point < raw_srccat->catsources[src].n_points; point++){
        //Check if point component above horizon
        if (raw_srccat->catsources[src].point_zas[point] < M_PI / 2.0){
          cropped_src->point_ras[point_crop_component_index] = raw_srccat->catsources[src].point_ras[point];
          cropped_src->point_decs[point_crop_component_index] = raw_srccat->catsources[src].point_decs[point];
          cropped_src->point_ref_stokesI[point_crop_component_index] = raw_srccat->catsources[src].point_ref_stokesI[point];
          cropped_src->point_ref_stokesQ[point_crop_component_index] = raw_srccat->catsources[src].point_ref_stokesQ[point];
          cropped_src->point_ref_stokesU[point_crop_component_index] = raw_srccat->catsources[src].point_ref_stokesU[point];
          cropped_src->point_ref_stokesV[point_crop_component_index] = raw_srccat->catsources[src].point_ref_stokesV[point];
          cropped_src->point_ref_freqs[point_crop_component_index] = raw_srccat->catsources[src].point_ref_freqs[point];
          cropped_src->point_SIs[point_crop_component_index] = raw_srccat->catsources[src].point_SIs[point];

          //Calculate az/za values for each point for all time steps
          for (int time_step = 0; time_step < num_time_steps; time_step++) {
            convert_radec2azza((double)raw_srccat->catsources[src].point_ras[point],
                               (double)raw_srccat->catsources[src].point_decs[point],
                               (double)lsts[time_step], &az, &za);
            cropped_src->point_azs[point_crop_component_index*num_time_steps + time_step] = (float)az;
            cropped_src->point_zas[point_crop_component_index*num_time_steps + time_step] = (float)za;
            // printf("AZ ZA %d %.7f %.7f\n",point_crop_component_index*num_time_steps + time_step,(float)az,(float)za );
          }//End az/za calculation loop

          point_crop_component_index ++;
        }//End point above horizon loop
      }//End point component loop


      //Begin gauss component loop
      for (size_t gauss = 0; gauss < raw_srccat->catsources[src].n_gauss; gauss++){
        //Check if gauss component above horizon
        if (raw_srccat->catsources[src].gauss_zas[gauss] < M_PI / 2.0){
          cropped_src->gauss_ras[gauss_crop_component_index] = raw_srccat->catsources[src].gauss_ras[gauss];
          cropped_src->gauss_decs[gauss_crop_component_index] = raw_srccat->catsources[src].gauss_decs[gauss];
          cropped_src->gauss_ref_stokesI[gauss_crop_component_index] = raw_srccat->catsources[src].gauss_ref_stokesI[gauss];
          cropped_src->gauss_ref_stokesQ[gauss_crop_component_index] = raw_srccat->catsources[src].gauss_ref_stokesQ[gauss];
          cropped_src->gauss_ref_stokesU[gauss_crop_component_index] = raw_srccat->catsources[src].gauss_ref_stokesU[gauss];
          cropped_src->gauss_ref_stokesV[gauss_crop_component_index] = raw_srccat->catsources[src].gauss_ref_stokesV[gauss];
          cropped_src->gauss_ref_freqs[gauss_crop_component_index] = raw_srccat->catsources[src].gauss_ref_freqs[gauss];
          cropped_src->gauss_SIs[gauss_crop_component_index] = raw_srccat->catsources[src].gauss_SIs[gauss];

          cropped_src->gauss_majors[gauss_crop_component_index] = raw_srccat->catsources[src].gauss_majors[gauss];
          cropped_src->gauss_minors[gauss_crop_component_index] = raw_srccat->catsources[src].gauss_minors[gauss];
          cropped_src->gauss_pas[gauss_crop_component_index] = raw_srccat->catsources[src].gauss_pas[gauss];

          //Calculate az/za values for each gauss for all time steps
          for (size_t time_step = 0; time_step < num_time_steps; time_step++) {
            convert_radec2azza((double)raw_srccat->catsources[src].gauss_ras[gauss],
                               (double)raw_srccat->catsources[src].gauss_decs[gauss],
                               (double)lsts[time_step], &az, &za);
            cropped_src->gauss_azs[gauss_crop_component_index*num_time_steps + time_step] = (float)az;
            cropped_src->gauss_zas[gauss_crop_component_index*num_time_steps + time_step] = (float)za;
          }//End az/za calculation loop

          gauss_crop_component_index ++;
        }//End gauss above horizon loop
      }//End gauss component loop

      //Loop over shapelet components
      for (size_t shape = 0; shape < raw_srccat->catsources[src].n_shapes; shape++){

        if (raw_srccat->catsources[src].shape_zas[shape] < M_PI / 2.0){
          cropped_src->shape_ras[shape_crop_component_index] = raw_srccat->catsources[src].shape_ras[shape];
          cropped_src->shape_decs[shape_crop_component_index] = raw_srccat->catsources[src].shape_decs[shape];
          cropped_src->shape_fluxes[shape_crop_component_index] = raw_srccat->catsources[src].shape_fluxes[shape];
          cropped_src->shape_freqs[shape_crop_component_index] = raw_srccat->catsources[src].shape_freqs[shape];

          cropped_src->shape_majors[shape_crop_component_index] = raw_srccat->catsources[src].shape_majors[shape];
          cropped_src->shape_minors[shape_crop_component_index] = raw_srccat->catsources[src].shape_minors[shape];
          cropped_src->shape_pas[shape_crop_component_index] = raw_srccat->catsources[src].shape_pas[shape];

          //Calculate az/za values for each shape for all time steps
          for (size_t time_step = 0; time_step < num_time_steps; time_step++) {
            convert_radec2azza((double)raw_srccat->catsources[src].shape_ras[shape],
                               (double)raw_srccat->catsources[src].shape_decs[shape],
                               (double)lsts[time_step], &az, &za);
            cropped_src->shape_azs[shape_crop_component_index*num_time_steps + time_step] = (float)az;
            cropped_src->shape_zas[shape_crop_component_index*num_time_steps + time_step] = (float)za;
          }//End az/za calculation loop

          //Loop through all shapelet coeffs,n1s,n2s for this source; these 1D arrays
          //can contain information from multiple shapelet components.
          //raw_srccat->catsources[src].shape_param_indexes is used to match
          //the coeff,n1,n2 to each component, so check those indexes and grab
          //information if correct.
          // Do all this work now as the 1D array goes nicely into a GPU kernel
          for (size_t coeff_ind = 0; coeff_ind < raw_srccat->catsources[src].n_shape_coeffs; coeff_ind++) {
            //Check if we are on the cofrect component
            if ( (int)raw_srccat->catsources[src].shape_param_indexes[coeff_ind] == shape ){
              cropped_src->shape_coeffs[shape_coeff_component_index] = raw_srccat->catsources[src].shape_coeffs[coeff_ind];
              cropped_src->shape_n1s[shape_coeff_component_index] = raw_srccat->catsources[src].shape_n1s[coeff_ind];
              cropped_src->shape_n2s[shape_coeff_component_index] = raw_srccat->catsources[src].shape_n2s[coeff_ind];
              cropped_src->shape_param_indexes[shape_coeff_component_index] = shape_crop_component_index;

              shape_coeff_component_index += 1;

            }//end test if correct shape coeff
          }// end shapelet coeff loop
          shape_crop_component_index += 1;
        }//End shapelet above horizon loop
      }//End shapelet component loop
    }//End raw_srccat source loop
  }//End if sky_crop_type == CROP_COMPONENTS

  // for (size_t comp = 0; comp < cropped_src->n_shapes; comp++) {
  //
  //   printf("%f %f %f %f %f %f %f %f %f %f %f\n", cropped_src->shape_ras[comp], cropped_src->shape_decs[comp],
  //         cropped_src->shape_fluxes[comp],cropped_src->shape_freqs[comp],
  //         cropped_src->shape_majors[comp],cropped_src->shape_minors[comp],
  //         cropped_src->shape_pas[comp],cropped_src->shape_n1s[comp],
  //         cropped_src->shape_n2s[comp],cropped_src->shape_coeffs[comp],
  //         cropped_src->shape_param_indexes[comp]);
  // }

  return cropped_src;

}
