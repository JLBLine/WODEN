/*******************************************************************************
*  Methods to read in and crop sky models to everything above the horizon.
*  @author J.L.B. Line
*
*  Please see documentation in ../include/create_sky_model.h or online at
*  https://woden.readthedocs.io/en/latest/index.html
*******************************************************************************/
#include <stdio.h>
#include <string.h>
#include <erfa.h>
#include <math.h>
#include "constants.h"
#include "create_sky_model.h"
#include "woden_precision_defs.h"

//Copied straight up from:
//https://stackoverflow.com/questions/744766/how-to-compare-ends-of-strings-in-c
//Checks whether a sting ends if a specific suffix
int EndsWith(const char *str, const char *suffix)
{
    if (!str || !suffix)
        return 0;
    size_t lenstr = strlen(str);
    size_t lensuffix = strlen(suffix);
    if (lensuffix >  lenstr)
        return 0;
    return strncmp(str + lenstr - lensuffix, suffix, lensuffix) == 0;
}

//Ends with .txt or .yaml
int EndsWithTxt(const char *str) { return EndsWith(str, ".txt"); }
int EndsWithYaml(const char *str) { return EndsWith(str, ".yaml"); }

//Read in either a WODEN style text file or a hyperdrive style yaml file
int read_skymodel(const char *srclist, source_catalogue_t *raw_srccat){

  int status = 0;

  int is_text = EndsWithTxt(srclist);
  int is_yaml = EndsWithYaml(srclist);

  if (is_text == 1) {
    status = read_text_skymodel(srclist, raw_srccat);
  }
  else if (is_yaml == 1) {
    status = read_yaml_skymodel(srclist, raw_srccat);
  }
  else {
    printf("The input sky model does not end in either '.txt' or '.yaml'. Unable to guess format\n");
    return 1;
  }

  return status;

}


/*********************************
// Takes an ra, dec, lst (all rad) and returns the azimuth and zenith angle
// assuming the latitude of the MWA. Uses the ERFA library to do the
// transformation using the ha.
**********************************/
void convert_radec2azza(double ra, double dec, double lst, double latitude,
     double * az, double * za){

  double erfa_az, el;
  double ha = lst - ra;

  eraHd2ae( ha, dec, latitude, &erfa_az, &el );

  * az = erfa_az;
  * za = M_PI / 2. - el;

}

/******************************************************************************

*******************************************************************************/
void horizon_test(e_sky_crop sky_crop_type, components_t * components,
     int num_comps, double *lsts, double latitude,
     e_horizon * all_comps_above_horizon){

  double az, za;

  for (int comp_ind = 0; comp_ind < num_comps; comp_ind++) {
    //Calculate az/za for all point components
    convert_radec2azza(components->ras[comp_ind],
                       components->decs[comp_ind],
                       lsts[0], latitude, &az, &za);

    components->azs[comp_ind] = (user_precision_t)az;
    components->zas[comp_ind] = (user_precision_t)za;


    //Check if component is above horizon, and flag the source if not
    if (sky_crop_type == CROP_SOURCES) {
      if (za >= M_PI / 2.0){
        //When cropping sources, if one component is below the horizon,
        //we're dumping the whole source
         * all_comps_above_horizon = BELOW;
      }
    }//end if CROP_SOURCES
  }//end loop over components
}//end function

void realloc_and_add_component(e_sky_crop croptype,
                               e_component_type comptype,
                               e_flux_type fluxtype,
                               source_t *cropped_src,
                               source_t *original_src,
                               int orig_flux_ind,
                               double *lsts, double latitude, int num_time_steps){

  components_t *components_original = NULL;
  components_t *components_cropped = NULL;
  if (comptype == POINT) {
    components_original = &original_src->point_components;
    components_cropped = &cropped_src->point_components;
  } else if (comptype == GAUSSIAN) {
    components_original = &original_src->gauss_components;
    components_cropped = &cropped_src->gauss_components;
  } else if (comptype == SHAPELET) {
    components_original = &original_src->shape_components;
    components_cropped = &cropped_src->shape_components;
  }


  int orig_comp_ind;

  if (fluxtype == POWER_LAW) {
    orig_comp_ind = components_original->power_comp_inds[orig_flux_ind];
  } else if (fluxtype == CURVED_POWER_LAW) {
    orig_comp_ind = components_original->curve_comp_inds[orig_flux_ind];
  } else if (fluxtype == LIST) {
    orig_comp_ind = components_original->list_comp_inds[orig_flux_ind];
  } else {
    printf("`realloc_and_add_component`: You've entered an unrecognised fluxtype, this is bad\n");
    orig_comp_ind = -1;
  }

  //If == 1, we should be adding this COMPONENT. If we're cropping by SOURCE,
  //we should only ever call this function when we want to add the whole SOURCE
  //so default this to 1
  int add_component = 1;

  //If cropping by COMPONENT, check if this COMPONENT is below the horizon
  //and set add_component = 0 if so, which means we skip adding the component
  if (croptype == CROP_COMPONENTS){
    if (components_original->zas[orig_comp_ind] >= M_PI / 2.0) {
      add_component = 0;
    }
  }

  if (add_component == 1) {

    int crop_comp_ind;
    int crop_power_ind;
    int crop_curve_ind;
    int crop_list_ind;

    if (comptype == POINT) {
      crop_comp_ind = cropped_src->n_points;
      crop_power_ind = cropped_src->n_point_powers;
      crop_curve_ind = cropped_src->n_point_curves;
      crop_list_ind = cropped_src->n_point_lists;
    }
    else if (comptype == GAUSSIAN) {
      crop_comp_ind = cropped_src->n_gauss;
      crop_power_ind = cropped_src->n_gauss_powers;
      crop_curve_ind = cropped_src->n_gauss_curves;
      crop_list_ind = cropped_src->n_gauss_lists;
    }
    else if (comptype == SHAPELET) {
      crop_comp_ind = cropped_src->n_shapes;
      crop_power_ind = cropped_src->n_shape_powers;
      crop_curve_ind = cropped_src->n_shape_curves;
      crop_list_ind = cropped_src->n_shape_lists;
    } else {
      printf("`realloc_and_add_component`:You've entered an unrecognised comptype, this is bad\n");
      crop_comp_ind = -1;
      crop_power_ind = -1;
      crop_curve_ind = -1;
      crop_list_ind = -1;
    }

    components_cropped->ras = realloc(components_cropped->ras,
                                              sizeof(double)*(crop_comp_ind + 1));
    components_cropped->ras[crop_comp_ind] = components_original->ras[orig_comp_ind];

    components_cropped->decs = realloc(components_cropped->decs,
                                              sizeof(double)*(crop_comp_ind + 1));
    components_cropped->decs[crop_comp_ind] = components_original->decs[orig_comp_ind];

    double az, za;
    //Calculate az/za values for  all time steps
    components_cropped->azs = realloc(components_cropped->azs,
                      sizeof(user_precision_t)*(crop_comp_ind + 1)*num_time_steps);
    components_cropped->zas = realloc(components_cropped->zas,
                      sizeof(user_precision_t)*(crop_comp_ind + 1)*num_time_steps);


    for (int time_step = 0; time_step < num_time_steps; time_step++) {
      convert_radec2azza(components_cropped->ras[crop_comp_ind],
                         components_cropped->decs[crop_comp_ind],
                         lsts[time_step], latitude, &az, &za);
      components_cropped->azs[crop_comp_ind*num_time_steps + time_step] = (user_precision_t)az;
      components_cropped->zas[crop_comp_ind*num_time_steps + time_step] = (user_precision_t)za;
    }


    if (comptype == GAUSSIAN || comptype == SHAPELET) {
      components_cropped->majors = realloc(components_cropped->majors,
                                                sizeof(user_precision_t)*(crop_comp_ind + 1));
      components_cropped->majors[crop_comp_ind] = components_original->majors[orig_comp_ind];

      components_cropped->minors = realloc(components_cropped->minors,
                                                sizeof(user_precision_t)*(crop_comp_ind + 1));
      components_cropped->minors[crop_comp_ind] = components_original->minors[orig_comp_ind];

      components_cropped->pas = realloc(components_cropped->pas,
                                                sizeof(user_precision_t)*(crop_comp_ind + 1));
      components_cropped->pas[crop_comp_ind] = components_original->pas[orig_comp_ind];

    }

    if (comptype == SHAPELET) {

      // printf("Start of new coeff loop %d================\n",orig_comp_ind );

      int cropped_coeff_ind = cropped_src->n_shape_coeffs;

      // Loop through all shapelet coeffs,n1s,n2s for this source; these 1D arrays
      //can contain information from multiple shapelet components.
      //raw_srccat->sources[src].shape_components.param_indexes is used to match
      //the coeff,n1,n2 to each component, so check those indexes and grab
      //information if correct.
      // Do all this work now as the 1D array goes nicely into a GPU kernel
      for (int orig_coeff_ind = 0; orig_coeff_ind < original_src->n_shape_coeffs; orig_coeff_ind++) {
        //Check if we are on the cofrect component
        if ( (int)components_original->param_indexes[orig_coeff_ind] == orig_comp_ind ){

          // printf("orig_coeff_ind, orig_comp_ind %d %d\n", orig_coeff_ind, orig_comp_ind );


          components_cropped->shape_coeffs = realloc(components_cropped->shape_coeffs,
                              (cropped_coeff_ind + 1)*sizeof(user_precision_t));
          components_cropped->n1s = realloc(components_cropped->n1s,
                              (cropped_coeff_ind + 1)*sizeof(user_precision_t));
          components_cropped->n2s = realloc(components_cropped->n2s,
                              (cropped_coeff_ind + 1)*sizeof(user_precision_t));
          components_cropped->param_indexes = realloc(components_cropped->param_indexes,
                              (cropped_coeff_ind + 1)*sizeof(user_precision_t));

          components_cropped->shape_coeffs[cropped_coeff_ind] = components_original->shape_coeffs[orig_coeff_ind];
          components_cropped->n1s[cropped_coeff_ind] = components_original->n1s[orig_coeff_ind];
          components_cropped->n2s[cropped_coeff_ind] = components_original->n2s[orig_coeff_ind];
          components_cropped->param_indexes[cropped_coeff_ind] = crop_comp_ind;

          cropped_coeff_ind += 1;

        }//end test if correct shape coeff
      }// end shapelet coeff loop

      cropped_src->n_shape_coeffs = cropped_coeff_ind;

    }//end if (comptype == SHAPELET)




    if (fluxtype == POWER_LAW) {
      components_cropped->power_comp_inds = realloc(components_cropped->power_comp_inds,
                                                 (crop_power_ind + 1)*sizeof(int));

      components_cropped->power_ref_freqs = realloc(components_cropped->power_ref_freqs,
                                         (crop_power_ind + 1)*sizeof(double));
      components_cropped->power_ref_stokesI = realloc(components_cropped->power_ref_stokesI,
                                    (crop_power_ind + 1)*sizeof(user_precision_t));
      components_cropped->power_ref_stokesQ = realloc(components_cropped->power_ref_stokesQ,
                                    (crop_power_ind + 1)*sizeof(user_precision_t));
      components_cropped->power_ref_stokesU = realloc(components_cropped->power_ref_stokesU,
                                    (crop_power_ind + 1)*sizeof(user_precision_t));
      components_cropped->power_ref_stokesV = realloc(components_cropped->power_ref_stokesV,
                                    (crop_power_ind + 1)*sizeof(user_precision_t));
      components_cropped->power_SIs = realloc(components_cropped->power_SIs,
                                    (crop_power_ind + 1)*sizeof(user_precision_t));

      components_cropped->power_comp_inds[crop_power_ind] = crop_comp_ind;
      components_cropped->power_ref_freqs[crop_power_ind] = components_original->power_ref_freqs[orig_flux_ind];
      components_cropped->power_ref_stokesI[crop_power_ind] = components_original->power_ref_stokesI[orig_flux_ind];
      components_cropped->power_ref_stokesQ[crop_power_ind] = components_original->power_ref_stokesQ[orig_flux_ind];
      components_cropped->power_ref_stokesU[crop_power_ind] = components_original->power_ref_stokesU[orig_flux_ind];
      components_cropped->power_ref_stokesV[crop_power_ind] = components_original->power_ref_stokesV[orig_flux_ind];
      components_cropped->power_SIs[crop_power_ind] = components_original->power_SIs[orig_flux_ind];

    }
    else if (fluxtype == CURVED_POWER_LAW) {
      components_cropped->curve_comp_inds = realloc(components_cropped->curve_comp_inds,
                                                 (crop_curve_ind + 1)*sizeof(int));

      components_cropped->curve_ref_freqs = realloc(components_cropped->curve_ref_freqs,
                                         (crop_curve_ind + 1)*sizeof(double));
      components_cropped->curve_ref_stokesI = realloc(components_cropped->curve_ref_stokesI,
                                    (crop_curve_ind + 1)*sizeof(user_precision_t));
      components_cropped->curve_ref_stokesQ = realloc(components_cropped->curve_ref_stokesQ,
                                    (crop_curve_ind + 1)*sizeof(user_precision_t));
      components_cropped->curve_ref_stokesU = realloc(components_cropped->curve_ref_stokesU,
                                    (crop_curve_ind + 1)*sizeof(user_precision_t));
      components_cropped->curve_ref_stokesV = realloc(components_cropped->curve_ref_stokesV,
                                    (crop_curve_ind + 1)*sizeof(user_precision_t));
      components_cropped->curve_SIs = realloc(components_cropped->curve_SIs,
                                    (crop_curve_ind + 1)*sizeof(user_precision_t));
      components_cropped->curve_qs = realloc(components_cropped->curve_qs,
                                    (crop_curve_ind + 1)*sizeof(user_precision_t));

      components_cropped->curve_comp_inds[crop_curve_ind] = crop_comp_ind;
      components_cropped->curve_ref_freqs[crop_curve_ind] = components_original->curve_ref_freqs[orig_flux_ind];
      components_cropped->curve_ref_stokesI[crop_curve_ind] = components_original->curve_ref_stokesI[orig_flux_ind];
      components_cropped->curve_ref_stokesQ[crop_curve_ind] = components_original->curve_ref_stokesQ[orig_flux_ind];
      components_cropped->curve_ref_stokesU[crop_curve_ind] = components_original->curve_ref_stokesU[orig_flux_ind];
      components_cropped->curve_ref_stokesV[crop_curve_ind] = components_original->curve_ref_stokesV[orig_flux_ind];
      components_cropped->curve_SIs[crop_curve_ind] = components_original->curve_SIs[orig_flux_ind];
      components_cropped->curve_qs[crop_curve_ind] = components_original->curve_qs[orig_flux_ind];
    }
    else if (fluxtype == LIST) {
      components_cropped->list_comp_inds = realloc(components_cropped->list_comp_inds,
                                                 (crop_list_ind + 1)*sizeof(int));

      //How much we have to realloc is going to depend on this number

      int extra_list_entries = components_original->num_list_values[orig_flux_ind];

      int current_list_entries;
      if (crop_list_ind == 0) {
        current_list_entries = 0;
      } else {
        current_list_entries = components_cropped->num_list_values[crop_list_ind - 1] + components_cropped->list_start_indexes[crop_list_ind - 1];
      }

      // printf("RIGHT HERE orig_flux_ind, crop_list_ind, current_list_entries %d %d %d\n",orig_flux_ind, crop_list_ind, current_list_entries);

      components_cropped->list_comp_inds = realloc(components_cropped->list_comp_inds,
                                                 (crop_list_ind + 1)*sizeof(int));
      components_cropped->num_list_values = realloc(components_cropped->num_list_values,
                                                 (crop_list_ind + 1)*sizeof(int));
      components_cropped->list_start_indexes = realloc(components_cropped->list_start_indexes,
                                                 (crop_list_ind + 1)*sizeof(int));

      components_cropped->list_comp_inds[crop_list_ind] = crop_comp_ind;
      components_cropped->num_list_values[crop_list_ind] = extra_list_entries;
      components_cropped->list_start_indexes[crop_list_ind] = current_list_entries;


      int realloc_size = extra_list_entries + current_list_entries;

      components_cropped->list_freqs = realloc(components_cropped->list_freqs,
                                               realloc_size*sizeof(double));
      components_cropped->list_stokesI = realloc(components_cropped->list_stokesI,
                                         realloc_size*sizeof(user_precision_t));
      components_cropped->list_stokesQ = realloc(components_cropped->list_stokesQ,
                                         realloc_size*sizeof(user_precision_t));
      components_cropped->list_stokesU = realloc(components_cropped->list_stokesU,
                                         realloc_size*sizeof(user_precision_t));
      components_cropped->list_stokesV = realloc(components_cropped->list_stokesV,
                                         realloc_size*sizeof(user_precision_t));

      int orig_start_ind = components_original->list_start_indexes[orig_flux_ind];

      for (int list_increment = 0; list_increment < extra_list_entries; list_increment++) {

        // printf("YUP %d %d %d %.3e\n",current_list_entries, list_increment, crop_list_ind,
        //      components_original->list_freqs[orig_start_ind + list_increment] );

        components_cropped->list_freqs[current_list_entries + list_increment] = components_original->list_freqs[orig_start_ind + list_increment];
        components_cropped->list_stokesI[current_list_entries + list_increment] = components_original->list_stokesI[orig_start_ind + list_increment];
        components_cropped->list_stokesQ[current_list_entries + list_increment] = components_original->list_stokesQ[orig_start_ind + list_increment];
        components_cropped->list_stokesU[current_list_entries + list_increment] = components_original->list_stokesU[orig_start_ind + list_increment];
        components_cropped->list_stokesV[current_list_entries + list_increment] = components_original->list_stokesV[orig_start_ind + list_increment];
      }

    }// end if (fluxtype == LIST)

    //Increment the counters as appropriate
    cropped_src->n_comps += 1;

    if (comptype == POINT) {
      cropped_src->n_points += 1;
      if (fluxtype == POWER_LAW) {
        cropped_src->n_point_powers += 1;
      } else if (fluxtype == CURVED_POWER_LAW) {
        cropped_src->n_point_curves += 1;
      } else if (fluxtype == LIST) {
        cropped_src->n_point_lists += 1;
      }
    }
    else if (comptype == GAUSSIAN) {
      cropped_src->n_gauss += 1;
      if (fluxtype == POWER_LAW) {
        cropped_src->n_gauss_powers += 1;
      } else if (fluxtype == CURVED_POWER_LAW) {
        cropped_src->n_gauss_curves += 1;
      } else if (fluxtype == LIST) {
        cropped_src->n_gauss_lists += 1;
      }
    }
    else if (comptype == SHAPELET) {
      cropped_src->n_shapes += 1;
      if (fluxtype == POWER_LAW) {
        cropped_src->n_shape_powers += 1;
      } else if (fluxtype == CURVED_POWER_LAW) {
        cropped_src->n_shape_curves += 1;
      } else if (fluxtype == LIST) {
        cropped_src->n_shape_lists += 1;
      }
    }
  }//end adding component

}

// if CROP_SOURCES
void add_source_to_cropped_source(e_sky_crop croptype,
                                  source_t *cropped_src, source_t *original_src,
                                  double *lsts, double latitude, int num_time_steps){
  //For everything in the selected original source, add to the cropped source

  // int orig_comp_ind;

  for (int orig_power_ind = 0; orig_power_ind < original_src->n_point_powers; orig_power_ind++) {

    realloc_and_add_component(croptype, POINT, POWER_LAW,
                              cropped_src, original_src,
                              orig_power_ind,
                              lsts, latitude, num_time_steps);
  } //end loop over POINT + POWER_LAW

  for (int orig_curve_ind = 0; orig_curve_ind < original_src->n_point_curves; orig_curve_ind++) {

    realloc_and_add_component(croptype, POINT, CURVED_POWER_LAW,
                              cropped_src, original_src,
                              orig_curve_ind,
                              lsts, latitude, num_time_steps);
  } //end loop over POINT + CURVED_POWER_LAW

  for (int orig_list_ind = 0; orig_list_ind < original_src->n_point_lists; orig_list_ind++) {

    realloc_and_add_component(croptype, POINT, LIST,
                              cropped_src, original_src,
                              orig_list_ind,
                              lsts, latitude, num_time_steps);
  } //end loop over POINT + LIST

  for (int orig_power_ind = 0; orig_power_ind < original_src->n_gauss_powers; orig_power_ind++) {

    realloc_and_add_component(croptype, GAUSSIAN, POWER_LAW,
                              cropped_src, original_src,
                              orig_power_ind,
                              lsts, latitude, num_time_steps);
  } //end loop over GAUSSIAN + POWER_LAW

  for (int orig_curve_ind = 0; orig_curve_ind < original_src->n_gauss_curves; orig_curve_ind++) {

    realloc_and_add_component(croptype, GAUSSIAN, CURVED_POWER_LAW,
                              cropped_src, original_src,
                              orig_curve_ind,
                              lsts, latitude, num_time_steps);
  } //end loop over GAUSSIAN + CURVED_POWER_LAW

  for (int orig_list_ind = 0; orig_list_ind < original_src->n_gauss_lists; orig_list_ind++) {

    realloc_and_add_component(croptype, GAUSSIAN, LIST,
                              cropped_src, original_src,
                              orig_list_ind,
                              lsts, latitude, num_time_steps);
  } //end loop over GAUSSIAN + LIST

  for (int orig_power_ind = 0; orig_power_ind < original_src->n_shape_powers; orig_power_ind++) {

    realloc_and_add_component(croptype, SHAPELET, POWER_LAW,
                              cropped_src, original_src,
                              orig_power_ind,
                              lsts, latitude, num_time_steps);
  } //end loop over SHAPELET + POWER_LAW

  for (int orig_curve_ind = 0; orig_curve_ind < original_src->n_shape_curves; orig_curve_ind++) {

    realloc_and_add_component(croptype, SHAPELET, CURVED_POWER_LAW,
                              cropped_src, original_src,
                              orig_curve_ind,
                              lsts, latitude, num_time_steps);
  } //end loop over SHAPELET + CURVED_POWER_LAW

  for (int orig_list_ind = 0; orig_list_ind < original_src->n_shape_lists; orig_list_ind++) {

    realloc_and_add_component(croptype, SHAPELET, LIST,
                              cropped_src, original_src,
                              orig_list_ind,
                              lsts, latitude, num_time_steps);
  } //end loop over SHAPELET + LIST


} //end of add_source_to_cropped_source function

void null_component_sky_model_arrays(components_t * components) {

  //Commonn variables===========================================================
  // printf("WE BE NULLING\n");
  components->ras = NULL;
  components->decs =  NULL;
  components->azs =  NULL;
  components->zas =  NULL;

  //Power law stuff
  components->power_ref_stokesI = NULL;
  components->power_ref_stokesQ = NULL;
  components->power_ref_stokesU = NULL;
  components->power_ref_stokesV = NULL;
  components->power_ref_freqs = NULL;
  components->power_SIs = NULL;
  components->power_comp_inds = NULL;

  //Curved power law stuff
  components->curve_ref_stokesI = NULL;
  components->curve_ref_stokesQ = NULL;
  components->curve_ref_stokesU = NULL;
  components->curve_ref_stokesV = NULL;
  components->curve_ref_freqs = NULL;
  components->curve_SIs = NULL;
  components->curve_qs = NULL;
  components->curve_comp_inds = NULL;

  //List flux things
  components->list_freqs = NULL;
  components->list_stokesI = NULL;
  components->list_stokesQ = NULL;
  components->list_stokesU = NULL;
  components->list_stokesV = NULL;
  components->num_list_values = NULL;
  components->list_start_indexes = NULL;
  components->list_comp_inds = NULL;

  //Component type specific stuff but doesn't hurt to NULL them even if they
  //aren't being used
  components->majors = NULL;
  components->minors = NULL;
  components->pas = NULL;
  components->shape_coeffs = NULL;
  components->n1s = NULL;
  components->n2s = NULL;
  components->param_indexes = NULL;

  components->azs = NULL;
  components->zas = NULL;
  components->beam_has = NULL;
  components->beam_decs = NULL;

}

//Initialise all the various component/flux counters as zero
//Will be reallocing a number of arrays so start them off as NULL
void source_zero_counters_and_null_components(source_t *source){

  //Set all the numbers to zero so we can increment from there
  source->n_comps = 0;
  source->n_points = 0;
  source->n_point_lists = 0;
  source->n_point_powers = 0;
  source->n_point_curves = 0;
  source->n_gauss = 0;
  source->n_gauss_lists = 0;
  source->n_gauss_powers = 0;
  source->n_gauss_curves = 0;
  source->n_shapes = 0;
  source->n_shape_lists = 0;
  source->n_shape_powers = 0;
  source->n_shape_curves = 0;
  source->n_shape_coeffs = 0;

  //NULL a number of pointers so we can realloc them
  null_component_sky_model_arrays(&source->point_components);
  null_component_sky_model_arrays(&source->gauss_components);
  null_component_sky_model_arrays(&source->shape_components);
}


/*********************************
// Crop a sky model contained in a source_catalogue_t struct (raw_srccat) and
// crop out all sources below the horizon (at the beginning of the observation).
// Return a new single source_t (cropped_srccat) that contains the full
// cropped sky model
// First part of the function calculates az/za for the initial time step, and
// counts how many components / sources are to be saved
// Second part mallocs a big enough source_t struct to contain all cropped
// components, and copies all relevant data across from raw_srccat

// Possible TODO update the sky model for every time step to account for sources that
// have risen and set during the observation?

**********************************/
source_t * crop_sky_model(source_catalogue_t *raw_srccat, double *lsts,
              double latitude, int num_time_steps, e_sky_crop sky_crop_type){

  for (int src = 0; src < raw_srccat->num_sources; src++){

    //Initialise containers for azimuth and zenith angles, used in cropping
    //later on

    raw_srccat->sources[src].point_components.azs = malloc( raw_srccat->sources[src].n_points * sizeof(user_precision_t) );
    raw_srccat->sources[src].point_components.zas = malloc( raw_srccat->sources[src].n_points * sizeof(user_precision_t) );

    raw_srccat->sources[src].gauss_components.azs = malloc( raw_srccat->sources[src].n_gauss * sizeof(user_precision_t) );
    raw_srccat->sources[src].gauss_components.zas = malloc( raw_srccat->sources[src].n_gauss * sizeof(user_precision_t) );

    raw_srccat->sources[src].shape_components.azs = malloc( raw_srccat->sources[src].n_shapes * sizeof(user_precision_t) );
    raw_srccat->sources[src].shape_components.zas = malloc( raw_srccat->sources[src].n_shapes * sizeof(user_precision_t) );

  }

  int num_sources_retained=0;
  // int *cropped_src_indexes=NULL;
  e_horizon all_comps_above_horizon;

  // //Make an empty source_t and malloc using the numbers gather above
  source_t *cropped_src=NULL;
  cropped_src = malloc(sizeof(source_t));

  source_zero_counters_and_null_components(cropped_src);

  //Go through all components, calculate the az/za, and shove into
  //the new cropped source if aove horizon
  for (int src = 0; src < raw_srccat->num_sources; src++){
    all_comps_above_horizon = ABOVE;

    //Loop over POINT components and do horizon check
    horizon_test(sky_crop_type, &raw_srccat->sources[src].point_components,
                    raw_srccat->sources[src].n_points, lsts, latitude,
                    &all_comps_above_horizon);

    //Loop over GAUSSIAN components and do horizon check
    horizon_test(sky_crop_type, &raw_srccat->sources[src].gauss_components,
                    raw_srccat->sources[src].n_gauss, lsts, latitude,
                    &all_comps_above_horizon);

    //Loop over SHAPELET components and do horizon check
    horizon_test(sky_crop_type, &raw_srccat->sources[src].shape_components,
                    raw_srccat->sources[src].n_shapes, lsts, latitude,
                    &all_comps_above_horizon);

    //After checking all components in a source, if cropping out sources,
    //check all components were above horizon. If so, update component type
    //counters, and mark down the index of the sources that survied in
    //cropped_src_indexes

    if (sky_crop_type == CROP_SOURCES) {
      if (all_comps_above_horizon == ABOVE) {
        num_sources_retained ++;

        add_source_to_cropped_source(sky_crop_type,
                                     cropped_src, &raw_srccat->sources[src],
                                     lsts, latitude, num_time_steps);

      }//end if all_comps_above_horizon == ABOVE
    }//end if sky_crop_type == CROP_SOURCES

    else if (sky_crop_type == CROP_COMPONENTS) {
      add_source_to_cropped_source(sky_crop_type,
                                   cropped_src, &raw_srccat->sources[src],
                                   lsts, latitude, num_time_steps);
    }

  }//Finish checking az/za loop  and adding components to cropped src

  printf("The full cropped catalogue contains: \n");

  // printf("\tcropped_src->n_comps %d\n", cropped_src->n_comps);
  printf("\tPoint components: %d\n", cropped_src->n_points);
  printf("\t\t(Power laws %d)\n", cropped_src->n_point_powers);
  printf("\t\t(Curved power laws %d)\n", cropped_src->n_point_curves);
  printf("\t\t(List fluxes %d)\n", cropped_src->n_point_lists);

  printf("\tGaussian components: %d\n", cropped_src->n_gauss);
  printf("\t\t(Power laws %d)\n", cropped_src->n_gauss_powers);
  printf("\t\t(Curved power laws %d)\n", cropped_src->n_gauss_curves);
  printf("\t\t(List fluxes %d)\n",  cropped_src->n_gauss_lists);

  printf("\tShapelet components: %d\n", cropped_src->n_shapes);
  printf("\tShapelet basis functions: %d\n", cropped_src->n_shape_coeffs);
  printf("\t\t(Power laws %d)\n", cropped_src->n_shape_powers);
  printf("\t\t(Curved power laws %d)\n", cropped_src->n_shape_curves);
  printf("\t\t(List fluxes %d)\n",  cropped_src->n_shape_lists);


  return cropped_src;
}


// void free_components(components_t *components, int num_lists, int num_powers,
//                      int num_curves, int num_shapes){
//   printf("FREEING RA<DEC\n");
//   //Instrinsic to COMPONENT values
//   free(components->ras); /*!< COMPONENT right ascensions (radians) */
//   free(components->decs); /*!< COMPONENT declinations (radians) */
//
//   printf("FREED RA<DEC\n");
//
//   if (num_powers > 0){
//     free(components->power_ref_freqs);
//     free(components->power_ref_stokesI);
//     free(components->power_ref_stokesQ);
//     free(components->power_ref_stokesU);
//     free(components->power_ref_stokesV);
//     free(components->power_SIs);
//     free(components->power_comp_inds);
//   }
//
//   if (num_curves > 0){
//     free(components->curve_ref_freqs);
//     free(components->curve_ref_stokesI);
//     free(components->curve_ref_stokesQ);
//     free(components->curve_ref_stokesU);
//     free(components->curve_ref_stokesV);
//     free(components->curve_SIs);
//     free(components->curve_qs);
//     free(components->curve_comp_inds);
//   }
//
//   if (num_lists > 0){
//     free(components->list_comp_inds);
//     free(components->list_freqs);
//     free(components->list_stokesI);
//     free(components->list_stokesQ);
//     free(components->list_stokesU);
//     free(components->list_stokesV);
//     free(components->num_list_values);
//     free(components->list_start_indexes);
//   }
//
//   // //something to store extrapolated output fluxes in
//   // user_precision_t **extrap_stokesI; /*!< extrapolated COMPONENT Stokes I flux densities (Jy) */
//   // user_precision_t **extrap_stokesQ; /*!< extrapolated COMPONENT Stokes I flux densities (Jy) */
//   // user_precision_t **extrap_stokesU; /*!< extrapolated COMPONENT Stokes I flux densities (Jy) */
//   // user_precision_t **extrap_stokesV; /*!< extrapolated COMPONENT Stokes I flux densities (Jy) */
//
//   if (num_shapes > 0) {
//     free(components->shape_coeffs);
//     free(components->n1s);
//     free(components->n2s);
//     free(components->majors);
//     free(components->minors);
//     free(components->pas);
//     free(components->param_indexes);
//   }
//   // //Specific to observation settings for these COMPONENTs
//   // user_precision_t *azs; /*!< SHAPELET source azimuth angles for all time steps */
//   // user_precision_t *zas; /*!< SHAPELET source zenith angles for all time steps */
//   // double *beam_has; /*!< Hour angle of COMPONENTs for all time steps, used for
//   //  beam calculations */
//   // double *beam_decs; /*!< Declinations of COMPONENTs for all time steps, used for
//   //  beam calculations */
//
//   // free(components);
//
// }
//
//
// void free_source(source_t *source){
//   if (source->n_points > 0){
//     free_components(&source->point_components, source->n_point_lists,
//                          source->n_point_powers, source->n_point_curves, 0);
//   }
//
//   if (source->n_gauss > 0){
//     free_components(&source->gauss_components, source->n_gauss_lists,
//                          source->n_gauss_powers, source->n_gauss_curves, 0);
//   }
//
//   if (source->n_shapes > 0){
//     free_components(&source->shape_components, source->n_shape_lists,
//                          source->n_shape_powers, source->n_shape_curves,
//                          source->n_shapes);
//   }
//
//   free(source);
//
// }
//
// void free_source_catalogue(source_catalogue_t *source_catalogue){
//   for (int i = 0; i < source_catalogue->num_sources; i++) {
//     free_source(&source_catalogue->sources[i]);
//   }
//   free(source_catalogue);
// }
