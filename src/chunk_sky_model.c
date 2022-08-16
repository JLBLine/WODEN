#include <stdio.h>
#include <math.h>
#include "chunk_sky_model.h"


//Switch through different POINT and GAUSSIAN chunking states - if P, only
//POINTs in a chunk, PG both POINT and GAUSSIAN, G just GAUSSIAN
enum component_case {P=0, PG, G};

void null_components(source_t *src, e_component_type component_type) {

  components_t *components;

  if (component_type == POINT) {
    components = &src->point_components;
    src->n_points = 0;
    src->n_point_lists = 0;
    src->n_point_powers = 0;
    src->n_point_curves = 0;
  }
  else if (component_type == GAUSSIAN){
    components = &src->gauss_components;
    src->n_gauss = 0;
    src->n_gauss_lists = 0;
    src->n_gauss_powers = 0;
    src->n_gauss_curves = 0;
  }
  else if (component_type == SHAPELET){
    components = &src->shape_components;
    src->n_shapes = 0;
    src->n_shape_coeffs = 0;
    src->n_shape_lists = 0;
    src->n_shape_powers = 0;
    src->n_shape_curves = 0;
  } else {
    //Shouldn't ever get here unless something that isn't a COMPONENT type
    //is entered
    return;
  }

  //Do the NULLing of arrays
  null_component_sky_model_arrays(components);

}



void increment_pointgauss(source_t *temp_cropped_src, source_t *cropped_src,
                          e_component_type component_type,
                          int power_iter, int curve_iter, int list_iter,
                          int num_chunk_power, int num_chunk_curve, int num_chunk_list){

  components_t *temp_components;
  components_t *components;
  if (component_type == POINT) {
    temp_components = &temp_cropped_src->point_components;
    components = &cropped_src->point_components;
    temp_cropped_src->n_points = num_chunk_power + num_chunk_curve + num_chunk_list;
    temp_cropped_src->n_point_powers = num_chunk_power;
    temp_cropped_src->n_point_curves = num_chunk_curve;
    temp_cropped_src->n_point_lists = num_chunk_list;
  }
  else if (component_type == GAUSSIAN){
    temp_components = &temp_cropped_src->gauss_components;
    components = &cropped_src->gauss_components;
    temp_cropped_src->n_gauss = num_chunk_power + num_chunk_curve + num_chunk_list;
    temp_cropped_src->n_gauss_powers = num_chunk_power;
    temp_cropped_src->n_gauss_curves = num_chunk_curve;
    temp_cropped_src->n_gauss_lists = num_chunk_list;

    //These are GAUSSIAN only params
    temp_components->majors = components->majors;
    temp_components->minors = components->minors;
    temp_components->pas = components->pas;
  }
  else {
    return;
  }

  //Some pointers we just copy, as we'll use indexing to access them later
  //When copying across to GPU memory, we'll only copy across what we need.
  //For now, just copy the pointer as that doesn't cost memory allocation
  temp_components->ras = components->ras;
  temp_components->decs = components->decs;
  temp_components->azs = components->azs;
  temp_components->zas = components->zas;
  temp_components->beam_has = components->beam_has;
  temp_components->beam_decs = components->beam_decs;

  //Others, we'll increment with pointer arithmatic so we can simply loop
  //over them (well be parallel over them on the GPU)

  temp_components->power_ref_stokesI = components->power_ref_stokesI + power_iter;
  temp_components->power_ref_stokesQ = components->power_ref_stokesQ + power_iter;
  temp_components->power_ref_stokesU = components->power_ref_stokesU + power_iter;
  temp_components->power_ref_stokesV = components->power_ref_stokesV + power_iter;
  temp_components->power_ref_freqs = components->power_ref_freqs + power_iter;
  temp_components->power_SIs = components->power_SIs + power_iter;
  temp_components->power_comp_inds = components->power_comp_inds + power_iter;

  temp_components->curve_ref_freqs = components->curve_ref_freqs + curve_iter;
  temp_components->curve_ref_stokesI = components->curve_ref_stokesI + curve_iter;
  temp_components->curve_ref_stokesQ = components->curve_ref_stokesQ + curve_iter;
  temp_components->curve_ref_stokesU = components->curve_ref_stokesU + curve_iter;
  temp_components->curve_ref_stokesV = components->curve_ref_stokesV + curve_iter;
  temp_components->curve_SIs = components->curve_SIs + curve_iter;
  temp_components->curve_qs = components->curve_qs + curve_iter;
  temp_components->curve_comp_inds = components->curve_comp_inds + curve_iter;

  //For the list flux stuff, each component can have any amount of
  //flux values, so we need to copy the full array pointers for that info
  //so we can easily access
  //What we will iterate are the arrays that list how many values each component
  //has, and where each component starts in the full freq/flux arrays
  //That will allow us to iterate over the correct components for this chunk
  temp_components->list_freqs = components->list_freqs;
  temp_components->list_stokesI = components->list_stokesI;
  temp_components->list_stokesQ = components->list_stokesQ;
  temp_components->list_stokesU = components->list_stokesU;
  temp_components->list_stokesV = components->list_stokesV;

  temp_components->num_list_values = components->num_list_values + list_iter;
  temp_components->list_start_indexes = components->list_start_indexes + list_iter;
  temp_components->list_comp_inds = components->list_comp_inds + list_iter;

}

/*
Here, given the overall lower and upper index in a given type of components,
work out how many of each flux type we have and increment the counters as
appropriate

Always order things as POWER_LAW, CURVED_POWER_LAW, LIST
*/
void increment_flux_type_counters(int * power_iter, int * curve_iter, int * list_iter,
                                  int * num_chunk_power, int * num_chunk_curve, int * num_chunk_list,
                                  int num_power, int num_curve, int num_list,
                                  int comps_per_chunk,
                                  int lower_comp_ind, int upper_comp_ind){

  int remainder = 0;
  int lower_flux_ind = 0;
  int upper_flux_ind = 0;

  // printf("LOWER COMP IND %d\n",lower_comp_ind );

  //Enough POWER_LAW to fill the whole chunk
  if (num_power > upper_comp_ind) {
    * num_chunk_power = comps_per_chunk;
    * power_iter = lower_comp_ind;
    * num_chunk_curve = 0;
    * num_chunk_list = 0;
  }
  //Not enough POWER_LAW to fill the whole chunk
  else {
    //There are enough POWER_LAW to partially fill chunk
    if (num_power >= lower_comp_ind) {
      * num_chunk_power = num_power - lower_comp_ind;
      * power_iter = lower_comp_ind;

      //How much is left to fill in this chunk
      remainder = comps_per_chunk - * num_chunk_power;

      //If there are enough CURVED_POWER_LAW to fill rest of the chunk
      if (num_curve >= remainder) {
        * num_chunk_curve = remainder;
        * curve_iter = 0;
        * num_chunk_list = 0;
      }
      //Not enough CURVED_POWER_LAW to fill rest of the chunk
      else {
        //There are some CURVED_POWER_LAW to add
        if (num_curve < remainder && num_curve != 0) {
          * num_chunk_curve = num_curve;
          * curve_iter = 0;
          remainder -= num_curve;
        }

        //There are enough LIST to fill the rest of the chunk
        if (num_list >= remainder ) {
          * num_chunk_list = remainder;
          * list_iter = 0;
        }
        //There are some LIST but not enough to fill the rest of the chunk
        else if (num_list != 0) {
          * num_chunk_list = num_list;
          * list_iter = 0;
        }

      }
    }
   //There aren't any POWER_LAW to put in this chunk
   //We may well have already chunked up a number of POWER_LAW so take
   //off the number of POWER_LAW from the lower_comp_ind, upper_comp_ind
    else {
      lower_flux_ind = lower_comp_ind - num_power;
      upper_flux_ind = upper_comp_ind - num_power;

      // printf("lower_comp_ind, upper_comp_ind %d %d\n",lower_comp_ind, upper_comp_ind );

      //There are enough CURVED_POWER_LAW to fill the rest of the chunk
      if (num_curve >= upper_flux_ind) {
        * num_chunk_curve = comps_per_chunk;
        * curve_iter = lower_flux_ind;
        * num_chunk_list = 0;
      } else {
        // printf("HERE BE WE\n");
        //There are some CURVED_POWER_LAW to add
        if (num_curve > lower_flux_ind && num_curve != 0) {
          * num_chunk_curve = num_curve - lower_flux_ind;
          * curve_iter = lower_flux_ind;
          remainder = comps_per_chunk - * num_chunk_curve;

          //There are enough LIST to fill rest of chunk
          if (num_list >= remainder) {
            * num_chunk_list = remainder;
            * list_iter = 0;
          }
          //There aren't enough LIST to fill chunk but there are some
          else {
            * num_chunk_list = num_list;
            * list_iter = 0;
          }

        }
        //There are no POWER_LAW or CURVED_POWER_LAW to add
        else {
        lower_flux_ind = lower_comp_ind - num_power - num_curve;
        upper_flux_ind = upper_comp_ind - num_power - num_curve;

          // printf("lower_comp_ind, upper_comp_ind %d %d\n",lower_comp_ind, upper_comp_ind );

          //There are enough LIST to fill the rest of the chunk
          if (num_list >= upper_flux_ind ) {
            * num_chunk_list = comps_per_chunk;
            * list_iter = lower_flux_ind;
          }
          //There are some LIST but not enough to fill the rest of the chunk
          else if (num_list > lower_flux_ind) {
            * num_chunk_list = num_list - lower_flux_ind;
            * list_iter = lower_flux_ind;
          }
        }
      }
    }
  }
}


void fill_chunk_src_with_pointgauss(source_t *temp_cropped_src,
     source_t *cropped_src, int chunk_ind, int comps_per_chunk,
     woden_settings_t *woden_settings, e_component_type comptype) {

  //Splitting POINTs and GAUSSIANS into lovely chunks that our GPU can chew
  //First we have to ascertain where in the chunking we are, and which type
  //of component we have to include

  //Lower and upper indexes of components covered in this chunk
  int lower_comp_ind = chunk_ind * comps_per_chunk;
  int upper_comp_ind = (chunk_ind + 1) * comps_per_chunk;

  //comp_case is used to ascertain what combo of POINT,GAUSSIAN,SHAPELET we have
  // int comp_case = -1;

  //These ints are used to do pointer arithmatic to grab the correct portions
  //of arrays out of `cropped_src` and into `temp_cropped_src`
  int power_iter = 0;
  int curve_iter = 0;
  int list_iter = 0;

  int num_chunk_power = 0;
  int num_chunk_curve = 0;
  int num_chunk_list = 0;

  int n_powers = 0;
  int n_curves = 0;
  int n_lists = 0;

  if (comptype == POINT) {
    n_powers = cropped_src->n_point_powers;
    n_curves = cropped_src->n_point_curves;
    n_lists = cropped_src->n_point_lists;
  } else if (comptype == GAUSSIAN){
    n_powers = cropped_src->n_gauss_powers;
    n_curves = cropped_src->n_gauss_curves;
    n_lists = cropped_src->n_gauss_lists;
  }

  increment_flux_type_counters(&power_iter, &curve_iter, &list_iter,
                &num_chunk_power, &num_chunk_curve, &num_chunk_list,
                n_powers, n_curves, n_lists,
                comps_per_chunk, lower_comp_ind, upper_comp_ind);

  if (comptype == POINT) {
    null_components(temp_cropped_src, GAUSSIAN);
    null_components(temp_cropped_src, SHAPELET);
    increment_pointgauss(temp_cropped_src, cropped_src, POINT,
                    power_iter, curve_iter, list_iter,
                    num_chunk_power, num_chunk_curve, num_chunk_list);
    temp_cropped_src->point_components.num_primarybeam_values = temp_cropped_src->n_points*woden_settings->num_freqs*woden_settings->num_time_steps;

  }
  else if (comptype == GAUSSIAN) {
    null_components(temp_cropped_src, POINT);
    null_components(temp_cropped_src, SHAPELET);
    increment_pointgauss(temp_cropped_src, cropped_src, GAUSSIAN,
                    power_iter, curve_iter, list_iter,
                    num_chunk_power, num_chunk_curve, num_chunk_list);
    temp_cropped_src->gauss_components.num_primarybeam_values = temp_cropped_src->n_gauss*woden_settings->num_freqs*woden_settings->num_time_steps;
  }

  temp_cropped_src->n_comps = temp_cropped_src->n_points + temp_cropped_src->n_gauss;

}



void increment_shapelet(source_t *temp_cropped_src, source_t *cropped_src,
                        int * shape_iter){
  //increment the required pointers to point at the beginning of the next chunk

  //for shapelets, we chunk over coeffs, not components, so need all the
  //ras,azs, minors, etc each time, just iterate the coeffs, n1s, n2s
  //do this because there should be many coeffs per components, and not
  //many components, due to the nature of shapelets
  temp_cropped_src->n_shapes = cropped_src->n_shapes;
  temp_cropped_src->n_comps = cropped_src->n_shapes;
  temp_cropped_src->n_shape_powers = cropped_src->n_shape_powers;
  temp_cropped_src->n_shape_curves = cropped_src->n_shape_curves;
  temp_cropped_src->n_shape_lists = cropped_src->n_shape_lists;

  temp_cropped_src->shape_components.ras = cropped_src->shape_components.ras;
  temp_cropped_src->shape_components.decs = cropped_src->shape_components.decs;
  temp_cropped_src->shape_components.majors = cropped_src->shape_components.majors;
  temp_cropped_src->shape_components.minors = cropped_src->shape_components.minors;
  temp_cropped_src->shape_components.pas = cropped_src->shape_components.pas;


  temp_cropped_src->shape_components.power_ref_stokesI = cropped_src->shape_components.power_ref_stokesI;
  temp_cropped_src->shape_components.power_ref_stokesQ = cropped_src->shape_components.power_ref_stokesQ;
  temp_cropped_src->shape_components.power_ref_stokesU = cropped_src->shape_components.power_ref_stokesU;
  temp_cropped_src->shape_components.power_ref_stokesV = cropped_src->shape_components.power_ref_stokesV;
  temp_cropped_src->shape_components.power_ref_freqs = cropped_src->shape_components.power_ref_freqs;
  temp_cropped_src->shape_components.power_SIs = cropped_src->shape_components.power_SIs;
  temp_cropped_src->shape_components.power_comp_inds = cropped_src->shape_components.power_comp_inds;

  temp_cropped_src->shape_components.curve_ref_freqs = cropped_src->shape_components.curve_ref_freqs;
  temp_cropped_src->shape_components.curve_ref_stokesI = cropped_src->shape_components.curve_ref_stokesI;
  temp_cropped_src->shape_components.curve_ref_stokesQ = cropped_src->shape_components.curve_ref_stokesQ;
  temp_cropped_src->shape_components.curve_ref_stokesU = cropped_src->shape_components.curve_ref_stokesU;
  temp_cropped_src->shape_components.curve_ref_stokesV = cropped_src->shape_components.curve_ref_stokesV;
  temp_cropped_src->shape_components.curve_SIs = cropped_src->shape_components.curve_SIs;
  temp_cropped_src->shape_components.curve_qs = cropped_src->shape_components.curve_qs;
  temp_cropped_src->shape_components.curve_comp_inds = cropped_src->shape_components.curve_comp_inds;

  temp_cropped_src->shape_components.list_freqs = cropped_src->shape_components.list_freqs;
  temp_cropped_src->shape_components.list_stokesI = cropped_src->shape_components.list_stokesI;
  temp_cropped_src->shape_components.list_stokesQ = cropped_src->shape_components.list_stokesQ;
  temp_cropped_src->shape_components.list_stokesU = cropped_src->shape_components.list_stokesU;
  temp_cropped_src->shape_components.list_stokesV = cropped_src->shape_components.list_stokesV;
  temp_cropped_src->shape_components.num_list_values = cropped_src->shape_components.num_list_values;
  temp_cropped_src->shape_components.list_start_indexes = cropped_src->shape_components.list_start_indexes;
  temp_cropped_src->shape_components.list_comp_inds = cropped_src->shape_components.list_comp_inds;


  //only chunk over coeffs, so we need all the az / za for every chunk,
  //so we don't iterate the pointer here
  temp_cropped_src->shape_components.azs = cropped_src->shape_components.azs;
  temp_cropped_src->shape_components.zas = cropped_src->shape_components.zas;
  temp_cropped_src->shape_components.beam_has = cropped_src->shape_components.beam_has;
  temp_cropped_src->shape_components.beam_decs = cropped_src->shape_components.beam_decs;


  temp_cropped_src->shape_components.shape_coeffs = cropped_src->shape_components.shape_coeffs + * shape_iter;
  temp_cropped_src->shape_components.n1s = cropped_src->shape_components.n1s + * shape_iter;
  temp_cropped_src->shape_components.n2s = cropped_src->shape_components.n2s + * shape_iter;
  temp_cropped_src->shape_components.param_indexes = cropped_src->shape_components.param_indexes + * shape_iter;

}



void fill_chunk_src_with_shapelets(source_t *temp_cropped_src,
     source_t *cropped_src, int chunk_ind, int coeffs_per_chunk,
     woden_settings_t *woden_settings) {

  //Upper indexes of components covered in this chunk
  int upper_comp_ind = (chunk_ind + 1) * coeffs_per_chunk;

  //These ints are used to do pointer arithmatic to grab the correct portions
  //of arrays out of `cropped_src` and into `temp_cropped_src`
  int shape_iter = chunk_ind * coeffs_per_chunk;

  //If there are enough coeffs to fill the chunk?
  if (cropped_src->n_shape_coeffs >= upper_comp_ind){
    temp_cropped_src->n_shape_coeffs = coeffs_per_chunk;

  }
  else { //Otherwise just set the remainder
    temp_cropped_src->n_shape_coeffs = cropped_src->n_shape_coeffs % coeffs_per_chunk;
  }

  null_components(temp_cropped_src, POINT);
  null_components(temp_cropped_src, GAUSSIAN);
  increment_shapelet(temp_cropped_src, cropped_src,
                    &shape_iter);
  temp_cropped_src->shape_components.num_primarybeam_values = temp_cropped_src->n_shapes*woden_settings->num_freqs*woden_settings->num_time_steps;
}

source_catalogue_t * create_chunked_sky_models(source_t *cropped_src,
                                               woden_settings_t *woden_settings) {

  //Ok, so we split POINT and GAUSSIANs up by whole COMPONENT into the chunked
  //sky models. We split SHAPELETs up by their bapower_SIs function information,
  //which I shorthand to 'coeff' here. I've had trouble with running SHAPELETs
  //at the same time as many POINT/GAUSSIANs, so split them up separately

  //Chunking size says how many visibility calculations we can do on the GPU
  //at once, where one calculation here means:
  // one component * one time step * one frequency step * one baseline
  //Here we work out how many components we can stick on the GPU and remain
  //below the `chunking_size`.

  long int chunking_size = woden_settings->chunking_size;
  int num_baselines = woden_settings->num_baselines;
  int num_freqs = woden_settings->num_freqs;
  int num_time_steps = woden_settings->num_time_steps;

  // int comps_per_chunk = (int)floorf((float)chunking_size / (float)(num_baselines * num_freqs * num_time_steps));

  // printf("\t ++ chunking_size, num_baselines, num_freqs, num_time_steps %ld %d %d %d\n", chunking_size, num_baselines, num_freqs, num_time_steps );

  int comps_per_chunk = (int)floorf((float)(chunking_size / (num_baselines * num_freqs * num_time_steps)));

  if (comps_per_chunk < 1) {
    comps_per_chunk = 1;
  }

  // printf("\t ++ comps_per_chunk %d\n", comps_per_chunk);

  //How many chunks we need for the POINT/GAUSS, and the SHAPELETs
  int num_point_chunks = (int)ceilf((float)cropped_src->n_points / (float)comps_per_chunk);
  int num_gauss_chunks = (int)ceilf((float)cropped_src->n_gauss / (float)comps_per_chunk);
  int num_shape_chunks = (int)ceilf((float)cropped_src->n_shape_coeffs / (float)comps_per_chunk);

  int num_chunks = num_point_chunks + num_gauss_chunks + num_shape_chunks;

  source_catalogue_t *chunked_sky_models;
  chunked_sky_models = malloc(sizeof(source_catalogue_t));

  chunked_sky_models->num_sources = num_chunks;
  chunked_sky_models->num_shapelets = 0;
  chunked_sky_models->sources = malloc(num_chunks*sizeof(source_t));

  source_t *temp_cropped_src = malloc(sizeof(source_t));

  printf("Number of chunks required is %d\n", num_chunks);
  printf("Chunking sky model.. ");
  //Chunk the POINT/GAUSS first
  for (int comp_chunk = 0; comp_chunk < num_point_chunks; comp_chunk++) {

    fill_chunk_src_with_pointgauss(temp_cropped_src, cropped_src, comp_chunk,
                                   comps_per_chunk, woden_settings, POINT);

    chunked_sky_models->sources[comp_chunk] = *temp_cropped_src;

  }

  for (int comp_chunk = 0; comp_chunk < num_gauss_chunks; comp_chunk++) {

    fill_chunk_src_with_pointgauss(temp_cropped_src, cropped_src, comp_chunk,
                                   comps_per_chunk, woden_settings, GAUSSIAN);

    chunked_sky_models->sources[num_point_chunks + comp_chunk] = *temp_cropped_src;

  }


  for (int coeff_chunk = 0; coeff_chunk < num_shape_chunks; coeff_chunk++) {

    fill_chunk_src_with_shapelets(temp_cropped_src, cropped_src, coeff_chunk,
                              comps_per_chunk, woden_settings);

    //Add the number of shapelets onto the full source catalogue value
    //so we know if we need to setup shapelet bapower_SIs functions in GPU memory
    //or not
    chunked_sky_models->num_shapelets += temp_cropped_src->n_shapes;
    chunked_sky_models->sources[num_point_chunks + num_gauss_chunks + coeff_chunk] = *temp_cropped_src;

  }

  printf("Sky model chunked.\n");
  // free(temp_cropped_src);
  return chunked_sky_models;
}

void remap_common_comp_values(components_t *remapped_components,
                              components_t *chunked_components,
                              int new_comp_index, int old_comp_index,
                              int num_time_steps, e_beamtype beamtype,
                              e_component_type comptype) {

  // printf("remap common %d %d\n",new_comp_index,old_comp_index );

  remapped_components->ras[new_comp_index] = chunked_components->ras[old_comp_index];
  remapped_components->decs[new_comp_index] = chunked_components->decs[old_comp_index];

  if (comptype == GAUSSIAN) {
    remapped_components->pas[new_comp_index] = chunked_components->pas[old_comp_index];
    remapped_components->majors[new_comp_index] = chunked_components->majors[old_comp_index];
    remapped_components->minors[new_comp_index] = chunked_components->minors[old_comp_index];
  }

  int old_time_ind, new_time_ind;

  for (int time_ind = 0; time_ind < num_time_steps; time_ind++) {
    old_time_ind = old_comp_index*num_time_steps + time_ind;
    new_time_ind = new_comp_index*num_time_steps + time_ind;

    // int n_comps = 76894;
    // if (new_time_ind > n_comps) {
      // printf("OH WE'RE IN HERE %d %d\n", new_time_ind, n_comps);
    // }

    // printf("OH WE'RE IN HERE %d %d\n", new_comp_index, old_comp_index);

    remapped_components->azs[new_time_ind] = chunked_components->azs[old_time_ind];
    remapped_components->zas[new_time_ind] = chunked_components->zas[old_time_ind];

    if (beamtype == MWA_ANALY || beamtype == GAUSS_BEAM) {

        remapped_components->beam_has[new_time_ind] = chunked_components->beam_has[old_time_ind];
        remapped_components->beam_decs[new_time_ind] = chunked_components->beam_decs[old_time_ind];
    }
  }

}

/*
As we chunk over POINT/GAUSS components, but over coefficients for SHAPELETS,
only have to do remapping for POINT/GAUSS in the below
*/
void remap_components(source_t *remapped_source, source_t *chunked_source,
                      int num_time_steps, e_beamtype beamtype,
                      e_component_type comptype) {
  // printf("DOING A REMAP\n");

  components_t *remapped_components = NULL;
  components_t *chunked_components = NULL;
  int n_powers = 0, n_curves = 0, n_lists = 0;

  if (comptype == POINT) {
    remapped_components = &remapped_source->point_components;
    chunked_components = &chunked_source->point_components;
    n_powers = chunked_source->n_point_powers;
    n_curves = chunked_source->n_point_curves;
    n_lists = chunked_source->n_point_lists;
  }
  else if (comptype == GAUSSIAN) {
    remapped_components = &remapped_source->gauss_components;
    chunked_components = &chunked_source->gauss_components;
    n_powers = chunked_source->n_gauss_powers;
    n_curves = chunked_source->n_gauss_curves;
    n_lists = chunked_source->n_gauss_lists;
  }

  int new_comp_index = 0;
  int old_comp_index;
  int n_comps = n_powers + n_curves + n_lists;

  remapped_components->num_primarybeam_values = chunked_components->num_primarybeam_values;

  //Malloc things that we have to remap
  remapped_components->ras = malloc(n_comps*sizeof(double));
  remapped_components->decs = malloc(n_comps*sizeof(double));
  remapped_components->azs = malloc(num_time_steps*n_comps*sizeof(user_precision_t));
  remapped_components->zas = malloc(num_time_steps*n_comps*sizeof(user_precision_t));

  if (comptype == GAUSSIAN) {
    remapped_components->pas = malloc(n_comps*sizeof(user_precision_t));
    remapped_components->majors = malloc(n_comps*sizeof(user_precision_t));
    remapped_components->minors = malloc(n_comps*sizeof(user_precision_t));
  }

  //Only need these arrays for certain types of beams
  if (beamtype == MWA_ANALY || beamtype == GAUSS_BEAM) {
    remapped_components->beam_has = malloc(num_time_steps*n_comps*sizeof(double));
    remapped_components->beam_decs = malloc(num_time_steps*n_comps*sizeof(double));
  }

  if (n_powers > 0) {
    // printf("POWER new_comp_index %d\n", new_comp_index);

    //These had pointer arithmatic done on them so we can copy the pointer
    remapped_components->power_ref_stokesI = chunked_components->power_ref_stokesI;
    remapped_components->power_ref_stokesQ = chunked_components->power_ref_stokesQ;
    remapped_components->power_ref_stokesU = chunked_components->power_ref_stokesU;
    remapped_components->power_ref_stokesV = chunked_components->power_ref_stokesV;
    remapped_components->power_ref_freqs = chunked_components->power_ref_freqs;
    remapped_components->power_SIs = chunked_components->power_SIs;

    //Need a new mapping so time for mallocing
    remapped_components->power_comp_inds = malloc(n_powers*sizeof(int));

    for (int pow_ind = 0; pow_ind < n_powers; pow_ind++) {
      old_comp_index = chunked_components->power_comp_inds[pow_ind];

      remapped_components->power_comp_inds[pow_ind] = new_comp_index;

      remap_common_comp_values(remapped_components, chunked_components,
                              new_comp_index, old_comp_index, num_time_steps,
                              beamtype, comptype);

      new_comp_index += 1;
    }
  }

  if (n_curves > 0) {
    // printf("CURVE new_comp_index %d\n", new_comp_index);

    //These had pointer arithmatic done on them so we can copy the pointer
    remapped_components->curve_ref_stokesI = chunked_components->curve_ref_stokesI;
    remapped_components->curve_ref_stokesQ = chunked_components->curve_ref_stokesQ;
    remapped_components->curve_ref_stokesU = chunked_components->curve_ref_stokesU;
    remapped_components->curve_ref_stokesV = chunked_components->curve_ref_stokesV;
    remapped_components->curve_ref_freqs = chunked_components->curve_ref_freqs;
    remapped_components->curve_SIs = chunked_components->curve_SIs;
    remapped_components->curve_qs = chunked_components->curve_qs;

    //Need a new mapping so time for mallocing
    remapped_components->curve_comp_inds = malloc(n_curves*sizeof(int));

    for (int cur_ind = 0; cur_ind < n_curves; cur_ind++) {
      old_comp_index = chunked_components->curve_comp_inds[cur_ind];

      remapped_components->curve_comp_inds[cur_ind] = new_comp_index;

      remap_common_comp_values(remapped_components, chunked_components,
                              new_comp_index, old_comp_index, num_time_steps,
                              beamtype, comptype);

      new_comp_index += 1;
    }
  }

  //Ok, we have to grab the correct flux list entries, so we are only
  //copying what we need into GPU memory. This will involve some iteration
  //and index matching siggggggggggggggh.
  if (n_lists > 0) {

    //These had pointer arithmatic done on them so we can copy the pointer
    remapped_components->num_list_values = chunked_components->num_list_values;

    //Sum up everything in num_list_values gives us how much we need to malloc
    int total_num_flux_entires = 0;
    for (int list_comp = 0; list_comp < n_lists; list_comp++) {
      total_num_flux_entires += remapped_components->num_list_values[list_comp];
      // printf("Adding %d\n",remapped_components->num_list_values[list_comp] );
    }

    remapped_components->total_num_flux_entires = total_num_flux_entires;

    remapped_components->list_freqs = malloc(total_num_flux_entires*sizeof(double));
    remapped_components->list_stokesI = malloc(total_num_flux_entires*sizeof(user_precision_t));
    remapped_components->list_stokesQ = malloc(total_num_flux_entires*sizeof(user_precision_t));
    remapped_components->list_stokesU = malloc(total_num_flux_entires*sizeof(user_precision_t));
    remapped_components->list_stokesV = malloc(total_num_flux_entires*sizeof(user_precision_t));

    remapped_components->list_start_indexes = malloc(n_lists*sizeof(int));
    remapped_components->list_comp_inds = malloc(n_lists*sizeof(int));

    int remap_list_start = 0;

    int new_list_entry_ind = 0;

    int old_list_entry_ind;

    //Loop over all list components - each has any number of flux entries,
    //given by remapped_components->num_list_values
    for (int list_comp = 0; list_comp < n_lists; list_comp++) {

      old_comp_index = chunked_components->list_comp_inds[list_comp];

      for (int list_ind = 0; list_ind < remapped_components->num_list_values[list_comp]; list_ind++) {

        old_list_entry_ind = chunked_components->list_start_indexes[list_comp] + list_ind;

        remapped_components->list_freqs[new_list_entry_ind] = chunked_components->list_freqs[old_list_entry_ind];
        remapped_components->list_stokesI[new_list_entry_ind] = chunked_components->list_stokesI[old_list_entry_ind];
        remapped_components->list_stokesQ[new_list_entry_ind] = chunked_components->list_stokesQ[old_list_entry_ind];
        remapped_components->list_stokesU[new_list_entry_ind] = chunked_components->list_stokesU[old_list_entry_ind];
        remapped_components->list_stokesV[new_list_entry_ind] = chunked_components->list_stokesV[old_list_entry_ind];

        new_list_entry_ind += 1;
      }

      remapped_components->list_start_indexes[list_comp] = remap_list_start;
      remap_list_start += remapped_components->num_list_values[list_comp];

      remapped_components->list_comp_inds[list_comp] = new_comp_index;

      remap_common_comp_values(remapped_components, chunked_components,
                              new_comp_index, old_comp_index, num_time_steps,
                              beamtype, comptype);

      new_comp_index += 1;
      // printf("new_comp_index %d\n",new_comp_index );

    }//END loop over LIST components
  }//END if n_lists > 0
}


//We don't want to have to copy across all of the component information
//onto the GPU at once, so we need to remap only what we need
//We have to do a remap, not just pointer arithmatic, so requires a memory
//allocation. Makes sense to only do this before copying across to the GPU
//and then releasing straight after, so don't do this during the initial
//chunking step. Just do this before each iteration over the chunks
//in calculate_visibilities::calculate_visibilities
void remap_source_for_gpu(source_t *remapped_source, source_t *chunked_source,
                          int num_time_steps, e_beamtype beamtype) {

  source_zero_counters_and_null_components(remapped_source);

  if (chunked_source->n_points > 0) {
    remap_components(remapped_source, chunked_source, num_time_steps, beamtype,
                     POINT);
  }
  if (chunked_source->n_gauss > 0) {
    remap_components(remapped_source, chunked_source, num_time_steps, beamtype,
                   GAUSSIAN);
  }

  //Ok, it's assumed that there will never be millions of SHAPELET components,
  //due to their nature (big on the sky). So we chunk over the coefficients of
  //all SHAPELETs combined, as those should be big arrays.  So we SHOULD be able
  //to just copy the chunked pointers here, and all should be able to fit on
  //the GPU

  remapped_source->shape_components = chunked_source->shape_components;

  if (chunked_source->n_shape_lists > 0) {

    //Sum up everything in num_list_values gives us how much we need to malloc
    int total_num_flux_entires = 0;
    for (int list_comp = 0; list_comp < chunked_source->n_shape_lists; list_comp++) {
    // printf("REMAP num list things %d\n",chunked_source->shape_components.num_list_values[list_comp] );
      total_num_flux_entires += chunked_source->shape_components.num_list_values[list_comp];

    }
    remapped_source->shape_components.total_num_flux_entires = total_num_flux_entires;
  }

  remapped_source->n_comps = chunked_source->n_comps;
  remapped_source->n_points = chunked_source->n_points;
  remapped_source->n_point_lists = chunked_source->n_point_lists;
  remapped_source->n_point_powers = chunked_source->n_point_powers;
  remapped_source->n_point_curves = chunked_source->n_point_curves;

  remapped_source->n_gauss = chunked_source->n_gauss;
  remapped_source->n_gauss_lists = chunked_source->n_gauss_lists;
  remapped_source->n_gauss_powers = chunked_source->n_gauss_powers;
  remapped_source->n_gauss_curves = chunked_source->n_gauss_curves;

  remapped_source->n_shapes = chunked_source->n_shapes;
  remapped_source->n_shape_lists = chunked_source->n_shape_lists;
  remapped_source->n_shape_powers = chunked_source->n_shape_powers;
  remapped_source->n_shape_curves = chunked_source->n_shape_curves;
  remapped_source->n_shape_coeffs = chunked_source->n_shape_coeffs;

}

/*
Only free things that would have been malloc'd by remap_source_for_gpu
A large number should have been pointer arithmatic, so don't want to
free the original memory - need it for other iterations
*/
void free_remapped_components_for_gpu(components_t *components,
      int num_lists, int num_powers,  int num_curves,
      e_beamtype beamtype, e_component_type comptype){

  //Instrinsic to COMPONENT values
  free(components->ras);
  free(components->decs);

  free(components->azs);
  free(components->zas);

  if(comptype == GAUSSIAN) {
    free(components->pas);
    free(components->majors);
    free(components->minors);
  }

  if (beamtype == MWA_ANALY || beamtype == GAUSS_BEAM) {
    free(components->beam_has);
    free(components->beam_decs);
  }

  if (num_powers > 0){
    free(components->power_comp_inds);
  }

  if (num_curves > 0){
    free(components->curve_comp_inds);
  }

  // printf("Trying to free the LIST stuff\n");
  if (num_lists > 0){
    free(components->list_freqs);
    free(components->list_stokesI);
    free(components->list_stokesQ);
    free(components->list_stokesU);
    free(components->list_stokesV);
    free(components->list_comp_inds);
    free(components->list_start_indexes);
    //This is handled by pointer arithmatic, so don't need to free
    // free(components->num_list_values);
  }
}


void free_remapped_source_for_gpu(source_t *source, e_beamtype beamtype){
  if (source->n_points > 0){
    free_remapped_components_for_gpu(&source->point_components, source->n_point_lists,
        source->n_point_powers, source->n_point_curves, beamtype, POINT);
  }

  if (source->n_gauss > 0){
    free_remapped_components_for_gpu(&source->gauss_components, source->n_gauss_lists,
        source->n_gauss_powers, source->n_gauss_curves, beamtype, GAUSSIAN);
  }

  free(source);

}
