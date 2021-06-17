#include <stdio.h>
#include "read_and_write.h"

//Switch through different POINT and GAUSSIAN chunking states - if P, only
//POINTs in a chunk, PG both POINT and GAUSSIAN, G just GAUSSIAN
enum component_case {P=0, PG, G};

void null_point_comps(catsource_t *temp_cropped_src){
  temp_cropped_src->point_ras = NULL;
  temp_cropped_src->point_decs = NULL;
  temp_cropped_src->point_ref_stokesI = NULL;
  temp_cropped_src->point_ref_stokesQ = NULL;
  temp_cropped_src->point_ref_stokesU = NULL;
  temp_cropped_src->point_ref_stokesV = NULL;
  temp_cropped_src->point_ref_freqs = NULL;
  temp_cropped_src->point_SIs = NULL;
  temp_cropped_src->point_azs = NULL;
  temp_cropped_src->point_zas = NULL;
  temp_cropped_src->cos_point_para_angs = NULL;
  temp_cropped_src->sin_point_para_angs = NULL;
  temp_cropped_src->point_gaussbeam_has = NULL;
  temp_cropped_src->point_gaussbeam_decs = NULL;

  temp_cropped_src->n_points = 0;
  temp_cropped_src->num_point_primarybeam_values = 0;

}

void null_gauss_comps(catsource_t *temp_cropped_src){
  temp_cropped_src->gauss_ras = NULL;
  temp_cropped_src->gauss_decs = NULL;
  temp_cropped_src->gauss_ref_stokesI = NULL;
  temp_cropped_src->gauss_ref_stokesQ = NULL;
  temp_cropped_src->gauss_ref_stokesU = NULL;
  temp_cropped_src->gauss_ref_stokesV = NULL;
  temp_cropped_src->gauss_ref_freqs = NULL;
  temp_cropped_src->gauss_SIs = NULL;
  temp_cropped_src->gauss_majors = NULL;
  temp_cropped_src->gauss_minors = NULL;
  temp_cropped_src->gauss_pas = NULL;
  temp_cropped_src->gauss_azs = NULL;
  temp_cropped_src->gauss_zas = NULL;
  temp_cropped_src->cos_gauss_para_angs = NULL;
  temp_cropped_src->sin_gauss_para_angs = NULL;
  temp_cropped_src->gauss_gaussbeam_has = NULL;
  temp_cropped_src->gauss_gaussbeam_decs = NULL;

  temp_cropped_src->n_gauss = 0;
  temp_cropped_src->num_gauss_primarybeam_values = 0;

}

void null_shapelet_comps(catsource_t *temp_cropped_src){
  temp_cropped_src->shape_ras = NULL;
  temp_cropped_src->shape_decs = NULL;
  temp_cropped_src->shape_ref_stokesI = NULL;
  temp_cropped_src->shape_ref_stokesQ = NULL;
  temp_cropped_src->shape_ref_stokesU = NULL;
  temp_cropped_src->shape_ref_stokesV = NULL;
  temp_cropped_src->shape_ref_freqs = NULL;
  temp_cropped_src->shape_SIs = NULL;
  temp_cropped_src->shape_majors = NULL;
  temp_cropped_src->shape_minors = NULL;
  temp_cropped_src->shape_pas = NULL;
  temp_cropped_src->shape_n1s = NULL;
  temp_cropped_src->shape_n2s = NULL;
  temp_cropped_src->shape_coeffs = NULL;
  temp_cropped_src->shape_param_indexes = NULL;
  temp_cropped_src->shape_azs = NULL;
  temp_cropped_src->shape_zas= NULL;
  temp_cropped_src->cos_shape_para_angs = NULL;
  temp_cropped_src->sin_shape_para_angs = NULL;
  temp_cropped_src->shape_gaussbeam_has = NULL;
  temp_cropped_src->shape_gaussbeam_decs = NULL;

  temp_cropped_src->n_shapes = 0;
  temp_cropped_src->n_shape_coeffs = 0;
  temp_cropped_src->num_shape_primarybeam_values = 0;

}

void increment_point(catsource_t *temp_cropped_src, catsource_t *cropped_src,
                     int * point_iter, int num_time_steps){
  //increment the required pointers to point at the beginning of the next chunk
  temp_cropped_src->point_ras = cropped_src->point_ras + * point_iter;
  temp_cropped_src->point_decs = cropped_src->point_decs + * point_iter;
  temp_cropped_src->point_ref_stokesI = cropped_src->point_ref_stokesI + * point_iter;
  temp_cropped_src->point_ref_stokesQ = cropped_src->point_ref_stokesQ + * point_iter;
  temp_cropped_src->point_ref_stokesU = cropped_src->point_ref_stokesU + * point_iter;
  temp_cropped_src->point_ref_stokesV = cropped_src->point_ref_stokesV + * point_iter;
  temp_cropped_src->point_ref_freqs = cropped_src->point_ref_freqs + * point_iter;
  temp_cropped_src->point_SIs = cropped_src->point_SIs + * point_iter;

  temp_cropped_src->point_azs = cropped_src->point_azs + (num_time_steps * * point_iter);
  temp_cropped_src->point_zas = cropped_src->point_zas + (num_time_steps * * point_iter);
  temp_cropped_src->sin_point_para_angs = cropped_src->sin_point_para_angs + (num_time_steps * * point_iter);
  temp_cropped_src->cos_point_para_angs = cropped_src->cos_point_para_angs + (num_time_steps * * point_iter);
  temp_cropped_src->point_gaussbeam_has = cropped_src->point_gaussbeam_has + (num_time_steps * * point_iter);
  temp_cropped_src->point_gaussbeam_decs = cropped_src->point_gaussbeam_decs + (num_time_steps * * point_iter);

}

void increment_gauss(catsource_t *temp_cropped_src, catsource_t *cropped_src,
                     int * gauss_iter, int num_time_steps){

  //increment the required pointers to point at the beginning of the next chunk
  temp_cropped_src->gauss_ras = cropped_src->gauss_ras + * gauss_iter;
  temp_cropped_src->gauss_decs = cropped_src->gauss_decs + * gauss_iter;
  temp_cropped_src->gauss_ref_stokesI = cropped_src->gauss_ref_stokesI + * gauss_iter;
  temp_cropped_src->gauss_ref_stokesQ = cropped_src->gauss_ref_stokesQ + * gauss_iter;
  temp_cropped_src->gauss_ref_stokesU = cropped_src->gauss_ref_stokesU + * gauss_iter;
  temp_cropped_src->gauss_ref_stokesV = cropped_src->gauss_ref_stokesV + * gauss_iter;
  temp_cropped_src->gauss_ref_freqs = cropped_src->gauss_ref_freqs + * gauss_iter;
  temp_cropped_src->gauss_SIs = cropped_src->gauss_SIs + * gauss_iter;
  temp_cropped_src->gauss_majors = cropped_src->gauss_majors + * gauss_iter;
  temp_cropped_src->gauss_minors = cropped_src->gauss_minors + * gauss_iter;
  temp_cropped_src->gauss_pas = cropped_src->gauss_pas + * gauss_iter;

  temp_cropped_src->gauss_azs = cropped_src->gauss_azs + (num_time_steps * * gauss_iter);
  temp_cropped_src->gauss_zas = cropped_src->gauss_zas + (num_time_steps * * gauss_iter);
  temp_cropped_src->sin_gauss_para_angs = cropped_src->sin_gauss_para_angs + (num_time_steps * * gauss_iter);
  temp_cropped_src->cos_gauss_para_angs = cropped_src->cos_gauss_para_angs + (num_time_steps * * gauss_iter);
  temp_cropped_src->gauss_gaussbeam_has = cropped_src->gauss_gaussbeam_has + (num_time_steps * * gauss_iter);
  temp_cropped_src->gauss_gaussbeam_decs = cropped_src->gauss_gaussbeam_decs + (num_time_steps * * gauss_iter);

}


void increment_shapelet(catsource_t *temp_cropped_src, catsource_t *cropped_src,
                        int * shape_iter, int num_time_steps){
  //increment the required pointers to point at the beginning of the next chunk

  //for shapelets, we chunk over coeffs, not components, so need all the
  //ras,azs, minors, etc each time, just iterate the coeffs, n1s, n2s
  //do this because there should be many coeffs per components, and not
  //many components, due to the nature of shapelets
  temp_cropped_src->n_shapes = cropped_src->n_shapes;

  temp_cropped_src->shape_ras = cropped_src->shape_ras;
  temp_cropped_src->shape_decs = cropped_src->shape_decs;
  temp_cropped_src->shape_ref_stokesI = cropped_src->shape_ref_stokesI;
  temp_cropped_src->shape_ref_stokesQ = cropped_src->shape_ref_stokesQ;
  temp_cropped_src->shape_ref_stokesU = cropped_src->shape_ref_stokesU;
  temp_cropped_src->shape_ref_stokesV = cropped_src->shape_ref_stokesV;
  temp_cropped_src->shape_ref_freqs = cropped_src->shape_ref_freqs;
  temp_cropped_src->shape_SIs = cropped_src->shape_SIs;

  temp_cropped_src->shape_majors = cropped_src->shape_majors;
  temp_cropped_src->shape_minors = cropped_src->shape_minors;
  temp_cropped_src->shape_pas = cropped_src->shape_pas;

  temp_cropped_src->shape_coeffs = cropped_src->shape_coeffs + * shape_iter;
  temp_cropped_src->shape_n1s = cropped_src->shape_n1s + * shape_iter;
  temp_cropped_src->shape_n2s = cropped_src->shape_n2s + * shape_iter;
  temp_cropped_src->shape_param_indexes = cropped_src->shape_param_indexes + * shape_iter;

  //only chunk over coeffs, so we need all the az / za for every chunk,
  //so we don't iterate the pointer here

  temp_cropped_src->shape_azs = cropped_src->shape_azs;
  temp_cropped_src->shape_zas = cropped_src->shape_zas;
  temp_cropped_src->sin_shape_para_angs = cropped_src->sin_shape_para_angs;
  temp_cropped_src->cos_shape_para_angs = cropped_src->cos_shape_para_angs;
  temp_cropped_src->shape_gaussbeam_has = cropped_src->shape_gaussbeam_has;
  temp_cropped_src->shape_gaussbeam_decs = cropped_src->shape_gaussbeam_decs;
}

void fill_chunk_src_with_pointgauss(catsource_t *temp_cropped_src,
     catsource_t *cropped_src, int chunk_ind, int comps_per_chunk,
     woden_settings_t *woden_settings) {

  //Splitting POINTs and GAUSSIANS into lovely chunks that our GPU can chew
  //First we have to ascertain where in the chunking we are, and which type
  //of component we have to include

  //Lower and upper indexes of components covered in this chunk
  int lower_comp_ind = chunk_ind * comps_per_chunk;
  int upper_comp_ind = (chunk_ind + 1) * comps_per_chunk;

  //comp_case is used to ascertain what combo of POINT,GAUSSIAN,SHAPELET we have
  int comp_case = -1;

  //These ints are used to do pointer arithmatic to grab the correct portions
  //of arrays out of `cropped_src` and into `temp_cropped_src`
  int point_iter;
  int gauss_iter;

  //BEGIN work out what component type and how many of each component
  //types we need to shove into the temp_cropped_src

  //If there are enough POINTs to fill the whole chunk
  if (cropped_src->n_points >= upper_comp_ind){
    comp_case = P;
    temp_cropped_src->n_points = comps_per_chunk;
    point_iter = chunk_ind * comps_per_chunk;
  }
  //If chunk contains POINT components, but there is space for GAUSSIANs
  else if ((cropped_src->n_points < upper_comp_ind) && (cropped_src->n_points > lower_comp_ind)){

    point_iter = chunk_ind * comps_per_chunk;
    //If there are no GAUSSIANs we only have POINT
    if (cropped_src->n_gauss == 0) {
      comp_case = P;
      temp_cropped_src->n_points = cropped_src->n_points % comps_per_chunk;
    }
    //Otherwise we have both POINT and GAUSS
    else {
      comp_case = PG;
      //Number of points left over must be modulus chunking size as we've only
      //chunked by POINTs before this point
      temp_cropped_src->n_points = cropped_src->n_points % comps_per_chunk;
      //This is the first time using GAUSS so no iterating the pointers
      gauss_iter = 0;
      //If the current number of POINTs in the chunk, plus all GAUSS are smaller
      //than a chunk size, simulate all GAUSS in this chunk
      if (temp_cropped_src->n_points + cropped_src->n_gauss < comps_per_chunk) {
        temp_cropped_src->n_gauss = cropped_src->n_gauss;
      }
      //Otherwise we need to work out how many GAUSS components will fit into
      //this current chunk, given the number of POINTs already present
      else {
        temp_cropped_src->n_gauss = comps_per_chunk - temp_cropped_src->n_points;
      }
    }//END if we need both POINT and GAUSS in this chunk
  }//END if there are POINT sources in this chunk
  //If we've gotten here, there are no POINT sources in this chunk, so work out
  //how many GAUSSIANs to add
  else {
    comp_case = G;
    //Need to work out if we used any GAUSS during a combined POINT and GAUSS
    //chunk to correctly iterate the GAUSS pointers
    //We can just use which chunk we're on and how many point sources there
    //were to work this out
    int gauss_remainder = cropped_src->n_points + cropped_src->n_gauss - chunk_ind*comps_per_chunk;
    gauss_iter = cropped_src->n_gauss - gauss_remainder;

    if (gauss_remainder > comps_per_chunk) {
      temp_cropped_src->n_gauss = comps_per_chunk;
    } else {
      temp_cropped_src->n_gauss = gauss_remainder;
    }
  }


  //FINISH work out what component type and how many of each component
  //types we need to shove into the temp_cropped_src

  //Now we know what kind of componen we have, fill the temp_cropped_src
  //Make sure to NULL out COMPONENTs that aren't present in case they are
  //hanging around in memory like sneaky ninjas
  //If just simulating POINT in this chunk
  if (comp_case == P) {
    //Null out the gauss and shapelet arrays in temp_cropped_src
    null_gauss_comps(temp_cropped_src);
    null_shapelet_comps(temp_cropped_src);
    //Increment the pointers to the correct indexes
    increment_point(temp_cropped_src, cropped_src,
                    &point_iter, woden_settings->num_time_steps);
  }
  else if (comp_case == PG) {
    null_shapelet_comps(temp_cropped_src);
    increment_point(temp_cropped_src, cropped_src,
                    &point_iter, woden_settings->num_time_steps);

    increment_gauss(temp_cropped_src, cropped_src,
                    &gauss_iter, woden_settings->num_time_steps);
  }
  else if (comp_case == G) {
    null_point_comps(temp_cropped_src);
    null_shapelet_comps(temp_cropped_src);
    increment_gauss(temp_cropped_src, cropped_src,
                    &gauss_iter, woden_settings->num_time_steps);
  }
  else {
    printf("Chunking failed to group POINT,GAUSSIAN components sensibly. Something terrible has happened\n");
  }

  temp_cropped_src->num_point_primarybeam_values = temp_cropped_src->n_points*woden_settings->num_freqs*woden_settings->num_time_steps;
  temp_cropped_src->num_gauss_primarybeam_values = temp_cropped_src->n_gauss*woden_settings->num_freqs*woden_settings->num_time_steps;
}

void fill_chunk_src_with_shapelets(catsource_t *temp_cropped_src,
     catsource_t *cropped_src, int chunk_ind, int coeffs_per_chunk,
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

  null_point_comps(temp_cropped_src);
  null_gauss_comps(temp_cropped_src);
  increment_shapelet(temp_cropped_src, cropped_src,
                    &shape_iter, woden_settings->num_time_steps);
  temp_cropped_src->num_shape_primarybeam_values = temp_cropped_src->n_shapes*woden_settings->num_freqs*woden_settings->num_time_steps;
}


source_catalogue_t * create_chunked_sky_models(catsource_t *cropped_src,
                                               woden_settings_t *woden_settings) {

  //Ok, so we split POINT and GAUSSIANs up by whole COMPONENT into the chunked
  //sky models. We split SHAPELETs up by their basis function information,
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


  int comps_per_chunk = (int)floorf((float)chunking_size / (float)(num_baselines * num_freqs * num_time_steps));

  if (comps_per_chunk < 1) {
    comps_per_chunk = 1;
  }
  //Numbers of things we are splitting up
  int num_comps_to_chunk = cropped_src->n_points + cropped_src->n_gauss;
  int num_coeffs_to_chunk = cropped_src->n_shape_coeffs;

  //How many chunks we need for the POINT/GAUSS, and the SHAPELETs
  int num_comp_chunks = (int)ceilf((float)num_comps_to_chunk / (float)comps_per_chunk);
  int num_coeff_chunks = (int)ceilf((float)num_coeffs_to_chunk / (float)comps_per_chunk);

  int num_chunks = num_comp_chunks + num_coeff_chunks;

  source_catalogue_t *chunked_sky_models;
  chunked_sky_models = malloc(sizeof(source_catalogue_t));

  chunked_sky_models->num_sources = num_chunks;
  chunked_sky_models->num_shapelets = 0;
  chunked_sky_models->catsources = malloc(num_chunks*sizeof(catsource_t));

  catsource_t *temp_cropped_src = malloc(sizeof(catsource_t));

  printf("Number of chunks required is %d\n", num_chunks);
  printf("Chunking sky model.. ");
  //Chunk the POINT/GAUSS first
  for (int comp_chunk = 0; comp_chunk < num_comp_chunks; comp_chunk++) {

    fill_chunk_src_with_pointgauss(temp_cropped_src, cropped_src, comp_chunk,
                               comps_per_chunk, woden_settings);

    chunked_sky_models->catsources[comp_chunk] = *temp_cropped_src;

  }

  for (int coeff_chunk = 0; coeff_chunk < num_coeff_chunks; coeff_chunk++) {

    fill_chunk_src_with_shapelets(temp_cropped_src, cropped_src, coeff_chunk,
                              comps_per_chunk, woden_settings);

    //Add the number of shapelets onto the full source catalogue value
    //so we know if we need to setup shapelet basis functions in GPU memory
    //or not
    chunked_sky_models->num_shapelets += temp_cropped_src->n_shapes;
    chunked_sky_models->catsources[num_comp_chunks + coeff_chunk] = *temp_cropped_src;

  }

  printf("Sky model chunked.\n");
  // free(temp_cropped_src);
  return chunked_sky_models;
}
