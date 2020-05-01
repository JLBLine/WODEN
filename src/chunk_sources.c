#include <stdio.h>
#include "read_and_write.h"

//Switch through different POINT, GAUSSIAN, SHAPELET simulation states, where
//we can either be simulating one type of component, or a combination
enum component_case {P=0, G, S, PG, PS, GS, PGS};

void null_point_comps(catsource_t *temp_cropped_src){
  temp_cropped_src->point_ras = NULL;
  temp_cropped_src->point_decs = NULL;
  temp_cropped_src->point_fluxes = NULL;
  temp_cropped_src->point_freqs = NULL;
  temp_cropped_src->point_azs = NULL;
  temp_cropped_src->point_zas = NULL;

  temp_cropped_src->n_points = 0;

}

void null_gauss_comps(catsource_t *temp_cropped_src){
  temp_cropped_src->gauss_ras = NULL;
  temp_cropped_src->gauss_decs = NULL;
  temp_cropped_src->gauss_fluxes = NULL;
  temp_cropped_src->gauss_freqs = NULL;
  temp_cropped_src->gauss_majors = NULL;
  temp_cropped_src->gauss_minors = NULL;
  temp_cropped_src->gauss_pas = NULL;
  temp_cropped_src->gauss_azs = NULL;
  temp_cropped_src->gauss_zas = NULL;

  temp_cropped_src->n_gauss = 0;

}

void null_shapelet_comps(catsource_t *temp_cropped_src){
  temp_cropped_src->shape_ras = NULL;
  temp_cropped_src->shape_decs = NULL;
  temp_cropped_src->shape_fluxes = NULL;
  temp_cropped_src->shape_freqs = NULL;
  temp_cropped_src->shape_majors = NULL;
  temp_cropped_src->shape_minors = NULL;
  temp_cropped_src->shape_pas = NULL;
  temp_cropped_src->shape_n1s = NULL;
  temp_cropped_src->shape_n2s = NULL;
  temp_cropped_src->shape_coeffs = NULL;
  temp_cropped_src->shape_param_indexes = NULL;
  temp_cropped_src->shape_azs = NULL;
  temp_cropped_src->shape_zas= NULL;

  temp_cropped_src->n_shapes = 0;
  temp_cropped_src->n_shape_coeffs = 0;

}

void increment_point(catsource_t *temp_cropped_src, catsource_t *cropped_src,
     int chunk, int chunking_size, int num_time_steps){
  //increment the required pointers to point at the beginning of the next chunk
  temp_cropped_src->point_ras = cropped_src->point_ras + (chunk * chunking_size);
  temp_cropped_src->point_decs = cropped_src->point_decs + (chunk * chunking_size);
  temp_cropped_src->point_fluxes = cropped_src->point_fluxes + (chunk * chunking_size);
  temp_cropped_src->point_freqs = cropped_src->point_freqs + (chunk * chunking_size);
  temp_cropped_src->point_azs = cropped_src->point_azs + (num_time_steps * chunk * chunking_size);
  temp_cropped_src->point_zas = cropped_src->point_zas + (num_time_steps * chunk * chunking_size);
}

void increment_gauss(catsource_t *temp_cropped_src, catsource_t *cropped_src,
     int chunk, int chunking_size, int gauss_iter, int num_time_steps){
  //increment the required pointers to point at the beginning of the next chunk
  temp_cropped_src->gauss_ras = cropped_src->gauss_ras + gauss_iter;
  temp_cropped_src->gauss_decs = cropped_src->gauss_decs + gauss_iter;
  temp_cropped_src->gauss_fluxes = cropped_src->gauss_fluxes + gauss_iter;
  temp_cropped_src->gauss_freqs = cropped_src->gauss_freqs + gauss_iter;

  temp_cropped_src->gauss_majors = cropped_src->gauss_majors + gauss_iter;
  temp_cropped_src->gauss_minors = cropped_src->gauss_minors + gauss_iter;
  temp_cropped_src->gauss_pas = cropped_src->gauss_pas + gauss_iter;

  temp_cropped_src->gauss_azs = cropped_src->gauss_azs + (num_time_steps * chunk * chunking_size);
  temp_cropped_src->gauss_zas = cropped_src->gauss_zas + (num_time_steps * chunk * chunking_size);
}

void increment_shapelet(catsource_t *temp_cropped_src, catsource_t *cropped_src,
     int chunk, int chunking_size, int shape_iter, int num_time_steps){
  //increment the required pointers to point at the beginning of the next chunk
  temp_cropped_src->shape_ras = cropped_src->shape_ras;
  temp_cropped_src->shape_decs = cropped_src->shape_decs;
  temp_cropped_src->shape_fluxes = cropped_src->shape_fluxes;
  temp_cropped_src->shape_freqs = cropped_src->shape_freqs;

  temp_cropped_src->shape_majors = cropped_src->shape_majors;
  temp_cropped_src->shape_minors = cropped_src->shape_minors;
  temp_cropped_src->shape_pas = cropped_src->shape_pas;

  temp_cropped_src->shape_coeffs = cropped_src->shape_coeffs + shape_iter;
  temp_cropped_src->shape_n1s = cropped_src->shape_n1s + shape_iter;
  temp_cropped_src->shape_n2s = cropped_src->shape_n2s + shape_iter;
  temp_cropped_src->shape_param_indexes = cropped_src->shape_param_indexes + shape_iter;

  temp_cropped_src->shape_azs = cropped_src->shape_azs + (num_time_steps * chunk * chunking_size);
  temp_cropped_src->shape_zas = cropped_src->shape_zas + (num_time_steps * chunk * chunking_size);
}

void fill_chunk_src(catsource_t *temp_cropped_src, catsource_t *cropped_src,
     int num_chunks, int chunk, int chunking_size, int num_time_steps ){

  //So want to chunk in order point, gaussian, shapelet
  //First we have to ascertain where in the chunking we are, and which type
  //of component we have to include

  //Lower and upper indexes of components covered in this chunk
  int lower_comp_ind = chunk * chunking_size;
  int upper_comp_ind = (chunk + 1) * chunking_size - 1;

  //comp_case is used to ascertain what combo of POINT,GAUSSIAN,SHAPELET we have
  int comp_case = -1;

  //Used to iterate the pointers of the GAUSS and SHAPELET parameters in the
  //temp_cropped_src. As POINT is done first, we don't need to make a record,
  //we can just use maths later on
  int gauss_iter = 0;
  int shape_iter = 0;

  //BEGIN work out what component type and how many of each component
  //types we need to shove into the temp_cropped_src

  //If chunk lies within the number of POINT components
  if (cropped_src->n_points > upper_comp_ind){
    comp_case = P;
    temp_cropped_src->n_points = chunking_size;
  }
  //If chunk contains POINT components, but extends beyond the number POINT components
  else if ((cropped_src->n_points < upper_comp_ind) && (cropped_src->n_points > lower_comp_ind)){
    //If there are no GAUSSIAN or SHAPELET, we only have POINT
    if ( (cropped_src->n_gauss == 0) && (cropped_src->n_shape_coeffs == 0) ) {
      comp_case = P;
      temp_cropped_src->n_points = cropped_src->n_points % chunking_size;
    }
    //If there are no SHAPELET sources, we need both POINT and GAUSS
    else if (cropped_src->n_shape_coeffs == 0) {
      comp_case = PG;
      temp_cropped_src->n_points = cropped_src->n_points % chunking_size;
      //If the current number of POINTs in the chunk, plus all GAUSS are smaller
      //than a chunk size, simulate all GAUSS in this chunk
      if (temp_cropped_src->n_points + cropped_src->n_gauss < chunking_size) {
        temp_cropped_src->n_gauss = cropped_src->n_gauss;
        //This is the first time using GAUSS so no iterating the pointers
        gauss_iter = 0;
      }
      //Otherwise we need to work out how many GAUSS components will fit into
      //this current chunk, given the number of POINTs already present
      else {
        temp_cropped_src->n_gauss = chunking_size - temp_cropped_src->n_points;
        //This is the first time using GAUSS so no iterating the pointers
        gauss_iter = 0;
      }//END if no SHAPELET components
    }//END if there are no SHAPELET sources, we need both POINT and GAUSS

    //If there are no GAUSS sources, we need both POINT and SHAPELET
    else if (cropped_src->n_gauss == 0) {
      comp_case = PS;
      temp_cropped_src->n_points = cropped_src->n_points % chunking_size;
      //If the current number of POINTs in the chunk, plus all SHAPELET fit
      //within a chunk size, simulate all SHAPELET in this chunk
      if (temp_cropped_src->n_points + cropped_src->n_gauss <= chunking_size) {
        temp_cropped_src->n_shapes = cropped_src->n_shapes;
        temp_cropped_src->n_shape_coeffs = cropped_src->n_shape_coeffs;
        //This is the first time using SHAPELET so no iterating the pointers
        shape_iter = 0;
      }
      //Otherwise we need to work out how many SHAPELET components will fit into
      //this current chunk, given the number of POINTs already present
      else {
        temp_cropped_src->n_shapes = cropped_src->n_shapes;
        temp_cropped_src->n_shape_coeffs = chunking_size - temp_cropped_src->n_points;
        //This is the first time using SHAPELET so no iterating the pointers
        shape_iter = 0;
      }
    }//END if there are no GAUSS sources, we need both POINT and SHAPELET

    //If we've gotten here, there are POINT sources in this chunk, and there
    //are both GAUSSIAN and SHAPELETS to fill the rest of chunk
    else {
      temp_cropped_src->n_points = cropped_src->n_points % chunking_size;

      //If we can fill the rest of the chunk with GAUSSIANs
      if (temp_cropped_src->n_points + cropped_src->n_gauss >= chunking_size) {
        comp_case = PG;
        temp_cropped_src->n_gauss = chunking_size - temp_cropped_src->n_points;
        //This is the first time using GAUSS so no iterating the pointers
        gauss_iter = 0;
      }
      //else we can't fill the chunk size with POINTS and GAUSS, need to add in
      //some SHAPELETs
      else {
        comp_case = PGS;
        temp_cropped_src->n_gauss = cropped_src->n_gauss;

        temp_cropped_src->n_shapes = cropped_src->n_shapes;
        temp_cropped_src->n_shape_coeffs = chunking_size - temp_cropped_src->n_points - temp_cropped_src->n_gauss;
        //This is the first time using SHAPELET so no iterating the pointers
        shape_iter = 0;
      }

    }//END else there are POINT sources in this chunk, and there are both GAUSSIAN and SHAPELETS to fill the rest of chunk
  }//END if chunk contains POINT components, but extends beyond the number POINT components

  //We have established there are no POINTs in this chunk
  //Here, if GAUSS extend beyond chunk, we're only simulating GAUSSIANS
  else if (cropped_src->n_gauss + cropped_src->n_points >= upper_comp_ind){
    comp_case = G;
    temp_cropped_src->n_gauss = chunking_size;

    //Need to work out if we used any GAUSS during a combined POINT and GAUSS
    //chunk to correctly iterate the GAUSS pointers
    int gauss_remainder = cropped_src->n_points + cropped_src->n_gauss - chunk*chunking_size;
    gauss_iter = cropped_src->n_gauss - gauss_remainder;
  }//END if GAUSS extend beyond chunk, we're only simulating GAUSSIANS

  //Here, there are no POINTs in the chunk, and not enough GAUSS to fill the chunk
  else if ((cropped_src->n_gauss + cropped_src->n_points < upper_comp_ind) && (cropped_src->n_gauss + cropped_src->n_points >= lower_comp_ind) ){
    int gauss_remainder = cropped_src->n_points + cropped_src->n_gauss - chunk*chunking_size;
    gauss_iter = cropped_src->n_gauss - gauss_remainder;
    temp_cropped_src->n_gauss = gauss_remainder;

    //Here there are no SHAPELET coeffs, so we're just doing GAUSSIAN
    if (cropped_src->n_shape_coeffs == 0) {
      comp_case = G;
    }
    else {
      comp_case = GS;
      temp_cropped_src->n_shapes = cropped_src->n_shapes;
      int shape_remainder = (chunk + 1)*chunking_size - cropped_src->n_points - cropped_src->n_gauss;
      temp_cropped_src->n_shape_coeffs = shape_remainder;
      //This is the first time using SHAPELET so no iterating the pointers
      shape_iter = 0;
    }
  }

  //if we get here, we should only have SHAPELETs left to fill the chunk with
  else {
    comp_case = S;
    //These are the number of SHAPELET coefficients yet to be simulated
    int shape_remainder = cropped_src->n_points + cropped_src->n_gauss + cropped_src->n_shape_coeffs - chunk*chunking_size;

    if (shape_remainder > chunking_size) {
      temp_cropped_src->n_shape_coeffs = chunking_size;
    }
    else {
      temp_cropped_src->n_shape_coeffs = shape_remainder;
    }

    shape_iter = cropped_src->n_shape_coeffs - shape_remainder;
  }

  //FINISH work out what component type and how many of each component
  //types we need to shove into the temp_cropped_src

  //Now we know what kind of componen we have, fill the temp_cropped_src
  //If just simulating POINT in this chunk
  if (comp_case == P) {
    //Null out the gauss and shapelet arrays in temp_cropped_src
    null_gauss_comps(temp_cropped_src);
    null_shapelet_comps(temp_cropped_src);
    //Increment the pointers to the correct indexes
    increment_point(temp_cropped_src, cropped_src, chunk, chunking_size, num_time_steps);
    // printf("COMP CASE is just doing P\n");
  }

  else if (comp_case == PG) {
    null_shapelet_comps(temp_cropped_src);
    increment_point(temp_cropped_src, cropped_src, chunk, chunking_size, num_time_steps);
    increment_gauss(temp_cropped_src, cropped_src,
         chunk, chunking_size, gauss_iter, num_time_steps);
    // printf("COMP CASE is doing PG, gauss_iter is %d\n",gauss_iter);
  }

  else if (comp_case == PS) {
    null_gauss_comps(temp_cropped_src);
    increment_point(temp_cropped_src, cropped_src, chunk, chunking_size, num_time_steps);
    increment_shapelet(temp_cropped_src, cropped_src,
         chunk, chunking_size, shape_iter, num_time_steps);
    // printf("COMP CASE is doing PS, shape_iter is %d\n",shape_iter);
  }

  else if (comp_case == G) {
    null_point_comps(temp_cropped_src);
    null_shapelet_comps(temp_cropped_src);
    increment_gauss(temp_cropped_src, cropped_src,
         chunk, chunking_size, gauss_iter, num_time_steps);
    // printf("COMP CASE is just doing G, gauss_iter is %d\n",gauss_iter);
  }

  else if (comp_case == GS) {
    null_point_comps(temp_cropped_src);
    increment_gauss(temp_cropped_src, cropped_src,
         chunk, chunking_size, gauss_iter, num_time_steps);
    increment_shapelet(temp_cropped_src, cropped_src,
         chunk, chunking_size, shape_iter, num_time_steps);
    // printf("COMP CASE is doing GS, gauss_iter %d, shape_iter %d\n", gauss_iter, shape_iter);
  }

  else if (comp_case == S) {
    null_point_comps(temp_cropped_src);
    null_gauss_comps(temp_cropped_src);
    increment_shapelet(temp_cropped_src, cropped_src,
         chunk, chunking_size, shape_iter, num_time_steps);
    // printf("COMP CASE is just doing S, shape_iter is %
  }

  else if (comp_case == PGS) {
    increment_point(temp_cropped_src, cropped_src, chunk, chunking_size, num_time_steps);
    increment_gauss(temp_cropped_src, cropped_src,
         chunk, chunking_size, gauss_iter, num_time_steps);
    increment_shapelet(temp_cropped_src, cropped_src,
         chunk, chunking_size, shape_iter, num_time_steps);
  }

  else {
    printf("Chunking failed to group POINT,GAUSSIAN,SHAPELET components sensibly. Something terrible has happened\n");
  }

}
