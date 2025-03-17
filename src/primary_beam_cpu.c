#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include "constants.h"
#include "fundamental_coords_cpu.h"
#include "primary_beam_cpu.h"

void twoD_Gaussian_cpu(user_precision_t x, user_precision_t y,
           user_precision_t xo, user_precision_t yo,
           user_precision_t sigma_x, user_precision_t sigma_y,
           user_precision_t cos_theta, user_precision_t sin_theta,
           user_precision_t sin_2theta,
           user_precision_t * beam_real, user_precision_t * beam_imag) {

  user_precision_t a, b, c;

  a = (cos_theta*cos_theta)/(2*sigma_x*sigma_x) + (sin_theta*sin_theta)/(2*sigma_y*sigma_y);
  b = -sin_2theta / (4*sigma_x*sigma_x) + sin_2theta / (4*sigma_y*sigma_y);
  c = (sin_theta*sin_theta)/(2*sigma_x*sigma_x) + (cos_theta*cos_theta)/(2*sigma_y*sigma_y);

  * beam_real = exp( -( a*(x-xo)*(x-xo) + 2*b*(x-xo)*(y-yo) + c*(y-yo)*(y-yo) ));
  * beam_imag = 0.0;

}

void gaussian_beam_from_lm_cpu(double *beam_ls, double *beam_ms,
                      double beam_ref_freq, double *freqs,
                      user_precision_t fwhm_lm, user_precision_t cos_theta,
                      user_precision_t sin_theta, user_precision_t sin_2theta,
                      int num_components, int num_time_steps, int num_freqs,
                      user_precision_complex_t *g1xs,
                      user_precision_complex_t *g1ys) {

  int component, time_ind, beam_ind;
  user_precision_t d_beam_real, d_beam_imag, std;
  user_precision_complex_t temp;
  
  for (int iFreq = 0; iFreq < num_freqs; iFreq++) {
    for (int iLMcoord = 0; iLMcoord < num_components * num_time_steps; iLMcoord++) {

      component = (int)floorf((float)iLMcoord / (float)num_time_steps);
      time_ind = iLMcoord - component*num_time_steps;
      beam_ind = num_freqs*time_ind*num_components + (num_components*iFreq) + component;

      //Convert FWHM into standard dev, and scale for frequency
      std = (fwhm_lm / FWHM_FACTOR) * (user_precision_t)(beam_ref_freq / freqs[iFreq]);

      twoD_Gaussian_cpu((user_precision_t)beam_ls[iLMcoord], (user_precision_t)beam_ms[iLMcoord], 0, 0,
                std, std, cos_theta, sin_theta, sin_2theta,
                &d_beam_real, &d_beam_imag);
      
      temp = d_beam_real + I*d_beam_imag;

      g1xs[beam_ind] = temp;
      g1ys[beam_ind] = temp;
    }
  }
}

void calculate_gaussian_beam_cpu(int num_components, int num_time_steps,
           int num_freqs, user_precision_t ha0,
           user_precision_t sdec0, user_precision_t cdec0,
           user_precision_t fwhm_lm, user_precision_t cos_theta,
           user_precision_t sin_theta, user_precision_t sin_2theta,
           double beam_ref_freq, double *freqs,
           double *beam_has, double *beam_decs,
           user_precision_complex_t *g1xs,
           user_precision_complex_t *g1ys){

  int num_beam_hadec = num_components * num_time_steps;

  double *beam_ls = malloc(num_beam_hadec*sizeof(double));
  double *beam_ms = malloc(num_beam_hadec*sizeof(double));
  double *beam_ns = malloc(num_beam_hadec*sizeof(double));

  calc_lmn_cpu(ha0, sdec0, cdec0, beam_has, beam_decs,
               beam_ls, beam_ms, beam_ns, num_beam_hadec);

  gaussian_beam_from_lm_cpu(beam_ls, beam_ms,
                      beam_ref_freq, freqs,
                      fwhm_lm, cos_theta,
                      sin_theta, sin_2theta,
                      num_components, num_time_steps, num_freqs,
                      g1xs, g1ys);
  free( beam_ns );
  free( beam_ms );
  free( beam_ls );
}

void analytic_dipole_cpu(user_precision_t az, user_precision_t za,
           user_precision_t wavelength,
           user_precision_complex_t * beam_X,
           user_precision_complex_t * beam_Y) {

  user_precision_t dipole_height_m = 0.3;

  //Here, X means a north-south aligned dipole
  //      Y means an east-west aligned dipole

  user_precision_t theta_parallel_X = acos(sin(za)*cos(az));
  user_precision_t theta_parallel_Y = acos(sin(za)*sin(az));

  user_precision_t d_in_lambda = (2. * dipole_height_m)/wavelength;
  user_precision_t gp_effect_array = 2. * sin(M_PI*d_in_lambda*cos(za));

  user_precision_t voltage_parallel_X = sin(theta_parallel_X) * gp_effect_array;
  user_precision_t voltage_parallel_Y = sin(theta_parallel_Y) * gp_effect_array;

  // user_precision_complex_t tempX;
  // user_precision_complex_t tempY;

  // tempX = voltage_parallel_X + I*0;
  // tempY = voltage_parallel_Y + I*0;

  * beam_X = voltage_parallel_X + I*0;
  * beam_Y = voltage_parallel_Y + I*0;

}

void calculate_analytic_dipole_beam_cpu(int num_components,
     int num_time_steps, int num_freqs,
     user_precision_t *azs, user_precision_t *zas, double *freqs,
     user_precision_complex_t *g1xs, user_precision_complex_t *g1ys){

  int num_beam_azza = num_components * num_time_steps;

  user_precision_complex_t beam_norm_X, beam_norm_Y, beam_X, beam_Y, normed_X, normed_Y;
  user_precision_t wavelength;
  int component, time_ind, beam_ind;

  for (int iFreq = 0; iFreq < num_freqs; iFreq++) {
    
    wavelength = VELC / freqs[iFreq];
    analytic_dipole_cpu(0.0, 0.0, wavelength,
               &beam_norm_X, &beam_norm_Y);

    for (int iCoord = 0; iCoord < num_beam_azza; iCoord++) {

      component = (int)floorf((float)iCoord / (float)num_time_steps);
      time_ind = iCoord - component*num_time_steps;
      beam_ind = num_freqs*time_ind*num_components + (num_components*iFreq) + component;

      analytic_dipole_cpu(azs[iCoord], zas[iCoord], wavelength,
                      &beam_X, &beam_Y);

      //Analytic beam is entirely real, so can just normalise by real values
      normed_X = creal(beam_X) / creal(beam_norm_X) + I*0;
      normed_Y = creal(beam_Y) / creal(beam_norm_Y) + I*0;

      g1xs[beam_ind] = normed_X;
      g1ys[beam_ind] = normed_Y;

    }
  }
}

void RTS_MWA_beam_cpu(user_precision_t az, user_precision_t za,
           double ha, double dec,
           double wavelength, double *metre_delays,
           double latitude, int norm,
           user_precision_complex_t * gx, user_precision_complex_t * Dx,
           user_precision_complex_t * Dy, user_precision_complex_t * gy) {

  // set elements of the look-dir vector
  double proj_e = sin(za)*sin(az);
  double proj_n = sin(za)*cos(az);
  double proj_z = cos(za);

  int n_cols = 4;
  int n_rows = 4;

  //Used in calculating the phase later on
  double multiplier = -2 * M_PI / wavelength;
  double dipl_e, dipl_n, dipl_z;

  user_precision_complex_t x_dip, y_dip;

  user_precision_complex_t gx_dip = 0.0 + I*0;
  user_precision_complex_t Dx_dip = 0.0 + I*0;
  user_precision_complex_t Dy_dip = 0.0 + I*0;
  user_precision_complex_t gy_dip = 0.0 + I*0;

  user_precision_t gainx = 1.0;
  user_precision_t gainy = 1.0;

  int k = 0;

  for (int i = 0; i < n_cols; i++) {
    for (int j = 0; j < n_rows; j++) {

      // set elements of the baseline vector
      dipl_e = (i - 1.5) * MWA_DIPOLE_SEP;
      dipl_n = (j - 1.5) * MWA_DIPOLE_SEP;
      dipl_z = 0.0;

      double phase = multiplier*(dipl_e*proj_e + dipl_n*proj_n + dipl_z*proj_z - metre_delays[k]);
      double phaseshift_re = cos(phase);
      double phaseshift_im = sin(phase);

      //TODO You could get gains here from some dipole flagging scheme in the future
      //For now, we've just set them to one above

      x_dip = gainx*phaseshift_re + I*gainx*phaseshift_im;
      y_dip = gainy*phaseshift_re + I*gainy*phaseshift_im;

      gx_dip += x_dip;
      Dx_dip += x_dip;
      Dy_dip += y_dip;
      gy_dip += y_dip;

      k += 1;
    }
  }

  //Calculate the effect of the ground plane
  user_precision_t ground_plane = 2.0*sin(2.0*M_PI*MWA_DIPOLE_HEIGHT/wavelength*cos(za));

  //Normalise the beam if requested
  if (norm == 1){
    ground_plane /= 2.0*sin(2.0*M_PI*MWA_DIPOLE_HEIGHT/wavelength);
  }

  //Used in some kind of parallatic rotation?
  double coslat = cos(latitude);
  double cosdec = cos(dec);
  double cosha = cos(ha);
  double sinlat = sin(latitude);
  double sindec = sin(dec);
  double sinha = sin(ha);

  //Some kind of parallatic rotation?
  user_precision_t rot0 = coslat*cosdec + sinlat*sindec*cosha;
  user_precision_t rot1 = -sinlat*sinha;
  user_precision_t rot2 = sindec*sinha;
  user_precision_t rot3 = cosha;

  //Normalise the ground plane to the number of dipoles??
  user_precision_t ground_plane_div_dipoles = ground_plane / NUM_DIPOLES;

  user_precision_complex_t pgx = gx_dip * rot0 * ground_plane_div_dipoles;
  user_precision_complex_t pDx = Dx_dip * rot1 * ground_plane_div_dipoles;
  user_precision_complex_t pDy = Dy_dip * rot2 * ground_plane_div_dipoles;
  user_precision_complex_t pgy = gy_dip * rot3 * ground_plane_div_dipoles;
  

  //   pgx = creal(pgx) + I*0;
  // pDx = creal(pDx) + I*0;
  // pDy = creal(pDy) + I*0;
  // pgy = creal(pgy) + I*0;

  // //Explicitly set the imag to zero, this beam is real only
  // * gx = pgx;
  // * Dx = pDx;
  // * Dy = pDy;
  // * gy = pgy;

  // printf("gx %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f \n",
  //        creal(pgx), cimag(pgx), creal(pDx), cimag(pDx),
  //        creal(pDy), cimag(pDy), creal(pgy), cimag(pgy));

  //Explicitly set the imag to zero, this beam is real only
  * gx = creal(pgx) + I*0;
  * Dx = creal(pDx) + I*0;
  * Dy = creal(pDy) + I*0;
  * gy = creal(pgy) + I*0;

}

void calculate_RTS_MWA_analytic_beam_cpu(int num_components,
     int num_time_steps, int num_freqs,
     user_precision_t *azs, user_precision_t *zas, int *delays,
     double latitude, int norm,
     double *beam_has, double *beam_decs, double *freqs,
     user_precision_complex_t *gxs, user_precision_complex_t *Dxs,
     user_precision_complex_t *Dys, user_precision_complex_t *gys){

  int num_coords = num_components * num_time_steps;
  int component, time_ind, beam_ind;
  user_precision_complex_t gx, Dx, Dy, gy;
  double wavelength;

  //Apply the actual delay length added by MWA circuitry here, in metres - saves
  //on computation inside the kernel
  //Make a copy array so we don't modify the original
  //Have to reorder the delays as listed in the metafits to match what
  //the RTS analytic code
  double *metre_delays = malloc(NUM_DIPOLES*sizeof(double));

  for (int i = 0; i < 4; i++) {

      metre_delays[3-i] = (double)delays[0+i*4];
      metre_delays[7-i] = (double)delays[1+i*4];
      metre_delays[11-i] = (double)delays[2+i*4];
      metre_delays[15-i] = (double)delays[3+i*4];

  }

  //I have NO IDEA what this is doing, blindy copying the RTS code
  //One change to the RTS code is I take out the division by speed
  //of light, as later on the delays are multipled by speed of light again
  double delay_0 = 0.0;
  for (int k=0; k<NUM_DIPOLES; k++ ){
    delay_0 += metre_delays[k] * DQ;
  }

  delay_0 /= (double)NUM_DIPOLES;

  for(int k=0; k<NUM_DIPOLES; k++) {
    metre_delays[k] = metre_delays[k] * DQ - delay_0;
  }

  for (int iFreq = 0; iFreq < num_freqs; iFreq++) {
    for (int iCoord = 0; iCoord < num_coords; iCoord++) {

      component = (int)floorf((float)iCoord / (float)num_time_steps);
      time_ind = iCoord - component*num_time_steps;
      beam_ind = num_freqs*time_ind*num_components + (num_components*iFreq) + component;

      wavelength = VELC / freqs[iFreq];

      RTS_MWA_beam_cpu(azs[iCoord], zas[iCoord],
              beam_has[iCoord], beam_decs[iCoord],
              wavelength, metre_delays,
              latitude, norm,
              &gx, &Dx, &Dy, &gy);

      gxs[beam_ind] = gx;
      Dxs[beam_ind] = Dx;
      Dys[beam_ind] = Dy;
      gys[beam_ind] = gy;

    }
  }
  free(metre_delays);
}


// // we have 4 dimensions to loop over here:
// //  - component
// //  - freqs
// //  - times
// //  - tiles
// //
// // it's likely that time will usually be the smallest, so I'm going to
// // loop over that. Time will tell whether that's a good idea or not LOL
// __global__ void kern_map_hyperbeam_gains(int num_components,
//            int num_times, int num_freqs, int num_tiles, int iTime, int num_unique_fee_freqs,
//            double *d_jones, const int *d_tile_map, const int *d_freq_map,
//            int parallactic,
//            user_precision_complex_t *d_gxs,
//            user_precision_complex_t *d_Dxs,
//            user_precision_complex_t *d_Dys,
//            user_precision_complex_t *d_gys) {

//   //All baselines at all freqs and all times
//   int iComponent = threadIdx.x + (blockDim.x*blockIdx.x);

//   //The tile to map - currently unused but hopefully in the future tis possible
//   int iTile = threadIdx.y + (blockDim.y*blockIdx.y);

//   //freq step
//   int iFreq = threadIdx.z + (blockDim.z*blockIdx.z);

//   if(iComponent < num_components && iTime < num_times && iFreq < num_freqs && iTile < num_tiles) {

//     user_precision_complex_t d_beam_J00;
//     user_precision_complex_t d_beam_J01;
//     user_precision_complex_t d_beam_J10;
//     user_precision_complex_t d_beam_J11;

//     // For *this tile* and *this frequency*, access the de-duplicated beam
//     // response.
//     int i_row = d_tile_map[iTile];
//     int i_col = d_freq_map[iFreq];

//     int hyper_ind, current_ind;

//     //If we are doing parallactic rotation corrections, the latitude needs
//     //to be fed in. As this change change with time, we call hyperbeam for
//     //each time step, and use pointer arithmatic when storing the outputs
//     //This means the striping for retrieving outputs is different.
//     if (parallactic == 1)
//     {
//       //Where the desired data sits within a single time step hyperdrive
//       //output
//       hyper_ind = ((num_components*num_unique_fee_freqs*i_row) + num_components * i_col);

//       //Add on how many previous time-steps-worth of outputs we've
//       //already remapped
//       current_ind = hyper_ind + iComponent;
//     } else {
//       hyper_ind = ((num_components*num_times*num_unique_fee_freqs*i_row) + num_components * i_col);

//       //There is no frequency used in first chunk of this remapping, as that's handled
//       //by the hyper_ind
//       current_ind = iTime*num_components*num_tiles + hyper_ind + iComponent;
//     }

//     d_beam_J00.x = (user_precision_t)d_jones[2*MAX_POLS*current_ind + 0];
//     d_beam_J00.y = (user_precision_t)d_jones[2*MAX_POLS*current_ind + 1];
//     d_beam_J01.x = (user_precision_t)d_jones[2*MAX_POLS*current_ind + 2];
//     d_beam_J01.y = (user_precision_t)d_jones[2*MAX_POLS*current_ind + 3];
//     d_beam_J10.x = (user_precision_t)d_jones[2*MAX_POLS*current_ind + 4];
//     d_beam_J10.y = (user_precision_t)d_jones[2*MAX_POLS*current_ind + 5];
//     d_beam_J11.x = (user_precision_t)d_jones[2*MAX_POLS*current_ind + 6];
//     d_beam_J11.y = (user_precision_t)d_jones[2*MAX_POLS*current_ind + 7];

//     //Get an index to split into the WODEN style containers
//     int new_ind =  num_times*num_freqs*num_components*iTile + num_freqs*iTime*num_components + (num_components*iFreq) + iComponent;

//     d_gxs[new_ind] = d_beam_J00;
//     d_Dxs[new_ind] = d_beam_J01;
//     d_Dys[new_ind] = d_beam_J10;
//     d_gys[new_ind] = d_beam_J11;

//   }
// }


// __global__ void fill_with_ones(int num_azza, double *d_jones) {

//   //All baselines at all freqs and all times
//   int iAzza = threadIdx.x + (blockDim.x*blockIdx.x);

//   if (iAzza < num_azza) {
//     d_jones[iAzza] = 1;
//   }
// }

// Helper function to check if a value already exists in the array
// Source: Adapted with assistance from ChatGPT (OpenAI), 2024.
bool is_present(double *arr, int size, double value, double epsilon) {
    for (int i = 0; i < size; i++) {
        if (fabs(arr[i] - value) < epsilon) {
            return true;
        }
    }
    return false;
}

// Function to find unique values in an array
// Source: Adapted with assistance from ChatGPT (OpenAI), 2024.
int find_unique(double *input, int n, double *output, double epsilon) {
    int unique_count = 0;
    for (int i = 0; i < n; i++) {
        if (!is_present(output, unique_count, input[i], epsilon)) {
            output[unique_count++] = input[i];
        }
    }
    return unique_count;
}


void map_hyperbeam_gains_cpu(double *jones, int num_beams, int beam_ind,
                             int num_time_steps, int time_ind,
                             int num_freqs, int unq_freq_ind,
                             int num_components, 
                             int *freq_map,
                             user_precision_complex_t *gxs,
                             user_precision_complex_t *Dxs,
                             user_precision_complex_t *Dys,
                             user_precision_complex_t *gys){

  //loop over all freq mapping indexes; if the freq_map matches the unq_freq_ind,
  //fill the outputs with the jones values

  for (int freq_map_ind = 0; freq_map_ind < num_freqs; freq_map_ind++) {
    if (freq_map[freq_map_ind] == unq_freq_ind) {
      for (int iComponent = 0; iComponent < num_components; iComponent++) {

        int gain_ind =  num_time_steps*num_freqs*num_components*beam_ind + num_freqs*time_ind*num_components + (num_components*freq_map_ind) + iComponent;

        gxs[gain_ind] = jones[2*MAX_POLS*iComponent + 0] + I*jones[2*MAX_POLS*iComponent + 1];
        Dxs[gain_ind] = jones[2*MAX_POLS*iComponent + 2] + I*jones[2*MAX_POLS*iComponent + 3];
        Dys[gain_ind] = jones[2*MAX_POLS*iComponent + 4] + I*jones[2*MAX_POLS*iComponent + 5];
        gys[gain_ind] = jones[2*MAX_POLS*iComponent + 6] + I*jones[2*MAX_POLS*iComponent + 7];
      }
    }
  }
  // current_ind = iTime*num_components*num_tiles + hyper_ind + iComponent;
  // new_ind =  num_times*num_freqs*num_components*iTile + num_freqs*iTime*num_components + (num_components*iFreq) + iComponent;
}


void run_hyperbeam_cpu(int num_components,
           int num_time_steps, int num_freqs,
           int num_beams, uint8_t parallactic,
           double *freqs, struct FEEBeam *fee_beam,
           uint32_t *hyper_delays, int num_amps, double *amps,
           double *azs, double *zas,
           double *latitudes,
           user_precision_complex_t *gxs,
           user_precision_complex_t *Dxs,
           user_precision_complex_t *Dys,
           user_precision_complex_t *gys){

  //Always do things in IAU order
  int iau_order = 1;

  //Always norm to zenith
  uint8_t norm_to_zenith = 1;

  // int num_azza = num_components * num_time_steps;
  // // int num_beam_values = num_azza * num_freqs * num_beams;
  // int num_beam_values = num_azza * num_beams;

  //First thing to do is find the unique frequencies in the input freqs
  //The GPU code has a helper function `get_num_unique_fee_freqs`, but the CPU
  //code does not. So do it ourselves here (there might be caching so we 
  //miiiiiight not be getting a compete speed up doing this, but still worth doing)

  //First, find the closest availble frequencies to the ones we want using
  //the hyperbeam `fee_closest_freq` function
  double *closest_freqs = malloc(num_freqs*sizeof(double));
  //Get the closest frequencies to the ones we want
  for (int iFreq = 0; iFreq < num_freqs; iFreq++) {
    closest_freqs[iFreq] = fee_closest_freq(fee_beam, freqs[iFreq]);
    // printf("Closest freq. to %.2f MHz: %.2f MHz\n", freqs[iFreq]/1e6, closest_freqs[iFreq]/1e6);
  }

  //Next, find all unique avaible frequencies
  double *unique_freqs = malloc(num_freqs*sizeof(double));
  int num_unique_fee_freqs = find_unique(closest_freqs, num_freqs, unique_freqs, 1);

  //Now make a map of all the frequencies we want to the unique ones
  //We can use this to map calculated gains to output arrays
  int *freq_map = malloc(num_freqs*sizeof(int));
  for (int iFreq = 0; iFreq < num_freqs; iFreq++) {
    for (int iUnique = 0; iUnique < num_unique_fee_freqs; iUnique++) {
      if (closest_freqs[iFreq] == unique_freqs[iUnique]) {
        freq_map[iFreq] = iUnique;
        // printf("We haz done the thing %d\n", freq_map[iFreq]);
        // break;
      }
    }
  }

  // for (int iFreq = 0; iFreq < num_freqs; iFreq++) {
  //   printf("Freq %.2f MHz maps to unique index %d freq %.2f\n", freqs[iFreq], freq_map[iFreq], unique_freqs[freq_map[iFreq]]);
  // }


  // printf("num_unique: %d\n", num_unique);

  //In CPU version, we'll be looping over tile, time (if parallactic rotating),
  //and frequencies, so we only need to pass enough memory to calculate all
  //directions. If not doing parallelactic rotation, we can just pass all
  //directions for all times at once (as the latitude isn't changing with time)
  //so we need more memory
  double *jones = NULL;
  if (parallactic == 1) {
    jones = malloc(2*MAX_POLS*num_components*sizeof(double));
  } else {
    jones = malloc(2*MAX_POLS*num_components*num_time_steps*sizeof(double));
  }

  int32_t status = 0;

  //To be as efficient as possible, we'd want to make something like this
  //hyperbeam GPU function, so we on'y unique combos of dipole amplitudes.
  //For now I think research applications will have unique dipole amplitudes
  //per tile, so probably not going to have any effect
  // tile_map = get_fee_tile_map(fee_beam);

  int amp_increment;
  int azza_increment;

  for (int tile_ind = 0; tile_ind < num_beams; tile_ind++) {
    // printf("HERE Tile %d, %d\n", tile_ind, parallactic);
    //Which tile we are calculating dictates which amplitudes to use
    amp_increment = tile_ind*num_amps;

    for (int unq_freq_ind = 0; unq_freq_ind < num_unique_fee_freqs; unq_freq_ind++) {
      // printf("HERE Freq %d\n", unq_freq_ind);

      if (parallactic == 1) {
        //The latitude of the array should change with time if we are precessing
        //it back to J2000. This means we have to call hyperbeam for as many time
        //steps as we have. Supply chunks of az, za, and jones via pointer
        //arithmatic
        
        for (int time_ind = 0; time_ind < num_time_steps; time_ind++) {
          // printf("HERE Time %d\n", time_ind);
          //How far through the time steps we are dictates which az,za to use
          azza_increment = time_ind*num_components;

          status = fee_calc_jones_array(fee_beam,(uint32_t)num_components,
                              azs + azza_increment, zas + azza_increment,
                              (uint32_t)unique_freqs[unq_freq_ind],
                              hyper_delays, amps + amp_increment,
                              (uint32_t)num_amps,
                              norm_to_zenith,
                              &latitudes[time_ind],
                              iau_order, jones);

          // printf("Time %d Freq %d %.3e %.3e\n", time_ind, unq_freq_ind, jones[0], jones[1]);

          map_hyperbeam_gains_cpu(jones, num_beams, tile_ind,
                             num_time_steps, time_ind,
                             num_freqs, unq_freq_ind,
                             num_components, freq_map,
                             gxs, Dxs, Dys, gys);

          
        
        }//end time loop

      } else {
        //By setting latitude to a NULL, it tells the hyperbeam code to not
        //do any parallactic rotation
        double *current_latitude = NULL;
        fee_calc_jones_array(fee_beam,(uint32_t)(num_components*num_time_steps),
                            azs, zas, (uint32_t)unique_freqs[unq_freq_ind],
                            hyper_delays, amps + amp_increment, (uint32_t)num_amps,
                            norm_to_zenith, current_latitude ,
                            iau_order, jones);
        //gain mapping function works per time step, so we need to loop over
        int gain_increment;
        for (int time_ind = 0; time_ind < num_time_steps; time_ind++) {
            //How far through the time steps we are dictates which az,za to use
            gain_increment = 2*MAX_POLS*time_ind*num_components;

            // printf("gain_increment %d\n", gain_increment);

            map_hyperbeam_gains_cpu(jones + gain_increment, num_beams, tile_ind,
                              num_time_steps, time_ind,
                              num_freqs, unq_freq_ind,
                              num_components, freq_map,
                              gxs, Dxs, Dys, gys);

        }//end time loop
      } //end parallactic if else
    } //end freq loop
  } //end tile loop

  if (status != 0) {
    handle_hyperbeam_error(__FILE__, __LINE__, "fee_calc_jones_array");
    // printf("Something went wrong running fee_calc_jones_array\n");
  }

  // ( gpuFree(jones) );

  free(closest_freqs);
  free(jones);

}