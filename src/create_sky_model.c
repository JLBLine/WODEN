/*******************************************************************************
*  Methods to read in and crop a WODEN style sky model to everything
*  above the horizon.
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

/*********************************
// Taken and edited from the RTS (Mitchell et al 2008)
// All credit to the original authors
// https://github.com/ICRAR/mwa-RTS.git
**********************************/
int read_source_catalogue(const char *filename, source_catalogue_t *srccat) {
  int result, n_src=0, n_comps=0;
  int src_found=0, comp_found=0;
  int src_key=0, src_end=0, comp_key=0, comp_end=0, freq_key=0;
  int linear_key=0, type_key=0, param_key=0, coeff_key=0;

  char comp_type[16];
  double ra, dec, freq;
  user_precision_t flux_I, flux_Q, flux_U, flux_V, SI;
  user_precision_t coeff1, coeff2, coeff3;

  int point_ind=0, gauss_ind=0, shape_ind=0;
  int coeff_ind=0;

  //Array that will end up being as long as the number of sources,
  //containing source_t structs
  source_t *srcs=NULL;

  //Total number of SHAPELET COMPONENTs in the entire source catalogue
  //Start it off at zero
  srccat->num_shapelets = 0;

  //Reading in things
  FILE *fp=NULL;
  char line[BUFSIZ];

  /* open file */
  if ((fp=fopen(filename,"r"))==NULL) {
    printf("read_source_catalogue error: failed to open source file <%s>\n", filename);
    // exit(1);
    return 1;
  }

  int line_number = 1;
  //Read in line by line and do "smart" things accordingly
  while(fgets(line,BUFSIZ,fp) != NULL) {
    /* skip blank lines and comments */
    if (line[0]=='\n' || line[0] == '#'|| line[0] == '\0') continue;

    /* lines must start with one of the following keywords... which one is it? */
    linear_key = src_key = src_end = comp_key = comp_end = freq_key = param_key = coeff_key = 0;
    if(strncmp(line, SRC_KEY, strlen(SRC_KEY))==0) src_key=1;
    if(strncmp(line, SRC_END, strlen(SRC_END))==0) src_end=1;
    if(strncmp(line, COMP_KEY, strlen(COMP_KEY))==0) comp_key=1;
    if(strncmp(line, COMP_END, strlen(COMP_END))==0) comp_end=1;
    if(strncmp(line, FREQ_KEY, strlen(FREQ_KEY))==0) freq_key=1;
    if(strncmp(line, LINEAR_KEY, strlen(LINEAR_KEY))==0) linear_key=1;
    if(strncmp(line, GPARAMS_KEY, strlen(GPARAMS_KEY))==0) param_key=1;
    if(strncmp(line, SPARAMS_KEY, strlen(SPARAMS_KEY))==0) param_key=1;
    if(strncmp(line, SCOEFF_KEY, strlen(SCOEFF_KEY))==0) coeff_key=1;

    /* if a source key hasn't been found, skip lines until one is found */
    if ( (src_key==0 && src_end==0 && comp_key==0 && comp_end==0 && freq_key==0 && param_key==0 && coeff_key==0 && linear_key==0) ||
         (src_found==0 && src_key==0) /* i.e., haven't started a new source */ ) {
        // printf("%c%c%c%c %s %d %d\n",line[0],line[1],line[2],line[3],FREQ_KEY,freq_key,strncmp(line, FREQ_KEY, strlen(FREQ_KEY)));
        printf("read_source_catalogue: line %d is bad, please fix this line:\n\t%s\n",line_number, line);
        return 1;
    }

    if (src_end) {
        if (srcs[n_src-1].n_comps == 0){
            printf("WARNING read_source_catalogue: source <%s> had no components\n", srcs[n_src-1].name );
            n_src--;
        }
        src_found=0;
        continue;
    }

    /* is this the start of a new source, or more information for the current source? */
    if (src_found && src_key) {
        printf("read_source_catalogue error: found source key in line: \n %s \n while already reading a source\n",
              line );
        // return NULL;
        return 1;
    }

    else if (src_key) {

      /* found a new source. make more room, get the details */
      src_found=1;
      n_src++;
      srcs = realloc(srcs,sizeof(source_t)*n_src);
      if (srcs == NULL) {
        printf("read_source_catalogue error: no realloc for cat\n");
        // return NULL;
        return 1;
      }

      comp_found=0;
      point_ind=0;
      gauss_ind=0;
      shape_ind=0;
      coeff_ind=0;

      // srcs[n_src-1].name = NULL;
      srcs[n_src-1].n_comps = 0;
      srcs[n_src-1].n_points = 0;
      srcs[n_src-1].n_gauss = 0;
      srcs[n_src-1].n_shapes = 0;
      srcs[n_src-1].n_shape_coeffs = 0;

      //Line should look like <SOURCE name P %d G %d S %d %d>
      //%*s ignores input so we don't have to read in useless input
      result = sscanf( line, "%*s %s %*s %d %*s %d %*s %d %d", srcs[n_src-1].name,
                        &(srcs[n_src-1].n_points), &(srcs[n_src-1].n_gauss),
                        &(srcs[n_src-1].n_shapes), &(srcs[n_src-1].n_shape_coeffs) );
      if (result != 5) {
          printf("read_source_catalogue error: problem reading line %d:\n\t %s \n", line_number,line );
          return 1;
          // return NULL;
      }

      srcs[n_src-1].n_comps = srcs[n_src-1].n_points + srcs[n_src-1].n_gauss + srcs[n_src-1].n_shapes;

      //Update the total count of SHAPELET COMPONENTs in the entire source catalogue
      srccat->num_shapelets += srcs[n_src-1].n_shapes;

      // printf("New source: name: <%s>, Comps:%d  P:%d  G:%d  S:%d-%d \n",
      //       srcs[n_src-1].name, srcs[n_src-1].n_comps, srcs[n_src-1].n_points, srcs[n_src-1].n_gauss,
      //       srcs[n_src-1].n_shapes, srcs[n_src-1].n_shape_coeffs );

      //Now we know the number of sources, do some mallocing
      //Pointsource params
      srcs[n_src-1].point_components.ras = malloc( srcs[n_src-1].n_points * sizeof(double) );
      srcs[n_src-1].point_components.decs = malloc( srcs[n_src-1].n_points * sizeof(double) );
      srcs[n_src-1].point_components.SIs = malloc( srcs[n_src-1].n_points * sizeof(user_precision_t) );
      srcs[n_src-1].point_components.ref_stokesI = malloc( srcs[n_src-1].n_points * sizeof(user_precision_t) );
      srcs[n_src-1].point_components.ref_stokesQ = malloc( srcs[n_src-1].n_points * sizeof(user_precision_t) );
      srcs[n_src-1].point_components.ref_stokesU = malloc( srcs[n_src-1].n_points * sizeof(user_precision_t) );
      srcs[n_src-1].point_components.ref_stokesV = malloc( srcs[n_src-1].n_points * sizeof(user_precision_t) );
      srcs[n_src-1].point_components.ref_freqs = malloc( srcs[n_src-1].n_points * sizeof(double) );

      //Gaussian params
      srcs[n_src-1].gauss_components.ras = malloc( srcs[n_src-1].n_gauss * sizeof(double) );
      srcs[n_src-1].gauss_components.decs = malloc( srcs[n_src-1].n_gauss * sizeof(double) );
      srcs[n_src-1].gauss_components.SIs = malloc( srcs[n_src-1].n_gauss * sizeof(user_precision_t) );
      srcs[n_src-1].gauss_components.ref_stokesI = malloc( srcs[n_src-1].n_gauss * sizeof(user_precision_t) );
      srcs[n_src-1].gauss_components.ref_stokesQ = malloc( srcs[n_src-1].n_gauss * sizeof(user_precision_t) );
      srcs[n_src-1].gauss_components.ref_stokesU = malloc( srcs[n_src-1].n_gauss * sizeof(user_precision_t) );
      srcs[n_src-1].gauss_components.ref_stokesV = malloc( srcs[n_src-1].n_gauss * sizeof(user_precision_t) );
      srcs[n_src-1].gauss_components.ref_freqs = malloc( srcs[n_src-1].n_gauss * sizeof(double) );
      srcs[n_src-1].gauss_components.majors = malloc( srcs[n_src-1].n_gauss * sizeof(user_precision_t) );
      srcs[n_src-1].gauss_components.minors = malloc( srcs[n_src-1].n_gauss * sizeof(user_precision_t) );
      srcs[n_src-1].gauss_components.pas = malloc( srcs[n_src-1].n_gauss * sizeof(user_precision_t) );
      //Shapelet params
      srcs[n_src-1].shape_components.ras = malloc( srcs[n_src-1].n_shapes * sizeof(double) );
      srcs[n_src-1].shape_components.decs = malloc( srcs[n_src-1].n_shapes * sizeof(double) );
      srcs[n_src-1].shape_components.ref_stokesI = malloc( srcs[n_src-1].n_shapes * sizeof(user_precision_t) );
      srcs[n_src-1].shape_components.SIs = malloc( srcs[n_src-1].n_shapes * sizeof(user_precision_t) );
      srcs[n_src-1].shape_components.ref_stokesQ = malloc( srcs[n_src-1].n_shapes * sizeof(user_precision_t) );
      srcs[n_src-1].shape_components.ref_stokesU = malloc( srcs[n_src-1].n_shapes * sizeof(user_precision_t) );
      srcs[n_src-1].shape_components.ref_stokesV = malloc( srcs[n_src-1].n_shapes * sizeof(user_precision_t) );
      srcs[n_src-1].shape_components.ref_freqs = malloc( srcs[n_src-1].n_shapes * sizeof(double) );
      srcs[n_src-1].shape_components.majors = malloc( srcs[n_src-1].n_shapes * sizeof(user_precision_t) );
      srcs[n_src-1].shape_components.minors = malloc( srcs[n_src-1].n_shapes * sizeof(user_precision_t) );
      srcs[n_src-1].shape_components.pas = malloc( srcs[n_src-1].n_shapes * sizeof(user_precision_t) );
      //Need bigger arrays to hold the coeffs - make a shape_param_indexes
      //array so we can map each n1.n2,coeff to the corresponding pa,major,minor
      //when we have multiple shapelet models within a source
      srcs[n_src-1].shape_components.n1s = malloc( srcs[n_src-1].n_shape_coeffs * sizeof(user_precision_t) );
      srcs[n_src-1].shape_components.n2s = malloc( srcs[n_src-1].n_shape_coeffs * sizeof(user_precision_t) );
      srcs[n_src-1].shape_components.shape_coeffs = malloc( srcs[n_src-1].n_shape_coeffs * sizeof(user_precision_t) );
      srcs[n_src-1].shape_components.param_indexes = malloc( srcs[n_src-1].n_shape_coeffs * sizeof(user_precision_t) );

    }

    else if (comp_key) {

      /* found a new source component. get the details */
      comp_found=1;
      n_comps++;

      result = sscanf( line, "%*s %s %lf %lf", comp_type, &ra, &dec );
      if (result != 3) {
          printf("read_source_catalogue error %d: problem reading cal component input line number %d\n\t %s \n", result, line_number, line );
          // return NULL;
          return 1;
      }

      if(strncmp(comp_type, POINT_KEY, strlen(POINT_KEY))==0) type_key=POINT;
      if(strncmp(comp_type, GAUSSIAN_KEY, strlen(GAUSSIAN_KEY))==0) type_key=GAUSSIAN;
      if(strncmp(comp_type, SHAPELET_KEY, strlen(SHAPELET_KEY))==0) type_key=SHAPELET;

      switch (type_key) {
        case POINT:
          point_ind++;
          /* convert to radian */
          srcs[n_src-1].point_components.ras[point_ind-1] = ra * DH2R;
          srcs[n_src-1].point_components.decs[point_ind-1] = dec * DD2R;
          // LOGV( LOG_LOW, "New source component: <%s> component %d, ra: %8.5f hrs, dec: %+9.5f deg\n",
          //       cat[n_src-1].components[n_comps-1].ra*DR2H, cat[n_src-1].components[n_comps-1].dec*DR2D );
        break;

        case GAUSSIAN:
          gauss_ind++;
          srcs[n_src-1].gauss_components.ras[gauss_ind-1] = ra * DH2R;
          srcs[n_src-1].gauss_components.decs[gauss_ind-1] = dec * DD2R;
        break;

        case SHAPELET:
          shape_ind++;
          srcs[n_src-1].shape_components.ras[shape_ind-1] = ra * DH2R;
          srcs[n_src-1].shape_components.decs[shape_ind-1] = dec * DD2R;
        break;

        default:
          printf("read_source_catalogue error: unknown source type: %s on line %d:\n\t %s \n", comp_type,line_number, line );
          // return NULL;
          return 1;


      }//switch (type_key)
    }//if (comp_key)

    else if (freq_key && comp_found==1) {

      //If doing double, read in as double, otherwise float
      #ifdef DOUBLE_PRECISION
        result = sscanf( line, "%*s %lf %lf %lf %lf %lf", &freq, &flux_I, &flux_Q, &flux_U, &flux_V );
      #else
        result = sscanf( line, "%*s %lf %f %f %f %f", &freq, &flux_I, &flux_Q, &flux_U, &flux_V );
      #endif
      if (result != 5) {
          printf("read_source_catalogue error %d: problem reading cal component input line number %d\n\t %s \n", result, line_number, line );
          // return NULL;
          return 1;
      }

      switch (type_key) {
        case POINT:

          srcs[n_src-1].point_components.ref_freqs[point_ind-1] = freq;
          srcs[n_src-1].point_components.ref_stokesI[point_ind-1] = flux_I;
          srcs[n_src-1].point_components.ref_stokesQ[point_ind-1] = flux_Q;
          srcs[n_src-1].point_components.ref_stokesU[point_ind-1] = flux_U;
          srcs[n_src-1].point_components.ref_stokesV[point_ind-1] = flux_V;
          srcs[n_src-1].point_components.SIs[point_ind-1] = DEFAULT_SI;
          // LOGV( LOG_LOW, "New source component: <%s> component %d, ra: %8.5f hrs, dec: %+9.5f deg\n",
          //       cat[n_src-1].components[n_comps-1].ra*DR2H, cat[n_src-1].components[n_comps-1].dec*DR2D );
        break;

        case GAUSSIAN:
          srcs[n_src-1].gauss_components.ref_freqs[gauss_ind-1] = freq;
          srcs[n_src-1].gauss_components.ref_stokesI[gauss_ind-1] = flux_I;
          srcs[n_src-1].gauss_components.ref_stokesQ[gauss_ind-1] = flux_Q;
          srcs[n_src-1].gauss_components.ref_stokesU[gauss_ind-1] = flux_U;
          srcs[n_src-1].gauss_components.ref_stokesV[gauss_ind-1] = flux_V;
          srcs[n_src-1].gauss_components.SIs[gauss_ind-1] = DEFAULT_SI;
        break;

        case SHAPELET:
          srcs[n_src-1].shape_components.ref_freqs[shape_ind-1] = freq;
          srcs[n_src-1].shape_components.ref_stokesI[shape_ind-1] = flux_I;
          srcs[n_src-1].shape_components.ref_stokesQ[shape_ind-1] = flux_Q;
          srcs[n_src-1].shape_components.ref_stokesU[shape_ind-1] = flux_U;
          srcs[n_src-1].shape_components.ref_stokesV[shape_ind-1] = flux_V;
          srcs[n_src-1].shape_components.SIs[shape_ind-1] = DEFAULT_SI;
        break;

      default:
        printf("%d read_source_catalogue error assigning info from line %d:\n\t<%s>\n", type_key,  line_number, line );
        // return NULL;
        return 1;

      }//switch (type_key)

    }//if (freq_key && comp_found==1)


    else if (linear_key && comp_found==1) {

      #ifdef DOUBLE_PRECISION
        result = sscanf( line, "%*s %lf %lf %lf %lf %lf %lf", &freq, &flux_I, &flux_Q, &flux_U, &flux_V, &SI );
      #else
        result = sscanf( line, "%*s %lf %f %f %f %f %f", &freq, &flux_I, &flux_Q, &flux_U, &flux_V, &SI );
      #endif

      if (result != 6) {
          printf("%d read_source_catalogue error: problem reading cal component input from line %d:\n\t<%s>\n", result, line_number, line );
          // return NULL;
          return 1;
      }

      switch (type_key) {
        case POINT:

          srcs[n_src-1].point_components.ref_freqs[point_ind-1] = freq;
          srcs[n_src-1].point_components.ref_stokesI[point_ind-1] = flux_I;
          srcs[n_src-1].point_components.ref_stokesQ[point_ind-1] = flux_Q;
          srcs[n_src-1].point_components.ref_stokesU[point_ind-1] = flux_U;
          srcs[n_src-1].point_components.ref_stokesV[point_ind-1] = flux_V;
          srcs[n_src-1].point_components.SIs[point_ind-1] = SI;
        break;

        case GAUSSIAN:
          srcs[n_src-1].gauss_components.ref_freqs[gauss_ind-1] = freq;
          srcs[n_src-1].gauss_components.ref_stokesI[gauss_ind-1] = flux_I;
          srcs[n_src-1].gauss_components.ref_stokesQ[gauss_ind-1] = flux_Q;
          srcs[n_src-1].gauss_components.ref_stokesU[gauss_ind-1] = flux_U;
          srcs[n_src-1].gauss_components.ref_stokesV[gauss_ind-1] = flux_V;
          srcs[n_src-1].gauss_components.SIs[gauss_ind-1] = SI;
        break;

        case SHAPELET:
          srcs[n_src-1].shape_components.ref_freqs[shape_ind-1] = freq;
          srcs[n_src-1].shape_components.ref_stokesI[shape_ind-1] = flux_I;
          srcs[n_src-1].shape_components.ref_stokesQ[shape_ind-1] = flux_Q;
          srcs[n_src-1].shape_components.ref_stokesU[shape_ind-1] = flux_U;
          srcs[n_src-1].shape_components.ref_stokesV[shape_ind-1] = flux_V;
          srcs[n_src-1].shape_components.SIs[shape_ind-1] = SI;
        break;

      default:
        printf("%d read_source_catalogue error assigning info from line %d:\n\t<%s>\n", type_key,  line_number, line );
        // return NULL;
        return 1;

      }//switch (type_key)

    }//if (freq_key && comp_found==1)

    else if (param_key && comp_found==1) {
      #ifdef DOUBLE_PRECISION
        result = sscanf( line, "%*s %lf %lf %lf", &coeff1, &coeff2, &coeff3 );
      #else
        result = sscanf( line, "%*s %f %f %f", &coeff1, &coeff2, &coeff3 );
      #endif

      if (result != 3) {
          printf("%d read_source_catalogue error: problem reading component input from line %d:\n\t<%s>\n", result, line_number, line );
          // return NULL;
          return 1;
      }

      switch (type_key) {
        case GAUSSIAN:
          srcs[n_src-1].gauss_components.pas[gauss_ind-1] =  coeff1 * DD2R;
          srcs[n_src-1].gauss_components.majors[gauss_ind-1] = coeff2 * (DD2R / 60.0);
          srcs[n_src-1].gauss_components.minors[gauss_ind-1] = coeff3 * (DD2R / 60.0);
        break;

        case SHAPELET:
          srcs[n_src-1].shape_components.pas[shape_ind-1] =  coeff1 * DD2R;
          srcs[n_src-1].shape_components.majors[shape_ind-1] = coeff2 * (DD2R / 60.0);
          srcs[n_src-1].shape_components.minors[shape_ind-1] = coeff3 * (DD2R / 60.0);
        break;

      default:
        printf("read_source_catalogue error: problem reading component input from line %d:\n\t %s \n", line_number, line );
        // return NULL;
        return 1;
      }//switch (type_key)
    }//if (param_key && comp_found==1)

    else if (coeff_key && comp_found==1) {
      #ifdef DOUBLE_PRECISION
        result = sscanf( line, "%*s %lf %lf %lf", &coeff1, &coeff2, &coeff3 );
      #else
        result = sscanf( line, "%*s %f %f %f", &coeff1, &coeff2, &coeff3 );
      #endif
      if (result != 3) {
          printf("read_source_catalogue error: problem reading coeffs component input from line %d:\n\t %s \n", line_number, line );
          // return NULL;
          return 1;
      }

      switch (type_key) {
        case SHAPELET:
          srcs[n_src-1].shape_components.n1s[coeff_ind] = coeff1;
          srcs[n_src-1].shape_components.n2s[coeff_ind] = coeff2;
          srcs[n_src-1].shape_components.shape_coeffs[coeff_ind] = coeff3;
          srcs[n_src-1].shape_components.param_indexes[coeff_ind] = shape_ind - 1;
        break;

      default:
        printf("read_source_catalogue error: problem reading input from line %d:\n\t<%s>\n", line_number, line );
        // return NULL;
        return 1;
      }//switch (type_key)
      coeff_ind++;
    }//if (param_key && comp_found==1)

    else if (comp_end==1) {
      comp_found = 0;
    }
    line_number += 1;
  }//while fgets(line,BUFSIZ,fp) != NULL)

  /* set pointers and counts in the final catalogue object */
  srccat->sources = srcs;
  srccat->num_sources = n_src;

  return 0;
} //read_sources


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
// Takes a zenith angle and checks if it's below the horizon. Depending on
// whether we are cropping all SOURCEs partially below the horizon, or
// retaining all COMPONENTs, update either all_comps_above_horizon or
// sky_crop_type. If the component is a shapelet, we also need to cycle
// through the catsource.shape_param_indexes, checking which ones match this
// shape index (int shape), and update num_shape_coeff_retained accordingly
*******************************************************************************/
void horizon_test(double za, e_sky_crop sky_crop_type,
     e_horizon * all_comps_above_horizon, int * num_comp_retained,
     int * num_shape_coeff_retained, int num_shape_coeff_component,
     user_precision_t *shape_param_indexes, int shape){
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
        for (int param_index = 0; param_index < num_shape_coeff_component; param_index++) {
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

  //Begin checking az/za loop here - for each SOURCE
  for (int src = 0; src < raw_srccat->num_sources; src++){
    all_comps_above_horizon = ABOVE;

    //Begin point source component loop
    for (int point = 0; point < raw_srccat->sources[src].n_points; point++) {
      //Calculate az/za for all point components
      convert_radec2azza(raw_srccat->sources[src].point_components.ras[point],
                         raw_srccat->sources[src].point_components.decs[point],
                         lsts[0], latitude, &az, &za);

      raw_srccat->sources[src].point_components.azs[point] = (user_precision_t)az;
      raw_srccat->sources[src].point_components.zas[point] = (user_precision_t)za;
      //Check if components are above the horizon, and count how many
      //components survive / flag a source if a component is below the horizon
      //Last three arguments only used for shapelets so pass 0
      horizon_test(za, sky_crop_type, &all_comps_above_horizon,
                  &num_point_comp_retained, &num_shape_coeff_retained,
                  0, raw_srccat->sources[src].shape_components.param_indexes, 0);

    }//End point source component loop

    //Begin gauss source component loop
    for (int gauss = 0; gauss < raw_srccat->sources[src].n_gauss; gauss++) {
        //Calculate az/za for all gauss components
      convert_radec2azza(raw_srccat->sources[src].gauss_components.ras[gauss],
                         raw_srccat->sources[src].gauss_components.decs[gauss],
                         lsts[0], latitude, &az, &za);
      raw_srccat->sources[src].gauss_components.azs[gauss] = (user_precision_t)az;
      raw_srccat->sources[src].gauss_components.zas[gauss] = (user_precision_t)za;
      //Check if components are above the horizon, and count how many
      //components survive / flag a source if a component is below the horizon
      //Last three arguments only used for shapelets so pass 0
      horizon_test(za, sky_crop_type, &all_comps_above_horizon,
                  &num_gauss_comp_retained, &num_shape_coeff_retained,
                  0, raw_srccat->sources[src].shape_components.param_indexes, 0);

    }//End gauss source component loop

    //Begin shape source component loop
    for (int shape = 0; shape < raw_srccat->sources[src].n_shapes; shape++) {
      //Calculate az/za for all gauss components
      convert_radec2azza(raw_srccat->sources[src].shape_components.ras[shape],
                         raw_srccat->sources[src].shape_components.decs[shape],
                         lsts[0], latitude, &az, &za);
      raw_srccat->sources[src].shape_components.azs[shape] = (user_precision_t)az;
      raw_srccat->sources[src].shape_components.zas[shape] = (user_precision_t)za;
      //Check if components are above the horizon, and count how many
      //components survive / flag a source if a component is below the horizon
      //Use last three arguments to correctly identify which shapelet coeffs
      //belong to which shapelet component so we can malloc them correctly later
      horizon_test(za, sky_crop_type, &all_comps_above_horizon,
                  &num_shape_comp_retained, &num_shape_coeff_retained,
                  raw_srccat->sources[src].n_shape_coeffs,
                  raw_srccat->sources[src].shape_components.param_indexes, shape);

    }//End shape source component loop

    //After checking all components in a source, if cropping out sources,
    //check all components were above horizon. If so, update component type
    //counters, and mark down the index of the sources that survied in
    //cropped_src_indexes
    if (sky_crop_type == CROP_SOURCES) {
      if (all_comps_above_horizon == ABOVE) {
        num_sources_retained ++;
        num_point_comp_retained += raw_srccat->sources[src].n_points;
        num_gauss_comp_retained += raw_srccat->sources[src].n_gauss;
        num_shape_comp_retained += raw_srccat->sources[src].n_shapes;
        num_shape_coeff_retained += raw_srccat->sources[src].n_shape_coeffs;
        cropped_src_indexes = realloc(cropped_src_indexes,sizeof(int)*num_sources_retained);
        cropped_src_indexes[num_sources_retained - 1] = src;
      }//end if all_comps_above_horizon == ABOVE
    }//end if sky_crop_type == CROP_SOURCES

  }//Finish checking az/za loop here

  //Make an empty source_t and malloc using the numbers gather above
  source_t *cropped_src=NULL;
  cropped_src = malloc(sizeof(source_t));

  cropped_src->n_points = num_point_comp_retained;
  cropped_src->n_gauss = num_gauss_comp_retained;
  cropped_src->n_shapes = num_shape_comp_retained;
  cropped_src->n_comps = num_point_comp_retained + num_gauss_comp_retained + num_shape_comp_retained;
  cropped_src->n_shape_coeffs = num_shape_coeff_retained;

  cropped_src->point_components.ras = malloc( num_point_comp_retained * sizeof(double) );
  cropped_src->point_components.decs = malloc( num_point_comp_retained * sizeof(double) );
  cropped_src->point_components.ref_stokesI = malloc( num_point_comp_retained * sizeof(user_precision_t) );
  cropped_src->point_components.ref_stokesQ = malloc( num_point_comp_retained * sizeof(user_precision_t) );
  cropped_src->point_components.ref_stokesU = malloc( num_point_comp_retained * sizeof(user_precision_t) );
  cropped_src->point_components.ref_stokesV = malloc( num_point_comp_retained * sizeof(user_precision_t) );
  cropped_src->point_components.ref_freqs = malloc( num_point_comp_retained * sizeof(double) );
  cropped_src->point_components.SIs = malloc( num_point_comp_retained * sizeof(user_precision_t) );

  cropped_src->point_components.azs = malloc( num_point_comp_retained * num_time_steps * sizeof(user_precision_t) );
  cropped_src->point_components.zas = malloc( num_point_comp_retained * num_time_steps * sizeof(user_precision_t) );

  cropped_src->gauss_components.ras = malloc( num_gauss_comp_retained * sizeof(double) );
  cropped_src->gauss_components.decs = malloc( num_gauss_comp_retained * sizeof(double) );
  cropped_src->gauss_components.ref_stokesI = malloc( num_gauss_comp_retained * sizeof(user_precision_t) );
  cropped_src->gauss_components.ref_stokesQ = malloc( num_gauss_comp_retained * sizeof(user_precision_t) );
  cropped_src->gauss_components.ref_stokesU = malloc( num_gauss_comp_retained * sizeof(user_precision_t) );
  cropped_src->gauss_components.ref_stokesV = malloc( num_gauss_comp_retained * sizeof(user_precision_t) );
  cropped_src->gauss_components.ref_freqs = malloc( num_gauss_comp_retained * sizeof(double) );
  cropped_src->gauss_components.SIs = malloc( num_gauss_comp_retained * sizeof(user_precision_t) );
  cropped_src->gauss_components.majors = malloc( num_gauss_comp_retained * sizeof(user_precision_t) );
  cropped_src->gauss_components.minors = malloc( num_gauss_comp_retained * sizeof(user_precision_t) );
  cropped_src->gauss_components.pas = malloc( num_gauss_comp_retained * sizeof(user_precision_t) );
  cropped_src->gauss_components.azs = malloc( num_gauss_comp_retained * num_time_steps * sizeof(user_precision_t) );
  cropped_src->gauss_components.zas = malloc( num_gauss_comp_retained * num_time_steps * sizeof(user_precision_t) );

  cropped_src->shape_components.ras = malloc( num_shape_comp_retained * sizeof(double) );
  cropped_src->shape_components.decs = malloc( num_shape_comp_retained * sizeof(double) );
  cropped_src->shape_components.ref_stokesI = malloc( num_shape_comp_retained * sizeof(user_precision_t) );
  cropped_src->shape_components.ref_stokesQ = malloc( num_shape_comp_retained * sizeof(user_precision_t) );
  cropped_src->shape_components.ref_stokesU = malloc( num_shape_comp_retained * sizeof(user_precision_t) );
  cropped_src->shape_components.ref_stokesV = malloc( num_shape_comp_retained * sizeof(user_precision_t) );
  cropped_src->shape_components.SIs = malloc( num_shape_comp_retained * sizeof(user_precision_t) );
  cropped_src->shape_components.ref_freqs = malloc( num_shape_comp_retained * sizeof(double) );
  cropped_src->shape_components.majors = malloc( num_shape_comp_retained * sizeof(user_precision_t) );
  cropped_src->shape_components.minors = malloc( num_shape_comp_retained * sizeof(user_precision_t) );
  cropped_src->shape_components.pas = malloc( num_shape_comp_retained * sizeof(user_precision_t) );

  cropped_src->shape_components.n1s = malloc( num_shape_coeff_retained * sizeof(user_precision_t) );
  cropped_src->shape_components.n2s = malloc( num_shape_coeff_retained * sizeof(user_precision_t) );
  cropped_src->shape_components.shape_coeffs = malloc( num_shape_coeff_retained * sizeof(user_precision_t) );
  cropped_src->shape_components.param_indexes = malloc( num_shape_coeff_retained * sizeof(user_precision_t) );
  cropped_src->shape_components.azs = malloc( num_shape_comp_retained * num_time_steps * sizeof(user_precision_t) );
  cropped_src->shape_components.zas = malloc( num_shape_comp_retained * num_time_steps * sizeof(user_precision_t) );

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
    for (int retained = 0; retained < num_sources_retained; retained++) {
      int src = cropped_src_indexes[retained];

      //Loop over point components
      for (int point = 0; point < raw_srccat->sources[src].n_points; point++){
        cropped_src->point_components.ras[point_crop_component_index] = raw_srccat->sources[src].point_components.ras[point];
        cropped_src->point_components.decs[point_crop_component_index] = raw_srccat->sources[src].point_components.decs[point];
        cropped_src->point_components.ref_stokesI[point_crop_component_index] = raw_srccat->sources[src].point_components.ref_stokesI[point];
        cropped_src->point_components.ref_stokesQ[point_crop_component_index] = raw_srccat->sources[src].point_components.ref_stokesQ[point];
        cropped_src->point_components.ref_stokesU[point_crop_component_index] = raw_srccat->sources[src].point_components.ref_stokesU[point];
        cropped_src->point_components.ref_stokesV[point_crop_component_index] = raw_srccat->sources[src].point_components.ref_stokesV[point];
        cropped_src->point_components.ref_freqs[point_crop_component_index] = raw_srccat->sources[src].point_components.ref_freqs[point];
        cropped_src->point_components.SIs[point_crop_component_index] = raw_srccat->sources[src].point_components.SIs[point];

        //Calculate az/za values for each point for all time steps
        for (int time_step = 0; time_step < num_time_steps; time_step++) {
          convert_radec2azza(raw_srccat->sources[src].point_components.ras[point],
                             raw_srccat->sources[src].point_components.decs[point],
                             lsts[time_step], latitude, &az, &za);
          cropped_src->point_components.azs[point_crop_component_index*num_time_steps + time_step] = (user_precision_t)az;
          cropped_src->point_components.zas[point_crop_component_index*num_time_steps + time_step] = (user_precision_t)za;
        }

        point_crop_component_index ++;

      }//End point component loop

      //Loop over gauss components
      for (int gauss = 0; gauss < raw_srccat->sources[src].n_gauss; gauss++){
        cropped_src->gauss_components.ras[gauss_crop_component_index] = raw_srccat->sources[src].gauss_components.ras[gauss];
        cropped_src->gauss_components.decs[gauss_crop_component_index] = raw_srccat->sources[src].gauss_components.decs[gauss];
        cropped_src->gauss_components.ref_stokesI[gauss_crop_component_index] = raw_srccat->sources[src].gauss_components.ref_stokesI[gauss];
        cropped_src->gauss_components.ref_stokesQ[gauss_crop_component_index] = raw_srccat->sources[src].gauss_components.ref_stokesQ[gauss];
        cropped_src->gauss_components.ref_stokesU[gauss_crop_component_index] = raw_srccat->sources[src].gauss_components.ref_stokesU[gauss];
        cropped_src->gauss_components.ref_stokesV[gauss_crop_component_index] = raw_srccat->sources[src].gauss_components.ref_stokesV[gauss];
        cropped_src->gauss_components.ref_freqs[gauss_crop_component_index] = raw_srccat->sources[src].gauss_components.ref_freqs[gauss];
        cropped_src->gauss_components.SIs[gauss_crop_component_index] = raw_srccat->sources[src].gauss_components.SIs[gauss];

        cropped_src->gauss_components.majors[gauss_crop_component_index] = raw_srccat->sources[src].gauss_components.majors[gauss];
        cropped_src->gauss_components.minors[gauss_crop_component_index] = raw_srccat->sources[src].gauss_components.minors[gauss];
        cropped_src->gauss_components.pas[gauss_crop_component_index] = raw_srccat->sources[src].gauss_components.pas[gauss];

        //Calculate az/za values for each gauss for all time steps
        for (int time_step = 0; time_step < num_time_steps; time_step++) {
          convert_radec2azza(raw_srccat->sources[src].gauss_components.ras[gauss],
                             raw_srccat->sources[src].gauss_components.decs[gauss],
                             lsts[time_step], latitude, &az, &za);
          cropped_src->gauss_components.azs[gauss_crop_component_index*num_time_steps + time_step] = (user_precision_t)az;
          cropped_src->gauss_components.zas[gauss_crop_component_index*num_time_steps + time_step] = (user_precision_t)za;
        }

        gauss_crop_component_index ++;

      }//End gauss component loop
      //
      //Loop over shapelet components
      for (int shape = 0; shape < raw_srccat->sources[src].n_shapes; shape++){
        cropped_src->shape_components.ras[shape_crop_component_index] = raw_srccat->sources[src].shape_components.ras[shape];
        cropped_src->shape_components.decs[shape_crop_component_index] = raw_srccat->sources[src].shape_components.decs[shape];

        cropped_src->shape_components.ref_stokesI[shape_crop_component_index] = raw_srccat->sources[src].shape_components.ref_stokesI[shape];
        cropped_src->shape_components.ref_stokesQ[shape_crop_component_index] = raw_srccat->sources[src].shape_components.ref_stokesQ[shape];
        cropped_src->shape_components.ref_stokesU[shape_crop_component_index] = raw_srccat->sources[src].shape_components.ref_stokesU[shape];
        cropped_src->shape_components.ref_stokesV[shape_crop_component_index] = raw_srccat->sources[src].shape_components.ref_stokesV[shape];
        cropped_src->shape_components.ref_freqs[shape_crop_component_index] = raw_srccat->sources[src].shape_components.ref_freqs[shape];
        cropped_src->shape_components.SIs[shape_crop_component_index] = raw_srccat->sources[src].shape_components.SIs[shape];

        cropped_src->shape_components.majors[shape_crop_component_index] = raw_srccat->sources[src].shape_components.majors[shape];
        cropped_src->shape_components.minors[shape_crop_component_index] = raw_srccat->sources[src].shape_components.minors[shape];
        cropped_src->shape_components.pas[shape_crop_component_index] = raw_srccat->sources[src].shape_components.pas[shape];

        //Calculate az/za values for each shapelet component for all time steps
        for (int time_step = 0; time_step < num_time_steps; time_step++) {
          convert_radec2azza(raw_srccat->sources[src].shape_components.ras[shape],
                             raw_srccat->sources[src].shape_components.decs[shape],
                             lsts[time_step], latitude, &az, &za);
          cropped_src->shape_components.azs[shape_crop_component_index*num_time_steps + time_step] = (user_precision_t)az;
          cropped_src->shape_components.zas[shape_crop_component_index*num_time_steps + time_step] = (user_precision_t)za;
        }

        if ((int)shape == 0) {
          // Loop over the coefficients for this shapelet source
          // Only do it once as all shapelet component coeffs, n1s, n2s for one
          // source are in 1D arrays in raw_srccat->sources[src]. As each
          // shapelet component can have any number of coeffs, n1s, n2s, we
          // relate the coeffs, n1s, n2s, to shapelet ra, dec, etc via the
          // cropped_src->shape_components.param_indexes array. So need to pull that
          // information out from raw_srccat->sources[src].shape_components.param_indexes
          // and keep track of how many shapelet coeff components are in
          // the new cropped_src using shape_coeff_component_index
          // Do all this work now as the 1D array goes nicely into a GPU kernel
          for (int coeff_ind = 0; coeff_ind < raw_srccat->sources[src].n_shape_coeffs; coeff_ind++) {
            cropped_src->shape_components.shape_coeffs[shape_coeff_component_index] = raw_srccat->sources[src].shape_components.shape_coeffs[coeff_ind];
            cropped_src->shape_components.n1s[shape_coeff_component_index] = raw_srccat->sources[src].shape_components.n1s[coeff_ind];
            cropped_src->shape_components.n2s[shape_coeff_component_index] = raw_srccat->sources[src].shape_components.n2s[coeff_ind];

            //We do the loop of shapelet coeffs, n1s, n2s only once per source, so we may well get ahead of the
            //shape_crop_component_index. So use to shape_param_index to account for that
            int shape_param_index = raw_srccat->sources[src].shape_components.param_indexes[coeff_ind] + shape_crop_component_index;

            cropped_src->shape_components.param_indexes[shape_coeff_component_index] = shape_param_index;
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
    for (int src = 0; src < raw_srccat->num_sources; src++){

      //Begin point component loop
      for (int point = 0; point < raw_srccat->sources[src].n_points; point++){
        //Check if point component above horizon
        if (raw_srccat->sources[src].point_components.zas[point] < M_PI / 2.0){
          cropped_src->point_components.ras[point_crop_component_index] = raw_srccat->sources[src].point_components.ras[point];
          cropped_src->point_components.decs[point_crop_component_index] = raw_srccat->sources[src].point_components.decs[point];
          cropped_src->point_components.ref_stokesI[point_crop_component_index] = raw_srccat->sources[src].point_components.ref_stokesI[point];
          cropped_src->point_components.ref_stokesQ[point_crop_component_index] = raw_srccat->sources[src].point_components.ref_stokesQ[point];
          cropped_src->point_components.ref_stokesU[point_crop_component_index] = raw_srccat->sources[src].point_components.ref_stokesU[point];
          cropped_src->point_components.ref_stokesV[point_crop_component_index] = raw_srccat->sources[src].point_components.ref_stokesV[point];
          cropped_src->point_components.ref_freqs[point_crop_component_index] = raw_srccat->sources[src].point_components.ref_freqs[point];
          cropped_src->point_components.SIs[point_crop_component_index] = raw_srccat->sources[src].point_components.SIs[point];

          //Calculate az/za values for each point for all time steps
          for (int time_step = 0; time_step < num_time_steps; time_step++) {
            convert_radec2azza(raw_srccat->sources[src].point_components.ras[point],
                               raw_srccat->sources[src].point_components.decs[point],
                               lsts[time_step], latitude, &az, &za);
            cropped_src->point_components.azs[point_crop_component_index*num_time_steps + time_step] = (user_precision_t)az;
            cropped_src->point_components.zas[point_crop_component_index*num_time_steps + time_step] = (user_precision_t)za;
            // printf("AZ ZA %d %.7f %.7f\n",point_crop_component_index*num_time_steps + time_step,(float)az,(float)za );
          }//End az/za calculation loop

          point_crop_component_index ++;
        }//End point above horizon loop
      }//End point component loop


      //Begin gauss component loop
      for (int gauss = 0; gauss < raw_srccat->sources[src].n_gauss; gauss++){
        //Check if gauss component above horizon
        if (raw_srccat->sources[src].gauss_components.zas[gauss] < M_PI / 2.0){
          cropped_src->gauss_components.ras[gauss_crop_component_index] = raw_srccat->sources[src].gauss_components.ras[gauss];
          cropped_src->gauss_components.decs[gauss_crop_component_index] = raw_srccat->sources[src].gauss_components.decs[gauss];
          cropped_src->gauss_components.ref_stokesI[gauss_crop_component_index] = raw_srccat->sources[src].gauss_components.ref_stokesI[gauss];
          cropped_src->gauss_components.ref_stokesQ[gauss_crop_component_index] = raw_srccat->sources[src].gauss_components.ref_stokesQ[gauss];
          cropped_src->gauss_components.ref_stokesU[gauss_crop_component_index] = raw_srccat->sources[src].gauss_components.ref_stokesU[gauss];
          cropped_src->gauss_components.ref_stokesV[gauss_crop_component_index] = raw_srccat->sources[src].gauss_components.ref_stokesV[gauss];
          cropped_src->gauss_components.ref_freqs[gauss_crop_component_index] = raw_srccat->sources[src].gauss_components.ref_freqs[gauss];
          cropped_src->gauss_components.SIs[gauss_crop_component_index] = raw_srccat->sources[src].gauss_components.SIs[gauss];

          cropped_src->gauss_components.majors[gauss_crop_component_index] = raw_srccat->sources[src].gauss_components.majors[gauss];
          cropped_src->gauss_components.minors[gauss_crop_component_index] = raw_srccat->sources[src].gauss_components.minors[gauss];
          cropped_src->gauss_components.pas[gauss_crop_component_index] = raw_srccat->sources[src].gauss_components.pas[gauss];

          //Calculate az/za values for each gauss for all time steps
          for (int time_step = 0; time_step < num_time_steps; time_step++) {
            convert_radec2azza(raw_srccat->sources[src].gauss_components.ras[gauss],
                               raw_srccat->sources[src].gauss_components.decs[gauss],
                               lsts[time_step], latitude, &az, &za);
            cropped_src->gauss_components.azs[gauss_crop_component_index*num_time_steps + time_step] = (user_precision_t)az;
            cropped_src->gauss_components.zas[gauss_crop_component_index*num_time_steps + time_step] = (user_precision_t)za;
          }//End az/za calculation loop

          gauss_crop_component_index ++;
        }//End gauss above horizon loop
      }//End gauss component loop

      //Loop over shapelet components
      for (int shape = 0; shape < raw_srccat->sources[src].n_shapes; shape++){

        if (raw_srccat->sources[src].shape_components.zas[shape] < M_PI / 2.0){
          cropped_src->shape_components.ras[shape_crop_component_index] = raw_srccat->sources[src].shape_components.ras[shape];
          cropped_src->shape_components.decs[shape_crop_component_index] = raw_srccat->sources[src].shape_components.decs[shape];
          cropped_src->shape_components.ref_stokesI[shape_crop_component_index] = raw_srccat->sources[src].shape_components.ref_stokesI[shape];
          cropped_src->shape_components.ref_stokesQ[shape_crop_component_index] = raw_srccat->sources[src].shape_components.ref_stokesQ[shape];
          cropped_src->shape_components.ref_stokesU[shape_crop_component_index] = raw_srccat->sources[src].shape_components.ref_stokesU[shape];
          cropped_src->shape_components.ref_stokesV[shape_crop_component_index] = raw_srccat->sources[src].shape_components.ref_stokesV[shape];
          cropped_src->shape_components.ref_freqs[shape_crop_component_index] = raw_srccat->sources[src].shape_components.ref_freqs[shape];
          cropped_src->shape_components.SIs[shape_crop_component_index] = raw_srccat->sources[src].shape_components.SIs[shape];

          cropped_src->shape_components.majors[shape_crop_component_index] = raw_srccat->sources[src].shape_components.majors[shape];
          cropped_src->shape_components.minors[shape_crop_component_index] = raw_srccat->sources[src].shape_components.minors[shape];
          cropped_src->shape_components.pas[shape_crop_component_index] = raw_srccat->sources[src].shape_components.pas[shape];

          //Calculate az/za values for each shape for all time steps
          for (int time_step = 0; time_step < num_time_steps; time_step++) {
            convert_radec2azza(raw_srccat->sources[src].shape_components.ras[shape],
                               raw_srccat->sources[src].shape_components.decs[shape],
                               lsts[time_step], latitude, &az, &za);
            cropped_src->shape_components.azs[shape_crop_component_index*num_time_steps + time_step] = (user_precision_t)az;
            cropped_src->shape_components.zas[shape_crop_component_index*num_time_steps + time_step] = (user_precision_t)za;
          }//End az/za calculation loop

          //Loop through all shapelet coeffs,n1s,n2s for this source; these 1D arrays
          //can contain information from multiple shapelet components.
          //raw_srccat->sources[src].shape_components.param_indexes is used to match
          //the coeff,n1,n2 to each component, so check those indexes and grab
          //information if correct.
          // Do all this work now as the 1D array goes nicely into a GPU kernel
          for (int coeff_ind = 0; coeff_ind < raw_srccat->sources[src].n_shape_coeffs; coeff_ind++) {
            //Check if we are on the cofrect component
            if ( (int)raw_srccat->sources[src].shape_components.param_indexes[coeff_ind] == shape ){
              cropped_src->shape_components.shape_coeffs[shape_coeff_component_index] = raw_srccat->sources[src].shape_components.shape_coeffs[coeff_ind];
              cropped_src->shape_components.n1s[shape_coeff_component_index] = raw_srccat->sources[src].shape_components.n1s[coeff_ind];
              cropped_src->shape_components.n2s[shape_coeff_component_index] = raw_srccat->sources[src].shape_components.n2s[coeff_ind];
              cropped_src->shape_components.param_indexes[shape_coeff_component_index] = shape_crop_component_index;

              shape_coeff_component_index += 1;

            }//end test if correct shape coeff
          }// end shapelet coeff loop
          shape_crop_component_index += 1;
        }//End shapelet above horizon loop
      }//End shapelet component loop
    }//End raw_srccat source loop
  }//End if sky_crop_type == CROP_COMPONENTS

  return cropped_src;

}
