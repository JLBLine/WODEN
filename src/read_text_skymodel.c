/*******************************************************************************
*  Methods to read in a WODEN style sky model to everything
*  above the horizon.
*  @author J.L.B. Line
*
*  Please see documentation in ../include/read_text_skymodel.h or online at
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
int read_text_skymodel(const char *filename, source_catalogue_t *srccat) {
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

      //We gonna be realloc a bunch of arrays, so set them to NULL to
      //start with. ALso set all the counter int numbers to 0.0
      source_zero_counters_and_null_components(&srcs[n_src-1]);

      comp_found=0;
      point_ind=0;
      gauss_ind=0;
      shape_ind=0;
      coeff_ind=0;

      //Line should look like <SOURCE name P %d G %d S %d %d>
      //%*s ignores input so we don't have to read in useless input
      result = sscanf( line, "%*s %s %*s %d %*s %d %*s %d %d", srcs[n_src-1].name,
                        &(srcs[n_src-1].n_points), &(srcs[n_src-1].n_gauss),
                        &(srcs[n_src-1].n_shapes), &(srcs[n_src-1].n_shape_coeffs) );

      //For this old sky model type, we only do power law flux behaviours,
      //so set the number of power law models of each COMPONENT equal
      srcs[n_src-1].n_point_powers = srcs[n_src-1].n_points;
      srcs[n_src-1].n_gauss_powers = srcs[n_src-1].n_gauss;
      srcs[n_src-1].n_shape_powers = srcs[n_src-1].n_shapes;

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
      srcs[n_src-1].point_components.power_SIs = malloc( srcs[n_src-1].n_points * sizeof(user_precision_t) );
      srcs[n_src-1].point_components.power_ref_stokesI = malloc( srcs[n_src-1].n_points * sizeof(user_precision_t) );
      srcs[n_src-1].point_components.power_ref_stokesQ = malloc( srcs[n_src-1].n_points * sizeof(user_precision_t) );
      srcs[n_src-1].point_components.power_ref_stokesU = malloc( srcs[n_src-1].n_points * sizeof(user_precision_t) );
      srcs[n_src-1].point_components.power_ref_stokesV = malloc( srcs[n_src-1].n_points * sizeof(user_precision_t) );
      srcs[n_src-1].point_components.power_ref_freqs = malloc( srcs[n_src-1].n_points * sizeof(double) );
      srcs[n_src-1].point_components.power_comp_inds = malloc( srcs[n_src-1].n_points * sizeof(int) );

      //Gaussian params
      srcs[n_src-1].gauss_components.ras = malloc( srcs[n_src-1].n_gauss * sizeof(double) );
      srcs[n_src-1].gauss_components.decs = malloc( srcs[n_src-1].n_gauss * sizeof(double) );
      srcs[n_src-1].gauss_components.power_SIs = malloc( srcs[n_src-1].n_gauss * sizeof(user_precision_t) );
      srcs[n_src-1].gauss_components.power_ref_stokesI = malloc( srcs[n_src-1].n_gauss * sizeof(user_precision_t) );
      srcs[n_src-1].gauss_components.power_ref_stokesQ = malloc( srcs[n_src-1].n_gauss * sizeof(user_precision_t) );
      srcs[n_src-1].gauss_components.power_ref_stokesU = malloc( srcs[n_src-1].n_gauss * sizeof(user_precision_t) );
      srcs[n_src-1].gauss_components.power_ref_stokesV = malloc( srcs[n_src-1].n_gauss * sizeof(user_precision_t) );
      srcs[n_src-1].gauss_components.power_ref_freqs = malloc( srcs[n_src-1].n_gauss * sizeof(double) );
      srcs[n_src-1].gauss_components.power_comp_inds = malloc( srcs[n_src-1].n_points * sizeof(int) );
      srcs[n_src-1].gauss_components.majors = malloc( srcs[n_src-1].n_gauss * sizeof(user_precision_t) );
      srcs[n_src-1].gauss_components.minors = malloc( srcs[n_src-1].n_gauss * sizeof(user_precision_t) );
      srcs[n_src-1].gauss_components.pas = malloc( srcs[n_src-1].n_gauss * sizeof(user_precision_t) );
      //Shapelet params
      srcs[n_src-1].shape_components.ras = malloc( srcs[n_src-1].n_shapes * sizeof(double) );
      srcs[n_src-1].shape_components.decs = malloc( srcs[n_src-1].n_shapes * sizeof(double) );
      srcs[n_src-1].shape_components.power_ref_stokesI = malloc( srcs[n_src-1].n_shapes * sizeof(user_precision_t) );
      srcs[n_src-1].shape_components.power_SIs = malloc( srcs[n_src-1].n_shapes * sizeof(user_precision_t) );
      srcs[n_src-1].shape_components.power_ref_stokesQ = malloc( srcs[n_src-1].n_shapes * sizeof(user_precision_t) );
      srcs[n_src-1].shape_components.power_ref_stokesU = malloc( srcs[n_src-1].n_shapes * sizeof(user_precision_t) );
      srcs[n_src-1].shape_components.power_ref_stokesV = malloc( srcs[n_src-1].n_shapes * sizeof(user_precision_t) );
      srcs[n_src-1].shape_components.power_ref_freqs = malloc( srcs[n_src-1].n_shapes * sizeof(double) );
      srcs[n_src-1].shape_components.power_comp_inds = malloc( srcs[n_src-1].n_points * sizeof(int) );
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
          //For text model only use power law flux behaviour, so just set to array position
          srcs[n_src-1].point_components.power_comp_inds[point_ind-1] = point_ind-1;
        break;

        case GAUSSIAN:
          gauss_ind++;
          srcs[n_src-1].gauss_components.ras[gauss_ind-1] = ra * DH2R;
          srcs[n_src-1].gauss_components.decs[gauss_ind-1] = dec * DD2R;
          srcs[n_src-1].gauss_components.power_comp_inds[gauss_ind-1] = gauss_ind-1;
        break;

        case SHAPELET:
          shape_ind++;
          srcs[n_src-1].shape_components.ras[shape_ind-1] = ra * DH2R;
          srcs[n_src-1].shape_components.decs[shape_ind-1] = dec * DD2R;
          srcs[n_src-1].shape_components.power_comp_inds[shape_ind-1] = shape_ind-1;
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

          srcs[n_src-1].point_components.power_ref_freqs[point_ind-1] = freq;
          srcs[n_src-1].point_components.power_ref_stokesI[point_ind-1] = flux_I;
          srcs[n_src-1].point_components.power_ref_stokesQ[point_ind-1] = flux_Q;
          srcs[n_src-1].point_components.power_ref_stokesU[point_ind-1] = flux_U;
          srcs[n_src-1].point_components.power_ref_stokesV[point_ind-1] = flux_V;
          srcs[n_src-1].point_components.power_SIs[point_ind-1] = DEFAULT_SI;
          // LOGV( LOG_LOW, "New source component: <%s> component %d, ra: %8.5f hrs, dec: %+9.5f deg\n",
          //       cat[n_src-1].components[n_comps-1].ra*DR2H, cat[n_src-1].components[n_comps-1].dec*DR2D );
        break;

        case GAUSSIAN:
          srcs[n_src-1].gauss_components.power_ref_freqs[gauss_ind-1] = freq;
          srcs[n_src-1].gauss_components.power_ref_stokesI[gauss_ind-1] = flux_I;
          srcs[n_src-1].gauss_components.power_ref_stokesQ[gauss_ind-1] = flux_Q;
          srcs[n_src-1].gauss_components.power_ref_stokesU[gauss_ind-1] = flux_U;
          srcs[n_src-1].gauss_components.power_ref_stokesV[gauss_ind-1] = flux_V;
          srcs[n_src-1].gauss_components.power_SIs[gauss_ind-1] = DEFAULT_SI;
        break;

        case SHAPELET:
          srcs[n_src-1].shape_components.power_ref_freqs[shape_ind-1] = freq;
          srcs[n_src-1].shape_components.power_ref_stokesI[shape_ind-1] = flux_I;
          srcs[n_src-1].shape_components.power_ref_stokesQ[shape_ind-1] = flux_Q;
          srcs[n_src-1].shape_components.power_ref_stokesU[shape_ind-1] = flux_U;
          srcs[n_src-1].shape_components.power_ref_stokesV[shape_ind-1] = flux_V;
          srcs[n_src-1].shape_components.power_SIs[shape_ind-1] = DEFAULT_SI;
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

          srcs[n_src-1].point_components.power_ref_freqs[point_ind-1] = freq;
          srcs[n_src-1].point_components.power_ref_stokesI[point_ind-1] = flux_I;
          srcs[n_src-1].point_components.power_ref_stokesQ[point_ind-1] = flux_Q;
          srcs[n_src-1].point_components.power_ref_stokesU[point_ind-1] = flux_U;
          srcs[n_src-1].point_components.power_ref_stokesV[point_ind-1] = flux_V;
          srcs[n_src-1].point_components.power_SIs[point_ind-1] = SI;
        break;

        case GAUSSIAN:
          srcs[n_src-1].gauss_components.power_ref_freqs[gauss_ind-1] = freq;
          srcs[n_src-1].gauss_components.power_ref_stokesI[gauss_ind-1] = flux_I;
          srcs[n_src-1].gauss_components.power_ref_stokesQ[gauss_ind-1] = flux_Q;
          srcs[n_src-1].gauss_components.power_ref_stokesU[gauss_ind-1] = flux_U;
          srcs[n_src-1].gauss_components.power_ref_stokesV[gauss_ind-1] = flux_V;
          srcs[n_src-1].gauss_components.power_SIs[gauss_ind-1] = SI;
        break;

        case SHAPELET:
          srcs[n_src-1].shape_components.power_ref_freqs[shape_ind-1] = freq;
          srcs[n_src-1].shape_components.power_ref_stokesI[shape_ind-1] = flux_I;
          srcs[n_src-1].shape_components.power_ref_stokesQ[shape_ind-1] = flux_Q;
          srcs[n_src-1].shape_components.power_ref_stokesU[shape_ind-1] = flux_U;
          srcs[n_src-1].shape_components.power_ref_stokesV[shape_ind-1] = flux_V;
          srcs[n_src-1].shape_components.power_SIs[shape_ind-1] = SI;
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
}
