#include "read_and_write.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include "/usr/local/mwa-RTS_dev/include/slamac.h"
// #include "point_source_lib.h"

source_catalogue_t * read_source_catalogue(char *filename) {
  int result, n_src=0, n_comps=0, n_freqs=0;
  int src_found=0, comp_found=0;
  int src_key, src_end, comp_key, comp_end, freq_key, type_key, param_key, coeff_key;

  char comp_type[16];
  float ra, dec;
  float freq, flux;
  float coeff1, coeff2, coeff3;

  int point_ind, gauss_ind, S1_ind, S2_ind;
  int coeff_ind;

  //So we know how many s_coeffs we have put into S1_coeffs for each source
  // int ind_scoeffs1;

  //Array that will end up being as long as the number of sources,
  //containing catsource_t structs
  catsource_t *srcs=NULL;

  //Will contain the array srcs, and the number of sources n_src
  source_catalogue_t *srccat=NULL;

  //Reading in things
  FILE *fp=NULL;
  char line[BUFSIZ];

  /* open file */
  if ((fp=fopen(filename,"r"))==NULL) {
    printf("read_source_catalogue error: failed to open source file <%s>\n", filename);
    exit(1);
  }

  /* allocate the catalogue object */
  srccat = malloc( sizeof(source_catalogue_t) );
  if (srccat==NULL) {
    printf("read_source_catalogue error: no malloc for srccat\n");
    return NULL;
  }

  //Read in line by line and do "smart" things accordingly
  while(fgets(line,BUFSIZ,fp) != NULL) {

    /* skip blank lines and comments */
    if (line[0]=='\n' || line[0] == '#'|| line[0] == '\0') continue;

    /* lines must start with one of the following keywords... which one is it? */
    src_key = src_end = comp_key = comp_end = freq_key = param_key = coeff_key = 0;
    if(strncmp(line, SRC_KEY, strlen(SRC_KEY))==0) src_key=1;
    if(strncmp(line, SRC_END, strlen(SRC_END))==0) src_end=1;
    if(strncmp(line, COMP_KEY, strlen(COMP_KEY))==0) comp_key=1;
    if(strncmp(line, COMP_END, strlen(COMP_END))==0) comp_end=1;
    if(strncmp(line, FREQ_KEY, strlen(FREQ_KEY))==0) freq_key=1;
    if(strncmp(line, GPARAMS_KEY, strlen(GPARAMS_KEY))==0) param_key=1;
    if(strncmp(line, S1PARAMS_KEY, strlen(S1PARAMS_KEY))==0) param_key=1;
    if(strncmp(line, S2PARAMS_KEY, strlen(S2PARAMS_KEY))==0) param_key=1;
    if(strncmp(line, S1COEFF_KEY, strlen(S1COEFF_KEY))==0) coeff_key=1;
    if(strncmp(line, S2COEFF_KEY, strlen(S2COEFF_KEY))==0) coeff_key=1;


    /* if a source key hasn't been found, skip lines until one is found */
    if ( (src_key==0 && src_end==0 && comp_key==0 && comp_end==0 && freq_key==0 && param_key==0 && coeff_key==0) ||
         (src_found==0 && src_key==0) /* i.e., haven't started a new source */ ) {
        // printf("%c%c%c%c %s %d %d\n",line[0],line[1],line[2],line[3],FREQ_KEY,freq_key,strncmp(line, FREQ_KEY, strlen(FREQ_KEY)));
        printf("read_source_catalogue: skipping bad line: %s\n",line);
        continue;
    }

    if (src_end) {
        if (srcs[n_src-1].n_comps == 0){
            printf("WARNING: source <%s> had no components\n", srcs[n_src-1].name );
            n_src--;
        }
        src_found=0;
        continue;
    }

    /* is this the start of a new source, or more information for the current source? */
    if (src_found && src_key) {
        printf("read_source_catalogue error: found source key in line <%s> while already reading a source\n",
              line );
        return NULL;
    }

    else if (src_key) {

      /* found a new source. make more room, get the details */
      src_found=1;
      n_src++;
      srcs = realloc(srcs,sizeof(catsource_t)*n_src);
      if (srcs == NULL) {
        printf("read_source_catalogue error: no realloc for cat\n");
        return NULL;
      }

      comp_found=0;
      point_ind=0;
      gauss_ind=0;
      S1_ind=0;
      S2_ind=0;
      coeff_ind=0;

      // srcs[n_src-1].name = NULL;
      srcs[n_src-1].n_comps = 0;
      srcs[n_src-1].n_points = 0;
      srcs[n_src-1].n_gauss = 0;
      srcs[n_src-1].n_S1s = 0;
      srcs[n_src-1].n_S1_coeffs = 0;
      srcs[n_src-1].n_S2s = 0;

      //Line should look like <SOURCE name P %d G %d S1 %d %d S2 %d %d>
      //%*s ignores input so we don't have to read in useless input
      result = sscanf( line, "%*s %s %*s %d %*s %d %*s %d %d %*s %d %d", srcs[n_src-1].name,
                        &(srcs[n_src-1].n_points), &(srcs[n_src-1].n_gauss),
                        &(srcs[n_src-1].n_S1s), &(srcs[n_src-1].n_S1_coeffs),
                        &(srcs[n_src-1].n_S2s), &(srcs[n_src-1].n_S2_coeffs) );
      if (result != 7) {
          printf("read_source_catalogue error %d: problem reading cal source input line <%s>\n", result,line );
          return NULL;
      }

      srcs[n_src-1].n_comps = srcs[n_src-1].n_points + srcs[n_src-1].n_gauss + srcs[n_src-1].n_S1s + srcs[n_src-1].n_S2s;

      printf("New source: name: <%s>, Comps:%d  P:%d  G:%d  S1:%d-%d  S2:%d-%d \n",
            srcs[n_src-1].name, srcs[n_src-1].n_comps, srcs[n_src-1].n_points, srcs[n_src-1].n_gauss,
            srcs[n_src-1].n_S1s, srcs[n_src-1].n_S1_coeffs, srcs[n_src-1].n_S2s, srcs[n_src-1].n_S2_coeffs );

      //Now we know the number of sources, do some mallocing
      //Pointsource params
      srcs[n_src-1].point_ras = malloc( srcs[n_src-1].n_points * sizeof(float) );
      srcs[n_src-1].point_decs = malloc( srcs[n_src-1].n_points * sizeof(float) );
      srcs[n_src-1].point_fluxes = malloc( srcs[n_src-1].n_points * sizeof(float) );
      srcs[n_src-1].point_freqs = malloc( srcs[n_src-1].n_points * sizeof(float) );
      //Gaussian params
      srcs[n_src-1].gauss_ras = malloc( srcs[n_src-1].n_gauss * sizeof(float) );
      srcs[n_src-1].gauss_decs = malloc( srcs[n_src-1].n_gauss * sizeof(float) );
      srcs[n_src-1].gauss_fluxes = malloc( srcs[n_src-1].n_gauss * sizeof(float) );
      srcs[n_src-1].gauss_freqs = malloc( srcs[n_src-1].n_gauss * sizeof(float) );
      srcs[n_src-1].gauss_majors = malloc( srcs[n_src-1].n_gauss * sizeof(float) );
      srcs[n_src-1].gauss_minors = malloc( srcs[n_src-1].n_gauss * sizeof(float) );
      srcs[n_src-1].gauss_pas = malloc( srcs[n_src-1].n_gauss * sizeof(float) );

      srcs[n_src-1].S2_ras = malloc( srcs[n_src-1].n_S2s * sizeof(float) );
      srcs[n_src-1].S2_decs = malloc( srcs[n_src-1].n_S2s * sizeof(float) );
      srcs[n_src-1].S2_fluxes = malloc( srcs[n_src-1].n_S2s * sizeof(float) );
      srcs[n_src-1].S2_freqs = malloc( srcs[n_src-1].n_S2s * sizeof(float) );
      srcs[n_src-1].S2_majors = malloc( srcs[n_src-1].n_S2s * sizeof(float) );
      srcs[n_src-1].S2_minors = malloc( srcs[n_src-1].n_S2s * sizeof(float) );
      srcs[n_src-1].S2_pas = malloc( srcs[n_src-1].n_S2s * sizeof(float) );

      //Need bigger arrays to hold the coeffs - make a S2_param_indexes
      //array so we can map each n1.n2,coeff to the corresponding pa,major,minor
      //when we have multiple shapelet models within a source
      srcs[n_src-1].S2_n1s = malloc( srcs[n_src-1].n_S2_coeffs * sizeof(float) );
      srcs[n_src-1].S2_n2s = malloc( srcs[n_src-1].n_S2_coeffs * sizeof(float) );
      srcs[n_src-1].S2_coeffs = malloc( srcs[n_src-1].n_S2_coeffs * sizeof(float) );
      srcs[n_src-1].S2_param_indexes = malloc( srcs[n_src-1].n_S2_coeffs * sizeof(float) );

    }

    else if (comp_key) {

      /* found a new source component. get the details */
      comp_found=1;
      n_comps++;

      result = sscanf( line, "%*s %s %f %f", comp_type, &ra, &dec );
      if (result != 3) {
          printf("CS_ReadSourceList error: problem reading cal component input line <%s>\n", line );
          return NULL;
      }

      if(strncmp(comp_type, POINT_KEY, strlen(POINT_KEY))==0) type_key=POINT;
      if(strncmp(comp_type, GAUSSIAN_KEY, strlen(GAUSSIAN_KEY))==0) type_key=GAUSSIAN;
      if(strncmp(comp_type, SHAPELET_KEY, strlen(SHAPELET_KEY))==0) type_key=SHAPELET;
      if(strncmp(comp_type, SHAPELET2_KEY, strlen(SHAPELET2_KEY))==0) type_key=SHAPELET2;

      switch (type_key) {
        case POINT:
          point_ind++;
          /* convert to radian */
          srcs[n_src-1].point_ras[point_ind-1] = ra * DH2R;
          srcs[n_src-1].point_decs[point_ind-1] = dec * DD2R;
          // LOGV( LOG_LOW, "New source component: <%s> component %d, ra: %8.5f hrs, dec: %+9.5f deg\n",
          //       cat[n_src-1].components[n_comps-1].ra*DR2H, cat[n_src-1].components[n_comps-1].dec*DR2D );
        break;

        case GAUSSIAN:
          gauss_ind++;
          srcs[n_src-1].gauss_ras[gauss_ind-1] = ra * DH2R;
          srcs[n_src-1].gauss_decs[gauss_ind-1] = dec * DD2R;
        break;

        case SHAPELET2:
          S2_ind++;
          srcs[n_src-1].S2_ras[S2_ind-1] = ra * DH2R;
          srcs[n_src-1].S2_decs[S2_ind-1] = dec * DD2R;
        break;

        default:
          printf("CS_ReadSourceList error: unknown source type: %s \n", comp_type );
          return NULL;


      }//switch (type_key)
    }//if (comp_key)

    else if (freq_key && comp_found==1) {

      //I'm ignoring Stokes Q,U,V for now
      result = sscanf( line, "%*s %f %f %*f %*f %*f", &freq, &flux );
      if (result != 2) {
          printf("%d CS_ReadSourceList error: problem reading cal component input line <%s>\n", result, line );
          return NULL;
      }

      switch (type_key) {
        case POINT:

          srcs[n_src-1].point_freqs[point_ind-1] = freq;
          srcs[n_src-1].point_fluxes[point_ind-1] = flux;
          // LOGV( LOG_LOW, "New source component: <%s> component %d, ra: %8.5f hrs, dec: %+9.5f deg\n",
          //       cat[n_src-1].components[n_comps-1].ra*DR2H, cat[n_src-1].components[n_comps-1].dec*DR2D );
        break;

        case GAUSSIAN:
          srcs[n_src-1].gauss_freqs[gauss_ind-1] = freq;
          srcs[n_src-1].gauss_fluxes[gauss_ind-1] = flux;
        break;

        case SHAPELET2:
          srcs[n_src-1].S2_freqs[S2_ind-1] = freq;
          srcs[n_src-1].S2_fluxes[S2_ind-1] = flux;
        break;

      default:
        printf("%d CS_ReadSourceList error assigning info from line: <%s>\n", type_key, line );
        return NULL;

      }//switch (type_key)

    }//if (freq_key && comp_found==1)

    else if (param_key && comp_found==1) {
      result = sscanf( line, "%*s %f %f %f", &coeff1, &coeff2, &coeff3 );
      if (result != 3) {
          printf("%d CS_ReadSourceList error: problem reading coeffs input line <%s>\n", result, line );
          return NULL;
      }

      switch (type_key) {
        case GAUSSIAN:
          srcs[n_src-1].gauss_pas[gauss_ind-1] =  coeff1 * DD2R;
          srcs[n_src-1].gauss_majors[gauss_ind-1] = coeff2 * (DD2R / 60.0);
          srcs[n_src-1].gauss_minors[gauss_ind-1] = coeff3 * (DD2R / 60.0);
        break;

        case SHAPELET2:
          srcs[n_src-1].S2_pas[S2_ind-1] =  coeff1 * DD2R;
          srcs[n_src-1].S2_majors[S2_ind-1] = coeff2 * (DD2R / 60.0);
          srcs[n_src-1].S2_minors[S2_ind-1] = coeff3 * (DD2R / 60.0);
        break;

      default:
        printf("CS_ReadSourceList error assigning info from line: <%s>\n", line );
        return NULL;
      }//switch (type_key)
    }//if (param_key && comp_found==1)

    else if (coeff_key && comp_found==1) {
      result = sscanf( line, "%*s %f %f %f", &coeff1, &coeff2, &coeff3 );
      if (result != 3) {
          printf("%d CS_ReadSourceList error: problem reading coeffs input line <%s>\n", result, line );
          return NULL;
      }

      switch (type_key) {
        case SHAPELET2:
          srcs[n_src-1].S2_n1s[coeff_ind] = coeff1;
          srcs[n_src-1].S2_n2s[coeff_ind] = coeff2;
          srcs[n_src-1].S2_coeffs[coeff_ind] = coeff3;
          srcs[n_src-1].S2_param_indexes[coeff_ind] = S2_ind - 1;
        break;

      default:
        printf("CS_ReadSourceList error assigning info from line: <%s>\n", line );
        return NULL;
      }//switch (type_key)
      coeff_ind++;
    }//if (param_key && comp_found==1)

    else if (comp_end==1) {
      comp_found = 0;
    }

  }//while fgets(line,BUFSIZ,fp) != NULL)

  /* set pointers and counts in the final catalogue object */
  srccat->catsources = srcs;
  srccat->num_sources = n_src;

  return srccat;
} //read_sources


// int main(void)
// {
//
//   float ra0 = 53.62297029690264*D2R;
//   float dec0 = -39.47871177225978*D2R;
//   float lst = 71.64335359166665*D2R;
//
//   // float ra0 = 71.64335359166665*D2R;
//   // float dec0 = -26.7*D2R;
//
//   int num_baselines = 8128;
//
//   float sdec0,cdec0;
//   sdec0 = sin(dec0); cdec0=cos(dec0);
//
//   int num_time_steps = 1;
//   float time_res = 2.0;
//   float ha0,sha0,cha0;
//
//
//   // char sourceinfo[] = "phase_cent_pointsource.txt";
//   // char sourceinfo[] = "dummy_data.txt";
//   // char sourceinfo[] = "simple_ForA.txt";
//   char cat_filename[] = "srclist_wsclean_ForA_phase1_cropped_percent100.txt";
//
//   // int num_components;
//   // num_components = get_num_components(sourceinfo);
//   //
//   // // struct pointsource *pointsources = (pointsource*)malloc(sizeof(struct pointsource) * num_components);
//   // float *comp_ras = (float*)malloc(sizeof(float)*num_components);
//   // float *comp_decs = (float*)malloc(sizeof(float)*num_components);
//   // float *comp_fluxes = (float*)malloc(sizeof(float)*num_components);
//   //
//   // float *comp_pa = (float*)malloc(sizeof(float)*num_components);
//   // float *comp_majors = (float*)malloc(sizeof(float)*num_components);
//   // float *comp_minors = (float*)malloc(sizeof(float)*num_components);
//   //
//   // int *comp_types = (float*)malloc(sizeof(char)*num_components);
//   //
//   // read_pointsources(sourceinfo, comp_ras, comp_decs, comp_fluxes, num_components);
//
//   source_catalogue_t *srccat;
//
//   srccat = read_source_catalogue(cat_filename);
//
//   for (int i = 0; i < 2; i++) {
//     printf("===================================\n");
//     printf("point_ra %d %f\n",i,srccat->catsource[0].point_ras[i]);
//     printf("gauss ra %d %f\n",i,srccat->catsource[0].gauss_ras[i]);
//     printf("point flux %d %f\n",i,srccat->catsource[0].point_fluxes[i]);
//     printf("gauss flux %d %f\n",i,srccat->catsource[0].gauss_fluxes[i]);
//     printf("gauss pa %d %f\n",i,srccat->catsource[0].gauss_pas[i]);
//     printf("point major %d %f\n",i,srccat->catsource[0].gauss_majors[i]);
//     printf("gauss minor %d %f\n",i,srccat->catsource[0].gauss_minors[i]);
//
//     /* code */
//   }
//
// }
