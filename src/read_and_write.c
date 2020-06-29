#include "read_and_write.h"
#include "constants.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <json.h>
#include <fitsio.h>

// enum component_type {POINT=0, GAUSSIAN, SHAPELET};

/*********************************
// Taken from the RTS (Mitchell et al 2008)
// All credit to the original authors
// https://github.com/ICRAR/mwa-RTS.git
  convert coords in local topocentric East, North, Height units to
  'local' XYZ units. Local means Z point north, X points through the equator from the geocenter
  along the local meridian and Y is East.
  This is like the absolute system except that zero lon is now
  the local meridian rather than prime meridian.
  Latitude is geodetic, in radian.
  This is what you want for constructing the local antenna positions in a UVFITS antenna table.
**********************************/
void ENH2XYZ_local(float E, float N, float H, float lat, float *X, float *Y, float *Z) {
  float sl,cl;

  sl = sinf(lat);
  cl = cosf(lat);
  *X = -N*sl + H*cl;
  *Y = E;
  *Z = N*cl + H*sl;
}

/*********************************
// Taken and edited from the RTS (Mitchell et al 2008)
// All credit to the original authors
// https://github.com/ICRAR/mwa-RTS.git
**********************************/
source_catalogue_t * read_source_catalogue(const char *filename) {
  int result, n_src=0, n_comps=0, n_freqs=0;
  int src_found=0, comp_found=0;
  int src_key, src_end, comp_key, comp_end, freq_key, type_key, param_key, coeff_key;

  char comp_type[16];
  float ra, dec;
  float freq, flux;
  float coeff1, coeff2, coeff3;

  int point_ind, gauss_ind, shape_ind;
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
    if(strncmp(line, SPARAMS_KEY, strlen(SPARAMS_KEY))==0) param_key=1;
    if(strncmp(line, SCOEFF_KEY, strlen(SCOEFF_KEY))==0) coeff_key=1;

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
          printf("read_source_catalogue error %d: problem reading cal source input line <%s>\n", result,line );
          return NULL;
      }

      srcs[n_src-1].n_comps = srcs[n_src-1].n_points + srcs[n_src-1].n_gauss + srcs[n_src-1].n_shapes;

      // printf("New source: name: <%s>, Comps:%d  P:%d  G:%d  S:%d-%d \n",
      //       srcs[n_src-1].name, srcs[n_src-1].n_comps, srcs[n_src-1].n_points, srcs[n_src-1].n_gauss,
      //       srcs[n_src-1].n_shapes, srcs[n_src-1].n_shape_coeffs );

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

      srcs[n_src-1].shape_ras = malloc( srcs[n_src-1].n_shapes * sizeof(float) );
      srcs[n_src-1].shape_decs = malloc( srcs[n_src-1].n_shapes * sizeof(float) );
      srcs[n_src-1].shape_fluxes = malloc( srcs[n_src-1].n_shapes * sizeof(float) );
      srcs[n_src-1].shape_freqs = malloc( srcs[n_src-1].n_shapes * sizeof(float) );
      srcs[n_src-1].shape_majors = malloc( srcs[n_src-1].n_shapes * sizeof(float) );
      srcs[n_src-1].shape_minors = malloc( srcs[n_src-1].n_shapes * sizeof(float) );
      srcs[n_src-1].shape_pas = malloc( srcs[n_src-1].n_shapes * sizeof(float) );

      //Need bigger arrays to hold the coeffs - make a shape_param_indexes
      //array so we can map each n1.n2,coeff to the corresponding pa,major,minor
      //when we have multiple shapelet models within a source
      srcs[n_src-1].shape_n1s = malloc( srcs[n_src-1].n_shape_coeffs * sizeof(float) );
      srcs[n_src-1].shape_n2s = malloc( srcs[n_src-1].n_shape_coeffs * sizeof(float) );
      srcs[n_src-1].shape_coeffs = malloc( srcs[n_src-1].n_shape_coeffs * sizeof(float) );
      srcs[n_src-1].shape_param_indexes = malloc( srcs[n_src-1].n_shape_coeffs * sizeof(float) );

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

        case SHAPELET:
          shape_ind++;
          srcs[n_src-1].shape_ras[shape_ind-1] = ra * DH2R;
          srcs[n_src-1].shape_decs[shape_ind-1] = dec * DD2R;
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

        case SHAPELET:
          srcs[n_src-1].shape_freqs[shape_ind-1] = freq;
          srcs[n_src-1].shape_fluxes[shape_ind-1] = flux;
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

        case SHAPELET:
          srcs[n_src-1].shape_pas[shape_ind-1] =  coeff1 * DD2R;
          srcs[n_src-1].shape_majors[shape_ind-1] = coeff2 * (DD2R / 60.0);
          srcs[n_src-1].shape_minors[shape_ind-1] = coeff3 * (DD2R / 60.0);
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
        case SHAPELET:
          srcs[n_src-1].shape_n1s[coeff_ind] = coeff1;
          srcs[n_src-1].shape_n2s[coeff_ind] = coeff2;
          srcs[n_src-1].shape_coeffs[coeff_ind] = coeff3;
          srcs[n_src-1].shape_param_indexes[coeff_ind] = shape_ind - 1;
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

woden_settings_t * read_json_settings(const char *filename){
  FILE *fp;
	char buffer[1024];

  struct json_object *parsed_json;
  struct json_object *lst_base;
  struct json_object *ra0;
  struct json_object *dec0;
  struct json_object *num_baselines;
  struct json_object *num_freqs;
  struct json_object *frequency_resolution;
  struct json_object *base_frequency;
  struct json_object *num_time_steps;
  struct json_object *time_res;
  struct json_object *cat_filename;
  struct json_object *metafits_filename;
  struct json_object *sky_crop_type;
  struct json_object *gaussian_beam;
  struct json_object *chunking_size;
  struct json_object *gauss_beam_FWHM;
  struct json_object *gauss_beam_ref_freq;
  struct json_object *FEE_beam;

	fp = fopen(filename,"r");
	fread(buffer, 1024, 1, fp);
	fclose(fp);

	parsed_json = json_tokener_parse(buffer);

  // json_object_object_get_ex(parsed_json, "lst_base", &lst_base);
  json_object_object_get_ex(parsed_json, "ra0", &ra0);
  json_object_object_get_ex(parsed_json, "dec0", &dec0);
  // json_object_object_get_ex(parsed_json, "num_baselines", &num_baselines);
  json_object_object_get_ex(parsed_json, "num_freqs", &num_freqs);
  json_object_object_get_ex(parsed_json, "frequency_resolution", &frequency_resolution);
  // json_object_object_get_ex(parsed_json, "base_frequency", &base_frequency);
  json_object_object_get_ex(parsed_json, "num_time_steps", &num_time_steps);
  json_object_object_get_ex(parsed_json, "time_res", &time_res);
  json_object_object_get_ex(parsed_json, "cat_filename", &cat_filename);
  json_object_object_get_ex(parsed_json, "metafits_filename", &metafits_filename);
  json_object_object_get_ex(parsed_json, "sky_crop_components", &sky_crop_type);
  json_object_object_get_ex(parsed_json, "use_gaussian_beam", &gaussian_beam);
  json_object_object_get_ex(parsed_json, "chunking_size", &chunking_size);

  json_object_object_get_ex(parsed_json, "gauss_beam_FWHM", &gauss_beam_FWHM);
  json_object_object_get_ex(parsed_json, "gauss_beam_ref_freq", &gauss_beam_ref_freq);

  json_object_object_get_ex(parsed_json, "use_FEE_beam", &FEE_beam);

  woden_settings_t * woden_settings;
  // woden_settings = NULL;
  woden_settings = malloc( sizeof(woden_settings_t) );

  // woden_settings->lst_base = (float)json_object_get_double(lst_base)*D2R;
  woden_settings->ra0 = (float)json_object_get_double(ra0)*DD2R;
  woden_settings->dec0 = (float)json_object_get_double(dec0)*DD2R;
  // woden_settings->num_baselines = json_object_get_int(num_baselines);
  woden_settings->num_freqs = json_object_get_int(num_freqs);
  woden_settings->frequency_resolution = (float)json_object_get_double(frequency_resolution);
  // woden_settings->base_frequency = (float)json_object_get_double(base_frequency);
  woden_settings->num_time_steps = json_object_get_int(num_time_steps);
  woden_settings->time_res = (float)json_object_get_double(time_res);
  woden_settings->cat_filename = json_object_get_string(cat_filename);
  woden_settings->metafits_filename = json_object_get_string(metafits_filename);

  //Boolean setting whether to crop sources by SOURCE or by COMPONENT
  woden_settings->sky_crop_type = json_object_get_boolean(sky_crop_type);

  //Boolean whether to use gaussian primary beam
  int gauss_beam = json_object_get_boolean(gaussian_beam);

  //Boolean whether to use gaussian primary beam
  int fee_beam = json_object_get_boolean(FEE_beam);

  if (gauss_beam) {
    woden_settings->beamtype = GAUSS_BEAM;

    float beam_FWHM = (float)json_object_get_double(gauss_beam_FWHM);
    //If the gauss_beam_FWHM has been set in the json file, use it
    //Otherwise, set the defult FWHM of 20 deg
    if (beam_FWHM > 0.0) {
      woden_settings->gauss_beam_FWHM = beam_FWHM;
    }
    else {
      woden_settings->gauss_beam_FWHM = 20.0;
    }

    float beam_ref_freq = (float)json_object_get_double(gauss_beam_ref_freq);
    //If gauss_beam_ref_freq has been set in the json file, use it
    //Otherwise, set the defult FWHM of 20 deg
    if (beam_ref_freq > 0.0) {
      woden_settings->gauss_beam_ref_freq = beam_ref_freq;
    }
    else {
      woden_settings->gauss_beam_ref_freq = 150e+6;
    }
  }
  //TODO once we implement MWA primary beam, check whether both gaussian beam and
  //MWA beam have been selected. If so, tell the user, and exit, get them to sort
  //themselves out

  //TODO make sure metafits contains correct beam info if using FEE_beam

  else if (fee_beam){
    woden_settings->beamtype = FEE_BEAM;
  }

  else {
    woden_settings->beamtype = NO_BEAM;
  }

  woden_settings->chunking_size = json_object_get_int(chunking_size);


  struct json_object *band_num;
  struct json_object *band_nums;
  size_t num_bands;
	size_t i;

  json_object_object_get_ex(parsed_json, "band_nums", &band_nums);
  num_bands = json_object_array_length(band_nums);

  woden_settings->num_bands = num_bands;
  woden_settings->band_nums = malloc( num_bands * sizeof(int));

	for(i=0;i<num_bands;i++) {
		band_num = json_object_array_get_idx(band_nums, i);
		woden_settings->band_nums[i] = json_object_get_int(band_num);
	}



  return woden_settings;

}

/*********************************
// Taken and edited from the RTS (Mitchell et al 2008)
// All credit to the original authors
// https://github.com/ICRAR/mwa-RTS.git
**********************************/
int RTS_init_meta_file(fitsfile *mfptr, MetaFfile_t *metafits, const char *nome){
    int status=0;
    int ncols, anynulls, colnum, nfound, i;
    long naxes[2], nrows, frow, felem;
    //char tflg[1];
    float nullval;
    // FILE *fplog=NULL;
    // extra variable for printing
    char card[FLEN_CARD], keyname[FLEN_KEYWORD], coltype[FLEN_VALUE], colname[FLEN_VALUE];
    int single = 0, nkeys, hdupos, ii;

    // fplog=LogGetFilehandle();
    // if (fplog==NULL) fplog=stderr;

    memset(metafits,'\0',sizeof(MetaFfile_t));

    // prints out metafits header contents
    int look = 0;
    if (look){
      fits_get_hdu_num(mfptr, &hdupos);
      for (; !status; hdupos++)  /* Main loop through each extension */
      {
        fits_get_hdrspace(mfptr, &nkeys, NULL, &status); /* get # of keywords */
        printf("Header listing for HDU #%d:\n", hdupos);
        for (ii = 1; ii <= nkeys; ii++) { /* Read and print each keywords */
           if (fits_read_record(mfptr, ii, card, &status))break;
           printf("%s\n", card);
        }
        printf("END\n\n");  /* terminate listing with END */
        if (single) break;  /* quit if only listing a single header */
        fits_movrel_hdu(mfptr, 1, NULL, &status);  /* try to move to next HDU */
      }
    } // ends look

    fits_read_key(mfptr,TSTRING, "VERSION", &(metafits->version), NULL, &status);
    if (status) {
        printf("No VERSION keyword in metafits. Continuing...\n");
        fits_clear_errmsg();
        status=0;
    }

    fits_read_key(mfptr, TSTRING, "MWAVER", &(metafits->mwaVersion), NULL, &status);
    if (status) {
        printf("No MWAVER keyword in metafits. Continuing...\n");
        fits_clear_errmsg();
        status=0;
    }
    fits_read_key(mfptr, TSTRING, "RECVRS", &(metafits->recvrs), NULL, &status);
    if (status) {
        // fits_report_error(fplog,status);
        printf("No RECVRS keyword in metafits.\n");
        return status;
    }
    fits_read_key(mfptr,TSTRING, "CALIBSRC", &(metafits->calsrc), NULL, &status);
    if (status) {
        printf("No CALIBSRC keyword in metafits. Continuing...\n");
        fits_clear_errmsg();
        status=0;
    }
    fits_read_key(mfptr,TSTRING, "TILEFLAG", &(metafits->tileflg), NULL, &status);
    if (status) {
       printf("No TILEFLAG keyword in metafits. Continuing...\n");
        fits_clear_errmsg();
        status=0;
    }

    fits_read_key(mfptr,TINT, "CENTCHAN", &(metafits->centchan), NULL, &status);
    if (status) {
       printf("No CENTCHAN keyword in metafits. Continuing...\n");
        fits_clear_errmsg();
        status=0;
    }

    fits_read_key(mfptr,TFLOAT, "LST", &(metafits->lst_base), NULL, &status);
    if (status) {
       printf("No LST keyword in metafits. Continuing...\n");
        fits_clear_errmsg();
        status=0;
    }

    fits_read_key(mfptr,TFLOAT, "RA", &(metafits->ra_point), NULL, &status);
    if (status) {
       printf("No RA keyword in metafits. Continuing...\n");
        fits_clear_errmsg();
        status=0;
    }

    fits_read_key(mfptr,TFLOAT, "DEC", &(metafits->dec_point), NULL, &status);
    if (status) {
       printf("No DEC keyword in metafits. Continuing...\n");
        fits_clear_errmsg();
        status=0;
    }

    fits_read_key(mfptr,TFLOAT, "FINECHAN", &(metafits->frequency_resolution), NULL, &status);
    if (status) {
       printf("No FINECHAN keyword in metafits. Continuing...\n");
        fits_clear_errmsg();
        status=0;
    }

    fits_read_key(mfptr,TFLOAT, "FREQCENT", &(metafits->frequency_cent), NULL, &status);
    if (status) {
       printf("No FREQCENT keyword in metafits. Continuing...\n");
        fits_clear_errmsg();
        status=0;
    }

    fits_read_key(mfptr,TFLOAT, "BANDWDTH", &(metafits->bandwidth), NULL, &status);
    if (status) {
       printf("No BANDWDTH keyword in metafits. Continuing...\n");
        fits_clear_errmsg();
        status=0;
    }

    fits_read_key(mfptr,TFLOAT, "INTTIME", &(metafits->time_res), NULL, &status);
    if (status) {
       printf("No INTTIME keyword in metafits. Continuing...\n");
        fits_clear_errmsg();
        status=0;
    }

    fits_read_key(mfptr,TINT, "NINPUTS", &(metafits->num_tiles), NULL, &status);
    if (status) {
       printf("No NINPUTS keyword in metafits. Continuing...\n");
        fits_clear_errmsg();
        status=0;
    }

    //NINPUTS has both XX and YY inputs into correlator so divide by 2 to
    //get the number of tiles
    metafits->num_tiles = metafits->num_tiles / 2;
    metafits->lst_base *= DD2R;
    metafits->ra_point *= DD2R;
    metafits->dec_point *= DD2R;
    metafits->frequency_resolution *= 1e+3;
    metafits->frequency_cent *= 1e+6;
    metafits->bandwidth *= 1e+6;
    metafits->base_low_freq = metafits->frequency_cent - (metafits->bandwidth / 2.0) - (metafits->frequency_resolution / 2.0);

    /* now move onto HDU2, which is where the binary table is */
  int hdutype=0;
  fits_movrel_hdu(mfptr,1 , &hdutype, &status);
  if (status) {
      // fits_report_error(fplog,status);
      printf("%s: Cannot move to binary table HDU.",__func__);
      return status;
  }
  // printf("%s: moved to next HDU. Type is: %d\n",__func__,hdutype);
  fits_read_keys_lng(mfptr, "NAXIS", 1, 2, naxes, &nfound, &status);
  (metafits->naxes[0]) = naxes[0]; (metafits->naxes[1]) = naxes[1];
  fits_get_num_rows(mfptr, &nrows, &status);
  fits_get_num_cols(mfptr, &ncols, &status);
  //print for checking stuff (default disabled)
  if (look) {
    for (ii = 1; ii <= ncols; ii++) {
      fits_make_keyn("TTYPE", ii, keyname, &status); /* make keyword */
      fits_read_key(mfptr, TSTRING, keyname, colname, NULL, &status);
      fits_make_keyn("TFORM", ii, keyname, &status); /* make keyword */
      fits_read_key(mfptr, TSTRING, keyname, coltype, NULL, &status);
      printf(" %3d %-16s %-16s\n", ii, colname, coltype);
    }
  }
  frow = 1;
  felem = 1;
  nullval = -99.;
  //read Input key
  fits_get_colnum(mfptr, CASEINSEN, "Input", &colnum, &status);
  fits_read_col(mfptr, TINT, colnum,frow, felem, metafits->naxes[1], &nullval, &(metafits->inps), &anynulls, &status);
  //read Antenna key
  fits_get_colnum(mfptr, CASEINSEN, "Antenna", &colnum, &status);
  fits_read_col(mfptr, TINT, colnum,frow, felem, metafits->naxes[1], &nullval, &(metafits->ants), &anynulls, &status);
  // read Tile key
  fits_get_colnum(mfptr, CASEINSEN, "Tile", &colnum, &status);
  fits_read_col(mfptr, TINT, colnum,frow, felem, metafits->naxes[1], &nullval, &(metafits->tile), &anynulls, &status);
  // read Dipole Delays
  fits_get_colnum(mfptr, CASEINSEN, "Delays", &colnum, &status);
  for (i=0; i<metafits->naxes[1]; i++){
     fits_read_col(mfptr, TINT, colnum,i+1, felem, 16, &nullval, &(metafits->FEE_delays[i]),
                   &anynulls, &status);
  }
// read Digital Gains
  fits_get_colnum(mfptr, CASEINSEN, "Gains", &colnum, &status);
  for (i=0; i<metafits->naxes[1]; i++){
     fits_read_col(mfptr, TINT, colnum,i+1, felem, 24, &nullval, &(metafits->dig_gains[i]),
                   &anynulls, &status);
  }

  // read Flags
  fits_get_colnum(mfptr, CASEINSEN, "Flag", &colnum, &status);
  fits_read_col(mfptr, TINT, colnum,frow, felem, metafits->naxes[1], &nullval, &(metafits->flags), &anynulls, &status);
  // read Antenna Pos E
  fits_get_colnum(mfptr, CASEINSEN, "East", &colnum, &status);
  fits_read_col(mfptr, TFLOAT, colnum,frow, felem, metafits->naxes[1], &nullval, &(metafits->E), &anynulls, &status);
  // read Antenna Pos N
  fits_get_colnum(mfptr, CASEINSEN, "North", &colnum, &status);
  fits_read_col(mfptr, TFLOAT, colnum,frow, felem, metafits->naxes[1], &nullval, &(metafits->N), &anynulls, &status);
  // read Antenna Pos H
  fits_get_colnum(mfptr, CASEINSEN, "Height", &colnum, &status);
  fits_read_col(mfptr, TFLOAT, colnum,frow, felem, metafits->naxes[1], &nullval, &(metafits->H), &anynulls, &status);
  // read cable Length
  fits_get_colnum(mfptr, CASEINSEN, "Length", &colnum, &status);
        //for (i=0; i<metafits->naxes[1]; i++){
	char *lengptr[metafits->naxes[1]];
	//following for loop written as a test solution for seg faults in reading lengths (James)
	for(i=0; i<metafits->naxes[1];i++){
            lengptr[i] = metafits->leng[i];
	}

  fits_read_col(mfptr, TSTRING, colnum, frow, felem, metafits->naxes[1], &nullval, &lengptr, &anynulls, &status);
        //}
	for(i=0;i<metafits->naxes[1];i++){
		metafits->leng[i][0] = 'E';
	}
  return status;
}

array_layout_t * calc_XYZ_diffs(MetaFfile_t *metafits, int num_tiles){

  array_layout_t * array_layout;
  array_layout = malloc( sizeof(array_layout_t) );

  array_layout->latitude = MWA_LAT*DD2R;
  array_layout->num_tiles = num_tiles;
  array_layout->num_baselines = (array_layout->num_tiles*(array_layout->num_tiles-1)) / 2;

  array_layout->ant_east = malloc( array_layout->num_tiles * sizeof(float) );
  array_layout->ant_north = malloc( array_layout->num_tiles * sizeof(float) );
  array_layout->ant_height = malloc( array_layout->num_tiles * sizeof(float) );

  array_layout->ant_X = malloc( array_layout->num_tiles * sizeof(float) );
  array_layout->ant_Y = malloc( array_layout->num_tiles * sizeof(float) );
  array_layout->ant_Z = malloc( array_layout->num_tiles * sizeof(float) );
  //
  for (int i = 0; i < array_layout->num_tiles; i++) {
    //Metafits e,n,h goes XX,YY,XX,YY so need to choose every other value
    array_layout->ant_east[i] = metafits->E[i*2];
    array_layout->ant_north[i] = metafits->N[i*2];
    array_layout->ant_height[i] = metafits->H[i*2];

    //Convert to local X,Y,Z
    ENH2XYZ_local(array_layout->ant_east[i], array_layout->ant_north[i], array_layout->ant_height[i],
                  array_layout->latitude,
                  &(array_layout->ant_X[i]), &(array_layout->ant_Y[i]), &(array_layout->ant_Z[i]));

  }

  array_layout->X_diff_metres = malloc( array_layout->num_baselines * sizeof(float) );
  array_layout->Y_diff_metres = malloc( array_layout->num_baselines * sizeof(float) );
  array_layout->Z_diff_metres = malloc( array_layout->num_baselines * sizeof(float) );

  int baseline_ind = 0;
  for (int ant1 = 0; ant1 < array_layout->num_tiles - 1; ant1++) {
    for (int ant2 = ant1 + 1; ant2 < array_layout->num_tiles; ant2++) {
      array_layout->X_diff_metres[baseline_ind] = array_layout->ant_X[ant1] - array_layout->ant_X[ant2];
      array_layout->Y_diff_metres[baseline_ind] = array_layout->ant_Y[ant1] - array_layout->ant_Y[ant2];
      array_layout->Z_diff_metres[baseline_ind] = array_layout->ant_Z[ant1] - array_layout->ant_Z[ant2];
      baseline_ind++;
    }
  }

  return array_layout;
//
}
