// #include <read_yaml_skymodel.h>
#include <math.h>
#include <stdio.h>
#include "create_sky_model.h"

user_precision_t zero_key_get_value_user_precision(char *event_value, int * key_int) {

  * key_int = 0;

  user_precision_t result;
  #ifdef DOUBLE_PRECISION
    result = strtod((char *)event_value, NULL);
  #else
    result = strtof((char *)event_value, NULL);
  #endif

  return (user_precision_t)result;
}

void check_key(char *event_value, char *key_char, int * key_int) {

  if (strcmp(event_value, key_char) == 0) {
    * key_int = 1;
  }
}

void update_realloc_num(int * num_malloc){
  if (* num_malloc == 0) {
    * num_malloc = INITIAL_NUM_COMPONENTS;
  }
  else {
    * num_malloc *= 2;
  }
}


int add_component_information(components_t *components, int num_srcs,
          int comp_ind, int power_flux_ind, int curve_flux_ind,
          double ra, double dec, double freq,
          e_flux_type fluxtype, e_component_type comptype,
          user_precision_t flux_I, user_precision_t flux_Q,
          user_precision_t flux_U, user_precision_t flux_V,
          user_precision_t curve_q,
          user_precision_t si, user_precision_t pa,
          user_precision_t maj, user_precision_t min,
          track_comp_malloc_t *track_comp_malloc){

  // printf("On component %d\n",comp_ind );

  if (comp_ind + 1 > track_comp_malloc->n_comps) {

    // printf("Before realloc comp %d\n", track_comp_malloc->n_comps);

    update_realloc_num(&track_comp_malloc->n_comps);
    // printf("After realloc comp %d\n", track_comp_malloc->n_comps);

    components->ras = realloc(components->ras,
                                     sizeof(double)*track_comp_malloc->n_comps);
    components->decs = realloc(components->decs,
                                     sizeof(double)*track_comp_malloc->n_comps);

    if (comptype == GAUSSIAN || comptype == SHAPELET) {

      components->majors = realloc(components->majors,
                           sizeof(user_precision_t)*track_comp_malloc->n_comps);
      components->minors = realloc(components->minors,
                           sizeof(user_precision_t)*track_comp_malloc->n_comps);
      components->pas = realloc(components->pas,
                           sizeof(user_precision_t)*track_comp_malloc->n_comps);
    }
  }

  components->ras[comp_ind] = ra*DD2R;
  components->decs[comp_ind] = dec*DD2R;

  if (comptype == GAUSSIAN || comptype == SHAPELET) {
    components->majors[comp_ind] = maj * (DD2R / 3600.0);
    components->minors[comp_ind] = min * (DD2R / 3600.0);
    components->pas[comp_ind] = pa * DD2R;
  }

  if (isnan(ra) != 0) {
    printf("Didn't find an RA, it's still NAN\n");
    return 1;
  }

  if (fluxtype == POWER_LAW) {

    if (isnan(flux_I) != 0 && isnan(flux_Q) != 0 && isnan(flux_U) != 0 && isnan(flux_V) != 0) {
      printf("Didn't find either I, Q, U, V, when reading a POWER_LAW\n");
      return 1;
    }
    if (isnan(freq) != 0) {
      printf("Didn't find reference freq when reading a POWER_LAW\n");
      return 1;
    }


    if (power_flux_ind + 1 > track_comp_malloc->n_powers) {

      update_realloc_num(&track_comp_malloc->n_powers);

      components->power_ref_freqs = realloc(components->power_ref_freqs,
                                     sizeof(double)*track_comp_malloc->n_powers);
      components->power_ref_stokesI = realloc(components->power_ref_stokesI,
                          sizeof(user_precision_t)*track_comp_malloc->n_powers);
      components->power_ref_stokesQ = realloc(components->power_ref_stokesQ,
                          sizeof(user_precision_t)*track_comp_malloc->n_powers);
      components->power_ref_stokesU = realloc(components->power_ref_stokesU,
                          sizeof(user_precision_t)*track_comp_malloc->n_powers);
      components->power_ref_stokesV = realloc(components->power_ref_stokesV,
                          sizeof(user_precision_t)*track_comp_malloc->n_powers);
      components->power_SIs = realloc(components->power_SIs,
                          sizeof(user_precision_t)*track_comp_malloc->n_powers);
      components->power_comp_inds = realloc(components->power_comp_inds,
                                       sizeof(int)*track_comp_malloc->n_powers);
    }

    components->power_ref_freqs[power_flux_ind] = freq;
    components->power_ref_stokesI[power_flux_ind] = flux_I;
    components->power_ref_stokesQ[power_flux_ind] = flux_Q;
    components->power_ref_stokesU[power_flux_ind] = flux_U;
    components->power_ref_stokesV[power_flux_ind] = flux_V;
    components->power_SIs[power_flux_ind] = si;
    components->power_comp_inds[power_flux_ind] = comp_ind;

  }

  else if ( fluxtype == CURVED_POWER_LAW) {

    if (isnan(flux_I) != 0 && isnan(flux_Q) != 0 && isnan(flux_U) != 0 && isnan(flux_V) != 0) {
      printf("Didn't find either I, Q, U, V, when reading a CURVED_POWER_LAW\n");
      return 1;
    }
    if (isnan(freq) != 0) {
      printf("Didn't find reference freq when reading a CURVED_POWER_LAW\n");
      return 1;
    }

    if (curve_flux_ind + 1 > track_comp_malloc->n_curves) {

      update_realloc_num(&track_comp_malloc->n_curves);

      components->curve_ref_freqs = realloc(components->curve_ref_freqs,
                                   sizeof(double)*track_comp_malloc->n_curves);
      components->curve_ref_stokesI = realloc(components->curve_ref_stokesI,
                          sizeof(user_precision_t)*track_comp_malloc->n_curves);
      components->curve_ref_stokesQ = realloc(components->curve_ref_stokesQ,
                          sizeof(user_precision_t)*track_comp_malloc->n_curves);
      components->curve_ref_stokesU = realloc(components->curve_ref_stokesU,
                          sizeof(user_precision_t)*track_comp_malloc->n_curves);
      components->curve_ref_stokesV = realloc(components->curve_ref_stokesV,
                          sizeof(user_precision_t)*track_comp_malloc->n_curves);
      components->curve_SIs = realloc(components->curve_SIs,
                          sizeof(user_precision_t)*track_comp_malloc->n_curves);
      components->curve_qs = realloc(components->curve_qs,
                          sizeof(user_precision_t)*track_comp_malloc->n_curves);
      components->curve_comp_inds = realloc(components->curve_comp_inds,
                                       sizeof(int)*track_comp_malloc->n_curves);
    }

    components->curve_ref_freqs[curve_flux_ind] = freq;
    components->curve_ref_stokesI[curve_flux_ind] = flux_I;
    components->curve_ref_stokesQ[curve_flux_ind] = flux_Q;
    components->curve_ref_stokesU[curve_flux_ind] = flux_U;
    components->curve_ref_stokesV[curve_flux_ind] = flux_V;
    components->curve_SIs[curve_flux_ind] = si;
    components->curve_qs[curve_flux_ind] = curve_q;
    components->curve_comp_inds[curve_flux_ind] = comp_ind;

  }

  return 0;
}

// Function to swap elements
void swap_int(int *a, int *b) {
  int temp = *a;
  *a = *b;
  *b = temp;
}

void swap_double(double *a, double *b) {
  double temp = *a;
  *a = *b;
  *b = temp;
}

// bubble sort function
void bubbleSort_freqs_index_array(int *index_array, double *freqs,  int num_freqs) {

  for (int index = 0; index < num_freqs; index++) {
    index_array[index] = index;
  }

  int i, j;
  for (i = 0; i < num_freqs-1; i++) {
    for (j = 0; j < num_freqs-i-1; j++) {
      if (freqs[j] > freqs[j+1]) {
          swap_double(&freqs[j], &freqs[j+1]);
          swap_int(&index_array[j], &index_array[j+1]);
      }
    }
  }
}

void add_list_flux_component_info(components_t *components, int num_list_entries,
                                 int comp_ind, int comp_list_ind,
                                 double *list_freqs_1D,
                                 user_precision_t *list_stokesI_1D,
                                 user_precision_t *list_stokesQ_1D,
                                 user_precision_t *list_stokesU_1D,
                                 user_precision_t *list_stokesV_1D,
                                 track_comp_malloc_t *track_comp_malloc){

  //sort the frequencies so they are in ascending order. This makes
  //exptrapolated frequencies easier down the line
  int *index_array = malloc(num_list_entries*sizeof(int));

  bubbleSort_freqs_index_array(index_array, list_freqs_1D, num_list_entries);

  if (comp_list_ind+1 > track_comp_malloc->n_lists) {

    update_realloc_num(&track_comp_malloc->n_lists);

    // printf("Reallocing a set of list, now have on component %d\n",
    //        track_comp_malloc->n_lists);

    components->num_list_values = realloc(components->num_list_values,
                                        track_comp_malloc->n_lists*sizeof(int));
    components->list_start_indexes = realloc(components->list_start_indexes,
                                        track_comp_malloc->n_lists*sizeof(int));
    components->list_comp_inds = realloc(components->list_comp_inds,
                                    track_comp_malloc->n_lists*sizeof(int));
  }

  int low_list_index;

  if (comp_list_ind == 0) {
    low_list_index = 0;
  } else {
    low_list_index = components->list_start_indexes[comp_list_ind - 1] + components->num_list_values[comp_list_ind - 1];
  }

  components->list_comp_inds[comp_list_ind] = comp_ind - 1;
  components->list_start_indexes[comp_list_ind] = low_list_index;
  components->num_list_values[comp_list_ind] = num_list_entries;

  int new_size = low_list_index + num_list_entries;

  if (new_size > track_comp_malloc->n_list_values){

    // printf("Reallocing a the big-ass list array, now %d\n", track_comp_malloc->n_list_values);

    update_realloc_num(&track_comp_malloc->n_list_values);

    components->list_freqs = realloc(components->list_freqs,
                       track_comp_malloc->n_list_values*sizeof(double));
    components->list_stokesI = realloc(components->list_stokesI,
                       track_comp_malloc->n_list_values*sizeof(user_precision_t));
    components->list_stokesQ = realloc(components->list_stokesQ,
                       track_comp_malloc->n_list_values*sizeof(user_precision_t));
    components->list_stokesU = realloc(components->list_stokesU,
                       track_comp_malloc->n_list_values*sizeof(user_precision_t));
    components->list_stokesV = realloc(components->list_stokesV,
                       track_comp_malloc->n_list_values*sizeof(user_precision_t));
  }

  for (int i = 0; i < num_list_entries; i++) {
    components->list_freqs[low_list_index + i] = list_freqs_1D[i];
    components->list_stokesI[low_list_index + i] = list_stokesI_1D[index_array[i]];
    components->list_stokesQ[low_list_index + i] = list_stokesQ_1D[index_array[i]];
    components->list_stokesU[low_list_index + i] = list_stokesU_1D[index_array[i]];
    components->list_stokesV[low_list_index + i] = list_stokesV_1D[index_array[i]];
  }
  free(index_array);
}

track_comp_malloc_t * initialise_track_comp_malloc(void){

  track_comp_malloc_t *track_comp_malloc = malloc(sizeof(track_comp_malloc_t));

  track_comp_malloc->n_comps = 0;
  track_comp_malloc->n_powers = 0;
  track_comp_malloc->n_curves = 0;
  track_comp_malloc->n_lists = 0;
  track_comp_malloc->n_list_values = 0;
  track_comp_malloc->n_shape_coeffs = 0;

  return track_comp_malloc;
}

void realloc_component_to_min_size(components_t *comps, int num_powers,
                                   int num_curves, int num_lists,
                                   int num_shape_coeffs,
                                   e_component_type comptype){

  int n_comps = num_powers +  num_curves + num_lists;

  comps->ras = realloc(comps->ras, sizeof(double)*n_comps);
  comps->decs = realloc(comps->decs, sizeof(double)*n_comps);

  if (comptype == GAUSSIAN || comptype == SHAPELET){
    comps->majors = realloc(comps->majors, sizeof(user_precision_t)*n_comps);
    comps->minors = realloc(comps->minors, sizeof(user_precision_t)*n_comps);
    comps->pas = realloc(comps->pas, sizeof(user_precision_t)*n_comps);
  }

  if (comptype == SHAPELET){
    comps->shape_coeffs = realloc(comps->shape_coeffs,
                                   num_shape_coeffs*sizeof(user_precision_t));
    comps->n1s = realloc(comps->n1s,
                                   num_shape_coeffs*sizeof(user_precision_t));
    comps->n2s = realloc(comps->n2s,
                                   num_shape_coeffs*sizeof(user_precision_t));
    comps->param_indexes = realloc(comps->param_indexes,
                                   num_shape_coeffs*sizeof(user_precision_t));
  }

  if (num_powers > 0){
      comps->power_ref_freqs = realloc(comps->power_ref_freqs,
                                          sizeof(double)*num_powers);
      comps->power_ref_stokesI = realloc(comps->power_ref_stokesI,
                                          sizeof(user_precision_t)*num_powers);
      comps->power_ref_stokesQ = realloc(comps->power_ref_stokesQ,
                                          sizeof(user_precision_t)*num_powers);
      comps->power_ref_stokesU = realloc(comps->power_ref_stokesU,
                                          sizeof(user_precision_t)*num_powers);
      comps->power_ref_stokesV = realloc(comps->power_ref_stokesV,
                                          sizeof(user_precision_t)*num_powers);
      comps->power_SIs = realloc(comps->power_SIs,
                                          sizeof(user_precision_t)*num_powers);
      comps->power_comp_inds = realloc(comps->power_comp_inds,
                                          sizeof(int)*num_powers);
  }

  if (num_curves > 0) {
    comps->curve_ref_freqs = realloc(comps->curve_ref_freqs,
                                          sizeof(double)*num_curves);
    comps->curve_ref_stokesI = realloc(comps->curve_ref_stokesI,
                                          sizeof(user_precision_t)*num_curves);
    comps->curve_ref_stokesQ = realloc(comps->curve_ref_stokesQ,
                                          sizeof(user_precision_t)*num_curves);
    comps->curve_ref_stokesU = realloc(comps->curve_ref_stokesU,
                                          sizeof(user_precision_t)*num_curves);
    comps->curve_ref_stokesV = realloc(comps->curve_ref_stokesV,
                                          sizeof(user_precision_t)*num_curves);
    comps->curve_SIs = realloc(comps->curve_SIs,
                                          sizeof(user_precision_t)*num_curves);
    comps->curve_qs = realloc(comps->curve_qs,
                                          sizeof(user_precision_t)*num_curves);
    comps->curve_comp_inds = realloc(comps->curve_comp_inds,
                                          sizeof(int)*num_curves);
  }

  if (num_lists > 0) {
    comps->list_comp_inds = realloc(comps->list_comp_inds,
      sizeof(int)*num_lists);
    comps->num_list_values = realloc(comps->num_list_values,
      sizeof(int)*num_lists);
    comps->list_start_indexes = realloc(comps->list_start_indexes,
      sizeof(int)*num_lists);

    comps->total_num_flux_entires = comps->list_start_indexes[num_lists-1] + comps->num_list_values[num_lists-1];

    comps->list_freqs = realloc(comps->list_freqs,
                        sizeof(double)*comps->total_num_flux_entires);
    comps->list_stokesI = realloc(comps->list_stokesI,
                        sizeof(user_precision_t)*comps->total_num_flux_entires);
    comps->list_stokesQ = realloc(comps->list_stokesQ,
                        sizeof(user_precision_t)*comps->total_num_flux_entires);
    comps->list_stokesU = realloc(comps->list_stokesU,
                        sizeof(user_precision_t)*comps->total_num_flux_entires);
    comps->list_stokesV = realloc(comps->list_stokesV,
                        sizeof(user_precision_t)*comps->total_num_flux_entires);
  }

}


int read_yaml_skymodel(const char *yaml_path, source_catalogue_t *srccat)
{
  FILE *fh = fopen(yaml_path, "r");
  yaml_parser_t parser;
  yaml_event_t  event;   /* New variable */

  /* Initialize parser */
  if(!yaml_parser_initialize(&parser)) {
    printf("Failed to initialize yaml parser\n");
    return 1;
  }
  if(fh == NULL) {
    printf("Failed to open file: %s\n", yaml_path);
    return 1;
  }

  /* Set input file */
  yaml_parser_set_input_file(&parser, fh);

  char* RA_KEY = "ra";
  char* DEC_KEY = "dec";
  char* COMP_HYP_KEY = "comp_type";

  //For power law / curved power law COMPONENTS
  char* FLUX_TYPE_KEY = "flux_type";
  char* SI_KEY = "si";
  char* REF_FLUX_KEY = "fd";
  char* REF_FREQ_KEY = "freq";
  char* I_KEY = "i";
  char* Q_KEY = "q";
  char* U_KEY = "u";
  char* V_KEY = "v";

  //For shapelet/gaussian COMPONENTs
  char* MAJ_KEY = "maj";
  char* MIN_KEY = "min";
  char* PA_KEY = "pa";

  char* COEFF_KEY = "coeffs";
  char* N1_KEY = "n1";
  char* N2_KEY = "n2";
  char* VALUE_KEY = "value";

  int ra_key = 0, dec_key = 0, comp_key = 0, flux_type_key = 0;
  int si_key = 0, freq_key = 0, ref_flux_key = 0;
  int i_key = 0, q_key = 0, u_key = 0, v_key = 0;

  int maj_key = 0, min_key = 0, pa_key = 0;
  int coeff_key = 0, n1_key = 0, n2_key = 0, value_key = 0;

  //Values we'll be grabbing out the catalogue
  double ra = NAN, dec=NAN, freq=NAN;
  user_precision_t si=NAN;
  user_precision_t flux_I=NAN;
  user_precision_t flux_Q=NAN, flux_U=NAN, flux_V=NAN;
  user_precision_t curve_q=NAN;

  user_precision_t maj = NAN, min = NAN, pa = NAN;
  user_precision_t n2 = NAN;
  user_precision_t n1 = NAN;
  user_precision_t shape_coeff_val = NAN;


  //Used to work out where we are in the process of reading things in
  int started_source = 0;
  //Everytime this a new field is read in, we get a YAML_MAPPING_START_EVENT
  // int started_mapping = 0;
  int num_mappings = 0;

  // int started_freq_list_entry = 0;
  int num_list_entries = 0;

  source_t *srcs=NULL;
  int num_srcs = 0;

  e_component_type comptype = -1;
  e_flux_type fluxtype = -1;
  int list_key = 0;

  //Various indexes because I'm mental - these things keep track of each
  //COMPONENT type and the flux model each uses
  int point_ind = 0;
  int point_list_ind = 0;

  int gauss_ind = 0;
  int gauss_list_ind = 0;

  int shape_ind = 0;
  int shape_list_ind = 0;
  int shape_coeff_ind = 0;

  // int n_shape_coeffs_incremented = 0;

  srccat->num_shapelets = 0;

  //Arrays to keep LIST style information in
  double *list_freqs_1D = malloc(sizeof(double)*INITIAL_NUM_FLUXES);
  user_precision_t *list_stokesI_1D = malloc(sizeof(user_precision_t)*INITIAL_NUM_FLUXES);
  user_precision_t *list_stokesQ_1D = malloc(sizeof(user_precision_t)*INITIAL_NUM_FLUXES);
  user_precision_t *list_stokesU_1D = malloc(sizeof(user_precision_t)*INITIAL_NUM_FLUXES);
  user_precision_t *list_stokesV_1D = malloc(sizeof(user_precision_t)*INITIAL_NUM_FLUXES);
  int num_list_malloc = INITIAL_NUM_FLUXES;


  track_comp_malloc_t *track_point_malloc = initialise_track_comp_malloc();
  track_comp_malloc_t *track_gauss_malloc = initialise_track_comp_malloc();
  track_comp_malloc_t *track_shape_malloc = initialise_track_comp_malloc();

  //Loop through the whole yaml object, and keep going until
  //you find the end of the YAML stream
  while(event.type != YAML_STREAM_END_EVENT)
  {
    //Wig out if the something fails
    if (!yaml_parser_parse(&parser, &event)) {
       printf("Parser error %d\n", parser.error);
       exit(EXIT_FAILURE);
    }

    switch(event.type)
    {
    case YAML_NO_EVENT: puts("No event!"); break;
    /* Stream start/end */
    case YAML_STREAM_START_EVENT: puts("Started reading yaml"); break;
    case YAML_STREAM_END_EVENT:   puts("Finsihed reading yaml");   break;
    /* Block delimeters */
    case YAML_DOCUMENT_START_EVENT: break; //puts("Start Document");
    case YAML_DOCUMENT_END_EVENT:   break; //puts("End Document");


    //WHAT is an alias event??
    case YAML_ALIAS_EVENT:
      printf("Got alias (anchor %s)\n", event.data.alias.anchor);
    break;

    //Alright, this means we have a new SOURCE
    case YAML_SEQUENCE_START_EVENT:
      // puts("Start Sequence");
      if (started_source == 0 && coeff_key == 0) {
        // printf("STARTED A SOURCE NOW\n");
        num_srcs++;
        srcs = realloc(srcs,sizeof(source_t)*num_srcs);

        //We gonna be realloc a bunch of arrays, so set them to NULL to
        //start with. ALso set all the counter int numbers to 0.0
        source_zero_counters_and_null_components(&srcs[num_srcs-1]);

        started_source = 1;

      }
    break;

    //If we've started SOURCE, then start counting the mappings, which will
    //tell us how many fields we've filled in (i.e. ra, dec, comp_type)
    case YAML_MAPPING_START_EVENT:
      if (started_source == 0) {
        num_mappings = 0;
      }
      else if (list_key == 1){
        //If this is a list frequency type thing, we don't want to
        //increment the mapping here
        ;
      }
      else if (coeff_key == 1){
        //If this is a shapelet coeff type thing, we don't want to
        //increment the mapping here
        ;
      }
      else {
        num_mappings += 1;
      }
      // printf("Start mapping, num maps %d\n",num_mappings );
    break;

    //Ok, here we've either found a header (i.e ra, dec, comp_type)
    //or the associated values. Either set a key so we know in the next
    //iteration to collect the data, or actually collect the data
    case YAML_SCALAR_EVENT:
    // printf("%s\n", event.data.scalar.value );

      //************************************************************************
      //These if statements mean we previously set a key, and are ready to
      //read in certain values. Reset the key once value is read
      if (ra_key == 1){
          ra = strtod((char *)event.data.scalar.value, NULL);
          // printf("The string read of RA is %s\n", event.data.scalar.value );
          // printf("THE RA IS %.32f\n", ra);
          //Reset the RA key, done with it now
          ra_key = 0;
      }
      else if (dec_key == 1){
          dec = strtod((char *)event.data.scalar.value, NULL);
          // printf("The string read of DEC is %s\n", event.data.scalar.value );
          // printf("THE DEC IS %.32f\n", dec);
          //Reset the RA key, done with it now
          dec_key = 0;
      }

      else if (comp_key == 1){

        if (strcmp((char *)event.data.scalar.value, "point") == 0) {
          comptype = POINT;
        }
        else if (strcmp((char *)event.data.scalar.value, "gaussian") == 0) {
          comptype = GAUSSIAN;
        }
        else if (strcmp((char *)event.data.scalar.value, "shapelet") == 0) {
          comptype = SHAPELET;
          // printf("WE FOUND A SHAPELET\n");
        }
        comp_key = 0;
      }

      else if (flux_type_key == 1){
        if (strcmp((char *)event.data.scalar.value, "power_law") == 0) {
          fluxtype = POWER_LAW;
        }
        else if (strcmp((char *)event.data.scalar.value, "curved_power_law") == 0) {
          fluxtype = CURVED_POWER_LAW;
        }
        else if (strcmp((char *)event.data.scalar.value, "list") == 0) {
          fluxtype = LIST;
          // printf("WE FOUND A LIST FLUX TYPE\n");
          list_key = 1;
        }
        else {
          printf("The `flux_type` key was neither 'power_law', 'curved_power_law', or 'list'.\n");
          return 1;
        }
        flux_type_key = 0;
      }

      else if (si_key == 1){
        si = zero_key_get_value_user_precision((char *)event.data.scalar.value, &si_key);
      }
      else if (freq_key == 1){
        freq = zero_key_get_value_user_precision((char *)event.data.scalar.value, &freq_key);
      }
      else if (i_key == 1){
        flux_I = zero_key_get_value_user_precision((char *)event.data.scalar.value, &i_key);
      }
      //Something special here - there is a 'q' for Stokes, and a 'q' for spectral
      //curvature (sigh). So update different values depending on whether the
      //flux related keys are set or not
      else if (q_key == 1){
        if (ref_flux_key == 1 || list_key == 1) {
          flux_Q = zero_key_get_value_user_precision((char *)event.data.scalar.value, &q_key);
        } else {
          curve_q = zero_key_get_value_user_precision((char *)event.data.scalar.value, &q_key);
        }


      }
      else if (u_key == 1){
        flux_U = zero_key_get_value_user_precision((char *)event.data.scalar.value, &u_key);
      }
      else if (v_key == 1){
        flux_V = zero_key_get_value_user_precision((char *)event.data.scalar.value, &v_key);
      }
      else if (maj_key == 1){
        maj = zero_key_get_value_user_precision((char *)event.data.scalar.value, &maj_key);
      }
      else if (min_key == 1){
        min = zero_key_get_value_user_precision((char *)event.data.scalar.value, &min_key);
      }
      else if (pa_key == 1){
        pa = zero_key_get_value_user_precision((char *)event.data.scalar.value, &pa_key);
      }
      else if (n2_key == 1){
        n2 = zero_key_get_value_user_precision((char *)event.data.scalar.value, &n2_key);
      }
      else if (n1_key == 1){
        n1 = zero_key_get_value_user_precision((char *)event.data.scalar.value, &n1_key);
      }
      else if (value_key == 1){
        // printf("The string read of value is %s\n", event.data.scalar.value );
        shape_coeff_val = zero_key_get_value_user_precision((char *)event.data.scalar.value, &value_key);
      }

      // //This is separate from the big "else if" above because we leave
      // //the coeff_key switched on until a "sequence" ends
      // if (coeff_key == 1){
      //   if (n_shape_coeffs_incremented == 0) {
      //     srcs[num_srcs-1].n_shape_coeffs += 1;
      //     n_shape_coeffs_incremented = 1;
      //   }
      // }
      //************************************************************************

      //************************************************************************
      //Here, we have no keys set, so read the event value to set a key
      //and work out what value we are grabbing next
      //Go through the possible events we are looking for, and setup
      //keys to instruct the next interation in the while loop to grab certain
      //values

      check_key((char *)event.data.scalar.value, RA_KEY, &ra_key);
      check_key((char *)event.data.scalar.value, DEC_KEY, &dec_key);
      check_key((char *)event.data.scalar.value, COMP_HYP_KEY, &comp_key);
      check_key((char *)event.data.scalar.value, FLUX_TYPE_KEY, &flux_type_key);
      check_key((char *)event.data.scalar.value, SI_KEY, &si_key);
      check_key((char *)event.data.scalar.value, REF_FREQ_KEY, &freq_key);
      check_key((char *)event.data.scalar.value, I_KEY, &i_key);
      check_key((char *)event.data.scalar.value, Q_KEY, &q_key);
      check_key((char *)event.data.scalar.value, U_KEY, &u_key);
      check_key((char *)event.data.scalar.value, V_KEY, &v_key);
      check_key((char *)event.data.scalar.value, MAJ_KEY, &maj_key);
      check_key((char *)event.data.scalar.value, MIN_KEY, &min_key);
      check_key((char *)event.data.scalar.value, PA_KEY, &pa_key);
      check_key((char *)event.data.scalar.value, COEFF_KEY, &coeff_key);
      check_key((char *)event.data.scalar.value, N1_KEY, &n1_key);
      check_key((char *)event.data.scalar.value, N2_KEY, &n2_key);
      check_key((char *)event.data.scalar.value, VALUE_KEY, &value_key);
      check_key((char *)event.data.scalar.value, REF_FLUX_KEY, &ref_flux_key);

      //************************************************************************

    break; //break for YAML_SCALAR_EVENT

    //Decrement the number of mappings if we finish a mapping
    //If we get to a mapping number of zero, it means we've finished this component
    case YAML_MAPPING_END_EVENT:

      if (started_source == 0) {
        num_mappings = 0;
      }
      else if (list_key == 1){

        if (num_list_entries + 1 > num_list_malloc){

          num_list_malloc *= 2;

          list_freqs_1D = realloc(list_freqs_1D, num_list_malloc*sizeof(double));
          list_stokesI_1D = realloc(list_stokesI_1D, num_list_malloc*sizeof(user_precision_t));
          list_stokesQ_1D = realloc(list_stokesQ_1D, num_list_malloc*sizeof(user_precision_t));
          list_stokesU_1D = realloc(list_stokesU_1D, num_list_malloc*sizeof(user_precision_t));
          list_stokesV_1D = realloc(list_stokesV_1D, num_list_malloc*sizeof(user_precision_t));

        }

        list_freqs_1D[num_list_entries] = freq;
        list_stokesI_1D[num_list_entries] = flux_I;

        // printf("MALLOCing this many things %d %.1e\n",(num_list_entries+1),freq, flux_Q );

        if (isnan(flux_Q) != 0) {
          list_stokesQ_1D[num_list_entries] = 0.0;
        } else {
          list_stokesQ_1D[num_list_entries] = flux_Q;
        }
        if (isnan(flux_U) != 0) {
          list_stokesU_1D[num_list_entries] = 0.0;
        } else {
          list_stokesU_1D[num_list_entries] = flux_U;
        }
        if (isnan(flux_V) != 0) {
          list_stokesV_1D[num_list_entries] = 0.0;
        } else {
          list_stokesV_1D[num_list_entries] = flux_V;
        }

        num_list_entries += 1;
        flux_I = NAN, flux_Q = NAN, flux_U = NAN, flux_V = NAN;
      } else if (coeff_key == 1){

        srcs[num_srcs - 1].n_shape_coeffs += 1;

        if (srcs[num_srcs - 1].n_shape_coeffs > track_shape_malloc->n_shape_coeffs) {

          update_realloc_num(&track_shape_malloc->n_shape_coeffs);

          srcs[num_srcs - 1].shape_components.n1s = realloc(srcs[num_srcs - 1].shape_components.n1s,
                     sizeof(user_precision_t)*track_shape_malloc->n_shape_coeffs);
          srcs[num_srcs - 1].shape_components.n2s = realloc(srcs[num_srcs - 1].shape_components.n2s,
                     sizeof(user_precision_t)*track_shape_malloc->n_shape_coeffs);
          srcs[num_srcs - 1].shape_components.shape_coeffs = realloc(srcs[num_srcs - 1].shape_components.shape_coeffs,
                     sizeof(user_precision_t)*track_shape_malloc->n_shape_coeffs);
          srcs[num_srcs - 1].shape_components.param_indexes = realloc(srcs[num_srcs - 1].shape_components.param_indexes,
                     sizeof(user_precision_t)*track_shape_malloc->n_shape_coeffs);
        }

        srcs[num_srcs - 1].shape_components.n1s[shape_coeff_ind] = n1;
        srcs[num_srcs - 1].shape_components.n2s[shape_coeff_ind] = n2;
        srcs[num_srcs - 1].shape_components.shape_coeffs[shape_coeff_ind] = shape_coeff_val;
        srcs[num_srcs - 1].shape_components.param_indexes[shape_coeff_ind] = shape_ind;

        n1 = NAN;
        n2 = NAN;
        shape_coeff_val = NAN;
        // printf("INSIDE READING\n %d %d\n",srcs[num_srcs-1].n_shape_coeffs, shape_coeff_ind);
        // printf("%.1e %.1e %.1e %d\n",n1, n2, shape_coeff_val,  shape_ind);
        // n_shape_coeffs_incremented = 0;
        shape_coeff_ind += 1;
      }
      else {
        num_mappings -= 1;
      }

      //If we were collecting reference flux values, this means we've finished
      //We can reset the ref_flux_key at this point
      if (ref_flux_key == 1){
        // printf("There was a mapping for the ref flux key ending?\n");
        ref_flux_key = 0;
      }

      //If we have completed all the mappings in this component, add
      //the information into the current src
      if (started_source == 1 && num_mappings == 0){
        if (comptype == POINT) {

          srcs[num_srcs-1].n_points += 1;

          add_component_information(&srcs[num_srcs - 1].point_components, num_srcs,
            point_ind, srcs[num_srcs - 1].n_point_powers, srcs[num_srcs - 1].n_point_curves,
            ra, dec, freq, fluxtype, comptype,
            flux_I, flux_Q, flux_U, flux_V, curve_q, si, pa, maj, min,
            track_point_malloc);

          point_ind += 1;
        }
        else if (comptype == GAUSSIAN) {
          srcs[num_srcs-1].n_gauss += 1;

          add_component_information(&srcs[num_srcs - 1].gauss_components, num_srcs,
            gauss_ind, srcs[num_srcs - 1].n_gauss_powers, srcs[num_srcs - 1].n_gauss_curves,
            ra, dec, freq, fluxtype, comptype,
            flux_I, flux_Q, flux_U, flux_V, curve_q, si, pa, maj, min,
            track_gauss_malloc);

          gauss_ind += 1;
        }
        else if (comptype == SHAPELET) {
          srcs[num_srcs-1].n_shapes += 1;

          add_component_information(&srcs[num_srcs - 1].shape_components, num_srcs,
            shape_ind, srcs[num_srcs - 1].n_shape_powers, srcs[num_srcs - 1].n_shape_curves,
            ra, dec, freq, fluxtype, comptype,
            flux_I, flux_Q, flux_U, flux_V, curve_q, si, pa, maj, min,
            track_shape_malloc);

          shape_ind += 1;
        } else{
          printf("The `comp_type` wasn't set. Check 'comp_type' spelling, and that it's value is either 'point', 'gaussian' or 'shapelet'\n");
          return 1;
        }

        if (fluxtype == LIST) {
          if (comptype == POINT) {

            add_list_flux_component_info(&srcs[num_srcs - 1].point_components,
                                             num_list_entries,
                                             point_ind, point_list_ind,
                                             list_freqs_1D,
                                             list_stokesI_1D,
                                             list_stokesQ_1D,
                                             list_stokesU_1D,
                                             list_stokesV_1D,
                                             track_point_malloc);

            srcs[num_srcs - 1].n_point_lists += 1;
            point_list_ind += 1;

          } else if (comptype == GAUSSIAN) {


            add_list_flux_component_info(&srcs[num_srcs - 1].gauss_components,
                                             num_list_entries,
                                             gauss_ind, gauss_list_ind,
                                             list_freqs_1D,
                                             list_stokesI_1D,
                                             list_stokesQ_1D,
                                             list_stokesU_1D,
                                             list_stokesV_1D,
                                             track_gauss_malloc);

            srcs[num_srcs - 1].n_gauss_lists += 1;
            gauss_list_ind += 1;
          } else if (comptype == SHAPELET) {

            add_list_flux_component_info(&srcs[num_srcs - 1].shape_components,
                                             num_list_entries,
                                             shape_ind, shape_list_ind,
                                             list_freqs_1D,
                                             list_stokesI_1D,
                                             list_stokesQ_1D,
                                             list_stokesU_1D,
                                             list_stokesV_1D,
                                             track_shape_malloc);

            srcs[num_srcs - 1].n_shape_lists += 1;
            shape_list_ind += 1;
          }
          //Reset for the next LIST flux type component
          num_list_entries = 0;
          // free(list_freqs_1D);
          // free(list_stokesI_1D);
          // free(list_stokesQ_1D);
          // free(list_stokesU_1D);
          // free(list_stokesV_1D);
          //realloc doesn't like it if I don't add a NULL in here. Good to make
          //sure we free and NULL so we don't copy old flux into new components
          // list_freqs_1D = NULL;
          // list_stokesI_1D = NULL;
          // list_stokesQ_1D = NULL;
          // list_stokesU_1D = NULL;
          // list_stokesV_1D = NULL;
        }

        else if (fluxtype == POWER_LAW) {
            if (comptype == POINT) {
              srcs[num_srcs - 1].n_point_powers += 1;
            } else if (comptype == GAUSSIAN) {
              srcs[num_srcs - 1].n_gauss_powers += 1;
            } else if (comptype == SHAPELET) {
              srcs[num_srcs - 1].n_shape_powers += 1;
            }

        }
        else if (fluxtype == CURVED_POWER_LAW) {
            if (comptype == POINT) {
              srcs[num_srcs - 1].n_point_curves += 1;
            } else if (comptype == GAUSSIAN) {
              srcs[num_srcs - 1].n_gauss_curves += 1;
            } else if (comptype == SHAPELET) {
              srcs[num_srcs - 1].n_shape_curves += 1;
            }

        }

        //Once new component information has been added, reset a number of values
        comptype = -1;
        fluxtype = -1;
        ra = NAN, dec = NAN, freq = NAN;
        flux_I = NAN, flux_Q = NAN, flux_U = NAN, flux_V = NAN;
        si = NAN;
        maj = NAN, min = NAN, pa = NAN;
        n2 = NAN, n1 = NAN, shape_coeff_val = NAN;
        curve_q = NAN;

      }
    break;

    case YAML_SEQUENCE_END_EVENT:
      // printf("End Sequence\n");

      //Ok in this event, our "sequence" was adding one set of shapelet
      //basis function information. So need to realloc to get correct space,
      //insert the new information, and then reset keys, coeff counting, and increment
      //the index of the coefficients
      if (coeff_key == 1) {

        coeff_key = 0;

      }

      else if (list_key == 1) {
        //The entries for a LIST style flux have ended, so switch key to zero

        list_key = 0;
      }

      //Otherwise, we have a SOURCE that is ending, so wrap up and reset SOURCE
      //related couting
      else {
        started_source = 0;
        point_ind = 0;
        point_list_ind = 0;

        gauss_ind = 0;
        gauss_list_ind = 0;

        shape_ind = 0;
        shape_list_ind = 0;
        shape_coeff_ind = 0;

        // n_shape_coeffs_incremented = 0;

        srccat->num_shapelets += srcs[num_srcs-1].n_shapes;
        // printf("THIS %d %d %d\n",srcs[num_srcs-1].n_points, srcs[num_srcs-1].n_gauss,  srcs[num_srcs - 1].n_shapes);
        srcs[num_srcs - 1].n_comps = srcs[num_srcs-1].n_points + srcs[num_srcs-1].n_gauss + srcs[num_srcs - 1].n_shapes;

        realloc_component_to_min_size(&srcs[num_srcs - 1].point_components,
                                              srcs[num_srcs - 1].n_point_powers,
                                              srcs[num_srcs - 1].n_point_curves,
                                              srcs[num_srcs - 1].n_point_lists,
                                              srcs[num_srcs - 1].n_shape_coeffs,
                                              POINT);

        realloc_component_to_min_size(&srcs[num_srcs - 1].gauss_components,
                                              srcs[num_srcs - 1].n_gauss_powers,
                                              srcs[num_srcs - 1].n_gauss_curves,
                                              srcs[num_srcs - 1].n_gauss_lists,
                                              srcs[num_srcs - 1].n_shape_coeffs,
                                              GAUSSIAN);

        realloc_component_to_min_size(&srcs[num_srcs - 1].shape_components,
                                              srcs[num_srcs - 1].n_shape_powers,
                                              srcs[num_srcs - 1].n_shape_curves,
                                              srcs[num_srcs - 1].n_shape_lists,
                                              srcs[num_srcs - 1].n_shape_coeffs,
                                              SHAPELET);
        //Reset components malloc counters for next SOURCE
        track_point_malloc->n_comps = 0;
        track_point_malloc->n_powers = 0;
        track_point_malloc->n_curves = 0;
        track_point_malloc->n_lists = 0;
        track_point_malloc->n_shape_coeffs = 0;
        track_point_malloc->n_list_values = 0;

        track_gauss_malloc->n_comps = 0;
        track_gauss_malloc->n_powers = 0;
        track_gauss_malloc->n_curves = 0;
        track_gauss_malloc->n_lists = 0;
        track_gauss_malloc->n_shape_coeffs = 0;
        track_gauss_malloc->n_list_values = 0;

        track_shape_malloc->n_comps = 0;
        track_shape_malloc->n_powers = 0;
        track_shape_malloc->n_curves = 0;
        track_shape_malloc->n_lists = 0;
        track_shape_malloc->n_shape_coeffs = 0;
        track_shape_malloc->n_list_values = 0;

      }

    break;

    }//end of switching over events
    //Unless this is the end of the STREAM, we have to delete the event to make
    //space for the next event
    if(event.type != YAML_STREAM_END_EVENT) {
      yaml_event_delete(&event);
    }
  }//end of while loop over events
  yaml_event_delete(&event);

  /* Cleanup */
  yaml_parser_delete(&parser);
  fclose(fh);

  //Have finished reading in everything now, so fill in some overall catalogue
  //information
  srccat->num_sources = num_srcs;
  srccat->sources = srcs;

  free(list_freqs_1D);
  free(list_stokesI_1D);
  free(list_stokesQ_1D);
  free(list_stokesU_1D);
  free(list_stokesV_1D);

  return 0;
}
