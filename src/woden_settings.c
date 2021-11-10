/*******************************************************************************
*  Methods to read in simulation parameters from a .json file and prepare
*  settings into a `woden_settings_t` struct
*  @author J.L.B. Line
*
*  Please see documentation in ../include/read_and_write.h or online at
*  https://woden.readthedocs.io/en/latest/index.html
*******************************************************************************/
#include "woden_settings.h"
#include "constants.h"
#include "woden_precision_defs.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <json.h>
#include <pal.h>

int read_json_settings(const char *filename,  woden_settings_t *woden_settings){
  FILE *fp;
	char buffer[1024];

  struct json_object *parsed_json;
  struct json_object *lst_base;
  struct json_object *latitude;
  struct json_object *ra0;
  struct json_object *dec0;
  struct json_object *num_freqs;
  struct json_object *frequency_resolution;
  struct json_object *base_low_freq;
  struct json_object *num_time_steps;
  struct json_object *time_res;
  struct json_object *cat_filename;
  struct json_object *coarse_band_width;
  struct json_object *sky_crop_type;
  struct json_object *gaussian_beam;
  struct json_object *chunking_size;
  struct json_object *gauss_beam_FWHM;
  struct json_object *gauss_beam_ref_freq;
  struct json_object *FEE_beam;
  struct json_object *hdf5_beam_path;
  struct json_object *jd_date;
  struct json_object *EDA2_beam;
  struct json_object *array_layout_file_path;
  struct json_object *no_precession;
  //
  /* open file */
  if ((fp=fopen(filename,"r"))==NULL) {
    printf("read_json_settings: failed to open json file:\n\t %s \n", filename);
    // exit(1);
    return 1;
  }
  //
  size_t status = fread(buffer, 1024, 1, fp);
	fclose(fp);

  if (status == 1) {
    printf("Failed to read data from the given .json file stream\n");
    return 1;
  }
  //
	parsed_json = json_tokener_parse(buffer);
  //
  json_object_object_get_ex(parsed_json, "LST", &lst_base);
  json_object_object_get_ex(parsed_json, "latitude", &latitude);
  json_object_object_get_ex(parsed_json, "ra0", &ra0);
  json_object_object_get_ex(parsed_json, "dec0", &dec0);
  json_object_object_get_ex(parsed_json, "num_freqs", &num_freqs);
  json_object_object_get_ex(parsed_json, "frequency_resolution", &frequency_resolution);
  json_object_object_get_ex(parsed_json, "lowest_channel_freq", &base_low_freq);
  json_object_object_get_ex(parsed_json, "coarse_band_width", &coarse_band_width);

  json_object_object_get_ex(parsed_json, "num_time_steps", &num_time_steps);
  json_object_object_get_ex(parsed_json, "time_res", &time_res);
  json_object_object_get_ex(parsed_json, "cat_filename", &cat_filename);
  json_object_object_get_ex(parsed_json, "sky_crop_components", &sky_crop_type);

  json_object_object_get_ex(parsed_json, "chunking_size", &chunking_size);
  json_object_object_get_ex(parsed_json, "array_layout", &array_layout_file_path);
  json_object_object_get_ex(parsed_json, "jd_date", &jd_date);

  json_object_object_get_ex(parsed_json, "use_gaussian_beam", &gaussian_beam);
  json_object_object_get_ex(parsed_json, "gauss_beam_FWHM", &gauss_beam_FWHM);
  json_object_object_get_ex(parsed_json, "gauss_beam_ref_freq", &gauss_beam_ref_freq);

  json_object_object_get_ex(parsed_json, "use_FEE_beam", &FEE_beam);
  json_object_object_get_ex(parsed_json, "hdf5_beam_path", &hdf5_beam_path);

  json_object_object_get_ex(parsed_json, "use_EDA2_beam", &EDA2_beam);
  json_object_object_get_ex(parsed_json, "no_precession", &no_precession);
  //
  //
  //Boolean whether to use gaussian primary beam
  // int lat_true = json_object_get_boolean(latitude);

  //See whether latitude has been set or not
  int lat_true = json_object_get_type(latitude);

  if (lat_true == 0) {
    woden_settings->latitude = MWA_LAT_RAD;
  }
  else {
    woden_settings->latitude = json_object_get_double(latitude)*DD2R;
  }

  woden_settings->lst_base = json_object_get_double(lst_base)*DD2R;
  woden_settings->ra0 = (user_precision_t)json_object_get_double(ra0)*DD2R;
  woden_settings->dec0 = (user_precision_t)json_object_get_double(dec0)*DD2R;
  woden_settings->num_freqs = json_object_get_int(num_freqs);
  woden_settings->frequency_resolution = (user_precision_t)json_object_get_double(frequency_resolution);

  woden_settings->base_low_freq = (user_precision_t)json_object_get_double(base_low_freq);
  woden_settings->coarse_band_width = (user_precision_t)json_object_get_double(coarse_band_width);

  woden_settings->num_time_steps = json_object_get_int(num_time_steps);
  woden_settings->time_res = (user_precision_t)json_object_get_double(time_res);
  woden_settings->cat_filename = json_object_get_string(cat_filename);
  woden_settings->hdf5_beam_path = json_object_get_string(hdf5_beam_path);
  woden_settings->jd_date = (user_precision_t)json_object_get_double(jd_date);

  //Boolean setting whether to crop sources by SOURCE or by COMPONENT
  woden_settings->sky_crop_type = json_object_get_boolean(sky_crop_type);

  //Boolean whether to use gaussian primary beam
  int gauss_beam = json_object_get_boolean(gaussian_beam);

  //Boolean whether to use the MWA FEE beam
  int fee_beam = json_object_get_boolean(FEE_beam);

  //Boolean whether to use the MWA FEE beam
  int eda2_beam = json_object_get_boolean(EDA2_beam);

  //Boolean whether to use apply precession to array or not
  int no_precess = json_object_get_boolean(no_precession);

  if (no_precess) {
    woden_settings->do_precession = 0;
  } else {
    woden_settings->do_precession = 1;
  }

  if (gauss_beam + fee_beam + eda2_beam > 1 ) {
    printf("You have selected more than one primary beam type in the .json file\n");
    printf("You can have only ONE of the following:\n");
    printf("\t\"use_gaussian_beam\": True\n");
    printf("\t\"use_FEE_beam\": True\n");
    printf("\t\"use_EDA2_beam\": True\n");
    return 1;
  }

  if (gauss_beam) {
    woden_settings->beamtype = GAUSS_BEAM;

    user_precision_t beam_FWHM = (user_precision_t)json_object_get_double(gauss_beam_FWHM);
    //If the gauss_beam_FWHM has been set in the json file, use it
    //Otherwise, set the defult FWHM of 20 deg
    if (beam_FWHM > 0.0) {
      woden_settings->gauss_beam_FWHM = beam_FWHM;
    } else {
      woden_settings->gauss_beam_FWHM = 20.0;
    }

    user_precision_t beam_ref_freq = (user_precision_t)json_object_get_double(gauss_beam_ref_freq);
    //If gauss_beam_ref_freq has been set in the json file, use it
    //Otherwise, set the reference to 150e+6
    if (beam_ref_freq > 0.0) {
      woden_settings->gauss_beam_ref_freq = beam_ref_freq;
    } else {
      woden_settings->gauss_beam_ref_freq = 150e+6;
    }

    struct json_object *gauss_ra_point;
    struct json_object *gauss_dec_point;
    json_object_object_get_ex(parsed_json, "gauss_ra_point", &gauss_ra_point);
    json_object_object_get_ex(parsed_json, "gauss_dec_point", &gauss_dec_point);
    woden_settings->gauss_ra_point = (user_precision_t)json_object_get_double(gauss_ra_point)*DD2R;
    woden_settings->gauss_dec_point = (user_precision_t)json_object_get_double(gauss_dec_point)*DD2R;

  }
  else if (fee_beam){
    woden_settings->beamtype = FEE_BEAM;

    struct json_object *delay;
    struct json_object *FEE_ideal_delays;
    int delays_length;

    json_object_object_get_ex(parsed_json, "FEE_delays", &FEE_ideal_delays);
    delays_length = json_object_array_length(FEE_ideal_delays);

    if (delays_length != 16) {
      printf("FEE_delays in json file must be an array of length 16");
      return 1;
    }

  	for(int i=0;i<delays_length;i++) {
  		delay = json_object_array_get_idx(FEE_ideal_delays, i);
  		woden_settings->FEE_ideal_delays[i] = (user_precision_t)json_object_get_double(delay);
  	}

    woden_settings->hdf5_beam_path = json_object_get_string(hdf5_beam_path);

    if (woden_settings->hdf5_beam_path == NULL) {
      printf("Must specify path the MWA FEE hdf5 file for MWA FEE Beam simulation \n");
      return 1;
    }

  }

  else if (EDA2_beam){
    woden_settings->beamtype = ANALY_DIPOLE;
  }

  else {
    woden_settings->beamtype = NO_BEAM;
  }

  woden_settings->chunking_size = json_object_get_int64(chunking_size);
  //If user selects an insanely large chunking size gonna have problems, so limit it
  if (woden_settings->chunking_size > MAX_CHUNKING_SIZE) {
    printf("Current maximum allowable chunk size is %ld. Anything above that your\nGPU is likely to catch fire. Defaulting to this value.\n", MAX_CHUNKING_SIZE);
    woden_settings->chunking_size = MAX_CHUNKING_SIZE;
  }
  else if (woden_settings->chunking_size < 1 ) {
    woden_settings->chunking_size = MAX_CHUNKING_SIZE;
  }

  woden_settings->array_layout_file = json_object_get_boolean(array_layout_file_path);
  if (woden_settings->array_layout_file) {
    woden_settings->array_layout_file_path = json_object_get_string(array_layout_file_path);
    // printf("Will use east,north,height coords from this file: %s\n", woden_settings->array_layout_file_path);

    FILE *fp_test_array=NULL;

    if ((fp_test_array=fopen(woden_settings->array_layout_file_path,"r"))==NULL) {
      printf("Reading of array_layout file:\n %s\nhas failed", woden_settings->array_layout_file_path);
      return 1;
    }


  } else {
    printf("WODEN needs an east,north,height coordinate list. Specify using \
    'array_layout' in setting file\n");
    return 1;
  }

  struct json_object *band_num;
  struct json_object *band_nums;
  size_t num_bands;

  json_object_object_get_ex(parsed_json, "band_nums", &band_nums);
  num_bands = json_object_array_length(band_nums);

  woden_settings->num_bands = num_bands;
  woden_settings->band_nums = malloc( num_bands * sizeof(int));

  for(int i=0;i<num_bands;i++) {
    band_num = json_object_array_get_idx(band_nums, i);
    woden_settings->band_nums[i] = json_object_get_int(band_num);
  }

  //Everything was fine
  return 0;

}
double * setup_lsts_and_phase_centre(woden_settings_t *woden_settings){
  //Useful number to have
  const int num_visis = woden_settings->num_baselines * woden_settings->num_time_steps * woden_settings->num_freqs;
  woden_settings->num_visis = num_visis;

  //Phase centre details
  double sdec0,cdec0;
  sdec0 = sin(woden_settings->dec0); cdec0=cos(woden_settings->dec0);

  printf("Setting initial LST to %.10fdeg\n",woden_settings->lst_base/DD2R );
  printf("Setting phase centre RA,DEC %.5fdeg %.5fdeg\n",woden_settings->ra0/DD2R, woden_settings->dec0/DD2R);

  //Used for calculating l,m,n for components
  woden_settings->sdec0 = sdec0;
  woden_settings->cdec0 = cdec0;

  //Calculate all lsts for this observation
  double *lsts = malloc(woden_settings->num_time_steps*sizeof(double));

  for ( int time_step = 0; time_step < woden_settings->num_time_steps; time_step++ ) {
    double lst = woden_settings->lst_base + time_step*woden_settings->time_res*SOLAR2SIDEREAL*DS2R;

    //Add half a time_res so we are sampling centre of each time step
    lst += 0.5*woden_settings->time_res*SOLAR2SIDEREAL*DS2R;
    lsts[time_step] = lst;
  }
  return lsts;
}
