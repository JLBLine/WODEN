#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fitsio.h>
#include <erfa.h>
#include <complex.h>

#include "create_sky_model.h"
#include "shapelet_basis.h"
#include "woden.h"
#include "constants.h"
#include "chunk_source.h"
#include "print_help.h"

// #include "FEE_primary_beam_cuda.h"

int main(int argc, char **argv) {

  //If not enough arguments, print help
  if (argc < 2) {
    print_cmdline_help();
    exit(1);
  }

  //If --help is passed, print help
  if (strcmp("--help", argv[1]) == 0) {
    print_cmdline_help();
    exit(1);

  }

  //Create the shapelet basis function array
  float *sbf;
  sbf = NULL;
  sbf = malloc( sbf_N * sbf_L * sizeof(float) );
  sbf = create_sbf(sbf);

  //Read in the settings from the controlling json file
  woden_settings_t *woden_settings;
  woden_settings = read_json_settings(argv[1]);

  if (woden_settings->chunking_size > MAX_CHUNKING_SIZE) {
    printf("Current maximum allowable chunk size is %d.  Defaulting to this value.", MAX_CHUNKING_SIZE);
    woden_settings->chunking_size = MAX_CHUNKING_SIZE;
  }
  else if (woden_settings->chunking_size < 1 ) {
    woden_settings->chunking_size = MAX_CHUNKING_SIZE;
  }

  //Read in information from the metafits file
  int status=0;
  static fitsfile *metaf_file=NULL;
  MetaFfile_t metafits;
  fits_open_file(&metaf_file, woden_settings->metafits_filename, READONLY, &status);
  status = RTS_init_meta_file(metaf_file, &metafits, woden_settings->metafits_filename);

  //Create the array layout in instrument-centric X,Y,Z using positions
  //from the metafits file
  array_layout_t * array_layout;
  array_layout = calc_XYZ_diffs(&metafits, metafits.num_tiles);

  //Propagate some of the metafits data into woden_settings
  woden_settings->lst_base = metafits.lst_base;
  woden_settings->base_low_freq = metafits.base_low_freq;
  woden_settings->num_baselines = array_layout->num_baselines;

  //Set some constants based on the settings
  float ha0, sha0, cha0;
  float wavelength;
  float frequency;
  const int num_visis = woden_settings->num_baselines * woden_settings->num_time_steps * woden_settings->num_freqs;
  float sdec0,cdec0;
  sdec0 = sin(woden_settings->dec0); cdec0=cos(woden_settings->dec0);

  printf("Setting phase centre (rad) to %f %f\n",woden_settings->ra0, woden_settings->dec0);
  printf("Obs pointing centre (rad) is %f %f\n",metafits.ra_point, metafits.dec_point);

  //Used for calculating l,m,n for components
  float angles_array[3] = {sdec0, cdec0, woden_settings->ra0};
  int num_time_steps = woden_settings->num_time_steps;

  //Calculate all lsts for this observation
  float lsts[num_time_steps];

  for ( int time_step = 0; time_step < woden_settings->num_time_steps; time_step++ ) {
    float lst = woden_settings->lst_base + time_step*woden_settings->time_res*SOLAR2SIDEREAL*DS2R;
    //TODO add half a time step is good? Add time decorrelation?
    lst += 0.5*woden_settings->time_res*SOLAR2SIDEREAL*DS2R;
    lsts[time_step] = lst;
  }

  //Read in the source catalogue
  source_catalogue_t *raw_srccat;
  raw_srccat = read_source_catalogue(woden_settings->cat_filename);

  //Crop emission below the horizon, and collapse all SOURCES from raw_srccat
  //into one single SOURCE
  printf("Horizon cropping sky model and calculating az/za for all components for observation\n");
  catsource_t *cropped_src;
  cropped_src = crop_sky_model(raw_srccat, lsts, num_time_steps, woden_settings->sky_crop_type);

  printf("Finished cropping and calculating az/za\n");

  //Setup primary beam settings for observation
  beam_settings_t beam_settings; //= malloc(sizeof(beam_settings_t));
  //Angles used in calculating beam style l,m,ns
  beam_settings.beam_angles_array = malloc(3*sizeof(float));
  beam_settings.beam_angles_array[0] = sinf(metafits.dec_point);
  beam_settings.beam_angles_array[1] = cosf(metafits.dec_point);
  beam_settings.beam_angles_array[2] = woden_settings->lst_base - metafits.ra_point;

  //Number of beam calculations needed for point components
  beam_settings.num_point_beam_values = cropped_src->n_points * woden_settings->num_time_steps * woden_settings->num_freqs;
  beam_settings.num_gausscomp_beam_values = cropped_src->n_gauss * woden_settings->num_time_steps * woden_settings->num_freqs;
  beam_settings.num_shape_beam_values = cropped_src->n_shapes * woden_settings->num_time_steps * woden_settings->num_freqs;

  //TODO pull the functionality of this GAUSS_BEAM loop into another file
  //If using a gaussian primary beam, do gaussian beam things
  if (woden_settings->beamtype == GAUSS_BEAM) {
    beam_settings.beamtype = GAUSS_BEAM;

    printf("Setting up Gaussian primary beam settings\n");
    printf("   setting beam FWHM to %.5fdeg and ref freq to %.3fMHz\n",woden_settings->gauss_beam_FWHM,woden_settings->gauss_beam_ref_freq / 1e+6  );

    //Set constants used in beam calculation
    beam_settings.beam_FWHM_rad = woden_settings->gauss_beam_FWHM * D2R;
    //TODO I cannot for the life of me work out how to cudaMalloc and Memcpy
    //a single float (argh) so put the ref freq in an array (embarrassment)
    float beam_ref_freq_array[1] = {woden_settings->gauss_beam_ref_freq};
    beam_settings.beam_ref_freq_array = malloc(sizeof(float));
    beam_settings.beam_ref_freq_array = beam_ref_freq_array;

    //Store all ha (which change with lst) that the beam needs to be calculated at.
    beam_settings.beam_point_has = malloc(woden_settings->num_time_steps * cropped_src->n_points * sizeof(float));
    beam_settings.beam_point_decs = malloc(woden_settings->num_time_steps * cropped_src->n_points * sizeof(float));

    beam_settings.beam_gausscomp_has = malloc(woden_settings->num_time_steps * cropped_src->n_gauss * sizeof(float));
    beam_settings.beam_gausscomp_decs = malloc(woden_settings->num_time_steps * cropped_src->n_gauss * sizeof(float));

    beam_settings.beam_shape_has = malloc(woden_settings->num_time_steps * cropped_src->n_shapes * sizeof(float));
    beam_settings.beam_shape_decs = malloc(woden_settings->num_time_steps * cropped_src->n_shapes * sizeof(float));

    //Loop over all time and point components and calculate ha
    for ( int time_step = 0; time_step < woden_settings->num_time_steps; time_step++ ) {
      for (int component = 0; component < cropped_src->n_points; component++) {
        int step = cropped_src->n_points*time_step + component;

        beam_settings.beam_point_has[step] = lsts[time_step] - cropped_src->point_ras[component];
        beam_settings.beam_point_decs[step] = cropped_src->point_decs[component];

        // printf("Original mothers step, ha %d %f \n",step,beam_settings.beam_point_has[step] );

      }//point loop

    //Loop over all time and gausscomp components and calculate ha
      for (int component = 0; component < cropped_src->n_gauss; component++) {
        int step = cropped_src->n_gauss*time_step + component;

        beam_settings.beam_gausscomp_has[step] = lsts[time_step] - cropped_src->gauss_ras[component];
        beam_settings.beam_gausscomp_decs[step] = cropped_src->gauss_decs[component];
      }//gausscomp loop

    //Loop over all time and shape components and calculate ha
      for (int component = 0; component < cropped_src->n_shapes; component++) {
        int step = cropped_src->n_shapes*time_step + component;

        beam_settings.beam_shape_has[step] = lsts[time_step] - cropped_src->shape_ras[component];
        beam_settings.beam_shape_decs[step] = cropped_src->shape_decs[component];
      }//shape loop

    }//gaussian beam time loop
  } // End if (woden_settings->gaussian_beam)

  else if (woden_settings->beamtype == FEE_BEAM) {
    beam_settings.beamtype = FEE_BEAM;

    //Get the parallactic angle of the beam pointing for every time step
    //Need to rotate the FEE model which is stored in theta/phi pols by the
    //parallactic angle to obtain XX/YY
    beam_settings.para_cosrot = malloc(woden_settings->num_time_steps*sizeof(float));
    beam_settings.para_sinrot = malloc(woden_settings->num_time_steps*sizeof(float));

    double para_angle;
    for ( int time_step = 0; time_step < woden_settings->num_time_steps; time_step++ ) {

      float FEE_HA = lsts[time_step] - metafits.ra_point;
      para_angle = eraHd2pa((double)FEE_HA, (double)metafits.dec_point, (double)MWA_LAT);

      beam_settings.para_cosrot[time_step] = cosf((float)para_angle);
      beam_settings.para_sinrot[time_step] = sinf((float)para_angle);

      printf("PARA ANGLEEEEE %.10f %.10f %.10f\n",para_angle,cosf((float)para_angle),sinf((float)para_angle) );
    }

  }

  else {
    printf("BEAM TYPE %d\n",(int)woden_settings->beamtype );
  }

  //TODO add in the FEE beam model as an option
  //else if (woden_settings->beamtype == MWA_BEAM) {}

  //MWA correlator data is split into 24 'coarse' bands of 1.28MHz bandwidth,
  //which is typically split into 10, 20, or 40kHz fine channels
  //Loop through each coarse frequency band, run the simulation and dump to
  //a binary file
  for (size_t band = 0; band < woden_settings->num_bands; band++) {
    //Set the lower frequency edge for this coarse band
    int band_num = woden_settings->band_nums[band];
    float base_band_freq = ((band_num - 1)*(metafits.bandwidth/24.0)) + woden_settings->base_low_freq;
    printf("Simulating band %02d with bottom freq %.8e\n",band_num,base_band_freq);

    if (woden_settings->beamtype == FEE_BEAM){

      //Just use one single tile beam for all for now - will need a certain
      //number in the future to include dipole flagging
      int st = 0;
      float base_middle_freq = base_band_freq + metafits.bandwidth/48.0;

      //TODO make this a path defined by run_woden.py
      char* HDFbeampath = "/home/jline/software/useful/MWA_embedded_element_pattern_V02.h5";

      printf("Middle freq is %f\n",base_middle_freq );

      // copy_primary_beam_t *FEE_beam;
      // FEE_beam = malloc(sizeof(copy_primary_beam_t));
      beam_settings.FEE_beam = malloc(sizeof(copy_primary_beam_t));
      printf("Setting up the FEE beam...");
      RTS_HDFBeamInit(HDFbeampath, base_middle_freq, beam_settings.FEE_beam, (float *)metafits.FEE_delays[st], st);
      printf(" done.\n");

      printf("Getting FEE beam normalisation...");
      get_HDFBeam_normalisation(beam_settings.FEE_beam);
      printf(" done.\n");
      // beam_settings.FEE_beam->norm_fac[0] = 0.25714464415545296 + 0*I;
      // beam_settings.FEE_beam->norm_fac[1] = 0.25714464415545296 + 0*I;
      // beam_settings.FEE_beam->norm_fac[2] = 0.25729902904652246 + 0*I;
      // beam_settings.FEE_beam->norm_fac[3] = 0.25729902904652246 + 0*I;

      // beam_settings.FEE_beam->norm_fac[0] = 0.46907393930481284 + 0*I;
      // beam_settings.FEE_beam->norm_fac[1] = 0.46907393930481284 + 0*I;
      // beam_settings.FEE_beam->norm_fac[2] = 0.4694292262227538 + 0*I;
      // beam_settings.FEE_beam->norm_fac[3] = 0.4694292262227538 + 0*I;


      // printf("NMAX in WODEN.C IS %d\n",beam_settings.FEE_beam->nmax );
      printf("Copying the FEE beam across to the GPU...");
      copy_FEE_primary_beam_to_GPU(beam_settings, woden_settings->num_time_steps);
      printf(" done.\n");
      // printf("Have sent the FEE beam to the GPU IS DIFFERENT\n");

    }

    visibility_set_t *visibility_set = malloc(sizeof(visibility_set_t));
    visibility_set->sum_visi_real = malloc( num_visis * sizeof(float) );
    visibility_set->sum_visi_imag = malloc( num_visis * sizeof(float) );
    visibility_set->us_metres = malloc( num_visis * sizeof(float) );
    visibility_set->vs_metres = malloc( num_visis * sizeof(float) );
    visibility_set->ws_metres = malloc( num_visis * sizeof(float) );
    visibility_set->sha0s = malloc( num_visis * sizeof(float) );
    visibility_set->cha0s = malloc( num_visis * sizeof(float) );
    visibility_set->lsts = malloc( num_visis * sizeof(float) );
    visibility_set->wavelengths = malloc( num_visis * sizeof(float) );
    visibility_set->channel_frequencies = malloc( (int)woden_settings->num_freqs * sizeof(float) );

    visibility_set->sum_visi_XX_real = malloc( num_visis * sizeof(float) );
    visibility_set->sum_visi_XX_imag = malloc( num_visis * sizeof(float) );
    visibility_set->sum_visi_XY_real = malloc( num_visis * sizeof(float) );
    visibility_set->sum_visi_XY_imag = malloc( num_visis * sizeof(float) );
    visibility_set->sum_visi_YX_real = malloc( num_visis * sizeof(float) );
    visibility_set->sum_visi_YX_imag = malloc( num_visis * sizeof(float) );
    visibility_set->sum_visi_YY_real = malloc( num_visis * sizeof(float) );
    visibility_set->sum_visi_YY_imag = malloc( num_visis * sizeof(float) );

    visibility_set->sum_visi_real = malloc( num_visis * sizeof(float) );
    visibility_set->sum_visi_imag = malloc( num_visis * sizeof(float) );

    //Useful for testing beam things
    // visibility_set->beam_has = malloc(woden_settings->num_time_steps * cropped_src->n_shapes * sizeof(float) );
    // visibility_set->beam_decs = malloc(woden_settings->num_time_steps * cropped_src->n_shapes * sizeof(float) );
    // visibility_set->beam_ls = malloc(woden_settings->num_time_steps * cropped_src->n_shapes * sizeof(float) );
    // visibility_set->beam_ms = malloc(woden_settings->num_time_steps * cropped_src->n_shapes * sizeof(float) );
    // visibility_set->beam_reals = malloc(woden_settings->num_time_steps * cropped_src->n_shapes * woden_settings->num_freqs * sizeof(float) );
    // visibility_set->beam_imags = malloc(woden_settings->num_time_steps * cropped_src->n_shapes * woden_settings->num_freqs * sizeof(float) );

    //Fill in the fine channel frequencies
    for (int freq_step = 0; freq_step < woden_settings->num_freqs; freq_step++) {
      frequency = base_band_freq + (woden_settings->frequency_resolution*freq_step);
      visibility_set->channel_frequencies[freq_step] = frequency;
    }

    //Fill in visibility settings in order of baseline,freq,time
    //Order matches that of a uvfits file (I live in the past)

    for ( int time_step = 0; time_step < woden_settings->num_time_steps; time_step++ ) {
      ha0 = lsts[time_step] - woden_settings->ra0;
      sha0 = sin(ha0); cha0=cos(ha0);

      for (int freq_step = 0; freq_step < woden_settings->num_freqs; freq_step++) {
        frequency = base_band_freq + (woden_settings->frequency_resolution*freq_step);
        wavelength = VELC / frequency;
        int step = woden_settings->num_baselines*(time_step*woden_settings->num_freqs + freq_step);

        for (int baseline = 0; baseline < woden_settings->num_baselines; baseline++) {
          visibility_set->cha0s[step + baseline] = cha0;
          visibility_set->sha0s[step + baseline] = sha0;
          visibility_set->lsts[step + baseline] = lsts[time_step];
          visibility_set->wavelengths[step + baseline] = wavelength;
        }//baseline loop
      }//freq loop
    }//time loop

    //Calculating a single shapelet coeff is equivalent to a point/gauss so treat as a
    //component here
    int num_components = cropped_src->n_points + cropped_src->n_gauss + cropped_src->n_shape_coeffs;

    //TODO should we chunk outside the band for-loop so that we can reuse the chunks for each band (should be the same)
    if (num_components > woden_settings->chunking_size) {
      printf("Chunking sky model\n");
      int num_chunks = num_components / woden_settings->chunking_size;

      if (num_components % woden_settings->chunking_size != 0) {
        num_chunks = num_chunks + 1;
      }

      printf("Number of chunks required is %d\n",  num_chunks);

      //setup a temporary visibility set that calculate_visibilities will populate
      visibility_set_t *temp_visibility_set = malloc(sizeof(visibility_set_t));
      // temp_visibility_set->sum_visi_real = malloc( num_visis * sizeof(float) );
      // temp_visibility_set->sum_visi_imag = malloc( num_visis * sizeof(float) );

      temp_visibility_set->us_metres = malloc( num_visis * sizeof(float) );
      temp_visibility_set->vs_metres = malloc( num_visis * sizeof(float) );
      temp_visibility_set->ws_metres = malloc( num_visis * sizeof(float) );

      temp_visibility_set->sha0s = visibility_set->sha0s;
      temp_visibility_set->cha0s = visibility_set->cha0s;
      temp_visibility_set->lsts = visibility_set->lsts;
      temp_visibility_set->wavelengths = visibility_set->wavelengths;
      temp_visibility_set->channel_frequencies = visibility_set->channel_frequencies;

      temp_visibility_set->sum_visi_XX_real = malloc( num_visis * sizeof(float) );
      temp_visibility_set->sum_visi_XX_imag = malloc( num_visis * sizeof(float) );
      temp_visibility_set->sum_visi_XY_real = malloc( num_visis * sizeof(float) );
      temp_visibility_set->sum_visi_XY_imag = malloc( num_visis * sizeof(float) );
      temp_visibility_set->sum_visi_YX_real = malloc( num_visis * sizeof(float) );
      temp_visibility_set->sum_visi_YX_imag = malloc( num_visis * sizeof(float) );
      temp_visibility_set->sum_visi_YY_real = malloc( num_visis * sizeof(float) );
      temp_visibility_set->sum_visi_YY_imag = malloc( num_visis * sizeof(float) );

      catsource_t *temp_cropped_src = malloc(sizeof(catsource_t));
      int point_iter = 0;
      int gauss_iter = 0;
      int shape_iter = 0;

      //For each chunk, calculate the visibilities for those components
      for (int chunk = 0; chunk < num_chunks; chunk++) {
        printf("Processing chunk %d\n", chunk);

        //ensure temp visi's are 0.0
        for (size_t visi = 0; visi < num_visis; visi++) {
          // temp_visibility_set->sum_visi_real[visi] = 0.0;
          // temp_visibility_set->sum_visi_imag[visi] = 0.0;

          temp_visibility_set->sum_visi_XX_real[visi] = 0.0;
          temp_visibility_set->sum_visi_XX_imag[visi] = 0.0;
          temp_visibility_set->sum_visi_XY_real[visi] = 0.0;
          temp_visibility_set->sum_visi_XY_imag[visi] = 0.0;
          temp_visibility_set->sum_visi_YX_real[visi] = 0.0;
          temp_visibility_set->sum_visi_YX_imag[visi] = 0.0;
          temp_visibility_set->sum_visi_YY_real[visi] = 0.0;
          temp_visibility_set->sum_visi_YY_imag[visi] = 0.0;

        }

        fill_chunk_src(temp_cropped_src, cropped_src, num_chunks, chunk,
                       woden_settings->chunking_size, woden_settings->num_time_steps,
                       &point_iter, &gauss_iter, &shape_iter);
        // printf("WOT %d %d %d\n", point_iter, gauss_iter, shape_iter);

        printf("\tNumber of components in chunk are: P %d G %d S_coeffs %d\n",temp_cropped_src->n_points,temp_cropped_src->n_gauss,temp_cropped_src->n_shape_coeffs );

        // for ( int time_step = 0; time_step < woden_settings->num_time_steps; time_step++ ) {
        //   for (int component = 0; component < temp_cropped_src->n_points; component++) {
        //     int chunk_step = temp_cropped_src->n_points*time_step + component;
        //     int crop_step = cropped_src->n_points*time_step + point_iter + component;
        //
        //     printf("\tHA val %f %f\n",beam_settings_chunk.beam_point_has[chunk_step],beam_settings.beam_point_has[crop_step] );
        //     printf("\tDec val %f %f\n",beam_settings_chunk.beam_point_decs[chunk_step],beam_settings.beam_point_decs[crop_step] );
        //   }
        // }

        if (woden_settings->beamtype == GAUSS_BEAM){

          beam_settings_t beam_settings_chunk;
          beam_settings_chunk = make_beam_settings_chunk(beam_settings, temp_cropped_src,
                                cropped_src, woden_settings, point_iter, gauss_iter, shape_iter);

          calculate_visibilities(array_layout->X_diff_metres, array_layout->Y_diff_metres, array_layout->Z_diff_metres,
                      *temp_cropped_src, angles_array, beam_settings_chunk,
                      woden_settings->num_baselines, woden_settings->num_time_steps,
                      num_visis, woden_settings->num_freqs, temp_visibility_set,
                      sbf);

          free(beam_settings_chunk.beam_point_has );
          free(beam_settings_chunk.beam_point_decs );
          free(beam_settings_chunk.beam_gausscomp_has );
          free(beam_settings_chunk.beam_gausscomp_decs );
          free(beam_settings_chunk.beam_shape_has );
          free(beam_settings_chunk.beam_shape_decs );
        }

        else {
          // printf("This be where we call num gauss %d\n", temp_cropped_src->n_gauss);

          calculate_visibilities(array_layout->X_diff_metres, array_layout->Y_diff_metres, array_layout->Z_diff_metres,
                      *temp_cropped_src, angles_array, beam_settings,
                      woden_settings->num_baselines, woden_settings->num_time_steps,
                      num_visis, woden_settings->num_freqs, temp_visibility_set,
                      sbf);
        }


        // printf("Adding temporary visibility set.\n");
	      //add to visiblity_set
        for (int visi = 0; visi < num_visis; visi++) {
          //if the first chunk then initialise our values, and copy across
          //the u,v,w coords
          if (chunk == 0) {
            // visibility_set->sum_visi_real[visi] = 0;
            // visibility_set->sum_visi_imag[visi] = 0;

            visibility_set->sum_visi_XX_real[visi] = 0;
            visibility_set->sum_visi_XX_imag[visi] = 0;
            visibility_set->sum_visi_XY_real[visi] = 0;
            visibility_set->sum_visi_XY_imag[visi] = 0;
            visibility_set->sum_visi_YX_real[visi] = 0;
            visibility_set->sum_visi_YX_imag[visi] = 0;
            visibility_set->sum_visi_YY_real[visi] = 0;
            visibility_set->sum_visi_YY_imag[visi] = 0;

            visibility_set->us_metres[visi] = temp_visibility_set->us_metres[visi];
            visibility_set->vs_metres[visi] = temp_visibility_set->vs_metres[visi];
            visibility_set->ws_metres[visi] = temp_visibility_set->ws_metres[visi];
          }

          //add each chunk of components to visibility set
          visibility_set->sum_visi_XX_real[visi] += temp_visibility_set->sum_visi_XX_real[visi];
          visibility_set->sum_visi_XX_imag[visi] += temp_visibility_set->sum_visi_XX_imag[visi];
          visibility_set->sum_visi_XY_real[visi] += temp_visibility_set->sum_visi_XY_real[visi];
          visibility_set->sum_visi_XY_imag[visi] += temp_visibility_set->sum_visi_XY_imag[visi];
          visibility_set->sum_visi_YX_real[visi] += temp_visibility_set->sum_visi_YX_real[visi];
          visibility_set->sum_visi_YX_imag[visi] += temp_visibility_set->sum_visi_YX_imag[visi];
          visibility_set->sum_visi_YY_real[visi] += temp_visibility_set->sum_visi_YY_real[visi];
          visibility_set->sum_visi_YY_imag[visi] += temp_visibility_set->sum_visi_YY_imag[visi];

          // if (visi == 0) {
          //   printf("Like this slappa da bass %.5f %.5f\n", temp_visibility_set->sum_visi_real[visi], temp_visibility_set->sum_visi_imag[visi]);
          // }



        }//visi loop
      }//chunk loop

      free(temp_cropped_src);

      // free( temp_visibility_set->sum_visi_real );
      // free( temp_visibility_set->sum_visi_imag );
      free( temp_visibility_set->us_metres );
      free( temp_visibility_set->vs_metres );
      free( temp_visibility_set->ws_metres );

      free(temp_visibility_set->sum_visi_XX_real);
      free(temp_visibility_set->sum_visi_XX_imag);
      free(temp_visibility_set->sum_visi_XY_real);
      free(temp_visibility_set->sum_visi_XY_imag);
      free(temp_visibility_set->sum_visi_YX_real);
      free(temp_visibility_set->sum_visi_YX_imag);
      free(temp_visibility_set->sum_visi_YY_real);
      free(temp_visibility_set->sum_visi_YY_imag);

      free( temp_visibility_set );

    }
    //If not chunking the components, just simulate all in one go
    else {
      //Throw all of the settings at the GPU and crank the handle on the simulation
      calculate_visibilities(array_layout->X_diff_metres, array_layout->Y_diff_metres, array_layout->Z_diff_metres,
                      *cropped_src, angles_array, beam_settings,
                      woden_settings->num_baselines, woden_settings->num_time_steps,
                      num_visis, woden_settings->num_freqs, visibility_set,
                      sbf);
    }

    if (woden_settings->beamtype == FEE_BEAM) {
      free_FEE_primary_beam_from_GPU(beam_settings.FEE_beam);
    }

    // free( beam_settings );

    // if (woden_settings->beamtype == FEE_BEAM){
    //   free( beam_settings.FEE_beam );
    // }

    //Dumps u,v,w (metres), Re(vis), Im(vis) to a binary file
    FILE *output_visi;
    char buf[0x100];
    snprintf(buf, sizeof(buf), "output_visi_band%02d.dat", band_num);

    output_visi = fopen(buf,"wb");

    if(output_visi == NULL)
    {
        printf("Could not open output_visi_band%02d.dat - exiting", band_num);
        exit(1);
    }

    fwrite(visibility_set->us_metres, num_visis*sizeof(float), 1, output_visi);
    fwrite(visibility_set->vs_metres, num_visis*sizeof(float), 1, output_visi);
    fwrite(visibility_set->ws_metres, num_visis*sizeof(float), 1, output_visi);
    // fwrite(visibility_set->sum_visi_real, num_visis*sizeof(float), 1, output_visi);
    // fwrite(visibility_set->sum_visi_imag, num_visis*sizeof(float), 1, output_visi);

    fwrite(visibility_set->sum_visi_XX_real, num_visis*sizeof(float), 1, output_visi);
    fwrite(visibility_set->sum_visi_XX_imag, num_visis*sizeof(float), 1, output_visi);
    fwrite(visibility_set->sum_visi_XY_real, num_visis*sizeof(float), 1, output_visi);
    fwrite(visibility_set->sum_visi_XY_imag, num_visis*sizeof(float), 1, output_visi);
    fwrite(visibility_set->sum_visi_YX_real, num_visis*sizeof(float), 1, output_visi);
    fwrite(visibility_set->sum_visi_YX_imag, num_visis*sizeof(float), 1, output_visi);
    fwrite(visibility_set->sum_visi_YY_real, num_visis*sizeof(float), 1, output_visi);
    fwrite(visibility_set->sum_visi_YY_imag, num_visis*sizeof(float), 1, output_visi);

    fflush(output_visi);
    fclose(output_visi);


    // Dumps u,v,w (metres), Re(vis), Im(vis) directly to text file - useful for
    // bug hunting with small outputs
    // FILE *output_visi_text;
    // char buff[0x100];
    // snprintf(buff, sizeof(buff), "output_visi_band%02d.txt", band_num);
    // output_visi_text = fopen(buff,"w");
    // for ( int time_step = 0; time_step < woden_settings->num_time_steps; time_step++ ) {
    //   for ( int freq_step = 0; freq_step < woden_settings->num_freqs; freq_step++ ) {
    //     for (int baseline = 0; baseline < woden_settings->num_baselines; baseline++) {
    //       int step = woden_settings->num_baselines*(time_step*woden_settings->num_freqs + freq_step);
    //       fprintf(output_visi_text,"%f %f %f %f %f\n",visibility_set->us_metres[step + baseline],
    //               visibility_set->vs_metres[step + baseline],visibility_set->ws_metres[step + baseline],
    //               visibility_set->sum_visi_XX_real[step + baseline],visibility_set->sum_visi_XX_imag[step + baseline]);
    //     }
    //   }
    // }
    // fflush(output_visi_text);
    // fclose(output_visi_text);
    //
    // // Dumps u,v,w (metres), Re(vis), Im(vis) directly to text file - useful for
    // // bug hunting with small outputs
    // FILE *output_visi_text;
    // char buff[0x100];
    // snprintf(buff, sizeof(buff), "output_visi_band%02d.txt", band_num);
    // output_visi_text = fopen(buff,"w");
    // for ( int time_step = 0; time_step < woden_settings->num_time_steps; time_step++ ) {
    //   for ( int freq_step = 0; freq_step < woden_settings->num_freqs; freq_step++ ) {
    //     for (int baseline = 0; baseline < woden_settings->num_baselines; baseline++) {
    //       int step = woden_settings->num_baselines*(time_step*woden_settings->num_freqs + freq_step);
    //       fprintf(output_visi_text,"%f %f %f %f %f\n",visibility_set->us_metres[step + baseline],
    //               visibility_set->vs_metres[step + baseline],visibility_set->ws_metres[step + baseline],
    //               visibility_set->sum_visi_real[step + baseline],visibility_set->sum_visi_imag[step + baseline]);
    //     }
    //   }
    // }
    // fflush(output_visi_text);
    // fclose(output_visi_text);
    // //
    // // // Beam testing text file
    // FILE *output_beamcoords;
    // char bufff[0x100];
    // snprintf(bufff, sizeof(bufff), "output_beam_coords%02d.txt", band_num);
    // output_beamcoords = fopen(bufff,"w");
    // // for ( int time_step = 0; time_step < woden_settings->num_time_steps; time_step++ ) {
    // //     for (int component = 0; component < cropped_src->n_points; component++) {
    // //       int step = cropped_src->n_points*time_step;
    // //       fprintf(output_beamcoords,"%f %f %f %f %f %f\n",
    // //               beam_settings.beam_point_has[step + component], beam_settings.beam_point_decs[step + component],
    // //               visibility_set->beam_has[step + component], visibility_set->beam_decs[step + component],
    // //               visibility_set->beam_ls[step + component], visibility_set->beam_ms[step + component]);
    // //     }
    // // }
    // for ( int time_step = 0; time_step < woden_settings->num_time_steps; time_step++ ) {
    //     for (int component = 0; component < cropped_src->n_shapes; component++) {
    //       int step = cropped_src->n_shapes*time_step;
    //       fprintf(output_beamcoords,"%f %f %f %f %f %f\n",
    //               beam_settings.beam_shape_has[step + component], beam_settings.beam_shape_decs[step + component],
    //               visibility_set->beam_has[step + component], visibility_set->beam_decs[step + component],
    //               visibility_set->beam_ls[step + component], visibility_set->beam_ms[step + component]);
    //     }
    // }
    //
    // fflush(output_beamcoords);
    // fclose(output_beamcoords);
    // // //
    // // //
    // // //// Beam testing text file
    // FILE *output_beamvals;
    // char buffff[0x100];
    // snprintf(buffff, sizeof(buffff), "output_beam_values%02d.txt", band_num);
    // output_beamvals = fopen(buffff,"w");
    // for ( int time_step = 0; time_step < woden_settings->num_time_steps; time_step++ ) {
    //   for ( int freq_step = 0; freq_step < woden_settings->num_freqs; freq_step++ ) {
    //     // for (int component = 0; component < cropped_src->n_points; component++) {
    //     //   int step = cropped_src->n_points*(time_step*woden_settings->num_freqs + freq_step);
    //     //   fprintf(output_beamvals,"%f %f\n",
    //     //           visibility_set->beam_reals[step + component],visibility_set->beam_imags[step + component]);
    //     // }
    //     // for (int component = 0; component < cropped_src->n_gauss; component++) {
    //     //   int step = cropped_src->n_gauss*(time_step*woden_settings->num_freqs + freq_step);
    //     //   fprintf(output_beamvals,"%f %f\n",
    //     //           visibility_set->beam_reals[step + component],visibility_set->beam_imags[step + component]);
    //     // }
    //     for (int component = 0; component < cropped_src->n_shapes; component++) {
    //       int step = cropped_src->n_shapes*(time_step*woden_settings->num_freqs + freq_step);
    //       fprintf(output_beamvals,"%f %f\n",
    //               visibility_set->beam_reals[step + component],visibility_set->beam_imags[step + component]);
    //     }
    //   }
    // }
    // fflush(output_beamvals);
    // fclose(output_beamvals);
    //
    // free( visibility_set->beam_imags );
    // free( visibility_set->beam_reals );
    // free( visibility_set->beam_ls );
    // free( visibility_set->beam_ms );
    // free( visibility_set->beam_has );
    // free( visibility_set->beam_decs );

    //Free up that memory
    // printf("Made it here 1?\n");
    // free( visibility_set->sum_visi_real );
    // printf("Made it here 2?\n");
    // free( visibility_set->sum_visi_imag );
    // printf("Made it here 3?\n");
    free( visibility_set->us_metres );
    // printf("Made it here 4?\n");
    free( visibility_set->vs_metres );
    // printf("Made it here 5?\n");
    free( visibility_set->ws_metres );
    // printf("Made it here 6?\n");
    free( visibility_set->sha0s );
    // printf("Made it here 7?\n");
    free( visibility_set->cha0s );
    // printf("Made it here 8?\n");
    free( visibility_set->lsts );
    // printf("Made it here 9?\n");
    free( visibility_set->wavelengths );
    // printf("Made it here 10?\n");

    free(visibility_set->sum_visi_XX_real);
    free(visibility_set->sum_visi_XX_imag);
    free(visibility_set->sum_visi_XY_real);
    free(visibility_set->sum_visi_XY_imag);
    free(visibility_set->sum_visi_YX_real);
    free(visibility_set->sum_visi_YX_imag);
    free(visibility_set->sum_visi_YY_real);
    free(visibility_set->sum_visi_YY_imag);




    free( visibility_set );
    // printf("Made it here 11?\n");

  }//band loop
  printf("THE END?\n");
}//main
