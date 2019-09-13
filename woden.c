#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "/usr/local/mwa-RTS_dev/include/slamac.h"
#include "read_and_write.h"
// #include "point_source_lib.h"
#include "shapelet_basis.h"
#include "woden.h"

int main(void)
{
  float *sbf2;
  sbf2 = NULL;
  sbf2 = malloc( sbf_N2 * sbf_L2 * sizeof(float) );
  sbf2 = create_sbf2(sbf2);

  // float lst_base = 71.64335359166665*D2R; //phase1
  float lst_base = 75.3035563889981*D2R; //phase2

  // float ra0 = 71.64335359166665*D2R;
  // float dec0 = -26.7*D2R;

  float ra0 = 50.67*D2R;
  float dec0 = -37.2*D2R;

  int num_baselines = 8128;
  int num_freqs = 128;
  // int num_freqs = 1;
  float frequency_resolution = 10.0e+3;
  // float base_frequency = 167.645e+6;
  // float base_frequency = 169.595e+6;
  //
  float base_frequency = 1.8496e+08;

  float sdec0,cdec0;
  sdec0 = sin(dec0); cdec0=cos(dec0);

  // int num_time_steps = 56;
  int num_time_steps = 20;
  float time_res = 0.5;
  float ha0,sha0,cha0;


  // char cat_filename[] = "srclist_wsclean_ForA_phase1_cropped_percent100.txt";
  char cat_filename[] = "srclist_wsclean_ForA_phase1+2_cropped_percent100.txt";
  // char cat_filename[] = "srclist_singlepoint.txt";
  // char cat_filename[] = "srclist_shape_ForA_phase1+2_percent100.txt";


  source_catalogue_t *srccat;

  srccat = read_source_catalogue(cat_filename);

  // for (int i = 0; i < 1; i++) {
  //   printf("===================================\n");
  //   printf("shape_ra %d %f\n",i,srccat->catsources[0].S2_ras[i]);
  //   printf("shape_dec %d %f\n",i,srccat->catsources[0].S2_decs[i]);
  //   // printf("gauss ra %d %f\n",i,srccat->catsource[0].gauss_ras[i]);
  //   printf("shape flux %d %f\n",i,srccat->catsources[0].S2_fluxes[i]);
  //
  //   for (int j = 0; j < 2; j++) {
  //     printf("------------------------------\n");
  //     printf("shape n1s %d %f\n",j,srccat->catsources[0].S2_n1s[j]);
  //     printf("shape n2s %d %f\n",j,srccat->catsources[0].S2_n2s[j]);
  //     printf("shape coeffs %d %f\n",j,srccat->catsources[0].S2_coeffs[j]);
  //     printf("shape params %d %f\n",j,srccat->catsources[0].S2_param_indexes[j]);
  //   }
  // }

  //Gets these arrays from telescope_XYZ.c
  extern float X_diff_metres[ ];
  extern float Y_diff_metres[ ];
  extern float Z_diff_metres[ ];

  float wavelength;
  float frequency;

  // assign_pointsource_on_GPU(srccat->catsource[0]);

  // d_XYZ_t d_XYZ_metres;
  // d_XYZ_metres.d_X_diff = NULL;
  // d_XYZ_metres.d_Y_diff = NULL;
  // d_XYZ_metres.d_Z_diff = NULL;
  // copy_XYZ_to_GPU(&(d_XYZ_metres.d_X_diff),&(d_XYZ_metres.d_Y_diff),&(d_XYZ_metres.d_Z_diff),
  //                   X_diff_metres, Y_diff_metres, Z_diff_metres,
  //                   num_baselines);

  for (int freq_step = 0; freq_step < num_freqs; freq_step++) {
    visibility_set_t visibility_set;

    // frequency = 167.035e+6 + (frequency_resolution*freq_step);
    frequency = base_frequency + (frequency_resolution*freq_step);
    wavelength = VELC / frequency;

    visibility_set.sum_visi_real = malloc( num_baselines * num_time_steps * sizeof(float) );
    visibility_set.sum_visi_imag = malloc( num_baselines * num_time_steps * sizeof(float) );
    visibility_set.us_metres = malloc( num_baselines * num_time_steps * sizeof(float) );
    visibility_set.vs_metres = malloc( num_baselines * num_time_steps * sizeof(float) );
    visibility_set.ws_metres = malloc( num_baselines * num_time_steps * sizeof(float) );
    visibility_set.sha0s = malloc( num_baselines * num_time_steps * sizeof(float) );
    visibility_set.cha0s = malloc( num_baselines * num_time_steps * sizeof(float) );
    visibility_set.lsts = malloc( num_baselines * num_time_steps * sizeof(float) );

    float wavelengths[num_baselines * num_time_steps];



    for ( int time_step = 0; time_step < num_time_steps; time_step++ ) {
      float lst = lst_base + (time_step*time_res*SOLAR2SIDEREAL*D2R)*(15.0/3600.0);
      ha0 = lst - ra0;
      sha0 = sin(ha0); cha0=cos(ha0);

      for (int baseline = 0; baseline < num_baselines; baseline++) {
        visibility_set.cha0s[time_step*num_baselines + baseline] = cha0;
        visibility_set.sha0s[time_step*num_baselines + baseline] = sha0;
        visibility_set.lsts[time_step*num_baselines + baseline] = lst;

        wavelengths[time_step*num_baselines + baseline] = wavelength;
      }

    }

    // for ( int time_step = 0; time_step < num_time_steps; time_step++ ) {
    //   for (int baseline = 0; baseline < num_baselines; baseline++) {
    //     printf("%f %f\n",visibility_set.cha0s[time_step*num_baselines + baseline],visibility_set.sha0s[time_step*num_baselines + baseline]);
    //   }
    //
    // }

      // printf("%d %f %f %f %f %f %f\n",num_components,lst,ha0,ra0,dec0,comp_ras[0],comp_decs[0] );

      // visibility_set.sum_visi_real = {0};
      // visibility_set.sum_visi_imag = {0};
      // visibility_set.us_metres = {0};
      // visibility_set.vs_metres = {0};
      // visibility_set.ws_metres = {0};


      // printf("%f\n", wavelength);
      // printf("%f %f %f %f %f\n",sdec0, cdec0, sha0, cha0, ra0 );
      // d_XYZ_metres.d_X_diff, d_XYZ_metres.d_Y_diff, d_XYZ_metres.d_Z_diff,

    float angles_array[4] = {sdec0, cdec0, ra0, wavelength};

    Atomic_time_step(X_diff_metres, Y_diff_metres, Z_diff_metres,
                    srccat->catsources[0],angles_array,
                    num_baselines, num_time_steps, visibility_set,
                    sbf2);
    //
    FILE *output_visi;


    output_visi = fopen("output_visi.dat","ab");

    if(output_visi == NULL)
    {
        printf("Could not open output_visi.txt - exiting");
        exit(1);
    }



    fwrite(visibility_set.us_metres, num_baselines*num_time_steps*sizeof(float), 1, output_visi);
    fwrite(visibility_set.vs_metres, num_baselines*num_time_steps*sizeof(float), 1, output_visi);
    fwrite(visibility_set.vs_metres, num_baselines*num_time_steps*sizeof(float), 1, output_visi);
    fwrite(visibility_set.sum_visi_real, num_baselines*num_time_steps*sizeof(float), 1, output_visi);
    fwrite(visibility_set.sum_visi_imag, num_baselines*num_time_steps*sizeof(float), 1, output_visi);

    fflush(output_visi);

    fclose(output_visi);

    // output_visi = fopen("output_visi.txt","a");
    //
    // int ii;
    // for ( int time_step = 0; time_step < num_time_steps; time_step++ ) {
    //   for (int baseline = 0; baseline < num_baselines; baseline++) {
    //     ii = time_step*num_baselines + baseline;
    //     fprintf(output_visi,"%f %f %f %f %f\n",visibility_set.us_metres[ii],
    //             visibility_set.vs_metres[ii],visibility_set.ws_metres[ii],visibility_set.sum_visi_real[ii],visibility_set.sum_visi_imag[ii]);
    //   }
    //
    // }
    // fclose(output_visi);

  }//freq loop

}
