#include <math.h>
#include <stdlib.h>
#include <complex.h>

#include "constants.h"
#include "woden_struct_defs.h"
#include "FEE_primary_beam.h"

//External CUDA code we're linking in
extern void test_RTS_CUDA_FEE_beam(int num_components,
           user_precision_t *azs, user_precision_t *zas, double latitude,
           RTS_MWA_FEE_beam_t *FEE_beam_zenith,
           RTS_MWA_FEE_beam_t *FEE_beam,
           int rotation, int scaling,
           user_precision_complex_t *FEE_beam_gains);

#define UNITY_INCLUDE_FLOAT

//Different delays settings, which control the pointing of the MWA beam
user_precision_t zenith_delays[16] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                           0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

void run_RTS_MWA_FEE_beam(user_precision_t *azs,  user_precision_t *zas, int num_azza,
                          double freq, user_precision_t *delays,
                          char* mwa_fee_hdf5,
                          char *outname) {

  //Call the C code to interrogate the hdf5 file and set beam things up
  RTS_MWA_FEE_beam_t *FEE_beam = malloc(sizeof(RTS_MWA_FEE_beam_t));

  //Get a zenith pointing beam for normalisation purposes
  RTS_MWA_FEE_beam_t *FEE_beam_zenith = malloc(sizeof(RTS_MWA_FEE_beam_t));

  printf("\n\tSetting up the zenith FEE beam...");
  RTS_MWAFEEInit(mwa_fee_hdf5, freq, FEE_beam_zenith, zenith_delays);
  printf(" done.\n");

  printf("\tSetting up the FEE beam...");
  RTS_MWAFEEInit(mwa_fee_hdf5, freq, FEE_beam, delays);
  printf(" done.\n");

  //Set up a bunch of az/za

  //Rotate by parallactic angles
  int rotation = 1;
  //Scale to zenith
  int scaling = 1;

  int num_pols = 4;
  user_precision_complex_t *FEE_beam_gains = malloc(num_pols*num_azza*sizeof(user_precision_complex_t));

  //Run the CUDA code
  test_RTS_CUDA_FEE_beam(num_azza,
             azs, zas, MWA_LAT,
             FEE_beam_zenith,
             FEE_beam,
             rotation, scaling,
             FEE_beam_gains);

    FILE *beam_values_out;
    // char buff[0x100];
    // snprintf(buff, sizeof(buff), outname);
    // beam_values_out = fopen(buff,"w");
    beam_values_out = fopen(outname,"w");

    for (size_t comp = 0; comp < num_azza; comp++) {
      fprintf(beam_values_out, "%.5f %.5f %.7f %.7f %.7f %.7f %.7f %.7f %.7f %.7f\n",
      azs[comp], zas[comp],
      creal(FEE_beam_gains[num_pols*comp+0]), cimag(FEE_beam_gains[num_pols*comp+0]),
      creal(FEE_beam_gains[num_pols*comp+1]), cimag(FEE_beam_gains[num_pols*comp+1]),
      creal(FEE_beam_gains[num_pols*comp+2]), cimag(FEE_beam_gains[num_pols*comp+2]),
      creal(FEE_beam_gains[num_pols*comp+3]), cimag(FEE_beam_gains[num_pols*comp+3]));
  }

  free(FEE_beam_gains);

}


void RTS_MWA_beam_by_freq(char* mwa_fee_hdf5){

  user_precision_t *single_az = malloc(sizeof(user_precision_t));
  user_precision_t *single_za = malloc(sizeof(user_precision_t));

  single_az[0] = 55.0 * (M_PI / 180.0);
  single_za[0] = 60.0 * (M_PI / 180.0);

  int num_azza = 1;

  double freq;

  int num_pols = 4;
  user_precision_complex_t *FEE_beam_gains = malloc(num_pols*num_azza*sizeof(user_precision_complex_t));

  // /Call the C code to interrogate the hdf5 file and set beam things up
  RTS_MWA_FEE_beam_t *FEE_beam = malloc(sizeof(RTS_MWA_FEE_beam_t));
  //Get a zenith pointing beam for normalisation purposes
  RTS_MWA_FEE_beam_t *FEE_beam_zenith = malloc(sizeof(RTS_MWA_FEE_beam_t));


  int num_freqs = 48;

  user_precision_t *out_gx_re = malloc(num_freqs*sizeof(user_precision_t));
  user_precision_t *out_gx_im = malloc(num_freqs*sizeof(user_precision_t));
  user_precision_t *out_Dx_re = malloc(num_freqs*sizeof(user_precision_t));
  user_precision_t *out_Dx_im = malloc(num_freqs*sizeof(user_precision_t));
  user_precision_t *out_Dy_re = malloc(num_freqs*sizeof(user_precision_t));
  user_precision_t *out_Dy_im = malloc(num_freqs*sizeof(user_precision_t));
  user_precision_t *out_gy_re = malloc(num_freqs*sizeof(user_precision_t));
  user_precision_t *out_gy_im = malloc(num_freqs*sizeof(user_precision_t));
  user_precision_t *freqs_used = malloc(num_freqs*sizeof(user_precision_t));

  for (int i = 0; i < num_freqs; i++) {
    freq = 168e+6 + i*240e+3;

    freqs_used[i] = freq;

    printf("\n\tSetting up the zenith FEE beam...");
    RTS_MWAFEEInit(mwa_fee_hdf5, freq, FEE_beam_zenith, zenith_delays);
    printf(" done.\n");

    printf("\tSetting up the FEE beam...");
    RTS_MWAFEEInit(mwa_fee_hdf5, freq, FEE_beam, zenith_delays);
    printf(" done.\n");

    //Set up a bunch of az/za

    //Rotate by parallactic angles
    int rotation = 1;
    //Scale to zenith
    int scaling = 1;

    //Run the CUDA code
    test_RTS_CUDA_FEE_beam(num_azza,
               single_az, single_za, MWA_LAT,
               FEE_beam_zenith,
               FEE_beam,
               rotation, scaling,
               FEE_beam_gains);

    printf("Freq %.3f %.8f\n",freq/1e+6, creal(FEE_beam_gains[0]));

    out_gx_re[i] = creal(FEE_beam_gains[0]);
    out_gx_im[i] = cimag(FEE_beam_gains[0]);
    out_Dx_re[i] = creal(FEE_beam_gains[1]);
    out_Dx_im[i] = cimag(FEE_beam_gains[1]);
    out_Dy_re[i] = creal(FEE_beam_gains[2]);
    out_Dy_im[i] = cimag(FEE_beam_gains[2]);
    out_gy_re[i] = creal(FEE_beam_gains[3]);
    out_gy_im[i] = cimag(FEE_beam_gains[3]);

    RTS_freeHDFBeam(FEE_beam);
    RTS_freeHDFBeam(FEE_beam_zenith);
  }
  free(FEE_beam_gains);

  FILE *beam_values_out;
  char buff[0x100];
  snprintf(buff, sizeof(buff), "MWAFEE_beamvalues-vs-freq.txt");
  beam_values_out = fopen(buff,"w");

  for (int i = 0; i < num_freqs; i++) {
    fprintf(beam_values_out,"%.1f %.7f %.7f %.7f %.7f %.7f %.7f %.7f %.7f\n",
                        freqs_used[i],
                        out_gx_re[i], out_gx_im[i], out_Dx_re[i], out_Dx_im[i],
                        out_Dy_re[i], out_Dy_im[i], out_gy_re[i], out_gy_im[i]);
  }
  free(out_gx_re);
  free(out_gx_im);
  free(out_Dx_re);
  free(out_Dx_im);
  free(out_Dy_re);
  free(out_Dy_im);
  free(out_gy_re);
  free(out_gy_im);
  free(freqs_used);
}



int main(void)
{
  FILE *fp=NULL;
  char line[BUFSIZ];

  if ((fp=fopen("azza_values.txt","r"))==NULL) {
    printf("Read of azs,zas failed <%s>\n", "azza_values.txt");
    exit(1);
  }

  int num_azzas = 0;
  user_precision_t *azs = malloc(sizeof(user_precision_t));
  user_precision_t *zas = malloc(sizeof(user_precision_t));

  // int result;

  while(fgets(line,BUFSIZ,fp) != NULL) {

    num_azzas += 1;

    azs = realloc(azs,sizeof(user_precision_t)*num_azzas);
    zas = realloc(zas,sizeof(user_precision_t)*num_azzas);

    sscanf( line, "%lf %lf", &azs[num_azzas-1], &zas[num_azzas-1] );

  }


  char* mwa_fee_hdf5 = getenv("MWA_FEE_HDF5");

  if (mwa_fee_hdf5) {
    printf("MWA_FEE_HDF5: %s", mwa_fee_hdf5 );
    run_RTS_MWA_FEE_beam(azs,  zas, num_azzas,
                              180e+6, zenith_delays,
                              mwa_fee_hdf5,
                              "MWAFEE_beamvalues_180MHz.txt");

    RTS_MWA_beam_by_freq(mwa_fee_hdf5);
  }
  else {
    printf("MWA_FEE_HDF5 not found - not running test_RTS_FEE_beam test");
  }

  free(azs);
  free(zas);

}
