#include <unity.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <erfa.h>

#include "constants.h"
#include "FEE_primary_beam.h"
#include "azza_radec_nside101.h"

#define FREQ 200000000

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

//External CUDA code we're linking in
extern void test_RTS_CUDA_FEE_beam(int num_components,
           user_precision_t *azs, user_precision_t *zas, double latitude,
           RTS_MWA_FEE_beam_t *FEE_beam_zenith,
           RTS_MWA_FEE_beam_t *FEE_beam,
           int rotation, int scaling,
           user_precision_complex_t *FEE_beam_gains);

#define UNITY_INCLUDE_FLOAT

void test_MWA_FEE_beam_nside101() {

  char* mwa_fee_hdf5 = getenv("MWA_FEE_HDF5");

  if (mwa_fee_hdf5) {
    printf("MWA_FEE_HDF5: %s\n", mwa_fee_hdf5 );

    int nside = 101;

    //just say the sky hasn't moved with time for this test
    int num_times = 1;
    // int num_freqs = 1;

    int num_components = nside*nside;
    int num_azza = num_components*num_times;

    user_precision_complex_t *FEE_beam_gains = malloc(4*num_azza*sizeof(user_precision_complex_t));

    // int norm = 1;
    user_precision_t delays_zen[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    user_precision_t delays[] = {0, 2, 4, 8, 2, 4, 8, 12, 4, 8, 12, 16, 8, 12, 16, 20};

    // user_precision_t delays[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                      // 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};


    int rotation=1;
    int scaling=1;
    double freq = FREQ;

    //Call the C code to interrogate the hdf5 file and set beam things up
    RTS_MWA_FEE_beam_t *FEE_beam = malloc(sizeof(RTS_MWA_FEE_beam_t));

    //Get a zenith pointing beam for normalisation purposes
    RTS_MWA_FEE_beam_t *FEE_beam_zenith = malloc(sizeof(RTS_MWA_FEE_beam_t));

    printf("\n\tSetting up the zenith FEE beam...");
    RTS_MWAFEEInit(mwa_fee_hdf5,
                   freq, FEE_beam_zenith, delays_zen);
    printf(" done.\n");

    printf("\tSetting up the FEE beam...");
    RTS_MWAFEEInit(mwa_fee_hdf5,
                   freq, FEE_beam, delays);
    printf(" done.\n");

    test_RTS_CUDA_FEE_beam(num_azza,
               nside101_azs, nside101_zas, MWA_LAT_RAD,
               FEE_beam_zenith,
               FEE_beam,
               rotation, scaling,
               FEE_beam_gains);

    FILE *output_text;

    char filename[100];
    sprintf(filename, "MWA_FEE_gains_azza%02d.txt", num_components);

    output_text = fopen(filename,"w");

    int count = 0;

    for (int comp = 0; comp < num_components; comp ++) {

      int beam_ind = 4*comp;

      fprintf(output_text,"%.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.1f\n",
       nside101_azs[comp], nside101_zas[comp],
       creal(FEE_beam_gains[beam_ind + 0]), cimag(FEE_beam_gains[beam_ind + 0]),
       creal(FEE_beam_gains[beam_ind + 1]), cimag(FEE_beam_gains[beam_ind + 1]),
       creal(FEE_beam_gains[beam_ind + 2]), cimag(FEE_beam_gains[beam_ind + 2]),
       creal(FEE_beam_gains[beam_ind + 3]), cimag(FEE_beam_gains[beam_ind + 3]),
       freq );

      // fprintf(output_text,"%.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.1f\n",
      //  nside101_azs[comp], nside101_zas[comp],
      //  cimag(FEE_beam_gains[beam_ind + 0]), creal(FEE_beam_gains[beam_ind + 0]),
      //  cimag(FEE_beam_gains[beam_ind + 1]), creal(FEE_beam_gains[beam_ind + 1]),
      //  cimag(FEE_beam_gains[beam_ind + 2]), creal(FEE_beam_gains[beam_ind + 2]),
      //  cimag(FEE_beam_gains[beam_ind + 3]), creal(FEE_beam_gains[beam_ind + 3]),
      //  freq );

       count ++;
    }

    fflush(output_text);
    fclose(output_text);

    RTS_freeHDFBeam(FEE_beam);
    RTS_freeHDFBeam(FEE_beam_zenith);

    free(FEE_beam);
    free(FEE_beam_zenith);

    free(FEE_beam_gains);
  }
  else {
    printf("MWA_FEE_HDF5 not found - not running comp_MWA_analytic_to_FEE\n");
  }

}


void test_MWA_FEE_beam_single_azza() {

  char* mwa_fee_hdf5 = getenv("MWA_FEE_HDF5");

  if (mwa_fee_hdf5) {
    printf("MWA_FEE_HDF5: %s\n", mwa_fee_hdf5 );

    int nside = 1;

    //just say the sky hasn't moved with time for this test
    int num_times = 1;
    // int num_freqs = 1;

    int num_components = nside*nside;
    int num_azza = num_components*num_times;

    user_precision_complex_t *FEE_beam_gains = malloc(4*num_azza*sizeof(user_precision_complex_t));

    // int norm = 1;
    user_precision_t delays_zen[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    user_precision_t delays[] = {6, 4, 2, 0, 8, 6, 4, 2,
                      10, 8, 6, 4, 12, 10, 8, 6};

    // user_precision_t delays[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                      // 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};


    int rotation=1;
    int scaling=1;
    double freq = FREQ;

    //Call the C code to interrogate the hdf5 file and set beam things up
    RTS_MWA_FEE_beam_t *FEE_beam = malloc(sizeof(RTS_MWA_FEE_beam_t));

    //Get a zenith pointing beam for normalisation purposes
    RTS_MWA_FEE_beam_t *FEE_beam_zenith = malloc(sizeof(RTS_MWA_FEE_beam_t));

    printf("\n\tSetting up the zenith FEE beam...");
    RTS_MWAFEEInit(mwa_fee_hdf5,
                   freq, FEE_beam_zenith, delays_zen);
    printf(" done.\n");

    printf("\tSetting up the FEE beam...");
    RTS_MWAFEEInit(mwa_fee_hdf5,
                   freq, FEE_beam, delays);
    printf(" done.\n");

    // int num_azza = 1;

    user_precision_t azs[1] = {3.7997899};
    user_precision_t zas[1] = {0.3080986};

    test_RTS_CUDA_FEE_beam(num_azza,
               azs, zas, MWA_LAT_RAD,
               FEE_beam_zenith,
               FEE_beam,
               rotation, scaling,
               FEE_beam_gains);

    printf("gx re, im %.8f %.8f\n",creal(FEE_beam_gains[0]), cimag(FEE_beam_gains[0]) );
    printf("Dx re, im %.8f %.8f\n",creal(FEE_beam_gains[1]), cimag(FEE_beam_gains[1]) );
    printf("Dy re, im %.8f %.8f\n",creal(FEE_beam_gains[2]), cimag(FEE_beam_gains[2]) );
    printf("gy re, im %.8f %.8f\n",creal(FEE_beam_gains[3]), cimag(FEE_beam_gains[3]) );

    RTS_freeHDFBeam(FEE_beam);
    RTS_freeHDFBeam(FEE_beam_zenith);

    free(FEE_beam);
    free(FEE_beam_zenith);

    free(FEE_beam_gains);
  }
  else {
    printf("MWA_FEE_HDF5 not found - not running comp_MWA_analytic_to_FEE\n");
  }

}

//Run the test with unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_MWA_FEE_beam_nside101);
    // RUN_TEST(test_MWA_FEE_beam_single_azza);

    return UNITY_END();
}
