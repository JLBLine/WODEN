#include <math.h>
#include <stdlib.h>
#include <complex.h>

#include "constants.h"
#include "woden_struct_defs.h"
#include "FEE_primary_beam.h"
#include "woden_precision_defs.h"

//CUDA code we are testing
extern void multifreq_get_MWAFEE_normalisation(beam_settings_t *beam_settings);

extern void free_FEE_primary_beam_from_GPU(RTS_MWA_FEE_beam_t *primary_beam);

//CUDA code we are testing
extern void test_run_and_map_multifreq_calc_CUDA_FEE_beam(beam_settings_t *beam_settings,
    user_precision_t *azs, user_precision_t *zas, double latitude,
    user_precision_complex_t *primay_beam_J00,
    user_precision_complex_t *primay_beam_J01,
    user_precision_complex_t *primay_beam_J10,
    user_precision_complex_t *primay_beam_J11,
    int num_freqs, int num_components, int num_times,
    int rotation, int scaling);

void compare_interp_MWA_FEE(char* mwa_fee_hdf5, char* mwa_fee_hdf5_interp){

  user_precision_t *single_az = (user_precision_t *)malloc(sizeof(user_precision_t));
  user_precision_t *single_za = (user_precision_t *)malloc(sizeof(user_precision_t));

  single_az[0] = 55.0 * (M_PI / 180.0);
  single_za[0] = 60.0 * (M_PI / 180.0);

  int num_freqs = 72;

  double freq_res = 160e+3;
  double base_low_freq = 168.0e+6;

  double *freqs = malloc(num_freqs*sizeof(double));

  for (int i = 0; i < num_freqs; i++) {
    freqs[i] = base_low_freq + i*freq_res;
  }

  beam_settings_t *beam_settings = malloc(sizeof(beam_settings_t));
  woden_settings_t *woden_settings = malloc(sizeof(woden_settings_t));

  woden_settings->base_low_freq = base_low_freq;
  woden_settings->num_freqs = num_freqs;
  woden_settings->hdf5_beam_path = mwa_fee_hdf5_interp;

  beam_settings->num_MWAFEE = num_freqs;

  for(int i=0;i<16;i++) {
      woden_settings->FEE_ideal_delays[i] = 0.0;
  }

  int status = multifreq_RTS_MWAFEEInit(beam_settings,  woden_settings, freqs);

  //Send them to GPU and calculate normalisations
  multifreq_get_MWAFEE_normalisation(beam_settings);

  //Test is setup for a single time step
  int num_times = 1;
  int num_components = 1;
  int rotation = 1;
  int scaling = 1;

  user_precision_complex_t *out_gx = malloc(num_freqs*sizeof(user_precision_complex_t));
  user_precision_complex_t *out_Dx = malloc(num_freqs*sizeof(user_precision_complex_t));
  user_precision_complex_t *out_Dy = malloc(num_freqs*sizeof(user_precision_complex_t));
  user_precision_complex_t *out_gy = malloc(num_freqs*sizeof(user_precision_complex_t));

  printf("WHAT %.6f\n",single_za[0] );


  //Test the CUDA code
  test_run_and_map_multifreq_calc_CUDA_FEE_beam(beam_settings,
    single_az, single_za, MWA_LAT_RAD,
    out_gx, out_Dx, out_Dy, out_gy,
    num_freqs, num_components, num_times,
    rotation, scaling);

  FILE *beam_values_out;
  char buff[0x100];
  snprintf(buff, sizeof(buff), "MWAFEEinterp_beamvalues-vs-freq.txt");
  beam_values_out = fopen(buff,"w");

  for (int i = 0; i < num_freqs; i++) {
    fprintf(beam_values_out,"%.8e %.7f %.7f %.7f %.7f %.7f %.7f %.7f %.7f\n",
                        freqs[i], creal(out_gx[i]), cimag(out_gx[i]),
                        creal(out_Dx[i]), cimag(out_Dx[i]),
                        creal(out_Dy[i]), cimag(out_Dy[i]),
                        creal(out_gy[i]), cimag(out_gy[i]));
  }


  free(out_gx);
  free(out_Dx);
  free(out_Dy);
  free(out_gy);

  for (int freq_ind = 0; freq_ind < beam_settings->num_MWAFEE; freq_ind++) {

    RTS_MWA_FEE_beam_t *FEE_beam = &beam_settings->FEE_beams[freq_ind];
    RTS_MWA_FEE_beam_t *FEE_beam_zenith = &beam_settings->FEE_beam_zeniths[freq_ind];

    RTS_freeHDFBeam(FEE_beam);
    RTS_freeHDFBeam(FEE_beam_zenith);

  }

}



int main(void)
{


  char* mwa_fee_hdf5 = getenv("MWA_FEE_HDF5");
  char* mwa_fee_hdf5_interp = getenv("MWA_FEE_HDF5_INTERP");

  if (mwa_fee_hdf5 && mwa_fee_hdf5_interp) {
    printf("MWA_FEE_HDF5: %s\n", mwa_fee_hdf5 );
    printf("MWA_FEE_HDF5_INTERP: %s\n", mwa_fee_hdf5_interp );

    compare_interp_MWA_FEE(mwa_fee_hdf5, mwa_fee_hdf5_interp);
  }
  else {
    printf("MWA_FEE_HDF5 not found - not running");
  }


}
