#include <math.h>
#include <stdio.h>
#include <complex.h>
#include <string.h>
#include "constants.h"
#include <stdlib.h>
// #include "read_and_write.h"

//External CUDA code being linked in
extern void test_calculate_gaussian_beam(int num_components, int num_time_steps,
     int num_freqs, float ha0, float sdec0, float cdec0,
     float fwhm_lm, float cos_theta, float sin_theta, float sin_2theta,
     float beam_ref_freq, float *freqs,
     float *beam_has, float *beam_decs,
     float _Complex *primay_beam_J00, float _Complex *primay_beam_J11);

void main(int argc, char *argv[]) {

  FILE *fp=NULL;
  char line[BUFSIZ];

  if ((fp=fopen("ha-dec_values.txt","r"))==NULL) {
    printf("Read of has,decs failed <%s>\n", "ha-dec_values.txt");
    exit(1);
  }

  int num_hadec = 0;
  float *beam_has;
  float *beam_decs;

  int result;

  while(fgets(line,BUFSIZ,fp) != NULL) {

    num_hadec += 1;

    beam_has = realloc(beam_has,sizeof(float)*num_hadec);
    beam_decs = realloc(beam_decs,sizeof(float)*num_hadec);

    result = sscanf( line, "%f %f", &beam_has[num_hadec-1], &beam_decs[num_hadec-1] );

  }

  float ha0 = 0;
  float dec0 = MWA_LAT_RAD;
  float sdec0, cdec0;
  sdec0 = sinf(MWA_LAT_RAD);
  cdec0 = cosf(MWA_LAT_RAD);

  const int num_freqs = 1;
  const int num_times = 1;

  int num_components = num_hadec;
  int num_beam_calcs = num_freqs * num_hadec;

  float freqs[1] = {100e+6};
  float beam_ref_freq = 150e+6;

  float cos_theta = 1.0;
  float sin_theta = 0.0;
  float sin_2theta = 0.0;

  float beam_FWHM_rad = 10*DD2R;
  float fwhm_lm = sinf(beam_FWHM_rad);

  float _Complex *primay_beam_J00 = malloc(num_freqs*num_times*num_hadec*sizeof(float _Complex));
  float _Complex *primay_beam_J11 = malloc(num_freqs*num_times*num_hadec*sizeof(float _Complex));

  test_calculate_gaussian_beam(num_hadec, num_times,
     num_freqs, ha0, sdec0, cdec0,
     fwhm_lm, cos_theta, sin_theta, sin_2theta,
     beam_ref_freq, freqs,
     beam_has, beam_decs,
     primay_beam_J00, primay_beam_J11);


  //Write out the ha/dec values we read in
  FILE *output;
  char buff[0x100];

  snprintf(buff, sizeof(buff), "Gaussian_beam_zenith_100MHz.txt");
  output = fopen(buff,"w");
  // //
  for (size_t coord = 0; coord < num_hadec; coord++) {
    fprintf(output, "%.8f %.8f\n",
            creal(primay_beam_J00[coord]), creal(primay_beam_J11[coord]));
  }

  fflush(output);
  fclose(output);

  ha0 = M_PI / 4;


  test_calculate_gaussian_beam(num_hadec, num_times,
     num_freqs, ha0, sdec0, cdec0,
     fwhm_lm, cos_theta, sin_theta, sin_2theta,
     beam_ref_freq, freqs,
     beam_has, beam_decs,
     primay_beam_J00, primay_beam_J11);



  snprintf(buff, sizeof(buff), "Gaussian_beam_offzenith_100MHz.txt");
  output = fopen(buff,"w");
  // //
  for (size_t coord = 0; coord < num_hadec; coord++) {
    fprintf(output, "%.8f %.8f\n",
            creal(primay_beam_J00[coord]), creal(primay_beam_J11[coord]));
  }

  fflush(output);
  fclose(output);

}//end main
