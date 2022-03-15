#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <erfa.h>

#include "constants.h"
#include "woden_struct_defs.h"
#include "FEE_primary_beam.h"

//External CUDA code we're linking in
extern void test_RTS_calculate_MWA_analytic_beam(int num_components,
     int num_time_steps, int num_freqs,
     user_precision_t *azs, user_precision_t *zas, user_precision_t *delays,
     double latitude, int norm,
     double *beam_has, double *beam_decs, double *freqs,
     user_precision_complex_t *gxs, user_precision_complex_t *Dxs,
     user_precision_complex_t *Dys, user_precision_complex_t *gys);

#define UNITY_INCLUDE_FLOAT

//Different delays settings, which control the pointing of the MWA beam
user_precision_t zenith_delays[16] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                           0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

void run_RTS_MWA_analy_beam(user_precision_t *azs,  user_precision_t *zas, int num_azza,
                          double freq, user_precision_t *delays,
                          char *outname) {

  double *beam_has = malloc(num_azza*sizeof(double));
  double *beam_decs = malloc(num_azza*sizeof(double));

  double ha, dec, el;
  for (int comp = 0; comp < num_azza; comp++) {
    el = M_PI/2.0 - (double)zas[comp];

    eraAe2hd((double)azs[comp], el, MWA_LAT_RAD, &ha, &dec);

    beam_has[comp] = ha;
    beam_decs[comp] = dec;

  }


  double *freqs = malloc(sizeof(double));
  freqs[0] = freq;

  user_precision_complex_t *gxs = malloc(num_azza*sizeof(user_precision_complex_t));
  user_precision_complex_t *Dxs = malloc(num_azza*sizeof(user_precision_complex_t));
  user_precision_complex_t *Dys = malloc(num_azza*sizeof(user_precision_complex_t));
  user_precision_complex_t *gys = malloc(num_azza*sizeof(user_precision_complex_t));

  test_RTS_calculate_MWA_analytic_beam(num_azza, 1, 1, azs, zas, zenith_delays,
       MWA_LAT_RAD, 1, beam_has, beam_decs, freqs, gxs, Dxs, Dys, gys);

  FILE *beam_values_out;
  beam_values_out = fopen(outname,"w");

  for (size_t comp = 0; comp < num_azza; comp++) {
    fprintf(beam_values_out, "%.5f %.5f %.7f %.7f %.7f %.7f %.7f %.7f %.7f %.7f\n",
    azs[comp], zas[comp],
    creal(gxs[comp]), cimag(gxs[comp]),
    creal(Dxs[comp]), cimag(Dxs[comp]),
    creal(Dys[comp]), cimag(Dys[comp]),
    creal(gys[comp]), cimag(gys[comp]));
  }

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

  run_RTS_MWA_analy_beam(azs,  zas, num_azzas,
                         180e+6, zenith_delays,
                         "MWAanaly_beamvalues_180MHz.txt");

  free(azs);
  free(zas);

}
