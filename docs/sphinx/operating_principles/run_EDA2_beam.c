#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <assert.h>

#include "constants.h"

void test_analytic_dipole_beam(int num_components,
     int num_time_steps, int num_freqs,
     float *azs, float *zas, float *freqs,
     float _Complex *analy_beam_X, float _Complex *analy_beam_Y);

int main(int argc, char *argv[]) {

  FILE *fp=NULL;
  char line[BUFSIZ];

  if ((fp=fopen("azza_values.txt","r"))==NULL) {
    printf("Read of azs,zas failed <%s>\n", "az-za_values.txt");
    exit(1);
  }

    int num_azza = 0;
    float *azs = malloc(sizeof(float));
    float *zas = malloc(sizeof(float));

    // int result;

    while(fgets(line,BUFSIZ,fp) != NULL) {

      num_azza += 1;

      azs = realloc(azs,sizeof(float)*num_azza);
      zas = realloc(zas,sizeof(float)*num_azza);

      sscanf( line, "%f %f", &azs[num_azza-1], &zas[num_azza-1] );

    }

  float freq = 70e+6;
  int num_freqs = 1;
  int num_time_steps = 1;

  float freqs[1] = {freq};

  float _Complex *analy_beam_X=NULL;
  analy_beam_X = malloc(num_freqs*num_time_steps*num_azza*sizeof(float _Complex));

  float _Complex *analy_beam_Y=NULL;
  analy_beam_Y = malloc(num_freqs*num_time_steps*num_azza*sizeof(float _Complex));

  test_analytic_dipole_beam(num_azza,
       num_time_steps, num_freqs,
       azs, zas, freqs,
       analy_beam_X, analy_beam_Y);


  FILE *output;
  char buff[0x100];

  snprintf(buff, sizeof(buff), "EDA2_beam_70MHz.txt");
  output = fopen(buff,"w");
  // //
  for (size_t coord = 0; coord < num_azza; coord++) {
    fprintf(output, "%.8f %.8f\n",
            creal(analy_beam_X[coord]), creal(analy_beam_Y[coord]));
  }

  fflush(output);
  fclose(output);

  // printf("Here be the end\n");
}
