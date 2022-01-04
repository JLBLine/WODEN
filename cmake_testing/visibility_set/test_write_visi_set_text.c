#include <unity.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "constants.h"
#include "woden_struct_defs.h"
#include "visibility_set.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

#define UNITY_INCLUDE_FLOAT

/*
Check the function that writes out visibilities to a text file works
*/
void test_write_visi_set_text() {

  //Settings struct
  woden_settings_t *woden_settings = malloc(sizeof(woden_settings_t));
  woden_settings->num_time_steps = 2;
  woden_settings->num_freqs = 3;
  woden_settings->num_baselines = 5;

  int num_visis = woden_settings->num_time_steps*woden_settings->num_freqs*woden_settings->num_baselines;

  //Output container
  visibility_set_t *visibility_set = malloc(sizeof(visibility_set_t));

  //Do so mallocing
  visibility_set->us_metres = malloc( num_visis * sizeof(user_precision_t) );
  visibility_set->vs_metres = malloc( num_visis * sizeof(user_precision_t) );
  visibility_set->ws_metres = malloc( num_visis * sizeof(user_precision_t) );

  visibility_set->sum_visi_XX_real = malloc( num_visis * sizeof(user_precision_t) );
  visibility_set->sum_visi_XX_imag = malloc( num_visis * sizeof(user_precision_t) );

  //Fill the outputs with some basic inputs
  //Make the five arrays slightly different do we know we aren't writing
  //out one array to all outputs
  for ( int time_step = 0; time_step < woden_settings->num_time_steps; time_step++ ) {
    for ( int freq_step = 0; freq_step < woden_settings->num_freqs; freq_step++ ) {
      for (int baseline = 0; baseline < woden_settings->num_baselines; baseline++) {
        int step = woden_settings->num_baselines*(time_step*woden_settings->num_freqs + freq_step);

        visibility_set->us_metres[step + baseline] = step + baseline + 1;
        visibility_set->vs_metres[step + baseline] = step + baseline + 2;
        visibility_set->ws_metres[step + baseline] = step + baseline + 3;
        visibility_set->sum_visi_XX_real[step + baseline] = step + baseline + 4;
        visibility_set->sum_visi_XX_imag[step + baseline] = step + baseline + 5;

      }
    }
  }

  //Try a random band number
  int band_num = 7;

  write_visi_set_text(visibility_set, band_num,
                      woden_settings);

  //Now try reading in the output we just made
  FILE *fp=NULL;
  char line[BUFSIZ];

  if ((fp=fopen("output_visi_band07.txt","r"))==NULL) {
    printf("Failed to read output_visi_band07.txt");
    exit(1);
  }

  //Things to store the outputs in
  user_precision_t *read_us_metres = malloc( num_visis * sizeof(user_precision_t) );
  user_precision_t *read_vs_metres = malloc( num_visis * sizeof(user_precision_t) );
  user_precision_t *read_ws_metres = malloc( num_visis * sizeof(user_precision_t) );

  user_precision_t *read_sum_visi_XX_real = malloc( num_visis * sizeof(user_precision_t) );
  user_precision_t *read_sum_visi_XX_imag = malloc( num_visis * sizeof(user_precision_t) );

  //Read it in, with the precision dictated at compile time
  int line_index = 0;
  while(fgets(line, BUFSIZ, fp) != NULL) {
    #ifdef DOUBLE_PRECISION

      sscanf( line, "%lf %lf %lf %lf %lf", &read_us_metres[line_index],
                                           &read_vs_metres[line_index],
                                           &read_ws_metres[line_index],
                                           &read_sum_visi_XX_real[line_index],
                                           &read_sum_visi_XX_imag[line_index]);
    #else
      sscanf( line, "%f %f %f %f %f", &read_us_metres[line_index],
                                           &read_vs_metres[line_index],
                                           &read_ws_metres[line_index],
                                           &read_sum_visi_XX_real[line_index],
                                           &read_sum_visi_XX_imag[line_index]);
    #endif

    line_index += 1;
  }

  //Check the outputs match expectation
  for (int visi_ind = 0; visi_ind < num_visis; visi_ind++) {
    TEST_ASSERT_EQUAL_FLOAT((float)(visi_ind + 1), (float)read_us_metres[visi_ind]);
    TEST_ASSERT_EQUAL_FLOAT((float)(visi_ind + 2), (float)read_vs_metres[visi_ind]);
    TEST_ASSERT_EQUAL_FLOAT((float)(visi_ind + 3), (float)read_ws_metres[visi_ind]);
    TEST_ASSERT_EQUAL_FLOAT((float)(visi_ind + 4), (float)read_sum_visi_XX_real[visi_ind]);
    TEST_ASSERT_EQUAL_FLOAT((float)(visi_ind + 5), (float)read_sum_visi_XX_imag[visi_ind]);
  }

  free(visibility_set->us_metres);
  free(visibility_set->vs_metres);
  free(visibility_set->ws_metres);
  free(visibility_set->sum_visi_XX_real);
  free(visibility_set->sum_visi_XX_imag);
}

//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_write_visi_set_text);

    return UNITY_END();
}
