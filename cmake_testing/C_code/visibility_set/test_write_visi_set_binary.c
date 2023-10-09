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
Check the function that writes out visibilities to a binary file works
*/
void test_write_visi_set_binary() {

  //Output container
  visibility_set_t *visibility_set = malloc(sizeof(visibility_set_t));

  int num_visis = 3;

  //Do so mallocing
  visibility_set->us_metres = malloc( num_visis * sizeof(user_precision_t) );
  visibility_set->vs_metres = malloc( num_visis * sizeof(user_precision_t) );
  visibility_set->ws_metres = malloc( num_visis * sizeof(user_precision_t) );

  visibility_set->sum_visi_XX_real = malloc( num_visis * sizeof(user_precision_t) );
  visibility_set->sum_visi_XX_imag = malloc( num_visis * sizeof(user_precision_t) );
  visibility_set->sum_visi_XY_real = malloc( num_visis * sizeof(user_precision_t) );
  visibility_set->sum_visi_XY_imag = malloc( num_visis * sizeof(user_precision_t) );
  visibility_set->sum_visi_YX_real = malloc( num_visis * sizeof(user_precision_t) );
  visibility_set->sum_visi_YX_imag = malloc( num_visis * sizeof(user_precision_t) );
  visibility_set->sum_visi_YY_real = malloc( num_visis * sizeof(user_precision_t) );
  visibility_set->sum_visi_YY_imag = malloc( num_visis * sizeof(user_precision_t) );


  //make some outputs
  for (int visi = 0; visi < num_visis; visi++) {
    visibility_set->us_metres[visi] = visi*11 + 0;
    visibility_set->vs_metres[visi] = visi*11 + 1;
    visibility_set->ws_metres[visi] = visi*11 + 2;
    visibility_set->sum_visi_XX_real[visi] = visi*11 + 3;
    visibility_set->sum_visi_XX_imag[visi] = visi*11 + 4;
    visibility_set->sum_visi_XY_real[visi] = visi*11 + 5;
    visibility_set->sum_visi_XY_imag[visi] = visi*11 + 6;
    visibility_set->sum_visi_YX_real[visi] = visi*11 + 7;
    visibility_set->sum_visi_YX_imag[visi] = visi*11 + 8;
    visibility_set->sum_visi_YY_real[visi] = visi*11 + 9;
    visibility_set->sum_visi_YY_imag[visi] = visi*11 + 10;

  }

  //Generate expected output array
  int num_outputs = 33;
  user_precision_t *expec_outputs = malloc(num_outputs*sizeof(user_precision_t));
  int expec_ind = 0;

  for (int out_ind = 0; out_ind < 11; out_ind++) {
    for (int visi = 0; visi < num_visis; visi++) {
      expec_outputs[expec_ind] = visi*11 + out_ind;
      expec_ind ++;
    }
  }

  //Try a random band number
  int band_num = 5;

  //Code to be tested
  write_visi_set_binary(visibility_set, band_num, num_visis);

  //Reading in things
  FILE *fp=NULL;
  // char line[BUFSIZ];
  fp = fopen("output_visi_band05.dat","rb");

  //Assert that we wrote something and so managed to read it in
  TEST_ASSERT_NOT_NULL(fp);

  int result = 0;

  //Decode the user_precision_ts and see if it matches expectations
  user_precision_t output;
  for (int out_ind = 0; out_ind < num_outputs; out_ind++) {
    result = fread((void*)(&output), sizeof(user_precision_t), 1, fp);
    TEST_ASSERT_EQUAL_INT(1, result);
    TEST_ASSERT_EQUAL_FLOAT((float)output, (float)expec_outputs[out_ind]);
  }
}

//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_write_visi_set_binary);

    return UNITY_END();
}
