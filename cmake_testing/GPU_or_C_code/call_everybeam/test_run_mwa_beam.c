#include <unity.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <math.h>
#include "call_everybeam_c.h"
// #include "mwa_jones_values.h"

#include "woden_precision_defs.h"
#include "woden_struct_defs.h"

// #include "run_hyperbeam_common.h"
// #include <mwa_hyperbeam.h>
#include "azza_para_nside051.h"
#include "test_run_mwa_beam.h"

// #include "erfa.h"

// nside051_azs

// extern void learn_cpp(double _Complex *jones, int num_coords);

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

#define TOL 1e-5

#define NUM_COORDS 51
#define NUM_DIRS NUM_COORDS*NUM_COORDS
#define NUM_TIMES 2
#define NUM_FREQS 3
#define NUM_STATIONS 1
#define D2R (M_PI/180.0)
#define PATCH_WIDTH 40.0*D2R
#define RA0_MWA 0.0*D2R
#define DEC0_MWA -26.703319405555554*D2R
#define LOW_FREQ 100e+6
#define FREQ_INC 50e+6
#define TIME_RES 3*3600.0
#define MJD_mwa 4875418129.0

#define MWA_LAT -26.703319405555554
#define MWA_LONG 116.67081523611111


void do_run_mwa_beam(const char *ms_path, bool apply_beam_norms,
                       bool rotate, bool element_only, double _Complex * jones) {

  const char coeff_path[] = "/home/jack-line/software/mwa_beam_files/mwa_full_embedded_element_pattern.h5";
  const char element_response_model[] = "MWA";
  bool iau_order = true;

  int eb_status = 0;

  double mjd_sec_times[NUM_TIMES]; 
  for (int timei = 0; timei < NUM_TIMES; timei++) {
    mjd_sec_times[timei] = MJD_mwa + timei*TIME_RES;
  }

  // int num_freqs = 1;
  double freqs[NUM_FREQS];
  
  for (int i = 0; i < NUM_FREQS; i++) {
    freqs[i] = LOW_FREQ + i*FREQ_INC;
  }

  // int num_stations = 1;
  int station_idxs[NUM_STATIONS] = {0};


  // double ras[NUM_DIRS];
  // double decs[NUM_DIRS];

  // make_radec(ras, decs);

  // double para_angles[NUM_TIMES*NUM_DIRS];

  // for (int timei = 0; timei < NUM_TIMES; timei++) {
  //   for (int diri = 0; diri < NUM_DIRS; diri++) {
  //     para_angles[diri*NUM_TIMES + timei] = 0.0;
  //   }
  // }

  //Iternally, the load_and_run_mwa_beam function expects za,az,para
  //to have time fastest changing, then direction. So do a cheeky reorder
  //of the input arrays
  double para_angles[NUM_TIMES*NUM_DIRS];
  double azs[NUM_TIMES*NUM_DIRS];
  double zas[NUM_TIMES*NUM_DIRS];
  for (int timei = 0; timei < NUM_TIMES; timei++) {
    for (int diri = 0; diri < NUM_DIRS; diri++) {
      para_angles[diri*NUM_TIMES + timei] = nside051_paras[timei*NUM_DIRS + diri];
      azs[diri*NUM_TIMES + timei] = nside051_azs[timei*NUM_DIRS + diri];
      zas[diri*NUM_TIMES + timei] = nside051_zas[timei*NUM_DIRS + diri];
    }
  }

  eb_status = load_and_run_mwa_beam(ms_path, element_response_model,
                          coeff_path,
                          NUM_STATIONS, station_idxs,
                          NUM_DIRS, RA0_MWA, DEC0_MWA,
                          azs, zas, 
                          para_angles,
                          NUM_TIMES, mjd_sec_times, NUM_FREQS, freqs,
                          apply_beam_norms, rotate, element_only, iau_order,
                          jones);

  TEST_ASSERT_EQUAL_INT(0, eb_status);

  FILE *beam_values_out;
  // char buff[0x100];

  // #ifdef DOUBLE_PRECISION
  // if (rotate == 1) {
  //   snprintf(buff, sizeof(buff), "%s_rot_double.txt", outname);
  // } else {
  //   snprintf(buff, sizeof(buff), "%s_double.txt", outname);
  // }

  // #else
  //     snprintf(buff, sizeof(buff), "%s_float.txt", outname);
  // #endif
  // //
  beam_values_out = fopen("mwa_everybeam_values.txt","w");

  int num_components = NUM_DIRS;
  int num_freqs = NUM_FREQS;
  int num_times = NUM_TIMES;

  for (int time = 0; time < num_times; time ++) {
    for (int freq = 0; freq < num_freqs; freq ++) {
      for (int comp = 0; comp < num_components; comp ++) {

        // int beam_ind = num_freqs*time*num_components + num_components*freq + comp;
        int beam_ind = 4*(num_freqs*time*num_components + num_components*freq + comp);
        int coord_ind = comp*num_times + time;

        fprintf(beam_values_out,"%.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.1f\n",
         nside051_azs[coord_ind], nside051_zas[coord_ind],
         creal(jones[beam_ind + 0]), cimag(jones[beam_ind + 0]),
         creal(jones[beam_ind + 1]), cimag(jones[beam_ind + 1]),
         creal(jones[beam_ind + 2]), cimag(jones[beam_ind + 2]),
         creal(jones[beam_ind + 3]), cimag(jones[beam_ind + 3]),
         freqs[freq] );

      }
    }
  }
}

void check_against_hyperbeam_values(double _Complex *jones,
                                    double *expected) {

  int num_components = NUM_DIRS;
  for (int comp = 0; comp < num_components; comp ++) {

      // int beam_ind = MAX_POLS*(num_freqs*time*num_components + num_components*freq_ind + comp);

      int beam_ind = 4*comp;

      int expected_base = 2*MAX_POLS*comp;

      TEST_ASSERT_DOUBLE_WITHIN(TOL, expected[expected_base+0], creal(jones[beam_ind+0]) );
      TEST_ASSERT_DOUBLE_WITHIN(TOL, expected[expected_base+1], cimag(jones[beam_ind+0]) );
      TEST_ASSERT_DOUBLE_WITHIN(TOL, expected[expected_base+2], creal(jones[beam_ind+1]) );
      TEST_ASSERT_DOUBLE_WITHIN(TOL, expected[expected_base+3], cimag(jones[beam_ind+1]) );
      TEST_ASSERT_DOUBLE_WITHIN(TOL, expected[expected_base+4], creal(jones[beam_ind+2]) );
      TEST_ASSERT_DOUBLE_WITHIN(TOL, expected[expected_base+5], cimag(jones[beam_ind+2]) );
      TEST_ASSERT_DOUBLE_WITHIN(TOL, expected[expected_base+6], creal(jones[beam_ind+3]) );
      TEST_ASSERT_DOUBLE_WITHIN(TOL, expected[expected_base+7], cimag(jones[beam_ind+3]) );

  }
}

void test_run_mwa_telescope(void) {

  const char ms_path[] = "../../../../test_installation/everybeam/MWA-single-timeslot.ms";
  bool apply_beam_norms = false;
  bool rotate = true;
  bool element_only = false;

  double _Complex *jones = malloc(4*NUM_DIRS*NUM_TIMES*NUM_FREQS*NUM_STATIONS*sizeof(double _Complex));

  //This calcs all the different freqs and times
  do_run_mwa_beam(ms_path, apply_beam_norms,
                    rotate, element_only, jones);

  int stripe = 0;
  int time_ind = 0;
  int freq_ind = 0;

  stripe = MAX_POLS*(NUM_FREQS*NUM_DIRS*time_ind + NUM_DIRS*freq_ind);
  check_against_hyperbeam_values(jones + stripe, hyper_jones_100_time0);

  freq_ind = 1;
  stripe = MAX_POLS*(NUM_FREQS*NUM_DIRS*time_ind + NUM_DIRS*freq_ind);
  check_against_hyperbeam_values(jones + stripe, hyper_jones_150_time0);

  freq_ind = 2;
  stripe = MAX_POLS*(NUM_FREQS*NUM_DIRS*time_ind + NUM_DIRS*freq_ind);
  check_against_hyperbeam_values(jones + stripe, hyper_jones_200_time0);

  time_ind = 1;
  freq_ind = 0;
  stripe = MAX_POLS*(NUM_FREQS*NUM_DIRS*time_ind + NUM_DIRS*freq_ind);
  check_against_hyperbeam_values(jones + stripe, hyper_jones_100_time1);

  freq_ind = 1;
  stripe = MAX_POLS*(NUM_FREQS*NUM_DIRS*time_ind + NUM_DIRS*freq_ind);
  check_against_hyperbeam_values(jones + stripe, hyper_jones_150_time1);

  freq_ind = 2;
  stripe = MAX_POLS*(NUM_FREQS*NUM_DIRS*time_ind + NUM_DIRS*freq_ind);
  check_against_hyperbeam_values(jones + stripe, hyper_jones_200_time1);

  free(jones);
}





//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_run_mwa_telescope);

    return UNITY_END();
}
