#include <unity.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <math.h>
#include "call_everybeam_c.h"
#include "test_run_lba_beam.h"
#include "woden_precision_defs.h"
#include "constants.h"

// extern void learn_cpp(double _Complex *jones, int num_coords);

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

#define TOL 1e-4

#define NUM_COORDS 31
#define NUM_DIRS NUM_COORDS*NUM_COORDS
#define NUM_TIMES 2
#define NUM_FREQS 2
#define NUM_STATIONS 2
#define D2R (M_PI/180.0)
#define RA0_LBA 277.3824204*D2R
#define DEC0_LBA 48.7461556*D2R
#define LOW_FREQ 50e+6
#define FREQ_INC 50e+6
#define TIME_RES 3*3600.0
#define MJD_LBA 4921833485.00695




void do_run_lofar_beam(const char *ms_path, bool apply_beam_norms,
                       bool rotate, bool element_only, bool iau_order,
                       double _Complex * jones) {

  const char coeff_path[] = "";
  const char element_response_model[] = "hamaker";

  int eb_status = 0;

  double mjd_sec_times[NUM_TIMES];
  for (int timei = 0; timei < NUM_TIMES; timei++) {
    mjd_sec_times[timei] = MJD_LBA + timei*TIME_RES;
  }

  double freqs[NUM_FREQS];
  
  for (int i = 0; i < NUM_FREQS; i++) {
    freqs[i] = LOW_FREQ + i*FREQ_INC;
  }

  int station_idxs[NUM_STATIONS] = {0, 9};

  eb_status = load_and_run_lofar_beam(ms_path, element_response_model,
                          coeff_path,
                          NUM_STATIONS, station_idxs,
                          NUM_DIRS, RA0_LBA, DEC0_LBA, ras, decs,
                          NUM_TIMES, mjd_sec_times, NUM_FREQS, freqs,
                          apply_beam_norms, rotate, element_only, iau_order,
                          jones);

  TEST_ASSERT_EQUAL_INT(0, eb_status);
}

void check_against_python_wrapper_values(double _Complex *jones,
                                    double *expected) {

  int num_components = NUM_DIRS;
  for (int comp = 0; comp < num_components; comp ++) {

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


void test_run_lba_telescope_rotate(void) {

  const char ms_path[] = "../../../../test_installation/everybeam/lba.MS";
  bool apply_beam_norms = false;
  bool rotate = true;
  bool element_only = false;
  bool iau_order = false;

  double _Complex *jones = malloc(4*NUM_DIRS*NUM_TIMES*NUM_FREQS*NUM_STATIONS*sizeof(double _Complex));

  do_run_lofar_beam(ms_path, apply_beam_norms,
                    rotate, element_only, iau_order, jones);

  check_against_python_wrapper_values(jones, jones_rotated);

  free(jones);
}

void test_run_lba_telescope_rotate_normed(void) {

  const char ms_path[] = "../../../../test_installation/everybeam/lba.MS";
  bool apply_beam_norms = true;
  bool rotate = true;
  bool element_only = false;
  bool iau_order = false;

  double _Complex *jones = malloc(4*NUM_DIRS*NUM_TIMES*NUM_FREQS*NUM_STATIONS*sizeof(double _Complex));

  do_run_lofar_beam(ms_path, apply_beam_norms,
                    rotate, element_only, iau_order, jones);

  check_against_python_wrapper_values(jones, jones_rotated_normed);

  free(jones);
}

void test_run_lba_telescope_reordered(void) {

  const char ms_path[] = "../../../../test_installation/everybeam/lba.MS";
  bool apply_beam_norms = false;
  bool rotate = false;
  bool element_only = false;
  bool iau_order = true;

  double _Complex *jones = malloc(4*NUM_DIRS*NUM_TIMES*NUM_FREQS*NUM_STATIONS*sizeof(double _Complex));

  do_run_lofar_beam(ms_path, apply_beam_norms,
                    rotate, element_only, iau_order, jones);

  check_against_python_wrapper_values(jones, jones_reordered);

  free(jones);
}



//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_run_lba_telescope_rotate);
    RUN_TEST(test_run_lba_telescope_rotate_normed);
    RUN_TEST(test_run_lba_telescope_reordered);

    return UNITY_END();
}
