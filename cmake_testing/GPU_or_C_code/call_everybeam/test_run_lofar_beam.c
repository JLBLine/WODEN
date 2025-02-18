#include <unity.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <math.h>
#include "call_everybeam_c.h"
#include "lofar_jones_values.h"

// extern void learn_cpp(double _Complex *jones, int num_coords);

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

#define TOL 5e-5

#define NUM_COORDS 5
#define NUM_DIRS NUM_COORDS*NUM_COORDS
#define NUM_TIMES 2
#define NUM_FREQS 2
#define NUM_STATIONS 2
#define D2R (M_PI/180.0)
#define PATCH_WIDTH 40.0*D2R
#define RA0_HBA 123.4002825*D2R
#define DEC0_HBA 48.21738361*D2R
#define LOW_FREQ 160e+6
#define FREQ_INC 50e+6
#define TIME_RES 3*3600.0
#define MJD_HBA 4929192872.0

void make_radec(double * ras, double * decs) {
  

  double ra_inc = PATCH_WIDTH / (NUM_COORDS - 1);
  double dec_inc = PATCH_WIDTH / (NUM_COORDS - 1);

  double ra, dec;

  for (int rai = 0; rai < NUM_COORDS; rai++) {
    for (int deci = 0; deci < NUM_COORDS; deci++) {
      ra = RA0_HBA - PATCH_WIDTH + rai * ra_inc;
      dec = DEC0_HBA - PATCH_WIDTH + deci * dec_inc;

      ras[rai * NUM_COORDS + deci] = ra;
      decs[rai * NUM_COORDS + deci] = dec;

      // printf("ra: %.2f, dec: %.2f\n", ra/D2R, dec/D2R);

    }
  }
}


void do_run_lofar_beam(const char *ms_path, bool apply_beam_norms,
                       bool rotate, bool element_only, double _Complex * jones) {

  const char coeff_path[] = "";
  const char element_response_model[] = "hamaker";

  int eb_status = 0;

  // double RA0 = 0.0;
  // double DEC0 = 90.0*D2R;

  // int num_coords = 3;
  // int num_dirs = num_coords*num_coords;

  // int num_times = 1;
  double mjd_sec_times[NUM_TIMES]; //= {4891507200.006634};
  for (int timei = 0; timei < NUM_TIMES; timei++) {
    mjd_sec_times[timei] = MJD_HBA + timei*TIME_RES;
  }

  // int num_freqs = 1;
  double freqs[NUM_FREQS];
  
  for (int i = 0; i < NUM_FREQS; i++) {
    freqs[i] = LOW_FREQ + i*FREQ_INC;
  }

  // int num_stations = 1;
  int station_idxs[NUM_STATIONS] = {0};


  double ras[NUM_DIRS];
  double decs[NUM_DIRS];

  make_radec(ras, decs);

  eb_status = load_and_run_lofar_beam(ms_path, element_response_model,
                          coeff_path,
                          NUM_STATIONS, station_idxs,
                          NUM_DIRS, RA0_HBA, DEC0_HBA, ras, decs,
                          NUM_TIMES, mjd_sec_times, NUM_FREQS, freqs,
                          apply_beam_norms, rotate, element_only, jones);

  TEST_ASSERT_EQUAL_INT(0, eb_status);
}

void check_beam_values(double _Complex *jones, const double *expected_values) {

  // for (int i = 0; i < NUM_DIRS*NUM_TIMES*NUM_FREQS*NUM_STATIONS; i++) {
  //   printf("jones[%d]: %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n",
  //          i, creal(jones[4*i]), cimag(jones[4*i]),
  //             creal(jones[4*i+1]), cimag(jones[4*i+1]),
  //             creal(jones[4*i+2]), cimag(jones[4*i+2]),
  //             creal(jones[4*i+3]), cimag(jones[4*i+3]));
  //   printf("expec[%d]: %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n",
  //          i, creal(expected_values[8*i+0]), cimag(expected_values[8*i+1]),
  //             creal(expected_values[8*i+2]), cimag(expected_values[8*i+3]),
  //             creal(expected_values[8*i+4]), cimag(expected_values[8*i+5]),
  //             creal(expected_values[8*i+6]), cimag(expected_values[8*i+7]));
  // }

  for (int i = 0; i < 4*NUM_DIRS*NUM_TIMES*NUM_FREQS*NUM_STATIONS; i++) {
    TEST_ASSERT_FLOAT_WITHIN(TOL, creal(jones[i]),  expected_values[2*i]);
    TEST_ASSERT_FLOAT_WITHIN(TOL, cimag(jones[i]), expected_values[2*i+1]);
  }
}


void test_run_hba_telescope(void) {

  const char ms_path[] = "../../../../test_installation/everybeam/LOFAR_HBA_MOCK.ms";
  bool apply_beam_norms = false;
  bool rotate = false;
  bool element_only = false;

  double _Complex *jones = malloc(4*NUM_DIRS*NUM_TIMES*NUM_FREQS*NUM_STATIONS*sizeof(double _Complex));

  do_run_lofar_beam(ms_path, apply_beam_norms,
                    rotate, element_only, jones);

  check_beam_values(jones, lofar_hba_jones);

  free(jones);
}

void test_run_hba_telescope_rotate(void) {

  const char ms_path[] = "../../../../test_installation/everybeam/LOFAR_HBA_MOCK.ms";
  bool apply_beam_norms = false;
  bool rotate = true;
  bool element_only = false;

  double _Complex *jones = malloc(4*NUM_DIRS*NUM_TIMES*NUM_FREQS*NUM_STATIONS*sizeof(double _Complex));

  do_run_lofar_beam(ms_path, apply_beam_norms,
                    rotate, element_only, jones);

  check_beam_values(jones, lofar_hba_jones_rotate);

  free(jones);
}

void test_run_hba_telescope_normed(void) {

  const char ms_path[] = "../../../../test_installation/everybeam/LOFAR_HBA_MOCK.ms";
  bool apply_beam_norms = true;
  bool rotate = false;
  bool element_only = false;

  double _Complex *jones = malloc(4*NUM_DIRS*NUM_TIMES*NUM_FREQS*NUM_STATIONS*sizeof(double _Complex));

  do_run_lofar_beam(ms_path, apply_beam_norms,
                    rotate, element_only, jones);

  check_beam_values(jones, lofar_hba_jones_normed);

  free(jones);
}

void test_run_hba_telescope_rotate_normed(void) {

  const char ms_path[] = "../../../../test_installation/everybeam/LOFAR_HBA_MOCK.ms";
  bool apply_beam_norms = true;
  bool rotate = true;
  bool element_only = false;

  double _Complex *jones = malloc(4*NUM_DIRS*NUM_TIMES*NUM_FREQS*NUM_STATIONS*sizeof(double _Complex));

  do_run_lofar_beam(ms_path, apply_beam_norms,
                    rotate, element_only, jones);

  check_beam_values(jones, lofar_hba_jones_rotate_normed);

  free(jones);
}

void test_run_hba_telescope_element_only(void) {

  const char ms_path[] = "../../../../test_installation/everybeam/LOFAR_HBA_MOCK.ms";
  bool apply_beam_norms = false;
  bool rotate = false;
  bool element_only = true;

  double _Complex *jones = malloc(4*NUM_DIRS*NUM_TIMES*NUM_FREQS*NUM_STATIONS*sizeof(double _Complex));

  do_run_lofar_beam(ms_path, apply_beam_norms,
                    rotate, element_only, jones);

  // check_beam_values(jones, lofar_hba_element);

  free(jones);
}



//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_run_hba_telescope);
    RUN_TEST(test_run_hba_telescope_rotate);
    RUN_TEST(test_run_hba_telescope_normed);
    RUN_TEST(test_run_hba_telescope_rotate_normed);
    // RUN_TEST(test_run_hba_telescope_element_only);

    return UNITY_END();
}
