#include <unity.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "call_everybeam_c.h"

// extern void learn_cpp(double _Complex *jones, int num_coords);

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

#ifdef DOUBLE_PRECISION
  double TOL = 1e-15;
#else
  double TOL = 1e-7;
#endif

// TEST_ASSERT_EQUAL_INT(NO_BEAM, beam_settings->beamtype );


void test_check_ms_telescope_type(const char *ms_path,
                                  const char *expected_telescope) {

  char *telescope = check_ms_telescope_type(ms_path);

  TEST_ASSERT_EQUAL_STRING(expected_telescope, telescope);

}


void test_check_lofar_telescope(void) {

  const char ms_path[] = "../../../../test_installation/everybeam/LOFAR_HBA_MOCK.ms";
  const char expected_telescope[] = "LOFAR";

  test_check_ms_telescope_type(ms_path, expected_telescope);

}

void test_check_mwa_telescope(void) {

  const char ms_path[] = "../../../../test_installation/everybeam/MWA-single-timeslot.ms";
  const char expected_telescope[] = "MWA";

  test_check_ms_telescope_type(ms_path, expected_telescope);

}

void test_check_oskar_telescope(void) {

  const char ms_path[] = "../../../../test_installation/everybeam/create_OSKAR-MWA_ms/OSKAR-MWA-layout.ms";
  const char expected_telescope[] = "OSKAR";

  test_check_ms_telescope_type(ms_path, expected_telescope);
}



//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_check_lofar_telescope);
    RUN_TEST(test_check_mwa_telescope);
    RUN_TEST(test_check_oskar_telescope);

    return UNITY_END();
}
