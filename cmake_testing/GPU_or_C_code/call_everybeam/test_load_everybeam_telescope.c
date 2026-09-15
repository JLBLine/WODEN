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


void test_load_telescope(const char *ms_path,
                          const char *element_response_model,
                          bool use_differential_beam,
                          bool use_channel_frequency,
                          const char *coeff_path) {

  int eb_status = 0;

  //Test that the correct arguments load up fine
  Telescope *telescope = load_everybeam_telescope(&eb_status, ms_path, element_response_model,
                                                  use_differential_beam, use_channel_frequency,
                                                coeff_path);
  TEST_ASSERT_EQUAL_INT(0, eb_status);
  destroy_everybeam_telescope(telescope);

  // Check it throws an error with incorrect model
  telescope = load_everybeam_telescope(&eb_status, ms_path, "jeeeezus",
                                        use_differential_beam, use_channel_frequency,
                                        coeff_path);
  TEST_ASSERT_EQUAL_INT(1, eb_status);
  destroy_everybeam_telescope(telescope);
}


void test_load_lofar_telescope(void) {

  const char ms_path[] = "../../../../test_installation/everybeam/LOFAR_HBA_MOCK.ms";
  bool use_differential_beam = false;
  bool use_channel_frequency = true;
  const char coeff_path[] = "";
  const char element_response_model[] = "hamaker";

  test_load_telescope(ms_path, element_response_model, use_differential_beam,
                      use_channel_frequency, coeff_path);

  // int num_coords = 3;
  // double _Complex *jones = malloc(4*num_coords * num_coords * sizeof(double _Complex));
  // run_lofar_beam(telescope, jones, num_coords);
  // // printf("telescope: %p\n", telescope);
  // printf("jones[0]: %.8f %.8f\n", creal(jones[0]), cimag(jones[0]));
  // printf("jones[1]: %.8f %.8f\n", creal(jones[1]), cimag(jones[1]));
  // free(jones);
}

void test_load_mwa_telescope(void) {

  const char ms_path[] = "../../../../test_installation/everybeam/MWA-single-timeslot.ms";
  bool use_differential_beam = false;
  bool use_channel_frequency = true;
  const char coeff_path[] = "";
  const char element_response_model[] = "MWA";

  int eb_status = 0;

  // Check it throws an error with MWA; we should be using a different function
  Telescope *telescope = load_everybeam_telescope(&eb_status, ms_path, element_response_model,
                                        use_differential_beam, use_channel_frequency,
                                        coeff_path);
  TEST_ASSERT_EQUAL_INT(1, eb_status);
  destroy_everybeam_telescope(telescope);
}

void test_load_oskar_telescope(void) {

  const char ms_path[] = "../../../../test_installation/everybeam/create_OSKAR-MWA_ms/OSKAR-MWA-layout.ms";
  bool use_differential_beam = false;
  bool use_channel_frequency = true;
  const char coeff_path[] = "";
  const char element_response_model[] = "SKALA40_WAVE";

  test_load_telescope(ms_path, element_response_model, use_differential_beam,
                      use_channel_frequency, coeff_path);
}



//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_load_lofar_telescope);
    RUN_TEST(test_load_mwa_telescope);
    RUN_TEST(test_load_oskar_telescope);

    return UNITY_END();
}
