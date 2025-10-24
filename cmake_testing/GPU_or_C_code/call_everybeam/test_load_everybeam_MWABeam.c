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


void test_load_mwa_beam(void) {

  // int eb_status = 0;

  char* coeff_path = getenv("MWA_FEE_HDF5");

  double delays[16] = {0.0,0.0,0.0,0.0,
                    0.0,0.0,0.0,0.0,
                    0.0,0.0,0.0,0.0,
                    0.0,0.0,0.0,0.0};

  double amps[16] = {1.0,1.0,1.0,1.0,
                     1.0,1.0,1.0,1.0,
                     1.0,1.0,1.0,1.0,
                     1.0,1.0,1.0,1.0};

  Beam2016Implementation *eb_mwa = load_everybeam_MWABeam(coeff_path, delays, amps);
  
  destroy_everybeam_MWABeam(eb_mwa);
}





//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_load_mwa_beam);

    return UNITY_END();
}
