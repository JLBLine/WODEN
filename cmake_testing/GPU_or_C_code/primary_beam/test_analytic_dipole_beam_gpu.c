#include <unity.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "constants.h"
#include "woden_precision_defs.h"
#include "analytic_dipole_beam_common.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

/*
Test that the analytic dipole beam code returns the correct values, in the
correct order, for two time and frequency steps, with 25 directions on the sky
*/
void test_analytic_dipole_beam_GivesCorrectValues_gpu(void) {
  test_analytic_dipole_beam_GivesCorrectValues(1);
}

//Run the test with unity
int main(void)
{
    UNITY_BEGIN();
    RUN_TEST(test_analytic_dipole_beam_GivesCorrectValues_gpu);
    return UNITY_END();
}
