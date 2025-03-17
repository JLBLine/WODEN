#include <unity.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "MWA_analytic_common.h"

#define FREQ 150000000

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

void test_MWA_analytic_beam_nside101_cpu(void) {
  test_MWA_analytic_beam_nside101(0);

}

//Same test but for a single az,za coord - good for diagnosing bugs
void test_single_azza_cpu(void) {
  test_single_azza(0);
}

//Run the test with unity
int main(void)
{
    UNITY_BEGIN();
    RUN_TEST(test_MWA_analytic_beam_nside101_cpu);
    // RUN_TEST(test_single_azza_cpu);

    return UNITY_END();
}
