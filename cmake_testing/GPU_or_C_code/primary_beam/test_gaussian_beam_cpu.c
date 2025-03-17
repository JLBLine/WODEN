#include <unity.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "constants.h"
#include "woden_precision_defs.h"
#include "gaussian_beam_common.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

void test_GaussBeam_GivesCorrectlValues_cpu(void) {
  test_GaussBeam_GivesCorrectlValues(0);
}

void test_GaussBeam_GivesCorrectmValues_cpu(void) {
  test_GaussBeam_GivesCorrectmValues(0);
}

void test_GaussBeam_GivesCorrectlValuesByFreq_cpu(void) {
  test_GaussBeam_GivesCorrectlValuesByFreq(0);
}

void test_GaussBeam_GivesCorrectmValuesByFreq_cpu(void) {
  test_GaussBeam_GivesCorrectmValuesByFreq(0);
}

int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_GaussBeam_GivesCorrectlValues_cpu);
    RUN_TEST(test_GaussBeam_GivesCorrectmValues_cpu);
    RUN_TEST(test_GaussBeam_GivesCorrectlValuesByFreq_cpu);
    RUN_TEST(test_GaussBeam_GivesCorrectmValuesByFreq_cpu);

    return UNITY_END();
}
