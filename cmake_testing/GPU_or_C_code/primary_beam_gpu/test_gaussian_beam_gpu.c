#include <unity.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "constants.h"
#include "woden_precision_defs.h"
#include "gaussian_beam_common.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

void test_GaussBeam_GivesCorrectlValues_gpu(void) {
  test_GaussBeam_GivesCorrectlValues(1);
}

void test_GaussBeam_GivesCorrectmValues_gpu(void) {
  test_GaussBeam_GivesCorrectmValues(1);
}

void test_GaussBeam_GivesCorrectlValuesByFreq_gpu(void) {
  test_GaussBeam_GivesCorrectlValuesByFreq(1);
}

void test_GaussBeam_GivesCorrectmValuesByFreq_gpu(void) {
  test_GaussBeam_GivesCorrectmValuesByFreq(1);
}

int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_GaussBeam_GivesCorrectlValues_gpu);
    RUN_TEST(test_GaussBeam_GivesCorrectmValues_gpu);
    RUN_TEST(test_GaussBeam_GivesCorrectlValuesByFreq_gpu);
    RUN_TEST(test_GaussBeam_GivesCorrectmValuesByFreq_gpu);

    return UNITY_END();
}
