#include "run_hyperbeam_multi_antennas_common.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */


/*
Run the test but vary the frequency and pointings. Compare to pre-calculated
values that are stored in test_RTS_FEE_beam.h
*/
void test_hyperbeam_100MHz_zenith_multiant_gpu(void) {
  test_hyperbeam_100MHz_zenith_multiant(1);
}

void test_hyperbeam_150MHz_zenith_multiant_gpu(void) {
  test_hyperbeam_150MHz_zenith_multiant(1);
}

void test_hyperbeam_200MHz_zenith_multiant_gpu(void) {
  test_hyperbeam_200MHz_zenith_multiant(1);
}

void test_hyperbeam_100MHz_off_zenith1_multiant_gpu(void) {
  test_hyperbeam_100MHz_off_zenith1_multiant(1);
}

void test_hyperbeam_150MHz_off_zenith1_multiant_gpu(void) {
  test_hyperbeam_150MHz_off_zenith1_multiant(1);
}

void test_hyperbeam_200MHz_off_zenith1_multiant_gpu(void) {
  test_hyperbeam_200MHz_off_zenith1_multiant(1);
}

void test_hyperbeam_100MHz_off_zenith2_multiant_gpu(void) {
  test_hyperbeam_100MHz_off_zenith2_multiant(1);
}

void test_hyperbeam_150MHz_off_zenith2_multiant_gpu(void) {
  test_hyperbeam_150MHz_off_zenith2_multiant(1);
}

void test_hyperbeam_200MHz_off_zenith2_multiant_gpu(void) {
  test_hyperbeam_200MHz_off_zenith2_multiant(1);
}


//Run the test with unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_hyperbeam_100MHz_zenith_multiant_gpu);
    RUN_TEST(test_hyperbeam_150MHz_zenith_multiant_gpu);
    RUN_TEST(test_hyperbeam_200MHz_zenith_multiant_gpu);
    
    RUN_TEST(test_hyperbeam_100MHz_off_zenith1_multiant_gpu);
    RUN_TEST(test_hyperbeam_150MHz_off_zenith1_multiant_gpu);
    RUN_TEST(test_hyperbeam_200MHz_off_zenith1_multiant_gpu);
    
    RUN_TEST(test_hyperbeam_100MHz_off_zenith2_multiant_gpu);
    RUN_TEST(test_hyperbeam_150MHz_off_zenith2_multiant_gpu);
    RUN_TEST(test_hyperbeam_200MHz_off_zenith2_multiant_gpu);

    return UNITY_END();
}
