#include "run_hyperbeam_interp_common.h"

// #define ROTATION 0

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

void test_delays1_freqs1_interpbeam_gpu(void){
  test_delays1_freqs1_interpbeam(1);
}

void test_delays2_freqs2_interpbeam_gpu(void){
  test_delays2_freqs2_interpbeam(1);
}

void test_delays3_freqs3_interpbeam_gpu(void){
  test_delays3_freqs3_interpbeam(1);
}

void test_delays4_freqs4_interpbeam_gpu(void){
  test_delays4_freqs4_interpbeam(1);
}

//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_delays1_freqs1_interpbeam_gpu);
    RUN_TEST(test_delays2_freqs2_interpbeam_gpu);
    RUN_TEST(test_delays3_freqs3_interpbeam_gpu);
    RUN_TEST(test_delays4_freqs4_interpbeam_gpu);

    return UNITY_END();
}
