#include "apply_beam_gains_common.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

#define UNITY_INCLUDE_FLOAT

/*
Test that the code applying beam gains to Stokes parameter fluxes works
*/
void test_apply_beam_gains_GiveCorrectValues_gpu(void) {
  test_apply_beam_gains_GiveCorrectValues(1);

}

//Run the test with unity
int main(void)
{
    UNITY_BEGIN();
    RUN_TEST(test_apply_beam_gains_GiveCorrectValues_gpu);
    return UNITY_END();
}
