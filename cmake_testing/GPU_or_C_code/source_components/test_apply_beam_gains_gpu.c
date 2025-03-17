#include "apply_beam_gains_common.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

#define UNITY_INCLUDE_FLOAT

/*
Test that the code applying beam gains to Stokes parameter fluxes works
Do it for on cardinal dipole case, full polarisation case
*/
void test_apply_beam_gains_GiveCorrectValues_on_cardinal_fullpol_gpu(void) {
  test_apply_beam_gains_GiveCorrectValues(1, 0, 1);
}

/*
Test that the code applying beam gains to Stokes parameter fluxes works
Do it for off cardinal dipole case, full polarisation case
*/
void test_apply_beam_gains_GiveCorrectValues_off_cardinal_fullpol_gpu(void) {
  test_apply_beam_gains_GiveCorrectValues(1, 1, 1);
}

/*
Test that the code applying beam gains to Stokes parameter fluxes works
Do it for on cardinal dipole case, just Stokes I case
*/
void test_apply_beam_gains_GiveCorrectValues_on_cardinal_stokesI_gpu(void) {
  test_apply_beam_gains_GiveCorrectValues(1, 0, 0);
}

/*
Test that the code applying beam gains to Stokes parameter fluxes works
Do it for off cardinal dipole case, just Stokes I case
*/
void test_apply_beam_gains_GiveCorrectValues_off_cardinal_stokesI_gpu(void) {
  test_apply_beam_gains_GiveCorrectValues(1, 1, 0);
}

//Run the test with unity
int main(void)
{
    UNITY_BEGIN();
    RUN_TEST(test_apply_beam_gains_GiveCorrectValues_on_cardinal_fullpol_gpu);
    RUN_TEST(test_apply_beam_gains_GiveCorrectValues_off_cardinal_fullpol_gpu);
    RUN_TEST(test_apply_beam_gains_GiveCorrectValues_on_cardinal_stokesI_gpu);
    RUN_TEST(test_apply_beam_gains_GiveCorrectValues_off_cardinal_stokesI_gpu);
    return UNITY_END();
}
