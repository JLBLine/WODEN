#include "extrap_stokes_common.h"
// #include "test_extrap_stokes.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

void test_extrap_stokes_GivesCorrectValues_gpu(void) {
  test_extrap_stokes_GivesCorrectValues(1);
}

//Run the test with unity
int main(void)
{
    UNITY_BEGIN();
    RUN_TEST(test_extrap_stokes_GivesCorrectValues_gpu);
    return UNITY_END();
}
