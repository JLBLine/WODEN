#include "calc_measurement_equation_common.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

void test_calc_measurement_equation_gpu(void) {
  test_calc_measurement_equation(1);
}

void test_calc_measurement_equation_GiveCorrectValues_gpu(void) {
  test_calc_measurement_equation_GiveCorrectValues(1);
}

//Run the test with unity
int main(void)
{
    UNITY_BEGIN();
    RUN_TEST(test_calc_measurement_equation_gpu);
    RUN_TEST(test_calc_measurement_equation_GiveCorrectValues_gpu);

    return UNITY_END();
}
