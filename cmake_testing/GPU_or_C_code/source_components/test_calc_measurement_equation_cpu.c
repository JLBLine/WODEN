#include "calc_measurement_equation_common.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

void test_calc_measurement_equation_cpu(void) {
  test_calc_measurement_equation(0);
}

void test_calc_measurement_equation_GiveCorrectValues_cpu(void) {
  test_calc_measurement_equation_GiveCorrectValues(0);
}

//Run the test with unity
int main(void)
{
    UNITY_BEGIN();
    RUN_TEST(test_calc_measurement_equation_cpu);
    RUN_TEST(test_calc_measurement_equation_GiveCorrectValues_cpu);

    return UNITY_END();
}
