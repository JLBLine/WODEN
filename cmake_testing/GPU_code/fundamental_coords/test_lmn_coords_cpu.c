#include <unity.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "woden_precision_defs.h"
#include "lmn_coords_common.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

void test_calc_lmn_GivesCorrectlCoords_cpu(void){
  test_calc_lmn_GivesCorrectlCoords(0);
}

void test_calc_lmn_GivesCorrectmCoords_cpu(void){
  test_calc_lmn_GivesCorrectmCoords(0);
}

int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_calc_lmn_GivesCorrectlCoords_cpu);
    RUN_TEST(test_calc_lmn_GivesCorrectmCoords_cpu);

    return UNITY_END();
}
