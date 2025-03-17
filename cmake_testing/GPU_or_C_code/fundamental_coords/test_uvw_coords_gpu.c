#include <unity.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "uvw_coords_common.h"

/*Unity needs these calls, stick them here and leave them alone */
void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

/*
Checking the function fundamental_coords.cu::kern_calc_uvw
Checks that the wavelength scaling of u,v,w is happening correctly. Set HA=0
to make checking easier
*/
void test_kern_calc_uvw_ScalesByWavelength(void){
    test_calc_uvw_ScalesByWavelength(1);
}

/*
Checking the function fundamental_coords.cu::kern_calc_uvw
Checks that u,v,w coords change with time as expected
Make checking easier by setting dec phase centre dec0=0.0
*/
void test_kern_calc_uvw_RotateWithTime(void){
     test_calc_uvw_RotateWithTime(1);
}

/*
Checking the function fundamental_coords.cu::kern_calc_uvw_shapelet
Checks that u,v,w coords change with time as expected
Make checking easier by setting dec phase centre dec0=0.0
Also checks that results are scaled by wavelength correctly
*/
void test_kern_calc_uvw_shapelet_RotateWithTime(void){
    test_calc_uvw_shapelet_RotateWithTime(1);
}

int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_kern_calc_uvw_ScalesByWavelength);
    RUN_TEST(test_kern_calc_uvw_RotateWithTime);
    RUN_TEST(test_kern_calc_uvw_shapelet_RotateWithTime);

    return UNITY_END();
}
