#include <unity.h>
#include <stdlib.h>
#include <math.h>

// #include <module_a.h>
#include "constants.h"
// #include "fundamental_coords.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

//External CUDA code we're linking in
extern void test_kern_calc_lmn(float ra0, float dec0,
                                  float *ras, float *decs, int num_coords,
                                  float * ls, float * ms, float * ns);

#define UNITY_INCLUDE_FLOAT

// void setup_lmn_params()

void test_kern_calc_lmn_GivesCorrectlCoords(void)
{
    float ra0 = 0.0*DD2R;
    float dec0 = 0.0*DD2R;

    int num_points = 5;
    float *decs = malloc(num_points*sizeof(float));
    float *zeroes = malloc(num_points*sizeof(float));

    float ras[5] = {0.0, M_PI/6, M_PI/4, M_PI/3, M_PI/2};

    for (int i = 0; i < num_points; i++) {
      decs[i] = dec0;
      zeroes[i] = 0.0;
    }

    float *ls = malloc(num_points*sizeof(float));
    float *ms = malloc(num_points*sizeof(float));
    float *ns = malloc(num_points*sizeof(float));

    test_kern_calc_lmn(ra0, dec0, ras, decs, num_points,
                                   ls, ms, ns);

    float l_expected[5] = {0.0, 0.5, sqrt(2)/2.0, sqrt(3)/2.0, 1.0};
    float n_expected[5] = {1.0, sqrt(3)/2.0, sqrt(2)/2.0, 0.5, 0.0};

    TEST_ASSERT_EQUAL_FLOAT_ARRAY(l_expected, ls, num_points);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(zeroes, ms, num_points);
    // TEST_ASSERT_EQUAL_FLOAT_ARRAY(n_expected, ns, num_points);
    //Here, when n should be 0.0, we get some floating point error from the
    //sin/cos functions
    for (int i = 0; i < num_points; i++) {
      TEST_ASSERT_FLOAT_WITHIN(1e-7, n_expected[i], ns[i]);
    }

    free(decs);
    free(zeroes);
    free(ls);
    free(ms);
    free(ns);

}

int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_kern_calc_lmn_GivesCorrectlCoords);

    return UNITY_END();
}
