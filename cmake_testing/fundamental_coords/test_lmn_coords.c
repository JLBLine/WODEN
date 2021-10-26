#include <unity.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

//External CUDA code we're linking in
extern void test_kern_calc_lmn(float ra0, float dec0,
                               double *ras, double *decs, int num_coords,
                               double * ls, double * ms, double * ns);

#define UNITY_INCLUDE_FLOAT

void test_kern_calc_lmn_GivesCorrectlCoords(void)
{
    float ra0 = 0.0*DD2R;
    float dec0 = 0.0*DD2R;

    int num_points = 9;
    double *decs = malloc(num_points*sizeof(double));
    double *zeroes = malloc(num_points*sizeof(double));

    //Keep RA between 0 and 2*pi here but enter RAs that should return
    //negative l values
    double ras[9] = {(3*M_PI)/2, (5*M_PI)/3, (7*M_PI)/4, (11*M_PI)/6,
                     0.0, M_PI/6, M_PI/4, M_PI/3, M_PI/2};

    for (int i = 0; i < num_points; i++) {
      decs[i] = dec0;
      zeroes[i] = 0.0;
    }

    double *ls = malloc(num_points*sizeof(double));
    double *ms = malloc(num_points*sizeof(double));
    double *ns = malloc(num_points*sizeof(double));

    test_kern_calc_lmn(ra0, dec0, ras, decs, num_points,
                                   ls, ms, ns);

    float l_expected[9] = {-1.0, -sqrt(3)/2.0, -sqrt(2)/2.0, -0.5,
                            0.0, 0.5, sqrt(2)/2.0, sqrt(3)/2.0, 1.0};
    float n_expected[9] = {0.0, 0.5, sqrt(2)/2.0, sqrt(3)/2.0,
                           1.0, sqrt(3)/2.0, sqrt(2)/2.0, 0.5, 0.0};

    // Here, when n should be 0.0, we get some floating point error from the
    // sin/cos functions
    for (int i = 0; i < num_points; i++) {
      TEST_ASSERT_EQUAL_FLOAT(l_expected[i], (float)ls[i]);
      TEST_ASSERT_EQUAL_FLOAT(0.0, (float)ms[i]);
      TEST_ASSERT_FLOAT_WITHIN(2e-7, n_expected[i], ns[i]);
    }    // Here, when n should be 0.0, we get some floating point error from the
        // sin/cos functions

    free(decs);
    free(zeroes);
    free(ls);
    free(ms);
    free(ns);

}

void test_kern_calc_lmn_GivesCorrectmCoords(void)
{
    float ra0 = 0.0*DD2R;
    float dec0 = 0.0*DD2R;

    int num_points = 9;
    double *ras = malloc(num_points*sizeof(double));
    double *zeroes = malloc(num_points*sizeof(double));

    //Test some know input/output values that should vary from -1.0 to 1.0
    //in m for a constant ra
    double decs[9] = {-M_PI/2, -M_PI/3, -M_PI/4, -M_PI/6,
                     0.0, M_PI/6, M_PI/4, M_PI/3, M_PI/2};

    for (int i = 0; i < num_points; i++) {
      ras[i] = ra0;
      zeroes[i] = 0.0;
    }

    double *ls = malloc(num_points*sizeof(double));
    double *ms = malloc(num_points*sizeof(double));
    double *ns = malloc(num_points*sizeof(double));

    test_kern_calc_lmn(ra0, dec0, ras, decs, num_points,
                                   ls, ms, ns);

    float m_expected[9] = {-1.0, -sqrt(3)/2.0, -sqrt(2)/2.0, -0.5,
                            0.0, 0.5, sqrt(2)/2.0, sqrt(3)/2.0, 1.0};
    float n_expected[9] = {0.0, 0.5, sqrt(2)/2.0, sqrt(3)/2.0,
                           1.0, sqrt(3)/2.0, sqrt(2)/2.0, 0.5, 0.0};

    // Here, when n should be 0.0, we get some floating point error from the
    // sin/cos functions
    for (int i = 0; i < num_points; i++) {
      TEST_ASSERT_EQUAL_FLOAT(0.0, (float)ls[i]);
      TEST_ASSERT_EQUAL_FLOAT(m_expected[i], (float)ms[i]);
      TEST_ASSERT_FLOAT_WITHIN(2e-7, n_expected[i], ns[i]);
    }

    free(ras);
    free(zeroes);
    free(ls);
    free(ms);
    free(ns);

}

int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_kern_calc_lmn_GivesCorrectlCoords);
    RUN_TEST(test_kern_calc_lmn_GivesCorrectmCoords);

    return UNITY_END();
}
