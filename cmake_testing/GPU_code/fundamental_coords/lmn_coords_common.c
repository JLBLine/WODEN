#include "lmn_coords_common.h"

// //Tolerance used for testing
#define TOL 1e-15

void test_calc_lmn_GivesCorrectlCoords(int do_gpu){
    double ra0 = 0.0*DD2R;
    double dec0 = 0.0*DD2R;

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

    woden_settings_t *woden_settings = malloc(sizeof(woden_settings_t));
    woden_settings->ra0 = ra0;
    woden_settings->dec0 = dec0;
    woden_settings->sdec0 = sin(dec0);
    woden_settings->cdec0 = cos(dec0);

    if (do_gpu == 1){
        // test_kern_calc_lmn(ra0, dec0, ras, decs, num_points,
        //                                ls, ms, ns);
        test_calc_lmn_for_components_gpu(ls, ms, ns, ras, decs,
                                         num_points, woden_settings);
    } else {
        calc_lmn_cpu(ra0, sin(dec0), cos(dec0), ras, decs, ls, ms, ns, num_points);
    }

    double l_expected[9] = {-1.0, -sqrt(3)/2.0, -sqrt(2)/2.0, -0.5,
                            0.0, 0.5, sqrt(2)/2.0, sqrt(3)/2.0, 1.0};
    double n_expected[9] = {0.0, 0.5, sqrt(2)/2.0, sqrt(3)/2.0,
                           1.0, sqrt(3)/2.0, sqrt(2)/2.0, 0.5, 0.0};

    // Here, when n should be 0.0, we get some user_precision_ting point error from the
    // sin/cos functions
    for (int i = 0; i < num_points; i++) {
      TEST_ASSERT_DOUBLE_WITHIN(TOL, l_expected[i], ls[i]);
      TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, ms[i]);
      TEST_ASSERT_DOUBLE_WITHIN(TOL, n_expected[i], ns[i]);
    }    // Here, when n should be 0.0, we get some user_precision_ting point error from the
        // sin/cos functions

    free(decs);
    free(zeroes);
    free(ls);
    free(ms);
    free(ns);

}

void test_calc_lmn_GivesCorrectmCoords(int do_gpu){

    double ra0 = 0.0*DD2R;
    double dec0 = 0.0*DD2R;

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

    woden_settings_t *woden_settings = malloc(sizeof(woden_settings_t));
    woden_settings->ra0 = ra0;
    woden_settings->dec0 = dec0;
    woden_settings->sdec0 = sin(dec0);
    woden_settings->cdec0 = cos(dec0);

    if (do_gpu == 1){
        // test_kern_calc_lmn(ra0, dec0, ras, decs, num_points,
        //                                ls, ms, ns);
        test_calc_lmn_for_components_gpu(ls, ms, ns, ras, decs,
                                         num_points, woden_settings);
    } else {
        calc_lmn_cpu(ra0, sin(dec0), cos(dec0), ras, decs, ls, ms, ns, num_points);
    }

    double m_expected[9] = {-1.0, -sqrt(3)/2.0, -sqrt(2)/2.0, -0.5,
                            0.0, 0.5, sqrt(2)/2.0, sqrt(3)/2.0, 1.0};
    double n_expected[9] = {0.0, 0.5, sqrt(2)/2.0, sqrt(3)/2.0,
                           1.0, sqrt(3)/2.0, sqrt(2)/2.0, 0.5, 0.0};

    // Here, when n should be 0.0, we get some user_precision_ting point error from the
    // sin/cos functions
    for (int i = 0; i < num_points; i++) {
      TEST_ASSERT_DOUBLE_WITHIN(TOL, 0.0, ls[i]);
      TEST_ASSERT_DOUBLE_WITHIN(TOL, m_expected[i], ms[i]);
      TEST_ASSERT_DOUBLE_WITHIN(TOL, n_expected[i], ns[i]);
    }

    free(ras);
    free(zeroes);
    free(ls);
    free(ms);
    free(ns);
}