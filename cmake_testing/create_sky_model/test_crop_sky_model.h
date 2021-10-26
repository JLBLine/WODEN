#include "constants.h"

#define LST0 0.0
#define LST1 M_PI / 2
#define LST2 M_PI
#define LST3 230.0*DD2R

//First POINT SOURCE information
double point0_ras[] = {0.0, DD2R, 2*DD2R};
double point0_decs[] = {MWA_LAT_RAD, MWA_LAT_RAD, MWA_LAT_RAD};
float point0_ref_freqs[] = {100e+6, 150e+6, 200e+6};
float point0_ref_stokesI[] = {1.0, 2.0, 3.0};
float point0_ref_stokesQ[] = {0.0, 0.0, 0.0};
float point0_ref_stokesU[] = {0.0, 0.0, 0.0};
float point0_ref_stokesV[] = {0.0, 0.0, 0.0};
float point0_SIs[] = {-0.8, -0.2, 0.3};

//Second POINT SOURCE information
double point1_ras[] = {106*DD2R, 108*DD2R, 110*DD2R};
double point1_decs[] = {-30*DD2R, -30*DD2R, -30*DD2R};
float point1_ref_freqs[] = {125e+6, 175e+6, 2225e+6};
float point1_ref_stokesI[] = {5.0, 6.0, 7.0};
float point1_ref_stokesQ[] = {0.0, 0.0, 0.0};
float point1_ref_stokesU[] = {0.0, 0.0, 0.0};
float point1_ref_stokesV[] = {0.0, 0.0, 0.0};
float point1_SIs[] = {1.0, -0.4, -1.1};

//First GAUSSIAN SOURCE information
//Doesn't matter if this is same as POINT, there are other
//unique fields that ensure wrong things aren't being tested
double gauss0_ras[] = {0.0, DD2R, 2*DD2R};
double gauss0_decs[] = {MWA_LAT_RAD, MWA_LAT_RAD, MWA_LAT_RAD};
float gauss0_ref_freqs[] = {100e+6, 150e+6, 200e+6};
float gauss0_ref_stokesI[] = {1.0, 2.0, 3.0};
float gauss0_ref_stokesQ[] = {0.0, 0.0, 0.0};
float gauss0_ref_stokesU[] = {0.0, 0.0, 0.0};
float gauss0_ref_stokesV[] = {0.0, 0.0, 0.0};
float gauss0_SIs[] = {-0.8, -0.2, 0.3};
float gauss0_majors[] = {6.0, 9.0, 12.0};
float gauss0_minors[] = {1.0, 2.0, 3.0};
float gauss0_pas[] = {0.0, 45.0, 90.0};

//Second GAUSSIAN SOURCE information
double gauss1_ras[] = {106*DD2R, 108*DD2R, 110*DD2R};
double gauss1_decs[] = {-30*DD2R, -30*DD2R, -30*DD2R};
float gauss1_ref_freqs[] = {125e+6, 175e+6, 2225e+6};
float gauss1_ref_stokesI[] = {5.0, 6.0, 7.0};
float gauss1_ref_stokesQ[] = {0.0, 0.0, 0.0};
float gauss1_ref_stokesU[] = {0.0, 0.0, 0.0};
float gauss1_ref_stokesV[] = {0.0, 0.0, 0.0};
float gauss1_SIs[] = {1.0, -0.4, -1.1};
float gauss1_majors[] = {2.0, 4.0, 12.0};
float gauss1_minors[] = {3.0, 2.0, 1.0};
float gauss1_pas[] = {180.0, 225.0, 270.0};

//First SHAPELET SOURCE information
//Doesn't matter if this is same as POINT, there are other
//unique fields that ensure wrong things aren't being tested
double shape0_ras[] = {0.0, DD2R, 2*DD2R};
double shape0_decs[] = {MWA_LAT_RAD, MWA_LAT_RAD, MWA_LAT_RAD};
float shape0_ref_freqs[] = {100e+6, 150e+6, 200e+6};
float shape0_ref_stokesI[] = {1.0, 2.0, 3.0};
float shape0_ref_stokesQ[] = {0.0, 0.0, 0.0};
float shape0_ref_stokesU[] = {0.0, 0.0, 0.0};
float shape0_ref_stokesV[] = {0.0, 0.0, 0.0};
float shape0_SIs[] = {-0.8, -0.2, 0.3};
float shape0_majors[] = {6.0, 9.0, 12.0};
float shape0_minors[] = {1.0, 2.0, 3.0};
float shape0_pas[] = {0.0, 45.0, 90.0};

//Three shapelet COMPONENTs,
//First has 2 basis functions
//Second has 3 basis functions
//Third has 5
float shape0_coeffs[] = {0.25, 0.1, 0.98, 0.66, 0.8, 0.62, 0.67, 0.36, 0.38, 0.63};
float shape0_n1s[] = {23, 69, 94, 28, 95, 23, 89, 82, 30, 93};
float shape0_n2s[] = {97, 91, 23, 14, 93, 18, 32, 14, 97, 11};
float shape0_param_indexes[] = {0, 0, 1, 1, 1, 2, 2, 2, 2, 2};

//Second SHAPELET SOURCE information
double shape1_ras[] = {106*DD2R, 108*DD2R, 110*DD2R};
double shape1_decs[] = {-30*DD2R, -30*DD2R, -30*DD2R};
float shape1_ref_freqs[] = {125e+6, 175e+6, 2225e+6};
float shape1_ref_stokesI[] = {5.0, 6.0, 7.0};
float shape1_ref_stokesQ[] = {0.0, 0.0, 0.0};
float shape1_ref_stokesU[] = {0.0, 0.0, 0.0};
float shape1_ref_stokesV[] = {0.0, 0.0, 0.0};
float shape1_SIs[] = {1.0, -0.4, -1.1};
float shape1_majors[] = {2.0, 4.0, 12.0};
float shape1_minors[] = {3.0, 2.0, 1.0};
float shape1_pas[] = {180.0, 225.0, 270.0};

//Three shapelet COMPONENTs,
//First has 4 basis functions
//Second has 1 basis functions
//Third has 2
float shape1_coeffs[] = {0.33, 0.63, 0.01, 0.95, 0.33, 0.2, 0.74};
float shape1_n1s[] = {1, 21, 80,  5, 92,  5, 98};
float shape1_n2s[] = {16, 4, 27, 94, 73, 55, 56};
float shape1_param_indexes[] = {0, 0, 0, 0, 1, 2, 2};

//Expected az/za for all time steps for the different LSTs
float expec_az_LST0[] = {0.000000, 4.711406, 4.710422, 1.574718,
                         1.573735, 1.572752, 1.578640, 1.577657,
                         1.576674, 2.157931, 2.155993, 2.154062 };
float expec_za_LST0[] = {0.000000, 0.003909, 0.007817, 0.015592,
                         0.011683, 0.007775, 0.031183, 0.027275,
                         0.023366, 1.559360, 1.556104, 1.552844 };

float expec_az_LST1[] = {4.290059, 4.288420, 4.286777, 4.296546, 4.294927,
                         4.293303, 4.302960, 4.301359, 4.299754, 1.865974,
                         1.868474, 1.871082, 1.849412, 1.851193, 1.853052,
                         1.837659, 1.838916, 1.840230 };
float expec_za_LST1[] = {1.367464, 1.371028, 1.374589, 1.353222, 1.356796,
                         1.360368, 1.338938, 1.342523, 1.346105, 0.252164,
                         0.248425, 0.244690, 0.282076, 0.278320, 0.274564,
                         0.312110, 0.308340, 0.304572 };

float expec_az_LST2[] = {4.325259, 4.324010, 4.322756, 4.335071, 4.333856,
                         4.332637, 4.344604, 4.343424, 4.342240 };
float expec_za_LST2[] = { 1.117496, 1.121114, 1.124731, 1.088563, 1.092195,
                          1.095826, 1.059518, 1.063164, 1.066809 };
