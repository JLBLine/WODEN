#include "constants.h"

#define LST0 0.0
#define LST1 M_PI / 2
#define LST2 M_PI
#define LST3 230.0*DD2R

//First POINT SOURCE information
float point0_ras[] = {0.0, DD2R, 2*DD2R};
float point0_decs[] = {MWA_LAT_RAD, MWA_LAT_RAD, MWA_LAT_RAD};
float point0_ref_freqs[] = {100e+6, 150e+6, 200e+6};
float point0_ref_stokesI[] = {1.0, 2.0, 3.0};
float point0_ref_stokesQ[] = {0.0, 0.0, 0.0};
float point0_ref_stokesU[] = {0.0, 0.0, 0.0};
float point0_ref_stokesV[] = {0.0, 0.0, 0.0};
float point0_SIs[] = {-0.8, -0.2, 0.3};

//Second POINT SOURCE information
float point1_ras[] = {106*DD2R, 108*DD2R, 110*DD2R};
float point1_decs[] = {-30*DD2R, -30*DD2R, -30*DD2R};
float point1_ref_freqs[] = {125e+6, 175e+6, 2225e+6};
float point1_ref_stokesI[] = {5.0, 6.0, 7.0};
float point1_ref_stokesQ[] = {0.0, 0.0, 0.0};
float point1_ref_stokesU[] = {0.0, 0.0, 0.0};
float point1_ref_stokesV[] = {0.0, 0.0, 0.0};
float point1_SIs[] = {1.0, -0.4, -1.1};

//First GAUSSIAN SOURCE information
//Doesn't matter if this is same as POINT, there are other
//unique fields that ensure wrong things aren't being tested
float gauss0_ras[] = {0.0, DD2R, 2*DD2R};
float gauss0_decs[] = {MWA_LAT_RAD, MWA_LAT_RAD, MWA_LAT_RAD};
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
float gauss1_ras[] = {106*DD2R, 108*DD2R, 110*DD2R};
float gauss1_decs[] = {-30*DD2R, -30*DD2R, -30*DD2R};
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
float shape0_ras[] = {0.0, DD2R, 2*DD2R};
float shape0_decs[] = {MWA_LAT_RAD, MWA_LAT_RAD, MWA_LAT_RAD};
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
float shape1_ras[] = {106*DD2R, 108*DD2R, 110*DD2R};
float shape1_decs[] = {-30*DD2R, -30*DD2R, -30*DD2R};
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
