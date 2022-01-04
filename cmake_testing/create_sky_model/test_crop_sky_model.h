#include "constants.h"

#define LST0 0.0
#define LST1 M_PI / 2
#define LST2 M_PI
#define LST3 230.0*DD2R

//First POINT SOURCE information
double point0_ras[] = {0.0, DD2R, 2*DD2R};
double point0_decs[] = {MWA_LAT_RAD, MWA_LAT_RAD, MWA_LAT_RAD};
double point0_ref_freqs[] = {100e+6, 150e+6, 200e+6};
user_precision_t point0_ref_stokesI[] = {1.0, 2.0, 3.0};
user_precision_t point0_ref_stokesQ[] = {0.0, 0.0, 0.0};
user_precision_t point0_ref_stokesU[] = {0.0, 0.0, 0.0};
user_precision_t point0_ref_stokesV[] = {0.0, 0.0, 0.0};
user_precision_t point0_SIs[] = {-0.8, -0.2, 0.3};

//Second POINT SOURCE information
double point1_ras[] = {106*DD2R, 108*DD2R, 110*DD2R};
double point1_decs[] = {-30*DD2R, -30*DD2R, -30*DD2R};
double point1_ref_freqs[] = {125e+6, 175e+6, 2225e+6};
user_precision_t point1_ref_stokesI[] = {5.0, 6.0, 7.0};
user_precision_t point1_ref_stokesQ[] = {0.0, 0.0, 0.0};
user_precision_t point1_ref_stokesU[] = {0.0, 0.0, 0.0};
user_precision_t point1_ref_stokesV[] = {0.0, 0.0, 0.0};
user_precision_t point1_SIs[] = {1.0, -0.4, -1.1};

//First GAUSSIAN SOURCE information
//Doesn't matter if this is same as POINT, there are other
//unique fields that ensure wrong things aren't being tested
double gauss0_ras[] = {0.0, DD2R, 2*DD2R};
double gauss0_decs[] = {MWA_LAT_RAD, MWA_LAT_RAD, MWA_LAT_RAD};
double gauss0_ref_freqs[] = {100e+6, 150e+6, 200e+6};
user_precision_t gauss0_ref_stokesI[] = {1.0, 2.0, 3.0};
user_precision_t gauss0_ref_stokesQ[] = {0.0, 0.0, 0.0};
user_precision_t gauss0_ref_stokesU[] = {0.0, 0.0, 0.0};
user_precision_t gauss0_ref_stokesV[] = {0.0, 0.0, 0.0};
user_precision_t gauss0_SIs[] = {-0.8, -0.2, 0.3};
user_precision_t gauss0_majors[] = {6.0, 9.0, 12.0};
user_precision_t gauss0_minors[] = {1.0, 2.0, 3.0};
user_precision_t gauss0_pas[] = {0.0, 45.0, 90.0};

//Second GAUSSIAN SOURCE information
double gauss1_ras[] = {106*DD2R, 108*DD2R, 110*DD2R};
double gauss1_decs[] = {-30*DD2R, -30*DD2R, -30*DD2R};
double gauss1_ref_freqs[] = {125e+6, 175e+6, 2225e+6};
user_precision_t gauss1_ref_stokesI[] = {5.0, 6.0, 7.0};
user_precision_t gauss1_ref_stokesQ[] = {0.0, 0.0, 0.0};
user_precision_t gauss1_ref_stokesU[] = {0.0, 0.0, 0.0};
user_precision_t gauss1_ref_stokesV[] = {0.0, 0.0, 0.0};
user_precision_t gauss1_SIs[] = {1.0, -0.4, -1.1};
user_precision_t gauss1_majors[] = {2.0, 4.0, 12.0};
user_precision_t gauss1_minors[] = {3.0, 2.0, 1.0};
user_precision_t gauss1_pas[] = {180.0, 225.0, 270.0};

//First SHAPELET SOURCE information
//Doesn't matter if this is same as POINT, there are other
//unique fields that ensure wrong things aren't being tested
double shape0_ras[] = {0.0, DD2R, 2*DD2R};
double shape0_decs[] = {MWA_LAT_RAD, MWA_LAT_RAD, MWA_LAT_RAD};
double shape0_ref_freqs[] = {100e+6, 150e+6, 200e+6};
user_precision_t shape0_ref_stokesI[] = {1.0, 2.0, 3.0};
user_precision_t shape0_ref_stokesQ[] = {0.0, 0.0, 0.0};
user_precision_t shape0_ref_stokesU[] = {0.0, 0.0, 0.0};
user_precision_t shape0_ref_stokesV[] = {0.0, 0.0, 0.0};
user_precision_t shape0_SIs[] = {-0.8, -0.2, 0.3};
user_precision_t shape0_majors[] = {6.0, 9.0, 12.0};
user_precision_t shape0_minors[] = {1.0, 2.0, 3.0};
user_precision_t shape0_pas[] = {0.0, 45.0, 90.0};

//Three shapelet COMPONENTs,
//First has 2 basis functions
//Second has 3 basis functions
//Third has 5
user_precision_t shape0_coeffs[] = {0.25, 0.1, 0.98, 0.66, 0.8, 0.62, 0.67, 0.36, 0.38, 0.63};
user_precision_t shape0_n1s[] = {23, 69, 94, 28, 95, 23, 89, 82, 30, 93};
user_precision_t shape0_n2s[] = {97, 91, 23, 14, 93, 18, 32, 14, 97, 11};
user_precision_t shape0_param_indexes[] = {0, 0, 1, 1, 1, 2, 2, 2, 2, 2};

//Second SHAPELET SOURCE information
double shape1_ras[] = {106*DD2R, 108*DD2R, 110*DD2R};
double shape1_decs[] = {-30*DD2R, -30*DD2R, -30*DD2R};
double shape1_ref_freqs[] = {125e+6, 175e+6, 2225e+6};
user_precision_t shape1_ref_stokesI[] = {5.0, 6.0, 7.0};
user_precision_t shape1_ref_stokesQ[] = {0.0, 0.0, 0.0};
user_precision_t shape1_ref_stokesU[] = {0.0, 0.0, 0.0};
user_precision_t shape1_ref_stokesV[] = {0.0, 0.0, 0.0};
user_precision_t shape1_SIs[] = {1.0, -0.4, -1.1};
user_precision_t shape1_majors[] = {2.0, 4.0, 12.0};
user_precision_t shape1_minors[] = {3.0, 2.0, 1.0};
user_precision_t shape1_pas[] = {180.0, 225.0, 270.0};

//Three shapelet COMPONENTs,
//First has 4 basis functions
//Second has 1 basis functions
//Third has 2
user_precision_t shape1_coeffs[] = {0.33, 0.63, 0.01, 0.95, 0.33, 0.2, 0.74};
user_precision_t shape1_n1s[] = {1, 21, 80,  5, 92,  5, 98};
user_precision_t shape1_n2s[] = {16, 4, 27, 94, 73, 55, 56};
user_precision_t shape1_param_indexes[] = {0, 0, 0, 0, 1, 2, 2};

//Expected az/za for all time steps for the different LSTs
double expec_az_LST0[] = {0.0000000000000, 4.7114059180153, 4.7104228481364,
                          1.5747179058039, 1.5737347986641, 1.5727517139711,
                          1.5786399615312, 1.5776566901038, 1.5766734710889,
                          2.1579308795782, 2.1559930030733, 2.1540620700128 };
double expec_za_LST0[] = {0.0000000000000, 0.0039086342333, 0.0078172646892,
                          0.0155917779073, 0.0116831661944, 0.0077745431904,
                          0.0311833160128, 0.0272747869414, 0.0233662315047,
                          1.5593602379436, 1.5561040790870, 1.5528437438724 };

double expec_az_LST1[] = {4.2900584634321, 4.2884202839903, 4.2867773246606,
                          4.2965462595342, 4.2949268733767, 4.2933028168844,
                          4.3029601663223, 4.3013591423068, 4.2997535551609,
                          1.8659737209395, 1.8684735052647, 1.8710817825899,
                          1.8494119996711, 1.8511930524957, 1.8530518584716,
                          1.8376594596641, 1.8389159857357, 1.8402301007951 };
double expec_za_LST1[] = {1.3674640552467, 1.3710279492557, 1.3745892053280,
                          1.3532215231587, 1.3567957699779, 1.3603674468138,
                          1.3389383598151, 1.3425226930563, 1.3461045223980,
                          0.2521635568594, 0.2484253832967, 0.2446901374855,
                          0.2820763899820, 0.2783194359053, 0.2745644502349,
                          0.3121100672983, 0.3083404306473, 0.3045721247998 };

double expec_az_LST2[] = {4.3252592918701, 4.3240097928876, 4.3227558775006,
                          4.3350706082457, 4.3338561058628, 4.3326372377189,
                          4.3446042801037, 4.3434244460185, 4.3422402797543 };
double expec_za_LST2[] = {1.1174960228303, 1.1211144827620, 1.1247310899676,
                          1.0885625249286, 1.0921953358046, 1.0958263892164,
                          1.0595178840782, 1.0631643040794, 1.0668090578716 };
