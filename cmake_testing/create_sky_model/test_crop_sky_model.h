#include "constants.h"

#define LST0 0.0
#define LST1 M_PI / 2
#define LST2 M_PI
#define LST3 230.0*DD2R

//-------------------------------------------------------------

//Give all the component types the same sets of RA/DEC so we
//can predict their cropping behaviours the same way
double comp0_ras[] = {0.0, DD2R, 2*DD2R, 0.0, DD2R, 2*DD2R, 0.0, DD2R, 2*DD2R};
double comp0_decs[] = {MWA_LAT_RAD, MWA_LAT_RAD, MWA_LAT_RAD,
                        MWA_LAT_RAD, MWA_LAT_RAD, MWA_LAT_RAD,
                        MWA_LAT_RAD, MWA_LAT_RAD, MWA_LAT_RAD};

double comp1_ras[] = {106*DD2R, 108*DD2R, 110*DD2R,
                       106*DD2R, 108*DD2R, 110*DD2R,
                       106*DD2R, 108*DD2R, 110*DD2R};
double comp1_decs[] = {-30*DD2R, -30*DD2R, -30*DD2R,
                        -30*DD2R, -30*DD2R, -30*DD2R,
                        -30*DD2R, -30*DD2R, -30*DD2R};


/*******************************************************************************
POINT SOURCE THINGIES
*******************************************************************************/

//First POINT SOURCE information
double point0_ref_freqs[] = {100e+6, 150e+6, 200e+6};
user_precision_t point0_ref_stokesI[] = {1.0, 2.0, 3.0};
user_precision_t point0_ref_stokesQ[] = {0.0, 0.0, 0.0};
user_precision_t point0_ref_stokesU[] = {0.0, 0.0, 0.0};
user_precision_t point0_ref_stokesV[] = {0.0, 0.0, 0.0};
user_precision_t point0_SIs[] = {-0.8, -0.2, 0.3};
user_precision_t point0_curve_qs[] = {0.4, 0.01, -0.5};
int point0_power_inds[] = {0, 1, 2};
int point0_curve_inds[] = {3, 4, 5};

//First list component - 3 entries
//Second list component - 2 entries
//Third list component - 4 entries
int point0_list_comp_inds[] = {6, 7, 8};
int point0_num_list_values[] = {3, 2, 4};
int point0_list_start_indexes[] = {0, 3, 5};
double point0_list_freqs[] = {217e+6, 224e+6, 246e+6, 55e+6, 86e+6,
                              106e+6, 118e+6, 129e+6, 141e+6};
user_precision_t point0_list_stokesI[] = {62.106, 76.393, 79.125, 10.446, 45.208,
                                   42.554, 5.386, 79.501, 72.694};
user_precision_t point0_list_stokesQ[] = {34.049, 43.548, 12.166, 97.499, 63.733,
                                   4.192, 51.157, 82.933, 63.51};
user_precision_t point0_list_stokesU[] = {0.475, 0.676, 0.829, 0.781, 0.694, 0.696,
                                   0.363, 0.904, 0.796};
user_precision_t point0_list_stokesV[] = {3.95 , 0.064, 8.741, 4.046, 8.169, 7.317,
                                   4.627, 8.497, 3.784};

//Second POINT SOURCE information

double point1_ref_freqs[] = {125e+6, 175e+6, 2225e+6};
user_precision_t point1_ref_stokesI[] = {5.0, 6.0, 7.0};
user_precision_t point1_ref_stokesQ[] = {0.0, 0.0, 0.0};
user_precision_t point1_ref_stokesU[] = {0.0, 0.0, 0.0};
user_precision_t point1_ref_stokesV[] = {0.0, 0.0, 0.0};
user_precision_t point1_SIs[] = {1.0, -0.4, -1.1};
user_precision_t point1_curve_qs[] = {0.1, 1, -0.0005};
int point1_power_inds[] = {0, 1, 2};
int point1_curve_inds[] = {3, 4, 5};

//First list component - 4 entries
//Second list component - 2 entries
//Third list component - 3 entries
int point1_list_comp_inds[] = {6, 7, 8};
int point1_num_list_values[] = {4, 2, 3};
int point1_list_start_indexes[] = {0, 4, 6};
double point1_list_freqs[] = {150e+6, 192e+6, 197e+6, 206e+6, 60e+6,  87e+6,
                              95e+6, 107e+6, 114e+6};
user_precision_t point1_list_stokesI[] = {224.154, 251.909, 368.004, 494.363,
                                      248.958, 62.811, 417.319, 131.881, 89.317};
user_precision_t point1_list_stokesQ[] = {0.013, 0.789, 0.41, 0.589, 0.035, 0.215,
                                          0.864, 0.227, 0.505};
user_precision_t point1_list_stokesU[] = {8.895, 8.861, 2.21 , 5.065, 2.164, 5.27,
                                          0.505, 5.613, 4.365};
user_precision_t point1_list_stokesV[] = {3.714, 5.898, 4.345, 5.359, 2.056,
                                                    0.586, 5.466, 9.743, 7.149};

/*******************************************************************************
GAUSSIAN SOURCE THINGIES
*******************************************************************************/

//First GAUSSIAN SOURCE information
//Doesn't matter if this is same as POINT, there are other
//unique fields that ensure wrong things aren't being tested
user_precision_t gauss0_majors[] = {6.0, 9.0, 12.0, 6.0, 9.0, 12.0, 6.0, 9.0, 12.0};
user_precision_t gauss0_minors[] = {1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0};
user_precision_t gauss0_pas[] = {0.0, 45.0, 90.0, 0.0, 45.0, 90.0, 0.0, 45.0, 90.0};

double gauss0_ref_freqs[] = {100e+6, 150e+6, 200e+6};
user_precision_t gauss0_ref_stokesI[] = {1.0, 2.0, 3.0};
user_precision_t gauss0_ref_stokesQ[] = {0.0, 0.0, 0.0};
user_precision_t gauss0_ref_stokesU[] = {0.0, 0.0, 0.0};
user_precision_t gauss0_ref_stokesV[] = {0.0, 0.0, 0.0};
user_precision_t gauss0_SIs[] = {-0.8, -0.2, 0.3};
user_precision_t gauss0_curve_qs[] = {0.3, -0.345, 0.928};
int gauss0_power_inds[] = {0, 1, 2};
int gauss0_curve_inds[] = {3, 4, 5};

//First list component - 3 entries
//Second list component - 3 entries
//Third list component - 4 entries
int gauss0_list_comp_inds[] = {6, 7, 8};
int gauss0_num_list_values[] = {3, 3, 4};
int gauss0_list_start_indexes[] = {0, 3, 6};
double gauss0_list_freqs[] = {55e+6,  61e+6, 107e+6, 135e+6, 149e+6, 189e+6,
                              197e+6, 239e+6, 241e+6, 245e+6};
user_precision_t gauss0_list_stokesI[] = {449.553,  14.159, 237.118, 101.845,
                                          234.062, 329.211, 160.568, 154.85,
                                          463.724, 443.669};
user_precision_t gauss0_list_stokesQ[] = {0.455, 0.295, 0.867, 0.473, 0.521,
                                              0.96 , 0.779, 0.812, 0.15 ,0.519};
user_precision_t gauss0_list_stokesU[] = {0.85 , 0.039, 0.012, 0.15 , 0.084,
                                              0.534, 0.817, 0.807, 0.147,  0.97};
user_precision_t gauss0_list_stokesV[] = {0.52 , 0.405, 0.681, 0.212, 0.58 ,
                                              0.544, 0.786, 0.344, 0.385, 0.047};



//Second GAUSSIAN SOURCE information
user_precision_t gauss1_majors[] = {2.0, 4.0, 12.0, 2.0, 4.0, 12.0, 2.0, 4.0, 12.0};
user_precision_t gauss1_minors[] = {3.0, 2.0, 1.0, 3.0, 2.0, 1.0, 3.0, 2.0, 1.0};
user_precision_t gauss1_pas[] = {180.0, 225.0, 270.0, 180.0, 225.0, 270.0,180.0, 225.0, 270.0};


double gauss1_ref_freqs[] = {125e+6, 175e+6, 2225e+6};
user_precision_t gauss1_ref_stokesI[] = {5.0, 6.0, 7.0};
user_precision_t gauss1_ref_stokesQ[] = {0.0, 0.0, 0.0};
user_precision_t gauss1_ref_stokesU[] = {0.0, 0.0, 0.0};
user_precision_t gauss1_ref_stokesV[] = {0.0, 0.0, 0.0};
user_precision_t gauss1_SIs[] = {1.0, -0.4, -1.1};
user_precision_t gauss1_curve_qs[] = {-3.3, -0.1234, 0.2345};
int gauss1_power_inds[] = {0, 1, 2};
int gauss1_curve_inds[] = {3, 4, 5};

//First list component - 5 entries
//Second list component - 3 entries
//Third list component - 4 entries
int gauss1_list_comp_inds[] = {6, 7, 8};
int gauss1_num_list_values[] = {5, 3, 4};
int gauss1_list_start_indexes[] = {0, 5, 8};
double gauss1_list_freqs[] = {160e+6, 204e+6, 209e+6, 236e+6,  247e+6, 52e+6,
                              55e+6,  74e+6, 100e+6, 105e+6, 136e+6, 148e+6};
user_precision_t gauss1_list_stokesI[] = {41.46, 25.88, 172.69, 208.48, 44.78,
                        162.81, 145.8 , 283.68, 120.13, 120.93, 143.61,  26.62};
user_precision_t gauss1_list_stokesQ[] = {0.844, 0.465, 0.946, 0.1, 0.841,
                                0.344, 0.955, 0.931, 0.27 ,0.788, 0.653, 0.693};
user_precision_t gauss1_list_stokesU[] = {0.783, 0.551, 0.046, 0.929, 0.649,
                                0.034, 0.064, 0.887, 0.203, 0.337, 0.507, 0.061};
user_precision_t gauss1_list_stokesV[] = {0.701, 0.862, 0.885, 0.924, 0.403,
                                  0.414, 0.91, 0.253, 0.543, 0.151, 0.066, 0.28};


/*******************************************************************************
SHAPELET SOURCE THINGIES
*******************************************************************************/

//First SHAPELET SOURCE information---------------------------------------------
//Doesn't matter if this is same as POINT, there are other
//unique fields that ensure wrong things aren't being tested
double shape0_ref_freqs[] = {100e+6, 150e+6, 200e+6};
user_precision_t shape0_majors[] = {6.0, 9.0, 12.0, 6.0, 9.0, 12.0, 6.0, 9.0, 12.0};
user_precision_t shape0_minors[] = {1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0};
user_precision_t shape0_pas[] = {0.0, 45.0, 90.0, 0.0, 45.0, 90.0, 0.0, 45.0, 90.0};

//power / curved flux stuff
user_precision_t shape0_ref_stokesI[] = {1.0, 2.0, 3.0};
user_precision_t shape0_ref_stokesQ[] = {0.0, 0.0, 0.0};
user_precision_t shape0_ref_stokesU[] = {0.0, 0.0, 0.0};
user_precision_t shape0_ref_stokesV[] = {0.0, 0.0, 0.0};
user_precision_t shape0_SIs[] = {-0.8, -0.2, 0.3};
user_precision_t shape0_curve_qs[] = {0.23, 0.1, -0.432};
int shape0_power_inds[] = {0, 1, 2};
int shape0_curve_inds[] = {3, 4, 5};


//First list component - 2 entries
//Second list component - 4 entries
//Third list component - 2 entries
int shape0_list_comp_inds[] = {6, 7, 8};
int shape0_num_list_values[] = {2, 4, 2};
int shape0_list_start_indexes[] = {0, 2, 6};
double shape0_list_freqs[] = {214e+6, 249e+6, 130e+6, 174e+6, 180e+6, 182e+6,
                              64e+6,  67e+6};
user_precision_t shape0_list_stokesI[] = {120.74 , 157.087,  69.962,  81.622,
                                            247.704,  74.592, 184.463, 137.092};
user_precision_t shape0_list_stokesQ[] = {0.644, 0.819, 0.544, 0.402, 0.404,
                                                           0.927, 0.295, 0.184};
user_precision_t shape0_list_stokesU[] = {0.41 , 0.774, 0.103, 0.813, 0.453,
                                                           0.442, 0.125, 0.431};
user_precision_t shape0_list_stokesV[] = {0.995, 0.646, 0.174, 0.592, 0.09,
                                                           0.314, 0.919, 0.485};

//We have 9 shapelet COMPONENTs, 3 of each flux type
//Need to assign different basis functions to each component. These all go
//into one dimensional arrays
//Number of coeffs in each component go as 4, 3, 2, 2, 5, 5, 5, 1, 2
int NUM_SHAPE0_COEFFS = 29;
user_precision_t shape0_param_indexes[] = {0, 0, 0, 0, 1, 1, 1, 2, 2, 3, 3,
    4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 8, 8};

user_precision_t shape0_coeffs[] = {0.752, 0.96 , 0.971, 0.337, 0.525, 0.19,
       0.284, 0.37 , 0.269, 0.584, 0.41 , 0.055, 0.168, 0.633, 0.028, 0.231,
       0.948, 0.261, 0.924, 0.162, 0.655, 0.391, 0.627, 0.719, 0.549, 0.562,
       0.306, 0.548, 0.361};
user_precision_t shape0_n1s[] = {7, 94, 85, 28,  4, 12, 17, 79, 41, 49,  8, 22,
             77, 45, 46, 42, 56, 52, 92, 40, 95, 69, 95, 78, 81, 58, 88, 91, 16};
user_precision_t shape0_n2s[] = {8, 89, 50, 60, 25, 24, 36, 52, 53, 26, 71,
          6, 34, 52, 46, 88, 68, 70,  9, 11, 18,  7, 46, 89, 51, 33, 22, 65, 35};



//First SHAPELET SOURCE information---------------------------------------------
user_precision_t shape1_majors[] = {2.0, 4.0, 12.0, 2.0, 4.0, 12.0,2.0, 4.0, 12.0};
user_precision_t shape1_minors[] = {3.0, 2.0, 1.0, 3.0, 2.0, 1.0,3.0, 2.0, 1.0};
user_precision_t shape1_pas[] = {180.0, 225.0, 270.0, 180.0, 225.0, 270.0,
                                                           180.0, 225.0, 270.0};

double shape1_ref_freqs[] = {125e+6, 175e+6, 225e+6};
user_precision_t shape1_ref_stokesI[] = {5.0, 6.0, 7.0};
user_precision_t shape1_ref_stokesQ[] = {0.0, 0.0, 0.0};
user_precision_t shape1_ref_stokesU[] = {0.0, 0.0, 0.0};
user_precision_t shape1_ref_stokesV[] = {0.0, 0.0, 0.0};
user_precision_t shape1_SIs[] = {1.0, -0.4, -1.1};
user_precision_t shape1_curve_qs[] = {0.23, 0.1, -0.432};
int shape1_power_inds[] = {0, 1, 2};
int shape1_curve_inds[] = {3, 4, 5};


//First list component - 3 entries
//Second list component - 5 entries
//Third list component - 1 entries
int shape1_list_comp_inds[] = {6, 7, 8};
int shape1_num_list_values[] = {3, 5, 1};
int shape1_list_start_indexes[] = {0, 3, 8};
double shape1_list_freqs[] = {194e+6, 199e+6, 227e+6, 94e+6, 96e+6, 171e+6,
                              183e+6, 186e+6, 88e+6};
user_precision_t shape1_list_stokesI[] = {82.61 , 107.258, 174.16 , 142.765,
                                    83.215,  81.279,  92.312, 102.499,  71.699};
user_precision_t shape1_list_stokesQ[] = {0.845, 0.644, 0.657, 0.981, 0.215,
                                                    0.568, 0.319, 0.348, 0.808};
user_precision_t shape1_list_stokesU[] = {0.291, 0.19 , 0.696, 0.011, 0.406,
                                                    0.86 , 0.173, 0.696, 0.649};
user_precision_t shape1_list_stokesV[] = {0.756, 0.908, 0.611, 0.793, 0.936,
                                                    0.479, 0.839, 0.363, 0.917};

//We have 9 shapelet COMPONENTs, 3 of each flux type
//Need to assign different basis functions to each component. These all go
//into one dimensional arrays
//Number of coeffs in each component go as 3, 3, 1, 5, 2, 4, 1, 2, 3

int NUM_SHAPE1_COEFFS = 24;
user_precision_t shape1_param_indexes[] = {0, 0, 0, 1, 1, 1, 2, 3, 3, 3, 3, 3,
             4, 4, 5, 5, 5, 5, 6, 7, 7, 8, 8, 8};

user_precision_t shape1_coeffs[] = {0.609, 0.695, 0.561, 0.489, 0.538, 0.793,
  0.501, 0.628, 0.06 ,0.53 , 0.772, 0.253, 0.303, 0.942, 0.536, 0.752, 0.451,
                            0.229, 0.805, 0.973, 0.165, 0.978, 0.747, 0.469};
user_precision_t shape1_n1s[] = {18, 24, 81, 69, 75, 50, 71, 50, 88, 95, 63, 33,
                                  20, 75, 74, 18, 96, 4, 89, 82, 12, 59, 91, 93};
user_precision_t shape1_n2s[] = {51, 47, 10, 75,  8, 15, 10, 45,  5,  8, 90, 44,
                                  48, 70, 59, 76, 32, 90, 55,  1, 42, 48, 93, 18};

/*******************************************************************************
Expected azimuth and zenith angle things
*******************************************************************************/


double expec_az_sourcecrop_LST0[] = {0.0000000000000, 4.7114059180153, 4.7104228481364,
                          1.5747179058039, 1.5737347986641, 1.5727517139711,
                          1.5786399615312, 1.5776566901038, 1.5766734710889,
                          0.0000000000000, 4.7114059180153, 4.7104228481364,
                          1.5747179058039, 1.5737347986641, 1.5727517139711,
                          1.5786399615312, 1.5776566901038, 1.5766734710889,
                          0.0000000000000, 4.7114059180153, 4.7104228481364,
                          1.5747179058039, 1.5737347986641, 1.5727517139711,
                          1.5786399615312, 1.5776566901038, 1.5766734710889};

double expec_za_sourcecrop_LST0[] = {0.0000000000000, 0.0039086342333, 0.0078172646892,
                          0.0155917779073, 0.0116831661944, 0.0077745431904,
                          0.0311833160128, 0.0272747869414, 0.0233662315047,
                          0.0000000000000, 0.0039086342333, 0.0078172646892,
                          0.0155917779073, 0.0116831661944, 0.0077745431904,
                          0.0311833160128, 0.0272747869414, 0.0233662315047,
                          0.0000000000000, 0.0039086342333, 0.0078172646892,
                          0.0155917779073, 0.0116831661944, 0.0077745431904,
                          0.0311833160128, 0.0272747869414, 0.0233662315047};

double expec_az_compcrop_LST0[] = {0.0000000000000, 4.7114059180153, 4.7104228481364,
                          1.5747179058039, 1.5737347986641, 1.5727517139711,
                          1.5786399615312, 1.5776566901038, 1.5766734710889,
                          0.0000000000000, 4.7114059180153, 4.7104228481364,
                          1.5747179058039, 1.5737347986641, 1.5727517139711,
                          1.5786399615312, 1.5776566901038, 1.5766734710889,
                          0.0000000000000, 4.7114059180153, 4.7104228481364,
                          1.5747179058039, 1.5737347986641, 1.5727517139711,
                          1.5786399615312, 1.5776566901038, 1.5766734710889,
                          2.1579308795782, 2.1559930030733, 2.1540620700128,
                          2.1579308795782, 2.1559930030733, 2.1540620700128,
                          2.1579308795782, 2.1559930030733, 2.1540620700128,};

double expec_za_compcrop_LST0[] = {0.0000000000000, 0.0039086342333, 0.0078172646892,
                          0.0155917779073, 0.0116831661944, 0.0077745431904,
                          0.0311833160128, 0.0272747869414, 0.0233662315047,
                          0.0000000000000, 0.0039086342333, 0.0078172646892,
                          0.0155917779073, 0.0116831661944, 0.0077745431904,
                          0.0311833160128, 0.0272747869414, 0.0233662315047,
                          0.0000000000000, 0.0039086342333, 0.0078172646892,
                          0.0155917779073, 0.0116831661944, 0.0077745431904,
                          0.0311833160128, 0.0272747869414, 0.0233662315047,
                          1.5593602379436, 1.5561040790870, 1.5528437438724,
                          1.5593602379436, 1.5561040790870, 1.5528437438724,
                          1.5593602379436, 1.5561040790870, 1.5528437438724};

//These need to be split up and duplicated agains - wait until doing all
//three flux types though

double expec_az_LST1[] = {4.2900584634321, 4.2884202839903, 4.2867773246606,
                          4.2965462595342, 4.2949268733767, 4.2933028168844,
                          4.3029601663223, 4.3013591423068, 4.2997535551609,
                          4.2900584634321, 4.2884202839903, 4.2867773246606,
                          4.2965462595342, 4.2949268733767, 4.2933028168844,
                          4.3029601663223, 4.3013591423068, 4.2997535551609,
                          4.2900584634321, 4.2884202839903, 4.2867773246606,
                          4.2965462595342, 4.2949268733767, 4.2933028168844,
                          4.3029601663223, 4.3013591423068, 4.2997535551609,
                          1.8659737209395, 1.8684735052647, 1.8710817825899,
                          1.8494119996711, 1.8511930524957, 1.8530518584716,
                          1.8376594596641, 1.8389159857357, 1.8402301007951,
                          1.8659737209395, 1.8684735052647, 1.8710817825899,
                          1.8494119996711, 1.8511930524957, 1.8530518584716,
                          1.8376594596641, 1.8389159857357, 1.8402301007951,
                          1.8659737209395, 1.8684735052647, 1.8710817825899,
                          1.8494119996711, 1.8511930524957, 1.8530518584716,
                          1.8376594596641, 1.8389159857357, 1.8402301007951 };

double expec_za_LST1[] = {1.3674640552467, 1.3710279492557, 1.3745892053280,
                          1.3532215231587, 1.3567957699779, 1.3603674468138,
                          1.3389383598151, 1.3425226930563, 1.3461045223980,
                          1.3674640552467, 1.3710279492557, 1.3745892053280,
                          1.3532215231587, 1.3567957699779, 1.3603674468138,
                          1.3389383598151, 1.3425226930563, 1.3461045223980,
                          1.3674640552467, 1.3710279492557, 1.3745892053280,
                          1.3532215231587, 1.3567957699779, 1.3603674468138,
                          1.3389383598151, 1.3425226930563, 1.3461045223980,
                          0.2521635568594, 0.2484253832967, 0.2446901374855,
                          0.2820763899820, 0.2783194359053, 0.2745644502349,
                          0.3121100672983, 0.3083404306473, 0.3045721247998,
                          0.2521635568594, 0.2484253832967, 0.2446901374855,
                          0.2820763899820, 0.2783194359053, 0.2745644502349,
                          0.3121100672983, 0.3083404306473, 0.3045721247998,
                          0.2521635568594, 0.2484253832967, 0.2446901374855,
                          0.2820763899820, 0.2783194359053, 0.2745644502349,
                          0.3121100672983, 0.3083404306473, 0.3045721247998 };

double expec_az_LST2[] = {4.3252592918701, 4.3240097928876, 4.3227558775006,
                          4.3350706082457, 4.3338561058628, 4.3326372377189,
                          4.3446042801037, 4.3434244460185, 4.3422402797543,
                          4.3252592918701, 4.3240097928876, 4.3227558775006,
                          4.3350706082457, 4.3338561058628, 4.3326372377189,
                          4.3446042801037, 4.3434244460185, 4.3422402797543,
                          4.3252592918701, 4.3240097928876, 4.3227558775006,
                          4.3350706082457, 4.3338561058628, 4.3326372377189,
                          4.3446042801037, 4.3434244460185, 4.3422402797543 };

double expec_za_LST2[] = {1.1174960228303, 1.1211144827620, 1.1247310899676,
                          1.0885625249286, 1.0921953358046, 1.0958263892164,
                          1.0595178840782, 1.0631643040794, 1.0668090578716,
                          1.1174960228303, 1.1211144827620, 1.1247310899676,
                          1.0885625249286, 1.0921953358046, 1.0958263892164,
                          1.0595178840782, 1.0631643040794, 1.0668090578716,
                          1.1174960228303, 1.1211144827620, 1.1247310899676,
                          1.0885625249286, 1.0921953358046, 1.0958263892164,
                          1.0595178840782, 1.0631643040794, 1.0668090578716 };
