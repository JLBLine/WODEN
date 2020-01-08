//Number of orders of basis functions
#define sbf_N 101
//Length of stored look up arrays
#define sbf_L 5001
//Zero pixel of the arrays
#define sbf_c 2500
//x resolution of those arrays
#define sbf_dx 0.01
//Constant for scaling
#define SQRT_M_PI_2_2_LN_2 2.668223128318498282851579 /* sqrt(pi^2 / (2 log_e(2)))*/

float * create_sbf(float *sbf);
