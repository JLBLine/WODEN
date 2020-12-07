#define _USE_MATH_DEFINES
//angle conversions
#define DH2R 0.26179938779914943653855361527329190701643078328126
#define DD2R 0.017453292519943295769236907684886127134428718885417
#define DS2R 7.2722052166430399038487115353692196393452995355905e-5
#define D2R M_PI/180.0

//astro constats
#define VELC 299792458.0
#define SOLAR2SIDEREAL 1.00274
#define DEFAULT_SI -0.8

//MWA defaults
#define MWA_LAT -26.703319
#define MWA_LAT_RAD -26.703319*DD2R

//shapelet related conversion
#define M_PI_2_2_LN_2 7.11941466249375271693034 /* pi^2 / (2 log_e(2)) */
#define SQRT_M_PI_2_2_LN_2 2.668223128318498282851579 /* sqrt(pi^2 / (2 log_e(2)))*/
#define FWHM_FACTOR 2.35482004503

//max number of visibilities to calculate simultaneously
#define MAX_CHUNKING_SIZE 260000

//RTS MWA FEE beam constants
#define N_LNA_FREQS 451
#define MAX_ZMATRIX_FREQS 256
#define NUM_DIPOLES 16
#define SIZ_Z_MATRIX (NUM_DIPOLES*2*NUM_DIPOLES*2)
#define DQ (435e-12*VELC)  // delay quantum of the MWA beamformer in meters.
#define MAX_POLS 4
#define N_COPOL 2
#define NUM_PARAMS_PER_DIPOLE 2
