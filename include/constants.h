/*! \file constants.h
  Numeric constants used throughout the code.
*/

//Somthing WEIRD is up with the way I'm using the documentation package
//'breathe', so having to include the definition value in the documentation
//string to get values to appear (sigh)

#pragma once
// #define _USE_MATH_DEFINES
//angle conversions-------------------------------------------------------------
/*! 0.26179938779914943653855361527329190701643078328126 \n
Convert hour angle to radians */
#define DH2R 0.26179938779914943653855361527329190701643078328126
/*! 0.017453292519943295769236907684886127134428718885417 \n
Covert degrees to radians */
#define DD2R 0.017453292519943295769236907684886127134428718885417
/*! 7.2722052166430399038487115353692196393452995355905e-5 \n
Convert sideral seconds to radians */
#define DS2R 7.2722052166430399038487115353692196393452995355905e-5

//astro constants---------------------------------------------------------------
/*! 299792458.0 \n
Speed of light (m/s) */
#define VELC 299792458.0
/*! 1.00274 \n
Convert solar angle to sidereal */
#define SOLAR2SIDEREAL 1.00274
/*! -0.8 \n
Default Spectral Index */
#define DEFAULT_SI -0.8

//MWA defaults------------------------------------------------------------------
/*! -26.703319405555554 \n
Latitude of the MWA (degrees) */
#define MWA_LAT -26.703319405555554
/*! -0.4660608448386394 \n
Latitude of the MWA (radians) */
#define MWA_LAT_RAD -0.4660608448386394

//shapelet related conversion---------------------------------------------------
/*! 7.11941466249375271693034 \n
\f$ \frac{\pi^2}{2 \ln(2)}  \f$ */
#define M_PI_2_2_LN_2 7.11941466249375271693034 /* pi^2 / (2 log_e(2)) */
/*! 2.668223128318498282851579 \n
\f$ \sqrt{ \frac{\pi^2}{2 \ln(2)} }  \f$ */
#define SQRT_M_PI_2_2_LN_2 2.668223128318498282851579 /* sqrt(pi^2 / (2 log_e(2)))*/
/*! 2.35482004503 \n
Convert standard deviation to FWHM for a Gaussian */
#define FWHM_FACTOR 2.35482004503
/*! 101 \n
Number of orders of basis functions (from 0 to 100 inclusive) */
#define sbf_N 101
/*! 10001 \n
Number of samples of each order of the basis functions */
#define sbf_L 10001
/*! 5000 \n
If shapelet basis function is B(x), this is the array index where x=0 */
#define sbf_c 5000
/*! 0.01 \n
If shapelet basis function is B(x), this is the x sampling resolution */
#define sbf_dx 0.01

//max number of COMPONENTs to calculate simultaneously--------------------------
/*! 260000 \n
Default number of COMPONENTs to process per chunk */
#define MAX_CHUNKING_SIZE 10000000000

//RTS MWA FEE beam constants----------------------------------------------------
/*! 16 \n
Number of dipoles in an MWA tile */
#define NUM_DIPOLES 16
/*! 435e-12*VELC \n
Delay quantum of the MWA beamformer in meters */
#define DQ (435e-12*VELC)  // delay quantum of the MWA beamformer in meters.
/*! 4 \n
Maximum number of polarisations after cross-correlation (XX, XY, YX, YY) */
#define MAX_POLS 4
/*! 2 \n
Number of polarisations before cross-correlation (X,Y) */
#define N_COPOL 2

/*! 0.29 \n
Height of MWA dipole (meters) */
#define MWA_DIPOLE_HEIGHT 0.3

/*! 1.1 \n
Separation of MWA dipoles (meters) */
#define MWA_DIPOLE_SEP 1.1


/*! 10000 \n
Initially make enough room to fit 1000000 components in a source*/
#define INITIAL_NUM_COMPONENTS 10000


/*! 100 \n
Initially make enough room to fit 100 flux entries in a list*/
#define INITIAL_NUM_FLUXES 100

/*! 200000000.0 \n
Frequency (Hz) that all power-law and curved power-law models are referenced to*/
#define REF_FREQ 200000000.0
