#pragma once
#include <stdlib.h>
#include <stdbool.h>
#include "constants.h"
#include "woden_precision_defs.h"

enum component_type {POINT=0, /*!< Point source type component */
                     GAUSSIAN, /*!< Gaussian type component */
                     SHAPELET, /*!< Shapelet type component */
                     };
typedef enum {NO_BEAM, /*!< Do not use a primary beam in the simulation */
              GAUSS_BEAM, /*!< Use a analytic Gaussian primary beam */
              FEE_BEAM, /*!< Use the RTS MWA FEE primary beam code */
              ANALY_DIPOLE, /*!< Use an analytic MWA dipole primary beam */
              }e_beamtype;

/*!
A struct to contain sky model values for a single catalogue SOURCE
*/
typedef struct _catsource_t {
  //General source info
  char name[32]; /*!< Source name */
  int n_comps; /*!< Total number of COMPONENTs in source  */
  int n_points; /*!< Number of POINT source COMPONENTs  */
  int n_gauss; /*!< Number of GAUSSIAN source COMPONENTs */
  int n_shapes; /*!< Number of SHAPELET source COMPONENTs */
  int n_shape_coeffs; /*!< Total number of SHAPELET coefficients */

  //Pointsource params
  double *point_ras; /*!< POINT source right ascensions (radians) */
  double *point_decs; /*!< POINT source declinations (radians) */
  user_precision_t *point_ref_freqs; /*!< POINT source Flux density reference frequencies (Hz) */
  user_precision_t *point_ref_stokesI; /*!< POINT source Stokes I reference flux density (Jy) */
  user_precision_t *point_ref_stokesQ; /*!< POINT source Stokes Q reference flux density (Jy) */
  user_precision_t *point_ref_stokesU; /*!< POINT source Stokes U reference flux density (Jy) */
  user_precision_t *point_ref_stokesV; /*!< POINT source Stokes V reference flux density (Jy) */
  user_precision_t *point_SIs; /*!<  POINT source spectral indexes */
  user_precision_t *point_azs; /*!< POINT source azimuth angles for all time steps */
  user_precision_t *point_zas; /*!< POINT source zenith angles for all time steps */
  user_precision_t *sin_point_para_angs; /*!< Sine of parallatic angle for all POINT source az,za */
  user_precision_t *cos_point_para_angs; /*!< Cosine of parallatic angle for all POINT source az,za */
  double *point_gaussbeam_has; /*!< Hour angle of POINT components used for Gaussian beam calculations */
  double *point_gaussbeam_decs; /*!< Declinations of POINT components used for Gaussian beam calculations */
  int num_point_primarybeam_values; /*!< Number of beam calculations needed for POINT components */

  //Gaussian params
  double *gauss_ras; /*!< GAUSSIAN source right ascensions (radians) */
  double *gauss_decs; /*!< GAUSSIAN source declinations (radians) */
  user_precision_t *gauss_ref_freqs; /*!< GAUSSIAN source Flux density reference frequencies (Hz) */
  user_precision_t *gauss_ref_stokesI; /*!< GAUSSIAN source Stokes I reference flux density (Jy) */
  user_precision_t *gauss_ref_stokesQ; /*!< GAUSSIAN source Stokes Q reference flux density (Jy) */
  user_precision_t *gauss_ref_stokesU; /*!< GAUSSIAN source Stokes U reference flux density (Jy) */
  user_precision_t *gauss_ref_stokesV; /*!< GAUSSIAN source Stokes V reference flux density (Jy) */
  user_precision_t *gauss_SIs; /*!<  GAUSSIAN source spectral indexes */
  user_precision_t *gauss_majors; /*!< GAUSSIAN major axis (FWHM, radians) */
  user_precision_t *gauss_minors; /*!< GAUSSIAN minor axis (FWHM, radians) */
  user_precision_t *gauss_pas; /*!< GAUSSIAN position angles (radians) */
  user_precision_t *gauss_azs; /*!< GAUSSIAN source azimuth angles for all time steps */
  user_precision_t *gauss_zas; /*!< GAUSSIAN source zenith angles for all time steps */
  user_precision_t *sin_gauss_para_angs; /*!< Sine of parallatic angle for all GAUSSIAN source az,za */
  user_precision_t *cos_gauss_para_angs; /*!< Cosine of parallatic angle for all GAUSSIAN source az,za */
  double *gauss_gaussbeam_has; /*!< Hour angle of GAUSSIAN components used for Gaussian beam calculations */
  double *gauss_gaussbeam_decs; /*!< Declinations of GAUSSIAN components used for Gaussian beam calculations */
  int num_gauss_primarybeam_values; /*!< Number of beam calculations needed for GAUSSIAN components */

  //Shapelet params
  double *shape_ras; /*!< SHAPELET source right ascensions (radians) */
  double *shape_decs; /*!< SHAPELET source declinations (radians) */
  user_precision_t *shape_ref_freqs; /*!< SHAPELET source Flux density reference frequencies (Hz) */
  user_precision_t *shape_ref_stokesI; /*!< SHAPELET source Stokes I reference flux density (Jy) */
  user_precision_t *shape_ref_stokesQ; /*!< SHAPELET source Stokes Q reference flux density (Jy) */
  user_precision_t *shape_ref_stokesU; /*!< SHAPELET source Stokes U reference flux density (Jy) */
  user_precision_t *shape_ref_stokesV; /*!< SHAPELET source Stokes V reference flux density (Jy) */
  user_precision_t *shape_SIs; /*!<  SHAPELET source spectral indexes */
  user_precision_t *shape_coeffs; /*!< Scaling coefficients for SHAPELET basis functions */
  user_precision_t *shape_n1s; /*!< 1st basis function order for SHAPELET basis functions */
  user_precision_t *shape_n2s; /*!< 2nd basis function order for SHAPELET basis functions */
  user_precision_t *shape_majors; /*!< SHAPELET major axis (beta1, radians) */
  user_precision_t *shape_minors; /*!< SHAPELET minor axis (beta2, radians) */
  user_precision_t *shape_pas; /*!< SHAPELET position angles (radians) */
  user_precision_t *shape_param_indexes; /*!< An index value to match each coeff, n1, and n2
  to the correct ra, dec, major, minor, pa for a SHAPELET */
  user_precision_t *shape_azs; /*!< SHAPELET source azimuth angles for all time steps */
  user_precision_t *shape_zas; /*!< SHAPELET source zenith angles for all time steps */
  user_precision_t *sin_shape_para_angs; /*!< Sine of parallatic angle for all SHAPELET source az,za */
  user_precision_t *cos_shape_para_angs; /*!< Cosine of parallatic angle for all SHAPELET source az,za */
  double *shape_gaussbeam_has; /*!< Hour angle of SHAPELET components used for Gaussian beam calculations */
  double *shape_gaussbeam_decs; /*!< Declinations of SHAPELET components used for Gaussian beam calculations */
  int num_shape_primarybeam_values; /*!< Number of beam calculations needed for SHAPELET components */

} catsource_t;

/*!
A struct to contain values for the MWA Fully Embbedded Element primary beam
*/
typedef struct _RTS_MWA_FEE_beam {
  double _Complex **Q1; /*!< Beam modes used for Spherical Harmonic model */
  double _Complex **Q2; /*!< Beam modes used for Spherical Harmonic model */
  double **M; /*!< First order of spherical harmonics */
  double **N; /*!< Second order of spherical harmonics */
  int nmax; /*!< Maximum order of spherical harmonic */
  int nMN; /*!< Total number of 1st and 2nd order harmnoic combinations */
  user_precision_complex_t norm_fac[MAX_POLS]; /*!< Zenith normalisation values */

  // BP 2019: All the Spherical Harmonic Beam data are double
  // so we will use them on the GPUs as well or there will be all kinds
  // of issues with copying

  user_precision_complex_t *d_Q1; /*!< Device copy of Q1 */
  user_precision_complex_t *d_Q2; /*!< Device copy of Q2 */
  user_precision_t *d_M; /*!< Device copy of M */
  user_precision_t *d_N; /*!< Device copy of N */

  user_precision_complex_t *emn_P; /*!< complex field values for phi polarisations
  separated by spherical harmonic ordering */
  user_precision_complex_t *emn_T; /*!< complex field values for theta polarisations
  separated by spherical harmonic ordering */

  user_precision_complex_t *d_emn_T_sum; /*!< complex field values for theta polarisations
  summed over spherical harmonics*/
  user_precision_complex_t *d_emn_P_sum; /*!< complex field values for phi polarisations
  summed over spherical harmonics*/

  user_precision_complex_t *rts_P1; /*!< calculated legendre polynomial values */
  user_precision_complex_t *rts_P_sin; /*!< calculated legendre polynomial / sin(theta) values*/

  user_precision_t *m_range; /*!< range of possible M spherical harmonic orders */

  user_precision_complex_t *d_FEE_beam_gain_matrices; /*!< output complex gains for all
  polarisation and dipole orientation combinations on the device*/

} RTS_MWA_FEE_beam_t;

/*!
A struct to contain settings pertaining to the primary beam
*/
typedef struct _beam_settings_t {
    user_precision_t gauss_sdec; /*!< Sine of the declination of the pointing for a Gaussian primary beam */
    user_precision_t gauss_cdec; /*!< Cosine of the declination of the pointing for a Gaussian primary beam */
    user_precision_t gauss_ha; /*!< Hour angle of the pointing for a Gaussian primary beam */

    user_precision_t beam_FWHM_rad; /*!< FWHM of requested Gaussian primary beam, at reference frequnecy */
    user_precision_t beam_ref_freq; /*!< Reference frequency for the given FWHM of Gaussian primary beam */
    int beamtype; /*!< What type of primary beam to simulate - see `e_beamtype` */

    RTS_MWA_FEE_beam_t *FEE_beam; /*!< Initialised MWA FEE beam model for desired pointing */
    RTS_MWA_FEE_beam_t *FEE_beam_zenith; /*!< Initialised MWA FEE beam model
    for zenith pointing, used for normalisation of the desired pointing */

} beam_settings_t;

/*!
A struct to contain multiple `catsource_t` type sky models and `beam_settings_t`
primary beam settings, to be iterated over by `calculate_visibilities`
*/
typedef struct _source_catalogue_t {
    int num_sources; /*!< Number of SOURCES in this `source_catalogue_t`*/
    int num_shapelets; /*!< Total number of SHAPELET components in this `source_catalogue_t` */
    catsource_t *catsources; /*!< Multiple sky models to simulate */
    // beam_settings_t *beam_settings; /*!< Primary beam settings corresponding to `catsources` */
} source_catalogue_t;

/**
Struct to contain simulation parameters and visibility outputs
*/
typedef struct _visibility_set_t {
  user_precision_t *us_metres; /*!< Output \f$u\f$ for all time steps, frequency steps,
  and baselines*/
  user_precision_t *vs_metres; /*!< Output \f$v\f$ for all time steps, frequency steps,
  and baselines*/
  user_precision_t *ws_metres; /*!< Output \f$w\f$ for all time steps, frequency steps,
  and baselines*/
  user_precision_t *allsteps_sha0s; /*!< Sine of hour angle of phase centre for all
  time steps, frequency steps, and baselines*/
  user_precision_t *allsteps_cha0s; /*!< Cosine of hour angle of phase centre for all
  time steps, frequency steps, and baselines*/
  user_precision_t *allsteps_lsts; /*!< Local sidereal time for all time steps,
  frequency steps, and baselines (radians)*/
  user_precision_t *allsteps_wavelengths; /*!< Wavelengths for all time steps,
  frequency steps, and baselines (metres)*/
  user_precision_t *channel_frequencies; /*!< Frequencies for a frequency steps (Hz)*/

  user_precision_t *sum_visi_XX_real; /*!< Real values for XX polarisation for all time
  steps, frequency steps, and baselines */
  user_precision_t *sum_visi_XX_imag; /*!< Imaginary values for XX polarisation for all time
  steps, frequency steps, and baselines */
  user_precision_t *sum_visi_XY_real; /*!< Real values for XY polarisation for all time
  steps, frequency steps, and baselines */
  user_precision_t *sum_visi_XY_imag; /*!< Imaginary values for XY polarisation for all time
  steps, frequency steps, and baselines */
  user_precision_t *sum_visi_YX_real; /*!< Real values for YX polarisation for all time
  steps, frequency steps, and baselines */
  user_precision_t *sum_visi_YX_imag; /*!< Imaginary values for YX polarisation for all time
  steps, frequency steps, and baselines */
  user_precision_t *sum_visi_YY_real; /*!< Real values for YY polarisation for all time
  steps, frequency steps, and baselines */
  user_precision_t *sum_visi_YY_imag; /*!< Imaginary values for YY polarisation for all time
  steps, frequency steps, and baselines */

} visibility_set_t;

/**
Struct to contain user defined settings for simulation
*/
typedef struct _woden_settings_t {
  user_precision_t lst_base; /*!< Local sidereal time for first time step (radians) */
  user_precision_t ra0;  /*!< Right ascension of phase centre (radians)*/
  user_precision_t dec0;  /*!< Declination of phase centre (radians)*/
  user_precision_t sdec0;  /*!< Sine of Declination of phase centre (radians)*/
  user_precision_t cdec0;  /*!< Cosine of Declination of phase centre (radians)*/
  int num_baselines;  /*!< Number of baselines this array layout has */
  int num_freqs;  /*!< Number of frequencies per coarse band*/
  user_precision_t frequency_resolution;  /*!< Frequency resolution of a fine channel (Hz)*/
  user_precision_t base_low_freq;  /*!< The lowest fine channel frequency of band 1*/
  int num_time_steps;  /*!< Number of time steps to simulate*/
  user_precision_t time_res;  /*!< Time resolution of simulation (seconds)*/
  const char* cat_filename;  /*!< Path to WODEN-style sky model*/
  int num_bands;  /*!< Number of coarse frequency bands to simulate */
  int *band_nums;  /*!< Which number coarse bands to simulate (e.g 1,4,6) */
  int sky_crop_type;  /*!< Whether to crop sky models by SOURCE or COMPONENT */
  e_beamtype beamtype;  /*!< What type of primary beam to simulate with */
  user_precision_t gauss_beam_FWHM;  /*!< FWHM of Gaussian primary beam (degrees)*/
  user_precision_t gauss_beam_ref_freq;  /*!< Reference frequency for given Gaussian primary beam FWHM*/
  long int chunking_size;  /*!< Maximum number of COMPONENTs to include in a single chunk*/
  const char* hdf5_beam_path;  /*!< Path to *.hf file containing MWA FEE beam
  spherical harmonic information*/
  double jd_date;  /*!< Julian date at beginning of simulation*/
  bool array_layout_file;  /*!< Do we have a path to the array layout or not */
  const char* array_layout_file_path;  /*!< Path to file containing E,N,H coords of array layout */
  double latitude;  /*!< Latitude of the array to simulate (radians) */
  user_precision_t longitude;  /*!< Longitude of the array to simulate (radians) */
  user_precision_t FEE_ideal_delays[16]; /*!< Delay values specifying the pointing for the MWA FEE beam model */
  user_precision_t coarse_band_width;  /*!< Frequency bandwidth of a single coarse band (Hz)*/
  user_precision_t gauss_ra_point;  /*!< The initial Right Ascension to point the Gaussian beam at (radians)*/
  user_precision_t gauss_dec_point;  /*!< The initial Declination to point the Gaussian beam at (radians)*/
  int num_visis;  /*!< Total number of visiblities to simulate, so freqs*times*baselines */
  user_precision_t base_band_freq;  /*!< The lowest fine channel frequency in the current band being simulated*/
  int do_precession; /*!< Boolean of whether to apply precession to the
  array layout or not*/

} woden_settings_t;

/**
Struct to contain array layout values. Here, a single receiving element is
sometimes called an antenna, sometimes called a 'tile' (MWA lingo). This is
equivalent to a 'station' in SKA_LOW talk.
*/
typedef struct _array_layout_t {
    user_precision_t *ant_X; /*!< Local \f$X\f$ location of all antenna/tiles*/
    user_precision_t *ant_Y; /*!< Local \f$Y\f$ location of all antenna/tiles*/
    user_precision_t *ant_Z; /*!< Local \f$Z\f$ location of all antenna/tiles*/
    user_precision_t *X_diff_metres; /*!< The length of all baselines in \f$X\f$ (metres)*/
    user_precision_t *Y_diff_metres; /*!< The length of all baselines in \f$Y\f$ (metres)*/
    user_precision_t *Z_diff_metres; /*!< The length of all baselines in \f$Z\f$ (metres)*/
    user_precision_t *ant_east; /*!< Local east location of all antenna/tiles */
    user_precision_t *ant_north; /*!< Local north location of all antenna/tiles */
    user_precision_t *ant_height; /*!< Local height location of all antenna/tiles */
    user_precision_t latitude; /*!< Latitude of the array (radians) */
    int num_baselines; /*!< Number of baselines in the array */
    int num_tiles; /*!< Number of antenna/tiles in the array*/
    user_precision_t lst_base; /*!< Local sidereal time of the first time step (radians)*/

} array_layout_t;
