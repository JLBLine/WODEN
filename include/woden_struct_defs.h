#pragma once
#include <stdlib.h>
#include <stdbool.h>
#include "constants.h"

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
  float *point_ras; /*!< POINT source right ascensions (radians) */
  float *point_decs; /*!< POINT source declinations (radians) */
  float *point_ref_freqs; /*!< POINT source Flux density reference frequencies (Hz) */
  float *point_ref_stokesI; /*!< POINT source Stokes I reference flux density (Jy) */
  float *point_ref_stokesQ; /*!< POINT source Stokes Q reference flux density (Jy) */
  float *point_ref_stokesU; /*!< POINT source Stokes U reference flux density (Jy) */
  float *point_ref_stokesV; /*!< POINT source Stokes V reference flux density (Jy) */
  float *point_SIs; /*!<  POINT source spectral indexes */
  float *point_azs; /*!< POINT source azimuth angles for all time steps */
  float *point_zas; /*!< POINT source zenith angles for all time steps */
  float *sin_point_para_angs; /*!< Sine of parallatic angle for all POINT source az,za */
  float *cos_point_para_angs; /*!< Cosine of parallatic angle for all POINT source az,za */
  float *point_gaussbeam_has; /*!< Hour angle of POINT components used for Gaussian beam calculations */
  float *point_gaussbeam_decs; /*!< Declinations of POINT components used for Gaussian beam calculations */
  int num_point_primarybeam_values; /*!< Number of beam calculations needed for POINT components */

  //Gaussian params
  float *gauss_ras; /*!< GAUSSIAN source right ascensions (radians) */
  float *gauss_decs; /*!< GAUSSIAN source declinations (radians) */
  float *gauss_ref_freqs; /*!< GAUSSIAN source Flux density reference frequencies (Hz) */
  float *gauss_ref_stokesI; /*!< GAUSSIAN source Stokes I reference flux density (Jy) */
  float *gauss_ref_stokesQ; /*!< GAUSSIAN source Stokes Q reference flux density (Jy) */
  float *gauss_ref_stokesU; /*!< GAUSSIAN source Stokes U reference flux density (Jy) */
  float *gauss_ref_stokesV; /*!< GAUSSIAN source Stokes V reference flux density (Jy) */
  float *gauss_SIs; /*!<  GAUSSIAN source spectral indexes */
  float *gauss_majors; /*!< GAUSSIAN major axis (FWHM, radians) */
  float *gauss_minors; /*!< GAUSSIAN minor axis (FWHM, radians) */
  float *gauss_pas; /*!< GAUSSIAN position angles (radians) */
  float *gauss_azs; /*!< GAUSSIAN source azimuth angles for all time steps */
  float *gauss_zas; /*!< GAUSSIAN source zenith angles for all time steps */
  float *sin_gauss_para_angs; /*!< Sine of parallatic angle for all GAUSSIAN source az,za */
  float *cos_gauss_para_angs; /*!< Cosine of parallatic angle for all GAUSSIAN source az,za */
  float *gauss_gaussbeam_has; /*!< Hour angle of GAUSSIAN components used for Gaussian beam calculations */
  float *gauss_gaussbeam_decs; /*!< Declinations of GAUSSIAN components used for Gaussian beam calculations */
  int num_gauss_primarybeam_values; /*!< Number of beam calculations needed for GAUSSIAN components */

  //Shapelet params
  float *shape_ras; /*!< SHAPELET source right ascensions (radians) */
  float *shape_decs; /*!< SHAPELET source declinations (radians) */
  float *shape_ref_freqs; /*!< SHAPELET source Flux density reference frequencies (Hz) */
  float *shape_ref_stokesI; /*!< SHAPELET source Stokes I reference flux density (Jy) */
  float *shape_ref_stokesQ; /*!< SHAPELET source Stokes Q reference flux density (Jy) */
  float *shape_ref_stokesU; /*!< SHAPELET source Stokes U reference flux density (Jy) */
  float *shape_ref_stokesV; /*!< SHAPELET source Stokes V reference flux density (Jy) */
  float *shape_SIs; /*!<  SHAPELET source spectral indexes */
  float *shape_coeffs; /*!< Scaling coefficients for SHAPELET basis functions */
  float *shape_n1s; /*!< 1st basis function order for SHAPELET basis functions */
  float *shape_n2s; /*!< 2nd basis function order for SHAPELET basis functions */
  float *shape_majors; /*!< SHAPELET major axis (beta1, radians) */
  float *shape_minors; /*!< SHAPELET minor axis (beta2, radians) */
  float *shape_pas; /*!< SHAPELET position angles (radians) */
  float *shape_param_indexes; /*!< An index value to match each coeff, n1, and n2
  to the correct ra, dec, major, minor, pa for a SHAPELET */
  float *shape_azs; /*!< SHAPELET source azimuth angles for all time steps */
  float *shape_zas; /*!< SHAPELET source zenith angles for all time steps */
  float *sin_shape_para_angs; /*!< Sine of parallatic angle for all SHAPELET source az,za */
  float *cos_shape_para_angs; /*!< Cosine of parallatic angle for all SHAPELET source az,za */
  float *shape_gaussbeam_has; /*!< Hour angle of SHAPELET components used for Gaussian beam calculations */
  float *shape_gaussbeam_decs; /*!< Declinations of SHAPELET components used for Gaussian beam calculations */
  int num_shape_primarybeam_values; /*!< Number of beam calculations needed for SHAPELET components */

} catsource_t;

/*!
A struct to contain values for the MWA Fully Embbedded Element primary beam
*/
typedef struct _RTS_MWA_FEE_beam {
  double _Complex **Q1; /*!< Beam modes used for Spherical Harmonic model */
  double _Complex **Q2; /*!< Beam modes used for Spherical Harmonic model */
  double _Complex **p_T; /*!< Some pre-computed theta related values used in tile response */
  double _Complex **p_P; /*!< Some pre-computed phi related values used in tile response */
  double **M; /*!< First order of spherical harmonics */
  double **N; /*!< Second order of spherical harmonics */
  int nmax; /*!< Maximum order of spherical harmonic */
  int nMN; /*!< Total number of 1st and 2nd order harmnoic combinations */
  float _Complex norm_fac[MAX_POLS]; /*!< Zenith normalisation values */

  // BP 2019: All the Spherical Harmonic Beam data are double
  // so we will use them on the GPUs as well or there will be all kinds
  // of issues with copying

  float _Complex *d_Q1; /*!< Device copy of Q1 */
  float _Complex *d_Q2; /*!< Device copy of Q2 */
  float *d_M; /*!< Device copy of M */
  float *d_N; /*!< Device copy of N */

  float _Complex *emn_P; /*!< complex field values for phi polarisations
  separated by spherical harmonic ordering */
  float _Complex *emn_T; /*!< complex field values for theta polarisations
  separated by spherical harmonic ordering */

  float _Complex *d_emn_T_sum; /*!< complex field values for theta polarisations
  summed over spherical harmonics*/
  float _Complex *d_emn_P_sum; /*!< complex field values for phi polarisations
  summed over spherical harmonics*/

  float _Complex *rts_P1; /*!< calculated legendre polynomial values */
  float _Complex *rts_P_sin; /*!< calculated legendre polynomial / sin(theta) values*/

  float *m_range; /*!< range of possible M spherical harmonic orders */

  float _Complex *d_FEE_beam_gain_matrices; /*!< output complex gains for all
  polarisation and dipole orientation combinations on the device*/

} RTS_MWA_FEE_beam_t;

/*!
A struct to contain settings pertaining to the primary beam
*/
typedef struct _beam_settings_t {
    float gauss_sdec; /*!< Sine of the declination of the pointing for a Gaussian primary beam */
    float gauss_cdec; /*!< Cosine of the declination of the pointing for a Gaussian primary beam */
    float gauss_ha; /*!< Hour angle of the pointing for a Gaussian primary beam */

    // float *beam_point_has; /*!< Hour angle of POINT components used for Gaussian beam calculations */
    // float *beam_point_decs; /*!< Declinations of POINT components used for Gaussian beam calculations */
    // int num_point_beam_values; /*!< Number of beam calculations needed for POINT components */
    //
    // float *beam_gausscomp_has; /*!< Hour angle of GAUSSIAN components used for Gaussian beam calculations */
    // float *beam_gausscomp_decs; /*!< Declinations of GAUSSIAN components used for Gaussian beam calculations */
    // int num_gauss_beam_values; /*!< Number of beam calculations needed for GAUSSIAN components */
    //
    // float *beam_shape_has; /*!< Hour angle of SHAPELET components used for Gaussian beam calculations */
    // float *beam_shape_decs; /*!< Declinations of SHAPELET components used for Gaussian beam calculations */
    // int num_shape_beam_values; /*!< Number of beam calculations needed for SHAPELET components */

    float beam_FWHM_rad; /*!< FWHM of requested Gaussian primary beam, at reference frequnecy */
    float beam_ref_freq; /*!< Reference frequency for the given FWHM of Gaussian primary beam */
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
    beam_settings_t *beam_settings; /*!< Primary beam settings corresponding to `catsources` */
} source_catalogue_t;

/**
Struct to contain simulation parameters and visibility outputs
*/
typedef struct _visibility_set_t {
  float *us_metres; /*!< Output \f$u\f$ for all time steps, frequency steps,
  and baselines*/
  float *vs_metres; /*!< Output \f$v\f$ for all time steps, frequency steps,
  and baselines*/
  float *ws_metres; /*!< Output \f$w\f$ for all time steps, frequency steps,
  and baselines*/
  float *allsteps_sha0s; /*!< Sine of hour angle of phase centre for all
  time steps, frequency steps, and baselines*/
  float *allsteps_cha0s; /*!< Cosine of hour angle of phase centre for all
  time steps, frequency steps, and baselines*/
  float *allsteps_lsts; /*!< Local sidereal time for all time steps,
  frequency steps, and baselines (radians)*/
  float *allsteps_wavelengths; /*!< Wavelengths for all time steps,
  frequency steps, and baselines (metres)*/
  float *channel_frequencies; /*!< Frequencies for a frequency steps (Hz)*/

  float *sum_visi_XX_real; /*!< Real values for XX polarisation for all time
  steps, frequency steps, and baselines */
  float *sum_visi_XX_imag; /*!< Imaginary values for XX polarisation for all time
  steps, frequency steps, and baselines */
  float *sum_visi_XY_real; /*!< Real values for XY polarisation for all time
  steps, frequency steps, and baselines */
  float *sum_visi_XY_imag; /*!< Imaginary values for XY polarisation for all time
  steps, frequency steps, and baselines */
  float *sum_visi_YX_real; /*!< Real values for YX polarisation for all time
  steps, frequency steps, and baselines */
  float *sum_visi_YX_imag; /*!< Imaginary values for YX polarisation for all time
  steps, frequency steps, and baselines */
  float *sum_visi_YY_real; /*!< Real values for YY polarisation for all time
  steps, frequency steps, and baselines */
  float *sum_visi_YY_imag; /*!< Imaginary values for YY polarisation for all time
  steps, frequency steps, and baselines */

} visibility_set_t;

/**
Struct to contain user defined settings for simulation
*/
typedef struct _woden_settings_t {
  float lst_base; /*!< Local sidereal time for first time step (radians) */
  float ra0;  /*!< Right ascension of phase centre (radians)*/
  float dec0;  /*!< Declination of phase centre (radians)*/
  float sdec0;  /*!< Sine of Declination of phase centre (radians)*/
  float cdec0;  /*!< Cosine of Declination of phase centre (radians)*/
  int num_baselines;  /*!< Number of baselines this array layout has */
  int num_freqs;  /*!< Number of frequencies per coarse band*/
  float frequency_resolution;  /*!< Frequency resolution of a fine channel (Hz)*/
  float base_low_freq;  /*!< The lowest fine channel frequency of band 1*/
  int num_time_steps;  /*!< Number of time steps to simulate*/
  float time_res;  /*!< Time resolution of simulation (seconds)*/
  const char* cat_filename;  /*!< Path to WODEN-style sky model*/
  int num_bands;  /*!< Number of coarse frequency bands to simulate */
  int *band_nums;  /*!< Which number coarse bands to simulate (e.g 1,4,6) */
  int sky_crop_type;  /*!< Whether to crop sky models by SOURCE or COMPONENT */
  e_beamtype beamtype;  /*!< What type of primary beam to simulate with */
  float gauss_beam_FWHM;  /*!< FWHM of Gaussian primary beam (degrees)*/
  float gauss_beam_ref_freq;  /*!< Reference frequency for given Gaussian primary beam FWHM*/
  int chunking_size;  /*!< Maximum number of COMPONENTs to include in a single chunk*/
  const char* hdf5_beam_path;  /*!< Path to *.hf file containing MWA FEE beam
  spherical harmonic information*/
  double jd_date;  /*!< Julian date at beginning of simulation*/
  bool array_layout_file;  /*!< Do we have a path to the array layout or not */
  const char* array_layout_file_path;  /*!< Path to file containing E,N,H coords of array layout */
  float latitude;  /*!< Latitude of the array to simulate (radians) */
  float longitude;  /*!< Longitude of the array to simulate (radians) */
  float FEE_ideal_delays[16]; /*!< Delay values specifying the pointing for the MWA FEE beam model */
  float coarse_band_width;  /*!< Frequency bandwidth of a single coarse band (Hz)*/
  float gauss_ra_point;  /*!< The initial Right Ascension to point the Gaussian beam at (radians)*/
  float gauss_dec_point;  /*!< The initial Declination to point the Gaussian beam at (radians)*/
  int num_visis;  /*!< Total number of visiblities to simulate, so freqs*times*baselines */
  float base_band_freq;  /*!< The lowest fine channel frequency in the current band being simulated*/

} woden_settings_t;

/**
Struct to contain array layout values. Here, a single receiving element is
sometimes called an antenna, sometimes called a 'tile' (MWA lingo). This is
equivalent to a 'station' in SKA_LOW talk.
*/
typedef struct _array_layout_t {
    float *ant_X; /*!< Local \f$X\f$ location of all antenna/tiles*/
    float *ant_Y; /*!< Local \f$Y\f$ location of all antenna/tiles*/
    float *ant_Z; /*!< Local \f$Z\f$ location of all antenna/tiles*/
    float *X_diff_metres; /*!< The length of all baselines in \f$X\f$ (metres)*/
    float *Y_diff_metres; /*!< The length of all baselines in \f$Y\f$ (metres)*/
    float *Z_diff_metres; /*!< The length of all baselines in \f$Z\f$ (metres)*/
    float *ant_east; /*!< Local east location of all antenna/tiles */
    float *ant_north; /*!< Local north location of all antenna/tiles */
    float *ant_height; /*!< Local height location of all antenna/tiles */
    float latitude; /*!< Latitude of the array (radians) */
    int num_baselines; /*!< Number of baselines in the array */
    int num_tiles; /*!< Number of antenna/tiles in the array*/
    float lst_base; /*!< Local sidereal time of the first time step (radians)*/

} array_layout_t;
