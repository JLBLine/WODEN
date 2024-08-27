#pragma once
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include "constants.h"
#include "woden_precision_defs.h"
#include <mwa_hyperbeam.h>

//Different
typedef enum {POINT=0, /*!< Point source type component */
              GAUSSIAN, /*!< Gaussian type component */
              SHAPELET, /*!< Shapelet type component */
              }e_component_type;

typedef enum {NO_BEAM, /*!< Do not use a primary beam in the simulation */
              GAUSS_BEAM, /*!< Use a analytic Gaussian primary beam */
              FEE_BEAM, /*!< Use the MWA FEE primary beam code */
              ANALY_DIPOLE, /*!< Use an analytic MWA dipole primary beam */
              FEE_BEAM_INTERP, /*!< Use the MWA FEE primary beam at all
              frequencies. Should be using an hdf5 file that has been frequency
              interpolated*/
              MWA_ANALY, /*!< Use an analytic MWA tile primary beam */
              EB_OSKAR, /*!< Use everybeam with the OSKAR primary beam */
              EB_LOFAR, /*!< Use everybeam with the LOFAR primary beam */
              }e_beamtype;

typedef enum {POWER_LAW=0, /*!< Power law flux model */
              CURVED_POWER_LAW, /*!< Curved power law flux model  */
              LIST, /*!< List of fluxes and freqs model */
              }e_flux_type;

/*!
A struct to contain COMPONENT information for either multiple POINT, GAUSSIAN,
or SHAPELET COMPONENTs. Three of these are included in `source_t` to specify
all COMPONENTs for a single source. This allows common values like RA,Dec,Flux
etc information to be stored the same way for different COMPONENT types
*/
typedef struct _components_t {

  //Instrinsic to COMPONENT values
  double *ras; /*!< COMPONENT right ascensions (radians) */
  double *decs; /*!< COMPONENT declinations (radians) */

  //power law params
  double *power_ref_freqs; /*!< COMPONENT Flux density reference frequencies (Hz) */
  user_precision_t *power_ref_stokesI; /*!< COMPONENT Stokes I reference flux density (Jy) */
  user_precision_t *power_SIs; /*!<  COMPONENT spectral indexes */

  // curved power law params
  double *curve_ref_freqs; /*!< COMPONENT Flux density reference frequencies (Hz) */
  user_precision_t *curve_ref_stokesI; /*!< COMPONENT Stokes I reference flux density (Jy) */
  user_precision_t *curve_SIs; /*!<  COMPONENT spectral indexes */
  user_precision_t *curve_qs; /*!<  COMPONENT curvature */

  int *power_comp_inds; /*!< The indexes of all power-law models w.r.t ra,dec */
  int *curve_comp_inds; /*!< The indexes of all curved power-law models w.r.t ra,dec */
  int *list_comp_inds; /*!< The indexes of all list models w.r.t ra,dec */

  //list flux params
  double *list_freqs; /*!< COMPONENT Flux density references frequencies (Hz) */
  user_precision_t *list_stokesI; /*!< COMPONENT Stokes I list flux density (Jy) */
  user_precision_t *list_stokesQ; /*!< COMPONENT Stokes Q list flux density (Jy) */
  user_precision_t *list_stokesU; /*!< COMPONENT Stokes U list flux density (Jy) */
  user_precision_t *list_stokesV; /*!< COMPONENT Stokes V list flux density (Jy) */
  int *num_list_values; /*!< How many freq/flux values are in each COMPONENT list*/
  int *list_start_indexes; /*!< How many freq/flux values are in each COMPONENT list*/

  int total_num_flux_entires; /*!< The total number of freq/flux values are in all lists combined*/

  //something to store extrapolated output fluxes in
  user_precision_t *extrap_stokesI; /*!< extrapolated COMPONENT Stokes I flux densities (Jy) */
  user_precision_t *extrap_stokesQ; /*!< extrapolated COMPONENT Stokes Q flux densities (Jy) */
  user_precision_t *extrap_stokesU; /*!< extrapolated COMPONENT Stokes U flux densities (Jy) */
  user_precision_t *extrap_stokesV; /*!< extrapolated COMPONENT Stokes V flux densities (Jy) */

  //SHAPELET / GAUSSIAN params
  user_precision_t *shape_coeffs; /*!< Scaling coefficients for SHAPELET basis functions */
  user_precision_t *n1s; /*!< 1st basis function order for SHAPELET basis functions */
  user_precision_t *n2s; /*!< 2nd basis function order for SHAPELET basis functions */
  user_precision_t *majors; /*!< GAUSSIAN/SHAPELET major axis (beta1, radians) */
  user_precision_t *minors; /*!< GAUSSIAN/SHAPELET minor axis (beta2, radians) */
  user_precision_t *pas; /*!< GAUSSIAN/SHAPELET position angles (radians) */
  user_precision_t *param_indexes; /*!< An index value to match each coeff, n1, and n2
  to the correct ra, dec, major, minor, pa for a SHAPELET */

  //Specific to observation settings for these COMPONENTs
  user_precision_t *azs; /*!< SHAPELET source azimuth angles for all time steps */
  user_precision_t *zas; /*!< SHAPELET source zenith angles for all time steps */
  double *beam_has; /*!< Hour angle of COMPONENTs for all time steps, used for
   beam calculations */
  double *beam_decs; /*!< Declinations of COMPONENTs for all time steps, used for
   beam calculations */
  int num_primarybeam_values; /*!< Number of beam calculations needed for
  COMPONENTs */

  /*
  Things to hold beam gain values
  */
  user_precision_complex_t *gxs; /*!< North-South Beam gain values for all directions,
  frequencies, and times for these COMPONENTS*/
  user_precision_complex_t *Dxs; /*!< North-South Beam leakage values for all directions,
  frequencies, and times for these COMPONENTS*/
  user_precision_complex_t *Dys; /*!< East-West Beam leakage values for all directions,
  frequencies, and times for these COMPONENTS*/
  user_precision_complex_t *gys; /*!< East-West Beam gain values for all directions,
  frequencies, and times for these COMPONENTS*/

  /*
  Things to hold beam gain values when you have different beams per antenna
  */
  user_precision_complex_t *gxs_ants; /*!< North-South Beam gain values for all directions,
  frequencies, times, and antennas for these COMPONENTS*/
  user_precision_complex_t *Dxs_ants; /*!< North-South Beam leakage values for all directions,
  frequencies, times, and antennas for these COMPONENTS*/
  user_precision_complex_t *Dys_ants; /*!< East-West Beam leakage values for all directions,
  frequencies, times, and antennas for these COMPONENTS*/
  user_precision_complex_t *gys_ants; /*!< East-West Beam gain values for all directions,
  frequencies, times, and antennas for these COMPONENTS*/

  //Leave off the d_ from these device values, as the components_t struct
  //itself will have the d_ label if doing things on the GPU
  double *ls; /*!< Device memory l cosine direction coords for these COMPONENTs*/
  double *ms; /*!< Device memory m cosine direction coords for these COMPONENTs*/
  double *ns; /*!< Device memory n cosine direction coords for these COMPONENTs*/

  user_precision_t *stokesV_pol_fracs; /*!< Stokes V polarisation fractions */
  int *stokesV_pol_frac_comp_inds; /*!< The indexes of all Stokes V polarisation fraction models w.r.t ra,dec */
  user_precision_t *stokesV_power_ref_flux; /*!< Stokes V reference flux for power-law */
  user_precision_t *stokesV_power_SIs; /*!< Stokes V spectral index for power-law */
  int *stokesV_power_comp_inds; /*!< The indexes of all Stokes V power-law models w.r.t ra,dec */
  user_precision_t *stokesV_curve_ref_flux; /*!< Stokes V reference flux for curved power-law */
  user_precision_t *stokesV_curve_SIs; /*!< Stokes V spectral index for curved power-law */
  user_precision_t *stokesV_curve_qs; /*!< Stokes V q param for curved power-law */
  int *stokesV_curve_comp_inds; /*!< The indexes of Stokes V curved power-law models w.r.t ra,dec */
  double *stokesV_list_ref_freqs; /*!< Stokes V reference frequencies for list type models */
  user_precision_t* stokesV_list_ref_flux; /*!< Stokes V reference fluxes for list type models */
  int *stokesV_list_comp_inds; /*!< The indexes of all Stokes V list type models w.r.t ra,dec */
  int *stokesV_num_list_values; /*!< How many freq/flux values there are for Stokes V list-type component*/
  int *stokesV_list_start_indexes; /*!< Where the Stokes V list-type component starts in the list flux arrays */
  user_precision_t *linpol_pol_fracs; /*!< Linear polarisation polarisation fractions */
  int *linpol_pol_frac_comp_inds; /*!< The indexes of all linear polarisation fraction models w.r.t ra,dec */
  user_precision_t *linpol_power_ref_flux; /*!< Linear polarisation reference flux for power-law */
  user_precision_t *linpol_power_SIs; /*!< Linear polarisation spectral index for power-law */
  int *linpol_power_comp_inds;  /*!< The indexes of all linear polarisation power-law models w.r.t ra,dec */
  user_precision_t *linpol_curve_ref_flux; /*!< Linear polarisation reference flux for curved power-law */
  user_precision_t *linpol_curve_SIs; /*!< Linear polarisation spectral index for curved power-law */
  user_precision_t *linpol_curve_qs; /*!< Linear polarisation q param for curved power-law */
  int *linpol_curve_comp_inds; /*!< The indexes of all linear polarisation curved power-law models w.r.t ra,dec */
  double *stokesQ_list_ref_freqs; /*!< Stokes Q reference frequencies for list type models */
  user_precision_t* stokesQ_list_ref_flux; /*!< Stokes Q reference fluxes for list type models */
  int *stokesQ_list_comp_inds; /*!< The indexes of all Stokes Q list type models w.r.t ra,dec */
  int *stokesQ_num_list_values; /*!< How many freq/flux values there are for Stokes Q list-type component*/
  int *stokesQ_list_start_indexes; /*!< Where the Stokes Q list-type component starts in the list flux arrays */
  double *stokesU_list_ref_freqs; /*!< Stokes U reference frequencies for list type models */
  user_precision_t* stokesU_list_ref_flux; /*!< Stokes U reference fluxes for list type models */
  int *stokesU_list_comp_inds; /*!< The indexes of all Stokes U list type models w.r.t ra,dec */
  int *stokesU_num_list_values; /*!< How many freq/flux values there are for Stokes U list-type component*/
  int *stokesU_list_start_indexes; /*!< Where the Stokes U list-type component starts in the list flux arrays */
  double *linpol_p_list_ref_freqs; /*!< Linear polarisation flux reference frequencies for list type models */
  user_precision_t* linpol_p_list_ref_flux; /*!< Linear polarisation flux reference fluxes for list type models */
  int *linpol_p_list_comp_inds; /*!< The indexes of all Linear polarisation flux list type models w.r.t ra,dec */
  int *linpol_p_num_list_values; /*!< How many freq/flux values there are for Linear polarisation flux list-type component*/
  int *linpol_p_list_start_indexes; /*!< Where the Linear polarisation flux list-type component starts in the list flux arrays */
  user_precision_t *rm_values; /*!< Linear polarisation rotation measures */
  user_precision_t *intr_pol_angle; /*!<  Linear polarisation instrinsic polarisation angles */
  int *linpol_angle_inds;  /*!< The indexes of all RM/intrinsic polarisation angles w.r.t ra,dec */
  int n_stokesV_pol_frac; /*!< The number of Stokes V polarisation fraction models */
  int n_stokesV_power; /*!< The number of Stokes V power-law models */
  int n_stokesV_curve; /*!< The number of Stokes V curved power-law models */
  int n_stokesV_list; /*!< The number of Stokes V list type models */
  int n_stokesV_list_flux_entries; /*!< The total number of Stokes V list type flux entries */
  int n_linpol_pol_frac;  /*!< The number of linear polarisation fraction models */
  int n_linpol_power; /*!< The number of linear polarisation power-law models */
  int n_linpol_curve; /*!< The number of linear polarisation curved power-law models */
  int n_linpol_list; /*!< The number of linear pol list type models */
  int n_stokesQ_list_flux_entries; /*!< The total number of Stokes Q list type flux entries */
  int n_stokesU_list_flux_entries;  /*!< The total number of Stokes U list type flux entries */
  int n_linpol_p_list; /*!< The number of linear polarised flux list type models */
  int n_linpol_p_list_flux_entries;  /*!< The total number of Polarised flux list type flux entries */
  int n_linpol_angles; /*!< The number of RM/intrinsic polarisation angles */
  int do_QUV; /*!< Set if doing any polarised information */

} components_t;



/*!
A struct to contain sky model values for a single SOURCE
*/
typedef struct _source_t {
  //General source info
  char name[32]; /*!< Source name */
  int n_comps; /*!< Total number of COMPONENTs in source  */

  int n_points; /*!< Number of POINT source COMPONENTs  */
  int n_point_lists; /*!< Number of POINTs with LIST type flux */
  int n_point_powers; /*!< Number of POINTs with POWER_LAW type flux */
  int n_point_curves; /*!< Number of POINTs with CURVED_POWER_LAW type flux */

  int n_gauss; /*!< Number of GAUSSIAN source COMPONENTs */
  int n_gauss_lists; /*!< Number of GAUSSIANs with LIST type flux */
  int n_gauss_powers; /*!< Number of GAUSSIANs with POWER_LAW type flux */
  int n_gauss_curves; /*!< Number of GAUSSIANs with CURVED_POWER_LAW type flux */

  int n_shapes; /*!< Number of SHAPELET source COMPONENTs */
  int n_shape_lists; /*!< Number of SHAPELETs with LIST type flux */
  int n_shape_powers; /*!< Number of SHAPELETs with POWER_LAW type flux */
  int n_shape_curves; /*!< Number of SHAPELETs with CURVED_POWER_LAW type flux */
  int n_shape_coeffs; /*!< Total number of SHAPELET coefficients */

  components_t point_components; /*!< `components_t` holding component
  information for all POINT COMPONENTs in this SOURCE.*/
  components_t gauss_components; /*!< `components_t` holding component
  information for all GAUSSIAN COMPONENTs in this SOURCE.*/
  components_t shape_components; /*!< `components_t` holding component
  information for all SHAPELET COMPONENTs in this SOURCE.*/

  //Device versions
  components_t d_point_components; /*!< `components_t` holding component
  information for all POINT COMPONENTs in this SOURCE.*/
  components_t d_gauss_components; /*!< `components_t` holding component
  information for all GAUSSIAN COMPONENTs in this SOURCE.*/
  components_t d_shape_components; /*!< `components_t` holding component
  information for all SHAPELET COMPONENTs in this SOURCE.*/

} source_t;

/*!
A struct to contain multiple `source_t` type sky models and `beam_settings_t`
primary beam settings, to be iterated over by `calculate_visibilities`
*/
typedef struct _source_catalogue_t {
    int num_sources; /*!< Number of SOURCES in this `source_catalogue_t`*/
    int num_shapelets; /*!< Total number of SHAPELET components in this `source_catalogue_t` */
    source_t *sources; /*!< Multiple sky models to simulate */
    // beam_settings_t *beam_settings; /*!< Primary beam settings corresponding to `sources` */
} source_catalogue_t;

/*!
A struct to contain settings pertaining to the primary beam
*/
typedef struct _beam_settings_t {
    user_precision_t gauss_sdec; /*!< Sine of the declination of the pointing for a Gaussian primary beam */
    user_precision_t gauss_cdec; /*!< Cosine of the declination of the pointing for a Gaussian primary beam */
    double gauss_ha; /*!< Hour angle of the pointing for a Gaussian primary beam */

    user_precision_t beam_FWHM_rad; /*!< FWHM of requested Gaussian primary beam, at reference frequnecy */
    double beam_ref_freq; /*!< Reference frequency for the given FWHM of Gaussian primary beam */
    e_beamtype beamtype; /*!< What type of primary beam to simulate - see `e_beamtype` */

    user_precision_t *MWAFEE_freqs; /*!< The frequencies of the initialised MWAFEE beams
    in FEE_beams */

    int num_MWAFEE; /*!< Number of MWAFEE beam instances to cover all desired frequencies */

    struct FEEBeamGpu *gpu_fee_beam; /*!< Single initialised hyperbeam device model for desired pointing */
    struct FEEBeam *fee_beam; /*!< Single initialised hyperbeam host model for desired pointing */

    char hyper_error_str[100]; /*!< Char array to hold error messages out of hyperbeam */

    double base_middle_freq; /*!< The frequency at the middle of the base coarse band */
    uint32_t *hyper_delays; /*!< MWA FEE delays in a format that hyperbeam likes */


} beam_settings_t;

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
  double *allsteps_sha0s; /*!< Sine of hour angle of phase centre for all
  time steps, frequency steps, and baselines*/
  double *allsteps_cha0s; /*!< Cosine of hour angle of phase centre for all
  time steps, frequency steps, and baselines*/
  double *allsteps_lsts; /*!< Local sidereal time for all time steps,
  frequency steps, and baselines (radians)*/
  user_precision_t *allsteps_wavelengths; /*!< Wavelengths for all time steps,
  frequency steps, and baselines (metres)*/
  double *channel_frequencies; /*!< Frequencies for a frequency steps (Hz)*/

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
  double lst_base; /*!< Local sidereal time for first time step (radians) */
  double lst_obs_epoch_base; /*!< Local sidereal time for first time step (radians)
  for the observation epoch (e.g. in 2020 for a 2020 obs) */
  double ra0;  /*!< Right ascension of phase centre (radians)*/
  double dec0;  /*!< Declination of phase centre (radians)*/
  double sdec0;  /*!< Sine of Declination of phase centre (radians)*/
  double cdec0;  /*!< Cosine of Declination of phase centre (radians)*/
  int num_baselines;  /*!< Number of baselines this array layout has */
  int num_ants;  /*!< Number of antennas this array layout has (MWA calls this
  number of tiles) */
  int num_freqs;  /*!< Number of frequencies per coarse band*/
  double frequency_resolution;  /*!< Frequency resolution of a fine channel (Hz)*/
  double base_low_freq;  /*!< The lowest fine channel frequency of band 1*/
  int num_time_steps;  /*!< Number of time steps to simulate*/
  double time_res;  /*!< Time resolution of simulation (seconds)*/
  const char* cat_filename;  /*!< Path to WODEN-style sky model*/
  int num_bands;  /*!< Number of coarse frequency bands to simulate */
  int *band_nums;  /*!< Which number coarse bands to simulate (e.g 1,4,6) */
  int sky_crop_type;  /*!< Whether to crop sky models by SOURCE or COMPONENT */
  e_beamtype beamtype;  /*!< What type of primary beam to simulate with */
  user_precision_t gauss_beam_FWHM;  /*!< FWHM of Gaussian primary beam (degrees)*/
  double gauss_beam_ref_freq;  /*!< Reference frequency for given Gaussian primary beam FWHM*/
  long int chunking_size;  /*!< Maximum number of COMPONENTs to include in a single chunk*/
  const char* hdf5_beam_path;  /*!< Path to *.hf file containing MWA FEE beam
  spherical harmonic information*/
  double jd_date;  /*!< Julian date at beginning of simulation*/
  bool array_layout_file;  /*!< Do we have a path to the array layout or not */
  const char* array_layout_file_path;  /*!< Path to file containing E,N,H coords of array layout */
  double latitude;  /*!< Latitude of the array to simulate (radians) */
  double latitude_obs_epoch_base;   /*!< Latitude of the array at the observation epoch (radians) */
  user_precision_t longitude;  /*!< Longitude of the array to simulate (radians) */
  int *FEE_ideal_delays; /*!< Delay values specifying the pointing for the MWA FEE beam model.
  Should be 16*num_beams in length*/
  double coarse_band_width;  /*!< Frequency bandwidth of a single coarse band (Hz)*/
  double gauss_ra_point;  /*!< The initial Right Ascension to point the Gaussian beam at (radians)*/
  double gauss_dec_point;  /*!< The initial Declination to point the Gaussian beam at (radians)*/
  int num_cross;  /*!< Total number of cross-correlations to simulate, so freqs*times*baselines */
  int num_autos;  /*!< Total number of auto-correlations to simulate, so freqs*times*baselines */
  int num_visis;  /*!< Total number of visibilities to simulate, so num_cross + num_autos */

  double base_band_freq;  /*!< The lowest fine channel frequency in the current band being simulated*/
  int do_precession; /*!< Boolean of whether to apply precession to the
  array layout or not*/
  double *lsts; /*!< Array to hold LSTs for all time centroids (these are
    different when precession is happening)*/
  double *latitudes; /*!< Array to hold latitudes for all time centroids (these
    are different when precession is happening)*/
  double *mjds; /*!< Array to hold modified julian dates for all time centroids*/
  int do_autos; /*!< Boolean of whether to simulate autos or not (0 False, 1 True)*/
  int use_dipamps; /*!< Boolean of whether to use dipole amplitudes, so have an individual beam per tile*/
  double *mwa_dipole_amps; /*!< Bespoke MWA dipole amplitudes for each antenna(tile). Should be 2*num_ants*16 long*/
  int single_everybeam_station; /*!< If using everybeam, add this to say we are only using a single station*/

} woden_settings_t;

/**
Struct to contain array layout values. Here, a single receiving element is
sometimes called an antenna, sometimes called a 'tile' (MWA lingo). This is
equivalent to a 'station' in SKA_LOW talk.
*/
typedef struct _array_layout_t {
    double *ant_X; /*!< Local \f$X\f$ location of all antenna/tiles*/
    double *ant_Y; /*!< Local \f$Y\f$ location of all antenna/tiles*/
    double *ant_Z; /*!< Local \f$Z\f$ location of all antenna/tiles*/
    double *X_diff_metres; /*!< The length of all baselines in \f$X\f$ (metres)*/
    double *Y_diff_metres; /*!< The length of all baselines in \f$Y\f$ (metres)*/
    double *Z_diff_metres; /*!< The length of all baselines in \f$Z\f$ (metres)*/
    double *ant_east; /*!< Local east location of all antenna/tiles */
    double *ant_north; /*!< Local north location of all antenna/tiles */
    double *ant_height; /*!< Local height location of all antenna/tiles */
    double latitude; /*!< Latitude of the array (radians) */
    int num_baselines; /*!< Number of baselines in the array */
    int num_tiles; /*!< Number of antenna/tiles in the array*/
    double lst_base; /*!< Local sidereal time of the first time step (radians)*/

} array_layout_t;
