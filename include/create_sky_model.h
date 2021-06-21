/*! \file
  Methods to read in and crop a WODEN style sky model to everything above
  the horizon.

  @author J.L.B. Line
*/
#include "woden_struct_defs.h"

//Somthing WEIRD is up with the way I'm using the documentation package
//'breathe', so having to include the definition value in the documentation
//string to get values to appear (sigh)

/**"SOURCE" - Line beginning / containing SOURCE information (includes number of COMPONENTS) */
#define SRC_KEY         "SOURCE"
/**"ENDSOURCE" - Line ending SOURCE information */
#define SRC_END         "ENDSOURCE"
/**"COMPONENT" - Line beginning / containing COMPONENT information */
#define COMP_KEY        "COMPONENT"
/**"ENDCOMPONENT" - Line ending COMPONENT information */
#define COMP_END        "ENDCOMPONENT"
/**"FREQ" (Deprecated) - Lines contains FREQ information */
#define FREQ_KEY        "FREQ"
/**"LINEAR" - Line contains simple exponential SED information */
#define LINEAR_KEY      "LINEAR"
/**"POINT" - Line contains POINT RA, Dec information */
#define POINT_KEY       "POINT"
/**"GAUSSIAN" - Line contains GAUSSIAN RA, Dec information */
#define GAUSSIAN_KEY    "GAUSSIAN"
/**"GPARAMS" - Line contains GAUSSIAN major, minor, PA information */
#define GPARAMS_KEY     "GPARAMS"
/**"SHAPELET" - Line containing SHAPELET RA, Dec information */
#define SHAPELET_KEY    "SHAPELET"
/**"SPARAMS" - Line contains SHAPELET beta1 (major), beta2 (minor), PA information  */
#define SPARAMS_KEY     "SPARAMS"
/**"SCOEFF" - Line contains SHAPELET basis numbers n1, n2m coefficient information */
#define SCOEFF_KEY      "SCOEFF"

/**
 @brief Takes a path to WODEN-style sky model and populates a
 `source_catalogue_t` struct with the contents of `filename`.

 @details The WODEN sourcelist at `filename` should contain a number of
 sources, with the basic structure:

            SOURCE source_name P 1 G 0 S 0 0
            COMPONENT POINT 4.0 -27.0
            LINEAR 1.8e+08 10.0 0 0 0 -0.8
            ENDCOMPONENT
            ENDSOURCE

For a more detailed explanation of the sky model, please see the
documentation at @todo Link the online documentation when there is a link.

Note that `srccat` should be memory intialised, so declare with something
like `source_catalogue_t *raw_srccat = malloc( sizeof(source_catalogue_t) );`
before feeding into this function.

@param[in] *srccat Struct to contain sky model information.
@param[in] *filename Path to a WODEN-style sky model
@return Integer where 0 if read was successful, 1 if failed
 */
int read_source_catalogue(const char *filename, source_catalogue_t *srccat);

/**
enum to describe if we are cropping the sky model by entire SOURCEs, or by
invidual COMPONENTs
*/
typedef enum {CROP_SOURCES, ///< cropping by SOURCE
             CROP_COMPONENTS, ///< cropping by COMPONENT
             }e_sky_crop;

/**
enum to describe if a SOURCE/COMPONENT is above or below the horizon
*/
typedef enum {BELOW, ///< above the horizon
              ABOVE, ///< below the horizon
              }e_horizon;

/**
@brief Convert Right Ascension, Declination into Azimuth and Zenith Angle
for the given LST and Latitude

@details All angles are in radians. Uses the `erfa` function `eraHd2ae` to
perform the calculation.

@param[in] ra Right Ascension (radians)
@param[in] dec Declination (radians)
@param[in] lst Local Sidereal Time (radians)
@param[in] latitude Latitude of the array (radians)
@param[in,out] *az Azimuth (radians)
@param[in,out] *za Zenith Angle (radians)
*/
void convert_radec2azza(double ra, double dec, double lst, double latitude,
     double * az, double * za);


/**
@brief Checks if a zenith angle for a COMPONENT is above the horizon, and
returns the correct indexes of COMPONENTS to retain if so

@details Depending on what type of sky cropping we are doing (see
`crop_sky_model` function below) as specified by `sky_crop_type`. If cropping
by SOURCE, and COMPONENT below horizon, sets `all_comps_above_horizon = BELOW`.
If cropping by COMPONENT, and above horizon, update the number of
retained COMPONENTs `num_comp_retained`

If the component is a SHAPELET, we also have to update `num_shape_coeff_retained`
this is how many shapelet coefficients are retained, as one SHAPELET COMPONENT
(with one az,za) can have multiple shapelet coefficients. The array of ints
`shape_param_indexes` maps the shapelet coefficients to each SHAPELET COMPONENT,
so the function uses the int `shape` to check which coefficients match this
particular COMPONENT by looping through `shape_param_indexes`.

@param[in] za Zenith Angle (radians)
@param[in] sky_crop_type `e_sky_crop` for SOURCE or COMPONENT cropping
@param[in,out] *all_comps_above_horizon Pointer to `e_horizon` type for whether
all COMPONENTs of the SOURCE this COMPONENT belongs to are above the horizon
@param[in,out] *num_comp_retained Number of COMPONENTs retained for this SOURCE
@param[in,out] *num_shape_coeff_retained Number of SHAPELET coefficients
retained for this SOURCE
@param[in] num_shape_coeff_component The total number of shapelet coefficients
in this SOURCE
@param[in] *shape_param_indexes Map of which shapelet coefficients match this
SHAPELET COMPONENT
@param[in] shape Index of this SHAPELET COMPONENT within this SOURCE
*/
void horizon_test(double za, e_sky_crop sky_crop_type,
     e_horizon * all_comps_above_horizon, int * num_comp_retained,
     int * num_shape_coeff_retained, int num_shape_coeff_component,
     float *shape_param_indexes, int shape);

/**
@brief Takes the WODEN sky model `raw_srccat`, and crops either all SOURCEs or
COMPONENTs that are below the horizon (at the initial time step),
returning the cropped sky model. Also calculates $az,za$ for all time steps
for the retained COMPONENTs.

@details A WODEN sky model can be made of multiple SOURCEs, each of which can
have any number of COMPONENTS. We can either crop an entire SOURCE if below
horizon, or retain all COMPONENTs that are above the horizon. Say if you put
an all-sky diffuse map into one SOURCE, you should crop by COMPONENT, as some
part of the SOURCE will always be below the horizon.

A complicating factor to the cropping is that each SHAPELET COMPONENT has any
number of associated shapelet coefficients, which must be mapped correctly in
the output cropped sky model `cropped_src`. The function `horizon_test` is
used to do this mapping.

@param[in] *raw_srccat Pointer to a populated `source_catalogue_t` struct of all input
sky model parameters
@param[in] *lsts Array of local sidereal times for the simulation
@param[in] latitude Latitude of the array (radians)
@param[in] num_time_steps Number of time steps for the simulation
@param[in] sky_crop_type `e_sky_crop` for SOURCE or COMPONENT cropping

@returns `cropped_src`, a `catsource_t` sky model with only SOURCE/COMPONENTS
above the horizon for simulation

@todo Consider cropping sources iteratively for each time step. If the user runs
a very long simulation, we may end up with large chunks of empty sky. This will
however require far more `erfa` az,za calculations, so for simulations with
millions of sources, more practical to just tell the user to run multiple
shorter simulations.
*/
catsource_t * crop_sky_model(source_catalogue_t *raw_srccat, float *lsts,
              double latitude, int num_time_steps, e_sky_crop sky_crop_type);
