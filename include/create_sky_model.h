/*! \file
  Methods to read in and crop a sky model to everything above the horizon.

  @author J.L.B. Line
*/
#pragma once
#include "woden_struct_defs.h"
#include "read_text_skymodel.h"
#include "read_yaml_skymodel.h"

//Somthing WEIRD is up with the way I'm using the documentation package
//'breathe', so having to include the definition value in the documentation
//string to get values to appear (sigh)

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
@brief Given an input char array `str`, check if the string ends with a given
sequence in `suffix`. Good for checking if a filename ends in '.txt' etc.

@details Copied straight up from:
https://stackoverflow.com/questions/744766/how-to-compare-ends-of-strings-in-c

@param[in] str String to be tested
@param[in] suffix String to test whether the end of `str` matches or not

@returns 1 if `str` ends in `suffix`, 0 otherwise
*/
int EndsWith(const char *str, const char *suffix);

/**
@brief Checks whether the string in `str` ends in ".txt"

@details Calls create_sky_model::EndsWith

@param[in] str String to be tested

@returns 1 if `str` ends in `suffix`, 0 otherwise
*/
int EndsWithTxt(const char *str);

/**
@brief Checks whether the string in `str` ends in ".yaml"

@details Calls create_sky_model::EndsWith

@param[in] str String to be tested

@returns 1 if `str` ends in `suffix`, 0 otherwise
*/
int EndsWithYaml(const char *str);

/**
@brief Reads in the skymodel file located at path `srclist` and populates
`raw_srccat` with the results

@details Checks whether the file ends in '.txt' or `.yaml`. If a text file
assumes the file is a WODEN-style format. If '.yaml', assumes file is a
`hyperdrive` style sky model. If `srclist` ends in neither `.txt` or `.yaml`,
returns a 1 and doesn't attempt to read in the file.

Calls either `read_text_skymodel::read_text_skymodel` or
`read_yaml_skymodel::read_yaml_skymodel` as appropriate.

@param[in] srclist Path to sky model to be read in
@param[in] *raw_srccat Pointer to a `source_catalogue_t` struct to hold outputs

@returns 0 if sky model read in successfully, 1 if else
*/
int read_skymodel(const char *srclist, source_catalogue_t *raw_srccat);


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
@brief Calculates the az/za for all ra,dec stored in `components->ra`,
`components->dec` for the first lst in `lsts`. If a single COMPONENT is
below the horizon, update `all_comps_above_horizon` to BELOW, so we know
to crop this SOURCE when sky_crop_type == SOURCE.

@details Calculated az/za are stored in `components->az`, `components->za`.

@param[in] sky_crop_type `e_sky_crop` either SOURCE or COMPONENT, sets how we are cropping the sky
@param[in] *components Pointer to a populated `components_t`, containing information to calculate az/za for
@param[in] num_comps Number of components in `components`
@param[in] lsts All LSTs to be used in simulation
@param[in] latitude latitude of array to use in az/za calculation
@param[in,out] *all_comps_above_horizon Pointer to `e_horizon` type for whether
all COMPONENTs of the SOURCE this COMPONENT belongs to are above the horizon
*/
void horizon_test(e_sky_crop sky_crop_type, components_t * components,
     int num_comps, double *lsts, double latitude,
     e_horizon * all_comps_above_horizon);

/**
@brief Takes sky model `raw_srccat`, and crops either all SOURCEs or
COMPONENTs that are below the horizon (at the initial time step),
returning the cropped sky model. Also calculates $az,za$ for all time steps
for the retained COMPONENTs.

@details The sky model can be made of multiple SOURCEs, each of which can
have any number of COMPONENTS. We can either crop an entire SOURCE if below
horizon, or retain all COMPONENTs that are above the horizon. Say if you put
an all-sky diffuse map into one SOURCE, you should crop by COMPONENT, as some
part of the SOURCE will always be below the horizon.

A complicating factor to the cropping is that each SHAPELET COMPONENT has any
number of associated shapelet coefficients, which must be mapped correctly in
the output cropped sky model `cropped_src`. Furthermore, LIST style flux
density catalogue entries can have any length, so logic must be done when cropping
to grab the correct information from the correct parts of certain arrays.

@param[in] *raw_srccat Pointer to a populated `source_catalogue_t` struct of all input
sky model parameters
@param[in] *lsts Array of local sidereal times for time centriods (radians)
@param[in] *latitudes Latitude of the array for all time centriods (radians).
These change with time when the array is precessed to J2000 for each time step
@param[in] num_time_steps Number of time steps for the simulation
@param[in] sky_crop_type `e_sky_crop` for SOURCE or COMPONENT cropping

@returns `cropped_src`, a `source_t` sky model with only SOURCE/COMPONENTS
above the horizon for simulation

@todo Consider cropping sources iteratively for each time step. If the user runs
a very long simulation, we may end up with large chunks of empty sky. This will
however require far more `erfa` az,za calculations, so for simulations with
millions of sources, more practical to just tell the user to run multiple
shorter simulations.
*/
source_t * crop_sky_model(source_catalogue_t *raw_srccat, double *lsts,
              double *latitudes, int num_time_steps, e_sky_crop sky_crop_type);


/**
@brief Given a pointer to a `source_t` struct, set all the int used
for counting to zero, and arrays for storing sky model to NULL.

@details This gets the model ready to have a number of `realloc`s performed
to store however much information is in the sky model

@param[in] *source Pointer to `source_t` struct

*/
void source_zero_counters_and_null_components(source_t *source);


/**
@brief Given a pointer to a `components` struct, set all arrays for
storing sky model to NULL.

@details

@param[in] *components Pointer to `components_t` struct

*/
void null_component_sky_model_arrays(components_t * components);


// void free_source_catalogue(source_catalogue_t *source_catalogue);
