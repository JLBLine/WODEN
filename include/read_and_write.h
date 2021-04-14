/*! \file
  Methods to read in simulation parameters from a .json file,
  read in data from source catalogues, read/generate array layouts, and move
  the array back to J2000.
  @author J.L.B. Line
*/
#pragma once
#include <math.h>
#include <stdint.h>
#include <fitsio.h>
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
 * @brief Convert coords in local topocentric East, North, Height units to
 * 'local' XYZ units.
 *
 * @details Local means Z point north, X points through the equator from the
 * geocenter along the local meridian and Y is East. This is like the absolute
 * system except that zero lon is now the local meridian rather than prim
 * meridian. Latitude is geodetic, in radian. This is what you want for
 * constructing the local antenna positions in a UVFITS antenna table.

 * @param[in] E         East coord of the array positions (metres)
 * @param[in] N         North coord of the array positions (metres)
 * @param[in] H         Height coord of the array positions (metres)
 * @param[in] lat       Latitude of the array (radians)
 * @param[in,out] *X        Pointer to local X coord to be filled (metres)
 * @param[in,out] *Y        Pointer to local Y coord to be filled (metres)
 * @param[in,out] *Z        Pointer to local Z coord to be filled (metres)
 */
void RTS_ENH2XYZ_local(float E, float N, float H, float lat,
                       float *X, float *Y, float *Z);

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
documentation at
@todo Link the online documentation when there is a link

@param[in] *filename Path to a WODEN-style sky model
@return A pointer to a populated `source_catalogue_t` struct
 */
source_catalogue_t * read_source_catalogue(const char *filename);


/**
 * @brief Takes a path to .json WODEN parameter file, and populates a
 * `woden_settings_t` struct with the contents of `filename`.
 *
 * @details For what can be included in the .json file, see the documentation
 * for
 * @todo Work out how to link the print_help function here
 *
 * @param[in] *filename Path to a WODEN *.json settings file
 * @return A pointer to a populated `woden_settings_t` struct
 */
woden_settings_t * read_json_settings(const char *filename);

/**
 * @brief Populates an `array_layout_t` struct with useful antenna coords and
 * baseline lengths for the array layout specified in `*woden_settings`
 *
 * @details Explicitly, tries to read the array layout stored at
 * `woden_settings->array_layout_file_path`, storeds the E,N,H coords,
 * converts to and stores X,Y,Z coords, and calculates the baseline lengths
 * in the X,Y,Z frame
 *
 * @param[in] *woden_settings Pointer to a `woden_settings_t` struct
 * @return A pointer to a populated `array_layout_t` struct
 */
array_layout_t * calc_XYZ_diffs(woden_settings_t *woden_settings);

/**
 * @brief Rotates the array coordinates in `x,y,z` from the current date to
 * the J2000 frame, returning the precessed coordinates in `xp,yp,zp`

 * @details wtf you do
 *
 * @param[in] rmat A 3D rotation matrix
 * @param[in] x X coord of the antenna in the current epoch
 * @param[in] y Y coord of the antenna in the current epoch
 * @param[in] z Z coord of the antenna in the current epoch
 * @param[in] lmst Current local mean sidereal time (radians)
 * @param[in,out] xp X coord of the antenna in J2000 epoch
 * @param[in,out] yp Y coord of the antenna in J2000 epoch
 * @param[in,out] zp Z coord of the antenna in J2000 epoch
 * @param[in] lmst2000 J2000 local mean sidereal time (radians)
 */
void RTS_precXYZ(double rmat[3][3], double x, double y, double z, double lmst,
         double *xp, double *yp, double *zp, double lmst2000);

/**
* @brief Performs calculations to rotate the array coordinates from the
* current date back to their positions in J2000
*
* @details For any simulation on a date other than 01/01/2000, a J2000 sky model
* must be precessed to the current date, or the sources are in an incorrect
* location. Rather than precess every source however, we can just rotate the
* array itself back to J2000. I think this ignores the intrinsic velocities
* of each source but covers the Earth based precessions and what not.
* Explicitly, these three arrays in `array_layout` are updated:
* array_layout->ant_X,
* array_layout->ant_Y,
* array_layout->ant_Z
*
* Many, many calls to `pal` are made.
*
* @param[in,out] *array_layout An `array_layout_t` struct containing the array
* layout
* @param[in] *woden_settings A `woden_settings_t` struct containing
* observational settings
*/
void RTS_PrecessXYZtoJ2000( array_layout_t *array_layout,
                       woden_settings_t *woden_settings);
