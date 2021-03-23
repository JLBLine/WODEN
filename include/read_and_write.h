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
// #include "FEE_primary_beam.h"
#include "woden_struct_defs.h"

// #define CS_LEN_SRC_NAME 16

#define SRC_KEY         "SOURCE"
#define SRC_END         "ENDSOURCE"
#define COMP_KEY        "COMPONENT"
#define COMP_END        "ENDCOMPONENT"
#define FREQ_KEY        "FREQ"
#define LINEAR_KEY      "LINEAR"
#define POINT_KEY       "POINT"
#define GAUSSIAN_KEY    "GAUSSIAN"
#define GPARAMS_KEY     "GPARAMS"
#define SHAPELET_KEY    "SHAPELET"
#define SPARAMS_KEY     "SPARAMS"
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
 * @brief Takes a path to WODEN-style sky model and populates a
 * `source_catalogue_t` struct with the contents of `filename`.
 * Adapted from the RTS (Mitchell et al 2008). All credit to the original authors.
 * RTS can be found here https://github.com/ICRAR/mwa-RTS.git
 *
 * @details The WODEN sourcelist at `filename` should contain a number of
 * sources, with the basic structure:\n
 *    SOURCE source_name P 1 G 0 S 0 0 \n
 *    COMPONENT POINT 4.0 -27.0 \n
 *    LINEAR 1.8e+08 10.0 0 0 0 -0.8\n
 *    ENDCOMPONENT \n
 *    ENDSOURCE
 *
 * For a more detailed explanation of the sky model, please see the
 * documentation at
 * @todo Link the online documentation when there is a link
 *
 * @param[in] *filename Path to a WODEN-style sky model
 * @return A pointer to a populated `source_catalogue_t` struct
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
