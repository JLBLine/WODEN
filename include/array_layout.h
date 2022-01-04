/*! \file
  Methods to read read/generate array layouts, and move the array back to J2000.
  @author J.L.B. Line
*/
#pragma once
#include <math.h>
#include <stdint.h>
#include "woden_struct_defs.h"

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
void RTS_ENH2XYZ_local(double E, double N, double H,
                       double lat,
                       double *X, double *Y, double *Z);

/**
@brief Populates an `array_layout_t` struct with useful antenna coords and
baseline lengths for the array layout specified in `*woden_settings`

@details Explicitly, tries to read the array layout stored at
`woden_settings->array_layout_file_path`, stores the E,N,H coords,
converts to and stores X,Y,Z coords, and calculates the baseline lengths
in the X,Y,Z frame. If `do_precession=1`, it will rotate the X,Y,Z to J2000
(based off of mjd in `woden_settings`). This should be equivalent to rotating
the sky to the present day to account for precession. Assumes the sky model is
in J2000 coords.

@param[in] *woden_settings Pointer to a `woden_settings_t` struct
@param[in] do_precession Whether to precess back to J2000 or not (0 False, 1 True)
@return A pointer to a populated `array_layout_t` struct
 */
array_layout_t * calc_XYZ_diffs(woden_settings_t *woden_settings,
                                int do_precession);

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
@brief Performs calculations to rotate the array coordinates from the
current date back to their positions in J2000

@details For any simulation on a date other than 01/01/2000, a J2000 sky model
must be precessed to the current date, or the sources are in an incorrect
location. Rather than precess every source however, we can just rotate the
array itself back to J2000. I think this ignores the intrinsic velocities
of each source but covers the Earth based precessions and what not.
Explicitly, these three arrays in `array_layout` are updated:
`array_layout->ant_X`,
`array_layout->ant_Y`,
`array_layout->ant_Z`, along with `woden_settings->lst_base`, as moving the
array effectively changes the LST. Many, many calls to `pal` are made.

@param[in,out] *array_layout An `array_layout_t` struct containing the array
layout
@param[in] *woden_settings A `woden_settings_t` struct containing
observational settings
*/
void RTS_PrecessXYZtoJ2000( array_layout_t *array_layout,
                       woden_settings_t *woden_settings);
