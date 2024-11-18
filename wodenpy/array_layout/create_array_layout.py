import numpy as np
import sys
import os
from typing import Union
import palpy
from ctypes import Structure
import argparse
from typing import Tuple
from erfa import gd2gc

##Constants
R2D = 180.0 / np.pi
D2R = np.pi / 180.0
MWA_LAT = -26.703319405555554
MWA_LONG = 116.67081523611111
MWA_HEIGHT = 377.827
VELC = 299792458.0
SOLAR2SIDEREAL = 1.00274
DS2R = 7.2722052166430399038487115353692196393452995355905e-5

from wodenpy.use_libwoden.create_woden_struct_classes import Woden_Struct_Classes
from wodenpy.use_libwoden.array_layout_struct import Array_Layout, setup_array_layout

##This call is so we can use it as a type annotation
woden_struct_classes = Woden_Struct_Classes()
Woden_Settings = woden_struct_classes.Woden_Settings


def convert_ecef_to_enh(ecef_X : np.ndarray, ecef_Y : np.ndarray, ecef_Z : np.ndarray,
                        lon : float, lat : float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    The array coords in a measurement set seem to be Earth Centred Earth Fixed
    coords
    
    I got the maths to go from ECEF to enh from this site:
    https://gssc.esa.int/navipedia/index.php/Transformations_between_ECEF_and_ENU_coordinates

    Parameters
    ----------
    ecef_X : np.ndarray
        ECEF X coordinates (metres)
    ecef_Y : np.ndarray
        ECEF Y coordinates (metres)
    ecef_Z : np.ndarray
        ECEF Z coordinates (metres)
    lon : float
        Longitude of the array (radians)
    lat : float
        Latitude of the array (radians

    Returns
    -------
    Tuple[np.ndarray, np.ndarray, np.ndarray]
        The east, north, and height coordinates
    """
                        
    ##Calculate the ecef xyz of the array position using the geodectic to
    ##geocentric function
    arrX, arrY, arrZ = gd2gc(1, lon, lat, 0)
    
    ##subtract it from the earth centred coords, as we want enh coords
    ##centred at the array
    X = ecef_X - arrX
    Y = ecef_Y - arrY
    Z = ecef_Z - arrZ
    
    ##convert xyz to enh
    east = -np.sin(lon)*X + np.cos(lon)*Y
    north = -np.cos(lon)*np.sin(lat)*X - np.sin(lon)*np.sin(lat)*Y + np.cos(lat)*Z
    height = np.cos(lon)*np.cos(lat)*X + np.sin(lon)*np.cos(lat)*Y + np.sin(lat)*Z
    
    return east, north, height

def enh2xyz(east : float, north : float, height : float, latitude : float) -> Tuple[float, float, float]:
    """
    Takes local east, north, height coords for a given latitude (radians)
    and returns local X,Y,Z coords to put in the uvfits antenna table

    Parameters
    ----------
    east : float
        Local east coorindate (metres)
    north : float
        Local north coorindate (metres)
    height : float
        Local height coorindate (metres)
    latitude : float
        Latitude of the array - defaults to MWA location (radians)

    Returns
    -------
    X : float
        Local X antenna location
    Y : float
        Local Y antenna location
    Z : float
        Local Z antenna location
    """

    sl = np.sin(latitude)
    cl = np.cos(latitude)
    X = -north*sl + height*cl
    Y = east
    Z = north*cl + height*sl

    return X,Y,Z

def RTS_precXYZ(rmat : np.ndarray, x : float, y : float, z : float, lmst : float,
                lmst2000 : float) -> Tuple[float, float, float]:
    """RTS magic for precessing local X,Y,Z from the current time frame back to
    J2000

    Parameters
    ----------
    rmat : np.ndarray
        2D rotation matrix as output by `palpy.prenut`
    x : float
        Local X coord of antenna (tile)
    y : float
        Local Y coord of antenna (tile)
    z : float
        Local X coord of antenna (tile)
    lmst : float
        LST of the array in the current frame
    lmst2000 : float
        LST of the array in the J2000 frame

    Returns
    -------
    Tuple[float, float, float]
        The precessed X,Y,Z coords
    """


    sep = np.sin(lmst)
    cep = np.cos(lmst)
    s2000 = np.sin(lmst2000)
    c2000 = np.cos(lmst2000)

    ##/* rotate to frame with x axis at zero RA */
    xpr = cep*x - sep*y
    ypr = sep*x + cep*y
    zpr = z

    ##/* apply the rotation matrix to account for precession/nutation */
    xpr2 = (rmat[0][0])*xpr + (rmat[0][1])*ypr + (rmat[0][2])*zpr
    ypr2 = (rmat[1][0])*xpr + (rmat[1][1])*ypr + (rmat[1][2])*zpr
    zpr2 = (rmat[2][0])*xpr + (rmat[2][1])*ypr + (rmat[2][2])*zpr

    ##/* rotate back to frame with xp pointing out at lmst2000 */
    xp = c2000*xpr2 + s2000*ypr2
    yp = -s2000*xpr2 + c2000*ypr2
    zp = zpr2

    return xp, yp, zp

def RTS_PrecessXYZtoJ2000( array_layout : Array_Layout,
        woden_settings : Woden_Settings) -> Array_Layout:
    """Given the populated `array_layout` and settings in `woden_settings`, use
    RTS functions to precess the array back to J2000, to account for the
    skymodel being in J2000.

    Parameters
    ----------
    array_layout : Array_Layout
        Populated array layout in the simulation epoch
    woden_settings : Woden_Settings
        Populated settings, where the `woden_settings.lsts` should have been
        precessed back to J2000 already.

    Returns
    -------
    Array_Layout
        The updated array layout with `array_layout.ant_X`, `array_layout.ant_Y`,
        and `array_layout.ant_Z` precessed back to J2000."""
        
    ##Rotate the array positions for each time step - they have different
    ##mjd dates and current epoch lsts so yield different XYZ over time
    for time_step in range(woden_settings.num_time_steps):

        lst_J2000 = woden_settings.lsts[time_step]
        mjd_current = woden_settings.mjds[time_step]

        ##Add on the angle accrued by current time step to the base LST
        lst_current = woden_settings.lst_obs_epoch_base + time_step*woden_settings.time_res*SOLAR2SIDEREAL*DS2R
        ##Add half a time_res so we are sampling centre of each time step
        lst_current += 0.5*woden_settings.time_res*SOLAR2SIDEREAL*DS2R

        rmatpn = palpy.prenut(2000, mjd_current)
        J2000_transformation = np.transpose(rmatpn)

        for st in range(array_layout.num_tiles):

            ##Where do this XYZ for this time step start
            st_offset = time_step*array_layout.num_tiles

            ##Calculate antennas positions in the J2000 frame
            X_epoch = array_layout.ant_X[st_offset + st]
            Y_epoch = array_layout.ant_Y[st_offset + st]
            Z_epoch = array_layout.ant_Z[st_offset + st]
            X_prec, Y_prec, Z_prec = RTS_precXYZ(J2000_transformation,
                                                 X_epoch, Y_epoch, Z_epoch,
                                                 lst_current, lst_J2000 )

            ##Update stored coordinates
            array_layout.ant_X[st_offset + st] = X_prec
            array_layout.ant_Y[st_offset + st] = Y_prec
            array_layout.ant_Z[st_offset + st] = Z_prec

    return array_layout



def calc_XYZ_diffs(woden_settings : Woden_Settings,
                   args : argparse.Namespace) -> Array_Layout:
    """
    Populates an Array_Layout class with the instrument layout, given the command
    line arguments. Calculates the differences in X, Y, and Z coordinates
    between all pairs of antennas in the array.

    Parameters
    ----------
    woden_settings: Woden_Settings
        An populated Woden_Settings instance
    args: argparse.Namespace
        The command line arguments checked by `wodenpy.wodenpy_setup.check_args`

    Returns
    -------
    array_layout
    - An instance of the Structure class containing the X, Y, and Z differences between all pairs of antennas in the array at each time step.
    """
    ##sets up the ctype array_layout struct equiv
    array_layout = setup_array_layout(woden_settings, args)

    # ##When precessing the array back to J2000, we apply a different precession
    # ##based on mjd of each time step. So store enough copies of XYZ for
    # ##all time steps
    for time_step in range(woden_settings.num_time_steps):
        for i in range(array_layout.num_tiles):

            ##Where do this XYZ for this time step start
            st_offset = time_step*array_layout.num_tiles

            ##Convert to local X,Y,Z
            x, y, z = enh2xyz(args.east[i], args.north[i],
                              args.height[i],
                              woden_settings.latitude_obs_epoch_base)
            
            array_layout.ant_X[st_offset + i] = x
            array_layout.ant_Y[st_offset + i] = y
            array_layout.ant_Z[st_offset + i] = z

    if woden_settings.do_precession:
        print("We are precessing the array")
        array_layout = RTS_PrecessXYZtoJ2000(array_layout, woden_settings)

    for time_step in range(woden_settings.num_time_steps):

        ##Where do these XYZ for this time step start
        ant_off = array_layout.num_tiles*time_step

        ##Where do these baselines for this time step start
        base_off = array_layout.num_baselines*time_step

        baseline_ind = 0
        for ant1 in range(array_layout.num_tiles - 1):
            for ant2 in range(ant1 + 1, array_layout.num_tiles):
                array_layout.X_diff_metres[base_off + baseline_ind] = array_layout.ant_X[ant_off + ant1] - array_layout.ant_X[ant_off + ant2]
                array_layout.Y_diff_metres[base_off + baseline_ind] = array_layout.ant_Y[ant_off + ant1] - array_layout.ant_Y[ant_off + ant2]
                array_layout.Z_diff_metres[base_off + baseline_ind] = array_layout.ant_Z[ant_off + ant1] - array_layout.ant_Z[ant_off + ant2]

                baseline_ind += 1

    return array_layout