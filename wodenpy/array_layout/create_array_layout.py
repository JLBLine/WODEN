"""Functions to create the array layout for the simulation, precess it, 
and convert it to the correct format for the C functions."""
import numpy as np
import sys
import os
from typing import Union
import palpy
from ctypes import Structure
import argparse
from typing import Tuple
from erfa import gd2gc
from ctypes import POINTER, c_double, c_int
from wodenpy.wodenpy_setup.woden_logger import simple_logger
from logging import Logger

##Constants
R2D = 180.0 / np.pi
D2R = np.pi / 180.0
MWA_LAT = -26.703319405555554
MWA_LONG = 116.67081523611111
MWA_HEIGHT = 377.827
VELC = 299792458.0
SOLAR2SIDEREAL = 1.00274
DS2R = 7.2722052166430399038487115353692196393452995355905e-5

from wodenpy.use_libwoden.woden_settings import Woden_Settings_Python
from wodenpy.use_libwoden.array_layout_struct import Array_Layout_Ctypes

##This call is so we can use it as a type annotation
# woden_struct_classes = Woden_Struct_Classes()
# Woden_Settings = woden_struct_classes.Woden_Settings

class Array_Layout_Python(object):
    """
    A class to represent the layout of an array in Python.

    Attributes
    ----------
    ant_X : None
        Antenna X coordinates (metres).
    ant_Y : None
        Antenna Y coordinates (metres).
    ant_Z : None
        Antenna Z coordinates (metres).
    X_diff_metres : None
        The difference in X coordinates (metres).
    Y_diff_metres : None
        The difference in Y coordinates (metres).
    Z_diff_metres : None
        The difference in Z coordinates (metres).
    ant_east : None
        Antenna east coordinates (metres)
    ant_north : None
        Antenna north coordinates (metres)
    ant_height : None
        Antenna height coordinates (metres)
    latitude : None
        The latitude of the array (radians).
    num_baselines : None
        The number of baselines.
    num_tiles : None
        The number of tiles.
    lst_base : None
        The base of the local sidereal time.
    """
    def __init__(self):
        self.ant_X = None
        self.ant_Y = None
        self.ant_Z = None
        self.X_diff_metres = None
        self.Y_diff_metres = None
        self.Z_diff_metres = None
        self.ant_east = None
        self.ant_north = None
        self.ant_height = None
        self.latitude = None
        self.num_baselines = None
        self.num_tiles = None
        self.lst_base = None


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


def convert_enh_to_ecef(east : np.ndarray, north : np.ndarray, height : np.ndarray,
                        lon : float, lat : float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Convert from array centred east, north, heigh (ENH) coordinates to
    Earth Centred, Earth Fixed (ECEF) coordinates.
    
    I got the maths from this site:
    https://gssc.esa.int/navipedia/index.php/Transformations_between_ECEF_and_ENU_coordinates

    Parameters
    ----------
    east : np.ndarray
        Local east (metres)
    north : np.ndarray
        Local north (metres)
    height : np.ndarray
        Local height (metres)
    lon : float
        Longitude of the array centre (radians)
    lat : float
        Latitude of the array centre (radians)

    Returns
    -------
    Tuple[np.ndarray, np.ndarray, np.ndarray]
        The ECEF X, Y, Z coordinates
    """
                        
    ##Calculate the ecef xyz of the array position using the geodectic to
    ##geocentric function
    arrX, arrY, arrZ = gd2gc(1, lon, lat, 0)
    

    
    ##convert enh to local earth centred coords (but not earth fixed??)
    
    X = -np.sin(lon)*east - np.sin(lat)*np.cos(lon)*north + np.cos(lat)*np.cos(lon)*height
    Y = np.cos(lon)*east - np.sin(lat) * np.sin(lon)*north + np.cos(lat)*np.sin(lon)*height
    Z = np.cos(lat)*north + np.sin(lat)*height
    
    ##add the earth centred coords of the array centre to get earth-fixed
    ecef_X = X + arrX
    ecef_Y = Y + arrY
    ecef_Z = Z + arrZ
    
    return ecef_X, ecef_Y, ecef_Z



def rotate_local_XYZ(X, Y, Z, curr_lat, curr_lon, new_lat, new_lon):
    
    east = -np.sin(curr_lon)*X + np.cos(curr_lon)*Y
    north = -np.cos(curr_lon)*np.sin(curr_lat)*X - np.sin(curr_lon)*np.sin(curr_lat)*Y + np.cos(curr_lat)*Z
    height = np.cos(curr_lon)*np.cos(curr_lat)*X + np.sin(curr_lon)*np.cos(curr_lat)*Y + np.sin(curr_lat)*Z
    
    X = -np.sin(new_lon)*east - np.sin(new_lat)*np.cos(new_lon)*north + np.cos(new_lat)*np.cos(new_lon)*height
    Y = np.cos(new_lon)*east - np.sin(new_lat) * np.sin(new_lon)*north + np.cos(new_lat)*np.sin(new_lon)*height
    Z = np.cos(new_lat)*north + np.sin(new_lat)*height
    
    return X, Y, Z

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

def RTS_PrecessXYZtoJ2000( array_layout : Array_Layout_Python,
        woden_settings : Woden_Settings_Python) -> Array_Layout_Python:
    """Given the populated `array_layout` and settings in `woden_settings`, use
    RTS functions to precess the array back to J2000, to account for the
    skymodel being in J2000.

    Parameters
    ----------
    array_layout : Array_Layout_Python
        Populated array layout in the simulation epoch
    woden_settings : Woden_Settings_Python
        Populated settings, where the `woden_settings.lsts` should have been
        precessed back to J2000 already.

    Returns
    -------
    Array_Layout_Ctypes
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

def setup_array_layout_python(woden_settings_python : Woden_Settings_Python,
                              args : argparse.Namespace) -> Array_Layout_Python: #type: ignore
    """Given the populated Woden_Settings_Python class and the command line arguments,
    set up the `Array_Layout_Python` class, and fill it with the correct values,
    as well as creating the correct size arrays:
     - array_layout.ant_X
     - array_layout.ant_Y
     - array_layout.ant_Z
     - array_layout.X_diff_metres
     - array_layout.Y_diff_metres
     - array_layout.Z_diff_metres

    Parameters
    ----------
    woden_settings_python : Woden_Settings_Python
        Populated Woden_Settings_Python class.
    args : argparse.Namespace
        Input arguments from the command line that have been checked using
        `wodenpy.use_libwoden.check_args.check_args`.

    Returns
    -------
    array_layout : Array_Layout_Python
        Initialised Array_Layout_Python class.
    """

    array_layout = Array_Layout_Python()

    array_layout.num_tiles = args.num_antennas
    woden_settings_python.num_ants = args.num_antennas

    array_layout.latitude = woden_settings_python.latitude
    array_layout.num_baselines = int((array_layout.num_tiles*(array_layout.num_tiles-1)) / 2)
    woden_settings_python.num_baselines = array_layout.num_baselines

    num_cross = woden_settings_python.num_baselines * woden_settings_python.num_time_steps * woden_settings_python.num_freqs

    woden_settings_python.num_cross = num_cross

    ##If user asked for auto-correlations, set this number to total autos to work out
    ##Otherwise stick it to zero
    if woden_settings_python.do_autos:
        num_autos = woden_settings_python.num_ants * woden_settings_python.num_time_steps * woden_settings_python.num_freqs
    else:
        num_autos = 0
    
    woden_settings_python.num_autos = num_autos
    woden_settings_python.num_visis = num_cross + num_autos
    
    array_layout.ant_X = np.empty(woden_settings_python.num_time_steps*args.num_antennas, dtype=np.float64)
    array_layout.ant_Y = np.empty(woden_settings_python.num_time_steps*args.num_antennas, dtype=np.float64)
    array_layout.ant_Z = np.empty(woden_settings_python.num_time_steps*args.num_antennas, dtype=np.float64)

    array_layout.X_diff_metres = np.empty(woden_settings_python.num_time_steps*array_layout.num_baselines, dtype=np.float64)
    array_layout.Y_diff_metres = np.empty(woden_settings_python.num_time_steps*array_layout.num_baselines, dtype=np.float64)
    array_layout.Z_diff_metres = np.empty(woden_settings_python.num_time_steps*array_layout.num_baselines, dtype=np.float64)
    
    return array_layout


def calc_XYZ_diffs(woden_settings_python : Woden_Settings_Python,
                   args : argparse.Namespace, 
                   logger : Logger = False) -> Array_Layout_Python:
    """
    Populates an Array_Layout_Ctypes class with the instrument layout, given the command
    line arguments. Calculates the differences in X, Y, and Z coordinates
    between all pairs of antennas in the array.

    Parameters
    ----------
    woden_settings_python: Woden_Settings_Python
        An populated Woden_Settings instance
    args: argparse.Namespace
        The command line arguments checked by `wodenpy.wodenpy_setup.check_args`
    logger: Logger
        A logger instance to log messages

    Returns
    -------
    array_layout
    - An instance of the Structure class containing the X, Y, and Z differences between all pairs of antennas in the array at each time step.
    """
    
    if logger == False:
        logger = simple_logger()

    array_layout = setup_array_layout_python(woden_settings_python, args)

    # ##When precessing the array back to J2000, we apply a different precession
    # ##based on mjd of each time step. So store enough copies of XYZ for
    # ##all time steps
    for time_step in range(woden_settings_python.num_time_steps):
        for i in range(array_layout.num_tiles):

            ##Where do this XYZ for this time step start
            st_offset = time_step*array_layout.num_tiles

            ##Convert to local X,Y,Z
            x, y, z = enh2xyz(args.east[i], args.north[i], args.height[i],
                              woden_settings_python.latitude_obs_epoch_base)
            
            array_layout.ant_X[st_offset + i] = x
            array_layout.ant_Y[st_offset + i] = y
            array_layout.ant_Z[st_offset + i] = z

    if woden_settings_python.do_precession:
        logger.info("Precessing array layout to J2000")
        array_layout = RTS_PrecessXYZtoJ2000(array_layout, woden_settings_python)
    else:
        logger.info("Not precessing the array layout to J2000")

    for time_step in range(woden_settings_python.num_time_steps):

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

def convert_array_layout_to_ctypes(array_layout_python : Array_Layout_Python,
                                   array_layout_ctypes : Array_Layout_Ctypes) -> Array_Layout_Ctypes:
    """Just convert things we need for C/GPU functions, not every field"""
    
    
    array_layout_ctypes.X_diff_metres = array_layout_python.X_diff_metres.ctypes.data_as(POINTER(c_double))
    array_layout_ctypes.Y_diff_metres = array_layout_python.Y_diff_metres.ctypes.data_as(POINTER(c_double))
    array_layout_ctypes.Z_diff_metres = array_layout_python.Z_diff_metres.ctypes.data_as(POINTER(c_double))
    
    array_layout_ctypes.latitude = array_layout_python.latitude
    array_layout_ctypes.num_baselines = array_layout_python.num_baselines
    array_layout_ctypes.num_tiles = array_layout_python.num_tiles
    # array_layout_ctypes.lst_base = array_layout_python.lst_base
    
    return array_layout_ctypes