import numpy as np
import sys
import os
from typing import Union
import palpy
from ctypes import Structure
import argparse

##Constants
R2D = 180.0 / np.pi
D2R = np.pi / 180.0
MWA_LAT = -26.703319405555554
MWA_LONG = 116.67081523611111
MWA_HEIGHT = 377.827
VELC = 299792458.0
SOLAR2SIDEREAL = 1.00274
DS2R = 7.2722052166430399038487115353692196393452995355905e-5

# ##If we are performing a ctest, this check means we use the code we are
# ##testing and NOT what has been pip or conda installed
# try:
#     testdir = os.environ['CMAKE_CURRENT_SOURCE_DIR']
#     sys.path.append('{:s}/../../../wodenpy/use_libwoden'.format(testdir))
#     from woden_settings import Woden_Settings_Float, Woden_Settings_Double
#     from array_layout_struct import Array_Layout, setup_array_layout
    
# except KeyError:
from wodenpy.use_libwoden.woden_settings import Woden_Settings_Float, Woden_Settings_Double
from wodenpy.use_libwoden.array_layout_struct import Array_Layout, setup_array_layout

def enh2xyz(east, north, height, latitude):
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
                lmst2000 : float):

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
        woden_settings : Union[Woden_Settings_Float, Woden_Settings_Double]):
    
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



def calc_XYZ_diffs(woden_settings : Union[Woden_Settings_Float, Woden_Settings_Double],
                   args : argparse.Namespace) -> Structure:

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