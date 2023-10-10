import numpy as np
import erfa
import palpy
from typing import Tuple


def RTS_Precess_LST_Lat_to_J2000(lst_current : float, latitude_current : float,
                                 mjd : float) -> Tuple[float, float]:
    """
    Convert the current local sidereal time and latitude to the J2000 mean system.

    Parameters
    ------------
    lst_current : float
        The current local sidereal time in radians.
    latitude_current : float
        The current latitude in radians.
    mjd : float
        The current modified Julian date.

    Returns
    ------------
    lst_J2000, latitude_J2000 : Tuple[float, float]
        A tuple containing the local sidereal time and latitude in the J2000 mean system, both in radians.
    """

    # Calculate a rotation matrix that accounts for precession and nutation
    # between the current modified julian date (mjd) and J2000
    #  palPrenut calls:
    #   - palPrec( 2000.0, palEpj(mjd), rmatp ); // form precession matrix: v_mean(mjd epoch) = rmatp * v_mean(J2000)
    #   - palNut( mjd, rmatn );                  // form nutation matrix: v_true(mjd epoch) = rmatn * v_mean(mjd epoch)
    #   - palDmxm( rmatn, rmatp, rmatpn );       // Combine the matrices:  pn = n x p
    
    rmatpn = palpy.prenut(2000, mjd)
    J2000_transformation = np.transpose(rmatpn)
    
    # /**
    # ****************************************************************************
    # * Change the various coordinates to the J2000 mean system
    # ****************************************************************************
    # * palDcs2c   - convert the apparent direction to direction cosines
    # * palDmxv    - perform the 3-d forward unitary transformation: v2 = tmatpn * v1
    # * palDcc2s   - convert cartesian coordinates back to spherical coordinates (i.e. zenith in the J2000 mean system).
    # * palDranrm  - normalize into range 0-2 pi.
    # */

    # // Change the coordinates of the initial zenith
    v1 = palpy.dcs2c(lst_current, latitude_current)
    v2 = palpy.dmxv(J2000_transformation, v1)
    lst_J2000, latitude_J2000 = palpy.dcc2s(v2)
    lst_J2000 = palpy.dranrm(lst_J2000)

    return lst_J2000, latitude_J2000