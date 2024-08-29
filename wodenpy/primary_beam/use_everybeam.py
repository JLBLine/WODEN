import numpy as np
from astropy.coordinates import ITRS, SkyCoord #AltAz, EarthLocation, 
from astropy.time import Time, TimeDelta
import astropy.units as u
import argparse
from wodenpy.use_libwoden.create_woden_struct_classes import Woden_Struct_Classes

import os
##Are we just making online documentation? If so, don't import everybeam
##Installing everybeam is non-trivial, so trying to get readthedocs to install
##it is a waste of time
read_the_docs_build = os.environ.get('READTHEDOCS', None) == 'True'

if read_the_docs_build:
    class EB:
        def __init__(self):
            """
            A fake `everybeam` class so we can build the documentation online
            in ReadTheDocs without installing `everybeam`, which is non-trivial
            """
            self.OSKAR = None
            self.LOFAR = None
            self.load_telescope = None
            self.Telescope = None
    eb = EB()
else:
    import everybeam as eb

##This call is so we can use it as a type annotation
woden_struct_classes = Woden_Struct_Classes()
Source_Catalogue = woden_struct_classes.Source_Catalogue
Woden_Settings = woden_struct_classes.Woden_Settings


def radec_to_xyz(ra : float, dec : float, time : Time):
    """
    Convert RA and Dec ICRS coordinates to ITRS cartesian coordinates.
    
    Taken from the everybeam documentation
    https://everybeam.readthedocs.io/en/latest/tree/demos/lofar-array-factor.html

    Args:
        ra (astropy.coordinates.Angle): Right ascension
        dec (astropy.coordinates.Angle): Declination
        time (float): astropy time instance

    Returns:
        pointing_xyz (ndarray): NumPy array containing the ITRS X, Y and Z coordinates
    """
    
    
    coord = SkyCoord(ra*u.rad, dec*u.rad, frame='icrs')
    coord_itrs = coord.transform_to(ITRS(obstime=time))
    
    return np.asarray(coord_itrs.cartesian.xyz.transpose())

def load_OSKAR_telescope(ms_path : str, response_model = "skala40_wave") -> eb.OSKAR:
    """Load an OSKAR telescope from a measurement set.

    Parameters
    ----------
    ms_path : str
        Path to the measurement set
    response_model : str, optional
        Response model to use, by default "skala40_wave"

    Returns
    -------
    eb.OSKAR
        Telescope object
    """

    ##TODO what does this mean, what is differential beam?
    use_differential_beam = False
    
    print("OSKAR response model", response_model)

    # Load the telescope
    telescope = eb.load_telescope(
        ms_path,
        use_differential_beam=use_differential_beam,
        element_response_model=response_model,
    )
    
    # assert type(telescope) == eb.OSKAR
    if type(telescope) != eb.OSKAR:
        print(f'WARNING: Telescope specified in {ms_path} is not an OSKAR telescope. Proceeding, but you might get nonsense results.')
    
    return telescope


def load_LOFAR_telescope(ms_path : str, response_model = "lobes") -> eb.LOFAR:
    """Load an LOFAR telescope from a measurement set. Settings lifted
    directly from https://everybeam.readthedocs.io/en/latest/tree/demos/lofar-lobes.html

    Parameters
    ----------
    ms_path : str
        Path to the measurement set
    response_model : str, optional
        Response model to use, by default "lobes"

    Returns
    -------
    eb.LOFAR
        Telescope object
    """

    ##TODO what does this mean, what is differential beam?
    use_differential_beam = False

    # Load the telescope
    telescope = eb.load_telescope(ms_path,
                                  use_differential_beam=use_differential_beam,
                                  element_response_model=response_model)
    
    # assert type(telescope) == eb.LOFAR
    if type(telescope) != eb.LOFAR:
        print(f'WARNING: Telescope specified in {ms_path} is not an OSKAR telescope. Proceeding, but you might get nonsense results.')
    
    return telescope


def get_everybeam_norm(phase_itrf : np.ndarray, time : Time, freq : float,
                       telescope : eb.Telescope, station_id = 0) -> np.ndarray:
    """Get a normalisation factor for the X and Y beams from everybeam for
    a given phase centre (ra0, dec), time, frequency, telescope, and station.

    Parameters
    ----------
    phase_itrf : np.ndarray
        XYZ itfr array (as output by `radec_to_xyz`) of beam phase centre
    time : Time
        Astropy Time object of observation
    freq : float
        Frequency of observation in Hz
    telescope : eb.Telescope
        An everybeam telescope object
    station_id : int, optional
        Integer index of station, by default 0

    Returns
    -------
    np.ndarray
        Normalisation factors for the X and Y beams (multiply by this number to apply the norm)
    """
    
    # phase_itrf = radec_to_xyz(ra0, dec0, time)
    
    # Full beam for station 0
    response = telescope.station_response(time.mjd*3600*24, station_id, freq,
                                          phase_itrf, phase_itrf,
                                          rotate=True)
    
    norm_x = 1 / np.abs(response[0,0])
    norm_y = 1 / np.abs(response[1,1])
    # print(station_id, norm_x, norm_y)
    
    return norm_x, norm_y

def run_everybeam(dir_itrf : np.ndarray, phase_itrf : np.ndarray, 
                  time : Time, freq : float, telescope: eb.Telescope,
                  station_id : int = 0, beam_norms : np.ndarray = np.ones(2),
                  reorder_jones : bool = True):
    """For a given direction, phase centre, time, frequency, telescope, station,
    calculate an everybeam jones matrix. Optionally normalise using the given
    beam norms [X_norm, Y_norm]. Explicitly, where
    `response` = [[j00, j01], [j10, j11]], norms are applied as
    
    response[0,:] *= beam_norms[0]
    response[1,:] *= beam_norms[1]
    
    By defauly, reorder the Jones matrix to be
    [-j11, j10, -j01, j00], as WODEN expects X = NS, Y = EW. The negatives
    may be down to a definition of Stokes V. This reordering is the same as
    done in `mwa_hyperbeam`. The reordering is done *after* the beam norms. As
    long as the beam norms have been calculated using `get_everybeam_norm` the
    ordering should be correct.

    Parameters
    ----------
    dir_itrf : np.ndarray
        XYZ itfr array (as output by `radec_to_xyz`) of direction of interest 
    phase_itrf : np.ndarray
        XYZ itfr array (as output by `radec_to_xyz`) of beam phase centre
    time : Time
        Astropy Time object of observation
    freq : float
        Frequency of observation in Hz
    telescope : eb.Telescope
        An everybeam telescope object
    station_id : int, optional
        Integer index of station, by default 0
    beam_norms : np.ndarray, optional
        Normalisation to apply, [X_norm, Y_norm], by default np.ones(2). Outputs are multiplied by these values
    reorder_jones : bool, optional
        If True, reorder the Jones matrix to be [-j11, j10, -j01, j00], which
        is the order expected by WODEN, by default True

    Returns
    -------
    np.ndarray
        2x2 array of complex beam jones matrix [[j00, j01], [j10, j11]]
    """
    
    ##Get the response
    response = telescope.station_response(time.mjd*3600*24, station_id, freq,
                                          dir_itrf, phase_itrf,
                                          rotate=True)
    
    ##normalise the beams using previously calculated norms
    response[0,:] *= beam_norms[0]
    response[1,:] *= beam_norms[1]
    
    if reorder_jones:
        # print("HERE")
        
        # print(response)
        
        ##Might be the case we have to do some different kind of reordering here
        if telescope.__class__ == eb.LOFAR:
            response = np.array([[-response[1,1], response[1,0]], [-response[0,1], response[0,0]]])
            
            ##diff stokes convention??
            response /= np.sqrt(2)
            
            ##rotate by 90 degrees??
            # response = np.dot(np.array([[0, 1], [-1, 0]]), response)
            
        elif telescope.__class__ == eb.OSKAR:
        # else:
        #     # Reorder the Jones matrix to be [-j11, j10, -j01, j00]
        #     # print('Doing this reorder')
            response = np.array([[-response[1,1], response[1,0]], [-response[0,1], response[0,0]]])
            # response = np.array([[response[1,1], response[1,0]], [-response[0,1], response[0,0]]])
        
    return response
