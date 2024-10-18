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
            self.MWA = None
            self.load_telescope = None
            self.Telescope = None
    eb = EB()
else:
    import everybeam as eb

##This call is so we can use it as a type annotation
woden_struct_classes = Woden_Struct_Classes()
Source_Catalogue = woden_struct_classes.Source_Catalogue
Woden_Settings = woden_struct_classes.Woden_Settings


USE_DIFFERENTIAL_BEAM = True

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

    print("OSKAR response model", response_model)

    # Load the telescope
    telescope = eb.load_telescope(
        ms_path,
        use_differential_beam=USE_DIFFERENTIAL_BEAM,
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

    # Load the telescope
    telescope = eb.load_telescope(ms_path,
                                  use_differential_beam=False,
                                  element_response_model=response_model)
    
    # assert type(telescope) == eb.LOFAR
    if type(telescope) != eb.LOFAR:
        print(f'WARNING: Telescope specified in {ms_path} is not an OSKAR telescope. Proceeding, but you might get nonsense results.')
    
    return telescope


def load_MWA_telescope(ms_path : str, coeff_path : str) -> eb.MWA:
    """Load an MWA telescope from a measurement set.

    Parameters
    ----------
    ms_path : str
        Path to the measurement set
    response_model : str, optional
        Response model to use, by default "lobes"

    Returns
    -------
    eb.MWA
        Telescope object
    """

    # Load the telescope
    telescope = eb.load_telescope(ms_path,
                                  use_differential_beam=USE_DIFFERENTIAL_BEAM,
                                  coeff_path=coeff_path)
                                #   element_response_model=response_model)
    
    # assert type(telescope) == eb.MWA
    if type(telescope) != eb.MWA:
        print(f'WARNING: Telescope specified in {ms_path} is not an MWA telescope. Proceeding, but you might get nonsense results.')
    
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
    
    response = telescope.station_response(time.mjd*3600*24, station_id, freq,
                                          phase_itrf, phase_itrf,
                                          rotate=True)
    
    norm_x = 1 / np.abs(response[0,0])
    norm_y = 1 / np.abs(response[1,1])
    
    return norm_x, norm_y

def run_everybeam(dir_itrf : np.ndarray, phase_itrf : np.ndarray, 
                  time : Time, freq : float, telescope: eb.Telescope,
                  station_id : int = 0, beam_norms : np.ndarray = np.ones(2),
                  reorder_jones : bool = True,
                  ra : float = False, dec : float = False,
                  parallactic_angle : float = 0,
                  para_angle_offset : float = -np.pi/2) -> np.ndarray:
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
        If True, reorder the Jones matrix to be [j11, j10, j01, j00], which
        is the order expected by WODEN, by default True
    ra : float, optional
        Right ascension of direction of interest (radians), by default False. Needed for MWA beam
    dec : float, optional
        Declination of direction of interest (radians), by default False. Needed for MWA beam
    parallactic_angle : float, optional
        Rotate by this parallactic angle (radians); by default 0
    para_angle_offset : float, optional
        Offset to apply to parallactic angle (radians); by default -np.pi/2

    Returns
    -------
    np.ndarray
        2x2 array of complex beam jones matrix [[j00, j01], [j10, j11]]
    """
    
    ##The MWA everybeam function takes ra, dec, not dir_itrf
    
    if type(telescope) == eb.MWA:
        ##Get the response
        response = telescope.station_response(time.mjd*3600*24, station_id, freq,
                                            ra, dec)
    else:
        
        parallactic_angle = 0
    
        ##Get the response
        response = telescope.station_response(time.mjd*3600*24, station_id, freq,
                                              dir_itrf, phase_itrf,
                                              rotate=True)
    ##normalise the beams using previously calculated norms
    response[0,:] *= beam_norms[0]
    response[1,:] *= beam_norms[1]
    
    if parallactic_angle:
        
        cospa = np.cos(parallactic_angle + para_angle_offset)
        sinpa = np.sin(parallactic_angle + para_angle_offset)
        
        rotated_response = np.zeros_like(response)
        
        rotated_response[0,0] = response[0,0]*cospa - response[0,1]*sinpa
        rotated_response[0,1] = response[0,0]*sinpa + response[0,1]*cospa
        rotated_response[1,0] = response[1,0]*cospa - response[1,1]*sinpa
        rotated_response[1,1] = response[1,0]*sinpa + response[1,1]*cospa
        
        response = rotated_response
    
    if reorder_jones:
        response = np.array([[response[1,1], -response[1,0]],
                             [response[0,1], -response[0,0]]])
        
    ##diff stokes convention??
    # if type(telescope) == eb.LOFAR:
    #     response /= np.sqrt(2)
    
    return response
