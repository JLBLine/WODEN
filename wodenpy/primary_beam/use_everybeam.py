"""Wrappers around the EveryBeam library to calculate primary beam values.
The core functionality is in the C++ library `libuse_everybeam.so`, which in
turn contains wrappers about the EveryBeam library. Functions here call the
C++ library to calculate the primary beam values, and return them to Python
usable data types, via `ctypes`.

See https://everybeam.readthedocs.io/en/latest/index.html for more information
on the EveryBeam library."""

from time import time
import numpy as np
from astropy.coordinates import ITRS, SkyCoord, AltAz, EarthLocation
from astropy.time import Time, TimeDelta
import astropy.units as u
import argparse
from wodenpy.use_libwoden.create_woden_struct_classes import Woden_Struct_Classes
import erfa
from typing import Union, Tuple
import concurrent.futures
from line_profiler import profile
import os
import astropy
import wodenpy
import importlib_resources
from wodenpy.use_libwoden.skymodel_structs import c_double_complex
from ctypes import c_char_p, c_int, c_double, POINTER, c_bool
import ctypes
from wodenpy.wodenpy_setup.woden_logger import simple_logger
from logging import Logger
from multiprocessing import Process, Queue
from wodenpy.use_libwoden.beam_settings import BeamTypes
import ctypes
from wodenpy.array_layout.create_array_layout import convert_ecef_to_enh, convert_enh_to_ecef
from erfa import gc2gd, gd2gc

##This call is so we can use it as a type annotation
woden_struct_classes = Woden_Struct_Classes()
Source_Catalogue = woden_struct_classes.Source_Catalogue
Woden_Settings = woden_struct_classes.Woden_Settings


def worker_get_num_stations(ms_path : str, q : Queue) -> int:
    """
    Worker function to get the number of stations in the measurement set
    given by `ms_path`, and put the in the queue `q`.
    
    This is done in a separate process because it uses `python-casacore`.
    Running `import casacore` creates a `c++` state with a number of library paths
    and casacore global variables. If `libuse_everybeam.so` has been built with
    against a different casacore library, this causes epic intermittent errors
    and segfaults. To void this, run the import in this separate thread process,
    which isolates the state of the `c++` library.
    
    TODO: All calls to `python-casacore` should be done via direct calls to
    casacore in c++ in `libuse_everybeam.so`. This removes all conflicts;
    however, this is a large task and will take time to implement. For now,
    this is a workaround to avoid the segfaults and errors.
    
    Parameters
    ----------
    ms_path : str
        Path to the measurement set.
    q : Queue
        Queue to put the number of stations in.
    """
    
    
    from casacore.tables import table
    
    with table(ms_path + '/ANTENNA') as t: 
        num_stations = len(t)
        
    q.put(num_stations)

def get_num_stations(ms_path : str) -> int:
    """
    Get the number of stations in a measurement set given by `ms_path`.
    
    Runs :func:`~worker_get_num_stations` in a separate process to
    avoid `python-casacore` clashing with WODEN-built `libuse_everybeam.so`.
    See :func:`~worker_get_num_stations` for more details.
    
    Parameters
    ----------
    ms_path : str
        Path to the measurement set.
        
    Returns
    -------
    int
        Number of stations in the measurement set.
    """
    
    q = Queue()
    p = Process(target=worker_get_num_stations, args=(ms_path, q,))
    p.start()
    p.join()

    num_stations = q.get()
        
    return num_stations

def enu_basis(lat : float, lon : float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Calculate the ENU basis vectors in terms of earth-centred earth-fixed coords,
    for a given latitude `lat` and longitude `lon`.
    
    Parameters
    ----------
    lat : float
        Latitude in radians.
    lon : float
        Longitude in radians.
        
    Returns
    -------
    Tuple[np.ndarray, np.ndarray, np.ndarray]
        Tuple of the east, north, and up basis vectors.
    """
    
    east = np.array([-np.sin(lon),  np.cos(lon), 0.0])
    north = np.array([-np.sin(lat) * np.cos(lon),
                      -np.sin(lat) * np.sin(lon),
                       np.cos(lat)])
    up = np.array([np.cos(lat) * np.cos(lon),
                   np.cos(lat) * np.sin(lon),
                   np.sin(lat)])
    return east, north, up


def calc_coordinate_axes(positions : np.ndarray,
                         central_lat : float, central_lon : float):
    """
    Internal to a LOFAR measurement set, inside the `LOFAR_ANTENNA_FIELD` table,
    there is a column called `COORDINATE_AXES`. This has something to do with
    defining what the normal direction of the incoming radation is relative
    to a particular station. If moving a measurement set from one lat/long to
    another, this must be recalculated. Below is a best guess by JLBLine at how
    to do this. See the following doc string copied from EveryBeam which
    describes the coordinate system, from the 
    `EveryBeam/cpp/antenna.h::Antenna::CoordinateSystem` function
    
    .. seealso:: from the EveryBeam documentation:
    
        Station coordinate system.
      
        A right handed, cartesian, local coordinate system with coordinate axes
        `p`, `q`, and `r` is associated with each antenna field.
      
        The r-axis is orthogonal to the antenna field, and points towards the
        local pseudo zenith.
      
        The q-axis is the northern bisector of the X and Y dipoles, i.e.
        it is the reference direction from which the orientation of the dual
        dipole antennae is determined. The q-axis points towards the North at
        the core. At remote sites it is defined as the intersection of the
        antenna field plane and a plane parallel to the meridian plane at the
        core. This ensures the reference directions at all sites are similar.
      
        The p-axis is orthogonal to both other axes, and points towards the East
        at the core.
      
        The axes and origin of the anntena field coordinate system are expressed
        as vectors in the geocentric, cartesian, ITRF coordinate system, in
        meters.
      
        "LOFAR Reference Plane and Reference Direction", M.A. Brentjens,
        LOFAR-ASTRON-MEM-248.
     
    .. todo::
        find out how to do this properly. As described above, everything is
        supposed to be calculated relative to the array centre, or for remote stations,
        "the intersection of the antenna field plane and a plane parallel to the 
        meridian plane at the core". Wtf does that mean? Below, I calculate things
        relative to the local ENU for each station. This means the calculations
        for remote stations are way off, as the curvature of the earth changes
        what local up is the further you get from the core. If people need accurate
        remote stations, this needs to be fixed (or just don't move the telescope).
    
    Parameters
    ----------
    positions : np.ndarray
        Array of station positions in ECEF coordinates. Shape (num_positions, 3).
    central_lat : float
        Latitude of the array centre in radians.
    central_lon : float
        Longitude of the array centre in radians.
        
    Returns
    -------
    coordinate_axes : np.ndarray
        Array of coordinate axes for each station. Shape (num_positions, 3, 3).
        coordinate_axes[:, 0, :] are the `p` vectors, coordinate_axes[:, 1, :]
        are the `q` vectors, and coordinate_axes[:, 2, :] are the `r` vectors, 
        as described in the docstring above.
        
    """
    
    
    arrX, arrY, arrZ = gd2gc(1, central_lon, central_lat, 0)
    arr_centre = np.array([arrX, arrY, arrZ])
    
    num_positions = positions.shape[0]
    coordinate_axes = np.zeros((num_positions, 3, 3), dtype=np.float64)

    for pos in range(num_positions):
        # delta = arr_centre - positions[pos, :]
        delta = positions[pos, :] - arr_centre
        lon, lat, height = gc2gd(1, positions[pos, :])

        east, north, up = enu_basis(lat, lon)
        
        ##IF doing everything relative to the array centre, then
        ##use the array centre as the ENU basis
        # east, north, up = enu_basis(central_lat, central_lon)

        east_comp = np.dot(delta, east)
        north_comp = np.dot(delta, north)
        up_comp = np.dot(delta, up)

        # Reconstruct component vectors
        east_vec = east_comp * east
        north_vec = north_comp * north
        up_vec = up_comp * up
        
        east_vec /= np.linalg.norm(east_vec)
        north_vec /= np.linalg.norm(north_vec)
        up_vec /= np.linalg.norm(up_vec)
        
        
        if lon < central_lon:
            east_vec *= -1
        if lat < central_lat:
            north_vec *= -1
        
        coordinate_axes[pos, 0, :] = east_vec
        coordinate_axes[pos, 1, :] = north_vec
        coordinate_axes[pos, 2, :] = up_vec
        

    return coordinate_axes

def worker_create_filtered_ms(ms_path : str, new_ms_path : str,
                              ra0 : float, dec0 : float,
                              recentre_array : bool = False,
                              current_latitude : float = False,
                              current_longitude : float = False,
                              new_latitude : float = False,
                              new_longitude : float = False) -> None:
    """
    Create a reduced measurement to path `new_ms_path` with only the first time 
    and frequency channel of the measurement set `ms_path`. Set the delay and 
    reference directions to the given RA/Dec `ra0,dec0`.
    
    The `EveryBeam` functions that `WODEN` wraps read the beam centre value
    directly from the measurement set. Rather than downloading a specific
    MS which exactly the pointing we want, we just make the smallest possible
    copy of the MS with the first time and frequency channel, and set
    the delay and reference direction to what we want.
    
    If `recentre_array` is True, also perform the necessary calculations and
    updates to move the array centre from `current_latitude` and `current_longitude`
    to `new_latitude` and `new_longitude`. 
    
    This function is run in a separate process to avoid `python-casacore`
    clashing with WODEN-built `libuse_everybeam.so`. Running `import casacore`
    creates a `c++` state with a number of library paths and casacore global variables.
    If `libuse_everybeam.so` has been built with against a different casacore library,
    this causes epic intermittent errors and segfaults. To void this, run the import
    in this separate thread process, which isolates the state of the `c++` library.
    
    
    Parameters
    ----------
    ms_path : str
        Path to the original measurement set.
    new_ms_path : str
        Path to the new measurement set.
    ra0 : float
        Right ascension of the beam centre in radians.
    dec0 : float
        Declination of the beam centre in radians.
    recentre_array : bool, optional
        Whether to recentre the array to the new latitude and longitude.
        Defaults to False. If True, `current_latitude`, `current_longitude`,
        `new_latitude`, and `new_longitude` must be set.
    current_latitude : float, optional
        Current latitude of the array in radians. Required if `recentre_array`
        is True.
    current_longitude : float, optional
        Current longitude of the array in radians. Required if `recentre_array`
        is True.
    new_latitude : float, optional
        New latitude of the array in radians. Required if `recentre_array`
        is True.
    new_longitude : float, optional
        New longitude of the array in radians. Required if `recentre_array`
        is True.
    """
    
    from casacore.tables import table, taql
    
    ##Find out telescope type
    with table(ms_path + "/OBSERVATION", ack=False) as tb:
        telescope_name = tb.getcol("TELESCOPE_NAME")[0]  # Get first entry
        
    ##First up, read in original MS, and filter it to only have 
    ## the first time and frequency channel
    with table(ms_path, readonly=True) as ms:
        time_col = ms.getcol("TIME")
        ddid_col = ms.getcol("DATA_DESC_ID")

        first_time = time_col[0]  # First timestamp
        first_ddid = ddid_col[0]  # First frequency channel

        # Use TaQL (Table Query Language) to select the subset efficiently
        query = f"SELECT * FROM {ms_path} WHERE TIME = {first_time} AND DATA_DESC_ID = {first_ddid}"
        filtered_ms = taql(query)
        
        filtered_ms.copy(new_ms_path, deep=True)
        
        filtered_ms.close()
        
    with table(new_ms_path+'::FIELD', readonly=False) as field_table:
        field_table.putcol('PHASE_DIR', np.array([[[ra0, dec0]]]))
        field_table.putcol('DELAY_DIR', np.array([[[ra0, dec0]]]))
        field_table.putcol('REFERENCE_DIR', np.array([[[ra0, dec0]]]))
        
        if telescope_name == "LOFAR":
            field_table.putcol('LOFAR_TILE_BEAM_DIR', np.array([[[ra0, dec0]]]))
            
    if recentre_array:
        
        with table(ms_path + '/ANTENNA', readonly=True) as t: 
            num_ants = len(t)
            num_antennas = num_ants
            ant_locations = np.array([t.getcell('POSITION', ant) for ant in range(num_ants)])
            ##convert from ECEF to ENH, as WODEN starts with enh coords
            east, north, height = convert_ecef_to_enh(ant_locations[:,0],
                                        ant_locations[:,1], ant_locations[:,2],
                                        current_longitude, current_latitude)
            
        ecef_X, ecef_Y, ecef_Z = convert_enh_to_ecef(east, north, height,
                                             new_longitude, new_latitude)
        
        with table(new_ms_path+'::ANTENNA', readonly=False) as antenna:
            
            new_locations = np.zeros((num_antennas, 3), dtype=np.float64)
            new_locations[:, 0] = ecef_X
            new_locations[:, 1] = ecef_Y
            new_locations[:, 2] = ecef_Z
            
            antenna.putcol('POSITION', new_locations)
            
            if telescope_name == "LOFAR":
                antenna.putcol('LOFAR_PHASE_REFERENCE', new_locations)
            
        if telescope_name == "LOFAR":
            
            with table(new_ms_path+'::LOFAR_ANTENNA_FIELD', readonly=False) as antenna:
                lof_ant_position = antenna.getcol('POSITION')
                
                east, north, height = convert_ecef_to_enh(lof_ant_position[:,0],
                                        lof_ant_position[:,1], lof_ant_position[:,2],
                                        current_longitude, current_latitude)
            
                ecef_X, ecef_Y, ecef_Z = convert_enh_to_ecef(east, north, height,
                                                    new_longitude, new_latitude)
                
                new_locations = np.zeros((num_antennas, 3), dtype=np.float64)
                new_locations[:, 0] = ecef_X
                new_locations[:, 1] = ecef_Y
                new_locations[:, 2] = ecef_Z
                
                antenna.putcol('POSITION', new_locations)
                
            new_coord_axes = calc_coordinate_axes(new_locations,
                                                        new_latitude,
                                                        new_longitude)
            
            with table(new_ms_path+'::LOFAR_ANTENNA_FIELD', readonly=False) as antenna:
                
                coord_axes = antenna.getcol('COORDINATE_AXES')
                antenna.putcol('COORDINATE_AXES', new_coord_axes)
                
                arrX, arrY, arrZ = gd2gc(1, new_longitude, new_latitude, 0)
                arr_centre = np.array([arrX, arrY, arrZ])
                
                for station in range(num_antennas):
                    offsets = antenna.getcell('ELEMENT_OFFSET', station)
                    offsets += ant_locations[station, :]
                    
                    cur_lon, cur_lat, height = gc2gd(1, ant_locations[station, :])
                    
                    east, north, height = convert_ecef_to_enh(offsets[:,0],
                                        offsets[:,1], offsets[:,2],
                                        cur_lon, cur_lat)
                    
                    new_lon, new_lat, height = gc2gd(1, new_locations[station, :])
            
                    ecef_X, ecef_Y, ecef_Z = convert_enh_to_ecef(east, north, height,
                                                        new_lon, new_lat)
                    
                    new_offsets = np.zeros_like(offsets)
                    new_offsets[:, 0] = ecef_X - new_locations[station, 0]
                    new_offsets[:, 1] = ecef_Y - new_locations[station, 1]
                    new_offsets[:, 2] = ecef_Z - new_locations[station, 2]
                        
                    antenna.putcell('ELEMENT_OFFSET', station, new_offsets)
                    
                for station in range(num_antennas):
                    offsets = antenna.getcell('TILE_ELEMENT_OFFSET', station)
                    
                    offsets += ant_locations[station, :]
                    
                    cur_lon, cur_lat, height = gc2gd(1, ant_locations[station, :])
                    
                    east, north, height = convert_ecef_to_enh(offsets[:,0],
                                        offsets[:,1], offsets[:,2],
                                        cur_lon, cur_lat)
                    
                    new_lon, new_lat, height = gc2gd(1, new_locations[station, :])
            
                    ecef_X, ecef_Y, ecef_Z = convert_enh_to_ecef(east, north, height,
                                                        new_lon, new_lat)
                    
                    new_offsets = np.zeros_like(offsets)
                    new_offsets[:, 0] = ecef_X - new_locations[station, 0]
                    new_offsets[:, 1] = ecef_Y - new_locations[station, 1]
                    new_offsets[:, 2] = ecef_Z - new_locations[station, 2]
                    
                    antenna.putcell('TILE_ELEMENT_OFFSET', station, new_offsets)


def create_filtered_ms(ms_path : str, new_ms_path : str,
                       ra0 : float, dec0 : float,
                       recentre_array : bool = False,
                       current_latitude : float = False,
                       current_longitude : float = False,
                       new_latitude : float = False,
                       new_longitude : float = False):
    """
    Create a reduced measurement to path `new_ms_path` with only the first time
    and frequency channel of the measurement set `ms_path`. Set the delay and
    reference directions to the given RA/Dec `ra0,dec0`. If `recentre_array`
    is True, also perform the necessary calculations and updates to move the
    array centre from `current_latitude` and `current_longitude` to 
    `new_latitude` and `new_longitude`.
    
    Runs :func:`~worker_create_filtered_ms` in a separate process to 
    avoid `python-casacore` clashing with WODEN-built `libuse_everybeam.so`.
    See :func:`~worker_create_filtered_ms` for more details.
    
    Parameters
    ----------
    ms_path : str
        Path to the original measurement set.
    new_ms_path : str
        Path to the new measurement set.
    ra0 : float
        Right ascension of the beam centre in radians.
    dec0 : float
        Declination of the beam centre in radians.
    recentre_array : bool, optional
        Whether to recentre the array to the new latitude and longitude.
        Defaults to False. If True, `current_latitude`, `current_longitude`,
        `new_latitude`, and `new_longitude` must be set.
    current_latitude : float, optional
        Current latitude of the array in radians. Required if `recentre_array`
        is True.
    current_longitude : float, optional
        Current longitude of the array in radians. Required if `recentre_array`
        is True.
    new_latitude : float, optional
        New latitude of the array in radians. Required if `recentre_array`
        is True.
    new_longitude : float, optional
        New longitude of the array in radians. Required if `recentre_array`
        is True.
    """
    

    p = Process(target=worker_create_filtered_ms,
                args=(ms_path, new_ms_path, ra0, dec0,
                      recentre_array, current_latitude,
                      current_longitude, new_latitude,
                      new_longitude))
    p.start()
    p.join()  # module is fully cleaned when the process exits
                
            
def check_ms_telescope_type_matches_element_response(ms_path : str,
                                                     element_response_model : str = 'default',
                                                     logger : Logger = False) -> Tuple[str, str]:
    """
    Check the requested `element_response_model` is compatible with the
    telescope type in the measurement set. If they don't match, set the element
    response model to the default for the telescope type. This function uses
    the `check_ms_telescope_type` function in `libuse_everybeam.so`.
    
    - LOFAR: 'hamaker', 'hamakerlba', 'lobes' or 'default'
    - OSKAR: 'skala40_wave' or 'default'
    
    Please note that 'hamakerlba', 'lobes' are not tested or supported in
    WODEN, and are not recommended for use. They are included for completeness.
    
    Furthermore, please note that WODEN doesn't use a base measurement set
    for the MWA primary beam, and so will throw an error here to make sure
    you're not doing something the c++ code can't handle.
    
    Parameters
    ----------
    ms_path : str
        Path to the measurement set.
    element_response_model : str, optional
        Requested Element response model to use. Defaults to 'default'.
    logger : Logger, optional
        Logger to use. Defaults to False; if False, create a new simple logger instance.
        
    Returns
    -------
    Tuple[str, str]
        Tuple of the telescope type and the element response model to use.
    """

    
    
    if not logger:
        logger = simple_logger()
    
    woden_path = importlib_resources.files(wodenpy).joinpath(f"libwoden_double.so")
    woden_lib = ctypes.cdll.LoadLibrary(woden_path)
    
    check_ms_telescope_type = woden_lib.check_ms_telescope_type
    check_ms_telescope_type.argtypes = [c_char_p]
    check_ms_telescope_type.restype = c_char_p
    
    # print(f"MS path: {type(ms_path)}, {ms_path}")
    
    ms_path_ctypes = ctypes.c_char_p(ms_path.encode('utf-8'))
    
    telescope_type = check_ms_telescope_type(ms_path_ctypes).decode('utf-8')
    
    use_element_response_model = False
    
    if telescope_type == 'MWA':
        exit_message = "WODEN does not support using a measurement set for the MWA primary beam model. "
        exit_message += "It instead calls the EveryBeam Beam2016Implementation code "
        exit_message += "directly to input az,za, and uses the array layout from either "
        exit_message += "an MWA metafits file, or a user-defined array layout. "
        
        logger.error(exit_message)
        exit(exit_message)
                
    elif telescope_type == 'LOFAR':
        if element_response_model == 'default':
            use_element_response_model = "hamaker"
        else:
            if element_response_model not in ['hamaker', 'hamakerlba', 'lobes']:
                logger.warning(f"Measurement set telescope type is LOFAR, but "
                               f"element_response_model was set to {element_response_model}, "
                               "which is not one of ['hamaker', 'hamakerlba','lobes']. "
                               "Defaulting to 'hamaker'")
                use_element_response_model = "hamaker"
            else:
                use_element_response_model = element_response_model
                
    elif telescope_type == 'OSKAR':
        if element_response_model == 'default':
            use_element_response_model = "skala40_wave"
        else:
            if element_response_model != 'skala40_wave':
                logger.warning(f"Measurement set telescope type is OSKAR, but element_response_model "
                               f"was set to {element_response_model}. Changing to 'skala40_wave'")
                use_element_response_model = "skala40_wave"
            else:
                use_element_response_model = element_response_model
    else:
        use_element_response_model = element_response_model
        exit_message = f"Measurement set telescope type is {telescope_type}. "
        exit_message += "WODEN currently supports LOFAR, OSKAR EveryBeam from a measurement set. "
        exit_message += "Cannot proceed as unknown behaviour will happen in C++ code. "
        exit_message += "Exiting now."
        logger.error(exit_message)
        exit(exit_message)
        
    return telescope_type, use_element_response_model


def convert_common_args_to_everybeam_args(ms_path : str, coeff_path : str,
                                          element_response_model : str,
                                          station_idxs : np.ndarray,
                                          freqs : np.ndarray,
                                          mjd_sec_times : np.ndarray) -> Tuple:
    """
    Convert input arugments common to all EveryBeam primary beam models
    into `ctypes` equivalents to pass to the C++ library.
    
    Parameters
    ----------
    ms_path : str
        Path to the measurement set.
    coeff_path : str
        Path to the coefficients file (only needed for MWA, pass the path
        to the hdf5 FEE file).
    element_response_model : str
        Element response model to use. Can be 'MWA', 'hamaker', 'skala40_wave',
        or 'default'.
    station_idxs : np.ndarray
        Array of station indices to use.
    freqs : np.ndarray
        Array of frequencies to use.
    mjd_sec_times : np.ndarray
        Array of times in MJD seconds.
        
    Returns
    -------
    Tuple[ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double)]
        Tuple of the converted arguments in the format required by the
        EveryBeam library. Order of outputs is identical to order of inputs arguments.
    """
    
    num_stations = len(station_idxs)
    num_times = len(mjd_sec_times)
    num_freqs = len(freqs)
                                          
    mjd_sec_times_ctypes = (ctypes.c_double * num_times)()
    for i in range(num_times):
        mjd_sec_times_ctypes[i] = mjd_sec_times[i]
        
    freqs_ctypes = (ctypes.c_double * num_freqs)()
    for i in range(num_freqs):
        freqs_ctypes[i] = freqs[i]
    
    station_idxs_ctypes = (ctypes.c_int * num_stations)()
    for i in range(num_stations):
        station_idxs_ctypes[i] = station_idxs[i]
    
    
    ms_path_ctypes = ctypes.c_char_p(ms_path.encode('utf-8'))
    element_response_model_ctypes = ctypes.c_char_p(element_response_model.encode('utf-8'))
    coeff_path_ctypes = ctypes.c_char_p(coeff_path.encode('utf-8'))
    
    return ms_path_ctypes, coeff_path_ctypes, element_response_model_ctypes, station_idxs_ctypes, freqs_ctypes, mjd_sec_times_ctypes

def run_everybeam(ras: np.ndarray, decs: np.ndarray,
                  freqs: np.ndarray,
                  
                  ms_path : str = False,
                  beam_ra0: float = np.nan, beam_dec0: float = np.nan,
                  times: np.ndarray[Time] = False,
                  station_ids: np.ndarray = False,
                  element_response_model='default',
                  apply_beam_norms: bool = True,
                  
                  mwa_coeff_path : str = False,
                  mwa_dipole_delays: np.ndarray = False,
                  mwa_dipole_amps: np.ndarray = np.ones(16),
                  j2000_latitudes: np.ndarray = False,
                  j2000_lsts: np.ndarray = False,
                  
                  iau_order: bool = False,
                  element_only: bool = False,
                  parallactic_rotate: bool = False,
                  logger : Logger = False) -> np.ndarray:
    """
    Calculate the Jones matrices for a given set of coordinates, times,
    frequencies, and station ids using the EveryBeam library. Currently,
    LOFAR and OSKAR beams require a measurement set as an input. The MWA beam
    does not need a measurement set, but does need the path to the FEE
    coefficients hdf5 file.
    
    For all beam types, `ras`, `decs`, `freqs` are required.
    
    For LOFAR/OSKAR beams, `ms_path`, `beam_ra0`, `beam_dec0`, `times`,
    and `station_ids` are required.

    For MWA beams, `mwa_coeff_path`, `mwa_dipole_delays`, `j2000_latitudes`,
    and `j2000_lsts` are required. `j2000_latitudes` should be the array
    latitude as precessed back to J2000, with `j2000_lsts` being the
    matching LST in J2000. 
    
    The returned Jones matrices will have shape
    (num_stations, num_times, num_freqs, num_radec, 2, 2). The MWA beam will
    always default to `num_stations=1` as dipole flagging has not been
    implemented in the WODEN wrapper. All other beams will have
    `num_stations=len(station_ids)`.
    
    Parameters
    ------------
    ras : np.ndarray
        Right ascensions of the coordinates in radians.
    decs : np.ndarray
         Declinations of the coordinates in radians.
    freqs : np.ndarray
        Array of frequencies (Hz)
        
    ms_path : str
        Path to the measurement set to load the EveryBeam telescope from.
    beam_ra0 : float
        Right ascension (radians) of the beam pointing centre; the beam will
        be centered on this position and therefore move with time. Does not
        apply to MWA beam as that is locked in az,za.
    beam_dec0 : float
        Declination (radians) of the beam pointing centre; the beam will
        be centered on this position and therefore move with time. Does not
        apply to MWA beam as that is locked in az,za.
    times : np.ndarray
        Array of `astropy.Time` observation times.
    station_ids : np.ndarray
        Array of station IDs. Will calculation n_stations=len(station_ids)
        worth of station beams (does not apply to MWA, which will only
        calculate a single beam).
    element_response_model : str, optional
        The Everybeam element response model to use. Defaults to 'hamaker'.
        Avaible options are 'hamaker' (LOFAR), 'skala40_wave' (OSKAR)
    apply_beam_norms : bool, optional
        Whether to apply beam normalisation. Defaults to True. Achieved by
        calculating the beam response at beam centre, and multiplying all
        Jones by the inverse of this central beam response (does not apply
        to MWA beam).
        
    mwa_coeff_path : str
        Path to the coefficients file (needed for MWA, pass the path
        to the hdf5 FEE file.)
    mwa_dipole_delays : np.ndarray
        Array of dipole delays for the MWA stations. Must be of length 16.
    mwa_dipole_amps : np.ndarray, optional
        Array of dipole amplitudes for the MWA stations. Must be of length 16.
    j2000_latitudes : np.ndarray
        Latitudes in J2000 coordinates (needed for MWA beam az,za calculations).
    j2000_lsts : np.ndarray
        Local sidereal times in J2000 coordinates (needed for MWA beam az,za calculations).
    
    iau_order : bool, optional
        If True, use IAU polarisation ordering, so set jones[0,0] to the NS dipole
        and jones[1,1] to EW. If False, jones[0,0] is EW. Defaults to False
    element_only : bool, optional
        Whether to use only the element response. Defaults to False. Use this to
        look at the dipole response only, not the beam formed response.
    parallactic_rotate : bool, optional
        Whether to apply parallactic angle rotation. Defaults to False.
    logger : Logger, optional
        Logger to use. Defaults to False; if False, create a new simple logger instance.
    
        
    Returns
    --------
    np.ndarray
        The calculated Jones matrices with shape (num_stations, num_times, num_freqs, num_coords, 2, 2).
    """

    if type(ms_path) == bool and type(mwa_coeff_path) == bool:
            logger.error("Either `ms_path` or `mwa_coeff_path` must be provided. Exiting now.")
            exit()
    
    if not logger:
        logger = simple_logger()
        

    if mwa_coeff_path:
        
        ##Check coeff path exists
        if not os.path.exists(mwa_coeff_path):
            logger.error(f"MWA coefficient path mwa_coeff_path={mwa_coeff_path} does not exist. Exiting now.")
            exit()
        
        if type(j2000_latitudes) == bool or type(j2000_lsts) == bool or type(mwa_dipole_delays) == bool:
            logger.error("j2000_latitudes, j2000_lsts, and mwa_dipole_delays must be provided for MWA beam calculations. Exiting now.")
            exit()
        
        jones = run_mwa_beam(mwa_dipole_delays,
                             mwa_coeff_path,
                             ras, decs, 
                             j2000_lsts, j2000_latitudes, freqs,
                             amps=mwa_dipole_amps,
                             apply_beam_norms=apply_beam_norms,
                             parallactic_rotate=parallactic_rotate,
                             iau_order=iau_order,
                             element_only=element_only)
    
    else:
        if type(ms_path) == bool:
            logger.error("ms_path must be provided for LOFAR and OSKAR beam calculations. Exiting now.")
            exit()
        
        telescope_type, checked_element_response_model = check_ms_telescope_type_matches_element_response(ms_path,
                                                                           element_response_model,
                                                                           logger)
        
        if type(station_ids) == bool or type(times) == bool or type(beam_ra0) == bool or type(beam_dec0) == bool:
            logger.error("station_ids, times, beam_ra0, and beam_dec0 must be provided for LOFAR and OSKAR beam calculations. Exiting now.")
            exit()
        
        mjd_sec_times = np.array([time.mjd * 86400.0 for time in times])
        

        if telescope_type == 'LOFAR':
            jones = run_lofar_beam(ms_path, checked_element_response_model,
                               station_ids,
                               beam_ra0, beam_dec0,
                               ras, decs, mjd_sec_times, freqs,
                               apply_beam_norms=apply_beam_norms,
                               parallactic_rotate=parallactic_rotate,
                               iau_order=iau_order,
                               element_only=element_only)
        
        elif telescope_type == 'OSKAR':
            jones = run_oskar_beam(ms_path, checked_element_response_model,
                                station_ids,
                                beam_ra0, beam_dec0,
                                ras, decs, mjd_sec_times, freqs,
                                apply_beam_norms=apply_beam_norms,
                                parallactic_rotate=parallactic_rotate,
                                iau_order=iau_order,
                                element_only=element_only)
            
        else:
            logger.error("Unknown telescope type. Exiting now.")
            exit()
    
    return jones


def run_everybeam_thread(num_threads : int, thread_id : int,
                         ras: np.ndarray, decs: np.ndarray,
                         freqs: np.ndarray,
                         
                         ms_path : str = False,
                         beam_ra0: float = np.nan, beam_dec0: float = np.nan,
                         times: np.ndarray[Time] = False,
                         station_ids: np.ndarray = False,
                         element_response_model='default',
                         apply_beam_norms: bool = True,
                         
                         mwa_coeff_path : str = False,
                         mwa_dipole_delays: np.ndarray = False,
                         mwa_dipole_amps: np.ndarray = np.ones(16),
                         j2000_latitudes: np.ndarray = False,
                         j2000_lsts: np.ndarray = False,
                         
                         iau_order: bool = False,
                         element_only: bool = False,
                         parallactic_rotate: bool = False,
                         logger : Logger = False) -> Tuple[np.ndarray, int]:
    """
    Thread function called by `run_everybeam_over_threads` to calculate the
    EveryBeam response in parrallel. Calls `run_everybeam` with a subset of
    the coordinates; see `run_everybeam` for more details of the parameters.
    
    Creates a new EveryBeam telescope object from `ms_path` or a new EveryBeam
    MWA Tile beam from `mwa_coeff_path` for each thread.
    This has to be done because `concurrent.futures.ProcessPoolExecutor` has
    to pickle the function and all it's arguments, and EveryBeam objects can't
    be pickled. This is somewhat wasteful but I can't work out a better way
    to make things parallel.
    
    Parameters
    ------------
    num_threads : int
        Number of threads being in call by `run_everybeam_over_threads`.
    thread_id : int
        ID of the current thread. Useds to work out what chunk of `ras` and `decs`
        to process.
    ras : np.ndarray
        Right ascensions of the coordinates in radians.
    decs : np.ndarray
         Declinations of the coordinates in radians.
    freqs : np.ndarray
        Array of frequencies (Hz)
        
    ms_path : str
        Path to the measurement set to load the EveryBeam telescope from.
    beam_ra0 : float
        Right ascension (radians) of the beam pointing centre; the beam will
        be centered on this position and therefore move with time. Does not
        apply to MWA beam as that is locked in az,za.
    beam_dec0 : float
        Declination (radians) of the beam pointing centre; the beam will
        be centered on this position and therefore move with time. Does not
        apply to MWA beam as that is locked in az,za.
    times : np.ndarray
        Array of `astropy.Time` observation times.
    station_ids : np.ndarray
        Array of station IDs. Will calculation n_stations=len(station_ids)
        worth of station beams (does not apply to MWA, which will only
        calculate a single beam).
    element_response_model : str, optional
        The Everybeam element response model to use. Defaults to 'hamaker'.
        Avaible options are 'hamaker' (LOFAR), 'skala40_wave' (OSKAR)
    apply_beam_norms : bool, optional
        Whether to apply beam normalisation. Defaults to True. Achieved by
        calculating the beam response at beam centre, and multiplying all
        Jones by the inverse of this central beam response (does not apply
        to MWA beam).
        
    mwa_coeff_path : str
        Path to the coefficients file (needed for MWA, pass the path
        to the hdf5 FEE file.)
    mwa_dipole_delays : np.ndarray
        Array of dipole delays for the MWA stations. Must be of length 16.
    mwa_dipole_amps : np.ndarray, optional
        Array of dipole amplitudes for the MWA stations. Must be of length 16.
    j2000_latitudes : np.ndarray
        Latitudes in J2000 coordinates (needed for MWA beam az,za calculations).
    j2000_lsts : np.ndarray
        Local sidereal times in J2000 coordinates (needed for MWA beam az,za calculations).
    
    iau_order : bool, optional
        If True, use IAU polarisation ordering, so set jones[0,0] to the NS dipole
        and jones[1,1] to EW. If False, jones[0,0] is EW. Defaults to False
    element_only : bool, optional
        Whether to use only the element response. Defaults to False. Use this to
        look at the dipole response only, not the beam formed response.
    parallactic_rotate : bool, optional
        Whether to apply parallactic angle rotation. Defaults to False.
    logger : Logger, optional
        Logger to use. Defaults to False; if False, create a new simple logger instance.
        
    Returns
    --------
    Tuple[np.ndarray, int]
        The calculated Jones matrices with shape
        (num_stations, num_times, num_freqs, num_coords_in_thread, 2, 2), as
        well as the thread ID. Use the thread ID to insert this thread output
        into the correct place in the final Jones matrix.
        
    """
    
    num_coords = len(ras)
    coords_per_thread = int(np.ceil(num_coords / num_threads))
    
    low_coord = thread_id * coords_per_thread
    high_coord = (thread_id + 1) * coords_per_thread
    
    # print(f"Thread {thread_id} processing coords {low_coord} to {high_coord}")
    
    jones = run_everybeam(ras[low_coord:high_coord],
                          decs[low_coord:high_coord],
                          freqs,
                          ms_path=ms_path,
                          times=times,
                          beam_ra0=beam_ra0, beam_dec0=beam_dec0,
                          mwa_coeff_path=mwa_coeff_path,
                          mwa_dipole_delays=mwa_dipole_delays,
                          mwa_dipole_amps=mwa_dipole_amps,
                          j2000_latitudes=j2000_latitudes,
                          j2000_lsts=j2000_lsts,
                          station_ids=station_ids,
                          element_response_model=element_response_model,
                          apply_beam_norms=apply_beam_norms,
                          iau_order=iau_order,
                          element_only=element_only,
                          parallactic_rotate=parallactic_rotate,
                          logger=logger)
    
    # print(f"Thread {thread_id} finished")
    
    return jones, thread_id

def run_everybeam_over_threads(num_threads : int,
                               ras: np.ndarray, decs: np.ndarray,
                               freqs: np.ndarray,
                               
                               ms_path : str = False,
                               beam_ra0: float = np.nan, beam_dec0: float = np.nan,
                               times: np.ndarray[Time] = False,
                               station_ids: np.ndarray = False,
                               element_response_model='default',
                               apply_beam_norms: bool = True,
                               
                               mwa_coeff_path : str = False,
                               mwa_dipole_delays: np.ndarray = False,
                               mwa_dipole_amps: np.ndarray = np.ones(16),
                               j2000_latitudes: np.ndarray = False,
                               j2000_lsts: np.ndarray = False,
                               
                               iau_order: bool = False,
                               element_only: bool = False,
                               parallactic_rotate: bool = False,
                               logger : Logger = False):
    """
    Runs `run_everybeam` in parallel over `num_threads` threads, using
    `concurrent.futures.ProcessPoolExecutor`. See `run_everybeam` for more
    details of what each parameter does.
    
    Creates a new EveryBeam telescope object from `ms_path` or a new EveryBeam
    MWA Tile beam from `mwa_coeff_path` for each thread.
    This has to be done because `concurrent.futures.ProcessPoolExecutor` has
    to pickle the function and all it's arguments, and EveryBeam objects can't
    be pickled. This is somewhat wasteful but I can't work out a better way
    to make things parallel (without using MPI on the other end of things).
    
    Parameters
    ------------
    num_threads : int
        Number of threads being in call by `run_everybeam_over_threads`.
        ras : np.ndarray
        Right ascensions of the coordinates in radians.
    ras : np.ndarray
        Right ascensions of the coordinates in radians.
    decs : np.ndarray
         Declinations of the coordinates in radians.
    freqs : np.ndarray
        Array of frequencies (Hz)
        
    ms_path : str
        Path to the measurement set to load the EveryBeam telescope from.
    beam_ra0 : float
        Right ascension (radians) of the beam pointing centre; the beam will
        be centered on this position and therefore move with time. Does not
        apply to MWA beam as that is locked in az,za.
    beam_dec0 : float
        Declination (radians) of the beam pointing centre; the beam will
        be centered on this position and therefore move with time. Does not
        apply to MWA beam as that is locked in az,za.
    times : np.ndarray
        Array of `astropy.Time` observation times.
    station_ids : np.ndarray
        Array of station IDs. Will calculation n_stations=len(station_ids)
        worth of station beams (does not apply to MWA, which will only
        calculate a single beam).
    element_response_model : str, optional
        The Everybeam element response model to use. Defaults to 'hamaker'.
        Avaible options are 'hamaker' (LOFAR), 'skala40_wave' (OSKAR)
    apply_beam_norms : bool, optional
        Whether to apply beam normalisation. Defaults to True. Achieved by
        calculating the beam response at beam centre, and multiplying all
        Jones by the inverse of this central beam response (does not apply
        to MWA beam).
        
    mwa_coeff_path : str
        Path to the coefficients file (needed for MWA, pass the path
        to the hdf5 FEE file.)
    mwa_dipole_delays : np.ndarray
        Array of dipole delays for the MWA stations. Must be of length 16.
    mwa_dipole_amps : np.ndarray, optional
        Array of dipole amplitudes for the MWA stations. Must be of length 16.
    j2000_latitudes : np.ndarray
        Latitudes in J2000 coordinates (needed for MWA beam az,za calculations).
    j2000_lsts : np.ndarray
        Local sidereal times in J2000 coordinates (needed for MWA beam az,za calculations).
    
    iau_order : bool, optional
        If True, use IAU polarisation ordering, so set jones[0,0] to the NS dipole
        and jones[1,1] to EW. If False, jones[0,0] is EW. Defaults to False
    element_only : bool, optional
        Whether to use only the element response. Defaults to False. Use this to
        look at the dipole response only, not the beam formed response.
    parallactic_rotate : bool, optional
        Whether to apply parallactic angle rotation. Defaults to False.
    logger : Logger, optional
        Logger to use. Defaults to False; if False, create a new simple logger instance.
        
    Returns
    --------
    np.ndarray
        The calculated Jones matrices with shape
        (num_stations, num_times, num_freqs, num_coord, 2, 2)
    """
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
        
        future_data = [executor.submit(run_everybeam_thread,
                                       num_threads, thread_id,
                                       ras, decs, freqs,
                                       ms_path=ms_path,
                                       times=times,
                                       beam_ra0=beam_ra0, beam_dec0=beam_dec0,
                                       mwa_coeff_path=mwa_coeff_path,
                                       mwa_dipole_delays=mwa_dipole_delays,
                                       mwa_dipole_amps=mwa_dipole_amps,
                                       j2000_latitudes=j2000_latitudes,
                                       j2000_lsts=j2000_lsts,
                                       station_ids=station_ids,
                                       element_response_model=element_response_model,
                                       apply_beam_norms=apply_beam_norms,
                                       iau_order=iau_order,
                                       element_only=element_only,
                                       parallactic_rotate=parallactic_rotate,
                                       logger=logger)
                                for thread_id in range(num_threads)]
            
        all_jones_chunks = []
        all_thread_ids = []
        for future in concurrent.futures.as_completed(future_data):
            jones_chunk, thread_id = future.result()
            all_jones_chunks.append(jones_chunk)
            all_thread_ids.append(thread_id)
    
    if type(mwa_coeff_path) != bool:
        num_stations = 1
        num_times = len(j2000_lsts)
    else:
        num_stations = len(station_ids)
        num_times = len(times)
        
    num_freqs = len(freqs)
    num_coords = len(ras)
    
    all_jones = np.zeros((num_stations, num_times, num_freqs, num_coords, 2, 2), dtype=np.complex128)*np.nan
    
    coords_per_thread = int(np.ceil(num_coords / num_threads))
    
    for jones_chunk, thread_id in zip(all_jones_chunks, all_thread_ids):
        
        low_coord = thread_id * coords_per_thread
        high_coord = (thread_id + 1) * coords_per_thread
        
        all_jones[:, :, :, low_coord:high_coord, :, :] = jones_chunk
    
    return all_jones

def run_lofar_beam(ms_path : str, element_response_model : bool,
                   station_idxs : np.ndarray,
                   beam_ra0 : float, beam_dec0 : float,
                   ras : np.ndarray, decs : np.ndarray,
                   mjd_sec_times : np.ndarray,
                   freqs : np.ndarray,
                   apply_beam_norms : bool = True,
                   parallactic_rotate : bool = False,
                   iau_order : bool = True,
                   element_only : bool = False):
    """
    Run the LOFAR beam model using the EveryBeam library.
    This function is a wrapper around the C++ library `libuse_everybeam.so`.
    It takes the input parameters and converts them to the appropriate
    ctypes types, then calls the C++ function to calculate the beam response.
    The output is a numpy array of the beam response, with shape
    (num_stations, num_times, num_freqs, num_coords, 2, 2).
    
    Parameters
    ----------
    ms_path : str
        Path to the measurement set.
    element_response_model : str
        Element response model to use. Can be 'hamaker', 'hamakerlba', 'lobes' or 'default'.
    station_idxs : np.ndarray
        Array of station indices to use.
    beam_ra0 : float
        Right ascension of the beam center in radians.
    beam_dec0 : float
        Declination of the beam center in radians.
    ras : np.ndarray
        Right ascensions of the coordinates in radians.
    decs : np.ndarray
        Declinations of the coordinates in radians.
    mjd_sec_times : np.ndarray
        Array of times in MJD seconds.
    freqs : np.ndarray
        Array of frequencies to use.
    apply_beam_norms : bool
        Whether to apply beam normalisation. Defaults to True. Achieved by
        calculating the beam response at beam centre, and multiplying all
        Jones by the inverse of this central beam response.
    parallactic_rotate : bool
        Whether to apply parallactic angle rotation. Defaults to False.
    iau_order : bool
        If True, use IAU polarisation ordering, so set jones[0,0] to the NS dipole and jones[1,1] to EW. If False, jones[0,0] is EW.
    element_only : bool
        Whether to use only the element response. Defaults to False. Use this to
        look at the dipole response only, not the beam formed response.
    Returns
    -------
    np.ndarray
        The calculated Jones matrices with shape
        (num_stations, num_times, num_freqs, num_coords, 2, 2).
    """
    
    woden_path = importlib_resources.files(wodenpy).joinpath(f"libuse_everybeam.so")
    woden_lib = ctypes.cdll.LoadLibrary(woden_path)
    
    load_and_run_lofar_beam = woden_lib.load_and_run_lofar_beam
    
    num_stations = len(station_idxs)
    num_dirs = len(ras)
    num_freqs = len(freqs)
    num_times = len(mjd_sec_times)
    
    ras_ctypes = (ctypes.c_double * num_dirs)()
    decs_ctypes = (ctypes.c_double * num_dirs)()
    for i in range(num_dirs):
        ras_ctypes[i] = ras[i]
        decs_ctypes[i] = decs[i]
        
    mjd_sec_times_ctypes = (ctypes.c_double * num_times)()
    for i in range(num_times):
        mjd_sec_times_ctypes[i] = mjd_sec_times[i]
        
    freqs_ctypes = (ctypes.c_double * num_freqs)()
    for i in range(num_freqs):
        freqs_ctypes[i] = freqs[i]
    
    station_idxs_ctypes = (ctypes.c_int * num_stations)()
    for i in range(num_stations):
        station_idxs_ctypes[i] = station_idxs[i]
    
    
    ms_path_ctypes = ctypes.c_char_p(ms_path.encode('utf-8'))
    element_response_model_ctypes = ctypes.c_char_p(element_response_model.encode('utf-8'))
    
    jones = ((num_stations*num_times*num_freqs*num_dirs*4)*c_double_complex)()
    
    load_and_run_lofar_beam.argtypes = [c_char_p, c_char_p,
                                        c_int, POINTER(c_int),
                                        c_int, c_double, c_double,
                                        POINTER(c_double), POINTER(c_double),
                                        c_int, POINTER(c_double),
                                        c_int, POINTER(c_double),
                                        c_bool, c_bool, c_bool, c_bool,
                                        POINTER(c_double_complex)]
    
    load_and_run_lofar_beam(ms_path_ctypes,
                            element_response_model_ctypes,
                            num_stations, station_idxs_ctypes,
                            num_dirs,
                            beam_ra0, beam_dec0,
                            ras_ctypes, decs_ctypes,
                            num_times, mjd_sec_times_ctypes,
                            num_freqs, freqs_ctypes,
                            apply_beam_norms, parallactic_rotate,
                            element_only, iau_order,
                            jones)
    
    # print(jones)
    
    jones_py = np.ctypeslib.as_array(jones, shape=(num_stations*num_times*num_freqs*num_dirs*4))
    jones_py = jones_py['real'] + 1j*jones_py['imag']
    
    jones_py = jones_py.reshape(num_stations, num_times, num_freqs, num_dirs, 2, 2)
    
    
    return jones_py

def run_oskar_beam(ms_path : str, element_response_model : bool,
                   station_idxs : np.ndarray,
                   beam_ra0 : float, beam_dec0 : float,
                   ras : np.ndarray, decs : np.ndarray,
                   mjd_sec_times : np.ndarray,
                   freqs : np.ndarray,
                   apply_beam_norms : bool,
                   parallactic_rotate : bool,
                   iau_order : bool = True,
                   element_only : bool = False):
    """
    Run the OSKAR beam model using the EveryBeam library.
    
    This function is a wrapper around the C++ library `libuse_everybeam.so`.
    It takes the input parameters and converts them to the appropriate
    ctypes types, then calls the C++ function to calculate the beam response.
    The output is a numpy array of the beam response, with shape
    (num_stations, num_times, num_freqs, num_coords, 2, 2).
    
    NOTE: This function is currently identical to the LOFAR beam function.
    Both LOFAR and OSKAR call the PhasedArrayBeam class in the C++ library.
    These functions have only been lightly tested for functionality. So I've
    kept them separate in case it's discovered they need different default
    values in the future.
    
    Parameters
    ----------
    ms_path : str
        Path to the measurement set.
    element_response_model : str
        Element response model to use. Can be 'hamaker', 'hamakerlba', 'lobes' or 'default'.
    station_idxs : np.ndarray
        Array of station indices to use.
    beam_ra0 : float
        Right ascension of the beam center in radians.
    beam_dec0 : float
        Declination of the beam center in radians.
    ras : np.ndarray
        Right ascensions of the coordinates in radians.
    decs : np.ndarray
        Declinations of the coordinates in radians.
    mjd_sec_times : np.ndarray
        Array of times in MJD seconds.
    freqs : np.ndarray
        Array of frequencies to use.
    apply_beam_norms : bool
        Whether to apply beam normalisation. Defaults to True. Achieved by
        calculating the beam response at beam centre, and multiplying all
        Jones by the inverse of this central beam response.
    parallactic_rotate : bool
        Whether to apply parallactic angle rotation. Defaults to False.
    iau_order : bool
        If True, use IAU polarisation ordering, so set jones[0,0] to the NS dipole and jones[1,1] to EW. If False, jones[0,0] is EW.
    element_only : bool
        Whether to use only the element response. Defaults to False. Use this to
        look at the dipole response only, not the beam formed response.
    Returns
    -------
    np.ndarray
        The calculated Jones matrices with shape
        (num_stations, num_times, num_freqs, num_coords, 2, 2).
    """
    
    woden_path = importlib_resources.files(wodenpy).joinpath(f"libuse_everybeam.so")
    woden_lib = ctypes.cdll.LoadLibrary(woden_path)
    
    load_and_run_oskar_beam = woden_lib.load_and_run_oskar_beam
    
    num_stations = len(station_idxs)
    num_dirs = len(ras)
    num_freqs = len(freqs)
    num_times = len(mjd_sec_times)
    
    ras_ctypes = (ctypes.c_double * num_dirs)()
    decs_ctypes = (ctypes.c_double * num_dirs)()
    for i in range(num_dirs):
        ras_ctypes[i] = ras[i]
        decs_ctypes[i] = decs[i]
        
    mjd_sec_times_ctypes = (ctypes.c_double * num_times)()
    for i in range(num_times):
        mjd_sec_times_ctypes[i] = mjd_sec_times[i]
        
    freqs_ctypes = (ctypes.c_double * num_freqs)()
    for i in range(num_freqs):
        freqs_ctypes[i] = freqs[i]
    
    station_idxs_ctypes = (ctypes.c_int * num_stations)()
    for i in range(num_stations):
        station_idxs_ctypes[i] = station_idxs[i]
    
    
    ms_path_ctypes = ctypes.c_char_p(ms_path.encode('utf-8'))
    element_response_model_ctypes = ctypes.c_char_p(element_response_model.encode('utf-8'))
    
    jones = ((num_stations*num_times*num_freqs*num_dirs*4)*c_double_complex)()
    
    load_and_run_oskar_beam.argtypes = [c_char_p, c_char_p, 
                                        c_int, POINTER(c_int),
                                        c_int, c_double, c_double,
                                        POINTER(c_double), POINTER(c_double),
                                        c_int, POINTER(c_double),
                                        c_int, POINTER(c_double),
                                        c_bool, c_bool, c_bool, c_bool,
                                        POINTER(c_double_complex)]
    
    load_and_run_oskar_beam(ms_path_ctypes,
                            element_response_model_ctypes,
                            num_stations, station_idxs_ctypes,
                            num_dirs,
                            beam_ra0, beam_dec0,
                            ras_ctypes, decs_ctypes,
                            num_times, mjd_sec_times_ctypes,
                            num_freqs, freqs_ctypes,
                            apply_beam_norms, parallactic_rotate,
                            element_only, iau_order,
                            jones)
    
    jones_py = np.ctypeslib.as_array(jones, shape=(num_stations*num_times*num_freqs*num_dirs*4))
    jones_py = jones_py['real'] + 1j*jones_py['imag']
    
    jones_py = jones_py.reshape(num_stations, num_times, num_freqs, num_dirs, 2, 2)
    
    return jones_py

def run_mwa_beam(delays : np.ndarray,
                 coeff_path : str,
                 ras : np.ndarray, decs : np.ndarray,
                 j2000_lsts : np.ndarray, j2000_latitudes : np.ndarray,
                 freqs : np.ndarray,
                 amps : np.ndarray = np.ones(16),
                 apply_beam_norms : bool = False,
                 parallactic_rotate : bool = True,
                 iau_order : bool = True,
                 element_only : bool = False):
    """
    Run the MWA beam model using the EveryBeam library.
    This function is a wrapper around the C++ library `libuse_everybeam.so`.
    It takes the input parameters and converts them to the appropriate
    ctypes types, then calls the C++ function to calculate the beam response.
    The output is a numpy array of the beam response, with shape
    (num_stations, num_times, num_freqs, num_coords, 2, 2).
    
    Parameters
    ----------
    delays : np.ndarray
        Array of length 16 of MWA pointing delays
    coeff_path : str
        Path to the hdf5 FEE file.
    ras : np.ndarray
        Right ascensions of the coordinates in radians.
    decs : np.ndarray
        Declinations of the coordinates in radians.
    mjd_sec_times : np.ndarray
        Array of times in MJD seconds.
    j2000_lsts : np.ndarray
        Local sidereal times in J2000 coordinates.
    j2000_latitudes : np.ndarray
        Latitudes in J2000 coordinates.
    freqs : np.ndarray
        Array of frequencies to use.
    amps : np.ndarray
        Array of length 16 of MWA dipole amplitudes. Defaults to all ones.
    apply_beam_norms : bool
        Whether to apply beam normalisation. Defaults to False. The beam seems
        to be normalised by EveryBeam by default, so this is not ne eded.
    parallactic_rotate : bool
        Whether to apply parallactic angle rotation. Defaults to True.
    iau_order : bool
        If True, use IAU polarisation ordering, so set jones[0,0] to the NS dipole and jones[1,1] to EW. If False, jones[0,0] is EW.
    element_only : bool
        Whether to use only the element response. Defaults to False. Use this to
        look at the dipole response only, not the beam formed response. Does this by settings
        all amplitudes bar the first dipole, so defaults to the first dipole only.
    Returns
    -------
    np.ndarray
        The calculated Jones matrices with shape
        (num_stations, num_times, num_freqs, num_coords, 2, 2).
        Note that this is to have a consistent shape with other EveryBeam
        functions which can have varying station patterns; at the moment,
        num_stations is hard-coded to 1 for MWA.
    """
    
    num_stations = 1
    station_idxs = np.array([0])
    
    num_dirs = len(ras)
    num_freqs = len(freqs)
    num_times = len(j2000_lsts)
    
    azs = np.empty(num_dirs*num_times)
    zas = np.empty(num_dirs*num_times)
    para_angles = np.empty(num_dirs*num_times)
    
    for comp_ind in range(num_dirs):
        comp_has = j2000_lsts - ras[comp_ind]
        these_azs, these_els = erfa.hd2ae(comp_has, decs[comp_ind], j2000_latitudes)
        these_zas = np.pi/2 - these_els
        
        these_para_angles = erfa.hd2pa(comp_has, decs[comp_ind], j2000_latitudes)
        azs[comp_ind*num_times:(comp_ind+1)*num_times] = these_azs
        zas[comp_ind*num_times:(comp_ind+1)*num_times] = these_zas
        para_angles[comp_ind*num_times:(comp_ind+1)*num_times] = these_para_angles
        
    woden_path = importlib_resources.files(wodenpy).joinpath(f"libuse_everybeam.so")
    woden_lib = ctypes.cdll.LoadLibrary(woden_path)
    
    load_and_run_mwa_beam = woden_lib.load_and_run_mwa_beam
    
    zas_ctypes = (ctypes.c_double*(num_dirs*num_times))()
    azs_ctypes = (ctypes.c_double*(num_dirs*num_times))()
    para_angles_ctypes = (ctypes.c_double*(num_dirs*num_times))()
    
    for i in range(num_dirs*num_times):
        zas_ctypes[i] = zas[i]
        azs_ctypes[i] = azs[i]
        para_angles_ctypes[i] = para_angles[i]
    
    freqs_ctypes = (ctypes.c_double * num_freqs)()
    for i in range(num_freqs):
        freqs_ctypes[i] = freqs[i]
    coeff_path_ctypes = ctypes.c_char_p(coeff_path.encode('utf-8'))
    
    jones = ((num_stations*num_times*num_freqs*num_dirs*4)*c_double_complex)()
    
    delays_ctypes = (ctypes.c_double * 16)()
    for i in range(16):
        delays_ctypes[i] = delays[i]
        
    amps_ctypes = (ctypes.c_double * 16)()
    for i in range(16):
        amps_ctypes[i] = amps[i]
    
    load_and_run_mwa_beam.argtypes = [POINTER(c_double), POINTER(c_double), 
                                      c_char_p, 
                                      c_int, POINTER(c_double), POINTER(c_double),
                                      POINTER(c_double),
                                      c_int, POINTER(c_double),
                                      c_int,
                                      c_bool, c_bool,
                                      POINTER(c_double_complex)]
    
    load_and_run_mwa_beam(delays_ctypes, amps_ctypes,
                          coeff_path_ctypes,
                          num_dirs, azs_ctypes, zas_ctypes,
                          para_angles_ctypes,
                          num_freqs, freqs_ctypes,
                          num_times, 
                          parallactic_rotate, iau_order,
                          jones)
    
    jones_py = np.ctypeslib.as_array(jones, shape=(num_stations*num_times*num_freqs*num_dirs*4))
    jones_py = jones_py['real'] + 1j*jones_py['imag']
    
    jones_py = jones_py.reshape(num_stations, num_times, num_freqs, num_dirs, 2, 2)
    
    return jones_py