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


from wodenpy.wodenpy_setup.run_setup import check_for_library
have_everybeam = check_for_library('everybeam')

class EB_fake:
    """
        A fake `everybeam` class so we can build the documentation online
        in ReadTheDocs without installing `everybeam`, which is non-trivial
    """
    def __init__(self):
        """
        Just set everything that is ever used to None
        """
        self.OSKAR = None
        self.LOFAR = None
        self.MWA = None
        self.MWALocal = None
        self.load_telescope = None
        self.Telescope = None
        
if have_everybeam:
    import everybeam as eb
else:
    eb = EB_fake()

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
        time float astropy time instance

    Returns:
        pointing_xyz (ndarray): NumPy array containing the ITRS X, Y and Z coordinates
    """
    
    coord = SkyCoord(ra*u.rad, dec*u.rad, frame='icrs')
    coord_itrs = coord.transform_to(ITRS(obstime=time))
    
    return np.asarray(coord_itrs.cartesian.xyz.transpose())

def load_OSKAR_telescope(ms_path : str, response_model = "skala40_wave",
                         use_differential_beam : bool = True) -> eb.OSKAR: # type: ignore
    """Load an OSKAR telescope from a measurement set.

    Parameters
    ----------
    ms_path : str
        Path to the measurement set
    response_model : str, optional
        Response model to use, by default "skala40_wave"
    use_differential_beam : bool, optional
        Use the differential beam a.k.a return a "normalised" beam, by default True

    Returns
    -------
    eb.OSKAR
        Telescope object
    """

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


def load_LOFAR_telescope(ms_path : str, response_model : str = "hamaker",
                         use_differential_beam : bool = False) -> eb.LOFAR: # type: ignore
    """Load an LOFAR telescope from a measurement set. Settings lifted
    directly from https://everybeam.readthedocs.io/en/latest/tree/demos/lofar-lobes.html

    Parameters
    ----------
    ms_path : str
        Path to the measurement set
    response_model : str, optional
        Response model to use, by default "hamaker"
    use_differential_beam : bool, optional
        Use the differential beam a.k.a return a "normalised" beam, by default False

    Returns
    -------
    eb.LOFAR
        Telescope object
    """

    telescope = eb.load_telescope(ms_path,
                                  use_differential_beam=use_differential_beam,
                                  element_response_model=response_model)
    
    if type(telescope) != eb.LOFAR:
        print(f'WARNING: Telescope specified in {ms_path} is not an OSKAR telescope. Proceeding, but you might get nonsense results.')
    
    return telescope


def load_MWA_telescope(ms_path : str, coeff_path : str,
                       use_local_mwa : bool = True) -> Union[eb.MWA, eb.MWALocal]: # type: ignore
    """Load an MWA telescope from a measurement set.

    Parameters
    ----------
    ms_path : str
        Path to the measurement set
    response_model : str, optional
        Response model to use, by default "lobes"
    use_local_mwa : bool, optional
        Use the local MWA model, which takes za/az instead of RA/Dec. Defaults
        to True

    Returns
    -------
    Union[eb.MWA, eb.MWALocal]
        Telescope object, either MWA (if use_local_mwa is False) or MWALocal
        (if use_local_mwa is True)
    """

    ## Load the telescope. Adding use_differential_beam seems to do nothing,
    ## so always leave it as False
    telescope = eb.load_telescope(ms_path,
                                  use_differential_beam=False,
                                  coeff_path=coeff_path,
                                  use_local_mwa=use_local_mwa)
    
    # assert type(telescope) == eb.MWA
    if type(telescope) != eb.MWA and type(telescope) != eb.MWALocal:
        print(f'WARNING: Telescope specified in {ms_path} is not an MWA telescope. Proceeding, but you might get nonsense results.')
    
    return telescope


def eb_local_xyz_from_radec(ra : float, dec : float,
                            altaz_frame : astropy.coordinates.AltAz,
                            delta_az : float = (1/2)*np.pi,
                            negative_azimuth=True):
    """
    Get the local cartesian coordinates used by EveryBeam from a given RA, Dec, and AltAz frame.
    This function transforms the given right ascension (RA) and declination (Dec) into local 
    cartesian coordinates based on the provided AltAz frame. The azimuth is adjusted by reversing 
    it and adding 90 degrees (or a specified delta) to match the expected coordinates.
    
    This function is used to emulate the way the EveryBeam does it's parallactic
    rotation. `delta_az` is set to π/2 by experiment to match outputs from EveryBeam.
    Not really needed, as when can use EveryBeam to do the parallactic rotation
    for everything aside the MWA beam.
    
    
    Parameters
    -----------
    ra : float
        Right ascension in radians.
    dec : float
        Declination in radians.
    altaz_frame : `astropy.coordinates.AltAz`
        The AltAz frame to transform the coordinates into.
    delta_az : float, optional
        The azimuth adjustment in radians. Default is π/2.
    negative_azimuth : bool, optional
        If True, the azimuth is reversed before adding the delta. Default is True.
        
    Returns
    --------
    numpy.ndarray
        A 3xN array of the local cartesian coordinates.
    """
    
    coord = SkyCoord(ra=ra*u.rad, dec=dec*u.rad, frame='icrs')
    coord = coord.transform_to(altaz_frame)
    
    delta_az = delta_az * u.rad
    
    if negative_azimuth:
        updated_coord = SkyCoord(az=-coord.az + delta_az,
                alt=coord.alt,
                distance=coord.distance,
                frame=coord.frame)
    else:
        updated_coord = SkyCoord(az=coord.az + delta_az,
                                 alt=coord.alt,
                                 distance=coord.distance,
                                 frame=coord.frame)
        
    return np.array(updated_coord.cartesian.xyz.transpose())

def eb_north_east(direction : np.ndarray, ncp_t : np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute the north and east vectors given a direction ITRF vector,
    and the ITRF vector towards the north celestial pole.
    This function calculates the east vector by taking the cross product of the normal vector (ncp_t) 
    and the direction vector, and then normalizes it. The north vector is then calculated as the cross 
    product of the direction vector and the east vector.
    
    Translated from EveryBeam station.cc Station::ComputeElementResponse
    
    const vector3r_t east = normalize(cross(ncp_t, direction));
    const vector3r_t north = cross(direction, east);
    options.east = east;
    options.north = north;
    
    Parameters
    ------------
    direction : np.ndarray
        A 3-element array representing the direction vector.
    ncp_t : np.ndarray
        A 3-element array representing the normal vector.
        
    Returns
    --------
    Tuple[np.ndarray, np.ndarray]
        A tuple containing the north and east vectors as 3-element arrays.
    """
    east = np.cross(ncp_t, direction)
    east = east/np.linalg.norm(east)
    north = np.cross(direction, east)
    
    return north, east

def calc_everybeam_rotation(direction : np.ndarray, north : np.ndarray,
                            east : np.ndarray) -> np.ndarray:
    """Given an ITRF 3-vector in the direction of interest `direction`,
    and the associated north and east vectors, calculate the 2x2 rotation 
    matrix to rotate by parallactic angle.
    
    Translated from EveryBeam beamformer.cc BeamFormer::LocalResponse
    const vector3r_t e_phi = normalize(cross(direction));
    const vector3r_t e_theta = cross(e_phi, direction);
    result *= {dot(e_theta, options.north), dot(e_theta, options.east),
               dot(e_phi, options.north), dot(e_phi, options.east)};
    
    Parameters
    ------------
    direction : np.ndarray
        A 3-element array representing the direction ITRF vector.
    north : np.ndarray
        A 3-element array representing the north vector.
    east : np.ndarray
        A 3-element array representing the east vector.
        
    Returns
    --------
    np.ndarray
        A 2x2 rotation matrix.
    """
    e_phi = np.cross([0.0, 0.0, 1.0], direction)
    e_phi = e_phi/np.linalg.norm(e_phi)
    
    e_theta = np.cross(e_phi, direction)
    e_theta = e_theta/np.linalg.norm(e_theta)
    
    rot_matrix = np.array([[np.dot(e_theta, north), np.dot(e_theta, east)],
                            [np.dot(e_phi, north), np.dot(e_phi, east)]])
    
    return rot_matrix

# @profile
def run_everybeam(ras: np.ndarray, decs: np.ndarray,
                  beam_ra0: float, beam_dec0: float,
                  j2000_latitudes: np.ndarray, j2000_lsts: np.ndarray,
                  current_latitude: float, current_longitude: float,
                  times: np.ndarray, freqs: np.ndarray,
                  telescope: eb.Telescope, # type: ignore
                  station_ids: np.ndarray,
                  full_accuracy: bool = True,
                  apply_beam_norms: bool = True,
                  reorder_jones: bool = False,
                  element_only: bool = False,
                  eb_rotate: bool = False,
                  parallactic_rotate: bool = False,
                  para_angle_offset: float = 0) -> np.ndarray:
    """
    Calculate the Jones matrices for a given set of coordinates, times,
    frequencies, and station ids using the EveryBeam library.
    `j2000_latitudes` should be the array latitude as precessed back to J2000,
    with `j2000_lsts` being the matching LST in J2000. `current_latitude` and
    `current_longitude` should be latitude and longitude of the array at the
    time of the observation. `telescope` should be an EveryBeam telescope object.
    
    
    Parameters
    ------------
    ras : np.ndarray
        Right ascensions of the coordinates in radians.
    decs : np.ndarray
         Declinations of the coordinates in radians.
    beam_ra0 : float
        Right ascension of the beam center in radians.
    beam_dec0 : float
        Declination of the beam center in radians.
    j2000_latitudes : np.ndarray
        Latitudes in J2000 coordinates.
    j2000_lsts : np.ndarray
        Local sidereal times in J2000 coordinates.
    current_latitude : float
        Current latitude in radians.
    current_longitude : float
        Current longitude in radians.
    times : np.ndarray
        Array of observation times.
    freqs : np.ndarray
        Array of frequencies.
    telescope : (eb.Telescope):
        Telescope object from the EveryBeam library.
    station_ids : np.ndarray
        Array of station IDs.
    full_accuracy : bool, optional
        Whether to use the full accuracy of the EveryBeam library. Defaults to True.
        If False, the array factor and element response are calculated separately
        (the two values that multiply to give the full response). The array factor
        is only calculated at the middle time and frequency, under the assumption
        that the range of frequencies and times is small enough that the array factor
        will not change significantly. Magnitude of the differences vary by
        frequency and direction so use with caution.
    apply_beam_norms : bool, optional
        Whether to apply beam normalisation. Defaults to True. Achieved by
        calculating the beam response at beam centre, and multiplying all
        Jones by the inverse of this central beam response.
    reorder_jones : bool, optional
        Whether to reorder the Jones matrices. Defaults to False. Just rearranges
        the Jones matrix from [[0,0, 0,1,], [1,0, 1,1]] to [[1,1, 1,0,], [0,1, 0,0]].
    element_only : bool, optional
        Whether to use only the element response. Defaults to False. Use this to
        look at the dipole response only, not the beam formed response.
    eb_rotate : bool, optional
        Whether to apply parallactic rotation using EveryBeam. Defaults to False.
        Should probably be used for everything apart from MWA beams.
    parallactic_rotate : bool, optional
        Whether to apply parallactic angle rotation using `wodenpy`. Defaults to False.
        Should be True for MWA beams if you want rotation. If True for a non-MWA beam,
        `wodenpy` should match the output as if `eb_rotate` was True.
    para_angle_offset : float, optional
        Offset to add to the parallactic angle. Defaults to 0.
        
    Returns
    --------
    np.ndarray
        The calculated Jones matrices with shape (num_stations, num_times, num_freqs, num_coords, 2, 2).
    """
    
    # print("YO MR WHITE", type(telescope))
    
    num_stations = len(station_ids)
    num_times = len(times)
    num_freqs = len(freqs)
    num_coords = len(ras)
    
    all_output_jones = np.zeros((num_stations, num_times, num_freqs, num_coords, 2, 2), dtype=np.complex128)*np.nan
    
    non_itrf_beams = [eb.MWA, eb.MWALocal]
    
    if not full_accuracy:
        mid_time = times[len(times)//2]
        mid_freq = freqs[len(freqs)//2]
        
        phase_itrf_mid = radec_to_xyz(beam_ra0, beam_dec0, mid_time)
        dir_itrfs_mid = radec_to_xyz(ras, decs, mid_time)
        
        array_factors = telescope.array_factor(mid_time.mjd*3600*24, station_ids,
                                            mid_freq, dir_itrfs_mid, phase_itrf_mid)
    
    if parallactic_rotate:
        if type(telescope) not in non_itrf_beams:
            coords = SkyCoord(ras*u.rad, decs*u.rad, frame='icrs')
            location = EarthLocation(lat=current_latitude*u.rad,
                                    lon=current_longitude*u.rad)
    
    for time_ind, time in enumerate(times):
        if type(telescope) in non_itrf_beams:
            comp_has = j2000_lsts[time_ind] - ras
            azs, els = erfa.hd2ae(comp_has, decs, j2000_latitudes[time_ind])
            zas = np.pi/2 - els
            
            if parallactic_rotate:
                beam_ha0 = j2000_lsts[time_ind] - beam_ra0
                beam_az0, beam_el0 = erfa.hd2ae(beam_ha0, beam_dec0,
                                                j2000_latitudes[time_ind])
                beam_za0 = np.pi/2 - beam_el0
        else:
            phase_itrf = radec_to_xyz(beam_ra0, beam_dec0, time)
            dir_itrfs = radec_to_xyz(ras, decs, time)
            
            if parallactic_rotate:
                altaz_frame = AltAz(obstime=time, location=location)
                ncp_t = eb_local_xyz_from_radec(0, np.radians(90), altaz_frame)
                dir_local = eb_local_xyz_from_radec(ras, decs, altaz_frame)
                
        
        time_mjd_secs = time.mjd*3600*24
        
        if parallactic_rotate:
            has = j2000_lsts[time_ind] - ras
            para_angles = erfa.hd2pa(has, decs, j2000_latitudes[time_ind])
            
            rot_matrix = np.empty((num_coords, 2,2))
            
            if type(telescope) in non_itrf_beams:
                rot_matrix[:,0,0] = np.sin(-para_angles)
                rot_matrix[:,0,1] = -np.cos(-para_angles)
                rot_matrix[:,1,0] = -np.cos(-para_angles)
                rot_matrix[:,1,1] = -np.sin(-para_angles)
                
            else:
                for dir_ind, dir_itrf in enumerate(dir_itrfs):
                    
                    dir_az = dir_local[dir_ind]
                    north, east = eb_north_east(dir_az, ncp_t)
                    rot = calc_everybeam_rotation(dir_az, north, east)
                    rot_matrix[dir_ind] = rot
                
        for station_ind, station_id in enumerate(station_ids):
            for freq_ind, freq in enumerate(freqs):
                
                if not full_accuracy:
                    element_responses = np.zeros((num_coords, 2, 2), dtype=np.complex128)*np.nan
                
                if apply_beam_norms:
                    if type(telescope) == eb.MWA:
                        ##Get the response
                        norm_jones = telescope.station_response(time_mjd_secs, station_id, freq,
                                                                beam_ra0, beam_dec0)
                    elif type(telescope) == eb.MWALocal:
                        norm_jones = telescope.station_response(time_mjd_secs, station_id, freq,
                                                                beam_az0, beam_za0)
                        print("norm_jones", norm_jones)
                        
                    else:
                        element_id = 0
                        if element_only:
                            norm_jones = telescope.element_response(time_mjd_secs, station_id, element_id, freq,
                                                                phase_itrf, rotate=eb_rotate)
                        else:
                            if full_accuracy:
                                norm_jones = telescope.station_response(time_mjd_secs, station_id, freq,
                                                            phase_itrf, phase_itrf, 
                                                            rotate=eb_rotate)
                            else:
                                array_factor_norm = telescope.array_factor(mid_time.mjd*3600*24, station_id,
                                            mid_freq, phase_itrf_mid, phase_itrf_mid)
                                element_norm = telescope.element_response(time_mjd_secs, station_id, element_id, freq,
                                                                phase_itrf, rotate=eb_rotate)
                                norm_jones = np.matmul(array_factor_norm, element_norm)
                        if parallactic_rotate:
                            dir_phase_local = eb_local_xyz_from_radec(beam_ra0, beam_dec0, altaz_frame)
                            north, east = eb_north_east(dir_phase_local, ncp_t)
                            rot = calc_everybeam_rotation(dir_phase_local, north, east)
               
                    if parallactic_rotate:
                        if type(telescope) in non_itrf_beams:
                            ha0 = j2000_lsts[time_ind] - beam_ra0
                            para_angles = erfa.hd2pa(ha0, beam_dec0, j2000_latitudes[time_ind])
                            rot = np.empty((2,2))
                            rot[0,0] = np.sin(-para_angles)
                            rot[0,1] = -np.cos(-para_angles)
                            rot[1,0] = -np.cos(-para_angles)
                            rot[1,1] = -np.sin(-para_angles)
                            # print("rot", rot)
                        
                        norm_jones = np.matmul(norm_jones, rot)
                    
                for coord_ind, (ra, dec) in enumerate(zip(ras, decs)):
                    ##Only MWA uses ra,dec as a direct input
                    if type(telescope) == eb.MWA:
                            response = telescope.station_response(time_mjd_secs, station_id, freq,
                                                                  ra, dec)
                            all_output_jones[station_ind, time_ind, freq_ind, coord_ind] = response
                    ##Only MWALocal uses az,za as a direct input
                    elif type(telescope) == eb.MWALocal:
                            response = telescope.station_response(time_mjd_secs, station_id, freq,
                                                                  zas[coord_ind], azs[coord_ind])
                            all_output_jones[station_ind, time_ind, freq_ind, coord_ind] = response
                    ##Everything else uses ITRF coordinates
                    else:
                        if element_only:
                            response = telescope.element_response(time_mjd_secs, station_id, freq,
                                                                dir_itrfs[coord_ind], rotate=eb_rotate)
                        else:
                            if full_accuracy:
                                response = telescope.station_response(time_mjd_secs, station_id, freq,
                                                                dir_itrfs[coord_ind], phase_itrf, 
                                                                rotate=eb_rotate)
                                all_output_jones[station_ind, time_ind, freq_ind, coord_ind] = response
                            else:
                                element_id = 0
                                element_response = telescope.element_response(time_mjd_secs,
                                                                station_id, element_id, freq,
                                                                dir_itrfs[coord_ind],
                                                                rotate=eb_rotate)
                                element_responses[coord_ind] = element_response
                                
                ##Parallactic angle doesn't change per station or freq, but
                ##if we are normalising the beam, we want to rotate before we normalise
                ##So do the rotation now. This also means if we are calculating array
                ## factor and element response separately, need to multiply them together here
                if not full_accuracy:
                    all_output_jones[station_ind, time_ind, freq_ind, :, :, :] = np.einsum('klm,kmn->kln', array_factors[station_ind, :, :, :], element_responses)
                
                if parallactic_rotate:
                    
                    rot_jones = np.einsum('klm,kmn->kln', all_output_jones[station_ind, time_ind, freq_ind, :, :, :], rot_matrix)
                    all_output_jones[station_ind, time_ind, freq_ind, :, :, :] = rot_jones
                    
                if apply_beam_norms:
                    ##Each station, time, and freq gets it's own normalisation
                    ##Same 2x2 normalisation for all directions
                    inv_beam_norms = np.linalg.inv(norm_jones)
                    output_jones = np.einsum('lm,kmn->kln', inv_beam_norms, all_output_jones[station_ind, time_ind, freq_ind, :, :, :])
                    all_output_jones[station_ind, time_ind, freq_ind, :, :, :] = output_jones
                    
    if reorder_jones:
        ##swap all_output_jones[:,:,:,:,0,0] with all_output_jones[:,:,:,:,1,1]
        all_output_jones[:, :, :, :, [0, 1], [0, 1]] = all_output_jones[:, :, :, :, [1, 0], [1, 0]]
        ##swap all_output_jones[:,:,:,:,0,1] with all_output_jones[:,:,:,:,1,0]
        all_output_jones[:, :, :, :, [0, 1], [1, 0]] = all_output_jones[:, :, :, :, [1, 0], [0, 1]]
                    
    return all_output_jones


def run_everybeam_thread(num_threads : int, thread_id : int,
                         ms_path : str, coeff_path : str,
                         ras : np.ndarray, decs : np.ndarray,
                         ra0 : float, dec0 : float,
                         j2000_latitudes : np.ndarray, j2000_lsts : np.ndarray,
                         current_latitude : float, current_longitude : float,
                         times : np.ndarray, freqs : np.ndarray,
                         station_ids : np.ndarray,
                         full_accuracy : bool = True,
                         use_differential_beam : bool = True,
                         apply_beam_norms : bool = True,
                         reorder_jones : bool = False,
                         element_only : bool = False,
                         eb_rotate : bool = False,
                         parallactic_rotate : bool = True,
                         para_angle_offset : float = 0,
                         element_response_model='hamaker',
                         use_local_mwa : bool = True) -> Tuple[np.ndarray, int]:
    """
    Thread function called by `run_everybeam_over_threads` to calculate the
    EveryBeam response in parrallel. Calls `run_everybeam` with a subset of
    the coordinates; see `run_everybeam` for more details of the parameters.
    
    Creates a new EveryBeam telescope object from `ms_path` for each thread.
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
    ms_path : str
        Path to the measurement set to load the EveryBeam telescope from.
    ras : np.ndarray
        Right ascensions of the coordinates in radians.
    decs : np.ndarray
         Declinations of the coordinates in radians.
    beam_ra0 : float
        Right ascension of the beam center in radians.
    beam_dec0 : float
        Declination of the beam center in radians.
    j2000_latitudes : np.ndarray
        Latitudes in J2000 coordinates.
    j2000_lsts : np.ndarray
        Local sidereal times in J2000 coordinates.
    current_latitude : float
        Current latitude in radians.
    current_longitude : float
        Current longitude in radians.
    times : np.ndarray
        Array of observation times.
    freqs : np.ndarray
        Array of frequencies.
    station_ids : np.ndarray
        Array of station IDs.
    full_accuracy : bool, optional
        Whether to use the full accuracy of the EveryBeam library. Defaults to True.
        If False, the array factor and element response are calculated separately
        (the two values that multiply to give the full response). The array factor
        is only calculated at the middle time and frequency, under the assumption
        that the range of frequencies and times is small enough that the array factor
        will not change significantly. Magnitude of the differences vary by
        frequency and direction so use with caution.
    apply_beam_norms : bool, optional
        Whether to apply beam normalisation. Defaults to True. Achieved by
        calculating the beam response at beam centre, and multiplying all
        Jones by the inverse of this central beam response.
    reorder_jones : bool, optional
        Whether to reorder the Jones matrices. Defaults to False. Just rearranges
        the Jones matrix from [[0,0, 0,1,], [1,0, 1,1]] to [[1,1, 1,0,], [0,1, 0,0]].
    element_only : bool, optional
        Whether to use only the element response. Defaults to False. Use this to
        look at the dipole response only, not the beam formed response.
    eb_rotate : bool, optional
        Whether to apply parallactic rotation using EveryBeam. Defaults to False.
        Should probably be used for everything apart from MWA beams.
    parallactic_rotate : bool, optional
        Whether to apply parallactic angle rotation using `wodenpy`. Defaults to False.
        Should be True for MWA beams if you want rotation. If True for a non-MWA beam,
        `wodenpy` should match the output as if `eb_rotate` was True.
    para_angle_offset : float, optional
        Offset to add to the parallactic angle. Defaults to 0.
    element_response_model : str, optional
        The Everybeam element response model to use. Defaults to 'hamaker'.
        Avaible options are 'hamaker' (LOFAR), 'skala40_wave' (OSKAR), and 'MWA' (MWA).
    use_local_mwa : bool, optional
        Whether to use the local MWA model. Defaults to True. The local MWA model
        takes za/az instead of RA/Dec.
        
    Returns
    --------
    Tuple[np.ndarray, int]
        The calculated Jones matrices with shape
        (num_stations, num_times, num_freqs, num_coords_in_thread, 2, 2), as
        well as the thread ID. Use the thread ID to insert this thread output
        into the correct place in the final Jones matrix.
    """
    
    telescope = eb.load_telescope(ms_path,
                                  use_differential_beam=use_differential_beam,
                                  coeff_path=coeff_path,
                                  element_response_model=element_response_model,
                                  use_local_mwa=use_local_mwa)
    
    num_coords = len(ras)
    coords_per_thread = int(np.ceil(num_coords / num_threads))
    
    low_coord = thread_id * coords_per_thread
    high_coord = (thread_id + 1) * coords_per_thread
    
    print(f"Thread {thread_id} processing coords {low_coord} to {high_coord}")
    
    jones = run_everybeam(ras[low_coord:high_coord],
                          decs[low_coord:high_coord],
                          ra0, dec0,
                          j2000_latitudes, j2000_lsts,
                          current_latitude, current_longitude,
                          times, freqs,
                          telescope, station_ids,
                          full_accuracy=full_accuracy,
                          apply_beam_norms=apply_beam_norms,
                          reorder_jones=reorder_jones,
                          element_only=element_only,
                          eb_rotate=eb_rotate,
                          parallactic_rotate=parallactic_rotate)
    
    print(f"Thread {thread_id} finished")
    
    return jones, thread_id

def run_everybeam_over_threads(num_threads : int,
                               ms_path : str,
                               coeff_path : str,
                               ras : np.ndarray, decs : np.ndarray,
                               ra0 : float, dec0 : float,
                               j2000_latitudes : np.ndarray, j2000_lsts : np.ndarray,
                               current_latitude : float, current_longitude : float,
                               times : np.ndarray, freqs : np.ndarray,
                               station_ids : np.ndarray,
                               full_accuracy : bool = True,
                               use_differential_beam : bool = False,
                               apply_beam_norms : bool = True,
                               reorder_jones : bool = False,
                               element_only : bool = False,
                               eb_rotate : bool = False,
                               parallactic_rotate : bool = True,
                               use_local_mwa : bool = True,
                               para_angle_offset : float = 0,
                               element_response_model='hamaker'):
    """
    Runs `run_everybeam` in parallel over `num_threads` threads, using
    `concurrent.futures.ProcessPoolExecutor`. See `run_everybeam` for more
    details of what each parameter does.
    
    Creates a new EveryBeam telescope object from `ms_path` for each thread.
    This has to be done because `concurrent.futures.ProcessPoolExecutor` has
    to pickle the function and all it's arguments, and EveryBeam objects can't
    be pickled. This is somewhat wasteful but I can't work out a better way
    to make things parallel.
    
    Parameters
    ------------
    num_threads : int
        Number of threads being in call by `run_everybeam_over_threads`.
    ms_path : str
        Path to the measurement set to load the EveryBeam telescope from.
    ras : np.ndarray
        Right ascensions of the coordinates in radians.
    decs : np.ndarray
         Declinations of the coordinates in radians.
    beam_ra0 : float
        Right ascension of the beam center in radians.
    beam_dec0 : float
        Declination of the beam center in radians.
    j2000_latitudes : np.ndarray
        Latitudes in J2000 coordinates.
    j2000_lsts : np.ndarray
        Local sidereal times in J2000 coordinates.
    current_latitude : float
        Current latitude in radians.
    current_longitude : float
        Current longitude in radians.
    times : np.ndarray
        Array of observation times.
    freqs : np.ndarray
        Array of frequencies.
    station_ids : np.ndarray
        Array of station IDs.
    full_accuracy : bool, optional
        Whether to use the full accuracy of the EveryBeam library. Defaults to True.
        If False, the array factor and element response are calculated separately
        (the two values that multiply to give the full response). The array factor
        is only calculated at the middle time and frequency, under the assumption
        that the range of frequencies and times is small enough that the array factor
        will not change significantly. Magnitude of the differences vary by
        frequency and direction so use with caution.
    use_differential_beam : bool, optional
        Whether to use the differential beam a.k.a return a "normalised" beam, 
        as normalised by EveryBeam. by default False
    apply_beam_norms : bool, optional
        Whether to apply beam normalisation. Defaults to True. Achieved by
        calculating the beam response at beam centre, and multiplying all
        Jones by the inverse of this central beam response.
    reorder_jones : bool, optional
        Whether to reorder the Jones matrices. Defaults to False. Just rearranges
        the Jones matrix from [[0,0, 0,1,], [1,0, 1,1]] to [[1,1, 1,0,], [0,1, 0,0]].
    element_only : bool, optional
        Whether to use only the element response. Defaults to False. Use this to
        look at the dipole response only, not the beam formed response.
    eb_rotate : bool, optional
        Whether to apply parallactic rotation using EveryBeam. Defaults to False.
        Should probably be used for everything apart from MWA beams.
    parallactic_rotate : bool, optional
        Whether to apply parallactic angle rotation using `wodenpy`. Defaults to False.
        Should be True for MWA beams if you want rotation. If True for a non-MWA beam,
        `wodenpy` should match the output as if `eb_rotate` was True.
    para_angle_offset : float, optional
        Offset to add to the parallactic angle. Defaults to 0.
    element_response_model : str, optional
        The Everybeam element response model to use. Defaults to 'hamaker'.
        Avaible options are 'hamaker' (LOFAR), 'skala40_wave' (OSKAR), and 'MWA' (MWA).
    use_local_mwa : bool, optional
        Whether to use the local MWA model. Defaults to True. The local MWA model
        takes za/az instead of RA/Dec.
        
    Returns
    --------
    np.ndarray
        The calculated Jones matrices with shape
        (num_stations, num_times, num_freqs, num_coord, 2, 2)
    """
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
        
            future_data = [executor.submit(run_everybeam_thread,
                                           num_threads, thread_id,
                                           ms_path, 
                                           coeff_path,
                                           ras, decs,
                                           ra0, dec0,
                                           j2000_latitudes, j2000_lsts,
                                           current_latitude, current_longitude,
                                           times, freqs,
                                           station_ids,
                                           full_accuracy=full_accuracy,
                                           use_differential_beam=use_differential_beam,
                                           apply_beam_norms=apply_beam_norms,
                                           reorder_jones=reorder_jones,
                                           element_only=element_only,
                                           eb_rotate=eb_rotate,
                                           parallactic_rotate=parallactic_rotate,
                                           para_angle_offset=para_angle_offset,
                                           element_response_model=element_response_model,
                                           use_local_mwa=use_local_mwa)
                                    for thread_id in range(num_threads)]
            
            all_jones_chunks = []
            all_thread_ids = []
            for future in concurrent.futures.as_completed(future_data):
                jones_chunk, thread_id = future.result()
                all_jones_chunks.append(jones_chunk)
                all_thread_ids.append(thread_id)
                
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