import numpy as np
from astropy.coordinates import ITRS, SkyCoord, AltAz, EarthLocation
from astropy.time import Time, TimeDelta
import astropy.units as u
import argparse
from wodenpy.use_libwoden.create_woden_struct_classes import Woden_Struct_Classes
import erfa
from typing import Union
import concurrent.futures
from line_profiler import profile
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
            self.MWALocal = None
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


def eb_local_xyz_from_radec(ra, dec, altaz_frame, delta_az=(1/2)*np.pi, 
                            negative_azimuth=True):
    """Get the local cartesian coords used by EveryBeam from a given RA, Dec,
    and AltAz frame. Reversing the azimuth and adding 90 degrees was found
    to match via trial and error"""
    
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

def eb_north_east(direction, ncp_t):
    ##taken from EveryBeam station.cc Station::ComputeElementResponse
    
    # const vector3r_t east = normalize(cross(ncp_t, direction));
    # const vector3r_t north = cross(direction, east);
    # options.east = east;
    # options.north = north;
    
    east = np.cross(ncp_t, direction)
    east = east/np.linalg.norm(east)
    north = np.cross(direction, east)
    
    
    return north, east

def calc_everybeam_rotation(direction, north, east):
    
    ##taken from EveryBeam beamformer.cc BeamFormer::LocalResponse
    # const vector3r_t e_phi = normalize(cross(direction));
    # const vector3r_t e_theta = cross(e_phi, direction);
    # result *= {dot(e_theta, options.north), dot(e_theta, options.east),
    #            dot(e_phi, options.north), dot(e_phi, options.east)};
    
    e_phi = np.cross([0.0, 0.0, 1.0], direction)
    e_phi = e_phi/np.linalg.norm(e_phi)
    
    e_theta = np.cross(e_phi, direction)
    e_theta = e_theta/np.linalg.norm(e_theta)
    
    rot_matrix = np.array([[np.dot(e_theta, north), np.dot(e_theta, east)],
                            [np.dot(e_phi, north), np.dot(e_phi, east)]])
    
    return rot_matrix

# @profile
def run_everybeam(ras : np.ndarray, decs : np.ndarray,
                  beam_ra0 : float, beam_dec0 : float,
                  j2000_latitudes : np.ndarray, j2000_lsts : np.ndarray,
                  current_latitude : float, current_longitude : float,
                  times : np.ndarray, freqs : np.ndarray,
                  telescope: eb.Telescope,  # type: ignore
                  station_ids : np.ndarray,
                  apply_beam_norms : bool = True,
                  reorder_jones : bool = False,
                  element_only : bool = False,
                  eb_rotate : bool = False,
                  parallactic_rotate : bool = False,
                  para_angle_offset : float = -np.pi/2) -> np.ndarray:
    
    
    num_stations = len(station_ids)
    num_times = len(times)
    num_freqs = len(freqs)
    num_coords = len(ras)
    
    all_output_jones = np.zeros((num_stations, num_times, num_freqs, num_coords, 2, 2), dtype=np.complex128)*np.nan
    
    if parallactic_rotate:
        if type(telescope) != eb.MWA or type(telescope) != eb.MWALocal:
            coords = SkyCoord(ras*u.rad, decs*u.rad, frame='icrs')
            location = EarthLocation(lat=current_latitude*u.rad,
                                    lon=current_longitude*u.rad)
    
    for time_ind, time in enumerate(times):
        
        if type(telescope) == eb.MWA or type(telescope) == eb.MWALocal:
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
                ncp_t = eb_local_xyz_from_radec(0, np.radians(90), altaz_frame)
                dir_local = eb_local_xyz_from_radec(ras, decs, altaz_frame)
                altaz_frame = AltAz(obstime=time, location=location)
        
        time_mjd_secs = time.mjd*3600*24
        
        if parallactic_rotate:
            has = j2000_lsts[time_ind] - ras
            para_angles = erfa.hd2pa(has, decs, j2000_latitudes[time_ind])
            
            rot_matrix = np.empty((num_coords, 2,2))
            
            if type(telescope) == eb.MWA or type(telescope) == eb.MWALocal:
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
                
                if apply_beam_norms:
                    if type(telescope) == eb.MWA:
                        ##Get the response
                        norm_jones = telescope.station_response(time_mjd_secs, station_id, freq,
                                                                beam_ra0, beam_dec0)
                    elif type(telescope) == eb.MWALocal:
                        norm_jones = telescope.station_response(time_mjd_secs, station_id, freq,
                                                                beam_az0, beam_za0)
                        
                    if type(telescope) == eb.MWA or type(telescope) == eb.MWALocal:
                        if parallactic_rotate:
                            ha0 = j2000_lsts[time_ind] - beam_ra0
                            para_angles = erfa.hd2pa(ha0, beam_dec0, j2000_latitudes[time_ind])
                            rot = np.empty((2,2))
                            rot[0,0] = np.sin(-para_angles)
                            rot[0,1] = -np.cos(-para_angles)
                            rot[1,0] = -np.cos(-para_angles)
                            rot[1,1] = -np.sin(-para_angles)
                        
                    else:
                        element_id = 0
                        if element_only:
                            norm_jones = telescope.element_response(time_mjd_secs, station_id, element_id, freq,
                                                                phase_itrf, rotate=eb_rotate)
                        else:
                            norm_jones = telescope.station_response(time_mjd_secs, station_id, freq,
                                                            phase_itrf, phase_itrf, 
                                                            rotate=eb_rotate)
                        if parallactic_rotate:
                            dir_phase_local = eb_local_xyz_from_radec(beam_ra0, beam_dec0, altaz_frame)
                            north, east = eb_north_east(dir_phase_local, ncp_t)
                            rot = calc_everybeam_rotation(dir_phase_local, north, east)
               
                    if parallactic_rotate:
                        norm_jones = np.matmul(norm_jones, rot)
                    
                for coord_ind, (ra, dec) in enumerate(zip(ras, decs)):
                    # if np.isnan(ra) or np.isnan(dec):
                    #     pass
                    # else:
                    
                    if type(telescope) == eb.MWA:
                            ##Get the response
                            print("WE BE DOING THE THINGS")
                            response = telescope.station_response(time_mjd_secs, station_id, freq,
                                                                  ra, dec)
                    elif type(telescope) == eb.MWALocal:
                            ##Get the response
                            response = telescope.station_response(time_mjd_secs, station_id, freq,
                                                                  zas[coord_ind], azs[coord_ind])
                            
                    else:
                        if element_only:
                            response = telescope.element_response(time_mjd_secs, station_id, 0, freq,
                                                                dir_itrfs[coord_ind], rotate=eb_rotate)
                        else:
                            response = telescope.station_response(time_mjd_secs, station_id, freq,
                                                            dir_itrfs[coord_ind], phase_itrf, 
                                                            rotate=eb_rotate)
                            
                    all_output_jones[station_ind, time_ind, freq_ind, coord_ind] = response
                
                if parallactic_rotate:
                    ##Parallactic angle doesn't change per station or freq, but
                    ##if we are normalising the beam, we want to rotate before we normalise
                    ##So do the rotation now
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
                         use_differential_beam : bool = True,
                         apply_beam_norms : bool = True,
                         reorder_jones : bool = False,
                         element_only : bool = False,
                         eb_rotate : bool = False,
                         parallactic_rotate : bool = True,
                         para_angle_offset : float = 0,
                         element_response_model='hamaker',
                         use_local_mwa : bool = True):
    """SIGH YOU CAN'T PICKLE AN EVERYBEAM TELESCOPE OBJECT"""
    
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
                               use_differential_beam : bool = True,
                               apply_beam_norms : bool = True,
                               reorder_jones : bool = False,
                               element_only : bool = False,
                               eb_rotate : bool = False,
                               parallactic_rotate : bool = True,
                               use_local_mwa : bool = True,
                               para_angle_offset : float = 0,
                               element_response_model='hamaker'):
    
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