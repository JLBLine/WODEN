import numpy as np
from astropy.coordinates import ITRS, SkyCoord, AltAz, EarthLocation
from astropy.time import Time, TimeDelta
import astropy.units as u
import argparse
from wodenpy.use_libwoden.create_woden_struct_classes import Woden_Struct_Classes
import erfa

# import mwa_hyperbeam

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
ROTATE = False
DO_ELEMENT_ONLY = False

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


def load_LOFAR_telescope(ms_path : str, response_model = "hamaker") -> eb.LOFAR:
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
    # telescope = eb.load_telescope(ms_path,
    #                               use_differential_beam=True,
    #                               element_response_model=response_model)
    
    telescope = eb.load_telescope(ms_path,
                                  use_differential_beam=USE_DIFFERENTIAL_BEAM,
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
    a given phase centre (beam_ra0, dec), time, frequency, telescope, and station.

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
                                          rotate=ROTATE)
    
    # if type(telescope) == eb.LOFAR:
    
    #     norm_x = 1 / np.abs(response[0,1])
    #     norm_y = 1 / np.abs(response[1,0])
        
    # else:
    norm_x = 1 / np.abs(response[0,0])
    norm_y = 1 / np.abs(response[1,1])
    
    return norm_x, norm_y

# def run_everybeam(dir_itrf : np.ndarray, phase_itrf : np.ndarray, 
#                   time : Time, freq : float, telescope: eb.Telescope,
#                   station_id : int = 0, beam_norms : np.ndarray = np.ones(2),
#                   reorder_jones : bool = False,
#                   ra : float = False, dec : float = False,
#                   parallactic_angle : float = 0,
#                   para_angle_offset : float = -np.pi/2) -> np.ndarray:
#     """For a given direction, phase centre, time, frequency, telescope, station,
#     calculate an everybeam jones matrix. Optionally normalise using the given
#     beam norms [X_norm, Y_norm]. Explicitly, where
#     `response` = [[j00, j01], [j10, j11]], norms are applied as
    
#     response[0,:] *= beam_norms[0]
#     response[1,:] *= beam_norms[1]
    
#     By defauly, reorder the Jones matrix to be
#     [-j11, j10, -j01, j00], as WODEN expects X = NS, Y = EW. The negatives
#     may be down to a definition of Stokes V. This reordering is the same as
#     done in `mwa_hyperbeam`. The reordering is done *after* the beam norms. As
#     long as the beam norms have been calculated using `get_everybeam_norm` the
#     ordering should be correct.

#     Parameters
#     ----------
#     dir_itrf : np.ndarray
#         XYZ itfr array (as output by `radec_to_xyz`) of direction of interest 
#     phase_itrf : np.ndarray
#         XYZ itfr array (as output by `radec_to_xyz`) of beam phase centre
#     time : Time
#         Astropy Time object of observation
#     freq : float
#         Frequency of observation in Hz
#     telescope : eb.Telescope
#         An everybeam telescope object
#     station_id : int, optional
#         Integer index of station, by default 0
#     beam_norms : np.ndarray, optional
#         Normalisation to apply, [X_norm, Y_norm], by default np.ones(2). Outputs are multiplied by these values
#     reorder_jones : bool, optional
#         If True, reorder the Jones matrix to be [j11, j10, j01, j00], which
#         is the order expected by WODEN, by default True
#     ra : float, optional
#         Right ascension of direction of interest (radians), by default False. Needed for MWA beam
#     dec : float, optional
#         Declination of direction of interest (radians), by default False. Needed for MWA beam
#     parallactic_angle : float, optional
#         Rotate by this parallactic angle (radians); by default 0
#     para_angle_offset : float, optional
#         Offset to apply to parallactic angle (radians); by default -np.pi/2

#     Returns
#     -------
#     np.ndarray
#         2x2 array of complex beam jones matrix [[j00, j01], [j10, j11]]
#     """
    
#     ##The MWA everybeam function takes ra, dec, not dir_itrf
    
#     if type(telescope) == eb.MWA:
#         ##Get the response
#         response = telescope.station_response(time.mjd*3600*24, station_id, freq,
#                                             ra, dec)
#     else:
        
#         # parallactic_angle = 0
#         if DO_ELEMENT_ONLY:
#             response = telescope.element_response(time.mjd*3600*24, station_id, 0, freq,
#                                                   dir_itrf,
#                                                   rotate=ROTATE)
#         else:
#             ##Get the response
#             response = telescope.station_response(time.mjd*3600*24, station_id, freq,
#                                               dir_itrf, phase_itrf, 
#                                               rotate=ROTATE)
        
#         # response = telescope.element_response(time.mjd*3600*24, station_id, 0, freq,
#         #                                       dir_itrf,
#         #                                       rotate=ROTATE)
        
#         # print(dir_itrf, phase_itrf)
        
#     ##normalise the beams using previously calculated norms
    
#     # print(beam_norms[0])
    
    
    
#     if parallactic_angle:
        
#         cospa = np.cos(parallactic_angle + para_angle_offset)
#         sinpa = np.sin(parallactic_angle + para_angle_offset)
        
#         rotated_response = np.zeros_like(response)
        
#         rotated_response[0,0] = response[0,0]*cospa - response[0,1]*sinpa
#         rotated_response[0,1] = response[0,0]*sinpa + response[0,1]*cospa
#         rotated_response[1,0] = response[1,0]*cospa - response[1,1]*sinpa
#         rotated_response[1,1] = response[1,0]*sinpa + response[1,1]*cospa
        
#         response = rotated_response
        
        
#     # response[0,:] *= beam_norms[0]
#     # response[1,:] *= beam_norms[1]
    
#     # if reorder_jones:
#     #     response = np.array([[response[1,1], -response[1,0]],
#     #                          [response[0,1], -response[0,0]]])
    
#     # if type(telescope) == eb.LOFAR:
#     #         response = np.array([[response[0,1], -response[0,0]],
#     #                              [response[1,1], -response[1,0]]])
        
#     ##diff stokes convention??
#     # if type(telescope) == eb.LOFAR:
#     #     response /= np.sqrt(2)
    
#     return response


def eb_local_xyz_from_radec(ra, dec, altaz_frame, delta_az=(1/2)*np.pi, 
                            negative_azimuth=True):
    """Get the local cartesian coords used by EveryBeam from a given RA, Dec,
    and AltAz frame. Reversing the azimuth and adding 90 degrees was found
    to match via trial and error"""
    
    coord = SkyCoord(ra=ra*u.rad, dec=dec*u.rad, frame='icrs')
    coord = coord.transform_to(altaz_frame)
    
    delta_az = delta_az * u.rad
    
    # print(coord.az, coord.az + delta_az)
    
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
        
    # print("BEFORE, AFTER", coord.az, updated_coord.az)
    
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


def run_everybeam(ras : np.ndarray, decs : np.ndarray,
                  beam_ra0 : float, beam_dec0 : float,
                  latitudes : np.ndarray, longitude : float,
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
    
    for time_ind, time in enumerate(times):
        if type(telescope) != eb.MWA:
            phase_itrf = radec_to_xyz(beam_ra0, beam_dec0, time)
        dir_itrfs = radec_to_xyz(ras, decs, time)
        
        time_mjd_secs = time.mjd*3600*24
        
        # ncp_t = radec_to_xyz(0, np.radians(90), time)
        
        coords = SkyCoord(ras*u.rad, decs*u.rad, frame='icrs')
        location = EarthLocation(lat=latitudes[time_ind]*u.rad,
                                 lon=longitude*u.rad)
        altaz_frame = AltAz(obstime=time, location=location)
        
        azza_coords = coords.transform_to(altaz_frame)
        azs = azza_coords.az.rad
        zas = np.pi/2 - azza_coords.alt.rad
        
        if type(telescope) == eb.MWA:
            # ncp_t = eb_local_xyz_from_radec(0, np.radians(90), altaz_frame,
            #                                 delta_az=0, negative_azimuth=False)
            # dir_local = eb_local_xyz_from_radec(ras, decs, altaz_frame,
            #                                 delta_az=0, negative_azimuth=False)
            
            ncp_t = eb_local_xyz_from_radec(0, np.radians(90), altaz_frame,
                                            negative_azimuth=False)
            dir_local = eb_local_xyz_from_radec(ras, decs, altaz_frame,
                                            negative_azimuth=False)
            
            beam = mwa_hyperbeam.FEEBeam()
            
        else:
            ncp_t = eb_local_xyz_from_radec(0, np.radians(90), altaz_frame)
            dir_local = eb_local_xyz_from_radec(ras, decs, altaz_frame)
        
        
        if parallactic_rotate:
            LST = time.sidereal_time('mean')
            LST_radians = np.radians(LST.value*15)
            has = LST_radians - ras
            para_angles = erfa.hd2pa(has, decs, latitudes[time_ind])
            # lat = latitudes[time_ind]
            
            rot_matrix = np.empty((num_coords, 2,2))
            
            if type(telescope) == eb.MWA:
            
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
                ##TODO get normalisation here
                beam_norms = np.ones(2)
                
                if apply_beam_norms:
                
                    if type(telescope) == eb.MWA:
                        ##Get the response
                        norm_jones = telescope.station_response(time_mjd_secs, station_id, freq,
                                                                beam_ra0, beam_dec0)
                        
                        # print("before rotate", norm_jones)
                        
                        # if parallactic_rotate:
                        #     LST = time.sidereal_time('mean')
                        #     LST_radians = np.radians(LST.value*15)
                        #     ha0 = LST_radians - beam_ra0
                        #     para_angles = erfa.hd2pa(ha0, beam_dec0, latitudes[time_ind])
                        #     # lat = latitudes[time_ind]
                        #     rot = np.empty((2,2))
                        #     rot[0,0] = np.sin(-para_angles)
                        #     rot[0,1] = -np.cos(-para_angles)
                        #     rot[1,0] = -np.cos(-para_angles)
                        #     rot[1,1] = -np.sin(-para_angles)
                        
                    else:
                        element_id = 0
                        if element_only:
                            norm_jones = telescope.element_response(time_mjd_secs, station_id, element_id, freq,
                                                                phase_itrf, rotate=eb_rotate)
                        else:
                            ##Get the response
                            # print("INPUT DIR, PHASE", dir_itrfs[coord_ind], phase_itrf)
                            norm_jones = telescope.station_response(time_mjd_secs, station_id, freq,
                                                            phase_itrf, phase_itrf, 
                                                            rotate=eb_rotate)
                
                    # if parallactic_rotate:
                    #     dir_phase_local = eb_local_xyz_from_radec(beam_ra0, beam_dec0, altaz_frame)
                    #     north, east = eb_north_east(dir_phase_local, ncp_t)
                    #     rot = calc_everybeam_rotation(dir_phase_local, north, east)
                    
                if parallactic_rotate and apply_beam_norms:
                    dir_phase_local = eb_local_xyz_from_radec(beam_ra0, beam_dec0, altaz_frame)
                    north, east = eb_north_east(dir_phase_local, ncp_t)
                    rot = calc_everybeam_rotation(dir_phase_local, north, east)
                    norm_jones = np.matmul(norm_jones, rot)
                    
                    # print("after rotate", norm_jones)
                    
                    
                # beam_norms = np.array([1/np.abs(norm_jones[0,0]), 1/np.abs(norm_jones[1,1])])
                
                for coord_ind, (ra, dec) in enumerate(zip(ras, decs)):
                    # if np.isnan(ra) or np.isnan(dec):
                    #     pass
                    # else:
                    
                    if type(telescope) == eb.MWA:
                            ##Get the response
                            response = telescope.station_response(time_mjd_secs, station_id, freq,
                                                                  ra, dec)
                            
                            # delays = [0,2,4,6,0,2,4,6,0,2,4,6,0,2,4,6]
                            # # response = beam.calc_jones(azs[coord_ind], zas[coord_ind], freq, delays, [1]*16, True,
                            # #                            np.radians(-26.7), True)
                            # response = beam.calc_jones(azs[coord_ind], zas[coord_ind], freq, delays, [1]*16, True)
                            # response.shape = (2,2)
                            
                            
                    else:
                        if element_only:
                            response = telescope.element_response(time_mjd_secs, station_id, 0, freq,
                                                                dir_itrfs[coord_ind], rotate=eb_rotate)
                        else:
                            ##Get the response
                            
                            # print("INPUT DIR, PHASE", dir_itrfs[coord_ind], phase_itrf)
                            response = telescope.station_response(time_mjd_secs, station_id, freq,
                                                            dir_itrfs[coord_ind], phase_itrf, 
                                                            rotate=eb_rotate)
                            
                    all_output_jones[station_ind, time_ind, freq_ind, coord_ind] = response
                
                if parallactic_rotate:
                    if type(telescope) == eb.MWA:
                        
                        # cospa = np.cos(para_angles + para_angle_offset)
                        # sinpa = np.sin(para_angles + para_angle_offset)
                        cospa = np.cos(para_angles)
                        sinpa = np.sin(para_angles)
                        
                        response = all_output_jones[:, time_ind, freq_ind, :, :, :]
                        
                        rotated_response = np.zeros_like(response)
                        
                        # rotated_response[:, :, 0,0] = response[:, :, 0,0]*cospa - response[:, :, 0,1]*sinpa
                        # rotated_response[:, :, 0,1] = response[:, :, 0,0]*sinpa + response[:, :, 0,1]*cospa
                        # rotated_response[:, :, 1,0] = response[:, :, 1,0]*cospa - response[:, :, 1,1]*sinpa
                        # rotated_response[:, :, 1,1] = response[:, :, 1,0]*sinpa + response[:, :, 1,1]*cospa
                        
                        rotated_response[:, :, 0,0] = -response[:, :, 1,0]*cospa + response[:, :, 1,1]*sinpa
                        rotated_response[:, :, 0,1] = -response[:, :, 1,0]*sinpa + -response[:, :, 1,1]*cospa
                        rotated_response[:, :, 1,0] = -response[:, :, 0,0]*cospa + response[:, :, 0,1]*sinpa
                        rotated_response[:, :, 1,1] = -response[:, :, 0,0]*sinpa + -response[:, :, 0,1]*cospa
                        
                        
                        # jones[2] * -c_rot + jones[3] * s_rot,
                        # jones[2] * -s_rot + jones[3] * -c_rot,
                        # jones[0] * -c_rot + jones[1] * s_rot,
                        # jones[0] * -s_rot + jones[1] * -c_rot,
                        
                        
                        all_output_jones[:, time_ind, freq_ind, :, :, :] = rotated_response
                        
                    else:  
                        rot_jones = np.einsum('jklm,kmn->jkln', all_output_jones[:, time_ind, freq_ind, :, :, :], rot_matrix)
                        all_output_jones[:, time_ind, freq_ind, :, :, :] = rot_jones
                    
                # if type(telescope) == eb.MWA:
                    
                #     reorder_jones = np.empty((num_stations, num_coords, 2, 2), dtype=np.complex128)
                    
                #     reorder_jones[:, :, 0, 0] = all_output_jones[:, time_ind, freq_ind, :, 1, 0]
                #     reorder_jones[:, :, 0, 1] = -all_output_jones[:, time_ind, freq_ind, :, 1, 1]
                #     reorder_jones[:, :, 1, 0] = all_output_jones[:, time_ind, freq_ind, :, 0, 0]
                #     reorder_jones[:, :, 1, 1] = -all_output_jones[:, time_ind, freq_ind, :, 0 ,1]
                    
                #     all_output_jones[:, time_ind, freq_ind, :, :, :] = reorder_jones
                    
                
                if apply_beam_norms:
                
                    ##Each station, time, and freq gets it's own normalisation
                    ##Same normalisation for all directions
                    
                    inv_beam_norms = np.linalg.inv(norm_jones)
                    output_jones = np.einsum('lm,kmn->kln', inv_beam_norms, all_output_jones[station_ind, time_ind, freq_ind, :, :, :])
                    all_output_jones[station_ind, time_ind, freq_ind, :, :, :] = output_jones
                    
    return all_output_jones