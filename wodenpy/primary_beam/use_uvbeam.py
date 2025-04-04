
from pyuvdata import ShortDipoleBeam, BeamInterface, UVBeam
from pyuvdata.analytic_beam import AnalyticBeam
import numpy as np
import erfa
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
import concurrent.futures
from typing import Tuple
from copy import deepcopy
import os
import shutil
import time
from wodenpy.use_libwoden.skymodel_structs import Components_Python
import sys

def setup_MWA_uvbeams(hdf5_path: str, freqs: np.ndarray,
                      delays: np.ndarray = None,
                      amplitudes: np.ndarray = None,
                      pixels_per_deg : float = 5) -> np.array:
    """
    Setup the MWA uvbeam objects for a given set of frequencies and delays.
    
    Parameters
    ------------
    freqs : np.ndarray
        Array of frequencies in Hz.
    delays : np.ndarray, optional
        Array of delays as integers as written in the metafits file.
    amplitudes : np.ndarray, optional
        Array of amplitudes. If None, amplitudes of 1.0 used everywhere.
        Should be in format output from wodenpy.run_setup.check_args a.k.a 16 amps for NS
        dipole for tile 1, 16 amps for EW dipole for tile 1, 16 amps for NS dipole for tile 2, etc.

    Returns
    --------
    np.array
        Array of MWA uvbeam objects.
    """
    if delays is None:
        use_delays = np.zeros((2, 16), dtype='int')
    else:
        check = delays.size / 16
        if check != 1:
            sys.exit(f"Delays array size {delays.size} is equal to 16 "
                     "WODEN only points all tiles/polarisations to the same delay ."
                     "Exiting now as can't proceed sensibly.")
        
        use_delays = np.empty((2, 16), dtype='int')
        
        use_delays[0, :] = delays
        use_delays[1, :] = delays
        
    if amplitudes is None:
        amplitudes = np.ones(32, dtype='float')
        
    else:
        check = amplitudes.size % 32
        if check != 0:
            sys.exit(f"Amplitudes array size {amplitudes.size} is not a multiple of 32. "
                     "Expected 32 delays per tile so exiting now.")
    
    ##need to go a certain distance outside the requested range
    ##as UVBeam does interpolation over frequencies
    freq_range = [freqs.min()-2*1.28e6, freqs.max()+2*1.28e+6]
    
    num_beams = int(amplitudes.size / 32)
    
    uvbeam_objs = np.empty(num_beams, dtype=object)
    for beam_ind in range(num_beams):
        
        use_amplitudes = np.empty((2, 16), dtype='float')
        
        pol1 = amplitudes[beam_ind*32:beam_ind*32+16]
        pol2 = amplitudes[beam_ind*32+16:(beam_ind+1)*32]
        
        use_amplitudes[0, :] = pol1
        use_amplitudes[1, :] = pol2
        
        # print(use_amplitudes)
    
        mwabeam = UVBeam.from_file(hdf5_path,
                                   pixels_per_deg=pixels_per_deg,
                                   delays=use_delays,
                                   amplitudes=use_amplitudes,
                                   freq_range=freq_range)
        
        mwabeam.peak_normalize()
        
        uvbeam_objs[beam_ind] = mwabeam
    
    return uvbeam_objs

def run_uvbeam(uvbeam_objs: np.array,
               ras: np.ndarray, decs: np.ndarray,
               j2000_latitudes: np.ndarray, j2000_lsts: np.ndarray,
               freqs: np.ndarray,
               iau_order: bool = False,
               parallactic_rotate: bool = True) -> np.ndarray:
    """
    Calculate the Jones matrices for a given set of coordinates, times,
    frequencies, and station ids using the pyuvdata module.
    `j2000_latitudes` should be the array latitude as precessed back to J2000,
    with `j2000_lsts` being the matching LST in J2000. `current_latitude` and
    `current_longitude` should be latitude and longitude of the array at the
    time of the observation. 
    
    
    Parameters
    ------------
    uvbeam_objs: np.array
        Array of initiliased UVBeam objects.
    ras : np.ndarray
        Right ascensions of the coordinates in radians.
    decs : np.ndarray
         Declinations of the coordinates in radians.
    j2000_latitudes : np.ndarray
        Latitudes in J2000 coordinates.
    j2000_lsts : np.ndarray
        Local sidereal times in J2000 coordinates.
    times : np.ndarray
        Array of observation times.
    freqs : np.ndarray
        Array of frequencies.
    parallactic_rotate : bool, optional
        Whether to apply parallactic angle rotation using `wodenpy`. Defaults to False.
        Should be True for MWA beams if you want rotation. If True for a non-MWA beam,
        `wodenpy` should match the output as if `eb_rotate` was True.
    iau_order : bool, optional
        Whether to return the Jones matrices in IAU order. Defaults to False.
        
    Returns
    --------
    np.ndarray
        The calculated Jones matrices with shape (num_stations, num_times, num_freqs, num_coords, 2, 2).
    """
    
    num_stations = len(uvbeam_objs)
    num_times = len(j2000_lsts)
    num_freqs = len(freqs)
    num_coords = len(ras)
    
        # use the faster interpolation method if appropriate
    beam = BeamInterface(uvbeam_objs[0], beam_type="efield")
    if beam._isuvbeam and beam.beam.pixel_coordinate_system == "az_za":
        interpol_fn = "az_za_map_coordinates"
    else:
        interpol_fn = None
        
    check_azza_domain = False
    spline_opts = None
    
    all_output_jones = np.zeros((num_stations, num_times, num_freqs, num_coords, 2, 2), dtype=np.complex128)*np.nan
    
    if parallactic_rotate:
        coords = SkyCoord(ras*u.rad, decs*u.rad, frame='icrs')
    
    for time_ind in range(num_times):
        
        comp_has = j2000_lsts[time_ind] - ras
        azs, els = erfa.hd2ae(comp_has, decs, j2000_latitudes[time_ind])
        zas = np.pi/2 - els
        
        use_azs = np.pi/2 - azs
        where_neg_az = np.nonzero(use_azs < 0)
        use_azs[where_neg_az] = use_azs[where_neg_az] + np.pi * 2.
        
        # import matplotlib.pyplot as plt
        
        # fig, ax = plt.subplots(1, 2)
        
        # ax[0].imshow(np.reshape(azs, (51, 51)), origin='lower')
        # ax[1].imshow(np.reshape(use_azs, (51, 51)), origin='lower')
        
        # plt.show()
        
        # use_azs = azs
            
        if parallactic_rotate:
            has = j2000_lsts[time_ind] - ras
            para_angles = erfa.hd2pa(has, decs, j2000_latitudes[time_ind])
            
            # para_angles = np.pi/2 - para_angles
            # para_angles
            
            rot_matrix = np.empty((num_coords, 2,2))
            
            rot_matrix[:,0,0] = np.sin(-para_angles)
            rot_matrix[:,0,1] = -np.cos(-para_angles)
            rot_matrix[:,1,0] = -np.cos(-para_angles)
            rot_matrix[:,1,1] = -np.sin(-para_angles)
            
            # rot_matrix[:,0,0] = np.cos(-para_angles)
            # rot_matrix[:,0,1] = -np.sin(-para_angles)
            # rot_matrix[:,1,0] = np.sin(-para_angles)
            # rot_matrix[:,1,1] = np.cos(-para_angles)
            
        for station_ind, beam_obj in enumerate(uvbeam_objs):
            # for freq_ind, freq in enumerate(freqs):
            beam = BeamInterface(beam_obj, beam_type="efield")
            
            
            
            response = beam.compute_response(az_array=use_azs,
                                             za_array=zas,
                                             freq_array=freqs,
                                             check_azza_domain=check_azza_domain
                                             )
            
            # freq_interp_kind = "nearest"
            # response = beam.compute_response(az_array=use_azs,
            #                                  za_array=zas,
            #                                  freq_array=freqs,
            #                                  interpolation_function=interpol_fn,
            #                                  freq_interp_kind=freq_interp_kind,
            #                                  spline_opts=spline_opts,
            #                                  check_azza_domain=check_azza_domain
            #                                  )
            
            response[:,:,:,np.isnan(zas)] = np.nan
            
            all_output_jones[station_ind, time_ind, :, :, 0, 0] = response[1,0,:,:]
            all_output_jones[station_ind, time_ind, :, :, 0, 1] = -response[0,0,:,:]
            all_output_jones[station_ind, time_ind, :, :, 1, 0] = response[1,1,:,:]
            all_output_jones[station_ind, time_ind, :, :, 1, 1] = -response[0,1,:,:]
            
        if parallactic_rotate:
            ##Parallactic angle doesn't change per station or freq, only by
            ##time and direction
            rot_jones = np.einsum('ijklm,kmn->ijkln', all_output_jones[:, time_ind, :, :, :, :], rot_matrix)
            all_output_jones[:, time_ind, :, :, :, :] = rot_jones
                    
                    
    if iau_order:
        ##swap all_output_jones[:,:,:,:,0,0] with all_output_jones[:,:,:,:,1,1]
        all_output_jones[:, :, :, :, [0, 1], [0, 1]] = all_output_jones[:, :, :, :, [1, 0], [1, 0]]
        ##swap all_output_jones[:,:,:,:,0,1] with all_output_jones[:,:,:,:,1,0]
        all_output_jones[:, :, :, :, [0, 1], [1, 0]] = all_output_jones[:, :, :, :, [1, 0], [0, 1]]
                    
    return all_output_jones



def calc_uvbeam_for_components(components : Components_Python, uvbeam_objs : np.array,
                               all_freqs : np.ndarray,
                               j2000_latitudes : np.ndarray, j2000_lsts : np.ndarray,
                               iau_order : bool = True,
                               parallactic_rotate : bool = True):
    """
    Given a set of components, calculate the Jones matrices for each component
    at each time and frequency.
    
    
    Parameters
    ----------
   
    """
    
    all_jones = run_uvbeam(uvbeam_objs, components.ras, components.decs,
                           j2000_latitudes, j2000_lsts,
                           all_freqs,
                           iau_order=iau_order,
                           parallactic_rotate=parallactic_rotate)
    
    ##all_jones comes out as shape (num_beams, num_times, num_freqs, num_components, 2, 2)
    ##WODEN wants this to be flat. I've intentionally made the shape so we can just
    ##use the flatten call
    
    gxs = all_jones[:,:,:,:, 0,0]
    Dxs = all_jones[:,:,:,:, 0,1]
    Dys = all_jones[:,:,:,:, 1,0]
    gys = all_jones[:,:,:,:, 1,1]
    
    components.gxs = gxs.flatten()
    components.Dxs = Dxs.flatten()
    components.Dys = Dys.flatten()
    components.gys = gys.flatten()
    
    
    

#DOESN'T WORK as something inside the UVBeam object is triggering
#the GIL (I think)
def run_uvbeam_thread(num_threads : int, thread_id : int,
                      uvbeam_objs: np.array,
                      ras: np.ndarray, decs: np.ndarray,
                      j2000_latitudes: np.ndarray, j2000_lsts: np.ndarray,
                      times: np.ndarray, freqs: np.ndarray,
                      iau_order: bool = False,
                      parallactic_rotate: bool = True) -> Tuple[np.ndarray, int]:
    """
    Thread function called by `run_uvbeam_over_threads` to calculate the
    EveryBeam response in parrallel. Calls `run_uvbeam` with a subset of
    the coordinates; see `run_uvbeam` for more details of the parameters.
    
    Creates a new EveryBeam telescope object from `ms_path` for each thread.
    This has to be done because `concurrent.futures.ProcessPoolExecutor` has
    to pickle the function and all it's arguments, and EveryBeam objects can't
    be pickled. This is somewhat wasteful but I can't work out a better way
    to make things parallel.
    
    Parameters
    ------------
    num_threads : int
        Number of threads being in call by `run_uvbeam_over_threads`.
    thread_id : int
        ID of the current thread. Useds to work out what chunk of `ras` and `decs`
        to process.
    
        
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
    
    print(f"Thread {thread_id} processing coords {low_coord} to {high_coord}")
    
    # these_objs = [deepcopy(uvbeam_obj) for uvbeam_obj in uvbeam_objs]
    
    # shutil.copyfile(os.environ['MWA_FEE_HDF5'], f'mwa_fee_{thread_id}.hdf5')
    
    # delays1 = np.zeros((2, 16), dtype='int')
    # start = time.time()
    # mwabeam1 = UVBeam.from_file(f'mwa_fee_{thread_id}.hdf5', pixels_per_deg=5, delays=delays1,
    #                             freq_range=[freqs.min()-2*1.28e6, freqs.max()+2*1.28e+6])
    # print("read parallel", time.time() - start)
    
    these_objs = uvbeam_objs
    
    # start = time.time()
    jones = run_uvbeam(these_objs, ras[low_coord:high_coord],
                              decs[low_coord:high_coord],
                              j2000_latitudes, j2000_lsts,
                              times, freqs,
                              iau_order=iau_order,
                              parallactic_rotate=parallactic_rotate)
    # print("run parallel", time.time() - start)
    
    print(f"Thread {thread_id} finished")
    
    return jones, thread_id

def run_uvbeam_over_threads(num_threads : int,
                            uvbeam_objs: np.array,
                            ras: np.ndarray, decs: np.ndarray,
                            j2000_latitudes: np.ndarray, j2000_lsts: np.ndarray,
                            times: np.ndarray, freqs: np.ndarray,
                            iau_order: bool = False,
                            parallactic_rotate: bool = True):
    """
    Runs `run_uvbeam` in parallel over `num_threads` threads, using
    `concurrent.futures.ProcessPoolExecutor`. See `run_uvbeam` for more
    details of what each parameter does.
    
    Creates a new EveryBeam telescope object from `ms_path` for each thread.
    This has to be done because `concurrent.futures.ProcessPoolExecutor` has
    to pickle the function and all it's arguments, and EveryBeam objects can't
    be pickled. This is somewhat wasteful but I can't work out a better way
    to make things parallel.
    
    Parameters
    ------------
    num_threads : int
        Number of threads being in call by `run_uvbeam_over_threads`.
    ms_path : str
        Path to the measurement set to load the EveryBeam telescope from.
    ras : np.ndarray
        Right ascensions of the coordinates in radians.
    decs : np.ndarray
         Declinations of the coordinates in radians.
    j2000_latitudes : np.ndarray
        Latitudes in J2000 coordinates.
    j2000_lsts : np.ndarray
        Local sidereal times in J2000 coordinates.
    times : np.ndarray
        Array of observation times.
    freqs : np.ndarray
        Array of frequencies.
    iau_order : bool, optional
        If True, use IAU polarisation ordering, so set jones[0,0] to the NS dipole and jones[1,1] to EW. If False, jones[0,0] is EW.
    parallactic_rotate : bool, optional
        Whether to apply parallactic angle. Defaults to False.
        
    Returns
    --------
    np.ndarray
        The calculated Jones matrices with shape
        (num_stations, num_times, num_freqs, num_coord, 2, 2)
    """
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
        
            future_data = [executor.submit(run_uvbeam_thread,
                                           num_threads, thread_id,
                                           uvbeam_objs, ras, decs,
                                           j2000_latitudes, j2000_lsts,
                                           times, freqs,
                                           iau_order=iau_order,
                                           parallactic_rotate=parallactic_rotate)
                                    for thread_id in range(num_threads)]
            
            all_jones_chunks = []
            all_thread_ids = []
            for future in concurrent.futures.as_completed(future_data):
                jones_chunk, thread_id = future.result()
                all_jones_chunks.append(jones_chunk)
                all_thread_ids.append(thread_id)
                
    num_stations = len(uvbeam_objs)
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