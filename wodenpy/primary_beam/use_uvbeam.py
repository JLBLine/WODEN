"""Wrappers around pyuvdata.UVBeam to calculate primary beam values.

See https://pyuvdata.readthedocs.io/en/latest/uvbeam.html for more information
on the UVBeam class.
"""
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
    
    Delays should be as appears in the metafits file, so 16 integers.
    
    Amplitudes should have a length which is a multiple of 32; 16 amplitudes
    for each dipole, for each tile. So if you want three unqiue tile beams,
    submit an amplitude array of length 96.
    
    Parameters
    ------------
    hdf5_path : str
        Path to the hdf5 file containing the MWA FEE beam model.
    freqs : np.ndarray
        Array of frequencies in Hz.
    delays : np.ndarray, optional
        Array of delays as integers as written in the metafits file. If not
        given, all delays set to 0.
    amplitudes : np.ndarray, optional
        Array of amplitudes. If None, amplitudes of 1.0 used everywhere.
        Should be in format output from wodenpy.run_setup.check_args a.k.a 16 amps for NS
        dipole for tile 1, 16 amps for EW dipole for tile 1, 16 amps for NS dipole for tile 2, etc.
        Should be a multiple of 32.
    pixels_per_deg : float, optional
        Number of pixels per degree for UVBeam to use. Defaults to 5.

    Returns
    --------
    np.array
        Array of MWA uvbeam objects, of length len(amplitudes)/32.
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
        
        mwabeam = UVBeam.from_file(hdf5_path,
                                   pixels_per_deg=pixels_per_deg,
                                   delays=use_delays,
                                   amplitudes=use_amplitudes,
                                   freq_range=freq_range)
        
        mwabeam.peak_normalize()
        
        uvbeam_objs[beam_ind] = mwabeam
    
    return uvbeam_objs

def run_uvbeam(uvbeam_objs: np.ndarray[UVBeam],
               ras: np.ndarray, decs: np.ndarray,
               j2000_latitudes: np.ndarray, j2000_lsts: np.ndarray,
               freqs: np.ndarray,
               iau_order: bool = False,
               parallactic_rotate: bool = True,
               freq_interp=True) -> np.ndarray:
    """
    Calculate the Jones matrices for a given set of directions, times (via `j2000_lsts`
    and `j2000_latitudes`), frequencies, and stations (a.k.a tiles for MWA),
    using an array of UVBeam objects.
    `j2000_latitudes` should be the array latitude as precessed back to J2000,
    with `j2000_lsts` being the matching LST in J2000. Uses these values to
    calculate the azimuth and zenith angles for each station at each time.
    
    
    Parameters
    ------------
    uvbeam_objs: np.ndarray[UVBeam]
        Array of initiliased UVBeam objects.
    ras : np.ndarray
        Right ascensions of the coordinates in radians.
    decs : np.ndarray
         Declinations of the coordinates in radians.
    j2000_latitudes : np.ndarray
        Latitudes in J2000 coordinates.
    j2000_lsts : np.ndarray
        Local sidereal times in J2000 coordinates.
    freqs : np.ndarray
        Array of frequencies.
    iau_order : bool, optional
        Whether to return the Jones matrices in IAU order (NS = X, EW = Y). Defaults to False.
    parallactic_rotate : bool, optional
        Whether to apply parallactic angle rotation. Defaults to True.
        
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
    
    all_output_jones = np.zeros((num_stations, num_times, num_freqs, num_coords, 2, 2), dtype=np.complex128)*np.nan
    
    if parallactic_rotate:
        coords = SkyCoord(ras*u.rad, decs*u.rad, frame='icrs')
    
    for time_ind in range(num_times):
        
        comp_has = j2000_lsts[time_ind] - ras
        azs, els = erfa.hd2ae(comp_has, decs, j2000_latitudes[time_ind])
        zas = np.pi/2 - els
        
        ##pyuvdata defines azimuth as 0 at East, increasing to 90 at North
        ##astropy, erfa, and wodenpy define azimuth as 0 at North, increasing to 90 at East
        ##So need to rotate azimuth by 90 degrees and invert rotation direction
        ##classic astronomy
        use_azs = np.pi/2 - azs
        where_neg_az = np.nonzero(use_azs < 0)
        use_azs[where_neg_az] = use_azs[where_neg_az] + np.pi * 2.
        
        # use_azs = azs
        
        if parallactic_rotate:
            has = j2000_lsts[time_ind] - ras
            para_angles = erfa.hd2pa(has, decs, j2000_latitudes[time_ind])
            
            rot_matrix = np.empty((num_coords, 2,2))
            
            rot_matrix[:,0,0] = np.sin(-para_angles)
            rot_matrix[:,0,1] = -np.cos(-para_angles)
            rot_matrix[:,1,0] = -np.cos(-para_angles)
            rot_matrix[:,1,1] = -np.sin(-para_angles)
            
        for station_ind, beam_obj in enumerate(uvbeam_objs):
            beam = BeamInterface(beam_obj, beam_type="efield")
            
            if freq_interp:
                response = beam.compute_response(az_array=use_azs,
                                                za_array=zas,
                                                freq_array=freqs,
                                                check_azza_domain=check_azza_domain
                                                )
            else:
                freq_interp_kind = "nearest"
                spline_opts = None
                response = beam.compute_response(az_array=use_azs,
                                                za_array=zas,
                                                freq_array=freqs,
                                                interpolation_function=interpol_fn,
                                                freq_interp_kind=freq_interp_kind,
                                                spline_opts=spline_opts,
                                                check_azza_domain=check_azza_domain
                                                )
            
            response[:,:,:,np.isnan(zas)] = np.nan
            
            ##Correct for MWA beam convention
            all_output_jones[station_ind, time_ind, :, :, 0, 0] = response[1,0,:,:]
            all_output_jones[station_ind, time_ind, :, :, 0, 1] = -response[0,0,:,:]
            all_output_jones[station_ind, time_ind, :, :, 1, 0] = response[1,1,:,:]
            all_output_jones[station_ind, time_ind, :, :, 1, 1] = -response[0,1,:,:]
            
            
            # all_output_jones[station_ind, time_ind, :, :, 0, 0] = response[0,1,:,:]
            # all_output_jones[station_ind, time_ind, :, :, 0, 1] = response[0,0,:,:]
            # all_output_jones[station_ind, time_ind, :, :, 1, 0] = response[1,1,:,:]
            # all_output_jones[station_ind, time_ind, :, :, 1, 1] = response[1,0,:,:]
            
            
        if parallactic_rotate:
            ##Parallactic angle doesn't change per station or freq, only by
            ##time and direction
            rot_jones = np.einsum('ijklm,kmn->ijkln', all_output_jones[:, time_ind, :, :, :, :], rot_matrix)
            all_output_jones[:, time_ind, :, :, :, :] = rot_jones
                    
                    
    ##different telescope beams seem to come out with different east-west/ north-south
    ##polarisation orders. So change whether we reorder based on telescope name??
    
    
    reorder = False
    if uvbeam_objs[0].telescope_name == "HERA":
        if not iau_order:
            reorder = True
              
    else:
        if iau_order:
            reorder = True      
                    
    if reorder:
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
    """For a set of `components` with filled `ras,decs` attributes,
    calculate the primary beam Jones matrices for all UVBeam objects
    in `uvbeam_objs` at the given frequencies `all_freqs`. Uses `j2000_latitudes`
    and `j2000_lsts` to calculate `az,za` for all time steps.
    Stores the results in the `components` `gxs,Dxs,Dys,gys`
    attributes. The beam values will eventually be used by the compiled code,
    either on the CPU or GPU.

    Parameters
    ----------
    components : Components_Python
        An initialised Components_Python object with the ras, decs attributes
        filled with the coordinates of the components. Should also have
        the gxs, Dxs, Dys, gys attributes setup as an array of np.float32 or
        np.float64 of length `len(uvbeam_objs)*len(all_freqs)*len(j2000_lsts)*len(components.ras)`.
        float32 or float64 depends on if running simulation in "float" or "double" mode.
    uvbeam_objs: np.ndarray[UVBeam]
        Array of initiliased UVBeam objects.
    all_freqs : np.ndarray
        Array of frequencies in Hz.
    j2000_latitudes : np.ndarray
        Latitudes in J2000 coordinates.
    j2000_lsts : np.ndarray
        Local sidereal times in J2000 coordinates.
    freqs : np.ndarray
        Array of frequencies.
    iau_order : bool, optional
        Whether to return the Jones matrices in IAU order (NS = X, EW = Y). Defaults to False.
        parallactic_rotate : bool, optional
    parallactic_rotate : bool, optional
        Whether to apply parallactic angle rotation. Defaults to True.
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