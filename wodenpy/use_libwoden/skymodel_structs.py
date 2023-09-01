import ctypes 
import importlib_resources
import numpy as np
from typing import Union
from ctypes import POINTER, c_float, c_double, c_int, c_char
import os
import sys
import erfa

from wodenpy.skymodel.woden_skymodel import Component_Type_Counter, CompTypes, Component_Info
from wodenpy.skymodel.chunk_sky_model import Skymodel_Chunk_Map
from wodenpy.use_libwoden.beam_settings import BeamTypes


Skymodel_Chunk_Map

VELC = 299792458.0
D2R = np.pi / 180.0

class Components_Float(ctypes.Structure):
    """A class structured equivalently to a `components_t` struct, used by 
    the C and CUDA code in libwoden_float.so
    
    :ivar POINTER(double) ras: COMPONENT right ascensions (radians)
    :cvar POINTER(double) decs: COMPONENT declinations (radians)
    :cvar POINTER(double) power_ref_freqs: COMPONENT Flux density reference frequencies (Hz)
    :cvar POINTER(user_precision_t) power_ref_stokesI: COMPONENT Stokes I reference flux density (Jy)
    :cvar POINTER(user_precision_t) power_ref_stokesQ: COMPONENT Stokes Q reference flux density (Jy)
    :cvar POINTER(user_precision_t) power_ref_stokesU: COMPONENT Stokes U reference flux density (Jy)
    :cvar POINTER(user_precision_t) power_ref_stokesV: COMPONENT Stokes V reference flux density (Jy)
    :cvar POINTER(user_precision_t) power_SIs:  COMPONENT spectral indexes
    :cvar POINTER(double) curve_ref_freqs: COMPONENT Flux density reference frequencies (Hz)
    :cvar POINTER(user_precision_t) curve_ref_stokesI: COMPONENT Stokes I reference flux density (Jy)
    :cvar POINTER(user_precision_t) curve_ref_stokesQ: COMPONENT Stokes Q reference flux density (Jy)
    :cvar POINTER(user_precision_t) curve_ref_stokesU: COMPONENT Stokes U reference flux density (Jy)
    :cvar POINTER(user_precision_t) curve_ref_stokesV: COMPONENT Stokes V reference flux density (Jy)
    :cvar POINTER(user_precision_t) curve_SIs:  COMPONENT spectral indexes
    :cvar POINTER(user_precision_t) curve_qs:  COMPONENT curvature
    :cvar POINTER(int) power_comp_inds: The indexes of all power-law models w.r.t ra,dec
    :cvar POINTER(int) curve_comp_inds: The indexes of all curved power-law models w.r.t ra,dec
    :cvar POINTER(int) list_comp_inds: The indexes of all list models w.r.t ra,dec
    :cvar POINTER(double) list_freqs: COMPONENT Flux density references frequencies (Hz)
    :cvar POINTER(user_precision_t) list_stokesI: COMPONENT Stokes I list flux density (Jy)
    :cvar POINTER(user_precision_t) list_stokesQ: COMPONENT Stokes Q list flux density (Jy)
    :cvar POINTER(user_precision_t) list_stokesU: COMPONENT Stokes U list flux density (Jy)
    :cvar POINTER(user_precision_t) list_stokesV: COMPONENT Stokes V list flux density (Jy)
    :cvar POINTER(int) num_list_values: How many freq/flux values are in each COMPONENT lis
    :cvar POINTER(int) list_start_indexes: How many freq/flux values are in each COMPONENT lis
    :cvar otal_num_flux_entires: The total number of freq/flux values are in all lists POINTER(int) combine
    :cvar POINTER(user_precision_t) extrap_stokesI: extrapolated COMPONENT Stokes I flux densities (Jy)
    :cvar POINTER(user_precision_t) extrap_stokesQ: extrapolated COMPONENT Stokes Q flux densities (Jy)
    :cvar POINTER(user_precision_t) extrap_stokesU: extrapolated COMPONENT Stokes U flux densities (Jy)
    :cvar POINTER(user_precision_t) extrap_stokesV: extrapolated COMPONENT Stokes V flux densities (Jy)
    :cvar POINTER(user_precision_t) shape_coeffs: Scaling coefficients for SHAPELET basis functions
    :cvar POINTER(user_precision_t) n1s: 1st basis function order for SHAPELET basis functions
    :cvar POINTER(user_precision_t) n2s: 2nd basis function order for SHAPELET basis functions
    :cvar POINTER(user_precision_t) majors: GAUSSIAN/SHAPELET major axis (beta1, radians)
    :cvar POINTER(user_precision_t) minors: GAUSSIAN/SHAPELET minor axis (beta2, radians)
    :cvar POINTER(user_precision_t) pas: GAUSSIAN/SHAPELET position angles (radians)
    :cvar POINTER(user_precision_t) param_indexes: An index value to match each coeff, n1, and n2 to the correct ra, dec, major, minor, pa for a SHAPELET
    :cvar POINTER(user_precision_t) azs: SHAPELET source azimuth angles for all time steps
    :cvar POINTER(user_precision_t) zas: SHAPELET source zenith angles for all time steps
    :cvar POINTER(double) beam_has: Hour angle of COMPONENTs for all time steps, used for beam calculations
    :cvar POINTER(double) beam_decs: Declinations of COMPONENTs for all time steps, used for beam calculations
    :cvar POINTER(int) num_primarybeam_values: Number of beam calculations needed for COMPONENTs
    :cvar POINTER(user_precision_complex_t) gxs: North-South Beam gain values for all directions, frequencies, and times for these COMPONENT
    :cvar POINTER(user_precision_complex_t) Dxs: North-South Beam leakage values for all directions, frequencies, and times for these COMPONENT
    :cvar POINTER(user_precision_complex_t) Dys: East-West Beam leakage values for all directions, frequencies, and times  for these COMPONENT
    :cvar POINTER(user_precision_complex_t) gys: East-West Beam gain values for all directions, frequencies, and times  for these COMPONENT
    :cvar POINTER(double) ls: Device memory l cosine direction coords for these COMPONENTs
    :cvar POINTER(double) ms: Device memory m cosine direction coords for these COMPONENTs
    :cvar POINTER(double) ns: Device memory n cosine direction coords for these COMPONENTs
    
    """
    
    _fields_ = [## Instrinsic to COMPONENT values
                ("ras", POINTER(c_double)),
                ("decs", POINTER(c_double)),
                ## power law params
                ("power_ref_freqs", POINTER(c_double)),
                ("power_ref_stokesI", POINTER(c_float)),
                ("power_ref_stokesQ", POINTER(c_float)),
                ("power_ref_stokesU", POINTER(c_float)),
                ("power_ref_stokesV", POINTER(c_float)),
                ("power_SIs", POINTER(c_float)),
                ##curved power law params
                ("curve_ref_freqs", POINTER(c_double)),
                ("curve_ref_stokesI", POINTER(c_float)),
                ("curve_ref_stokesQ", POINTER(c_float)),
                ("curve_ref_stokesU", POINTER(c_float)),
                ("curve_ref_stokesV", POINTER(c_float)),
                ("curve_SIs", POINTER(c_float)),
                ("curve_qs", POINTER(c_float)),
                ##indexes of types
                ("power_comp_inds", POINTER(c_int)),
                ("curve_comp_inds", POINTER(c_int)),
                ("list_comp_inds", POINTER(c_int)),
                ##list flux params
                ("list_freqs", POINTER(c_double)),
                ("list_stokesI", POINTER(c_float)),
                ("list_stokesQ", POINTER(c_float)),
                ("list_stokesU", POINTER(c_float)),
                ("list_stokesV", POINTER(c_float)),
                ("num_list_values", POINTER(c_int)),
                ("list_start_indexes", POINTER(c_int)),
                ("total_num_flux_entires", c_int),
                ##something to store extrapolated output fluxes in
                ("extrap_stokesI", POINTER(c_float)),
                ("extrap_stokesQ", POINTER(c_float)),
                ("extrap_stokesU", POINTER(c_float)),
                ("extrap_stokesV", POINTER(c_float)),
                ##SHAPELET params
                ("shape_coeffs", POINTER(c_float)),
                ("n1s", POINTER(c_float)),
                ("n2s", POINTER(c_float)),
                ##SHAPELET / GAUSSIAN params
                ("majors", POINTER(c_float)),
                ("minors", POINTER(c_float)),
                ("pas", POINTER(c_float)),
                ("param_indexes", POINTER(c_float)),
                ##Specific to observation settings for these COMPONENTs
                ("azs", POINTER(c_float)),
                ("zas", POINTER(c_float)),
                ("beam_has", POINTER(c_double)),
                ("beam_decs", POINTER(c_double)),
                ("num_primarybeam_values", c_int),
                ##things to hold the beam gain
                ##these are _complex in the C struct; we are not allocating
                ##memory on the python side, so I think we can just call them
                ##a double here?
                ("gxs", POINTER(c_float)),
                ("Dxs", POINTER(c_float)),
                ("Dys", POINTER(c_float)),
                ("gys", POINTER(c_float)),
                ##used to hold l,m,n coords on the GPU
                ("ls", POINTER(c_double)),
                ("ms", POINTER(c_double)),
                ("ns", POINTER(c_double)),
                ]
    
class Components_Double(ctypes.Structure):
    """A class structured equivalently to a `components_t` struct, used by 
    the C and CUDA code in libwoden_double.so
    
    :ivar POINTER(double) ras: COMPONENT right ascensions (radians)
    :cvar POINTER(double) decs: COMPONENT declinations (radians)
    :cvar POINTER(double) power_ref_freqs: COMPONENT Flux density reference frequencies (Hz)
    :cvar POINTER(user_precision_t) power_ref_stokesI: COMPONENT Stokes I reference flux density (Jy)
    :cvar POINTER(user_precision_t) power_ref_stokesQ: COMPONENT Stokes Q reference flux density (Jy)
    :cvar POINTER(user_precision_t) power_ref_stokesU: COMPONENT Stokes U reference flux density (Jy)
    :cvar POINTER(user_precision_t) power_ref_stokesV: COMPONENT Stokes V reference flux density (Jy)
    :cvar POINTER(user_precision_t) power_SIs:  COMPONENT spectral indexes
    :cvar POINTER(double) curve_ref_freqs: COMPONENT Flux density reference frequencies (Hz)
    :cvar POINTER(user_precision_t) curve_ref_stokesI: COMPONENT Stokes I reference flux density (Jy)
    :cvar POINTER(user_precision_t) curve_ref_stokesQ: COMPONENT Stokes Q reference flux density (Jy)
    :cvar POINTER(user_precision_t) curve_ref_stokesU: COMPONENT Stokes U reference flux density (Jy)
    :cvar POINTER(user_precision_t) curve_ref_stokesV: COMPONENT Stokes V reference flux density (Jy)
    :cvar POINTER(user_precision_t) curve_SIs:  COMPONENT spectral indexes
    :cvar POINTER(user_precision_t) curve_qs:  COMPONENT curvature
    :cvar POINTER(int) power_comp_inds: The indexes of all power-law models w.r.t ra,dec
    :cvar POINTER(int) curve_comp_inds: The indexes of all curved power-law models w.r.t ra,dec
    :cvar POINTER(int) list_comp_inds: The indexes of all list models w.r.t ra,dec
    :cvar POINTER(double) list_freqs: COMPONENT Flux density references frequencies (Hz)
    :cvar POINTER(user_precision_t) list_stokesI: COMPONENT Stokes I list flux density (Jy)
    :cvar POINTER(user_precision_t) list_stokesQ: COMPONENT Stokes Q list flux density (Jy)
    :cvar POINTER(user_precision_t) list_stokesU: COMPONENT Stokes U list flux density (Jy)
    :cvar POINTER(user_precision_t) list_stokesV: COMPONENT Stokes V list flux density (Jy)
    :cvar POINTER(int) num_list_values: How many freq/flux values are in each COMPONENT lis
    :cvar POINTER(int) list_start_indexes: How many freq/flux values are in each COMPONENT lis
    :cvar otal_num_flux_entires: The total number of freq/flux values are in all lists POINTER(int) combine
    :cvar POINTER(user_precision_t) extrap_stokesI: extrapolated COMPONENT Stokes I flux densities (Jy)
    :cvar POINTER(user_precision_t) extrap_stokesQ: extrapolated COMPONENT Stokes Q flux densities (Jy)
    :cvar POINTER(user_precision_t) extrap_stokesU: extrapolated COMPONENT Stokes U flux densities (Jy)
    :cvar POINTER(user_precision_t) extrap_stokesV: extrapolated COMPONENT Stokes V flux densities (Jy)
    :cvar POINTER(user_precision_t) shape_coeffs: Scaling coefficients for SHAPELET basis functions
    :cvar POINTER(user_precision_t) n1s: 1st basis function order for SHAPELET basis functions
    :cvar POINTER(user_precision_t) n2s: 2nd basis function order for SHAPELET basis functions
    :cvar POINTER(user_precision_t) majors: GAUSSIAN/SHAPELET major axis (beta1, radians)
    :cvar POINTER(user_precision_t) minors: GAUSSIAN/SHAPELET minor axis (beta2, radians)
    :cvar POINTER(user_precision_t) pas: GAUSSIAN/SHAPELET position angles (radians)
    :cvar POINTER(user_precision_t) param_indexes: An index value to match each coeff, n1, and n2 to the correct ra, dec, major, minor, pa for a SHAPELET
    :cvar POINTER(user_precision_t) azs: SHAPELET source azimuth angles for all time steps
    :cvar POINTER(user_precision_t) zas: SHAPELET source zenith angles for all time steps
    :cvar beam_has: Hour angle of COMPONENTs for all time steps, used for beam calculPOINTER(double) ations
    :cvar beam_decs: Declinations of COMPONENTs for all time steps, used for beam calculPOINTER(double) ations
    :cvar POINTER(int) num_primarybeam_values: Number of beam calculations needed for COMPONENTs
    :cvar POINTER(user_precision_complex_t) gxs: North-South Beam gain values for all directions, frequencies, and times for these COMPONENT
    :cvar POINTER(user_precision_complex_t) Dxs: North-South Beam leakage values for all directions, frequencies, and times for these COMPONENT
    :cvar POINTER(user_precision_complex_t) Dys: East-West Beam leakage values for all directions, frequencies, and times  for these COMPONENT
    :cvar POINTER(user_precision_complex_t) gys: East-West Beam gain values for all directions, frequencies, and times  for these COMPONENT
    :cvar POINTER(double) ls: Device memory l cosine direction coords for these COMPONENTs
    :cvar POINTER(double) ms: Device memory m cosine direction coords for these COMPONENTs
    :cvar POINTER(double) ns: Device memory n cosine direction coords for these COMPONENTs
    
    """
    
    _fields_ = [## Instrinsic to COMPONENT values
                ("ras", POINTER(c_double)),
                ("decs", POINTER(c_double)),
                ## power law params
                ("power_ref_freqs", POINTER(c_double)),
                ("power_ref_stokesI", POINTER(c_double)),
                ("power_ref_stokesQ", POINTER(c_double)),
                ("power_ref_stokesU", POINTER(c_double)),
                ("power_ref_stokesV", POINTER(c_double)),
                ("power_SIs", POINTER(c_double)),
                ##curved power law params
                ("curve_ref_freqs", POINTER(c_double)),
                ("curve_ref_stokesI", POINTER(c_double)),
                ("curve_ref_stokesQ", POINTER(c_double)),
                ("curve_ref_stokesU", POINTER(c_double)),
                ("curve_ref_stokesV", POINTER(c_double)),
                ("curve_SIs", POINTER(c_double)),
                ("curve_qs", POINTER(c_double)),
                ##indexes of types
                ("power_comp_inds", POINTER(c_int)),
                ("curve_comp_inds", POINTER(c_int)),
                ("list_comp_inds", POINTER(c_int)),
                ##list flux params
                ("list_freqs", POINTER(c_double)),
                ("list_stokesI", POINTER(c_double)),
                ("list_stokesQ", POINTER(c_double)),
                ("list_stokesU", POINTER(c_double)),
                ("list_stokesV", POINTER(c_double)),
                ("num_list_values", POINTER(c_int)),
                ("list_start_indexes", POINTER(c_int)),
                ("total_num_flux_entires", c_int),
                ##something to store extrapolated output fluxes in
                ("extrap_stokesI", POINTER(c_double)),
                ("extrap_stokesQ", POINTER(c_double)),
                ("extrap_stokesU", POINTER(c_double)),
                ("extrap_stokesV", POINTER(c_double)),
                ##SHAPELET params
                ("shape_coeffs", POINTER(c_double)),
                ("n1s", POINTER(c_double)),
                ("n2s", POINTER(c_double)),
                ##SHAPELET / GAUSSIAN params
                ("majors", POINTER(c_double)),
                ("minors", POINTER(c_double)),
                ("pas", POINTER(c_double)),
                ("param_indexes", POINTER(c_double)),
                ##Specific to observation settings for these COMPONENTs
                ("azs", POINTER(c_double)),
                ("zas", POINTER(c_double)),
                ("beam_has", POINTER(c_double)),
                ("beam_decs", POINTER(c_double)),
                ("num_primarybeam_values", c_int),
                ##things to hold the beam gain
                ##these are _complex in the C struct; we are not allocating
                ##memory on the python side, so I think we can just call them
                ##a double here?
                ("gxs", POINTER(c_double)),
                ("Dxs", POINTER(c_double)),
                ("Dys", POINTER(c_double)),
                ("gys", POINTER(c_double)),
                ##used to hold l,m,n coords on the GPU
                ("ls", POINTER(c_double)),
                ("ms", POINTER(c_double)),
                ("ns", POINTER(c_double)),
                ]
    
class Source_Float(ctypes.Structure):
    """A class structured equivalent to a `source_t` struct, used by 
    the C and CUDA code in libwoden_float.so
    """
    
    _fields_ = [
        ('name', (c_char*32)),
        ("n_comps", c_int),
        ("n_points", c_int),
        ("n_point_lists", c_int),
        ("n_point_powers", c_int),
        ("n_point_curves", c_int),
        ("n_gauss", c_int),
        ("n_gauss_lists", c_int),
        ("n_gauss_powers", c_int),
        ("n_gauss_curves", c_int),
        ("n_shapes", c_int),
        ("n_shape_lists", c_int),
        ("n_shape_powers", c_int),
        ("n_shape_curves", c_int),
        ("n_shape_coeffs", c_int),
        ("point_components", Components_Float),
        ("gauss_components", Components_Float),
        ("shape_components", Components_Float),
        ("d_point_components", Components_Float),
        ("d_gauss_components", Components_Float),
        ("d_shape_components", Components_Float),
    ]
    
class Source_Double(ctypes.Structure):
    """A class structured equivalent to a `source_t` struct, used by 
    the C and CUDA code in libwoden_double.so
    
    :cvar c_char*32 name: Source name
    :cvar c_int n_comps: Total number of COMPONENTs in source 
    :cvar c_int n_points: Number of POINT source COMPONENTs 
    :cvar c_int n_point_lists: Number of POINTs with LIST type flux
    :cvar c_int n_point_powers: Number of POINTs with POWER_LAW type flux
    :cvar c_int n_point_curves: Number of POINTs with CURVED_POWER_LAW type flux
    :cvar c_int n_gauss: Number of GAUSSIAN source COMPONENTs
    :cvar c_int n_gauss_lists: Number of GAUSSIANs with LIST type flux
    :cvar c_int n_gauss_powers: Number of GAUSSIANs with POWER_LAW type flux
    :cvar c_int n_gauss_curves: Number of GAUSSIANs with CURVED_POWER_LAW type flux
    :cvar c_int n_shapes: Number of SHAPELET source COMPONENTs
    :cvar c_int n_shape_lists: Number of SHAPELETs with LIST type flux
    :cvar c_int n_shape_powers: Number of SHAPELETs with POWER_LAW type flux
    :cvar c_int n_shape_curves: Number of SHAPELETs with CURVED_POWER_LAW type flux
    :cvar c_int n_shape_coeffs: Total number of SHAPELET coefficients
    :cvar Components_Double point_components: `Components_Double` holding component information for all POINT COMPONENTs in this SOURCE
    :cvar Components_Double gauss_components: `Components_Double` holding component information for all GAUSSIAN COMPONENTs in this SOURCE
    :cvar Components_Double shape_components: `Components_Double` holding component information for all SHAPELET COMPONENTs in this SOURCE
    :cvar Components_Double d_point_components: `Components_Double` holding component information on the device for all POINT COMPONENTs in this SOURCE
    :cvar Components_Double d_gauss_components: `Components_Double` holding component information on the device for all GAUSSIAN COMPONENTs in this SOURCE
    :cvar Components_Double d_shape_components: `Components_Double` holding component information on the device for all SHAPELET COMPONENTs in this SOURCE
    
    """
    
    _fields_ = [
        ('name', (c_char*32)),
        ("n_comps", c_int),
        ("n_points", c_int),
        ("n_point_lists", c_int),
        ("n_point_powers", c_int),
        ("n_point_curves", c_int),
        ("n_gauss", c_int),
        ("n_gauss_lists", c_int),
        ("n_gauss_powers", c_int),
        ("n_gauss_curves", c_int),
        ("n_shapes", c_int),
        ("n_shape_lists", c_int),
        ("n_shape_powers", c_int),
        ("n_shape_curves", c_int),
        ("n_shape_coeffs", c_int),
        ("point_components", Components_Double),
        ("gauss_components", Components_Double),
        ("shape_components", Components_Double),
        ("d_point_components", Components_Double),
        ("d_gauss_components", Components_Double),
        ("d_shape_components", Components_Double),
    ]
    
class Source_Catalogue_Float(ctypes.Structure):
    """A class structured equivalent to a `source_t` struct, used by 
    the C and CUDA code in libwoden_float.so
    """
    
    _fields_ = [("num_sources", c_int),
                ("num_shapelets", c_int),
                ("sources", POINTER(Source_Float))]
    
class Source_Catalogue_Double(ctypes.Structure):
    """A class structured equivalent to a `source_t` struct, used by 
    the C and CUDA code in libwoden_float.so
    """
    
    _fields_ = [("num_sources", c_int),
                ("num_shapelets", c_int),
                ("sources", POINTER(Source_Double))]
    
    
def setup_source_catalogue(num_sources : int, num_shapelets : int,
                           precision = "double"):
    
    if precision == 'float':
        source_catalogue = Source_Catalogue_Float()
        # source_catalogue = POINTER(Source_Catalogue_Float)
        source_array = num_sources*Source_Float
        
    else:
        source_catalogue = Source_Catalogue_Double()
        # source_catalogue = POINTER(Source_Catalogue_Double)
        source_array = num_sources*Source_Double
        
    source_catalogue.sources = source_array()
    source_catalogue.num_sources = num_sources
    source_catalogue.num_shapelets = num_shapelets
    
    return source_catalogue
    
def setup_components(chunk_map : Skymodel_Chunk_Map,
                     chunked_source : Union[Source_Float, Source_Double],
                     num_freqs : int,
                     num_times : int, comp_type : CompTypes,
                     beamtype : int, 
                     c_user_precision : Union[c_float, c_double]):
    """choose which components and do the ctypes "malloc" thing
    """
    
    if comp_type == CompTypes.POINT:
        
        n_comps = chunk_map.n_points
        components = chunked_source.point_components
        
        ##flux type specific things - need arrays of certain lengths, in
        ##in 'user' precision, int, and double
        power_user_ncomps_arr = c_user_precision*chunk_map.n_point_powers
        curve_user_ncomps_arr = c_user_precision*chunk_map.n_point_curves
        
        power_int_ncomps_arr = c_int*chunk_map.n_point_powers
        curve_int_ncomps_arr = c_int*chunk_map.n_point_curves
        
        power_double_ncomps_arr = c_double*chunk_map.n_point_powers
        curve_double_ncomps_arr = c_double*chunk_map.n_point_curves
        
        ##the number of entires for list flux types isn't the number
        ##of components, are any component can have multiple list entries
        list_user_nflux_arr = c_user_precision*chunk_map.point_components.total_num_flux_entires
        list_double_nflux_arr = c_double*chunk_map.point_components.total_num_flux_entires
        
        list_int_ncomps_arr = c_int*chunk_map.n_point_lists
        n_lists = chunk_map.point_components.total_num_flux_entires
        
    elif comp_type == CompTypes.GAUSSIAN:
        
        n_comps = chunk_map.n_gauss
        components = chunked_source.gauss_components
        
        ##flux type specific things
        power_user_ncomps_arr = c_user_precision*chunk_map.n_gauss_powers
        curve_user_ncomps_arr = c_user_precision*chunk_map.n_gauss_curves
        
        power_int_ncomps_arr = c_int*chunk_map.n_gauss_powers
        curve_int_ncomps_arr = c_int*chunk_map.n_gauss_curves
        
        power_double_ncomps_arr = c_double*chunk_map.n_gauss_powers
        curve_double_ncomps_arr = c_double*chunk_map.n_gauss_curves
        
        ##the number of entires for list flux types isn't the number
        ##of components, are any component can have multiple list entries
        list_user_nflux_arr = c_user_precision*chunk_map.gauss_components.total_num_flux_entires
        list_double_nflux_arr = c_double*chunk_map.gauss_components.total_num_flux_entires
        
        list_int_ncomps_arr = c_int*chunk_map.n_gauss_lists
        n_lists = chunk_map.gauss_components.total_num_flux_entires
        
    elif comp_type == CompTypes.SHAPELET:
        
        n_comps = chunk_map.n_shapes
        components = chunked_source.shape_components
        
        ##flux type specific things
        power_user_ncomps_arr = c_user_precision*chunk_map.n_shape_powers
        curve_user_ncomps_arr = c_user_precision*chunk_map.n_shape_curves
        
        power_int_ncomps_arr = c_int*chunk_map.n_shape_powers
        curve_int_ncomps_arr = c_int*chunk_map.n_shape_curves
        
        power_double_ncomps_arr = c_double*chunk_map.n_shape_powers
        curve_double_ncomps_arr = c_double*chunk_map.n_shape_curves
        
        ##the number of entires for list flux types isn't the number
        ##of components, are any component can have multiple list entries
        list_user_nflux_arr = c_user_precision*chunk_map.shape_components.total_num_flux_entires
        list_double_nflux_arr = c_double*chunk_map.shape_components.total_num_flux_entires
        
        list_int_ncomps_arr = c_int*chunk_map.n_shape_lists
        n_lists = chunk_map.shape_components.total_num_flux_entires
        
    ##component type specific things
    user_ncomps_arr = c_user_precision*n_comps
    double_ncomps_arr = c_double*n_comps
    num_primarybeam_values = n_comps*num_freqs*num_times
        
    components.ras = double_ncomps_arr()
    components.decs = double_ncomps_arr()
    components.num_primarybeam_values = num_primarybeam_values
    
    ##power-law flux things
    components.power_ref_freqs = power_double_ncomps_arr() ##this is always double
    components.power_ref_stokesI = power_user_ncomps_arr()
    components.power_ref_stokesQ = power_user_ncomps_arr()
    components.power_ref_stokesU = power_user_ncomps_arr()
    components.power_ref_stokesV = power_user_ncomps_arr()
    components.power_SIs = power_user_ncomps_arr()
    
    ##curved power-law flux things
    components.curve_ref_freqs = curve_double_ncomps_arr() ##this is always double
    components.curve_ref_stokesI = curve_user_ncomps_arr()
    components.curve_ref_stokesQ = curve_user_ncomps_arr()
    components.curve_ref_stokesU = curve_user_ncomps_arr()
    components.curve_ref_stokesV = curve_user_ncomps_arr()
    components.curve_SIs = curve_user_ncomps_arr()
    components.curve_qs = curve_user_ncomps_arr()
    
    ##list flux things
    components.list_freqs = list_double_nflux_arr() ##this is always double
    components.list_stokesI = list_user_nflux_arr()
    components.list_stokesQ = list_user_nflux_arr()
    components.list_stokesU = list_user_nflux_arr()
    components.list_stokesV = list_user_nflux_arr()
    
    components.num_list_values = list_int_ncomps_arr()
    components.list_start_indexes = list_int_ncomps_arr()
    components.total_num_flux_entires = n_lists
    
    # print("WE BE ASSIGNED THIS NUM", components.total_num_flux_entires)
    
    ##indexes of types
    components.power_comp_inds = power_int_ncomps_arr()
    components.curve_comp_inds = curve_int_ncomps_arr()
    components.list_comp_inds = list_int_ncomps_arr()
    
    if comp_type == CompTypes.GAUSSIAN or comp_type == CompTypes.SHAPELET:
        components.majors = user_ncomps_arr()
        components.minors = user_ncomps_arr()
        components.pas = user_ncomps_arr()
        
        
    if comp_type == CompTypes.SHAPELET:
        user_ncoeffs_arr = c_user_precision*chunk_map.shape_components.total_shape_coeffs
        components.shape_coeffs = user_ncoeffs_arr()
        components.n1s = user_ncoeffs_arr()
        components.n2s = user_ncoeffs_arr()
        components.param_indexes = user_ncoeffs_arr()
        
    ##----------------------------------------------------------------
    ##now we make space for coordinates that are need for the primary beam
    
    if beamtype == BeamTypes.GAUSS_BEAM.value or beamtype == BeamTypes.MWA_ANALY.value:
        hadec_arr = c_double*(n_comps*num_times)
        components.beam_has = hadec_arr()
        components.beam_decs = hadec_arr()
        
    ##only the NO_BEAM and GAUSS_BEAM options don't need az,za coords
    if beamtype == BeamTypes.GAUSS_BEAM.value or beamtype == BeamTypes.NO_BEAM.value:
        pass
    else:
        azza_arr = c_user_precision*(n_comps*num_times)
        components.azs = azza_arr()
        components.zas = azza_arr()
        
    
def setup_chunked_source(chunk_map : Skymodel_Chunk_Map, num_freqs : int,
                         num_times : int, beamtype : int,
                         precision='double') -> ctypes.Structure:
    """Sets up a ctypes structure class to contain a chunked sky model.
    This class is compatible with the C/CUDA code, and will allocate the
    correct amount of memory, based on whether the precision is either
    'double' or 'float'.
    
    """
    
    if precision == 'float':
        
        c_user_precision = c_float
        chunked_source = Source_Float()
        
    else:
        c_user_precision = c_double
        chunked_source = Source_Double()
        
    chunked_source.n_points = chunk_map.n_points
    chunked_source.n_point_lists = chunk_map.n_point_lists
    chunked_source.n_point_powers = chunk_map.n_point_powers
    chunked_source.n_point_curves = chunk_map.n_point_curves
    chunked_source.n_gauss = chunk_map.n_gauss
    chunked_source.n_gauss_lists = chunk_map.n_gauss_lists
    chunked_source.n_gauss_powers = chunk_map.n_gauss_powers
    chunked_source.n_gauss_curves = chunk_map.n_gauss_curves
    chunked_source.n_shapes = chunk_map.n_shapes
    chunked_source.n_shape_lists = chunk_map.n_shape_lists
    chunked_source.n_shape_powers = chunk_map.n_shape_powers
    chunked_source.n_shape_curves = chunk_map.n_shape_curves
    chunked_source.n_shape_coeffs = chunk_map.n_shape_coeffs
    chunked_source.n_comps = chunk_map.n_comps
        
    if chunk_map.n_points > 0:
        setup_components(chunk_map, chunked_source, num_freqs, num_times,
                         CompTypes.POINT, beamtype, c_user_precision)
        
    if chunk_map.n_gauss > 0:
        setup_components(chunk_map, chunked_source, num_freqs, num_times,
                         CompTypes.GAUSSIAN, beamtype, c_user_precision)
        
    if chunk_map.n_shapes > 0:
        setup_components(chunk_map, chunked_source, num_freqs, num_times,
                         CompTypes.SHAPELET, beamtype, c_user_precision)
    
    return chunked_source


def add_info_to_source_catalogue(chunked_skymodel_maps : list,
                                 source_catalogue : Union[Source_Catalogue_Float, Source_Catalogue_Double],
                                 orig_comp_ind : int, comp_info : Component_Info,
                                 map_comp_to_chunk : np.ndarray,
                                 all_chunk_comp_indexes : np.ndarray,
                                 beamtype : int, lsts : np.ndarray,
                                 latitude : float,
                                 collected_comps : int):
    
    ##Finalise all the info inside `comp_info`, and warn in 
    empty_fluxes = comp_info.finalise_comp()
                            
    if len(empty_fluxes) > 0:
        print(f"WARNING: {len(empty_fluxes)} components in source '{comp_info.source_name}' have no flux values. Setting to zero")
    
    
    chunk_indexes = np.unique(map_comp_to_chunk[np.where(all_chunk_comp_indexes == orig_comp_ind)]).astype(int)
    
    ##you can have the same component in multiple chunks for a shapelet
    ##as it's split over shapelet basis function, not component
    ##this iteration will only every be length one for point/gaussian
    for chunk_ind in chunk_indexes:
        
        chunk_map = chunked_skymodel_maps[chunk_ind]
        chunked_source = source_catalogue.sources[chunk_ind]
        
        ##TODO - gotta work out some way count the index of each component
        ##type to make sure we're adding the right things to
        # chunked_source.*_shapelet.param_indexes
        
        ##select the component type we need to populate
        if comp_info.point:
            
            source_components = chunked_source.point_components
            map_components = chunk_map.point_components
            
            n_powers = chunk_map.n_point_powers
            n_curves = chunk_map.n_point_curves
            
        elif comp_info.gaussian:
            
            source_components = chunked_source.gauss_components
            map_components = chunk_map.gauss_components
            
            n_powers = chunk_map.n_gauss_powers
            n_curves = chunk_map.n_gauss_curves
            
        elif comp_info.shapelet:
            
            source_components = chunked_source.shape_components
            map_components = chunk_map.shape_components
            
            n_powers = chunk_map.n_shape_powers
            n_curves = chunk_map.n_shape_curves
            
        ##always shove things into the source as power, curve, list
        ## - chunk_flux_type_index is the index to access with etc
        ##   source_components.power_ref_freqs
        ##   source_components.curve_ref_stokesQ
        ## - chunk_comp_index is the index to access within e.g
        ##   source_components.ras
        
        # print(np.isnan(comp_info.fluxes))
        
        if comp_info.flux_power:
            ##this is the index
            chunk_flux_type_index = int(np.where(map_components.power_orig_inds == orig_comp_ind)[0])
            chunk_comp_index = chunk_flux_type_index
            
            source_components.power_comp_inds[chunk_flux_type_index] = chunk_comp_index
            
            source_components.power_ref_freqs[chunk_flux_type_index] = comp_info.freqs[0]
            source_components.power_ref_stokesI[chunk_flux_type_index] = comp_info.fluxes[0][0]
            source_components.power_ref_stokesQ[chunk_flux_type_index] = comp_info.fluxes[0][1]
            source_components.power_ref_stokesU[chunk_flux_type_index] = comp_info.fluxes[0][2]
            source_components.power_ref_stokesV[chunk_flux_type_index] = comp_info.fluxes[0][3]
            source_components.power_SIs[chunk_flux_type_index] = comp_info.si
            
        if comp_info.flux_curve:
            ##this is the index
            chunk_flux_type_index = int(np.where(map_components.curve_orig_inds == orig_comp_ind)[0])
            chunk_comp_index = n_powers + chunk_flux_type_index
            
            source_components.curve_comp_inds[chunk_flux_type_index] = chunk_comp_index
            
            source_components.curve_ref_freqs[chunk_flux_type_index] = comp_info.freqs[0]
            source_components.curve_ref_stokesI[chunk_flux_type_index] = comp_info.fluxes[0][0]
            source_components.curve_ref_stokesQ[chunk_flux_type_index] = comp_info.fluxes[0][1]
            source_components.curve_ref_stokesU[chunk_flux_type_index] = comp_info.fluxes[0][2]
            source_components.curve_ref_stokesV[chunk_flux_type_index] = comp_info.fluxes[0][3]
            source_components.curve_SIs[chunk_flux_type_index] = comp_info.si
            source_components.curve_qs[chunk_flux_type_index] = comp_info.curve_q
            
        if comp_info.flux_list:
            ##this is the index
            chunk_flux_type_index = int(np.where(map_components.list_orig_inds == orig_comp_ind)[0])
            chunk_comp_index = n_powers + n_curves + chunk_flux_type_index
            
            source_components.list_comp_inds[chunk_flux_type_index] = int(chunk_comp_index)
            
            ##things get complicated with the flux list stuff
            
            if chunk_flux_type_index == 0:
                list_start_index = 0
            else:
                # print("HERE GO", source_components.num_list_values[:chunk_flux_type_index])

                list_start_index = int(np.sum(source_components.num_list_values[:chunk_flux_type_index]))
                
            # print("WHY", comp_info.num_fluxes, list_start_index, chunk_flux_type_index)

            for f_ind in range(comp_info.num_fluxes):

                # print(list_start_index, f_ind, type(list_start_index), type(f_ind))
                
                source_components.list_freqs[list_start_index + f_ind] = comp_info.freqs[f_ind]
                source_components.list_stokesI[list_start_index + f_ind] = comp_info.fluxes[f_ind][0]
                source_components.list_stokesQ[list_start_index + f_ind] = comp_info.fluxes[f_ind][1]
                source_components.list_stokesU[list_start_index + f_ind] = comp_info.fluxes[f_ind][2]
                source_components.list_stokesV[list_start_index + f_ind] = comp_info.fluxes[f_ind][3]

            # print(list_start_index, chunk_flux_type_index, comp_info.num_fluxes)
            source_components.num_list_values[chunk_flux_type_index] = int(comp_info.num_fluxes)
            source_components.list_start_indexes[chunk_flux_type_index] = int(list_start_index)
            
            # print(chunk_flux_type_index,
            #       source_components.num_list_values[chunk_flux_type_index],
            #       source_components.list_start_indexes[chunk_flux_type_index])
            
        ##everybody needs an ra/dec
        source_components.ras[chunk_comp_index] = comp_info.ra
        source_components.decs[chunk_comp_index] = comp_info.dec
        
        num_time_steps = len(lsts)
        
        if beamtype == BeamTypes.GAUSS_BEAM.value or beamtype == BeamTypes.MWA_ANALY.value:
            comp_has = lsts - comp_info.ra
            
            ##OK these ctype arrays cannot be sliced, so let's increment
            ##over them at a snail's pace
            for time_ind in range(num_time_steps):
                hadec_low = chunk_comp_index*num_time_steps
                source_components.beam_has[hadec_low + time_ind] = comp_has[time_ind]
                source_components.beam_decs[hadec_low + time_ind] = comp_info.dec
            
        ##only the NO_BEAM and GAUSS_BEAM options don't need az,za coords
        if beamtype == BeamTypes.GAUSS_BEAM.value or beamtype == BeamTypes.NO_BEAM.value:
            pass
        else:
            ##Calculate lst, and then azimuth/elevation
            comp_has = lsts - comp_info.ra
            comp_azs, comp_els = erfa.hd2ae(comp_has, comp_info.dec,
                                            latitude)
            
            ##OK these ctype arrays cannot be sliced, so let's increment
            ##over them at a snail's pace
            for time_ind in range(num_time_steps):
                azza_low = chunk_comp_index*num_time_steps
                source_components.azs[azza_low + time_ind] = comp_azs[time_ind]
                source_components.zas[azza_low + time_ind] = np.pi/2 - comp_els[time_ind]
        
        ##only some people need major, minor, pas
        if comp_info.gaussian or comp_info.shapelet:
            source_components.majors[chunk_comp_index] = comp_info.major
            source_components.minors[chunk_comp_index] = comp_info.minor
            source_components.pas[chunk_comp_index] = comp_info.pa
            
        ##only shapelets need these other scary things
        if comp_info.shapelet:
        
            if comp_info.flux_power:
                coeff_indexes_to_add = map_components.power_shape_basis_inds[np.where(map_components.power_shape_orig_inds == orig_comp_ind)].astype(int)
                
            elif comp_info.flux_curve:
                coeff_indexes_to_add = map_components.curve_shape_basis_inds[np.where(map_components.curve_shape_orig_inds == orig_comp_ind)].astype(int)
                
            elif comp_info.flux_list:
                coeff_indexes_to_add = map_components.list_shape_basis_inds[np.where(map_components.list_shape_orig_inds == orig_comp_ind)].astype(int)
                
            ##OK these ctype arrays cannot be sliced, so let's increment
            ##over them at a snail's pace
            
            low = chunk_map.current_shape_basis_index
            high = int(chunk_map.current_shape_basis_index + len(coeff_indexes_to_add))
            chunk_basis_indexes = range(low, high)
            
            for chunk_basis_index, coeff_index in zip(chunk_basis_indexes, coeff_indexes_to_add):
            
                source_components.n1s[chunk_basis_index] = comp_info.n1s[coeff_index]
                source_components.n2s[chunk_basis_index] = comp_info.n2s[coeff_index]
                source_components.shape_coeffs[chunk_basis_index] = comp_info.coeffs[coeff_index]
                source_components.param_indexes[chunk_basis_index] = chunk_comp_index
        
            chunk_map.current_shape_basis_index += len(coeff_indexes_to_add)
            
        collected_comps += 1
    
    return collected_comps


class _Components_Python(object):
    """python equivalent to Components_Float or Components_Double"""
    
    def __init__(self, components : Union[Components_Float, Components_Double],
                       comp_type : CompTypes, n_powers : int,
                       n_curves : int, n_lists : int,
                       n_shape_coeffs = 0):
        """
        docstring
        
        """
        
        n_comps = n_powers + n_curves + n_lists
        
        total_num_flux_entires = int(components.total_num_flux_entires)
        
        self.ras = np.ctypeslib.as_array(components.ras, shape=(n_comps, ))
        self.decs = np.ctypeslib.as_array(components.decs, shape=(n_comps, ))
        
        self.power_ref_freqs = np.ctypeslib.as_array(components.power_ref_freqs, shape=(n_powers, ))
        self.power_ref_stokesI = np.ctypeslib.as_array(components.power_ref_stokesI, shape=(n_powers, ))
        self.power_ref_stokesQ = np.ctypeslib.as_array(components.power_ref_stokesQ, shape=(n_powers, ))
        self.power_ref_stokesU = np.ctypeslib.as_array(components.power_ref_stokesU, shape=(n_powers, ))
        self.power_ref_stokesV = np.ctypeslib.as_array(components.power_ref_stokesV, shape=(n_powers, ))
        self.power_SIs = np.ctypeslib.as_array(components.power_SIs, shape=(n_powers, ))
        self.power_comp_inds = np.ctypeslib.as_array(components.power_comp_inds, shape=(n_powers, ))
        
        self.curve_ref_freqs = np.ctypeslib.as_array(components.curve_ref_freqs, shape=(n_curves, ))
        self.curve_ref_stokesI = np.ctypeslib.as_array(components.curve_ref_stokesI, shape=(n_curves, ))
        self.curve_ref_stokesQ = np.ctypeslib.as_array(components.curve_ref_stokesQ, shape=(n_curves, ))
        self.curve_ref_stokesU = np.ctypeslib.as_array(components.curve_ref_stokesU, shape=(n_curves, ))
        self.curve_ref_stokesV = np.ctypeslib.as_array(components.curve_ref_stokesV, shape=(n_curves, ))
        self.curve_SIs = np.ctypeslib.as_array(components.curve_SIs, shape=(n_curves, ))
        self.curve_qs = np.ctypeslib.as_array(components.curve_SIs, shape=(n_curves, ))
        self.curve_comp_inds = np.ctypeslib.as_array(components.curve_comp_inds, shape=(n_curves, ))
        
        # self.list_freqs = np.ctypeslib.as_array(components.list_freqs, shape=(total_num_flux_entires, ))
        # self.list_stokesI = np.ctypeslib.as_array(components.list_stokesI, shape=(total_num_flux_entires, ))
        # self.list_stokesQ = np.ctypeslib.as_array(components.list_stokesQ, shape=(total_num_flux_entires, ))
        # self.list_stokesU = np.ctypeslib.as_array(components.list_stokesU, shape=(total_num_flux_entires, ))
        # self.list_stokesV = np.ctypeslib.as_array(components.list_stokesV, shape=(total_num_flux_entires, ))
        # self.list_comp_inds = np.ctypeslib.as_array(components.list_comp_inds, shape=(n_lists, ))
        
        self.list_freqs = np.empty(total_num_flux_entires, dtype=float)
        self.list_stokesI = np.empty(total_num_flux_entires, dtype=float)
        self.list_stokesQ = np.empty(total_num_flux_entires, dtype=float)
        self.list_stokesU = np.empty(total_num_flux_entires, dtype=float)
        self.list_stokesV = np.empty(total_num_flux_entires, dtype=float)
        
        for ind in range(total_num_flux_entires):
        
            self.list_freqs[ind] = components.list_freqs[ind]
            self.list_stokesI[ind] = components.list_stokesI[ind]
            self.list_stokesQ[ind] = components.list_stokesQ[ind]
            self.list_stokesU[ind] = components.list_stokesU[ind]
            self.list_stokesV[ind] = components.list_stokesV[ind]
        
        self.list_comp_inds = np.empty(n_lists, dtype=int)
        self.num_list_values = np.empty(n_lists, dtype=int)
        self.list_start_indexes = np.empty(n_lists, dtype=int)
        
        for ind in range(n_lists):
            self.num_list_values[ind] = components.num_list_values[ind]
            self.list_start_indexes[ind] = components.list_start_indexes[ind]
            self.list_comp_inds[ind] = components.list_comp_inds[ind]
        
        self.num_list_values = np.ctypeslib.as_array(components.num_list_values, shape=(int(n_lists), ))
        self.list_start_indexes = np.ctypeslib.as_array(components.list_start_indexes, shape=(n_lists, ))
        
        if comp_type == CompTypes.GAUSSIAN or comp_type == CompTypes.SHAPELET:
            self.minors = np.ctypeslib.as_array(components.minors, shape=(n_comps, ))
            self.majors = np.ctypeslib.as_array(components.majors, shape=(n_comps, ))
            self.pas = np.ctypeslib.as_array(components.pas, shape=(n_comps, ))
            
        if comp_type == CompTypes.SHAPELET:
            self.param_indexes = np.ctypeslib.as_array(components.param_indexes, shape=(n_shape_coeffs, ))
            self.n1s = np.ctypeslib.as_array(components.n1s, shape=(n_shape_coeffs, ))
            self.n2s = np.ctypeslib.as_array(components.n2s, shape=(n_shape_coeffs, ))
            self.shape_coeffs = np.ctypeslib.as_array(components.shape_coeffs, shape=(n_shape_coeffs, ))
        
    
    

class _Ctype_Source_Into_Python(object):
    """
    Class to convert a ctype Source model into a pythonic version
    """
    
    def __init__(self, source : Union[Source_Float, Source_Double]):
        """
        docstring
        """
        
        self.n_points = source.n_points
        self.n_point_lists = source.n_point_lists
        self.n_point_powers = source.n_point_powers
        self.n_point_curves = source.n_point_curves
        self.n_gauss = source.n_gauss
        self.n_gauss_lists = source.n_gauss_lists
        self.n_gauss_powers = source.n_gauss_powers
        self.n_gauss_curves = source.n_gauss_curves
        self.n_shapes = source.n_shapes
        self.n_shape_lists = source.n_shape_lists
        self.n_shape_powers = source.n_shape_powers
        self.n_shape_curves = source.n_shape_curves
        self.n_shape_coeffs = source.n_shape_coeffs
        self.n_comps = source.n_comps
        
        if source.n_points:
            self.point_components = _Components_Python(source.point_components,
                                                       CompTypes.POINT,
                                                       self.n_point_powers,
                                                       self.n_point_curves,
                                                       self.n_point_lists)
            
        if source.n_gauss:
            self.gauss_components = _Components_Python(source.gauss_components,
                                                       CompTypes.GAUSSIAN,
                                                       self.n_gauss_powers,
                                                       self.n_gauss_curves,
                                                       self.n_gauss_lists)
            
        if source.n_shapes:
            self.shape_components = _Components_Python(source.shape_components,
                                                       CompTypes.SHAPELET,
                                                       self.n_shape_powers,
                                                       self.n_shape_curves,
                                                       self.n_shape_lists,
                                                       self.n_shape_coeffs)
            
            
    # num_cross = num_time_steps * num_freq_channels * num_baselines

    ##Righto, this converts from the ctype POINTER into a numpy array
    ##This is grabbing all the lovely things calculated by the GPU
    # us_metres = np.ctypeslib.as_array(visibility_set.us_metres, shape=(num_visi,))
 