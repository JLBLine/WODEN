import ctypes 
import importlib_resources
import numpy as np
from typing import Union
from ctypes import POINTER, c_float, c_double, c_int, c_char, pointer
import os
import sys
import erfa

from wodenpy.skymodel.woden_skymodel import Component_Type_Counter, CompTypes, Component_Info
from wodenpy.skymodel.chunk_sky_model import Skymodel_Chunk_Map
from wodenpy.use_libwoden.beam_settings import BeamTypes

VELC = 299792458.0
D2R = np.pi / 180.0


def create_components_struct(precision="double"):
    
    if precision == "float":
        c_user_precision = c_float
    else:
        c_user_precision = c_double
    
    class Components(ctypes.Structure):
        """A class structured equivalently to a `components_t` struct, used by 
        the C and CUDA code in libwoden_float.so
        
        :cvar POINTER(double) ras: COMPONENT right ascensions (radians)
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
        :cvar POINTER(int) total_num_flux_entires: The total number of freq/flux values are in all lists combine
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
        :cvar POINTER(user_precision_complex_t) gxs_ants: North-South Beam gain values for all directions, frequencies, times, and ants for these COMPONENT
        :cvar POINTER(user_precision_complex_t) Dxs_ants: North-South Beam leakage values for all directions, frequencies, times, and ants for these COMPONENT
        :cvar POINTER(user_precision_complex_t) Dys_ants: East-West Beam leakage values for all directions, frequencies, times, and ants  for these COMPONENT
        :cvar POINTER(user_precision_complex_t) gys_ants: East-West Beam gain values for all directions, frequencies, times, and ants  for these COMPONENT
        :cvar POINTER(double) ls: Device memory l cosine direction coords for these COMPONENTs
        :cvar POINTER(double) ms: Device memory m cosine direction coords for these COMPONENTs
        :cvar POINTER(double) ns: Device memory n cosine direction coords for these COMPONENTs
        :cvar POINTER(c_float) stokesV_pol_fracs:  Stokes V polarisation fractions
        :cvar POINTER(c_int) stokesV_pol_frac_comp_inds: The indexes of all Stokes V polarisation fraction models w.r.t ra,dec
        :cvar POINTER(c_float) stokesV_power_ref_flux: Stokes V reference flux for power-law
        :cvar POINTER(c_float) stokesV_power_SIs: Stokes V spectral index for power-law
        :cvar POINTER(c_int) stokesV_power_comp_inds: The indexes of all Stokes V power-law models w.r.t ra,dec
        :cvar POINTER(c_float) stokesV_curve_ref_flux: Stokes V reference flux for curved power-law
        :cvar POINTER(c_float) stokesV_curve_SIs: Stokes V spectral index for curved power-law
        :cvar POINTER(c_float) stokesV_curve_qs: Stokes V q param for curved power-law
        :cvar POINTER(c_int) stokesV_curve_comp_inds: The indexes of Stokes V curved power-law models w.r.t ra,dec
        :cvar POINTER(c_float) linpol_pol_fracs: Linear polarisation polarisation fractions
        :cvar POINTER(c_int) linpol_pol_frac_comp_inds: The indexes of all linear polarisation fraction models w.r.t ra,dec
        :cvar POINTER(c_float) linpol_power_ref_flux: Linear polarisation reference flux for power-law
        :cvar POINTER(c_float) linpol_power_SIs: Linear polarisation spectral index for power-law
        :cvar POINTER(c_int) linpol_power_comp_inds: The indexes of all linear polarisation power-law models w.r.t ra,dec
        :cvar POINTER(c_float) linpol_curve_ref_flux: Linear polarisation reference flux for curved power-law
        :cvar POINTER(c_float) linpol_curve_SIs: Linear polarisation spectral index for curved power-law
        :cvar POINTER(c_float) linpol_curve_qs: Linear polarisation q param for curved power-law
        :cvar POINTER(c_int) linpol_curve_comp_inds: The indexes of all linear polarisation curved power-law models w.r.t ra,dec
        :cvar POINTER(c_float) rm_values: Linear polarisation rotation measures
        :cvar POINTER(c_float) intr_pol_angle: Linear polarisation instrinsic polarisation angles
        :cvar POINTER(c_int) linpol_angle_inds: The indexes of all RM/intrinsic polarisation angles w.r.t ra,dec
        :cvar c_int n_stokesV_pol_frac: The number of Stokes V polarisation fraction models
        :cvar c_int n_stokesV_power: The number of Stokes V power-law models
        :cvar c_int n_stokesV_curve: The number of V curved power-law models
        :cvar c_int n_linpol_pol_frac: The number of linear polarisation fraction models
        :cvar c_int n_linpol_power: The number of linear polarisation power-law models
        :cvar c_int n_linpol_curve: The number of linear polarisation curved power-law models
        :cvar c_int n_linpol_angles: The number of RM/intrinsic polarisation angles
        :cvar c_int do_QUV: Set if doing any polarised information
        
        """
        
        _fields_ = [## Instrinsic to COMPONENT values
                    ("ras", POINTER(c_double)),
                    ("decs", POINTER(c_double)),
                    ## power law params
                    ("power_ref_freqs", POINTER(c_double)),
                    ("power_ref_stokesI", POINTER(c_user_precision)),
                    ("power_SIs", POINTER(c_user_precision)),
                    ##curved power law params
                    ("curve_ref_freqs", POINTER(c_double)),
                    ("curve_ref_stokesI", POINTER(c_user_precision)),
                    ("curve_SIs", POINTER(c_user_precision)),
                    ("curve_qs", POINTER(c_user_precision)),
                    ##indexes of types
                    ("power_comp_inds", POINTER(c_int)),
                    ("curve_comp_inds", POINTER(c_int)),
                    ("list_comp_inds", POINTER(c_int)),
                    ##list flux params
                    ("list_freqs", POINTER(c_double)),
                    ("list_stokesI", POINTER(c_user_precision)),
                    ("list_stokesQ", POINTER(c_user_precision)),
                    ("list_stokesU", POINTER(c_user_precision)),
                    ("list_stokesV", POINTER(c_user_precision)),
                    ("num_list_values", POINTER(c_int)),
                    ("list_start_indexes", POINTER(c_int)),
                    ("total_num_flux_entires", c_int),
                    ##something to store extrapolated output fluxes in
                    ("extrap_stokesI", POINTER(c_user_precision)),
                    ("extrap_stokesQ", POINTER(c_user_precision)),
                    ("extrap_stokesU", POINTER(c_user_precision)),
                    ("extrap_stokesV", POINTER(c_user_precision)),
                    ##SHAPELET params
                    ("shape_coeffs", POINTER(c_user_precision)),
                    ("n1s", POINTER(c_user_precision)),
                    ("n2s", POINTER(c_user_precision)),
                    ##SHAPELET / GAUSSIAN params
                    ("majors", POINTER(c_user_precision)),
                    ("minors", POINTER(c_user_precision)),
                    ("pas", POINTER(c_user_precision)),
                    ("param_indexes", POINTER(c_user_precision)),
                    ##Specific to observation settings for these COMPONENTs
                    ("azs", POINTER(c_user_precision)),
                    ("zas", POINTER(c_user_precision)),
                    ("beam_has", POINTER(c_double)),
                    ("beam_decs", POINTER(c_double)),
                    ("num_primarybeam_values", c_int),
                    ##things to hold the beam gain
                    ##these are _complex in the C struct; we are not allocating
                    ##memory on the python side, so I think we can just call them
                    ##a double here?
                    ("gxs", POINTER(c_user_precision)),
                    ("Dxs", POINTER(c_user_precision)),
                    ("Dys", POINTER(c_user_precision)),
                    ("gys", POINTER(c_user_precision)),
                    ("gxs_ants", POINTER(c_user_precision)),
                    ("Dxs_ants", POINTER(c_user_precision)),
                    ("Dys_ants", POINTER(c_user_precision)),
                    ("gys_ants", POINTER(c_user_precision)),
                    ##used to hold l,m,n coords on the GPU
                    ("ls", POINTER(c_double)),
                    ("ms", POINTER(c_double)),
                    ("ns", POINTER(c_double)),
                    ###polarisation information
                    ("stokesV_pol_fracs", POINTER(c_user_precision)),
                    ("stokesV_pol_frac_comp_inds", POINTER(c_int)),
                    ("stokesV_power_ref_flux", POINTER(c_user_precision)),
                    ("stokesV_power_SIs", POINTER(c_user_precision)),
                    ("stokesV_power_comp_inds", POINTER(c_int)),
                    ("stokesV_curve_ref_flux", POINTER(c_user_precision)),
                    ("stokesV_curve_SIs", POINTER(c_user_precision)),
                    ("stokesV_curve_qs", POINTER(c_user_precision)),
                    ("stokesV_curve_comp_inds", POINTER(c_int)),
                    ("linpol_pol_fracs", POINTER(c_user_precision)),
                    ("linpol_pol_frac_comp_inds", POINTER(c_int)),
                    ("linpol_power_ref_flux", POINTER(c_user_precision)),
                    ("linpol_power_SIs", POINTER(c_user_precision)),
                    ("linpol_power_comp_inds", POINTER(c_int)),
                    ("linpol_curve_ref_flux", POINTER(c_user_precision)),
                    ("linpol_curve_SIs", POINTER(c_user_precision)),
                    ("linpol_curve_qs", POINTER(c_user_precision)),
                    ("linpol_curve_comp_inds", POINTER(c_int)),
                    ("rm_values", POINTER(c_user_precision)),
                    ("intr_pol_angle", POINTER(c_user_precision)),
                    ("linpol_angle_inds", POINTER(c_int)),
                    ("n_stokesV_pol_frac", c_int),
                    ("n_stokesV_power", c_int),
                    ("n_stokesV_curve", c_int),
                    ("n_stokesV_list", c_int),
                    ("n_stokesV_list_flux_entries", c_int),
                    ("n_linpol_pol_frac", c_int),
                    ("n_linpol_power", c_int),
                    ("n_linpol_curve", c_int),
                    ("n_linpol_list", c_int),
                    ("n_stokesQ_list_flux_entries", c_int),
                    ("n_stokesU_list_flux_entries", c_int),
                    ("n_linpol_p_list", c_int),
                    ("n_linpol_p_list_flux_entries", c_int),
                    ("n_linpol_angles", c_int),
                    ("do_QUV", c_int),
                    ]
        
    return Components
    
    
def create_source_struct(Components):
    
    class Source(ctypes.Structure):
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
            ("point_components",  Components),
            ("gauss_components",  Components),
            ("shape_components",  Components),
            ("d_point_components",  Components),
            ("d_gauss_components",  Components),
            ("d_shape_components",  Components),
        ]
        
    return Source
    
# class Source_Catalogue_Double(ctypes.Structure):
#     """
#     A class structured equivalent to a `source_t` struct, used by 
#     the C and CUDA code in libwoden_double.so
    
#     Attributes
#     -----------
#     num_sources : int
#         The number of sources in the catalogue
#     num_shapelets : int
#         The total number of shapelets components in the catalogue
#     sources : POINTER(Source_Double)
#         A pointer to an array of `Source_Double` objects representing the sources in the catalogue
#     """
    
#     _fields_ = [("num_sources", c_int),
#                 ("num_shapelets", c_int),
#                 ("sources", POINTER(Source_Double))]


def create_source_catalogue_struct(Source):
    
    # Source = create_source_struct(precision)
    
    class Source_Catalogue(ctypes.Structure):
        """
        A class structured equivalent to a `source_t` struct, used by 
        the C and CUDA code in libwoden_float.so
        
        Attributes
        -----------
        num_sources : int
            The number of sources in the catalogue
        num_shapelets : int
            The total number of shapelets components in the catalogue
        sources : POINTER(Source_Float)
            A pointer to an array of `Source_Float` objects representing the sources in the catalogue
        """
        
        # _fields_ = [("num_sources", c_int),
        #             ("num_shapelets", c_int),
        #             ("sources", POINTER(Source))]
        
        _fields_ = [("num_sources", c_int),
                    ("num_shapelets", c_int),
                    ("sources", POINTER(Source))]
        
    return Source_Catalogue
    
    
def setup_source_catalogue(Source, Source_Catalogue, num_sources : int, num_shapelets : int,
                           precision = "double"):
    """
    Creates a source catalogue object with the specified number of sources and shapelets.
    
    Parameters
    ------------
    num_sources: int
        The number of sources in the catalogue.
    num_shapelets: int
        The number of shapelets for each source.
    precision: str
        The precision of the source catalogue. Can be "float" or "double". Default is "double".
    
    Returns
    ---------
    source_catalogue : (Source_Catalogue_Float or Source_Catalogue_Double)
        The initialised source_catalogue object.
    """
    
    source_catalogue = Source_Catalogue()
    
    source_array = num_sources*Source
        
    source_catalogue.sources = source_array()
    source_catalogue.num_sources = num_sources
    source_catalogue.num_shapelets = num_shapelets
    
    return source_catalogue
    
def setup_components(woden_struct_classes,
                     chunk_map : Skymodel_Chunk_Map,
                     chunked_source,
                     num_freqs : int,
                     num_times : int, comp_type : CompTypes,
                     beamtype : int, 
                     c_user_precision : Union[c_float, c_double]):
    """
    Given the mapping information in `chunk_map`, initialise the necessary
    components in `chunked_source` for the given `comp_type`, this being one
    of
     - chunked_source.point_components
     - chunked_source.gaussion_components
     - chunked_source.shape_components
    This setup includes allocating memory.

    Parameters
    ------------
    chunk_map: Skymodel_Chunk_Map
        Object containing information about the chunk of the sky model.
    chunked_source: Union
       [Source_Float, Source_Double] object containing the chunked source information.
    num_freqs: int
        representing the number of frequencies.
    num_times: int
        representing the number of times.
    comp_type: CompTypes
        enum representing the type of component.
    beamtype: int
        representing the type of beam.
    c_user_precision: Union
       [c_float, c_double] representing the user precision.


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
        
        map_components = chunk_map.point_components
        
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
        
        map_components = chunk_map.gauss_components
        
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
        
        map_components = chunk_map.shape_components
        
    ##component type specific things
    user_ncomps_arr = c_user_precision*n_comps
    double_ncomps_arr = c_double*n_comps
    num_primarybeam_values = n_comps*num_freqs*num_times
        
    components.ras = double_ncomps_arr()
    components.decs = double_ncomps_arr()
    components.num_primarybeam_values = num_primarybeam_values
    
    ##power-law flux things
    components.power_ref_freqs = power_double_ncomps_arr() ##this is always double
    
    # print("DIS", type(components.power_ref_stokesI))
    
    components.power_ref_stokesI = power_user_ncomps_arr()
    components.power_SIs = power_user_ncomps_arr()
    
    ##curved power-law flux things
    components.curve_ref_freqs = curve_double_ncomps_arr() ##this is always double
    components.curve_ref_stokesI = curve_user_ncomps_arr()
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
        
    ##now do the polarisation thingies------------------------------------------
    
    power_user_ncomps_arr = c_user_precision*chunk_map.n_shape_powers
    
    components.stokesV_pol_fracs = (c_user_precision*map_components.num_v_pol_fracs)()
    components.stokesV_pol_frac_comp_inds = (c_int*map_components.num_v_pol_fracs)()
    components.stokesV_power_ref_flux = (c_user_precision*map_components.num_v_powers)()
    components.stokesV_power_SIs = (c_user_precision*map_components.num_v_powers)()
    components.stokesV_power_comp_inds = (c_int*map_components.num_v_powers)()
    components.stokesV_curve_ref_flux = (c_user_precision*map_components.num_v_curves)()
    components.stokesV_curve_SIs = (c_user_precision*map_components.num_v_curves)()
    components.stokesV_curve_qs = (c_user_precision*map_components.num_v_curves)()
    components.stokesV_curve_comp_inds = (c_int*map_components.num_v_curves)()
    components.linpol_pol_fracs = (c_user_precision*map_components.num_lin_pol_fracs)()
    components.linpol_pol_frac_comp_inds = (c_int*map_components.num_lin_pol_fracs)()
    components.linpol_power_ref_flux = (c_user_precision*map_components.num_lin_powers)()
    components.linpol_power_SIs = (c_user_precision*map_components.num_lin_powers)()
    components.linpol_power_comp_inds = (c_int*map_components.num_lin_powers)()
    components.linpol_curve_ref_flux = (c_user_precision*map_components.num_lin_curves)()
    components.linpol_curve_SIs = (c_user_precision*map_components.num_lin_curves)()
    components.linpol_curve_qs = (c_user_precision*map_components.num_lin_curves)()
    components.linpol_curve_comp_inds = (c_int*map_components.num_lin_curves)()
    components.rm_values = (c_user_precision*map_components.num_lin_angles)()
    components.intr_pol_angle = (c_user_precision*map_components.num_lin_angles)()
    components.linpol_angle_inds = (c_int*map_components.num_lin_angles)()
    
    components.n_stokesV_pol_frac = map_components.num_v_pol_fracs
    components.n_stokesV_power = map_components.num_v_powers
    components.n_stokesV_curve = map_components.num_v_curves
    components.n_linpol_pol_frac = map_components.num_lin_pol_fracs
    components.n_linpol_power = map_components.num_lin_powers
    components.n_linpol_curve = map_components.num_lin_curves
    components.n_linpol_angles = map_components.num_lin_angles
    
    # print(components.n_stokesV_pol_frac)
    # print(components.n_stokesV_power)
    # print(components.n_stokesV_curve)
    # print(components.n_linpol_pol_frac)
    # print(components.n_linpol_power)
    # print(components.n_linpol_curve)
    # print(components.n_linpol_angles)
        
    
def setup_chunked_source(woden_struct_classes, chunked_source, chunk_map : Skymodel_Chunk_Map, num_freqs : int,
                         num_times : int, beamtype : int,
                         precision='double'):
    """
    Sets up a ctypes structure class to contain a chunked sky model.
    This class is compatible with the C/CUDA code, and will allocate the
    correct amount of memory, based on whether the precision is either
    'double' or 'float'.

    Parameters
    ------------
    chunk_map : Skymodel_Chunk_Map
        Object containing information about the chunk of the sky model.
    num_freqs : int
        The number of frequency channels.
    num_times : int
        The number of time samples.
    beamtype : int
        The type of primary beam.
    precision : str, optional
        The precision to use. Can be either "float" or "double. Defaults to "double".

    Returns
    ---------
      Source_Float or Source_Double
          The ctypes structure class containing the chunked sky model.
    """
    
    # set_precision_global(precision)
    
    if precision == 'float':
        c_user_precision = c_float
    else:
        c_user_precision = c_double
    
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
        setup_components(woden_struct_classes, chunk_map, chunked_source, num_freqs, num_times,
                         CompTypes.POINT, beamtype, c_user_precision)
        
    if chunk_map.n_gauss > 0:
        setup_components(woden_struct_classes, chunk_map, chunked_source, num_freqs, num_times,
                         CompTypes.GAUSSIAN, beamtype, c_user_precision)
        
    if chunk_map.n_shapes > 0:
        setup_components(woden_struct_classes, chunk_map, chunked_source, num_freqs, num_times,
                         CompTypes.SHAPELET, beamtype, c_user_precision)
    
    return chunked_source


class _Components_Python(object):
    """python equivalent to Components_Float or Components_Double"""
    
    def __init__(self, components,
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
        self.power_SIs = np.ctypeslib.as_array(components.power_SIs, shape=(n_powers, ))
        self.power_comp_inds = np.ctypeslib.as_array(components.power_comp_inds, shape=(n_powers, ))
        
        self.curve_ref_freqs = np.ctypeslib.as_array(components.curve_ref_freqs, shape=(n_curves, ))
        self.curve_ref_stokesI = np.ctypeslib.as_array(components.curve_ref_stokesI, shape=(n_curves, ))
        self.curve_SIs = np.ctypeslib.as_array(components.curve_SIs, shape=(n_curves, ))
        self.curve_qs = np.ctypeslib.as_array(components.curve_qs, shape=(n_curves, ))
        self.curve_comp_inds = np.ctypeslib.as_array(components.curve_comp_inds, shape=(n_curves, ))
        
        # self.list_freqs = np.ctypeslib.as_array(components.list_freqs, shape=(total_num_flux_entires, ))
        
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
            
            
        # self.curve_comp_inds = np.ctypeslib.as_array(components.curve_comp_inds, shape=(n_curves, ))
        
        self.n_stokesV_pol_frac = 0
        self.n_stokesV_power = 0
        self.n_stokesV_curve = 0
        self.n_linpol_pol_frac = 0
        self.n_linpol_power = 0
        self.n_linpol_curve = 0
        self.n_linpol_angles = 0
        
        if components.n_stokesV_pol_frac > 0:
            self.n_stokesV_pol_frac = components.n_stokesV_pol_frac
            self.stokesV_pol_fracs = np.ctypeslib.as_array(components.stokesV_pol_fracs, shape=(components.n_stokesV_pol_frac, ))
            self.stokesV_pol_frac_comp_inds = np.ctypeslib.as_array(components.stokesV_pol_frac_comp_inds, shape=(components.n_stokesV_pol_frac, ))
            
            
        if components.n_stokesV_power > 0:
            self.n_stokesV_power = components.n_stokesV_power
            self.stokesV_power_ref_flux = np.ctypeslib.as_array(components.stokesV_power_ref_flux, shape=(components.n_stokesV_power, ))
            self.stokesV_power_SIs = np.ctypeslib.as_array(components.stokesV_power_SIs, shape=(components.n_stokesV_power, ))
            self.stokesV_power_comp_inds = np.ctypeslib.as_array(components.stokesV_power_comp_inds, shape=(components.n_stokesV_power, ))
            
        if components.n_stokesV_curve > 0:
            self.n_stokesV_curve = components.n_stokesV_curve
            self.stokesV_curve_ref_flux = np.ctypeslib.as_array(components.stokesV_curve_ref_flux, shape=(components.n_stokesV_curve, ))
            self.stokesV_curve_SIs = np.ctypeslib.as_array(components.stokesV_curve_SIs, shape=(components.n_stokesV_curve, ))
            self.stokesV_curve_qs = np.ctypeslib.as_array(components.stokesV_curve_qs, shape=(components.n_stokesV_curve, ))
            self.stokesV_curve_comp_inds = np.ctypeslib.as_array(components.stokesV_curve_comp_inds, shape=(components.n_stokesV_curve, ))
            
        if components.n_linpol_pol_frac > 0:
            self.n_linpol_pol_frac = components.n_linpol_pol_frac
            self.linpol_pol_fracs = np.ctypeslib.as_array(components.linpol_pol_fracs, shape=(components.n_linpol_pol_frac, ))
            self.linpol_pol_frac_comp_inds = np.ctypeslib.as_array(components.linpol_pol_frac_comp_inds, shape=(components.n_linpol_pol_frac, ))
            
            
        if components.n_linpol_power > 0:
            self.n_linpol_power = components.n_linpol_power
            self.linpol_power_ref_flux = np.ctypeslib.as_array(components.linpol_power_ref_flux, shape=(components.n_linpol_power, ))
            self.linpol_power_SIs = np.ctypeslib.as_array(components.linpol_power_SIs, shape=(components.n_linpol_power, ))
            self.linpol_power_comp_inds = np.ctypeslib.as_array(components.linpol_power_comp_inds, shape=(components.n_linpol_power, ))
            
            
        if components.n_linpol_curve > 0:
            self.n_linpol_curve = components.n_linpol_curve
            self.linpol_curve_ref_flux = np.ctypeslib.as_array(components.linpol_curve_ref_flux, shape=(components.n_linpol_curve, ))
            self.linpol_curve_SIs = np.ctypeslib.as_array(components.linpol_curve_SIs, shape=(components.n_linpol_curve, ))
            self.linpol_curve_qs = np.ctypeslib.as_array(components.linpol_curve_qs, shape=(components.n_linpol_curve, ))
            self.linpol_curve_comp_inds = np.ctypeslib.as_array(components.linpol_curve_comp_inds, shape=(components.n_linpol_curve, ))
            
        if components.n_linpol_angles > 0:
            self.n_linpol_angles = components.n_linpol_angles
            self.rm_values = np.ctypeslib.as_array(components.rm_values, shape=(components.n_linpol_angles, ))
            self.intr_pol_angle = np.ctypeslib.as_array(components.intr_pol_angle, shape=(components.n_linpol_angles, ))
            self.linpol_angle_inds = np.ctypeslib.as_array(components.linpol_angle_inds, shape=(components.n_linpol_angles, ))
            
        self.do_QUV = components.do_QUV
            

class _Ctype_Source_Into_Python(object):
    """
    Class to convert a ctype Source model into a pythonic version
    """
    
    # def __init__(self, source : Union[Source_Float, Source_Double]):
    def __init__(self, source):
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
 