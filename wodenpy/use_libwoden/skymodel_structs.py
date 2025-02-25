import ctypes 
import importlib_resources
import numpy as np
from typing import Union
from ctypes import POINTER, c_float, c_double, c_int, c_char, pointer, Structure
import os
import sys
import erfa

from wodenpy.skymodel.woden_skymodel import Component_Type_Counter, CompTypes, Component_Info
from wodenpy.skymodel.chunk_sky_model import Skymodel_Chunk_Map
from wodenpy.use_libwoden.beam_settings import BeamTypes, BeamGroups

VELC = 299792458.0
D2R = np.pi / 180.0

class c_double_complex(Structure): 
    """ctypes doesn't have a built-in complex type, define one here.
    Thanks to random internet person for this code.
    https://stackoverflow.com/questions/13373291/complex-number-in-ctypes
    
    Note you could use python complex types, and some wrapper functions
    inside C, but I'd keep C code as simple as possible
    
    """
    _fields_ = [("real", c_double),("imag", c_double)]
    # @property
    # def value(self):
    #     return self.real+1j*self.imag # fields declared above
    
class c_float_complex(Structure): 
    """ctypes doesn't have a built-in complex type, define one here.
    Thanks to random internet person for this code.
    https://stackoverflow.com/questions/13373291/complex-number-in-ctypes
    
    Note you could use python complex types, and some wrapper functions
    inside C, but I'd keep C code as simple as possible
    
    """
    _fields_ = [("real", c_float),("imag", c_float)]

def create_components_struct(precision="double"):
    """Creates a `Components_Ctypes` class structured equivalently to a `components_t`
    struct in the C/CUDA code. Created dynamically based on the `precision`,
    to match the compile time precision flag `-DUSE_DOUBLE` in the C code.

    Parameters
    ----------
    precision : str, optional
        Either "float" or "double:, by default "double"

    Returns
    -------
    Components_Ctypes
        The Components_Ctypes class structured equivalently to a `components_t` struct
    """
    
    if precision == "float":
        c_user_precision = c_float
        c_user_precision_complex = c_float_complex
    else:
        c_user_precision = c_double
        c_user_precision_complex = c_double_complex
    
    class Components_Ctypes(Structure):
        """A class structured equivalently to a `components_t` struct, used by 
        the C and CUDA code in libwoden_float.so or libwoden_double.so.
        
        Created by the function `create_components_struct`, which sets
        `user_precision_t` to either `c_float` or `c_double`.
        
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
        :cvar POINTER(user_precision_t) azs: Azimuth angles for all time steps
        :cvar POINTER(user_precision_t) zas: Zenith angles for all time steps
        :cvar POINTER(user_precision_t) para_angles: Parallactic angles for all time steps
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
                    ("para_angles", POINTER(c_user_precision)),
                    ("beam_has", POINTER(c_double)),
                    ("beam_decs", POINTER(c_double)),
                    ("num_primarybeam_values", c_int),
                    ##things to hold the beam gain
                    ("gxs", POINTER(c_user_precision_complex)),
                    ("Dxs", POINTER(c_user_precision_complex)),
                    ("Dys", POINTER(c_user_precision_complex)),
                    ("gys", POINTER(c_user_precision_complex)),
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
                    ("stokesV_list_ref_freqs", POINTER(c_double)),
                    ("stokesV_list_ref_flux", POINTER(c_user_precision)),
                    ("stokesV_list_comp_inds", POINTER(c_int)),
                    ("stokesV_num_list_values", POINTER(c_int)),
                    ("stokesV_list_start_indexes", POINTER(c_int)),
                    ("linpol_pol_fracs", POINTER(c_user_precision)),
                    ("linpol_pol_frac_comp_inds", POINTER(c_int)),
                    ("linpol_power_ref_flux", POINTER(c_user_precision)),
                    ("linpol_power_SIs", POINTER(c_user_precision)),
                    ("linpol_power_comp_inds", POINTER(c_int)),
                    ("linpol_curve_ref_flux", POINTER(c_user_precision)),
                    ("linpol_curve_SIs", POINTER(c_user_precision)),
                    ("linpol_curve_qs", POINTER(c_user_precision)),
                    ("linpol_curve_comp_inds", POINTER(c_int)),
                    ("stokesQ_list_ref_freqs", POINTER(c_double)),
                    ("stokesQ_list_ref_flux", POINTER(c_user_precision)),
                    ("stokesQ_list_comp_inds", POINTER(c_int)),
                    ("stokesQ_num_list_values", POINTER(c_int)),
                    ("stokesQ_list_start_indexes", POINTER(c_int)),
                    ("stokesU_list_ref_freqs", POINTER(c_double)),
                    ("stokesU_list_ref_flux", POINTER(c_user_precision)),
                    ("stokesU_list_comp_inds", POINTER(c_int)),
                    ("stokesU_num_list_values", POINTER(c_int)),
                    ("stokesU_list_start_indexes", POINTER(c_int)),
                    ("linpol_p_list_ref_freqs", POINTER(c_double)),
                    ("linpol_p_list_ref_flux", POINTER(c_user_precision)),
                    ("linpol_p_list_comp_inds", POINTER(c_int)),
                    ("linpol_p_num_list_values", POINTER(c_int)),
                    ("linpol_p_list_start_indexes", POINTER(c_int)),
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
        
    return Components_Ctypes

class Components_Python(object):
    """A class structured equivalently to a `components_t` struct but using
    Python native data types. This can be pickled and therefore fed to Python
    multiprocessing functions. Needed as `ctypes` pointers cannot be pickled,
    and therefore cannot be passed to Python multiprocessing functions.
    
    """
    
    def __init__(self) -> None:
        self.ras = None
        self.decs = None
        self.power_ref_freqs = None
        self.power_ref_stokesI = None
        self.power_SIs = None
        self.curved = None
        self.curve_ref_freqs = None
        self.curve_ref_stokesI = None
        self.curve_SIs = None
        self.curve_qs = None
        self.indexes = None
        self.power_comp_inds = None
        self.curve_comp_inds = None
        self.list_comp_inds = None
        self.list = None
        self.list_freqs = None
        self.list_stokesI = None
        self.list_stokesQ = None
        self.list_stokesU = None
        self.list_stokesV = None
        self.num_list_values = None
        self.list_start_indexes = None
        self.total_num_flux_entires = None
        self.something = None
        self.extrap_stokesI = None
        self.extrap_stokesQ = None
        self.extrap_stokesU = None
        self.extrap_stokesV = None
        self.shape_coeffs = None
        self.n1s = None
        self.n2s = None
        self.majors = None
        self.minors = None
        self.pas = None
        self.param_indexes = None
        self.azs = None
        self.zas = None
        self.para_angles = None
        self.beam_has = None
        self.beam_decs = None
        self.num_primarybeam_values = None
        self.gxs = None
        self.Dxs = None
        self.Dys = None
        self.gys = None
        self.stokesV_pol_fracs = None
        self.stokesV_pol_frac_comp_inds = None
        self.stokesV_power_ref_flux = None
        self.stokesV_power_SIs = None
        self.stokesV_power_comp_inds = None
        self.stokesV_curve_ref_flux = None
        self.stokesV_curve_SIs = None
        self.stokesV_curve_qs = None
        self.stokesV_curve_comp_inds = None
        self.stokesV_list_ref_freqs = None
        self.stokesV_list_ref_flux = None
        self.stokesV_list_comp_inds = None
        self.stokesV_num_list_values = None
        self.stokesV_list_start_indexes = None
        self.linpol_pol_fracs = None
        self.linpol_pol_frac_comp_inds = None
        self.linpol_power_ref_flux = None
        self.linpol_power_SIs = None
        self.linpol_power_comp_inds = None
        self.linpol_curve_ref_flux = None
        self.linpol_curve_SIs = None
        self.linpol_curve_qs = None
        self.linpol_curve_comp_inds = None
        self.stokesQ_list_ref_freqs = None
        self.stokesQ_list_ref_flux = None
        self.stokesQ_list_comp_inds = None
        self.stokesQ_num_list_values = None
        self.stokesQ_list_start_indexes = None
        self.stokesU_list_ref_freqs = None
        self.stokesU_list_ref_flux = None
        self.stokesU_list_comp_inds = None
        self.stokesU_num_list_values = None
        self.stokesU_list_start_indexes = None
        self.linpol_p_list_ref_freqs = None
        self.linpol_p_list_ref_flux = None
        self.linpol_p_list_comp_inds = None
        self.linpol_p_num_list_values = None
        self.linpol_p_list_start_indexes = None
        self.rm_values = None
        self.intr_pol_angle = None
        self.linpol_angle_inds = None
        self.n_stokesV_pol_frac = None
        self.n_stokesV_power = None
        self.n_stokesV_curve = None
        self.n_stokesV_list = None
        self.n_stokesV_list_flux_entries = None
        self.n_linpol_pol_frac = None
        self.n_linpol_power = None
        self.n_linpol_curve = None
        self.n_linpol_list = None
        self.n_stokesQ_list_flux_entries = None
        self.n_stokesU_list_flux_entries = None
        self.n_linpol_p_list = None
        self.n_linpol_p_list_flux_entries = None
        self.n_linpol_angles = None
        self.do_QUV = None


##This call is so we can use it as a type annotation, and so sphinx can document the class
Components_Ctypes = create_components_struct("double")
    
def create_source_struct(Components_Ctypes : Components_Ctypes): # type: ignore
    """Creates a `Source_Ctypes` class structured equivalent to a `source_t` struct,
    used by the C and CUDA code in libwoden_double.so. The `Components_Ctypes` class
    is passed in as the `Components_Ctypes` class is created dynamically depending on
    the required precision.
    
    Parameters
    ----------
    Components_Ctypes : Components_Ctypes
        The Components_Ctypes class structured equivalently to a `components_t` struct
        
    Returns
    -------
    Source_Ctypes
        The Source_Ctypes class structured equivalent to a `source_t` struct
    """
    
    class Source_Ctypes(Structure):
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
        :cvar Components_Ctypes_Double point_components: `Components_Ctypes` holding component information for all POINT COMPONENTs in this SOURCE
        :cvar Components_Ctypes_Double gauss_components: `Components_Ctypes` holding component information for all GAUSSIAN COMPONENTs in this SOURCE
        :cvar Components_Ctypes_Double shape_components: `Components_Ctypes` holding component information for all SHAPELET COMPONENTs in this SOURCE
        
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
            ("point_components",  Components_Ctypes),
            ("gauss_components",  Components_Ctypes),
            ("shape_components",  Components_Ctypes),
        ]
        
    return Source_Ctypes

class Source_Python(Structure):
    """A class structured equivalent to a `source_t` struct,  but using
    Python native data types. This can be pickled and therefore fed to Python
    multiprocessing functions. Needed as `ctypes` pointers cannot be pickled,
    and therefore cannot be passed to Python multiprocessing functions.
    """
        
    def __init__(self) -> None:
        self.name = None
        self.n_comps = None
        self.n_points = None
        self.n_point_lists = None
        self.n_point_powers = None
        self.n_point_curves = None
        self.n_gauss = None
        self.n_gauss_lists = None
        self.n_gauss_powers = None
        self.n_gauss_curves = None
        self.n_shapes = None
        self.n_shape_lists = None
        self.n_shape_powers = None
        self.n_shape_curves = None
        self.n_shape_coeffs = None
        self.point_components = Components_Python()
        self.gauss_components = Components_Python()
        self.shape_components = Components_Python()
        

##This call is so we can use it as a type annotation, and so sphinx can document the class
Source_Ctypes = create_source_struct(Components_Ctypes)
    
def create_source_catalogue_struct(Source_Ctypes : Source_Ctypes): # type: ignore
    """Creates a `Source_Catalogue` class structured equivalent to a 
    `source_catalogue_t` struct, used by the C and CUDA code. The `Source_Ctypes` class
    is passed in as the `Source_Ctypes` class is created dynamically depending on
    the required precision.
    
    Parameters
    ----------
    Source_Ctypes : Source_Ctypes
        The Source_Ctypes class structured equivalent to a `source_t` struct
        
    Returns
    -------
    Source_Catalogue
        The Source_Catalogue class structured equivalent to a `source_catalogue_t` struct
    """
    
    class Source_Catalogue(Structure):
        """
        A class structured equivalent to a `source_t` struct, used by 
        the C/GPU code
        
        Attributes
        -----------
        num_sources : int
            The number of sources in the catalogue
        num_shapelets : int
            The total number of shapelets components in the catalogue
        sources : POINTER(Source_Ctypes)
            A pointer to an array of `Source_Float` objects representing the sources in the catalogue
        """
        
        _fields_ = [("num_sources", c_int),
                    ("num_shapelets", c_int),
                    ("sources", POINTER(Source_Ctypes))]
        
    return Source_Catalogue

##This call is so we can use it as a type annotation, and so sphinx can document the class
Source_Catalogue = create_source_catalogue_struct(Source_Ctypes)
    
    
def setup_source_catalogue(Source_Ctypes : Source_Ctypes, Source_Catalogue : Source_Catalogue, # type: ignore
                           num_sources : int, num_shapelets : int,
                           precision = "double") -> Source_Catalogue: # type: ignore
    """
    Creates a `Source_Catalogue` with the specified number of sources and shapelets.
    Sets source_catalogue.sources an array of `Source_Ctypes` objects of length `num_sources`
    
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
    source_catalogue : Source_Catalogue
        The initialised source_catalogue object.
    """
    
    source_catalogue = Source_Catalogue()
    
    source_catalogue.sources = (num_sources*Source_Ctypes)()
    source_catalogue.num_sources = num_sources
    source_catalogue.num_shapelets = num_shapelets
    
    return source_catalogue
    
def setup_components(chunk_map: Skymodel_Chunk_Map,
                     chunked_source: Source_Python, # type: ignore
                     num_freqs: int, num_times: int, num_beams: int,
                     comp_type: CompTypes, beamtype: int,
                     precision: str = "double") -> None:
    """
    Given the mapping information in `chunk_map`, initialise the necessary
    components in `chunked_source` for the given `comp_type`, this being one
    of
     - chunked_source.point_components
     - chunked_source.gaussion_components
     - chunked_source.shape_components
    This setup includes allocating numpy arrays of specific precision for
    the components, depending on what type of `chunked_source` is passed in.

    Parameters
    ------------
    chunk_map: Skymodel_Chunk_Map
        Object containing information about the chunk of the sky model.
    chunked_source: Union[Source_Ctypes, Source_Python]
        Object containing the chunked source information.
    num_freqs: int
        Number of frequencies
    num_times: int
        Number of times
    num_beams : int
        Number of beams to be calculated (only used for certain beam types).
        If using one beam for all station/tiles, num_beams=1, otherwise
        num_beams=num_ants.
    comp_type: CompTypes
        enum representing the type of component.
    beamtype: int
        representing the type of beam.
    precision: str
       Either "float" or "double", default "double".
    """
    
    if comp_type == CompTypes.POINT:
        components = chunked_source.point_components
        n_comps = chunk_map.n_points
        map_components = chunk_map.point_components
        n_powers = chunk_map.n_point_powers
        n_curves = chunk_map.n_point_curves
        n_lists = chunk_map.n_point_lists
        total_num_flux_entires = chunk_map.point_components.total_num_flux_entires
        
    elif comp_type == CompTypes.GAUSSIAN:
        components = chunked_source.gauss_components
        n_comps = chunk_map.n_gauss
        map_components = chunk_map.gauss_components
        n_powers = chunk_map.n_gauss_powers
        n_curves = chunk_map.n_gauss_curves
        n_lists = chunk_map.n_gauss_lists
        total_num_flux_entires = chunk_map.gauss_components.total_num_flux_entires
        
    # elif comp_type == CompTypes.SHAPELET:
    else:
        components = chunked_source.shape_components
        n_comps = chunk_map.n_shapes
        map_components = chunk_map.shape_components
        n_powers = chunk_map.n_shape_powers
        n_curves = chunk_map.n_shape_curves
        n_lists = chunk_map.n_shape_lists
        total_num_flux_entires = chunk_map.shape_components.total_num_flux_entires
    
    ##Assign empty numpy arrays to everything. As these aren't being
    ##sent to the GPU (they will be copied into a ctypes array later),
    ##Make sure we match the float/double precision that will be
    ##going into the ctypes equivalent arrays
        
    if precision == 'float':
        np_precision = np.float32
    else:
        np_precision = np.float64
    
    components.ras = np.empty(n_comps, dtype=np.float64)
    components.decs = np.empty(n_comps, dtype=np.float64)
    
    ##power-law flux things
    components.power_ref_freqs = np.empty(n_powers, dtype=np.float64)
    components.power_ref_stokesI = np.empty(n_powers, dtype=np_precision)
    components.power_SIs = np.empty(n_powers, dtype=np_precision)
    
    ##curved power-law flux things
    components.curve_ref_freqs = np.empty(n_curves, dtype=np.float64)
    components.curve_ref_stokesI = np.empty(n_curves, dtype=np_precision)
    components.curve_SIs = np.empty(n_curves, dtype=np_precision)
    components.curve_qs = np.empty(n_curves, dtype=np_precision)
    
    ##list flux things
    components.list_freqs = np.empty(total_num_flux_entires, dtype=np.float64)
    components.list_stokesI = np.empty(total_num_flux_entires, dtype=np_precision)
    components.list_stokesQ = np.empty(total_num_flux_entires, dtype=np_precision)
    components.list_stokesU = np.empty(total_num_flux_entires, dtype=np_precision)
    components.list_stokesV = np.empty(total_num_flux_entires, dtype=np_precision)
    
    components.num_list_values = np.empty(n_lists, dtype=np.int32)
    components.list_start_indexes = np.empty(n_lists, dtype=np.int32)
    
    ##indexes of types
    components.power_comp_inds = np.empty(n_powers, dtype=np.int32)
    components.curve_comp_inds = np.empty(n_curves, dtype=np.int32)
    components.list_comp_inds = np.empty(n_lists, dtype=np.int32)
    
    if comp_type == CompTypes.GAUSSIAN or comp_type == CompTypes.SHAPELET:
        components.majors = np.empty(n_comps, dtype=np_precision)
        components.minors = np.empty(n_comps, dtype=np_precision)
        components.pas = np.empty(n_comps, dtype=np_precision)
        
        
    if comp_type == CompTypes.SHAPELET:
        total_shape_coeffs = chunk_map.shape_components.total_shape_coeffs
        components.shape_coeffs = np.empty(total_shape_coeffs, dtype=np_precision)
        components.n1s = np.empty(total_shape_coeffs, dtype=np_precision)
        components.n2s = np.empty(total_shape_coeffs, dtype=np_precision)
        components.param_indexes = np.empty(total_shape_coeffs, dtype=np_precision)
        
    ##----------------------------------------------------------------
    ##now we make space for coordinates that are need for the primary beam
    if beamtype in BeamGroups.hadec_beam_values:
        components.beam_has = np.empty(n_comps*num_times, dtype=np.float64)
        components.beam_decs = np.empty(n_comps*num_times, dtype=np.float64)
        
    if beamtype in BeamGroups.azza_beam_values:
        components.azs = np.empty(n_comps*num_times, dtype=np_precision)
        components.zas = np.empty(n_comps*num_times, dtype=np_precision)
        
    if beamtype == BeamTypes.EB_MWA.value:
        components.para_angles = np.empty(n_comps*num_times, dtype=np_precision)
        
    ##NOTE: use this to calculate beam values in Python
    ##Needs to be done in future for pyuvbeam
    # if beamtype in BeamGroups.eb_beam_values:
    #     components.gxs = np.empty(n_comps*num_beams*num_freqs*num_times, dtype=np.complex128)
    #     components.Dxs = np.empty(n_comps*num_beams*num_freqs*num_times, dtype=np.complex128)
    #     components.Dys = np.empty(n_comps*num_beams*num_freqs*num_times, dtype=np.complex128)
    #     components.gys = np.empty(n_comps*num_beams*num_freqs*num_times, dtype=np.complex128)
        
    ##now do the polarisation thingies------------------------------------------
    
    components.stokesV_pol_fracs = np.empty(map_components.num_v_pol_fracs, dtype=np_precision)
    components.stokesV_pol_frac_comp_inds = np.empty(map_components.num_v_pol_fracs, dtype=np.int32)
    components.stokesV_power_ref_flux = np.empty(map_components.num_v_powers, dtype=np_precision)
    components.stokesV_power_SIs = np.empty(map_components.num_v_powers, dtype=np_precision)
    components.stokesV_power_comp_inds = np.empty(map_components.num_v_powers, dtype=np.int32)
    components.stokesV_curve_ref_flux = np.empty(map_components.num_v_curves, dtype=np_precision)
    components.stokesV_curve_SIs = np.empty(map_components.num_v_curves, dtype=np_precision)
    components.stokesV_curve_qs = np.empty(map_components.num_v_curves, dtype=np_precision)
    components.stokesV_curve_comp_inds = np.empty(map_components.num_v_curves, dtype=np.int32)
    components.linpol_pol_fracs = np.empty(map_components.num_lin_pol_fracs, dtype=np_precision)
    components.linpol_pol_frac_comp_inds = np.empty(map_components.num_lin_pol_fracs, dtype=np.int32)
    components.linpol_power_ref_flux =np.empty(map_components.num_lin_powers, dtype=np_precision)
    components.linpol_power_SIs = np.empty(map_components.num_lin_powers, dtype=np_precision)
    components.linpol_power_comp_inds = np.empty(map_components.num_lin_powers, dtype=np.int32)
    components.linpol_curve_ref_flux = np.empty(map_components.num_lin_curves, dtype=np_precision)
    components.linpol_curve_SIs = np.empty(map_components.num_lin_curves, dtype=np_precision)
    components.linpol_curve_qs = np.empty(map_components.num_lin_curves, dtype=np_precision)
    components.linpol_curve_comp_inds = np.empty(map_components.num_lin_curves, dtype=np.int32)
    components.rm_values = np.empty(map_components.num_lin_angles, dtype=np_precision)
    components.intr_pol_angle = np.empty(map_components.num_lin_angles, dtype=np_precision)
    components.linpol_angle_inds = np.empty(map_components.num_lin_angles, dtype=np.int32)
    
    components.stokesV_list_ref_freqs = np.empty(map_components.total_num_v_flux_entires, dtype=np_precision)
    components.stokesV_list_ref_flux = np.empty(map_components.total_num_v_flux_entires, dtype=np_precision)
    components.stokesV_list_comp_inds = np.empty(map_components.num_v_lists, dtype=np.int32)
    components.stokesV_num_list_values = np.empty(map_components.num_v_lists, dtype=np.int32)
    components.stokesV_list_start_indexes = np.empty(map_components.num_v_lists, dtype=np.int32)
    
    components.stokesQ_list_ref_freqs = np.empty(map_components.total_num_q_flux_entires, dtype=np_precision)
    components.stokesQ_list_ref_flux = np.empty(map_components.total_num_q_flux_entires, dtype=np_precision)
    components.stokesQ_list_comp_inds = np.empty(map_components.num_lin_lists, dtype=np.int32)
    components.stokesQ_num_list_values = np.empty(map_components.num_lin_lists, dtype=np.int32)
    components.stokesQ_list_start_indexes = np.empty(map_components.num_lin_lists, dtype=np.int32)
    
    components.stokesU_list_ref_freqs = np.empty(map_components.total_num_u_flux_entires, dtype=np_precision)
    components.stokesU_list_ref_flux = np.empty(map_components.total_num_u_flux_entires, dtype=np_precision)
    components.stokesU_list_comp_inds = np.empty(map_components.num_lin_lists, dtype=np.int32)
    components.stokesU_num_list_values = np.empty(map_components.num_lin_lists, dtype=np.int32)
    components.stokesU_list_start_indexes = np.empty(map_components.num_lin_lists, dtype=np.int32)
    
    components.linpol_p_list_ref_freqs = np.empty(map_components.total_num_lin_p_flux_entires, dtype=np_precision)
    components.linpol_p_list_ref_flux = np.empty(map_components.total_num_lin_p_flux_entires, dtype=np_precision)
    components.linpol_p_list_comp_inds = np.empty(map_components.num_lin_p_lists, dtype=np.int32)
    components.linpol_p_num_list_values = np.empty(map_components.num_lin_p_lists, dtype=np.int32)
    components.linpol_p_list_start_indexes = np.empty(map_components.num_lin_p_lists, dtype=np.int32)
    
    
    num_primarybeam_values = n_comps*num_freqs*num_times
    components.num_primarybeam_values = num_primarybeam_values
    components.total_num_flux_entires = total_num_flux_entires
    
    components.n_stokesV_pol_frac = map_components.num_v_pol_fracs
    components.n_stokesV_power = map_components.num_v_powers
    components.n_stokesV_curve = map_components.num_v_curves
    components.n_linpol_pol_frac = map_components.num_lin_pol_fracs
    components.n_stokesV_curve = map_components.num_v_curves
    components.n_linpol_pol_frac = map_components.num_lin_pol_fracs
    components.n_linpol_power = map_components.num_lin_powers
    components.n_linpol_curve = map_components.num_lin_curves
    components.n_linpol_angles = map_components.num_lin_angles
    
    components.n_stokesV_list = map_components.num_v_lists
    components.n_linpol_list = map_components.num_lin_lists
    components.n_linpol_p_list = map_components.num_lin_p_lists
    
    components.n_stokesV_list_flux_entries = map_components.total_num_v_flux_entires
    components.n_stokesQ_list_flux_entries = map_components.total_num_q_flux_entires
    components.n_stokesU_list_flux_entries = map_components.total_num_u_flux_entires
    components.n_linpol_p_list_flux_entries = map_components.total_num_lin_p_flux_entires
    
        
    
def setup_chunked_source(chunked_source : Source_Python,
                         chunk_map : Skymodel_Chunk_Map, 
                         num_freqs : int, num_times : int, 
                         num_beams : int, beamtype : int,
                         precision='double') -> Source_Python:
    """
    Sets up `chunked_source`, which will be entually be copied to a ctypes
    struct and then passed to the C/GPU functions.  Allocates the
    correct amount of memory, based on whether the precision is either
    'double' or 'float', as well as what is detailed in the `chunk_map`.

    Parameters
    ------------
    chunked_source : Source_Python
        Object to contain the chunked source information.
    chunk_map : Skymodel_Chunk_Map
        Object containing information about the chunk of the sky model.
    num_freqs : int
        The number of frequency channels.
    num_times : int
        The number of time samples.
    num_beams : int
        Number of beams to be calculated (only used for certain beam types).
        If using one beam for all station/tiles, num_beams=1, otherwise
        num_beams=num_ants.
    beamtype : int
        The type of primary beam.
    precision : str, optional
        The precision to use. Can be either "float" or "double. Defaults to "double".

    Returns
    ---------
      chunked_source : Source_Python
          The chunked_source that was passed in, but now with all the necessary
          components allocated.
    """
    
    # set_precision_global(precision)
    
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
                         num_beams, CompTypes.POINT, beamtype, precision)
        
    if chunk_map.n_gauss > 0:
        setup_components(chunk_map, chunked_source, num_freqs, num_times,
                         num_beams, CompTypes.GAUSSIAN, beamtype, precision)
        
    if chunk_map.n_shapes > 0:
        setup_components(chunk_map, chunked_source, num_freqs, num_times,
                         num_beams, CompTypes.SHAPELET, beamtype, precision)
    
    return chunked_source



def copy_python_components_to_ctypes(python_comps: Components_Python,
                                     ctypes_comps: Components_Ctypes,  # type: ignore
                                     comp_type: CompTypes, n_powers: int,
                                     n_curves: int, n_lists: int,
                                     beamtype: int,
                                     precision='double'):
    """
    Copies data from `python_comps` to `ctypes_comps` ctypes structures for use
    in C libraries. Handles memory allocation itself via the
    np.ndarray.ctypes.data_as() method, combined with the precision as
    specified by `precision``
    
    Parameters
    -----------
    python_comps : Components_Python
        The Python components containing the data to be copied.
    ctypes_comps : Components_Ctypes
        The ctypes components where the data will be copied to.
    comp_type : CompTypes
        The type of components being copied (e.g., POINT, GAUSSIAN, SHAPELET).
    n_powers : int
        Number of power components.
    n_curves : int
        Number of curve components.
    n_lists : int
        Number of list components.
    beamtype : int
        The type of beam being used.
    precision : str, optional
        Precision of the data ('float' or 'double', default is 'double').
    """

    if precision == 'float':
        c_user_precision = c_float
        c_user_precision_complex = c_float_complex
    else:
        c_user_precision = c_double
        c_user_precision_complex = c_double_complex
    
    n_comps = n_powers + n_curves + n_lists
    total_num_flux_entires = int(python_comps.total_num_flux_entires)
    
    ctypes_comps.ras = python_comps.ras.ctypes.data_as(POINTER(c_double))
    ctypes_comps.decs = python_comps.decs.ctypes.data_as(POINTER(c_double))
    ctypes_comps.power_ref_freqs = python_comps.power_ref_freqs.ctypes.data_as(POINTER(c_double))
    ctypes_comps.power_ref_stokesI = python_comps.power_ref_stokesI.ctypes.data_as(POINTER(c_user_precision))
    ctypes_comps.power_SIs = python_comps.power_SIs.ctypes.data_as(POINTER(c_user_precision))
    ctypes_comps.curve_ref_freqs = python_comps.curve_ref_freqs.ctypes.data_as(POINTER(c_double))
    ctypes_comps.curve_ref_stokesI = python_comps.curve_ref_stokesI.ctypes.data_as(POINTER(c_user_precision))
    ctypes_comps.curve_SIs = python_comps.curve_SIs.ctypes.data_as(POINTER(c_user_precision))
    ctypes_comps.curve_qs = python_comps.curve_qs.ctypes.data_as(POINTER(c_user_precision))
    ctypes_comps.power_comp_inds = python_comps.power_comp_inds.ctypes.data_as(POINTER(c_int))
    ctypes_comps.curve_comp_inds = python_comps.curve_comp_inds.ctypes.data_as(POINTER(c_int))
    ctypes_comps.list_comp_inds = python_comps.list_comp_inds.ctypes.data_as(POINTER(c_int))
    ctypes_comps.list_freqs = python_comps.list_freqs.ctypes.data_as(POINTER(c_double))
    ctypes_comps.list_stokesI = python_comps.list_stokesI.ctypes.data_as(POINTER(c_user_precision))
    ctypes_comps.list_stokesQ = python_comps.list_stokesQ.ctypes.data_as(POINTER(c_user_precision))
    ctypes_comps.list_stokesU = python_comps.list_stokesU.ctypes.data_as(POINTER(c_user_precision))
    ctypes_comps.list_stokesV = python_comps.list_stokesV.ctypes.data_as(POINTER(c_user_precision))
    
    ctypes_comps.num_list_values = python_comps.num_list_values.ctypes.data_as(POINTER(c_int))
    ctypes_comps.list_start_indexes = python_comps.list_start_indexes.ctypes.data_as(POINTER(c_int))
    ctypes_comps.total_num_flux_entires = python_comps.total_num_flux_entires
    
    if comp_type == CompTypes.GAUSSIAN or comp_type == CompTypes.SHAPELET:
        ctypes_comps.majors = python_comps.majors.ctypes.data_as(POINTER(c_user_precision))
        ctypes_comps.minors = python_comps.minors.ctypes.data_as(POINTER(c_user_precision))
        ctypes_comps.pas = python_comps.pas.ctypes.data_as(POINTER(c_user_precision))
    
    if comp_type == CompTypes.SHAPELET:
        ctypes_comps.shape_coeffs = python_comps.shape_coeffs.ctypes.data_as(POINTER(c_user_precision))
        ctypes_comps.n1s = python_comps.n1s.ctypes.data_as(POINTER(c_user_precision))
        ctypes_comps.n2s = python_comps.n2s.ctypes.data_as(POINTER(c_user_precision))
        ctypes_comps.param_indexes = python_comps.param_indexes.ctypes.data_as(POINTER(c_user_precision))
    
    
    if beamtype in BeamGroups.azza_beam_values:
        ctypes_comps.azs = python_comps.azs.ctypes.data_as(POINTER(c_user_precision))
        ctypes_comps.zas = python_comps.zas.ctypes.data_as(POINTER(c_user_precision))
        
    if beamtype == BeamTypes.EB_MWA.value:
        ctypes_comps.para_angles = python_comps.para_angles.ctypes.data_as(POINTER(c_user_precision))
        # print("WE DID THIS MOTHER HUCKER", python_comps.para_angles, python_comps.azs, python_comps.zas)
        
    if beamtype in BeamGroups.hadec_beam_values:
        ctypes_comps.beam_has = python_comps.beam_has.ctypes.data_as(POINTER(c_double))
        ctypes_comps.beam_decs = python_comps.beam_decs.ctypes.data_as(POINTER(c_double))
    
    ctypes_comps.num_primarybeam_values = python_comps.num_primarybeam_values
    ##things to hold the beam gain

    ##NOTE: this needs to happen when calculating beam values in Python    
    # if beamtype in BeamGroups.eb_beam_values:
    #     num_gains = len(python_comps.gxs)
        
    #     complex_num_beams = c_user_precision_complex*num_gains
            
    #     ctypes_comps.gxs = complex_num_beams()
    #     ctypes_comps.Dxs = complex_num_beams()
    #     ctypes_comps.Dys = complex_num_beams()
    #     ctypes_comps.gys = complex_num_beams()
        
    #     ##Actually iterate over the complex beam gains, as the ctypes complex
    #     ##object is actually a bespoke class
    #     for beam_ind in range(num_gains):
    #         ctypes_comps.gxs[beam_ind].real = python_comps.gxs[beam_ind].real
    #         ctypes_comps.gxs[beam_ind].imag = python_comps.gxs[beam_ind].imag
    #         ctypes_comps.Dxs[beam_ind].real = python_comps.Dxs[beam_ind].real
    #         ctypes_comps.Dxs[beam_ind].imag = python_comps.Dxs[beam_ind].imag
    #         ctypes_comps.Dys[beam_ind].real = python_comps.Dys[beam_ind].real
    #         ctypes_comps.Dys[beam_ind].imag = python_comps.Dys[beam_ind].imag
    #         ctypes_comps.gys[beam_ind].real = python_comps.gys[beam_ind].real
    #         ctypes_comps.gys[beam_ind].imag = python_comps.gys[beam_ind].imag
        
    ctypes_comps.n_stokesV_pol_frac = 0
    ctypes_comps.n_stokesV_power = 0
    ctypes_comps.n_stokesV_curve = 0
    ctypes_comps.n_stokesV_list = 0
    ctypes_comps.n_linpol_pol_frac = 0
    ctypes_comps.n_linpol_power = 0
    ctypes_comps.n_linpol_curve = 0
    ctypes_comps.n_linpol_list = 0
    ctypes_comps.n_linpol_p_list = 0
    ctypes_comps.n_linpol_angles = 0
    ctypes_comps.n_stokesV_list_flux_entries = 0
    ctypes_comps.n_stokesQ_list_flux_entries = 0
    ctypes_comps.n_stokesU_list_flux_entries = 0
    ctypes_comps.n_linpol_p_list_flux_entries = 0
    
    if python_comps.n_stokesV_pol_frac > 0:
        ctypes_comps.n_stokesV_pol_frac = python_comps.n_stokesV_pol_frac
        ctypes_comps.stokesV_pol_fracs = python_comps.stokesV_pol_fracs.ctypes.data_as(POINTER(c_user_precision))
        ctypes_comps.stokesV_pol_frac_comp_inds = python_comps.stokesV_pol_frac_comp_inds.ctypes.data_as(POINTER(c_int))
        
        
    if python_comps.n_stokesV_power > 0:
        ctypes_comps.n_stokesV_power = python_comps.n_stokesV_power
        ctypes_comps.stokesV_power_ref_flux = python_comps.stokesV_power_ref_flux.ctypes.data_as(POINTER(c_user_precision))
        ctypes_comps.stokesV_power_SIs = python_comps.stokesV_power_SIs.ctypes.data_as(POINTER(c_user_precision))
        ctypes_comps.stokesV_power_comp_inds = python_comps.stokesV_power_comp_inds.ctypes.data_as(POINTER(c_int))
        
    if python_comps.n_stokesV_curve > 0:
        ctypes_comps.n_stokesV_curve = python_comps.n_stokesV_curve
        ctypes_comps.stokesV_curve_ref_flux = python_comps.stokesV_curve_ref_flux.ctypes.data_as(POINTER(c_user_precision))
        ctypes_comps.stokesV_curve_SIs = python_comps.stokesV_curve_SIs.ctypes.data_as(POINTER(c_user_precision))
        ctypes_comps.stokesV_curve_qs = python_comps.stokesV_curve_qs.ctypes.data_as(POINTER(c_user_precision))
        ctypes_comps.stokesV_curve_comp_inds = python_comps.stokesV_curve_comp_inds.ctypes.data_as(POINTER(c_int))
        
    if python_comps.n_stokesV_list > 0:
        ctypes_comps.n_stokesV_list = python_comps.n_stokesV_list
        ctypes_comps.n_stokesV_list_flux_entries = python_comps.n_stokesV_list_flux_entries
        ctypes_comps.stokesV_list_ref_freqs = python_comps.stokesV_list_ref_freqs.ctypes.data_as(POINTER(c_double))
        ctypes_comps.stokesV_list_ref_flux = python_comps.stokesV_list_ref_flux.ctypes.data_as(POINTER(c_user_precision))
        ctypes_comps.stokesV_list_comp_inds = python_comps.stokesV_list_comp_inds.ctypes.data_as(POINTER(c_int))
        ctypes_comps.stokesV_num_list_values = python_comps.stokesV_num_list_values.ctypes.data_as(POINTER(c_int))
        ctypes_comps.stokesV_list_start_indexes = python_comps.stokesV_list_start_indexes.ctypes.data_as(POINTER(c_int))
        
    if python_comps.n_linpol_pol_frac > 0:
        ctypes_comps.n_linpol_pol_frac = python_comps.n_linpol_pol_frac
        ctypes_comps.linpol_pol_fracs = python_comps.linpol_pol_fracs.ctypes.data_as(POINTER(c_user_precision))
        ctypes_comps.linpol_pol_frac_comp_inds = python_comps.linpol_pol_frac_comp_inds.ctypes.data_as(POINTER(c_int))
        
        
    if python_comps.n_linpol_power > 0:
        ctypes_comps.n_linpol_power = python_comps.n_linpol_power
        ctypes_comps.linpol_power_ref_flux = python_comps.linpol_power_ref_flux.ctypes.data_as(POINTER(c_user_precision))
        ctypes_comps.linpol_power_SIs = python_comps.linpol_power_SIs.ctypes.data_as(POINTER(c_user_precision))
        ctypes_comps.linpol_power_comp_inds = python_comps.linpol_power_comp_inds.ctypes.data_as(POINTER(c_int))
        
        
    if python_comps.n_linpol_curve > 0:
        ctypes_comps.n_linpol_curve = python_comps.n_linpol_curve
        ctypes_comps.linpol_curve_ref_flux = python_comps.linpol_curve_ref_flux.ctypes.data_as(POINTER(c_user_precision))
        ctypes_comps.linpol_curve_SIs = python_comps.linpol_curve_SIs.ctypes.data_as(POINTER(c_user_precision))
        ctypes_comps.linpol_curve_qs = python_comps.linpol_curve_qs.ctypes.data_as(POINTER(c_user_precision))
        ctypes_comps.linpol_curve_comp_inds = python_comps.linpol_curve_comp_inds.ctypes.data_as(POINTER(c_int))
        
    if python_comps.n_linpol_angles > 0:
        ctypes_comps.n_linpol_angles = python_comps.n_linpol_angles
        ctypes_comps.rm_values = python_comps.rm_values.ctypes.data_as(POINTER(c_user_precision))
        ctypes_comps.intr_pol_angle = python_comps.intr_pol_angle.ctypes.data_as(POINTER(c_user_precision))
        ctypes_comps.linpol_angle_inds = python_comps.linpol_angle_inds.ctypes.data_as(POINTER(c_int))
        
    if python_comps.n_linpol_list > 0:
        ctypes_comps.n_linpol_list = python_comps.n_linpol_list
        
        ctypes_comps.n_stokesQ_list_flux_entries = python_comps.n_stokesQ_list_flux_entries
        ctypes_comps.stokesQ_list_ref_freqs = python_comps.stokesQ_list_ref_freqs.ctypes.data_as(POINTER(c_double))
        ctypes_comps.stokesQ_list_ref_flux = python_comps.stokesQ_list_ref_flux.ctypes.data_as(POINTER(c_user_precision))
        ctypes_comps.stokesQ_list_comp_inds = python_comps.stokesQ_list_comp_inds.ctypes.data_as(POINTER(c_int))
        ctypes_comps.stokesQ_num_list_values = python_comps.stokesQ_num_list_values.ctypes.data_as(POINTER(c_int))
        ctypes_comps.stokesQ_list_start_indexes = python_comps.stokesQ_list_start_indexes.ctypes.data_as(POINTER(c_int))
        
        ctypes_comps.n_stokesU_list_flux_entries = python_comps.n_stokesU_list_flux_entries
        ctypes_comps.stokesU_list_ref_freqs = python_comps.stokesU_list_ref_freqs.ctypes.data_as(POINTER(c_double))
        ctypes_comps.stokesU_list_ref_flux = python_comps.stokesU_list_ref_flux.ctypes.data_as(POINTER(c_user_precision))
        ctypes_comps.stokesU_list_comp_inds = python_comps.stokesU_list_comp_inds.ctypes.data_as(POINTER(c_int))
        ctypes_comps.stokesU_num_list_values = python_comps.stokesU_num_list_values.ctypes.data_as(POINTER(c_int))
        ctypes_comps.stokesU_list_start_indexes = python_comps.stokesU_list_start_indexes.ctypes.data_as(POINTER(c_int))
        
    if python_comps.n_linpol_p_list > 0:
        ctypes_comps.n_linpol_p_list = python_comps.n_linpol_p_list
        ctypes_comps.n_linpol_p_list_flux_entries = python_comps.n_linpol_p_list_flux_entries
        ctypes_comps.linpol_p_list_ref_freqs = python_comps.linpol_p_list_ref_freqs.ctypes.data_as(POINTER(c_double))
        ctypes_comps.linpol_p_list_ref_flux = python_comps.linpol_p_list_ref_flux.ctypes.data_as(POINTER(c_user_precision))
        ctypes_comps.linpol_p_list_comp_inds = python_comps.linpol_p_list_comp_inds.ctypes.data_as(POINTER(c_int))
        ctypes_comps.linpol_p_num_list_values = python_comps.linpol_p_num_list_values.ctypes.data_as(POINTER(c_int))
        ctypes_comps.linpol_p_list_start_indexes = python_comps.linpol_p_list_start_indexes.ctypes.data_as(POINTER(c_int))
        
    ctypes_comps.do_QUV = python_comps.do_QUV
    


def copy_python_source_to_ctypes(python_source: Source_Python,
                                 ctypes_source: Source_Ctypes, # type: ignore
                                 beamtype: int, precision: str = 'double'):
    """
    Copies data from a Python source object to a ctypes source object, which
    can be fed to the C/GPU code. Handles memory allocation itself via the
    np.ndarray.ctypes.data_as() method, combined with the precision as
    specified by `precision``
    
    Parameters
    -----------
    python_source : Source_Python
        The source object containing data in Python format.
    ctypes_source : Source_Ctypes
        The source object to be populated with data in ctypes format.
    beamtype : int
        The type of beam to be used.
    precision : str, optional
        The precision of the data to be copied, either 'float' or 'double'; 
        default is 'double'.
    Notes
    ------
    This function transfers various attributes from the Python source object
    to the ctypes source object, including the number of points, Gaussian
    components, and shapelet components. It also calls `copy_python_components_to_ctypes`
    to copy the components of each type if they exist in the Python source object.
    """
    
    ctypes_source.n_points = python_source.n_points
    ctypes_source.n_point_lists = python_source.n_point_lists
    ctypes_source.n_point_powers = python_source.n_point_powers
    ctypes_source.n_point_curves = python_source.n_point_curves
    ctypes_source.n_gauss = python_source.n_gauss
    ctypes_source.n_gauss_lists = python_source.n_gauss_lists
    ctypes_source.n_gauss_powers = python_source.n_gauss_powers
    ctypes_source.n_gauss_curves = python_source.n_gauss_curves
    ctypes_source.n_shapes = python_source.n_shapes
    ctypes_source.n_shape_lists = python_source.n_shape_lists
    ctypes_source.n_shape_powers = python_source.n_shape_powers
    ctypes_source.n_shape_curves = python_source.n_shape_curves
    ctypes_source.n_shape_coeffs = python_source.n_shape_coeffs
    ctypes_source.n_comps = python_source.n_comps
    
    if python_source.n_points:
        copy_python_components_to_ctypes(python_source.point_components,
                                         ctypes_source.point_components,
                                         CompTypes.POINT,
                                         python_source.n_point_powers,
                                         python_source.n_point_curves,
                                         python_source.n_point_lists,
                                         beamtype,
                                         precision=precision)
        
    if python_source.n_gauss:
        copy_python_components_to_ctypes(python_source.gauss_components,
                                         ctypes_source.gauss_components,
                                         CompTypes.GAUSSIAN,
                                         ctypes_source.n_gauss_powers,
                                         ctypes_source.n_gauss_curves,
                                         ctypes_source.n_gauss_lists,
                                         beamtype,
                                         precision=precision)
        
    if python_source.n_shapes:
        copy_python_components_to_ctypes(python_source.shape_components,
                                         ctypes_source.shape_components,
                                         CompTypes.SHAPELET,
                                         ctypes_source.n_shape_powers,
                                         ctypes_source.n_shape_curves,
                                         ctypes_source.n_shape_lists,
                                         beamtype,
                                         precision=precision)
        
##TODO we now already have a Components_Python class, so we should use that
##here. This was only written for testing so todo not a high priority
class _Components_Python(object):
    """python equivalent to Components_Ctypes_Float or Components_Ctypes_Double"""
    
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
        self.n_stokesV_list = 0
        self.n_linpol_pol_frac = 0
        self.n_linpol_power = 0
        self.n_linpol_curve = 0
        self.n_linpol_list = 0
        self.n_linpol_p_list = 0
        self.n_linpol_angles = 0
        self.n_stokesV_list_flux_entries = 0
        self.n_stokesQ_list_flux_entries = 0
        self.n_stokesU_list_flux_entries = 0
        self.n_linpol_p_list_flux_entries = 0
        
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
            
        if components.n_stokesV_list > 0:
            self.n_stokesV_list = components.n_stokesV_list
            self.n_stokesV_list_flux_entries = components.n_stokesV_list_flux_entries
            self.stokesV_list_ref_freqs = np.ctypeslib.as_array(components.stokesV_list_ref_freqs, shape=(components.n_stokesV_list_flux_entries, ))
            self.stokesV_list_ref_flux = np.ctypeslib.as_array(components.stokesV_list_ref_flux, shape=(components.n_stokesV_list_flux_entries, ))
            self.stokesV_list_comp_inds = np.ctypeslib.as_array(components.stokesV_list_comp_inds, shape=(components.n_stokesV_list, ))
            self.stokesV_num_list_values = np.ctypeslib.as_array(components.stokesV_num_list_values, shape=(components.n_stokesV_list, ))
            self.stokesV_list_start_indexes = np.ctypeslib.as_array(components.stokesV_list_start_indexes, shape=(components.n_stokesV_list, ))
            
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
            
        if components.n_linpol_list > 0:
            self.n_linpol_list = components.n_linpol_list
            
            self.n_stokesQ_list_flux_entries = components.n_stokesQ_list_flux_entries
            self.stokesQ_list_ref_freqs = np.ctypeslib.as_array(components.stokesQ_list_ref_freqs, shape=(components.n_stokesQ_list_flux_entries, ))
            self.stokesQ_list_ref_flux = np.ctypeslib.as_array(components.stokesQ_list_ref_flux, shape=(components.n_stokesQ_list_flux_entries, ))
            self.stokesQ_list_comp_inds = np.ctypeslib.as_array(components.stokesQ_list_comp_inds, shape=(components.n_linpol_list, ))
            self.stokesQ_num_list_values = np.ctypeslib.as_array(components.stokesQ_num_list_values, shape=(components.n_linpol_list, ))
            self.stokesQ_list_start_indexes = np.ctypeslib.as_array(components.stokesQ_list_start_indexes, shape=(components.n_linpol_list, ))
            
            self.n_stokesU_list_flux_entries = components.n_stokesU_list_flux_entries
            self.stokesU_list_ref_freqs = np.ctypeslib.as_array(components.stokesU_list_ref_freqs, shape=(components.n_stokesU_list_flux_entries, ))
            self.stokesU_list_ref_flux = np.ctypeslib.as_array(components.stokesU_list_ref_flux, shape=(components.n_stokesU_list_flux_entries, ))
            self.stokesU_list_comp_inds = np.ctypeslib.as_array(components.stokesU_list_comp_inds, shape=(components.n_linpol_list, ))
            self.stokesU_num_list_values = np.ctypeslib.as_array(components.stokesU_num_list_values, shape=(components.n_linpol_list, ))
            self.stokesU_list_start_indexes = np.ctypeslib.as_array(components.stokesU_list_start_indexes, shape=(components.n_linpol_list, ))
            
        if components.n_linpol_p_list > 0:
            self.n_linpol_p_list = components.n_linpol_p_list
            self.n_linpol_p_list_flux_entries = components.n_linpol_p_list_flux_entries
            
            self.linpol_p_list_ref_freqs = np.ctypeslib.as_array(components.linpol_p_list_ref_freqs, shape=(components.n_linpol_p_list_flux_entries, ))
            self.linpol_p_list_ref_flux = np.ctypeslib.as_array(components.linpol_p_list_ref_flux, shape=(components.n_linpol_p_list_flux_entries, ))
            self.linpol_p_list_comp_inds = np.ctypeslib.as_array(components.linpol_p_list_comp_inds, shape=(components.n_linpol_p_list, ))
            self.linpol_p_num_list_values = np.ctypeslib.as_array(components.linpol_p_num_list_values, shape=(components.n_linpol_p_list, ))
            self.linpol_p_list_start_indexes = np.ctypeslib.as_array(components.linpol_p_list_start_indexes, shape=(components.n_linpol_p_list, ))
            
        self.do_QUV = components.do_QUV
            
##TODO we now already have a Source_Python class, so we should use that
##here. This was only written for testing so todo not a high priority
class _Ctype_Source_Into_Python(object):
    """
    Class to convert a ctype Source_Ctypes model into a pythonic version
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

 