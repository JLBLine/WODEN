from sys import path
import os
import unittest
import numpy as np
import numpy.testing as npt

# ##Code we are testing
from wodenpy.skymodel import read_yaml_skymodel
# import wodenpy
from wodenpy.skymodel.woden_skymodel import Component_Type_Counter, CompTypes
from wodenpy.skymodel.chunk_sky_model import map_chunk_pointgauss, Skymodel_Chunk_Map
# from read_skymodel_common import Skymodel_Settings

D2R = np.pi/180.0

##for now, WODEN has three flux types: power law, curved power law, and list
NUM_FLUX_TYPES = 3

##limits for "all sky" sky model
LOW_DEC = -90.0*D2R
HIGH_DEC = 30.0*D2R

MWA_LATITUDE = -26.7*D2R

NUM_STOKESV_MODELS = 4
NUM_LINPOL_MODELS = 5

class Args:
    """reduced argparse.Namespace class for testing everybeams"""
    def __init__(self):
        self.date = "2000-01-01T12:00:00"
        self.beam_ms_path = "lba.MS"
        self.station_id = 0
        
        self.ra0 = 0
        self.dec0 = 0
        self.num_time_steps = 0
        self.band_nums = [1]
        self.coarse_band_width = 1.28e+6
        self.lowest_channel_freq = 167.035e+6
        self.num_freq_channels = 0
        self.freq_res = 80e+3
        self.hdf5_beam_path = ""
        self.fast_everybeam = False

class Skymodel_Settings:
    """Something to hold all the various settings and pass around between
    functions"""
    def __init__(self, deg_between_comps : float = 0,
                 num_coeff_per_shape : int = 0,
                 num_list_values : int = 0,
                 comps_per_source : int = 0,
                 stokesV_cadence : int = 0,
                 stokesV_frac_cadence : int = 0,
                 stokesV_pl_cadence : int = 0,
                 stokesV_cpl_cadence : int = 0,
                 stokesV_list_cadence : int = 0,
                 stokesV_num_list : int = 0,
                 linpol_cadence : int = 0,
                 linpol_frac_cadence : int = 0,
                 linpol_pl_cadence : int = 0,
                 linpol_cpl_cadence : int = 0,
                 linpol_list_cadence : int = 0,
                 linpol_p_list_cadence : int = 0,
                 linpol_num_list : int = 0,
                 linpol_num_p_list : int = 0,
                 max_num_visibilities : int = 0,
                 num_points : int = 0,
                 num_gauss : int = 0,
                 num_shapes : int = 0,
                 num_time_steps : int = 0,
                 num_baselines : int = 0,
                 num_freqs : int= 0):
        
        self.deg_between_comps = deg_between_comps
        self.num_coeff_per_shape = num_coeff_per_shape
        self.num_list_values = num_list_values
        self.comps_per_source = comps_per_source
        self.stokesV_frac_cadence = stokesV_frac_cadence
        self.stokesV_pl_cadence = stokesV_pl_cadence
        self.stokesV_cpl_cadence = stokesV_cpl_cadence
        self.linpol_frac_cadence = linpol_frac_cadence
        self.linpol_pl_cadence = linpol_pl_cadence
        self.linpol_cpl_cadence = linpol_cpl_cadence
        self.stokesV_list_cadence = stokesV_list_cadence
        self.stokesV_num_list = stokesV_num_list
        self.linpol_list_cadence = linpol_list_cadence
        self.linpol_p_list_cadence = linpol_p_list_cadence
        self.linpol_num_list = linpol_num_list
        self.linpol_num_p_list = linpol_num_p_list
        self.stokesV_cadence = stokesV_cadence
        self.linpol_cadence = linpol_cadence
        self.max_num_visibilities = max_num_visibilities
        self.num_points = num_points
        self.num_gauss = num_gauss
        self.num_time_steps = num_time_steps
        self.num_baselines = num_baselines
        self.num_freqs = num_freqs
        self.num_shapes = num_shapes
        
        self.before_crop_num_coords = 0

def fill_comp_counter_for_chunking(settings) -> Component_Type_Counter:
    """Fill up a Component_Type_Counter based on given numbers. Fill
    all flux model types for each given point, gaussian, and shapelet models,
    e.g. if you set num_points = 5, you get 5 power law points, 5 curved points,
    5 list point types.

    Parameters
    ----------
    num_points : int
        Number of POINT type components of each flux type to add
    num_gauss : int
        Number of GAUSSIAN type components of each flux type to add
    num_shapes : int
        Number of SHAPELET type components of each flux type to add
    num_coeff_per_shape : int
        How many coefficients per shapelet there should be
    num_list_values : int
        How many list entires for the list-type flux components
    num_time_steps : int
        How many time steps in the simulation

    Returns
    -------
    Component_Type_Counter
        A filled Component_Type_Counter class
    """
    
    num_points = settings.num_points
    num_gauss = settings.num_gauss
    num_shapes = settings.num_shapes
    num_coeff_per_shape = settings.num_coeff_per_shape
    num_list_values = settings.num_list_values
    stokesV_cadence = settings.stokesV_cadence
    linpol_cadence = settings.linpol_cadence

    full_comp_counter = Component_Type_Counter()

    ##Overall total number of components
    full_comp_counter.total_comps = NUM_FLUX_TYPES*(num_points + num_gauss + num_shapes);
    
    ##POINT related number
    full_comp_counter.total_point_comps = NUM_FLUX_TYPES*num_points;
    full_comp_counter.num_point_flux_powers = num_points;
    full_comp_counter.num_point_flux_curves = num_points;
    full_comp_counter.num_point_flux_lists = num_points;
    ##GAUSSIAN related number
    full_comp_counter.total_gauss_comps = NUM_FLUX_TYPES*num_gauss;
    full_comp_counter.num_gauss_flux_powers = num_gauss;
    full_comp_counter.num_gauss_flux_curves = num_gauss;
    full_comp_counter.num_gauss_flux_lists = num_gauss;
    ##SHAPELET related number
    full_comp_counter.total_shape_comps = NUM_FLUX_TYPES*num_shapes;
    full_comp_counter.num_shape_flux_powers = num_shapes;
    full_comp_counter.num_shape_flux_curves = num_shapes;
    full_comp_counter.num_shape_flux_lists = num_shapes;
    full_comp_counter.total_shape_basis = NUM_FLUX_TYPES*num_shapes*num_coeff_per_shape;
    
    ##Make the RA/Dec, file_line_nums values equal to index
    total_comps = full_comp_counter.total_comps
    full_comp_counter.comp_ras = np.arange(total_comps)
    full_comp_counter.comp_decs = np.arange(total_comps)
    full_comp_counter.file_line_nums = np.arange(total_comps)
    
    full_comp_counter.comp_types = np.zeros(total_comps, dtype=int)
    full_comp_counter.num_list_fluxes = np.zeros(total_comps, dtype=int)
    full_comp_counter.num_shape_coeffs = np.zeros(total_comps, dtype=int)
    
    ##Don't need source_indexes in testing, as source index used in horizon
    ##cropping, which should have been done before chunking
    # full_comp_counter.source_indexes = np.zeros(total_comps, dtype=int)
    
    ##Put the requested number of components into the counting arrays
    ##Do it in order of
    # POINT_POWER
    # POINT_CURVE
    # POINT_LIST
    # GAUSS_POWER
    # GAUSS_CURVE
    # GAUSS_LIST
    # SHAPE_POWER
    # SHAPE_CURVE
    # SHAPE_LIST
    
    ##POINT related things
    full_comp_counter.comp_types[:num_points] = CompTypes.POINT_POWER.value
    full_comp_counter.comp_types[num_points:2*num_points] = CompTypes.POINT_CURVE.value
    full_comp_counter.comp_types[2*num_points:3*num_points] = CompTypes.POINT_LIST.value
    full_comp_counter.num_list_fluxes[2*num_points:3*num_points] = num_list_values
    
    ##GAUSS related things - go after point sources (if there any)
    offset = NUM_FLUX_TYPES*num_points
    full_comp_counter.comp_types[offset:offset+num_gauss] = CompTypes.GAUSS_POWER.value
    full_comp_counter.comp_types[offset+num_gauss:offset+2*num_gauss] = CompTypes.GAUSS_CURVE.value
    full_comp_counter.comp_types[offset+2*num_gauss:offset+3*num_gauss] = CompTypes.GAUSS_LIST.value
    full_comp_counter.num_list_fluxes[offset+2*num_gauss:offset+3*num_gauss] = num_list_values
    
    ##SHAPE related things go after point and gaussian sources (if there any)
    offset = NUM_FLUX_TYPES*num_points + NUM_FLUX_TYPES*num_gauss
    full_comp_counter.comp_types[offset:offset+num_shapes] = CompTypes.SHAPE_POWER.value
    full_comp_counter.comp_types[offset+num_shapes:offset+2*num_shapes] = CompTypes.SHAPE_CURVE.value
    full_comp_counter.comp_types[offset+2*num_shapes:offset+3*num_shapes] = CompTypes.SHAPE_LIST.value
    full_comp_counter.num_list_fluxes[offset+2*num_shapes:offset+3*num_shapes] = num_list_values
    full_comp_counter.num_shape_coeffs[offset:offset+3*num_shapes] = num_coeff_per_shape

    ##OK, do ALL the types of polarisation models at the same time. Keep count
    ##of how many instances we've added, so we can keep swapping between types
    ##Each pol model type needs to be linked to point,gaussian,shapelet for
    ##chunking mapping purposes
    
    num_v = 0
    if stokesV_cadence:
        full_comp_counter.v_comp_types = np.full(total_comps, np.nan)
        full_comp_counter.num_v_list_fluxes = np.full(total_comps, np.nan)
        for comp in range(total_comps):
            ##We've hit the cadence number, so add a component
            if comp % stokesV_cadence == 0:
                comp_type = full_comp_counter.comp_types[comp]
                
                ##Stick in a polarisation fraction component
                if num_v % NUM_STOKESV_MODELS == 0:
                    if comp_type == CompTypes.POINT_POWER.value or comp_type == CompTypes.POINT_CURVE.value or comp_type == CompTypes.POINT_LIST.value:
                        full_comp_counter.v_comp_types[comp] = CompTypes.V_POINT_POL_FRAC.value
                    elif comp_type == CompTypes.GAUSS_POWER.value or comp_type == CompTypes.GAUSS_CURVE.value or comp_type == CompTypes.GAUSS_LIST.value:
                        full_comp_counter.v_comp_types[comp] = CompTypes.V_GAUSS_POL_FRAC.value
                    else: 
                        full_comp_counter.v_comp_types[comp] = CompTypes.V_SHAPE_POL_FRAC.value
                ##Stick in a power-law component
                elif num_v % NUM_STOKESV_MODELS == 1:
                    if comp_type == CompTypes.POINT_POWER.value or comp_type == CompTypes.POINT_CURVE.value or comp_type == CompTypes.POINT_LIST.value:
                        full_comp_counter.v_comp_types[comp] = CompTypes.V_POINT_POWER.value
                    elif comp_type == CompTypes.GAUSS_POWER.value or comp_type == CompTypes.GAUSS_CURVE.value or comp_type == CompTypes.GAUSS_LIST.value:
                        full_comp_counter.v_comp_types[comp] = CompTypes.V_GAUSS_POWER.value
                    else: 
                        full_comp_counter.v_comp_types[comp] = CompTypes.V_SHAPE_POWER.value
                ##Stick in a curved power-law component
                elif num_v % NUM_STOKESV_MODELS == 2:
                    if comp_type == CompTypes.POINT_POWER.value or comp_type == CompTypes.POINT_CURVE.value or comp_type == CompTypes.POINT_LIST.value:
                        full_comp_counter.v_comp_types[comp] = CompTypes.V_POINT_CURVE.value
                    elif comp_type == CompTypes.GAUSS_POWER.value or comp_type == CompTypes.GAUSS_CURVE.value or comp_type == CompTypes.GAUSS_LIST.value:
                        full_comp_counter.v_comp_types[comp] = CompTypes.V_GAUSS_CURVE.value
                    else: 
                        full_comp_counter.v_comp_types[comp] = CompTypes.V_SHAPE_CURVE.value
                else:
                    full_comp_counter.num_v_list_fluxes[comp] = settings.stokesV_num_list
                    if comp_type == CompTypes.POINT_POWER.value or comp_type == CompTypes.POINT_CURVE.value or comp_type == CompTypes.POINT_LIST.value:
                        full_comp_counter.v_comp_types[comp] = CompTypes.V_POINT_LIST.value
                    elif comp_type == CompTypes.GAUSS_POWER.value or comp_type == CompTypes.GAUSS_CURVE.value or comp_type == CompTypes.GAUSS_LIST.value:
                        full_comp_counter.v_comp_types[comp] = CompTypes.V_GAUSS_LIST.value
                    else: 
                        full_comp_counter.v_comp_types[comp] = CompTypes.V_SHAPE_LIST.value
                        
                num_v += 1
                
    num_lin = 0
    if linpol_cadence:
        full_comp_counter.num_q_list_fluxes = np.full(total_comps, np.nan)
        full_comp_counter.num_u_list_fluxes = np.full(total_comps, np.nan)
        full_comp_counter.num_p_list_fluxes = np.full(total_comps, np.nan)
        full_comp_counter.lin_comp_types = np.full(total_comps, np.nan)
        for comp in range(total_comps):
            ##We've hit the cadence number, so add a component
            if comp % linpol_cadence == 0:
                comp_type = full_comp_counter.comp_types[comp]
                
                ##Stick in a polarisation fraction component
                if num_lin % NUM_LINPOL_MODELS == 0:
                    if comp_type == CompTypes.POINT_POWER.value or comp_type == CompTypes.POINT_CURVE.value or comp_type == CompTypes.POINT_LIST.value:
                        full_comp_counter.lin_comp_types[comp] = CompTypes.LIN_POINT_POL_FRAC.value
                    elif comp_type == CompTypes.GAUSS_POWER.value or comp_type == CompTypes.GAUSS_CURVE.value or comp_type == CompTypes.GAUSS_LIST.value:
                        full_comp_counter.lin_comp_types[comp] = CompTypes.LIN_GAUSS_POL_FRAC.value
                    else: 
                        full_comp_counter.lin_comp_types[comp] = CompTypes.LIN_SHAPE_POL_FRAC.value
                ##Stick in a power-law component
                elif num_lin % NUM_LINPOL_MODELS == 1:
                    if comp_type == CompTypes.POINT_POWER.value or comp_type == CompTypes.POINT_CURVE.value or comp_type == CompTypes.POINT_LIST.value:
                        full_comp_counter.lin_comp_types[comp] = CompTypes.LIN_POINT_POWER.value
                    elif comp_type == CompTypes.GAUSS_POWER.value or comp_type == CompTypes.GAUSS_CURVE.value or comp_type == CompTypes.GAUSS_LIST.value:
                        full_comp_counter.lin_comp_types[comp] = CompTypes.LIN_GAUSS_POWER.value
                    else: 
                        full_comp_counter.lin_comp_types[comp] = CompTypes.LIN_SHAPE_POWER.value
                ##Stick in a curved power-law component
                elif num_lin % NUM_LINPOL_MODELS == 2:
                    if comp_type == CompTypes.POINT_POWER.value or comp_type == CompTypes.POINT_CURVE.value or comp_type == CompTypes.POINT_LIST.value:
                        full_comp_counter.lin_comp_types[comp] = CompTypes.LIN_POINT_CURVE.value
                    elif comp_type == CompTypes.GAUSS_POWER.value or comp_type == CompTypes.GAUSS_CURVE.value or comp_type == CompTypes.GAUSS_LIST.value:
                        full_comp_counter.lin_comp_types[comp] = CompTypes.LIN_GAUSS_CURVE.value
                    else: 
                        full_comp_counter.lin_comp_types[comp] = CompTypes.LIN_SHAPE_CURVE.value
                        
                elif num_lin % NUM_LINPOL_MODELS == 3:
                    full_comp_counter.num_q_list_fluxes[comp] = settings.linpol_num_list
                    full_comp_counter.num_u_list_fluxes[comp] = settings.linpol_num_list
                    if comp_type == CompTypes.POINT_POWER.value or comp_type == CompTypes.POINT_CURVE.value or comp_type == CompTypes.POINT_LIST.value:
                        full_comp_counter.lin_comp_types[comp] = CompTypes.LIN_POINT_LIST.value
                    elif comp_type == CompTypes.GAUSS_POWER.value or comp_type == CompTypes.GAUSS_CURVE.value or comp_type == CompTypes.GAUSS_LIST.value:
                        full_comp_counter.lin_comp_types[comp] = CompTypes.LIN_GAUSS_LIST.value
                    else: 
                        full_comp_counter.lin_comp_types[comp] = CompTypes.LIN_SHAPE_LIST.value
                        
                else:
                    full_comp_counter.num_p_list_fluxes[comp] = settings.linpol_num_p_list
                    if comp_type == CompTypes.POINT_POWER.value or comp_type == CompTypes.POINT_CURVE.value or comp_type == CompTypes.POINT_LIST.value:
                        full_comp_counter.lin_comp_types[comp] = CompTypes.LIN_POINT_P_LIST.value
                    elif comp_type == CompTypes.GAUSS_POWER.value or comp_type == CompTypes.GAUSS_CURVE.value or comp_type == CompTypes.GAUSS_LIST.value:
                        full_comp_counter.lin_comp_types[comp] = CompTypes.LIN_GAUSS_P_LIST.value
                    else: 
                        full_comp_counter.lin_comp_types[comp] = CompTypes.LIN_SHAPE_P_LIST.value
                        
                num_lin += 1
                
    return full_comp_counter

class Expec_Counter(object):
    """something to count expected values"""
    
    def __init__(self):
        self.point_comp_accum = 0
        self.point_power_accum = 0
        self.point_curve_accum = 0
        self.point_list_accum = 0
        
        self.gauss_comp_accum = 0
        self.gauss_power_accum = 0
        self.gauss_curve_accum = 0
        self.gauss_list_accum = 0
        
##Vehicle for running tests
class BaseChunkTest(unittest.TestCase):
    """Test `wodenpy.skymodel.chunk_sky_model.map_chunk_pointgauss`
    works correctly"""
    
    def check_polarisation_chunking(self, full_comp_counter,
                                    expec_power_inds, expec_curve_inds,
                                    expec_list_inds, components,
                                    settings):
        
        if settings.linpol_cadence:
            ##this is every orig index that would have had a linear comp
            orig_lin_inds = np.arange(0, full_comp_counter.total_comps, settings.linpol_cadence)
            ##we cycled through pol frac, power, curve, so select subsets
            lin_pol_frac_inds = orig_lin_inds[np.arange(0, len(orig_lin_inds), NUM_LINPOL_MODELS)]
            lin_power_inds = orig_lin_inds[np.arange(1, len(orig_lin_inds), NUM_LINPOL_MODELS)]
            lin_curve_inds = orig_lin_inds[np.arange(2, len(orig_lin_inds), NUM_LINPOL_MODELS)]
            lin_list_inds = orig_lin_inds[np.arange(3, len(orig_lin_inds), NUM_LINPOL_MODELS)]
            lin_p_list_inds = orig_lin_inds[np.arange(4, len(orig_lin_inds), NUM_LINPOL_MODELS)]
            
            ##As we know which orig components we expect, we can just check if
            ##they had polarisation info by matching the subset arrays here
            
            all_expec_inds = np.concatenate((expec_power_inds, expec_curve_inds,
                                            expec_list_inds))
            
            expec_lin_pol_frac_inds = np.intersect1d(lin_pol_frac_inds, all_expec_inds)
            expec_lin_power_inds = np.intersect1d(lin_power_inds, all_expec_inds)
            expec_lin_curve_inds = np.intersect1d(lin_curve_inds, all_expec_inds)
            expec_lin_list_inds = np.intersect1d(lin_list_inds, all_expec_inds)
            expec_lin_p_list_inds = np.intersect1d(lin_p_list_inds, all_expec_inds)
            
            npt.assert_array_equal(expec_lin_pol_frac_inds, components.lin_pol_frac_orig_inds)
            npt.assert_array_equal(expec_lin_power_inds, components.lin_power_orig_inds)
            npt.assert_array_equal(expec_lin_curve_inds, components.lin_curve_orig_inds)
            npt.assert_array_equal(expec_lin_list_inds, components.lin_list_orig_inds)
            npt.assert_array_equal(expec_lin_p_list_inds, components.lin_p_list_orig_inds)
            
            expec_num_q_flux_entires = len(expec_lin_list_inds) * settings.linpol_num_list
            self.assertEqual(expec_num_q_flux_entires, components.total_num_q_flux_entires)
            
            expec_num_u_flux_entires = len(expec_lin_list_inds) * settings.linpol_num_list
            self.assertEqual(expec_num_u_flux_entires, components.total_num_u_flux_entires)
            
            expec_num_lin_p_flux_entires = len(expec_lin_p_list_inds) * settings.linpol_num_p_list
            self.assertEqual(expec_num_lin_p_flux_entires, components.total_num_lin_p_flux_entires)
            
            
        if settings.stokesV_cadence:
            ##this is every orig index that would have had a linear comp
            orig_v_inds = np.arange(0, full_comp_counter.total_comps, settings.stokesV_cadence)
            ##we cycled through pol frac, power, curve, so select subsets
            v_pol_frac_inds = orig_v_inds[np.arange(0, len(orig_v_inds), NUM_STOKESV_MODELS)]
            v_power_inds = orig_v_inds[np.arange(1, len(orig_v_inds), NUM_STOKESV_MODELS)]
            v_curve_inds = orig_v_inds[np.arange(2, len(orig_v_inds), NUM_STOKESV_MODELS)]
            v_list_inds = orig_v_inds[np.arange(3, len(orig_v_inds), NUM_STOKESV_MODELS)]
            
            ##As we know which orig components we expect, we can just check if
            ##they had polarisation info by matching the subset arrays here
            
            all_expec_inds = np.concatenate((expec_power_inds, expec_curve_inds,
                                            expec_list_inds))
            
            expec_v_pol_frac_inds = np.intersect1d(v_pol_frac_inds, all_expec_inds)
            expec_v_power_inds = np.intersect1d(v_power_inds, all_expec_inds)
            expec_v_curve_inds = np.intersect1d(v_curve_inds, all_expec_inds)
            expec_v_list_inds = np.intersect1d(v_list_inds, all_expec_inds)
            
            npt.assert_array_equal(expec_v_pol_frac_inds, components.v_pol_frac_orig_inds)
            npt.assert_array_equal(expec_v_power_inds, components.v_power_orig_inds)
            npt.assert_array_equal(expec_v_curve_inds, components.v_curve_orig_inds)
            npt.assert_array_equal(expec_v_list_inds, components.v_list_orig_inds)
            
            expec_num_flux_entries = len(expec_v_list_inds) * settings.stokesV_num_list
            self.assertEqual(expec_num_flux_entries, components.total_num_v_flux_entires)
    
    def check_pointgauss_chunking(self, chunk_ind : int, comps_per_chunk : int,
                              num_chunks : int, orig_n_comps : int,
                              settings : Skymodel_Settings,
                              comp_type : CompTypes,
                              full_comp_counter : Component_Type_Counter,
                              chunk_map : Skymodel_Chunk_Map,
                              expec_counter : Expec_Counter) -> Expec_Counter:
        """check that we've got the right chunking outcomes after doing the
        chunking
        
        orig_n_comps is how many each each of power, curve, and list type we
        originally put into the sky model"""
        
        linpol_cadence = settings.linpol_cadence
        stokesV_cadence = settings.stokesV_cadence
        num_list_values = settings.num_list_values
        
        if comp_type == CompTypes.POINT:
            total_n_powers = full_comp_counter.num_point_flux_powers
            total_n_curves = full_comp_counter.num_point_flux_curves
            total_n_lists = full_comp_counter.num_point_flux_lists
            total_n_comps = full_comp_counter.total_point_comps

            found_n_powers = chunk_map.n_point_powers
            found_n_curves = chunk_map.n_point_curves
            found_n_lists = chunk_map.n_point_lists
            
        elif comp_type == CompTypes.GAUSSIAN:
            total_n_powers = full_comp_counter.num_gauss_flux_powers
            total_n_curves = full_comp_counter.num_gauss_flux_curves
            total_n_lists = full_comp_counter.num_gauss_flux_lists
            total_n_comps = full_comp_counter.total_gauss_comps

            found_n_powers = chunk_map.n_gauss_powers
            found_n_curves = chunk_map.n_gauss_curves
            found_n_lists = chunk_map.n_gauss_lists
            
        chunk_comp_ind = chunk_ind*comps_per_chunk
        
        expec_n_comps = 0
        expec_n_powers = 0
        expec_n_curves = 0
        expec_n_lists = 0
        
        ##Go through a bunch of logic to work out what we would expect;
        ##this is mirroring what should be happening inside the function
        
        if (total_n_powers >= chunk_comp_ind + comps_per_chunk):
            expec_n_powers = comps_per_chunk
        
        else:
            ##if chunk_comp_ind >= total_n_powers, we should already have chunked
            ##all the power-law components
            if (chunk_comp_ind >= total_n_powers ):
                expec_n_powers = 0
            else:
                expec_n_powers = total_n_powers - chunk_comp_ind

            ##update the index where this comp begins with however many power-law
            ##components we found

            chunk_comp_ind += expec_n_powers

            ##If there are enough curved power laws to fill the rest of this chunk
            ##take off how ever many power-law components we already have
            if (total_n_curves >= comps_per_chunk + chunk_comp_ind - total_n_powers ):
                expec_n_curves = comps_per_chunk - expec_n_powers
            
            else:
                ##There are some curve components left, but not enough to fill the chunk
                if (total_n_curves > chunk_comp_ind - total_n_powers) and (total_n_curves != 0):

                    expec_n_curves = total_n_curves - chunk_comp_ind + total_n_powers
                else:
                    expec_n_curves = 0

                ##Update the lower index with how many curved power law sources we added
                chunk_comp_ind += expec_n_curves

                if (total_n_lists >= comps_per_chunk + chunk_comp_ind - total_n_powers - total_n_curves ):
                    expec_n_lists = comps_per_chunk - expec_n_powers - expec_n_curves
                
                ##There are some curve components left, but not enough to fill the chunk
                elif (total_n_lists > chunk_comp_ind - total_n_powers - total_n_curves) and (total_n_lists != 0):

                    expec_n_lists = total_n_lists - chunk_comp_ind + total_n_powers + total_n_curves
                else:
                    expec_n_lists = 0
                            
        expec_n_comps = expec_n_powers + expec_n_curves + expec_n_lists
        
        ##the gaussians go after the point sources, which effects their
        ##indexing. When considering just POINTs, this doesn't matter
        point_offset = 0
        
        ##check specific numbers for POINT and grab generic numbers to check
        if comp_type == CompTypes.POINT:
            
            self.assertEqual(expec_n_comps, chunk_map.n_points)
            self.assertEqual(expec_n_powers, chunk_map.n_point_powers)
            self.assertEqual(expec_n_curves, chunk_map.n_point_curves)
            self.assertEqual(expec_n_lists, chunk_map.n_point_lists)
            
            n_power_accum = expec_counter.point_power_accum
            n_curve_accum = expec_counter.point_curve_accum
            n_list_accum = expec_counter.point_list_accum
            
            components = chunk_map.point_components
            
        ##check specific numbers for GAUSSIAN and grab generic numbers to check
        elif comp_type == CompTypes.GAUSSIAN:
            
            self.assertEqual(expec_n_comps, chunk_map.n_gauss)
            self.assertEqual(expec_n_powers, chunk_map.n_gauss_powers)
            self.assertEqual(expec_n_curves, chunk_map.n_gauss_curves)
            self.assertEqual(expec_n_lists, chunk_map.n_gauss_lists)
            
            n_power_accum = expec_counter.gauss_power_accum
            n_curve_accum = expec_counter.gauss_curve_accum
            n_list_accum = expec_counter.gauss_list_accum
            
            components = chunk_map.gauss_components
            
            ##how many point sources come before the gaussians, needed to offset
            ##indexes and line numbers
            point_offset = full_comp_counter.total_point_comps
            
        ##Ok, in the original, we put the power first, curve second, and
        ##then list types third. So given how many of each type we have
        ##accumulated, we can predict what the indexes of sources are
        
        expec_power_inds = point_offset + np.arange(expec_n_powers) + n_power_accum
        expec_curve_inds = point_offset + orig_n_comps + np.arange(expec_n_curves) + n_curve_accum
        expec_list_inds = point_offset + 2*orig_n_comps + np.arange(expec_n_lists) + n_list_accum
        
        ##how many flux list entries we should have in total
        expec_num_flux_entries = expec_n_lists * num_list_values
        
        ##again, depending on how many components we have gone through during
        ##chunking, we can extrapolate what the lowest file number should be
        ##as we set that equal to component index
        expec_lowest_file_num = np.nan
        
        if n_power_accum < total_n_powers:
            expec_lowest_file_num = point_offset + n_power_accum
        elif n_curve_accum < total_n_curves:
            expec_lowest_file_num = point_offset + total_n_powers + n_curve_accum
        else:
            expec_lowest_file_num = point_offset + 2*total_n_powers + n_list_accum
        
        ##check everything is correct  
        self.assertTrue(np.array_equal(expec_power_inds, components.power_orig_inds))
        self.assertTrue(np.array_equal(expec_curve_inds, components.curve_orig_inds))
        self.assertTrue(np.array_equal(expec_list_inds, components.list_orig_inds))
        self.assertEqual(expec_num_flux_entries, components.total_num_flux_entires)
        # self.assertEqual(expec_lowest_file_num, components.lowest_file_num)
        
        self.check_polarisation_chunking(full_comp_counter,
                                         expec_power_inds, expec_curve_inds,
                                         expec_list_inds, components,
                                         settings)
            
        ##now update the accumulated components
        if comp_type == CompTypes.POINT:
            
            expec_counter.point_comp_accum += expec_n_comps
            expec_counter.point_power_accum += expec_n_powers
            expec_counter.point_curve_accum += expec_n_curves
            expec_counter.point_list_accum += expec_n_lists
            
            components = chunk_map.point_components
            
        elif comp_type == CompTypes.GAUSSIAN:
            
            expec_counter.gauss_comp_accum += expec_n_comps
            expec_counter.gauss_power_accum += expec_n_powers
            expec_counter.gauss_curve_accum += expec_n_curves
            expec_counter.gauss_list_accum += expec_n_lists
        
        return expec_counter
    
    
    def check_shapelet_chunking(self, chunk_ind : int, num_coeff_per_shape : int,
                                coeffs_per_chunk : int, orig_n_comps : int,
                                settings : Skymodel_Settings,
                                full_comp_counter : Component_Type_Counter,
                                chunk_map : Skymodel_Chunk_Map,
                                total_point_comps = 0, 
                                total_gauss_comps = 0,
                                do_check = True):
        """check that we've got the right chunking outcomes after doing the
        chunking
        
        orig_n_comps is how many each each of power, curve, and list type we
        originally put into the sky model"""
        
        
        ##always jam orig_n_comps of all three flux types into the test
        ##in order of power, curve, and list. Each of those components is given
        ##coeffs per chunk. So we should be able to logic what we expect given
        ##those numbers
        
        low_coeff_ind = chunk_ind * coeffs_per_chunk
        high_coeff_ind = (chunk_ind + 1) * coeffs_per_chunk - 1
        
        ##if the high coeff is greater than the total number of coeffs, just
        ##means that we don't have a full final chunk. So reduce high coeff if
        ##needed
        
        if high_coeff_ind >= full_comp_counter.total_shape_basis:
            high_coeff_ind = full_comp_counter.total_shape_basis - 1
        
        low_comp_ind = int(np.floor(low_coeff_ind / num_coeff_per_shape))
        high_comp_ind = int(np.floor(high_coeff_ind / num_coeff_per_shape))
        
        expec_power_inds = np.array([])
        expec_curve_inds = np.array([])
        expec_list_inds = np.array([])
        
        ##only power laws
        if high_comp_ind < orig_n_comps:
            # expec_n_powers = (high_comp_ind - low_comp_ind) + 1
            expec_power_inds = np.arange(low_comp_ind, high_comp_ind+1)
            expec_n_powers = len(expec_power_inds)
            
        ##both power and curve
        elif low_comp_ind < orig_n_comps and high_comp_ind < 2*orig_n_comps:
            # expec_n_powers = (orig_n_comps - low_comp_ind)
            # expec_n_curves = (high_comp_ind - orig_n_comps) + 1
            expec_power_inds = np.arange(low_comp_ind, orig_n_comps)
            expec_curve_inds = np.arange(orig_n_comps, high_comp_ind + 1)
            
            
        ##all three exist
        elif low_comp_ind < orig_n_comps and high_comp_ind >= 2*orig_n_comps:
            expec_power_inds = np.arange(low_comp_ind, orig_n_comps)
            expec_curve_inds = np.arange(orig_n_comps, 2*orig_n_comps)
            expec_list_inds = np.arange(2*orig_n_comps, 2*orig_n_comps + 1)
            
        ##defo have curves
        elif low_comp_ind < 2*orig_n_comps:
            ##also have lists
            if high_comp_ind >= 2*orig_n_comps:
                expec_curve_inds = np.arange(low_comp_ind, 2*orig_n_comps)
                expec_list_inds = np.arange(2*orig_n_comps, high_comp_ind + 1)
            ##only have curves
            else:
                expec_curve_inds = np.arange(low_comp_ind, high_comp_ind + 1)
        
        ##should ONLY get here if we just have lists
        else:
            expec_list_inds = np.arange(low_comp_ind, high_comp_ind + 1)
            
        ##If there were point or shapelet models in the original model, during
        ##testing we will have put them first. So offset the predictions by
        ##however many we have
        
        expec_power_inds += total_point_comps + total_gauss_comps
        expec_curve_inds += total_point_comps + total_gauss_comps
        expec_list_inds += total_point_comps + total_gauss_comps
            
        expec_n_powers = len(expec_power_inds)
        expec_n_curves = len(expec_curve_inds)
        expec_n_lists = len(expec_list_inds)
            
        expec_n_comp = expec_n_powers + expec_n_curves + expec_n_lists
        
        expec_num_flux_entries = expec_n_lists*settings.num_list_values
        # # ##check everything is correct  
        
        if do_check:
            components = chunk_map.shape_components
            self.assertEqual(expec_n_powers, chunk_map.n_shape_powers)
            self.assertEqual(expec_n_curves, chunk_map.n_shape_curves)
            self.assertEqual(expec_n_lists, chunk_map.n_shape_lists)
            self.assertTrue(np.array_equal(expec_power_inds, components.power_orig_inds))
            self.assertTrue(np.array_equal(expec_curve_inds, components.curve_orig_inds))
            self.assertTrue(np.array_equal(expec_list_inds, components.list_orig_inds))
            self.assertEqual(expec_num_flux_entries, components.total_num_flux_entires)
        
        if expec_n_powers:
            min_power = expec_power_inds.min()
        else:
            min_power = np.nan
            
        if expec_n_curves:
            min_curve = expec_curve_inds.min()
        else:
            min_curve = np.nan
            
        if expec_n_lists:
            min_list = expec_list_inds.min()
        else:
            min_list = np.nan
        
        expec_lowest_file_num = np.nanmin([min_power, min_curve, min_list])
        # self.assertEqual(expec_lowest_file_num, components.lowest_file_num)
        
        if do_check:
            self.check_polarisation_chunking(full_comp_counter,
                                         expec_power_inds, expec_curve_inds,
                                         expec_list_inds, components,
                                         settings)
            
        return expec_n_comp
            
        
class classname(object):
    """
    something to hold expect nums of things for a test sky model
    """
    
    def __init__(self, deg_between_comps : float,
                       num_coeff_per_shape : int,
                       num_list_values : int,
                       comps_per_source : int):
        """
        count le things
        """
        self.num_comps

class Expected_Components(object):
    
    def __init__(self, comp_type : CompTypes, num_chunk_power = 0,
                 num_chunk_curve = 0, num_chunk_list = 0,
                 num_list_values = 0, comps_per_chunk = 0,
                 num_v_pol_frac = 0, num_v_power = 0, num_v_curve = 0,
                 num_v_list = 0, num_v_list_values = 0,
                 num_lin_pol_frac = 0, num_lin_power = 0, num_lin_curve = 0,
                 num_lin_list = 0, num_lin_list_values = 0, 
                 num_lin_p_list = 0, num_lin_p_list_values = 0):
        """
        docstring
        """
        
        num_comps = num_chunk_power + num_chunk_curve + num_chunk_list
    
        self.ras = np.empty(num_comps)
        self.decs = np.empty(num_comps)
        
        self.power_ref_freqs = np.empty(num_chunk_power)
        self.power_ref_stokesI = np.empty(num_chunk_power)
        self.power_SIs = np.empty(num_chunk_power)
        self.power_comp_inds = np.empty(num_chunk_power)
        
        self.curve_ref_freqs = np.empty(num_chunk_curve)
        self.curve_ref_stokesI = np.empty(num_chunk_curve)
        self.curve_SIs = np.empty(num_chunk_curve)
        self.curve_qs = np.empty(num_chunk_curve)
        self.curve_comp_inds = np.empty(num_chunk_curve)
        
        self.list_freqs = np.empty(num_list_values*num_chunk_list)
        self.list_stokesI = np.empty(num_list_values*num_chunk_list)
        self.list_stokesQ = np.empty(num_list_values*num_chunk_list)
        self.list_stokesU = np.empty(num_list_values*num_chunk_list)
        self.list_stokesV = np.empty(num_list_values*num_chunk_list)
        self.num_list_values = np.empty(num_chunk_list)
        self.list_start_indexes = np.empty(num_chunk_list)
        self.list_comp_inds = np.empty(num_chunk_list)
        
        if comp_type == CompTypes.GAUSSIAN or comp_type == CompTypes.SHAPELET:
            self.minors = np.empty(num_comps)
            self.majors = np.empty(num_comps)
            self.pas = np.empty(num_comps)
        if comp_type == CompTypes.SHAPELET:
            self.param_indexes = np.empty(comps_per_chunk)
            self.n1s = np.empty(comps_per_chunk)
            self.n2s = np.empty(comps_per_chunk)
            self.shape_coeffs = np.empty(comps_per_chunk)
        
        
        ##circular polarisation stuff===========================================
        
        self.stokesV_pol_fracs = np.empty(num_v_pol_frac)
        self.stokesV_pol_frac_comp_inds = np.empty(num_v_pol_frac, dtype=int)
        
        self.stokesV_power_ref_flux = np.empty(num_v_power)
        self.stokesV_power_SIs = np.empty(num_v_power)
        self.stokesV_power_comp_inds = np.empty(num_v_power, dtype=int)
        
        self.stokesV_curve_ref_flux = np.empty(num_v_curve)
        self.stokesV_curve_SIs = np.empty(num_v_curve)
        self.stokesV_curve_qs = np.empty(num_v_curve)
        self.stokesV_curve_comp_inds = np.empty(num_v_curve, dtype=int)
        
        self.stokesV_list_ref_flux = np.empty(num_v_list*num_v_list_values)
        self.stokesV_list_ref_freqs = np.empty(num_v_list*num_v_list_values)
        self.stokesV_list_comp_inds = np.empty(num_v_list, dtype=int)
        
         ##linear polarisation stuff===========================================
         
        ##These values are needed for any linear polarisation model
        ##These get used after the polarised flux has been calculated, so need
        ##their own indexes to refer to the polarised flux
        self.rm_values = np.empty(num_lin_pol_frac + num_lin_power + num_lin_curve + num_lin_p_list)
        self.intr_pol_angle = np.empty(num_lin_pol_frac + num_lin_power + num_lin_curve + num_lin_p_list)
        # self.linpol_angle_inds = np.empty(num_lin_pol_frac + num_lin_power + num_lin_curve + num_lin_p_list, dtype=int)
        self.linpol_angle_inds = np.full(num_lin_pol_frac + num_lin_power + num_lin_curve + num_lin_p_list, -1, dtype=int)
        
        ##The following are different models for linear polarisation
        self.linpol_pol_fracs = np.empty(num_lin_pol_frac)
        self.linpol_pol_frac_comp_inds = np.empty(num_lin_pol_frac, dtype=int)
        
        self.linpol_power_ref_flux = np.empty(num_lin_power)
        self.linpol_power_SIs = np.empty(num_lin_power)
        self.linpol_power_comp_inds = np.empty(num_lin_power, dtype=int)
        
        self.linpol_curve_ref_flux = np.empty(num_lin_curve)
        self.linpol_curve_SIs = np.empty(num_lin_curve)
        self.linpol_curve_qs = np.empty(num_lin_curve)
        self.linpol_curve_comp_inds = np.empty(num_lin_curve, dtype=int)
        
        self.linpol_list_ref_flux = np.empty(num_lin_list*num_lin_list_values)
        self.linpol_list_ref_freqs = np.empty(num_lin_list*num_lin_list_values)
        self.linpol_list_comp_inds = np.empty(num_lin_list, dtype=int)
        
        self.linpol_p_list_ref_flux = np.empty(num_lin_p_list*num_lin_p_list_values)
        self.linpol_p_list_ref_freqs = np.empty(num_lin_p_list*num_lin_p_list_values)
        self.linpol_p_list_comp_inds = np.empty(num_lin_p_list, dtype=int)
        
        self.n_v_pol_frac = num_v_pol_frac
        self.n_v_power = num_v_power
        self.n_v_curve = num_v_curve
        self.n_v_list = num_v_list
        self.n_lin_pol_frac = num_lin_pol_frac
        self.n_lin_power = num_lin_power
        self.n_lin_curve = num_lin_curve
        self.n_lin_list = num_lin_list
        self.n_lin_p_list = num_lin_p_list
                        
class Expected_Sky_Chunk(object):
    
    def __init__(self):
        """
        laksjfdkljshadloasd
        """
        self.n_points = 0
        self.n_point_lists = 0
        self.n_point_powers = 0
        self.n_point_curves = 0
        self.n_gauss = 0
        self.n_gauss_lists = 0
        self.n_gauss_powers = 0
        self.n_gauss_curves = 0
        self.n_shapes = 0
        self.n_shape_lists = 0
        self.n_shape_powers = 0
        self.n_shape_curves = 0
        self.n_shape_coeffs = 0
        self.n_comps = 0
        
        self.point_components = None
        self.gauss_components = None
        self.shape_components = None
        
        self.n_v_point_pol_frac = 0
        self.n_v_point_power = 0
        self.n_v_point_curve = 0
        self.n_lin_point_pol_frac = 0
        self.n_lin_point_power = 0
        self.n_lin_point_curve = 0
        
        self.n_v_gauss_pol_frac = 0
        self.n_v_gauss_power = 0
        self.n_v_gauss_curve = 0
        self.n_lin_gauss_pol_frac = 0
        self.n_lin_gauss_power = 0
        self.n_lin_gauss_curve = 0
        
        self.n_v_shape_pol_frac = 0
        self.n_v_shape_power = 0
        self.n_v_shape_curve = 0
        self.n_lin_shape_pol_frac = 0
        self.n_lin_shape_power = 0
        self.n_lin_shape_curve = 0
        
    def init_point_components(self, num_chunk_power : int,
                 num_chunk_curve : int, num_chunk_list : int,
                 num_list_values : int, comps_per_chunk : int,
                 num_v_pol_frac = 0, num_v_power = 0, num_v_curve = 0,
                 num_v_list = 0, num_v_list_values = 0,
                 num_lin_pol_frac = 0, num_lin_power = 0, num_lin_curve = 0,
                 num_lin_list = 0, num_lin_list_values = 0, 
                 num_lin_p_list = 0, num_lin_p_list_values = 0):
        
        self.n_point_powers = num_chunk_power
        self.n_point_curves = num_chunk_curve
        self.n_point_lists = num_chunk_list
        self.n_points = num_chunk_power + num_chunk_curve + num_chunk_list
        
        self.point_components = Expected_Components(CompTypes.POINT,
                                 num_chunk_power, num_chunk_curve, num_chunk_list,
                                 num_list_values, comps_per_chunk,
                                 num_v_pol_frac, num_v_power, num_v_curve,
                                 num_v_list, num_v_list_values,
                                 num_lin_pol_frac, num_lin_power, num_lin_curve,
                                 num_lin_list, num_lin_list_values,
                                 num_lin_p_list, num_lin_p_list_values)
        
    def init_gauss_components(self, num_chunk_power : int,
                 num_chunk_curve : int, num_chunk_list : int,
                 num_list_values : int, comps_per_chunk : int,
                 num_v_pol_frac = 0, num_v_power = 0, num_v_curve = 0,
                 num_v_list = 0, num_v_list_values = 0,
                 num_lin_pol_frac = 0, num_lin_power = 0, num_lin_curve = 0,
                 num_lin_list = 0, num_lin_list_values = 0, 
                 num_lin_p_list = 0, num_lin_p_list_values = 0):
        
        self.gauss_components = Expected_Components(CompTypes.GAUSSIAN,
                                 num_chunk_power, num_chunk_curve, num_chunk_list,
                                 num_list_values, comps_per_chunk,
                                 num_v_pol_frac, num_v_power, num_v_curve,
                                 num_v_list, num_v_list_values,
                                 num_lin_pol_frac, num_lin_power, num_lin_curve,
                                 num_lin_list, num_lin_list_values,
                                 num_lin_p_list, num_lin_p_list_values)
        
        self.n_gauss_powers = num_chunk_power
        self.n_gauss_curves = num_chunk_curve
        self.n_gauss_lists = num_chunk_list
        self.n_gauss = num_chunk_power + num_chunk_curve + num_chunk_list
        
    def init_shape_components(self, num_chunk_power : int,
                 num_chunk_curve : int, num_chunk_list : int,
                 num_list_values : int, comps_per_chunk : int):
        
        self.shape_components = Expected_Components(CompTypes.SHAPELET,
                                 num_chunk_power, num_chunk_curve, num_chunk_list,
                                 num_list_values, comps_per_chunk)
        
        
        ##with shapelets, easier to start things as lists and append
        ##polarisation information
        
        self.shape_components.stokesV_pol_fracs = []
        self.shape_components.stokesV_pol_frac_comp_inds = []
        self.shape_components.stokesV_power_ref_flux = []
        self.shape_components.stokesV_power_SIs = []
        self.shape_components.stokesV_power_comp_inds = []
        self.shape_components.stokesV_curve_ref_flux = []
        self.shape_components.stokesV_curve_SIs = []
        self.shape_components.stokesV_curve_qs = []
        self.shape_components.stokesV_curve_comp_inds = []
        self.shape_components.rm_values = []
        self.shape_components.intr_pol_angle = []
        self.shape_components.linpol_angle_inds = []
        self.shape_components.linpol_pol_fracs = []
        self.shape_components.linpol_pol_frac_comp_inds = []
        self.shape_components.linpol_power_ref_flux = []
        self.shape_components.linpol_power_SIs = []
        self.shape_components.linpol_power_comp_inds = []
        self.shape_components.linpol_curve_ref_flux = []
        self.shape_components.linpol_curve_SIs = []
        self.shape_components.linpol_curve_qs = []
        self.shape_components.linpol_curve_comp_inds = []
        self.shape_components.stokesV_list_comp_inds = []
        self.shape_components.stokesV_list_ref_flux = []
        self.shape_components.stokesV_list_ref_freqs = []
        self.shape_components.linpol_list_comp_inds = []
        self.shape_components.linpol_list_ref_flux = []
        self.shape_components.linpol_list_ref_freqs = []
        self.shape_components.linpol_p_list_comp_inds = []
        self.shape_components.linpol_p_list_ref_flux = []
        self.shape_components.linpol_p_list_ref_freqs = []
        
        
        self.n_shape_powers = num_chunk_power
        self.n_shape_curves = num_chunk_curve
        self.n_shape_lists = num_chunk_list
        self.n_shapes = num_chunk_power + num_chunk_curve + num_chunk_list
        
    def _print_things(self):
        print("n_points", self.n_points)
        print("n_point_lists", self.n_point_lists)
        print("n_point_powers", self.n_point_powers)
        print("n_point_curves", self.n_point_curves)
        print("n_gauss", self.n_gauss)
        print("n_gauss_lists", self.n_gauss_lists)
        print("n_gauss_powers", self.n_gauss_powers)
        print("n_gauss_curves", self.n_gauss_curves)
        print("n_shapes", self.n_shapes)
        print("n_shape_lists", self.n_shape_lists)
        print("n_shape_powers", self.n_shape_powers)
        print("n_shape_curves", self.n_shape_curves)
        print("n_shape_coeffs", self.n_shape_coeffs)
        print("n_comps", self.n_comps)