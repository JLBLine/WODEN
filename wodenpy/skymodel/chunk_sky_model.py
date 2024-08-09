
import numpy as np
import sys
import os
from enum import Enum

from wodenpy.skymodel.woden_skymodel import Component_Type_Counter, CompTypes
    
NUM_FLUX_TYPES = 3
    
class Components_Map(object):
    """
    Mapping information for a set of components, of either POINT,
    GAUSSIAN or SHAPELET type.
    
    
    :cvar bool power_orig_inds: Relative indexes for all power-law components w.r.t to the original sky model.
    :cvar bool curve_orig_inds: Relative indexes for all curved power-law components w.r.t to the original sky model.
    :cvar bool list_orig_inds: Relative indexes for all list-type components w.r.t to the original sky model.
    :cvar bool power_shape_basis_inds: Index of a basis function entry relative to its component for power components.
    :cvar bool curve_shape_basis_inds: Index of a basis function entry relative to its component for curve components.
    :cvar bool list_shape_basis_inds: Index of a basis function entry relative to its component for list components.
    :cvar int lowest_file_num: The line in the original sky model file that each component appears in. Ignore all lines before the smallest line number for all components in this chunk, makes reading faster.
    :cvar int total_num_flux_entires: Number of flux list entries there are in total so we can allocate correct amount when reading in full information.
    :cvar int total_shape_coeffs: Number of shapelet basis functions that are in total so we can allocate correct amount when reading in full information.
    :cvar int num_v_pol_fracs: Number of Stokes V fractional polarisation types.
    :cvar int num_v_powers: Number of Stokes V power-law components.
    :cvar int num_v_curves: Number of Stokes V curved power-law components.
    :cvar bool v_pol_frac_orig_inds: Relative indexes for all Stokes V fractional polarisation components w.r.t the original sky model.
    :cvar bool v_power_orig_inds: Relative indexes for all Stokes V power-law components w.r.t the original sky model
    :cvar bool v_curve_orig_inds: Relative indexes for all Stokes V curved power-law components w.r.t the original sky model
    :cvar int num_lin_pol_fracs: Number of Linear polarisation fractional polarisation types.
    :cvar int num_lin_powers: Number of Linear polarisation power-law components.
    :cvar int num_lin_curves: Number of Linear polarisation curved power-law components.
    :cvar int num_lin_angles: Total number of Linear polarisation angle components (num_lin_pol_fracs + num_lin_powers + num_lin_curves).
    :cvar bool lin_pol_frac_orig_inds: Relative indexes for all Linear polarisation fractional polarisation components w.r.t the original sky model.
    :cvar bool lin_power_orig_inds: Relative indexes for all Linear polarisation power-law components w.r.t the original sky model
    :cvar bool lin_curve_orig_inds: Relative indexes for all Linear polarisation curved power-law components w.r.t the original sky model
    
    """
    ##Mapping information for a set of components, of either POINT,
    ##GAUSSIAN or SHAPELET type
    
    def __init__(self):
        """
        Setup required fields
        """
        
        ##These are relative indexes w.r.t all components in the original
        ##sky model
        self.power_orig_inds = False
        self.curve_orig_inds = False
        self.list_orig_inds = False
        
        ##These are only used for a SHAPELET component
        ##They map the index of a basis function entry relative to it's component
        self.power_shape_basis_inds = False
        self.curve_shape_basis_inds = False
        self.list_shape_basis_inds = False
        
        ##when reading file back, quickest to stick into one single
        ##array with
        
        ##the line in the original sky model file that each component
        ##appears in. Ignore all lines before the smallest line number
        ##for all components in this chunk, makes reading faster
        ##only applied to txt/yaml models
        self.lowest_file_num = np.nan
        
        ##use this to count how many flux list entries there are in total
        ##so we can allocate correct amount when reading in full information
        self.total_num_flux_entires = 0
        
        #use this to count how many shapelet basis functions that are in total
        #so we can allocate correct amount when reading in full information
        self.total_shape_coeffs = 0
        
        ##count polarised flux types
        #so we can allocate correct amount when reading in full information
        self.num_v_pol_fracs = 0
        self.num_v_powers = 0
        self.num_v_curves = 0
        self.num_v_lists = 0
        self.v_pol_frac_orig_inds = False
        self.v_power_orig_inds = False
        self.v_curve_orig_inds = False
        self.v_list_orig_inds = False
        
        self.total_num_v_flux_entires = 0
        
        self.num_lin_pol_fracs = 0
        self.num_lin_powers = 0
        self.num_lin_curves = 0
        self.num_lin_lists = 0
        self.num_lin_p_lists = 0
        self.num_lin_angles = 0
        self.lin_pol_frac_orig_inds = False
        self.lin_power_orig_inds = False
        self.lin_curve_orig_inds = False
        self.lin_list_orig_inds = False
        self.lin_p_list_orig_inds = False
        
        self.total_num_u_flux_entires = 0
        self.total_num_q_flux_entires = 0
        self.total_num_lin_p_flux_entires = 0
        
    
    
class Skymodel_Chunk_Map(object):
    """A class representing a chunk of a sky model, containing information about
    the number of point, Gaussian, and shape components, as well as their
    respective power-law, curved power-law, or list type flux info. This class
    also provides methods for consolidating component and flux types into one
    array of original component indexes, and for printing information
    about the chunk.
    
    :cvar int n_point_lists:
        Number of POINT list-type components.
    :cvar int n_point_powers:
        Number of POINT power-law components.
    :cvar int n_point_curves:
        Number of POINT curved power-law components.
    :cvar int n_gauss_lists:
        Number of GAUSSIAN list-type components.
    :cvar int n_gauss_powers:
        Number of GAUSSIAN power-law components.
    :cvar int n_gauss_curves:
        Number of GAUSSIAN curved power-law components.
    :cvar int n_shape_lists:
        Number of SHAPELET list-type components.
    :cvar int n_shape_powers:
        Number of SHAPELET power-law components.
    :cvar int n_shape_curves:
        Number of SHAPELET curved power-law components.
    :cvar int n_shape_coeffs:
        Number of SHAPELET coefficients.
    :cvar int n_points:
        Number of POINT components
    :cvar int n_gauss:
        Number of GAUSSIAN components
    :cvar int n_shapes:
        Number of SHAPELET components
    :cvar int n_comps:
        Number of all components
    :cvar Components_Map point_components:
        Mapping object for POINT components.
    :cvar Components_Map gauss_components:
        Mapping object for GAUSSIAN components.
    :cvar Components_Map shape_components:
        Mapping object for SHAPELET components.
    
    """

    def __init__(self, n_point_powers = 0, n_point_curves = 0, n_point_lists = 0,
                       n_gauss_powers = 0, n_gauss_curves = 0, n_gauss_lists = 0,
                       n_shape_powers = 0, n_shape_curves = 0, n_shape_lists = 0,
                       n_shape_coeffs = 0):
        """Setup everything with zeros as default"""
        
        self.n_point_lists = n_point_lists
        self.n_point_powers = n_point_powers
        self.n_point_curves = n_point_curves
        
        self.n_gauss_lists = n_gauss_lists
        self.n_gauss_powers = n_gauss_powers
        self.n_gauss_curves = n_gauss_curves
        
        self.n_shape_lists = n_shape_lists
        self.n_shape_powers = n_shape_powers
        self.n_shape_curves = n_shape_curves
        self.n_shape_coeffs = n_shape_coeffs
        
        self.n_points = n_point_lists + n_point_powers + n_point_curves
        self.n_gauss =  n_gauss_lists + n_gauss_powers + n_gauss_curves
        self.n_shapes =  n_shape_lists + n_shape_powers + n_shape_curves
        
        self.n_comps = self.n_points + self.n_gauss + self.n_shapes
        
        ##Setup the POINT, GAUSS, and SHAPE classes
        ##TODO set these up regardless of size as empty things take up
        ##small RAM??
        if self.n_points > 0:
            self.point_components = Components_Map()
        if self.n_gauss > 0:
            self.gauss_components = Components_Map()
        if self.n_shapes > 0:
            self.shape_components = Components_Map()
            
        self.lowest_file_number = np.nan
        
        ##Used to count how many basis function values have already been
        ##added, gets updated by `use_libwoden.add_info_to_source_catalogue`
        ##when reading in the full model from the catalogue file
        self.current_shape_basis_index = 0
        
        ##Does the catalogue have the INTR_POL_ANGLE column?
        self.has_intr_pol_angle = False
            
    def print_info(self):
        """
        Prints information about the ChunkSkyModel object, including the number of points, Gaussians, and shapes,
        as well as the number of powers, curves, lists, and coefficients associated with each type of object.
        """
        
        print("n_points", self.n_points)
        print("\tn_point_powers", self.n_point_powers)
        print("\tn_point_curves", self.n_point_curves)
        print("\tn_point_lists", self.n_point_lists)
        
        print("n_gauss", self.n_gauss)
        print("\tn_gauss_powers", self.n_gauss_powers)
        print("\tn_gauss_curves", self.n_gauss_curves)
        print("\tn_gauss_lists", self.n_gauss_lists)
        
        print("n_shapes", self.n_shapes)
        print("\tn_shape_powers", self.n_shape_powers)
        print("\tn_shape_curves", self.n_shape_curves)
        print("\tn_shape_lists", self.n_shape_lists)
        print("\tn_shape_coeffs", self.n_shape_coeffs)
    

def increment_flux_type_counters(power_iter : int, curve_iter : int,
                                 list_iter : int, num_chunk_power : int,
                                 num_chunk_curve : int, num_chunk_list : int,
                                 num_power : int, num_curve : int, num_list : int,
                                 comps_per_chunk : int,
                                 lower_comp_ind : int, upper_comp_ind : int):
    """
    Here, given the overall lower and upper index in a given type of components,
    work out how many of each flux type we have and increment the counters as appropriate

    Always order things as POWER_LAW, CURVED_POWER_LAW, LIST
    
    Parameters
    -----------
    power_iter : int
        The current iteration of the power law flux type.
    curve_iter : int
        The current iteration of the curved power law flux type.
    list_iter : int
        The current iteration of the list flux type.
    num_chunk_power : int
        The number of power law flux types in the current chunk.
    num_chunk_curve : int
        The number of curved power law flux types in the current chunk.
    num_chunk_list : int
        The number of list flux types in the current chunk.
    num_power : int
        The total number of power law flux types.
    num_curve : int
        The total number of curved power law flux types.
    num_list : int
        The total number of list flux types.
    comps_per_chunk : int
        The number of components per chunk.
    lower_comp_ind : int
        The lower index of the current chunk.
    upper_comp_ind : int
        The upper index of the current chunk.
    
    Returns
    --------
    Tuple of integers
        The updated values of power_iter, curve_iter, list_iter, num_chunk_power, num_chunk_curve, num_chunk_list.
    """

    remainder = 0
    lower_flux_ind = 0
    upper_flux_ind = 0

    ## Enough POWER_LAW to fill the whole chunk
    if (num_power > upper_comp_ind):
        num_chunk_power = comps_per_chunk
        power_iter = lower_comp_ind
        num_chunk_curve = 0
        num_chunk_list = 0
    ##Not enough POWER_LAW to fill the whole chunk
    else:
        ##There are enough POWER_LAW to partially fill chunk
        if (num_power >= lower_comp_ind):
            num_chunk_power = num_power - lower_comp_ind
            power_iter = lower_comp_ind

            ##How much is left to fill in this chunk
            remainder = comps_per_chunk - num_chunk_power

            ##If there are enough CURVED_POWER_LAW to fill rest of the chunk
            if (num_curve >= remainder):
                num_chunk_curve = remainder
                curve_iter = 0
                num_chunk_list = 0
            ##Not enough CURVED_POWER_LAW to fill rest of the chunk
            else:
            ##There are some CURVED_POWER_LAW to add
                if (num_curve < remainder) and (num_curve != 0):
                    num_chunk_curve = num_curve
                    curve_iter = 0
                    remainder -= num_curve

                ##There are enough LIST to fill the rest of the chunk
                if (num_list >= remainder):
                    num_chunk_list = remainder
                    list_iter = 0
                ##There are some LIST but not enough to fill the rest of the chunk
                elif (num_list != 0):
                    num_chunk_list = num_list
                    list_iter = 0

        ##There aren't any POWER_LAW to put in this chunk
        ##We may well have already chunked up a number of POWER_LAW so take
        ##off the number of POWER_LAW from the lower_comp_ind, upper_comp_ind
        else:
            lower_flux_ind = lower_comp_ind - num_power
            upper_flux_ind = upper_comp_ind - num_power

            ##There are enough CURVED_POWER_LAW to fill the rest of the chunk
            if (num_curve >= upper_flux_ind):
                num_chunk_curve = comps_per_chunk
                curve_iter = lower_flux_ind
                num_chunk_list = 0
            else:
            ##There are some CURVED_POWER_LAW to add
                if (num_curve > lower_flux_ind) and (num_curve != 0):
                    num_chunk_curve = num_curve - lower_flux_ind
                    curve_iter = lower_flux_ind
                    remainder = comps_per_chunk - num_chunk_curve

                    ##There are enough LIST to fill rest of chunk
                    if (num_list >= remainder):
                        num_chunk_list = remainder
                        list_iter = 0
                    ##There aren't enough LIST to fill chunk but there are some
                    else:
                        num_chunk_list = num_list
                        list_iter = 0
                ##There are no POWER_LAW or CURVED_POWER_LAW to add
                else:
                    lower_flux_ind = lower_comp_ind - num_power - num_curve
                    upper_flux_ind = upper_comp_ind - num_power - num_curve

                    ##There are enough LIST to fill the rest of the chunk
                    if (num_list > upper_flux_ind):
                        num_chunk_list = comps_per_chunk
                        list_iter = lower_flux_ind
                    
                    ##There are some LIST but not enough to fill the rest of the chunk
                    elif (num_list > lower_flux_ind):
                        num_chunk_list = num_list - lower_flux_ind
                        list_iter = lower_flux_ind

    return power_iter, curve_iter, list_iter, num_chunk_power, num_chunk_curve, num_chunk_list


def fill_chunk_map_polarised_info(comp_type : CompTypes,
        chunk_map : Skymodel_Chunk_Map,
        cropped_comp_counter : Component_Type_Counter) -> Skymodel_Chunk_Map:
    """This should be called within map_chunk_pointgauss or map_chunk_shapelet,
    as a bunch of things should be filled"""
    
    if comp_type == CompTypes.POINT:
    
        components = chunk_map.point_components
        v_pol_frac_inds = cropped_comp_counter.orig_v_point_pol_frac_inds
        v_power_inds = cropped_comp_counter.orig_v_point_power_inds
        v_curve_inds = cropped_comp_counter.orig_v_point_curve_inds
        v_list_inds = cropped_comp_counter.orig_v_point_list_inds
        lin_pol_frac_inds = cropped_comp_counter.orig_lin_point_pol_frac_inds
        lin_power_inds = cropped_comp_counter.orig_lin_point_power_inds
        lin_curve_inds = cropped_comp_counter.orig_lin_point_curve_inds
        lin_list_inds = cropped_comp_counter.orig_lin_point_list_inds
        lin_p_list_inds = cropped_comp_counter.orig_lin_point_p_list_inds
        
    elif comp_type == CompTypes.GAUSSIAN:
        
        components = chunk_map.gauss_components
        v_pol_frac_inds = cropped_comp_counter.orig_v_gauss_pol_frac_inds
        v_power_inds = cropped_comp_counter.orig_v_gauss_power_inds
        v_curve_inds = cropped_comp_counter.orig_v_gauss_curve_inds
        v_list_inds = cropped_comp_counter.orig_v_gauss_list_inds
        lin_pol_frac_inds = cropped_comp_counter.orig_lin_gauss_pol_frac_inds
        lin_power_inds = cropped_comp_counter.orig_lin_gauss_power_inds
        lin_curve_inds = cropped_comp_counter.orig_lin_gauss_curve_inds
        lin_list_inds = cropped_comp_counter.orig_lin_gauss_list_inds
        lin_p_list_inds = cropped_comp_counter.orig_lin_gauss_p_list_inds
        
    else:
        components = chunk_map.shape_components
        v_pol_frac_inds = cropped_comp_counter.orig_v_shape_pol_frac_inds
        v_power_inds = cropped_comp_counter.orig_v_shape_power_inds
        v_curve_inds = cropped_comp_counter.orig_v_shape_curve_inds
        v_list_inds = cropped_comp_counter.orig_v_shape_list_inds
        lin_pol_frac_inds = cropped_comp_counter.orig_lin_shape_pol_frac_inds
        lin_power_inds = cropped_comp_counter.orig_lin_shape_power_inds
        lin_curve_inds = cropped_comp_counter.orig_lin_shape_curve_inds
        lin_list_inds = cropped_comp_counter.orig_lin_shape_list_inds
        lin_p_list_inds = cropped_comp_counter.orig_lin_shape_p_list_inds
        
    ##Default things to being zero
    components.num_v_pol_fracs = 0
    components.num_v_powers = 0
    components.num_v_curves = 0
    components.num_v_lists = 0
    components.num_lin_pol_fracs = 0
    components.num_lin_powers = 0
    components.num_lin_curves = 0
    components.num_lin_lists = 0
    components.num_lin_p_lists = 0
    
    ##OK, we want to find if there was any polarisation information in this
    ##chunk, so we use np.intersect1d, which finds the common elements between
    ##two arrays
    ##Doesn't matter what Stokes I flux type, if any of them have polarisation
    ##info we want to know, so combine them into one array
    
    all_orig_inds = np.concatenate((components.power_orig_inds,
                                   components.curve_orig_inds,
                                   components.list_orig_inds))
    ##These index arrays will be False by default, so test if they are a
    ##numpy array, and proceed if they are
    if type(v_pol_frac_inds) == np.ndarray:
        components.v_pol_frac_orig_inds = np.intersect1d(v_pol_frac_inds, all_orig_inds)
        components.num_v_pol_fracs = len(components.v_pol_frac_orig_inds)
        
    if type(v_power_inds) == np.ndarray:
        components.v_power_orig_inds = np.intersect1d(v_power_inds, all_orig_inds)    
        components.num_v_powers = len(components.v_power_orig_inds)
        
    if type(v_curve_inds) == np.ndarray:
        components.v_curve_orig_inds = np.intersect1d(v_curve_inds, all_orig_inds)
        components.num_v_curves = len(components.v_curve_orig_inds)
        
    if type(v_list_inds) == np.ndarray:
        components.v_list_orig_inds = np.intersect1d(v_list_inds, all_orig_inds)
        components.num_v_lists = len(components.v_list_orig_inds)
        
        ##Each list-type flux entry can have a different number of entires
        ##as sometimes you have NaNs in a column etc. Count the total here as
        ##we need it for mallocing later
        list_inds = np.where(np.isin(cropped_comp_counter.orig_comp_indexes, components.v_list_orig_inds) == True)[0]
        components.total_num_v_flux_entires = np.sum(cropped_comp_counter.num_v_list_fluxes[list_inds])
        
        
    if type(lin_pol_frac_inds) == np.ndarray:
        components.lin_pol_frac_orig_inds = np.intersect1d(lin_pol_frac_inds, all_orig_inds)
        components.num_lin_pol_fracs = len(components.lin_pol_frac_orig_inds)
        
    if type(lin_power_inds) == np.ndarray:
        components.lin_power_orig_inds = np.intersect1d(lin_power_inds, all_orig_inds)
        components.num_lin_powers = len(components.lin_power_orig_inds)
        
    if type(lin_curve_inds) == np.ndarray:
        components.lin_curve_orig_inds = np.intersect1d(lin_curve_inds, all_orig_inds)
        components.num_lin_curves = len(components.lin_curve_orig_inds)
        
    if type(lin_list_inds) == np.ndarray:
        components.lin_list_orig_inds = np.intersect1d(lin_list_inds, all_orig_inds)
        components.num_lin_lists = len(components.lin_list_orig_inds)
        
        list_inds = np.where(np.isin(cropped_comp_counter.orig_comp_indexes, components.lin_list_orig_inds) == True)[0]
        components.total_num_q_flux_entires = np.sum(cropped_comp_counter.num_q_list_fluxes[list_inds])
        components.total_num_u_flux_entires = np.sum(cropped_comp_counter.num_u_list_fluxes[list_inds])
        
        
    if type(lin_p_list_inds) == np.ndarray:
        components.lin_p_list_orig_inds = np.intersect1d(lin_p_list_inds, all_orig_inds)
        components.num_lin_p_lists = len(components.lin_p_list_orig_inds)
        
        list_inds = np.where(np.isin(cropped_comp_counter.orig_comp_indexes, components.lin_p_list_orig_inds) == True)[0]
        components.total_num_lin_p_flux_entires = np.sum(cropped_comp_counter.num_lin_p_list_fluxes[list_inds])
        
    
    ##Don't count num_lin_lists in the angles, because that model doesn't
    ##use the rotation measure
    components.num_lin_angles = components.num_lin_pol_fracs + components.num_lin_powers \
                              + components.num_lin_curves + components.num_lin_p_lists
    
    return chunk_map

def fill_chunk_component(comp_type : CompTypes,
                         cropped_comp_counter : Component_Type_Counter,
                         power_iter: int, num_chunk_power : int,
                         curve_iter: int, num_chunk_curve : int,
                         list_iter: int, num_chunk_list : int) -> Skymodel_Chunk_Map:
    """
    Fills in the relevant fields inside a `Skymodel_Chunk_Map` based on the given component type and
    the number of components of each type that are required in the chunk.
    The function returns a `Skymodel_Chunk_Map` in prepartion to read a number
    of chunks from the skymodel
    
    Parameters
    -------------
    comp_type : CompTypes
        The type of component to be filled in the chunk.
    cropped_comp_counter : Component_Type_Counter
        The counter object that contains the indices of the cropped components.
    power_iter : int
        The starting index of the power components in the cropped component counter.
    num_chunk_power : int
        The number of power components to be included in the chunk.
    curve_iter : int
        The starting index of the curve components in the cropped component counter.
    num_chunk_curve : int
        The number of curve components to be included in the chunk.
    list_iter : int
        The starting index of the list components in the cropped component counter.
    num_chunk_list : int
        The number of list components to be included in the chunk.

    Returns
    --------
    chunk_map: Skymodel_Chunk_Map
        A `Skymodel_Chunk_Map` object that contains the filled-in `Components_Map`.
    """
    
    # 
    
    if comp_type == CompTypes.POINT:
    
        chunk_map = Skymodel_Chunk_Map(n_point_powers = num_chunk_power,
                                       n_point_curves = num_chunk_curve,
                                       n_point_lists = num_chunk_list)
        
        components = chunk_map.point_components
        power_inds = cropped_comp_counter.orig_point_power_inds
        curve_inds = cropped_comp_counter.orig_point_curve_inds
        list_inds = cropped_comp_counter.orig_point_list_inds
        ##Need this to count the total number of list type flux entries
        cropped_list_inds = np.where(cropped_comp_counter.comp_types == CompTypes.POINT_LIST.value)[0]
        
    elif comp_type == CompTypes.GAUSSIAN:
        
        chunk_map = Skymodel_Chunk_Map(n_gauss_powers = num_chunk_power,
                                       n_gauss_curves = num_chunk_curve,
                                       n_gauss_lists = num_chunk_list)
        
        components = chunk_map.gauss_components
        power_inds = cropped_comp_counter.orig_gauss_power_inds
        curve_inds = cropped_comp_counter.orig_gauss_curve_inds
        list_inds = cropped_comp_counter.orig_gauss_list_inds
        ##Need this to count the total number of list type flux entries
        cropped_list_inds = np.where(cropped_comp_counter.comp_types == CompTypes.GAUSS_LIST.value)[0]
        
    ##chop everything down to what we want in this chunk
    components.power_orig_inds = power_inds[power_iter:power_iter+num_chunk_power]
    components.curve_orig_inds = curve_inds[curve_iter:curve_iter+num_chunk_curve]
    components.list_orig_inds = list_inds[list_iter:list_iter+num_chunk_list]
    
    # do this one to count how many flux entries there are in total
    cropped_list_inds = cropped_list_inds[list_iter:list_iter+num_chunk_list]
    components.total_num_flux_entires = np.sum(cropped_comp_counter.num_list_fluxes[cropped_list_inds])
    
    ##Now add in any polarisation mapping as needed
    fill_chunk_map_polarised_info(comp_type, chunk_map, cropped_comp_counter)
    
    # print(components.num_v_pol_fracs,components.num_v_powers,components.num_v_curves,components.num_lin_pol_fracs,components.num_lin_powers,components.num_lin_curves)
    
    return chunk_map
    
def map_chunk_pointgauss(cropped_comp_counter : Component_Type_Counter,
                         chunk_ind : int, comps_per_chunk : int,
                         point_source = False, gaussian_source = False) -> Components_Map:
    """
    For a given chunk index `chunk_ind`, work out how many of each
    type of COMPONENT to fit in the chunk, and then map specifics across
    to that one chunk

    Parameters
    -----------
    cropped_comp_counter : Component_Type_Counter
        A cropped component counter object that contains information about the components to use in the chunk.
    chunk_ind : int
        The index of the chunk.
    comps_per_chunk : int
        The number of components per chunk.
    point_source : bool
        Whether to use point sources. Default is False.
    gaussian_source : bool
        Whether to use gaussian sources. Default is False.

    Returns
    ---------
    chunk_map : Components_Map
        A Components_Map object that contains mapping information about the components in the chunk.
    """
    
    if not point_source and not gaussian_source:
        print("You must set either `point_source` or `gaussian_source` to True")
        return 1
    
    elif point_source and gaussian_source:
        print("You must set one `point_source` or `gaussian_source` to True, " 
              "but both are set to True")
        return 1
    else:
    
        ##Splitting POINTs and GAUSSIANS into lovely chunks that our GPU can chew
        ##First we have to ascertain where in the chunking we are, and which type
        ##of component we have to include

        ##Lower and upper indexes of components covered in this chunk
        lower_comp_ind = chunk_ind * comps_per_chunk
        upper_comp_ind = (chunk_ind + 1) * comps_per_chunk

        ##These ints are used to do pointer arithmatic to grab the correct portions
        ##of arrays out of `cropped_src` and into `temp_cropped_src`
        power_iter = 0
        curve_iter = 0
        list_iter = 0
        
        num_chunk_power = 0
        num_chunk_curve = 0
        num_chunk_list = 0
        
        n_powers = 0
        n_curves = 0
        n_lists = 0

        if point_source:
            n_powers = cropped_comp_counter.num_point_flux_powers
            n_curves = cropped_comp_counter.num_point_flux_curves
            n_lists = cropped_comp_counter.num_point_flux_lists
            
        elif gaussian_source:
            n_powers = cropped_comp_counter.num_gauss_flux_powers
            n_curves = cropped_comp_counter.num_gauss_flux_curves
            n_lists = cropped_comp_counter.num_gauss_flux_lists
            
        ##Given the information about either point of gaussian components,
        ##spit out numbers of where in all the components we have reached
        ##for this chunk (*_iter), and how many of them there are
        ##(num_chunk*)
        
        power_iter, curve_iter, list_iter, num_chunk_power, num_chunk_curve, num_chunk_list = increment_flux_type_counters(power_iter, curve_iter, list_iter, num_chunk_power, num_chunk_curve, num_chunk_list, n_powers, n_curves, n_lists, comps_per_chunk, lower_comp_ind, upper_comp_ind)
        
        ##using these numbers, populate the chunk_comp_counter, which is
        ##our map for reading in the full information from the sky model for
        ##this chunk
        if point_source:
            chunk_map = fill_chunk_component(CompTypes.POINT,
                                             cropped_comp_counter,
                                             power_iter, num_chunk_power,
                                             curve_iter, num_chunk_curve,
                                             list_iter, num_chunk_list)
            
        elif gaussian_source:
            chunk_map = fill_chunk_component(CompTypes.GAUSSIAN,
                                             cropped_comp_counter,
                                             power_iter, num_chunk_power,
                                             curve_iter, num_chunk_curve,
                                             list_iter, num_chunk_list)
            
        chunk_map.has_intr_pol_angle = cropped_comp_counter.has_intr_pol_angle
    
        return chunk_map


def create_shape_basis_maps(cropped_comp_counter : Component_Type_Counter):
    """
    Creates maps that associate each shape basis function with its corresponding original component index, 
    component type, and parameter index.
    
    Parameters:
    -----------
    cropped_comp_counter : Component_Type_Counter
        An instance of the Component_Type_Counter class that contains information about the number of shape 
        coefficients for each component type and the corresponding indices of the components in the original 
        component list.
        
    Returns:
    --------
    shape_basis_to_orig_comp_index_map : numpy.ndarray
        An array that maps each shape basis function to its corresponding original component index.
    shape_basis_to_comp_type_map : numpy.ndarray
        An array that maps each shape basis function to its corresponding component type.
    shape_basis_param_index : numpy.ndarray
        An array that maps each shape basis function to its corresponding parameter index within its component.
    """
    
    shape_basis_to_orig_comp_index_map = np.empty(cropped_comp_counter.total_shape_basis)
    shape_basis_to_comp_type_map = np.empty(cropped_comp_counter.total_shape_basis)
    
    ##this holds the index of each basis function within a shapelet component
    shape_basis_param_index = np.empty(cropped_comp_counter.total_shape_basis)
    
    coeff_iter = 0
    
    ##Go through power law flux indexes first, then curved, then list
    for comp_ind, orig_comp_ind in zip(cropped_comp_counter.shape_power_inds, cropped_comp_counter.orig_shape_power_inds):
        
        num_basis = cropped_comp_counter.num_shape_coeffs[comp_ind]
        shape_basis_to_orig_comp_index_map[coeff_iter:coeff_iter+num_basis] = orig_comp_ind
        shape_basis_to_comp_type_map[coeff_iter:coeff_iter+num_basis] = CompTypes.SHAPE_POWER.value
        shape_basis_param_index[coeff_iter:coeff_iter+num_basis] = np.arange(num_basis)
        coeff_iter += num_basis
        
    for comp_ind, orig_comp_ind in zip(cropped_comp_counter.shape_curve_inds, cropped_comp_counter.orig_shape_curve_inds):
        
        num_basis = cropped_comp_counter.num_shape_coeffs[comp_ind]
        shape_basis_to_orig_comp_index_map[coeff_iter:coeff_iter+num_basis] = orig_comp_ind
        shape_basis_to_comp_type_map[coeff_iter:coeff_iter+num_basis] = CompTypes.SHAPE_CURVE.value
        shape_basis_param_index[coeff_iter:coeff_iter+num_basis] = np.arange(num_basis)
        coeff_iter += num_basis
        
    for comp_ind, orig_comp_ind in zip(cropped_comp_counter.shape_list_inds, cropped_comp_counter.orig_shape_list_inds):
        
        num_basis = cropped_comp_counter.num_shape_coeffs[comp_ind]
        shape_basis_to_orig_comp_index_map[coeff_iter:coeff_iter+num_basis] = orig_comp_ind
        shape_basis_to_comp_type_map[coeff_iter:coeff_iter+num_basis] = CompTypes.SHAPE_LIST.value
        shape_basis_param_index[coeff_iter:coeff_iter+num_basis] = np.arange(num_basis)
        coeff_iter += num_basis
        
    return shape_basis_to_orig_comp_index_map, shape_basis_to_comp_type_map, shape_basis_param_index


def map_chunk_shapelets(cropped_comp_counter : Component_Type_Counter,
                        shape_basis_to_orig_comp_index_map : np.ndarray,
                        shape_basis_to_orig_type_map : np.ndarray,
                        shape_basis_param_index : np.ndarray,
                        chunk_ind : int,
                        coeffs_per_chunk : int):
    """
    Maps the shapelet components in a chunk of the sky model to their corresponding
    indices in the original sky model. This function is used to create a mapping
    between the shapelet components in the cropped sky model and their corresponding
    components in the original sky model. This mapping is used to extract the correct
    shapelet coefficients from the original sky model when creating a chunked sky model.

    Parameters
    -----------
    cropped_comp_counter : Component_Type_Counter
        A Component_Type_Counter object containing information about the components
        in the cropped sky model.
    shape_basis_to_orig_comp_index_map : np.ndarray
        An array mapping the indices of the shapelet basis functions in the cropped
        sky model to their corresponding indices in the original sky model.
    shape_basis_to_orig_type_map : np.ndarray
        An array mapping the indices of the shapelet basis functions in the cropped
        sky model to their corresponding component types in the original sky model.
    shape_basis_param_index : np.ndarray
        An array containing the indices of the shapelet basis functions in the
        original sky model.
    chunk_ind : int
        The index of the chunk being mapped.
    coeffs_per_chunk : int
        The number of shapelet coefficients in each chunk.

    Returns
    --------
    None
    """
    
    ##Upper indexes of components covered in this chunk
    upper_coeff_ind = (chunk_ind + 1) * coeffs_per_chunk

    ##These ints are used to do pointer arithmatic to grab the correct portions
    ##of arrays out of `cropped_src` and into `temp_cropped_src`
    lower_coeff_ind = chunk_ind * coeffs_per_chunk

    ##If there are enough coeffs to fill the chunk?
    if (cropped_comp_counter.total_shape_basis >= upper_coeff_ind):
        n_shape_coeffs = coeffs_per_chunk
        
    else:
        n_shape_coeffs = cropped_comp_counter.total_shape_basis % coeffs_per_chunk
        
    ##the ranges of comp types being sampled depends on which basis function
    ##coeffs we are sampling, so work out that range from the mapping arrays
    orig_index_chunk = shape_basis_to_orig_comp_index_map[lower_coeff_ind:upper_coeff_ind]
    orig_type_chunk = shape_basis_to_orig_type_map[lower_coeff_ind:upper_coeff_ind]
    
    shape_basis_param_index_chunk = shape_basis_param_index[lower_coeff_ind:upper_coeff_ind]
    
    ##cop that for an annoyingly complicated piece of logic
    ##this selects the subset of original component indexes that we want
    power_orig_inds = np.unique(orig_index_chunk[orig_type_chunk == CompTypes.SHAPE_POWER.value]).astype(int)
    curve_orig_inds = np.unique(orig_index_chunk[orig_type_chunk == CompTypes.SHAPE_CURVE.value]).astype(int)
    list_orig_inds = np.unique(orig_index_chunk[orig_type_chunk == CompTypes.SHAPE_LIST.value]).astype(int)
    
    power_shape_orig_inds = orig_index_chunk[orig_type_chunk == CompTypes.SHAPE_POWER.value]
    curve_shape_orig_inds = orig_index_chunk[orig_type_chunk == CompTypes.SHAPE_CURVE.value]
    list_shape_orig_inds = orig_index_chunk[orig_type_chunk == CompTypes.SHAPE_LIST.value]
    
    power_shape_basis_inds = shape_basis_param_index_chunk[orig_type_chunk == CompTypes.SHAPE_POWER.value]
    curve_shape_basis_inds = shape_basis_param_index_chunk[orig_type_chunk == CompTypes.SHAPE_CURVE.value]
    list_shape_basis_inds = shape_basis_param_index_chunk[orig_type_chunk == CompTypes.SHAPE_LIST.value]
    
    num_chunk_power = len(power_orig_inds)
    num_chunk_curve = len(curve_orig_inds)
    num_chunk_list = len(list_orig_inds)
    
    chunk_map = Skymodel_Chunk_Map(n_shape_powers = num_chunk_power,
                                   n_shape_curves = num_chunk_curve,
                                   n_shape_lists = num_chunk_list,
                                   n_shape_coeffs = n_shape_coeffs)
    
    ##shorthand so we're not typing as many things
    components = chunk_map.shape_components
    
    ##need some way to know what shapelet basis function indexes we
    ##want; similar to orig_comp_ind but for the basis functions
    components.power_shape_orig_inds = power_shape_orig_inds
    components.curve_shape_orig_inds = curve_shape_orig_inds
    components.list_shape_orig_inds = list_shape_orig_inds
    components.power_shape_basis_inds = power_shape_basis_inds
    components.curve_shape_basis_inds = curve_shape_basis_inds
    components.list_shape_basis_inds = list_shape_basis_inds
    
    ##Indexes of the shapelet components in the original sky model
    components.power_orig_inds = power_orig_inds
    components.curve_orig_inds = curve_orig_inds
    components.list_orig_inds = list_orig_inds
    
    ##how many shapelet coeffs we have
    components.total_shape_coeffs = n_shape_coeffs
    
    
    ##these are the indexes of each included component, within the cropped
    ##sky model itself
    cropped_list_inds = np.where(np.isin(cropped_comp_counter.orig_comp_indexes, list_orig_inds) == True)[0]
    
    ##if we have list type fluxes, count have many entries in total there are
    if num_chunk_list > 0:
        ##how many flux list entries in total are shared by these components
        
        components.total_num_flux_entires = np.sum(cropped_comp_counter.num_list_fluxes[cropped_list_inds])
    
    ##chuck in any polarisation information if needed
    fill_chunk_map_polarised_info(CompTypes.SHAPELET, chunk_map, cropped_comp_counter)
    
    chunk_map.has_intr_pol_angle = cropped_comp_counter.has_intr_pol_angle
    
    return chunk_map


def create_skymodel_chunk_map(comp_counter : Component_Type_Counter,
                              max_num_visibilities : int, num_baselines : int,
                              num_freqs : int, num_time_steps : int,
                              text_file=False) -> list:
                              
    """
    Given all the information in `comp_counter`, make a map of how to split
    the whole sky model up into manageable chunks to fit in memory. The
    purpose of this function is to record what to 'malloc' in each
    `Components_t` and `Source_t` ctype class before we lazy-load all the 
    values into them directly from the skymodel.

    Parameters
    ----------
    comp_counter: Component_Type_Counter 
         object that contains information about the number of components of each type in the sky model.
    max_num_visibilities: int
        The maximum number of visibilities that can be loaded into memory at once.
    num_baselines: int
        The number of baselines in the observation.
    num_freqs: int
        The number of frequency channels in the observation.
    num_time_steps: int
        The number of time steps in the observation.
    text_file: Boolean
        A boolean flag indicating whether to we are reading in from text file
        or not (default False)
    

    Returns
    -------
    list:
        A list of dictionaries containing information about the chunked sky model.
    """
        
    ##The number of components per chunk is set by how many visibilities
    ##we have
    comps_per_chunk = int(np.floor(max_num_visibilities / (num_baselines * num_freqs * num_time_steps)))
    
    ##pray this never happens, probably means we're going to run out of
    ##GPU memory TODO don't pray, submit a warning?
    if comps_per_chunk < 1: comps_per_chunk = 1
    
    ##chunks numbers for each type of component
    num_point_chunks = int(np.ceil(comp_counter.total_point_comps / comps_per_chunk))
    num_gauss_chunks = int(np.ceil(comp_counter.total_gauss_comps / comps_per_chunk))
    
    ##we split SHAPELET by the basis components (number of coeffs)
    num_coeff_chunks = int(np.ceil(comp_counter.total_shape_basis / comps_per_chunk))
    
    ##total number of chunks the sky model is splitted into
    num_chunks = num_point_chunks + num_gauss_chunks + num_coeff_chunks

    ##TODO maybe more efficient to set an array and shove in
    ##elements rather than appending?
    chunked_skymodel_maps = []
    
    ##Go through the point sources and add chunked maps
    for chunk_ind in range(num_point_chunks):
        chunk_map = map_chunk_pointgauss(comp_counter, chunk_ind,
                                         comps_per_chunk,
                                         point_source = True)
        chunked_skymodel_maps.append(chunk_map)
    
    ##Go through the gaussian sources and add chunked maps
    for chunk_ind in range(num_gauss_chunks):
        chunk_map = map_chunk_pointgauss(comp_counter, chunk_ind,
                                         comps_per_chunk,
                                         gaussian_source = True)
        chunked_skymodel_maps.append(chunk_map)
        
    ##need some extra mapping arrays to be able to grab the SHAPELET component
    ##that matches each basis function
    shape_basis_to_orig_comp_index_map, shape_basis_to_comp_type_map, shape_basis_param_index = create_shape_basis_maps(comp_counter)
        
    for chunk_ind in range(num_coeff_chunks):
        chunk_map = map_chunk_shapelets(comp_counter,
                                        shape_basis_to_orig_comp_index_map,
                                        shape_basis_to_comp_type_map,
                                        shape_basis_param_index,
                                        chunk_ind, comps_per_chunk)
        
        chunked_skymodel_maps.append(chunk_map)
        
    print(f"After chunking there are {len(chunked_skymodel_maps)} chunks")
        
    return chunked_skymodel_maps