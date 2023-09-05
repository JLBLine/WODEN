
import numpy as np
import sys
import os
from enum import Enum

from wodenpy.skymodel.woden_skymodel import Component_Type_Counter, CompTypes
    
NUM_FLUX_TYPES = 3
    
class Components_Map(object):
    """Mapping information for a set of components, of either POINT,
    GAUSSIAN or SHAPELET type"""
    
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
        self.lowest_file_num = np.nan
        
        ##use this to count how many flux list entries there are in total
        ##so we can allocate correct amount when reading in full information
        self.total_num_flux_entires = 0
        
        #use this to count how many shapelet basis functions that are in total
        #so we can allocate correct amount when reading in full information
        self.total_shape_coeffs = 0
    
    
class Skymodel_Chunk_Map(object):
    """
    Something to hold all the chunked information in a format that mirrors
    `components_t` in the C code. Should make the logic of reading in chunks
    of sky model in the ctypes class easier
    """
    
    def __init__(self, n_point_powers = 0, n_point_curves = 0, n_point_lists = 0,
                       n_gauss_powers = 0, n_gauss_curves = 0, n_gauss_lists = 0,
                       n_shape_powers = 0, n_shape_curves = 0, n_shape_lists = 0,
                       n_shape_coeffs = 0):
        """Setup everything with zeros as default, and total everything up"""
        
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
        
        self.n_shape_coeffs = n_shape_coeffs
        
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
            
    def make_all_orig_inds_array(self):
        """Look through all component and flux types and consolidate into one
        array of original component indexes. Use this when reading in full
        information from the sky model"""
        
        self.all_orig_inds = np.empty(self.n_points + self.n_gauss + self.n_shape_coeffs, dtype=int)
        
        lowest_file_lines = []
        
        if self.n_points > 0:
            lowest_file_lines.append(self.point_components.lowest_file_num)
            low_ind = 0
            if self.n_point_powers > 0:
                self.all_orig_inds[low_ind:low_ind+self.n_point_powers] = self.point_components.power_orig_inds
                low_ind += self.n_point_powers
                
            if self.n_point_curves > 0:
                self.all_orig_inds[low_ind:low_ind+self.n_point_curves] = self.point_components.curve_orig_inds
                low_ind += self.n_point_curves
                
            if self.n_point_lists > 0:
                self.all_orig_inds[low_ind:low_ind+self.n_point_lists] = self.point_components.list_orig_inds
                low_ind += self.n_point_lists
                
        if self.n_gauss > 0:
            lowest_file_lines.append(self.gauss_components.lowest_file_num)
            low_ind = self.n_points
            if self.n_gauss_powers > 0:
                self.all_orig_inds[low_ind:low_ind+self.n_gauss_powers] = self.gauss_components.power_orig_inds
                low_ind += self.n_gauss_powers
                
            if self.n_gauss_curves > 0:
                self.all_orig_inds[low_ind:low_ind+self.n_gauss_curves] = self.gauss_components.curve_orig_inds
                low_ind += self.n_gauss_curves
                
            if self.n_gauss_lists > 0:
                self.all_orig_inds[low_ind:low_ind+self.n_gauss_lists] = self.gauss_components.list_orig_inds
                low_ind += self.n_gauss_lists
        
        ##The component
        
        
        if self.n_shapes > 0:
            
            lowest_file_lines.append(self.shape_components.lowest_file_num)
            low_ind = self.n_points + self.n_gauss
            if self.n_shape_powers > 0:
                
                power_indexes = self.shape_components.power_shape_orig_inds
                self.all_orig_inds[low_ind:low_ind+len(power_indexes)] = power_indexes
                low_ind += len(power_indexes)
                
            if self.n_shape_curves > 0:
                curve_indexes = self.shape_components.curve_shape_orig_inds
                self.all_orig_inds[low_ind:low_ind+len(curve_indexes)] = curve_indexes
                low_ind += len(curve_indexes)
                
            if self.n_shape_lists > 0:
                list_indexes = self.shape_components.list_shape_orig_inds
                self.all_orig_inds[low_ind:low_ind+len(list_indexes)] = list_indexes
                low_ind += len(list_indexes)
                
        self.lowest_file_number = min(lowest_file_lines)
            
    def print_info(self):
        
        
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
    """Here, given the overall lower and upper index in a given type of components,
    work out how many of each flux type we have and increment the counters as appropriate

    Always order things as POWER_LAW, CURVED_POWER_LAW, LIST"""
    

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

def find_lowest_file_line(cropped_comp_counter : Component_Type_Counter,
                          components : Components_Map,
                          cropped_power_inds : np.ndarray, cropped_curve_inds : np.ndarray,
                          cropped_list_inds : np.ndarray):
    
    if len(cropped_power_inds) > 0:
        min_power = np.nanmin(cropped_comp_counter.file_line_nums[cropped_power_inds])
    else:
        min_power = np.nan
        
    if len(cropped_curve_inds) > 0:
        min_curve = np.nanmin(cropped_comp_counter.file_line_nums[cropped_curve_inds])
    else:
        min_curve = np.nan
    
    if len(cropped_list_inds) > 0:
        min_list = np.nanmin(cropped_comp_counter.file_line_nums[cropped_list_inds])
    else:
        min_list = np.nan
        
    lowest_file_num = int(np.nanmin([min_power, min_curve, min_list]))
    
    return lowest_file_num
    
def fill_chunk_component(comp_type : CompTypes,
                         cropped_comp_counter : Component_Type_Counter,
                         power_iter: int, num_chunk_power : int,
                         curve_iter: int, num_chunk_curve : int,
                         list_iter: int, num_chunk_list : int) -> Skymodel_Chunk_Map:
    """Having worked out what numbers of which components we want in the
    chunk, fill in the relevant fields inside a `Components_Map` inside a chunk"""
    
    if comp_type == CompTypes.POINT:
    
        chunk_map = Skymodel_Chunk_Map(n_point_powers = num_chunk_power,
                                       n_point_curves = num_chunk_curve,
                                       n_point_lists = num_chunk_list)
        
        components = chunk_map.point_components
        # power_inds = cropped_comp_counter.point_power_inds
        # curve_inds = cropped_comp_counter.point_curve_inds
        # list_inds = cropped_comp_counter.point_list_inds
        power_inds = cropped_comp_counter.orig_point_power_inds
        curve_inds = cropped_comp_counter.orig_point_curve_inds
        list_inds = cropped_comp_counter.orig_point_list_inds
        
        cropped_power_inds = np.where(cropped_comp_counter.comp_types == CompTypes.POINT_POWER.value)[0]
        cropped_curve_inds = np.where(cropped_comp_counter.comp_types == CompTypes.POINT_CURVE.value)[0]
        cropped_list_inds = np.where(cropped_comp_counter.comp_types == CompTypes.POINT_LIST.value)[0]
        
    elif comp_type == CompTypes.GAUSSIAN:
        
        chunk_map = Skymodel_Chunk_Map(n_gauss_powers = num_chunk_power,
                                       n_gauss_curves = num_chunk_curve,
                                       n_gauss_lists = num_chunk_list)
        
        components = chunk_map.gauss_components
        # power_inds = cropped_comp_counter.gauss_power_inds
        # curve_inds = cropped_comp_counter.gauss_curve_inds
        # list_inds = cropped_comp_counter.gauss_list_inds
        power_inds = cropped_comp_counter.orig_gauss_power_inds
        curve_inds = cropped_comp_counter.orig_gauss_curve_inds
        list_inds = cropped_comp_counter.orig_gauss_list_inds
        
        cropped_power_inds = np.where(cropped_comp_counter.comp_types == CompTypes.GAUSS_POWER.value)[0]
        cropped_curve_inds = np.where(cropped_comp_counter.comp_types == CompTypes.GAUSS_CURVE.value)[0]
        cropped_list_inds = np.where(cropped_comp_counter.comp_types == CompTypes.GAUSS_LIST.value)[0]
        
        
    components.power_orig_inds = power_inds[power_iter:power_iter+num_chunk_power]
    components.curve_orig_inds = curve_inds[curve_iter:curve_iter+num_chunk_curve]
    components.list_orig_inds = list_inds[list_iter:list_iter+num_chunk_list]
    
    cropped_power_inds = cropped_power_inds[power_iter:power_iter+num_chunk_power]
    cropped_curve_inds = cropped_curve_inds[curve_iter:curve_iter+num_chunk_curve]
    cropped_list_inds = cropped_list_inds[list_iter:list_iter+num_chunk_list]
    
    components.total_num_flux_entires = np.sum(cropped_comp_counter.num_list_fluxes[cropped_list_inds])
    
    ##We don't have to start reading the information for this chunk until
    ##this line of the 
    components.lowest_file_num = find_lowest_file_line(cropped_comp_counter,
                                                               components,
                                                               cropped_power_inds,
                                                               cropped_curve_inds,
                                                               cropped_list_inds)
    
    min_comp_inds = []
    if len(cropped_power_inds) > 0: min_comp_inds.append(components.power_orig_inds.min())
    if len(cropped_curve_inds) > 0: min_comp_inds.append(components.curve_orig_inds.min())
    if len(cropped_list_inds) > 0: min_comp_inds.append(components.list_orig_inds.min())
    
    components.min_comp_ind = np.min(min_comp_inds)
    
    return chunk_map
    
def map_chunk_pointgauss(cropped_comp_counter : Component_Type_Counter,
                         chunk_ind : int, comps_per_chunk : int,
                         point_source = False, gaussian_source = False) -> Components_Map:
    """For a given chunk index `chunk_ind`, work out how many of each
    type of COMPONENT to fit in the chunk, and then map specifics across
    to that one chunk"""
    
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
    
        chunk_map.make_all_orig_inds_array()

        return chunk_map


def create_shape_basis_maps(cropped_comp_counter : Component_Type_Counter):

    
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
                        coeffs_per_chunk : int,
                        text_file = False):
    
    ##Upper indexes of components covered in this chunk
    upper_coeff_ind = (chunk_ind + 1) * coeffs_per_chunk

    ##These ints are used to do pointer arithmatic to grab the correct portions
    ##of arrays out of `cropped_src` and into `temp_cropped_src`
    lower_coeff_ind = chunk_ind * coeffs_per_chunk

    ##If there are enough coeffs to fill the chunk?
    if (cropped_comp_counter.total_shape_basis >= upper_coeff_ind):
        n_shape_coeffs = coeffs_per_chunk;
        
    else:
        n_shape_coeffs = cropped_comp_counter.total_shape_basis % coeffs_per_chunk
        
    ##the ranges of comp types being sampled depends on which basis function
    ##coeffs we are sampling, so work out that range from the mapping arrays
    orig_index_chunk = shape_basis_to_orig_comp_index_map[lower_coeff_ind:upper_coeff_ind]
    orig_type_chunk = shape_basis_to_orig_type_map[lower_coeff_ind:upper_coeff_ind]
    
    shape_basis_param_index_chunk = shape_basis_param_index[lower_coeff_ind:upper_coeff_ind]
    
    ##TODO get this reordered based on component type??
    
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
    
    ##TODO need some way to know what shapelet basis function indexes we
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
    
    cropped_power_inds = np.where(np.isin(cropped_comp_counter.orig_comp_indexes, power_orig_inds) == True)[0]
    cropped_curve_inds = np.where(np.isin(cropped_comp_counter.orig_comp_indexes, curve_orig_inds) == True)[0]
    cropped_list_inds = np.where(np.isin(cropped_comp_counter.orig_comp_indexes, list_orig_inds) == True)[0]
    
    ##if we have list type fluxes, count have many entries in total there are
    if num_chunk_list > 0:
        ##how many flux list entries in total are shared by these components
        
        components.total_num_flux_entires = np.sum(cropped_comp_counter.num_list_fluxes[cropped_list_inds])
    
    # print('--map_chunk_shapelets-----------------------')
    # # print(cropped_comp_counter.comp_types)
    # print(n_shape_coeffs, cropped_comp_counter.num_shape_flux_powers)
    # print(power_orig_inds, curve_orig_inds, list_orig_inds)
    # print('----------------------------------------')
    
    
    
    if text_file:
        ##lowest line we want to read from the 
        components.lowest_file_num = find_lowest_file_line(cropped_comp_counter,
                                                           components,
                                                           cropped_power_inds,
                                                           cropped_curve_inds,
                                                           cropped_list_inds)
    
    chunk_map.make_all_orig_inds_array()
    
    return chunk_map


def create_skymodel_chunk_map(comp_counter : Component_Type_Counter,
                              max_num_visibilities : int, num_baselines : int,
                              num_freqs : int, num_time_steps : int,
                              text_file=False) -> list:
    """Given all the information in `comp_counter`, make a map of how to split
    the whole sky model up into managable chunks to fit in memory. The
    purpose of this function is to record what to 'malloc' in each
    `Components_t` and `Source_t` ctype class before we lazy-load all the 
    values into them directly from the skymodel
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
                                        chunk_ind, comps_per_chunk,
                                        text_file = text_file)
        
        chunked_skymodel_maps.append(chunk_map)
        
    print(f"After chunking there are {len(chunked_skymodel_maps)} chunks")
        
    return chunked_skymodel_maps