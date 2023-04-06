import numpy as np
import erfa
from enum import Enum, auto
from typing import Union
import os

D2R = np.pi / 180.0

class CompTypes(Enum):
    """That's right, C-style enum inside Python
    This Class let's us label the component/flux type combinations with a unique
    name, but as it's an enum each label only takes 8 bytes of memory, so we can
    stack loads of them into an array. We can also do numpy operations
    on them like np.where"""
    
    ##broad difference between component type
    POINT = auto()
    GAUSSIAN = auto()
    SHAPELET = auto()
    
    ##finer difference component type + flux type
    POINT_POWER = auto()
    POINT_CURVE = auto()
    POINT_LIST = auto()
    GAUSS_POWER = auto()
    GAUSS_CURVE = auto()
    GAUSS_LIST = auto()
    SHAPE_POWER = auto()
    SHAPE_CURVE = auto()
    SHAPE_LIST = auto()
    
class Component_Type_Counter():
    """Holds counts for the various types of components in a source list,
    including type (point, gaaussian, shapelet) and flux model (power law,
    curved power law, list). Contains methods to count properties of current
    source being read in, and then to add that component"""
    def __init__(self):
        """Setup all the parameters and intialise them to 0"""
        
        
        ##Things we need to collect to be able to do a sky crop and then
        ##malloc chunks of sky models. Start off with a certain size array
        ##and make things bigger when we need it - appending to lists is
        ##memory scary
        self.array_size = 100000
        # self.basis_array_size = 1000
        
        ##Use source_indexes, comp_types, and file_line_nums for tests against 
        ##an index, so initialise with -1, as anything real will be assigned 0 upward
        
        ##index of SOURCE that each COMPONENT belongs to
        self.source_indexes = np.full(self.array_size, -1, dtype=int)
        
        ##what type of COMPONENT - given by valyes in CompTypes()
        self.comp_types = np.full(self.array_size, -1, dtype=int)
        
        ##what number line in original sky model file each component starts
        ##at. Let's us start reading the file at a given point when reading in
        ##chunks of data; makes things faster
        self.file_line_nums = np.full(self.array_size, -1, dtype=int)
        
        ##these are used to calculate az/za, so use np.nan to make sure
        ##we don't calculate on anything we don't want to
        self.comp_ras = np.full(self.array_size, np.nan, dtype=np.float64)
        self.comp_decs = np.full(self.array_size, np.nan, dtype=np.float64)
        
        ##These are only filled if necessary so just use zeros
        self.num_list_fluxes = np.zeros(self.array_size, dtype=int)
        self.num_shape_coeffs = np.zeros(self.array_size, dtype=int)
        
        # ##OK, this will be used to record which basis function matches
        # ##what component
        # self.shape_basis_param_index = np.full(self.basis_array_size, np.nan, dtype=np.float64)
        
        ##used to index components we read things in
        self.comp_index = -1
        
        ##used to index the current source as we read things in
        self.source_index = -1
        
        ##used to keep track of what line we've got to in the file we're reading
        self.file_line_num = -1
        
        ##store current ra and dec
        self.ra = np.nan
        self.dec = np.nan

        ##used to set comp type during reading
        self.point = 0
        self.gaussian = 0
        self.shapelet = 0

        ##used to select flux type during reading
        self.flux_power = 0
        self.flux_curve = 0
        self.flux_list = 0
        
        ##used to count list fluxes during reading
        self.num_list_flux = 0
        
        #used to count shapelet basis functions
        self.num_shape_basis = 0

        ##used later during cropping for below the horizon        
        self.include_flags = False
        
        ## user later to reference the original component index in the mother
        ## sky model when chunking; need this index to pull out the correct
        ## component information from the full list
        self.orig_comp_indexes = None

        ##used later on to nail down exactly where each type of component
        ##can be indexed from the parent sky mode        
        self.point_power_inds = False
        self.point_curve_inds = False
        self.point_list_inds = False
        self.gauss_power_inds = False
        self.gauss_curve_inds = False
        self.gauss_list_inds = False
        self.shape_power_inds = False
        self.shape_curve_inds = False
        self.shape_list_inds = False
        
    def new_file_line(self):
        """Iterates the line in the file that we are on"""
        self.file_line_num += 1
        
    def new_source(self):
        """Iterates the source index if we are in a new source"""
        self.source_index += 1
        
    def count_list_flux(self):
        """Count an entry to the current number of flux entries for a list
        type flux, if the component has a list type flux model"""

        if self.flux_list:
            self.num_list_flux += 1
            
    def initiate_component(self):
        """Starting a new component, so need to iterate the component index
        and capture the current source index. Make arrays bigger if
        necessary"""
        
        ##make arrays bigger if we need more space
        if self.array_size <= self.comp_index + 1:
            self.source_indexes = np.concatenate((self.source_indexes,
                                    np.full(self.array_size, -1, dtype=int)))
            self.comp_types = np.concatenate((self.comp_types,
                                    np.full(self.array_size, -1, dtype=int)))
            self.file_line_nums = np.concatenate((self.file_line_nums,
                                    np.full(self.array_size, -1, dtype=int)))
            
            self.comp_ras = np.concatenate((self.comp_ras,
                                    np.full(self.array_size, np.nan, dtype=np.float64)))
            self.comp_decs = np.concatenate((self.comp_decs,
                                    np.full(self.array_size, np.nan, dtype=np.float64)))
            
            self.num_list_fluxes = np.concatenate((self.num_list_fluxes,
                                    np.zeros(self.array_size, dtype=int)))
            self.num_shape_coeffs = np.concatenate((self.num_shape_coeffs,
                                    np.zeros(self.array_size, dtype=int)))
            
            ##array size has now doubled
            self.array_size *= 2
            
        ##iterate the component index up by one
        self.comp_index += 1
        ##capture the current source index    
        self.source_indexes[self.comp_index] = self.source_index
        
        self.file_line_nums[self.comp_index] = self.file_line_num
        

    def add_component_reset_counters(self):
        """Add the current component to overall counts, and reset the current
        component values for the next component"""
        
        self.comp_ras[self.comp_index] = self.ra
        self.comp_decs[self.comp_index] = self.dec
        
        if self.point:
            if self.flux_power:
                self.comp_types[self.comp_index] = CompTypes.POINT_POWER.value
            elif self.flux_curve:
                self.comp_types[self.comp_index] = CompTypes.POINT_CURVE.value
            elif self.flux_list:
                self.comp_types[self.comp_index] = CompTypes.POINT_LIST.value

        elif self.gaussian:
            if self.flux_power:
                self.comp_types[self.comp_index] = CompTypes.GAUSS_POWER.value
            elif self.flux_curve:
                self.comp_types[self.comp_index] = CompTypes.GAUSS_CURVE.value
            elif self.flux_list:
                self.comp_types[self.comp_index] = CompTypes.GAUSS_LIST.value

        elif self.shapelet:
            if self.flux_power:
                self.comp_types[self.comp_index] = CompTypes.SHAPE_POWER.value
            elif self.flux_curve:
                self.comp_types[self.comp_index] = CompTypes.SHAPE_CURVE.value
            elif self.flux_list:
                self.comp_types[self.comp_index] = CompTypes.SHAPE_LIST.value
                
        if self.flux_list:
            self.num_list_fluxes[self.comp_index] = self.num_list_flux
            
        if self.shapelet:
            self.num_shape_coeffs[self.comp_index] = self.num_shape_basis


        ##reset the things
        self.ra = np.nan
        self.dec = np.nan
        self.point = 0
        self.gaussian = 0
        self.shapelet = 0
        self.flux_power = 0
        self.flux_curve = 0
        self.flux_list = 0
        self.num_list_flux = 0
        self.num_shape_basis = 0
        
    def remove_array_padding(self):
        """Once everything is read in, we can shrink down the arrays to however
        much information was read in (we keep doubling array size as needed)"""
        
        include = self.source_indexes > -1
        
        self.source_indexes = self.source_indexes[include]
        self.comp_types = self.comp_types[include]
        self.file_line_nums = self.file_line_nums[include]
        self.comp_ras = self.comp_ras[include]
        self.comp_decs = self.comp_decs[include]
        self.num_list_fluxes = self.num_list_fluxes[include]
        self.num_shape_coeffs = self.num_shape_coeffs[include]
        
    def total_components(self):
        """Create final total counts once all information has been read in"""
        # self.num_sources = self.source_index + 1
        
        self.num_sources = len(np.unique(self.source_indexes))
        
        if not type(self.orig_comp_indexes) == np.ndarray:
            self.orig_comp_indexes = np.arange(len(self.comp_types), dtype=int)
            
        ##Count point source related things
        ##Grab the indexes of all the point source types. Used later when
        self.point_power_inds = np.where(self.comp_types == CompTypes.POINT_POWER.value)[0]
        self.point_curve_inds = np.where(self.comp_types == CompTypes.POINT_CURVE.value)[0]
        self.point_list_inds = np.where(self.comp_types == CompTypes.POINT_LIST.value)[0]
        
        self.orig_point_power_inds = self.orig_comp_indexes[self.point_power_inds]
        self.orig_point_curve_inds = self.orig_comp_indexes[self.point_curve_inds]
        self.orig_point_list_inds = self.orig_comp_indexes[self.point_list_inds]
        
        self.num_point_flux_powers = len(self.point_power_inds)
        self.num_point_flux_curves = len(self.point_curve_inds)
        self.num_point_flux_lists = len(self.point_list_inds)
        self.total_point_comps = self.num_point_flux_powers + self.num_point_flux_curves + self.num_point_flux_lists
        
        ##Count gaussian source related things
        ##Grab the indexes of all the gaussian source types. Used later when
        self.gauss_power_inds = np.where(self.comp_types == CompTypes.GAUSS_POWER.value)[0]
        self.gauss_curve_inds = np.where(self.comp_types == CompTypes.GAUSS_CURVE.value)[0]
        self.gauss_list_inds = np.where(self.comp_types == CompTypes.GAUSS_LIST.value)[0]
        
        self.orig_gauss_power_inds = self.orig_comp_indexes[self.gauss_power_inds]
        self.orig_gauss_curve_inds = self.orig_comp_indexes[self.gauss_curve_inds]
        self.orig_gauss_list_inds = self.orig_comp_indexes[self.gauss_list_inds]
        
        self.num_gauss_flux_powers = len(self.gauss_power_inds)
        self.num_gauss_flux_curves = len(self.gauss_curve_inds)
        self.num_gauss_flux_lists = len(self.gauss_list_inds)
        self.total_gauss_comps = self.num_gauss_flux_powers + self.num_gauss_flux_curves + self.num_gauss_flux_lists
        
        ##Count shaplet related things
        ##Grab the indexes of all the shape source types. Used later when
        self.shape_power_inds = np.where(self.comp_types == CompTypes.SHAPE_POWER.value)[0]
        self.shape_curve_inds = np.where(self.comp_types == CompTypes.SHAPE_CURVE.value)[0]
        self.shape_list_inds = np.where(self.comp_types == CompTypes.SHAPE_LIST.value)[0]
        
        self.orig_shape_power_inds = self.orig_comp_indexes[self.shape_power_inds]
        self.orig_shape_curve_inds = self.orig_comp_indexes[self.shape_curve_inds]
        self.orig_shape_list_inds = self.orig_comp_indexes[self.shape_list_inds]
        
        self.num_shape_flux_powers = len(self.shape_power_inds)
        self.num_shape_flux_curves = len(self.shape_curve_inds)
        self.num_shape_flux_lists = len(self.shape_list_inds)
        self.total_shape_comps = self.num_shape_flux_powers + self.num_shape_flux_curves + self.num_shape_flux_lists
        
        self.total_shape_basis = np.sum(self.num_shape_coeffs)
        
        self.total_comps =  self.total_point_comps + self.total_gauss_comps + self.total_shape_comps
        
        ##TODO can we do delete self.comp_types ? As in
        ##del self.comp_types
        
    def print_info(self):
        
        
        print("total_point_comps", self.total_point_comps)
        print("\tnum_point_flux_powers", self.num_point_flux_powers)
        print("\tnum_point_flux_curves", self.num_point_flux_curves)
        print("\tnum_point_flux_lists", self.num_point_flux_lists)
        
        print("total_gauss_comps", self.total_gauss_comps)
        print("\tnum_gauss_flux_powers", self.num_gauss_flux_powers)
        print("\tnum_gauss_flux_curves", self.num_gauss_flux_curves)
        print("\tnum_gauss_flux_lists", self.num_gauss_flux_lists)
        
        print("total_shape_comps", self.total_shape_comps)
        print("\tnum_shape_flux_powers", self.num_shape_flux_powers)
        print("\tnum_shape_flux_curves", self.num_shape_flux_curves)
        print("\tnum_shape_flux_lists", self.num_shape_flux_lists)
        print("\ttotal_shape_basis", self.total_shape_basis)
        
# @profile
def crop_below_horizon(lst : float, latitude : float,
                       comp_counter : Component_Type_Counter, 
                       crop_by_component=True) -> np.ndarray:
    """Crop anything below the horizon - just work out cropping flags"""
    
    ##Want numpy arrays to do maths with
    if type(comp_counter.comp_ras) == list:
        comp_counter.comp_ras = np.array(comp_counter.comp_ras)
        
    if type(comp_counter.comp_decs) == list:
        comp_counter.comp_decs = np.array(comp_counter.comp_decs)
    
    if type(comp_counter.source_indexes) == list:
        comp_counter.source_indexes = np.array(comp_counter.source_indexes)
    
    ##Calculate lst, and then azimuth/elevation
    comp_has = lst - comp_counter.comp_ras
    comp_azs, comp_els = erfa.hd2ae(comp_has, comp_counter.comp_decs,
                                    latitude)
    
    ##Crop below the horizon based on COMPONENT/SOURCE
    include_flags = np.zeros(comp_counter.total_comps)
    above_horizon = np.where(comp_els >= 0)[0]
    
    print(f"Have read in {comp_counter.total_comps} components")
    
    ##If crop by COMPONENT, just include everything above the horizon
    if crop_by_component:
        include_flags[above_horizon] = 1.0
    else:
        ##TODO this is super slow, sort that shit out
        for source_ind in range(comp_counter.num_sources):
            ##These are the indexes of all COMPONENTS that belong to this SOURCE
            component_indexes = np.where(comp_counter.source_indexes == source_ind)
            
            ##If all the COMPONENTs for this SOURCE are above horizon,
            ##add the include flag for all of them
            if (comp_els[component_indexes] >= 0).all():
                include_flags[component_indexes] = 1.0
      
    ##original component indexes in the mother sky model - needed for retrieving
    ##full details when loading in chunks of sky model          
    comp_counter.orig_comp_indexes = np.where(include_flags == 1)[0]
    
    ##Shrink down sky model to only things above the horizon
    comp_counter.comp_types = comp_counter.comp_types[comp_counter.orig_comp_indexes]
    comp_counter.comp_ras = comp_counter.comp_ras[comp_counter.orig_comp_indexes]
    comp_counter.comp_decs = comp_counter.comp_decs[comp_counter.orig_comp_indexes]
    comp_counter.num_list_fluxes = comp_counter.num_list_fluxes[comp_counter.orig_comp_indexes]
    comp_counter.num_shape_coeffs = comp_counter.num_shape_coeffs[comp_counter.orig_comp_indexes]
    comp_counter.file_line_nums = comp_counter.file_line_nums[comp_counter.orig_comp_indexes]

    ##re-total everything now we have changed the sizes
    comp_counter.total_components()

    print(f"After cropping there are {comp_counter.total_comps} components")
    
    ##free some memory explicitly as these can be big arrays
    del include_flags
    
    return comp_counter


class Component_Info(object):
    """Holds all the information for the current component when a sourcelist
    is being read in"""
    
    def __init__(self):
        """Setup all the parameters and intialise them to 0"""
        
        ##what type of COMPONENT - given by valyes in CompTypes()
        self.comp_type = -1
        
        ##store current ra and dec
        self.ra = np.nan
        self.dec = np.nan

        ##used to set comp type during reading
        self.point = 0
        self.gaussian = 0
        self.shapelet = 0

        ##used to select flux type during reading
        self.flux_power = 0
        self.flux_curve = 0
        self.flux_list = 0
        
        ##GAUSSIAN/SHAPELET type things
        self.major = np.nan
        self.minor = np.nan
        self.pa = np.nan
        
        self.n1s = []
        self.n2s = []
        self.coeffs = []
        
        ##power/curve power law type things
        self.si = np.nan
        self.curve_q = np.nan        
        
        ##frequency information. this will always be of length one for 
        ##power and curved power laws, and any length for a list-style flux
        self.freqs = []
        
        ##this will be a list of Stokes vectors so shape = (num_refs, 4)
        ##num refs is always 1 for power / curved power laws
        self.fluxes = []
        
        ##use to keep track of how many freqs have been added
        self.num_fluxes = 0
        
        #used to count shapelet basis functions
        self.num_shape_basis = 0
        
        ##parent source name; useful for debugging
        self.source_name = ''
    
    
    ##bunch of methods to get information is determining component
    ##type
    def set_point(self):
        self.point = 1
        self.gaussian = 0
        self.shapelet = 0
        
    def set_gaussian(self):
        self.point = 0
        self.gaussian = 1
        self.shapelet = 0
        
    def set_shapelet(self):
        self.point = 0
        self.gaussian = 0
        self.shapelet = 1
        
    def set_flux_power(self):
        self.flux_power = 1
        
    def set_flux_curve(self):
        self.flux_curve = 1
        
    def set_flux_list(self):
        self.flux_list = 1
        
    
    ##bunch of methods for gathering component information
    def add_ra(self, ra : float):
        self.ra = ra
        
    def add_dec(self, dec : float):
        self.dec = dec
    
    ##gaussian/shapelet things
    def add_major(self, major : float):
        self.major = major
        
    def add_minor(self, minor : float):
        self.minor = minor
        
    def add_pa(self, pa : float):
        self.pa = pa
        
    ##shapelet things
    def add_n1(self, n1 : float):
        self.n1s.append(n1)
        
    def add_n2(self, n2 : float):
        self.n2s.append(n2)
        
    def add_coeff(self, coeffs : float):
        self.coeffs.append(coeffs)
        
    ##power/curved law related things
    def add_si(self, si : float):
        self.si = si
        
    def add_curve_q(self, curve_q : float):
        self.curve_q = curve_q
        
        
    def add_ref_freq(self, freq : float):
        """Setup another frequency entry, and add a corresponding
        empty Stokes flux vector"""
        self.num_fluxes += 1
        self.freqs.append(freq)
        self.fluxes.append(np.array([np.nan, np.nan, np.nan, np.nan]))
        
    def add_stokesI(self, stokesI):
        """Uses the current self.num_fluxes to add Stokes I to the
        correct reference frequency"""
        
        self.fluxes[self.num_fluxes - 1][0] = stokesI
        
    def add_stokesQ(self, stokesQ):
        """Uses the current self.num_fluxes to add Stokes Q to the
        correct reference frequency"""
        
        self.fluxes[self.num_fluxes - 1][1] = stokesQ
        
    def add_stokesU(self, stokesU):
        """Uses the current self.num_fluxes to add Stokes U to the
        correct reference frequency"""
        
        self.fluxes[self.num_fluxes - 1][2] = stokesU
        
    def add_stokesV(self, stokesV):
        """Uses the current self.num_fluxes to add Stokes V to the
        correct reference frequency"""
        
        self.fluxes[self.num_fluxes - 1][3] = stokesV
        
    def finalise_comp(self):
        
        if self.point:
            if self.flux_power:
                self.comp_type = CompTypes.POINT_POWER
            elif self.flux_curve:
                self.comp_type = CompTypes.POINT_CURVE
            elif self.flux_list:
                self.comp_type = CompTypes.POINT_LIST
        elif self.gaussian:
            if self.flux_power:
                self.comp_type = CompTypes.GAUSS_POWER
            elif self.flux_curve:
                self.comp_type = CompTypes.GAUSS_CURVE
            elif self.flux_list:
                self.comp_type = CompTypes.GAUSS_LIST
        if self.shapelet:
            if self.flux_power:
                self.comp_type = CompTypes.SHAPE_POWER
            elif self.flux_curve:
                self.comp_type = CompTypes.SHAPE_CURVE
            elif self.flux_list:
                self.comp_type = CompTypes.SHAPE_LIST
                
            self.n1s = np.array(self.n1s)
            self.n2s = np.array(self.n2s)
            self.coeffs = np.array(self.coeffs)
        
        
        self.num_fluxes = len(self.fluxes)
        self.fluxes = np.array(self.fluxes)
        
        empty_fluxes = []
        for flux_ind in range(self.num_fluxes):
            if np.all(np.isnan(self.fluxes[flux_ind])):
                empty_fluxes.append(flux_ind)
                
        ##set missing fluxes to zero
        self.fluxes[np.isnan(self.fluxes)] = 0.0
        
        return empty_fluxes