import numpy as np
import erfa
from enum import Enum, auto
from typing import Union, Tuple
import os

D2R = np.pi / 180.0

from enum import Enum, auto

class CompTypes(Enum):
    """
    That's right, C-style enum inside Python
    This Class let's us label the component/flux type combinations with a unique
    name, but as it's an enum each label only takes 8 bytes of memory, so we can
    stack loads of them into an array. We can also do numpy operations
    on them like np.where
    
    :cvar auto() POINT: point source
    :cvar auto() GAUSSIAN: gaussian source
    :cvar auto() SHAPELET: shapelet source
    :cvar auto() POINT_POWER: point source + power law
    :cvar auto() POINT_CURVE: point source + curved power law
    :cvar auto() POINT_LIST: point source + list type law
    :cvar auto() GAUSS_POWER: gaussian source + power law
    :cvar auto() GAUSS_CURVE: gaussian source + curved power law
    :cvar auto() GAUSS_LIST: gaussian source + list type law
    :cvar auto() SHAPE_POWER: shapelet source + power law
    :cvar auto() SHAPE_CURVE: shapelet source + curved power law
    :cvar auto() SHAPE_LIST: shapelet source + list type law
    """
    
    POINT = auto()
    GAUSSIAN = auto()
    SHAPELET = auto()
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
    source being read in, and then to add that component
    
    Does a large amount of book keeping on indexes of components within the
    original sky model, so that we can pull out the correct information when
    lazy loading chunks of sky model.
    
    This is used by `read_radec_count_components` and the functions within.
    
    """
    def __init__(self, initial_size=100000):
        """Setup all the parameters and intialise them to 0"""
        
        
        ##Things we need to collect to be able to do a sky crop and then
        ##malloc chunks of sky models. Start off with a certain size array
        ##and make things bigger when we need it - appending to lists is
        ##memory scary
        self.array_size = initial_size
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
        """Print out the totals of everything
        """
        
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
                       crop_by_component=True) -> Component_Type_Counter:
    """
    Crop a `Component_Type_Counter` to only include components above the horizon,
    given the local sidereal time and latitude of the array

    Parameters
    ----------
    lst : float
        The local sidereal time (in radians).
    latitude : float
        The array latitude (in radians).
    comp_counter : Component_Type_Counter
        An object containing information about the components in the sky model.
    crop_by_component : bool, optional
        If True, crop components individually. If False, crop sources as a whole.
        Default is True.

    Returns
    -------
    comp_counter : Component_Type_Counter
        An updated version of the input `comp_counter` object, with only the
        components above the horizon included.

    Notes
    -----
    This function crops a sky model to only include components above the horizon.
    It uses the input `lst` and `latitude` to calculate the azimuth and elevation
    of each component in the sky model, and then includes only those components
    with elevation greater than or equal to zero.

    If `crop_by_component` is True, each component is cropped individually. If
    False, sources are cropped as a whole. In the latter case, all components
    belonging to a source are included if at least one of them is above the horizon.

    The input `comp_counter` object is updated in place, and the updated object
    is returned.

    """
    
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
    """
    A class that stores information about a sky model component.

    :cvar int comp_type: The type of component, given by values in CompTypes().
    :cvar float ra: The right ascension of the component in radians.
    :cvar float dec: The declination of the component in radians.
    :cvar int point: Used to set the component type during reading.
    :cvar int gaussian: Used to set the component type during reading.
    :cvar int shapelet: Used to set the component type during reading.
    :cvar int flux_power: Used to select flux type during reading.
    :cvar int flux_curve: Used to select flux type during reading.
    :cvar int flux_list: Used to select flux type during reading.
    :cvar float major: The major axis of the component in radians.
    :cvar float minor: The minor axis of the component in radians.
    :cvar float pa: The position angle of the component in radians.
    :cvar list n1s: A list of shapelet basis function n1 values.
    :cvar list n2s: A list of shapelet basis function n2 values.
    :cvar list coeffs: A list of shapelet basis function coefficients.
    :cvar float si: The spectral index of the component.
    :cvar float curve_q: The curvature of the component spectrum.
    :cvar list freqs: A list of frequencies at which the flux is defined (Hz).
    :cvar list fluxes: A list of Stokes vectors defining the flux (Jy).
    :cvar int num_fluxes: The number of frequencies that have been added.
    :cvar int num_shape_basis: The number of shapelet basis functions.
    :cvar str source_name: The name of the parent source.
    :cvar float norm_comp_pl: The reference flux of the power-law component at 200MHz (Jy)
    :cvar float norm_comp_cpl: The reference flux of the curved power-law component at 200MHz (Jy)
    """
    
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
        
        ##used when converting other skymodel formats into a FITS style
        self.norm_comp_pl = np.nan
        self.norm_comp_cpl = np.nan
    
    
def calc_pl_norm_at_200MHz(component : Component_Info) -> Component_Info:
    """The FITS style sky model references everything to 200MHz, so extrap
    a power-law model reference flux to 200MHz, and set frequency to 200MHz

    Parameters
    ----------
    component : Component_Info
        A populated Component_Info instance with power law info

    Returns
    -------
    Component_Info
        Updated component info with power law extrapolated to 200MHz
    """

    ##There are four stokes params, just extrap them all
    for ref_flux_ind in range(4):

        ref_flux = component.fluxes[0][ref_flux_ind]
        ref_freq = component.freqs[0]

        new_ref_flux = ref_flux*(200e+6 / ref_freq)**component.si
        
        component.fluxes[0][ref_flux_ind] = new_ref_flux

    component.norm_comp_pl = component.fluxes[0][0]
    component.freqs[0] = 200e+6

    return component

def calc_cpl_norm_at_200MHz(component : Component_Info) -> Component_Info:
    """The FITS style sky model references everything to 200MHz, so extrap
    a curved power-law model reference flux to 200MHz, and set frequency to 200MHz

    Parameters
    ----------
    component : Component_Info
        A populated Component_Info instance with curved power law info

    Returns
    -------
    Component_Info
        Updated component info with curved power law extrapolated to 200MHz
    """
    
    if component.freqs[0] == 200e+6:
        component.norm_comp_cpl = component.fluxes[0][0]

    else:
        ##There are four stokes params, just extrap them all
        for ref_flux_ind in range(4):
            ref_freq = component.freqs[0]
            ref_flux = component.fluxes[0][ref_flux_ind]
            
            si_ratio = (200e+6 / ref_freq)**component.si
            exp_bit = np.exp(component.curve_q*(np.log(200e+6 / ref_freq))**2)

            new_ref_flux = ref_flux*si_ratio*exp_bit
            component.fluxes[0][ref_flux_ind] = new_ref_flux
            
            ##TODO for now, only use the Stokes I flux to do the spectral index
            ##extrapolation as everything is Stokes I
            if ref_flux_ind == 0:
                logratio = np.log(ref_freq / 200e+6)
                q = component.curve_q
                new_si = (np.log(ref_flux / new_ref_flux) - q*logratio**2) / logratio
                component.si = new_si
            
        component.norm_comp_cpl = component.fluxes[0][0]
        component.freqs[0] = 200e+6

    return component