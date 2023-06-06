import numpy as np
import sys
import os
from typing import Union

##If we are performing a ctest, this check means we use the code we are
##testing and NOT what has been pip or conda installed
try:
    testdir = os.environ['CMAKE_CURRENT_SOURCE_DIR']
    sys.path.append('{:s}/../../../wodenpy'.format(testdir))
    sys.path.append('{:s}/../../../wodenpy/use_libwoden'.format(testdir))
    sys.path.append('{:s}/../../../wodenpy/skymodel'.format(testdir))
    sys.path.append('{:s}/../../../wodenpy/skymodel'.format(testdir))
    
    from woden_skymodel import Component_Type_Counter, Component_Info, CompTypes
    from chunk_sky_model import Skymodel_Chunk_Map
    from skymodel_structs import setup_chunked_source, setup_source_catalogue, Source_Catalogue_Float, Source_Catalogue_Double, add_info_to_source_catalogue, _Ctype_Source_Into_Python, Components_Float, Components_Double, Source_Float, Source_Double
    from beam_settings import BeamTypes
    
except KeyError:
    from wodenpy.skymodel.woden_skymodel import Component_Type_Counter, Component_Info, CompTypes
    from wodenpy.skymodel.chunk_sky_model import Skymodel_Chunk_Map
    from wodenpy.use_libwoden.skymodel_structs import setup_chunked_source, setup_source_catalogue, Source_Catalogue_Float, Source_Catalogue_Double, add_info_to_source_catalogue, _Ctype_Source_Into_Python, Components_Float, Components_Double, Source_Float, Source_Double
    from wodenpy.use_libwoden.beam_settings import BeamTypes


import erfa

D2R = np.pi/180.0

# @profile
def read_yaml_radec_count_components(yaml_path : str):
    """Read just the  ra, dec, and count how many POINT/GAUSS/SHAPE and
    POWER/CURVE/LIST entries there are"""
    
    if not os.path.isfile(yaml_path):
        sys.exit(f"Cannot read sky model from {yaml_path}. Please check your paths, exiting now.")

    comp_counter = Component_Type_Counter()

    with open(yaml_path) as file:
        
        current_source = -1
        first_component = True
        
        ##These are all things we have to count to be able to malloc correct
        ##amount of memory down the line
        for line in file:
            comp_counter.new_file_line()
            if line != '---\n' and '#' not in line and line != ''  and line != ' ' and line != '\n':
                
                if line[0] != ' ':
                    current_source += 1
                    source_name = line[:-1]
                    comp_counter.new_source()
                
                elif 'ra:' in line:
                    ##If this is the first component, no need to reset
                    ##and count things up
                    if first_component:
                        first_component = False
                    else:
                        
                        ##ra should be the first thing in a component, so we need
                        ##to append all the previously found values and reset the
                        ##counters
                        comp_counter.add_component_reset_counters()
                        
                    ##Start a new component
                    comp_counter.initiate_component()
                        
                    ra = float(line.split()[-1])*D2R
                    comp_counter.ra = ra
                    
                elif 'dec:' in line:
                    dec = float(line.split()[-1])*D2R
                    comp_counter.dec = dec

                ##component type related things
                elif 'point' in line:
                    comp_counter.point = 1
                elif 'gaussian:' in line:
                    comp_counter.gaussian = 1
                elif 'shapelet:' in line:
                    comp_counter.shapelet = 1
                elif 'n1:' in line:
                    comp_counter.num_shape_basis += 1
                
                ##flux behaviou related things
                elif 'power_law:' in line and 'curved' not in line:
                    comp_counter.flux_power = 1
                elif 'curved_power_law:' in line:
                    comp_counter.flux_curve = 1
                elif 'list:' in line:
                    comp_counter.flux_list = 1
                elif 'freq:' in line:
                    comp_counter.count_list_flux()

        ##The final component needs to be counted as we used a trigger of the
        ##next to count the pervious components
        comp_counter.add_component_reset_counters()
    
    ##shrink the array length of all info to the amount actually in the model
    comp_counter.remove_array_padding()
    
    ##total up some component counts
    comp_counter.total_components()
        
    return comp_counter


# @profile
def read_yaml_skymodel_chunks(yaml_path : str,
                              chunked_skymodel_maps : list,
                              num_freqs : int, num_time_steps : int,
                              beamtype : int,
                              lsts : np.ndarray, latitude : float,
                              precision = "double") -> Union[Source_Catalogue_Float, Source_Catalogue_Double]:
    
    ##want to know how many shapelets are in all the chunks (used later
    # by "calculate_visiblities.cu")
    num_shapelets = 0
    for chunk_map in chunked_skymodel_maps:
        num_shapelets += chunk_map.n_shapes
        
    ##setup the source catalogue, which is going to store all of the information
    ##of each source and be fed straight into C/CUDA
    source_catalogue = setup_source_catalogue(len(chunked_skymodel_maps), num_shapelets,
                                precision = precision)

    num_comps_all_chunks = 0
    ##for each chunk map, create a Source_Float or Source_Double ctype
    ##struct, and "malloc" the right amount of arrays to store required infor
    for chunk_ind, chunk_map in enumerate(chunked_skymodel_maps):
        chunked_source = setup_chunked_source(chunk_map, num_freqs, num_time_steps,
                                              beamtype, precision=precision)
        
        # chunked_sources.append(chunked_source)
        source_catalogue.sources[chunk_ind] = chunked_source
        
        ##count up the total number of components across all chunks
        ##annoyingly, beacuse Jack sucks, we split shapelet us by basis 
        ##and not component, so this number is actually a combination of
        ##component and basis numbers
        num_comps_all_chunks += chunk_map.n_points + chunk_map.n_gauss + chunk_map.n_shape_coeffs
        
    ##this hold the original index in the sky model file of every
    ##component in this set of chunks
    all_chunk_comp_indexes = np.empty(num_comps_all_chunks)
    
    ##this maps those component indexes to each chunk index
    # map_comp_to_chunk = np.empty(num_comps_all_chunks)
    map_comp_to_chunk = np.full(num_comps_all_chunks, -1)
    lowest_file_numbers = []
    
    low_ind = 0
    for chunk_ind, chunk_map in enumerate(chunked_skymodel_maps):
        
        all_chunk_comp_indexes[low_ind:low_ind+len(chunk_map.all_orig_inds)] = chunk_map.all_orig_inds
        
        map_comp_to_chunk[low_ind:low_ind+len(chunk_map.all_orig_inds)] = chunk_ind
        
        low_ind += chunk_map.n_points + chunk_map.n_gauss + chunk_map.n_shape_coeffs
        
        lowest_file_numbers.append(chunk_map.lowest_file_number)
        
    ##as we iterate through, as soon as we find a component, we add one to
    ##the component index. So start counter one below the
    comp_ind = all_chunk_comp_indexes.min() - 1
    
    min_line_number = min(lowest_file_numbers)
    
    with open(yaml_path) as file:
        
        current_source = -1
        
        use_component = False
        comp_info = False
        
        collected_comps = 0
        
        source_name = "first in read"
        
        line_ind = -1
        
        # for line_ind, line in enumerate(file):
        for line in file:
            line_ind += 1
        
            # comp_counter.new_file_line()
            
            ##Stop iterating if we've collected everything we want
            if collected_comps >= num_comps_all_chunks:
                break
            
            # if line_ind < min(lowest_file_numbers):
            #     pass
            # else:
                
            if line_ind >= min_line_number:
            
                if line != '---\n' and '#' not in line and line != ''  and line != ' ' and line != '\n':
                    
                    if line[0] != ' ':
                        current_source += 1
                        source_name = line[:-1]
                        # comp_counter.new_source()
                    
                    elif 'ra:' in line:
                        
                        ##ra is the first bit of information
                        if comp_info and use_component:
                            
                            comp_info.source_name = source_name
                            collected_comps = add_info_to_source_catalogue(chunked_skymodel_maps,
                                     source_catalogue, comp_ind, comp_info,
                                     map_comp_to_chunk, all_chunk_comp_indexes,
                                     beamtype, lsts, latitude,
                                     collected_comps)
                            
                        ##ra is the first thing in the component, so we know we've
                        ##got a new component
                        comp_ind += 1
                        
                        ##reset all the temporary counters
                        use_component = False
                        
                        ##work out if we want to save information for this
                        ##particular component
                        if comp_ind in all_chunk_comp_indexes:
                            use_component = True
                        
                        comp_info = Component_Info()
                        
                        if use_component:
                            ra = float(line.split()[-1])*D2R
                            comp_info.add_ra(ra)
                            
                    else:
                        if use_component:
                        
                            if 'dec:' in line and use_component:
                                dec = float(line.split()[-1])*D2R
                                comp_info.add_dec(dec)

                            ##component type related things
                            elif 'point' in line and use_component:
                                comp_info.set_point()
                            elif 'gaussian:' in line and use_component:
                                comp_info.set_gaussian()
                            elif 'shapelet:' in line and use_component:
                                comp_info.set_shapelet()
                            
                            ##gaussian/shapelet things
                            elif "maj:" in line and use_component:
                                major = float(line.split()[-1])*(D2R / 3600.0)
                                comp_info.add_major(major)
                            elif "min:" in line and use_component:
                                minor = float(line.split()[-1])*(D2R / 3600.0)
                                comp_info.add_minor(minor)
                            elif "pa:" in line and use_component:
                                pa = float(line.split()[-1])*D2R
                                comp_info.add_pa(pa)

                            ##shapelet things
                            elif 'n1:' in line and use_component:
                                n1 = float(line.split()[-1])
                                comp_info.add_n1(n1)
                            elif 'n2:' in line and use_component:
                                n2 = float(line.split()[-1])
                                comp_info.add_n2(n2)
                            elif 'value:' in line and use_component:
                                coeff = float(line.split()[-1])
                                comp_info.add_coeff(coeff)
                            
                            ##power/curved law related things
                            elif 'si:' in line and use_component:
                                si = float(line.split()[-1])
                                comp_info.add_si(si)
                            
                            ##flux behaviour related things
                            elif 'power_law:' in line and 'curved' not in line and use_component:
                                comp_info.set_flux_power()
                            elif 'curved_power_law:' in line and use_component:
                                comp_info.set_flux_curve()
                            elif 'list:' in line and use_component:
                                comp_info.set_flux_list()
                            
                            elif 'freq:' in line and use_component:
                                freq = float(line.split()[-1])
                                
                                comp_info.add_ref_freq(freq)
                            
                                ##See what indent this freq entry starts at - used to
                                ##line up following freq entries, as `q` can either mean
                                ##stokes Q or q curvature param
                                freq_indent = line.index('f')
                                
                            elif ' i:' in line and 'si' not in line and use_component:
                                stokesI = float(line.split()[-1])
                                comp_info.add_stokesI(stokesI)
                            
                            ##Gotta be fancy here to work out if this is a Stokes Q or a 
                            ##curved power law 'q' param
                            elif ' q:' in line and use_component:
                                q = float(line.split()[-1])
                                if line.index('q') == freq_indent:
                                    comp_info.add_stokesQ(q)
                                else:
                                    if comp_info.flux_curve:
                                        comp_info.add_curve_q(q)
                                        
                            elif ' u:' in line and use_component:
                                stokesU = float(line.split()[-1])
                                comp_info.add_stokesU(stokesU)
                                
                            elif ' v:' in line and use_component:
                                stokesV = float(line.split()[-1])
                                comp_info.add_stokesV(stokesV)
        
        ##We use a new component appearing to trigger collecting information
        ##from the previous component. So for the very last component, need
        ##to check if we need to add the information, and add if needed
        if use_component:
            comp_info.source_name = source_name
            collected_comps = add_info_to_source_catalogue(chunked_skymodel_maps,
                                     source_catalogue, comp_ind, comp_info,
                                     map_comp_to_chunk, all_chunk_comp_indexes,
                                     beamtype, lsts, latitude,
                                     collected_comps)

    ##TODO some kind of consistency check between the chunk_maps and the
    ##sources in the catalogue - make sure we read in the correct information
    
    return source_catalogue

class Components_Noctype(object):
    
    def __init__(self) -> None:
    
        self.ras = None
        self.decs = None
        self.power_ref_freqs = None
        self.power_ref_stokesI = None
        self.power_ref_stokesQ = None
        self.power_ref_stokesU = None
        self.power_ref_stokesV = None
        self.power_SIs = None
        self.curve_ref_freqs = None
        self.curve_ref_stokesI = None
        self.curve_ref_stokesQ = None
        self.curve_ref_stokesU = None
        self.curve_ref_stokesV = None
        self.curve_SIs = None
        self.curve_qs = None
        self.power_comp_inds = None
        self.curve_comp_inds = None
        self.list_comp_inds = None
        self.list_freqs = None
        self.list_stokesI = None
        self.list_stokesQ = None
        self.list_stokesU = None
        self.list_stokesV = None
        self.num_list_values = None
        self.list_start_indexes = None
        self.total_num_flux_entires = None
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
        self.beam_has = None
        self.beam_decs = None
        self.num_primarybeam_values = None

class Source_Noctype(object):
    
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
        self.point_components = Components_Noctype()
        self.gauss_components = Components_Noctype()
        self.shape_components = Components_Noctype()
        
        
class Source_Catalogue_Noctype(object):
    
    def __init__(self) -> None:
    
        self.num_sources = None
        self.num_shapelets = None
        self.sources = []
        
def setup_source_catalogue_noctype(num_sources : int, num_shapelets : int):

    source_catalogue = Source_Catalogue_Noctype()
    
    source_catalogue.sources = [Source_Noctype() for source in range(num_sources)]
    source_catalogue.num_sources = num_sources
    source_catalogue.num_shapelets = num_shapelets

    return source_catalogue
    
    
def setup_components_noctype(chunk_map : Skymodel_Chunk_Map,
                     chunked_source : Source_Noctype,
                     num_freqs : int,
                     num_times : int, comp_type : CompTypes,
                     beamtype : int):
    """choose which components and do the ctypes "malloc" thing
    """
    
    if comp_type == CompTypes.POINT:
        
        n_comps = chunk_map.n_points
        components = chunked_source.point_components
        
        n_comp_powers = chunk_map.n_point_powers
        n_comp_curves = chunk_map.n_point_curves
        n_comp_lists = chunk_map.n_point_lists
        n_lists = chunk_map.point_components.total_num_flux_entires
        
    elif comp_type == CompTypes.GAUSSIAN:
        
        n_comps = chunk_map.n_gauss
        components = chunked_source.gauss_components
        
        n_comp_powers = chunk_map.n_gauss_powers
        n_comp_curves = chunk_map.n_gauss_curves
        n_comp_lists = chunk_map.n_gauss_lists
        n_lists = chunk_map.gauss_components.total_num_flux_entires
        
    elif comp_type == CompTypes.SHAPELET:
        
        n_comps = chunk_map.n_shapes
        components = chunked_source.shape_components
        
        n_comp_powers = chunk_map.n_shape_powers
        n_comp_curves = chunk_map.n_shape_curves
        n_comp_lists = chunk_map.n_shape_lists
        n_lists = chunk_map.shape_components.total_num_flux_entires
        
    ##component type specific things
    num_primarybeam_values = n_comps*num_freqs*num_times
        
    components.ras = np.empty(n_comps, dtype=float)
    components.decs = np.empty(n_comps, dtype=float)
    components.num_primarybeam_values = num_primarybeam_values
    
    ##power-law flux things
    components.power_ref_freqs = np.empty(n_comp_powers)
    components.power_ref_stokesI = np.empty(n_comp_powers)
    components.power_ref_stokesQ = np.empty(n_comp_powers)
    components.power_ref_stokesU = np.empty(n_comp_powers)
    components.power_ref_stokesV = np.empty(n_comp_powers)
    components.power_SIs = np.empty(n_comp_powers)
    
    ##curved power-law flux things
    components.curve_ref_freqs = np.empty(n_comp_curves)
    components.curve_ref_stokesI = np.empty(n_comp_curves)
    components.curve_ref_stokesQ = np.empty(n_comp_curves)
    components.curve_ref_stokesU = np.empty(n_comp_curves)
    components.curve_ref_stokesV = np.empty(n_comp_curves)
    components.curve_SIs = np.empty(n_comp_curves)
    components.curve_qs = np.empty(n_comp_curves)
    
    ##list flux things
    components.list_freqs = np.empty(n_lists) ##this is always double
    components.list_stokesI = np.empty(n_lists)
    components.list_stokesQ = np.empty(n_lists)
    components.list_stokesU = np.empty(n_lists)
    components.list_stokesV = np.empty(n_lists)
    
    components.num_list_values = np.empty(n_comp_lists, dtype=int)
    components.list_start_indexes = np.empty(n_comp_lists, dtype=int)
    components.total_num_flux_entires = n_lists
    
    # print("WE BE ASSIGNED THIS NUM", components.total_num_flux_entires)
    
    ##indexes of types
    components.power_comp_inds = np.empty(n_comp_powers, dtype=int)
    components.curve_comp_inds = np.empty(n_comp_curves, dtype=int)
    components.list_comp_inds = np.empty(n_comp_lists, dtype=int)
    
    if comp_type == CompTypes.GAUSSIAN or comp_type == CompTypes.SHAPELET:
        components.majors = np.empty(n_comps)
        components.minors = np.empty(n_comps)
        components.pas = np.empty(n_comps)
        
        
    if comp_type == CompTypes.SHAPELET:
        ncoeffs = chunk_map.shape_components.total_shape_coeffs
        components.shape_coeffs = np.empty(ncoeffs)
        components.n1s = np.empty(ncoeffs)
        components.n2s = np.empty(ncoeffs)
        components.param_indexes = np.empty(ncoeffs)
        
    ##----------------------------------------------------------------
    ##now we make space for coordinates that are need for the primary beam
    
    if beamtype == BeamTypes.GAUSS_BEAM.value or beamtype == BeamTypes.MWA_ANALY.value:
        components.beam_has = np.empty(n_comps*num_times)
        components.beam_decs = np.empty(n_comps*num_times)
        
    ##only the NO_BEAM and GAUSS_BEAM options don't need az,za coords
    if beamtype == BeamTypes.GAUSS_BEAM.value or beamtype == BeamTypes.NO_BEAM.value:
        pass
    else:
        components.azs = np.empty(n_comps*num_times)
        components.zas = np.empty(n_comps*num_times)
    
    
def setup_chunked_source_noctype(chunk_map : Skymodel_Chunk_Map, num_freqs : int,
                         num_times : int, beamtype : int):
    
    chunked_source = Source_Noctype()
    
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
        setup_components_noctype(chunk_map, chunked_source, num_freqs, num_times,
                         CompTypes.POINT, beamtype)
        
    if chunk_map.n_gauss > 0:
        setup_components_noctype(chunk_map, chunked_source, num_freqs, num_times,
                         CompTypes.GAUSSIAN, beamtype)
        
    if chunk_map.n_shapes > 0:
        setup_components_noctype(chunk_map, chunked_source, num_freqs, num_times,
                         CompTypes.SHAPELET, beamtype)
    
    
    return chunked_source
    

def read_yaml_skymodel_chunks_noctype(yaml_path : str,
                              chunked_skymodel_maps : list,
                              num_freqs : int, num_time_steps : int,
                              beamtype : int, lsts : np.ndarray,
                              latitude : float) -> Source_Catalogue_Noctype:
    
    ##want to know how many shapelets are in all the chunks (used later
    # by "calculate_visiblities.cu")
    num_shapelets = 0
    for chunk_map in chunked_skymodel_maps:
        num_shapelets += chunk_map.n_shapes
        
    ##setup the source catalogue, which is going to store all of the information
    ##of each source and be fed straight into C/CUDA
    source_catalogue = setup_source_catalogue_noctype(len(chunked_skymodel_maps),
                                                      num_shapelets)

    num_comps_all_chunks = 0
    ##for each chunk map, create a Source_Float or Source_Double ctype
    ##struct, and "malloc" the right amount of arrays to store required infor
    for chunk_ind, chunk_map in enumerate(chunked_skymodel_maps):
        chunked_source = setup_chunked_source_noctype(chunk_map, num_freqs, num_time_steps,
                                              beamtype)
        
        # chunked_sources.append(chunked_source)
        source_catalogue.sources[chunk_ind] = chunked_source
        
        ##count up the total number of components across all chunks
        ##annoyingly, beacuse Jack sucks, we split shapelet us by basis 
        ##and not component, so this number is actually a combination of
        ##component and basis numbers
        num_comps_all_chunks += chunk_map.n_points + chunk_map.n_gauss + chunk_map.n_shape_coeffs
        
    ##this hold the original index in the sky model file of every
    ##component in this set of chunks
    all_chunk_comp_indexes = np.empty(num_comps_all_chunks)
    
    ##this maps those component indexes to each chunk index
    # map_comp_to_chunk = np.empty(num_comps_all_chunks)
    map_comp_to_chunk = np.full(num_comps_all_chunks, -1)
    lowest_file_numbers = []
    
    low_ind = 0
    for chunk_ind, chunk_map in enumerate(chunked_skymodel_maps):
        
        all_chunk_comp_indexes[low_ind:low_ind+len(chunk_map.all_orig_inds)] = chunk_map.all_orig_inds
        
        map_comp_to_chunk[low_ind:low_ind+len(chunk_map.all_orig_inds)] = chunk_ind
        
        low_ind += chunk_map.n_points + chunk_map.n_gauss + chunk_map.n_shape_coeffs
        
        lowest_file_numbers.append(chunk_map.lowest_file_number)
        
    ##as we iterate through, as soon as we find a component, we add one to
    ##the component index. So start counter one below the
    comp_ind = all_chunk_comp_indexes.min() - 1
    
    min_line_number = min(lowest_file_numbers)
    
    with open(yaml_path) as file:
        
        current_source = -1
        
        use_component = False
        comp_info = False
        
        collected_comps = 0
        
        source_name = "first in read"
        
        line_ind = -1
        
        # for line_ind, line in enumerate(file):
        for line in file:
            line_ind += 1
        
            # comp_counter.new_file_line()
            
            ##Stop iterating if we've collected everything we want
            if collected_comps >= num_comps_all_chunks:
                break
            
            # if line_ind < min(lowest_file_numbers):
            #     pass
            # else:
                
            if line_ind >= min_line_number:
            
                if line != '---\n' and '#' not in line and line != ''  and line != ' ' and line != '\n':
                    
                    if line[0] != ' ':
                        current_source += 1
                        source_name = line[:-1]
                        # comp_counter.new_source()
                    
                    elif 'ra:' in line:
                        
                        ##ra is the first bit of information
                        if comp_info and use_component:
                            
                            comp_info.source_name = source_name
                            collected_comps = add_info_to_source_catalogue(chunked_skymodel_maps,
                                     source_catalogue, comp_ind, comp_info,
                                     map_comp_to_chunk, all_chunk_comp_indexes,
                                     beamtype, lsts, latitude,
                                     collected_comps)
                            
                        ##ra is the first thing in the component, so we know we've
                        ##got a new component
                        comp_ind += 1
                        
                        ##reset all the temporary counters
                        use_component = False
                        
                        ##work out if we want to save information for this
                        ##particular component
                        if comp_ind in all_chunk_comp_indexes:
                            use_component = True
                        
                        comp_info = Component_Info()
                        
                        if use_component:
                            ra = float(line.split()[-1])*D2R
                            comp_info.add_ra(ra)
                            
                    else:
                        if use_component:
                        
                            if 'dec:' in line and use_component:
                                dec = float(line.split()[-1])*D2R
                                comp_info.add_dec(dec)

                            ##component type related things
                            elif 'point' in line and use_component:
                                comp_info.set_point()
                            elif 'gaussian:' in line and use_component:
                                comp_info.set_gaussian()
                            elif 'shapelet:' in line and use_component:
                                comp_info.set_shapelet()
                            
                            ##gaussian/shapelet things
                            elif "maj:" in line and use_component:
                                major = float(line.split()[-1])*(D2R / 3600.0)
                                comp_info.add_major(major)
                            elif "min:" in line and use_component:
                                minor = float(line.split()[-1])*(D2R / 3600.0)
                                comp_info.add_minor(minor)
                            elif "pa:" in line and use_component:
                                pa = float(line.split()[-1])*D2R
                                comp_info.add_pa(pa)

                            ##shapelet things
                            elif 'n1:' in line and use_component:
                                n1 = float(line.split()[-1])
                                comp_info.add_n1(n1)
                            elif 'n2:' in line and use_component:
                                n2 = float(line.split()[-1])
                                comp_info.add_n2(n2)
                            elif 'value:' in line and use_component:
                                coeff = float(line.split()[-1])
                                comp_info.add_coeff(coeff)
                            
                            ##power/curved law related things
                            elif 'si:' in line and use_component:
                                si = float(line.split()[-1])
                                comp_info.add_si(si)
                            
                            ##flux behaviour related things
                            elif 'power_law:' in line and 'curved' not in line and use_component:
                                comp_info.set_flux_power()
                            elif 'curved_power_law:' in line and use_component:
                                comp_info.set_flux_curve()
                            elif 'list:' in line and use_component:
                                comp_info.set_flux_list()
                            
                            elif 'freq:' in line and use_component:
                                freq = float(line.split()[-1])
                                
                                comp_info.add_ref_freq(freq)
                            
                                ##See what indent this freq entry starts at - used to
                                ##line up following freq entries, as `q` can either mean
                                ##stokes Q or q curvature param
                                freq_indent = line.index('f')
                                
                            elif ' i:' in line and 'si' not in line and use_component:
                                stokesI = float(line.split()[-1])
                                comp_info.add_stokesI(stokesI)
                            
                            ##Gotta be fancy here to work out if this is a Stokes Q or a 
                            ##curved power law 'q' param
                            elif ' q:' in line and use_component:
                                q = float(line.split()[-1])
                                if line.index('q') == freq_indent:
                                    comp_info.add_stokesQ(q)
                                else:
                                    if comp_info.flux_curve:
                                        comp_info.add_curve_q(q)
                                        
                            elif ' u:' in line and use_component:
                                stokesU = float(line.split()[-1])
                                comp_info.add_stokesU(stokesU)
                                
                            elif ' v:' in line and use_component:
                                stokesV = float(line.split()[-1])
                                comp_info.add_stokesV(stokesV)
        
        ##We use a new component appearing to trigger collecting information
        ##from the previous component. So for the very last component, need
        ##to check if we need to add the information, and add if needed
        if use_component:
            comp_info.source_name = source_name
            collected_comps = add_info_to_source_catalogue(chunked_skymodel_maps,
                                     source_catalogue, comp_ind, comp_info,
                                     map_comp_to_chunk, all_chunk_comp_indexes,
                                     beamtype, lsts, latitude,
                                     collected_comps)

    ##TODO some kind of consistency check between the chunk_maps and the
    ##sources in the catalogue - make sure we read in the correct information

    return source_catalogue


def _convert_components_to_ctype(source_py : Source_Noctype,
                source : Union[Source_Float, Source_Double],
                comp_type : CompTypes,
                num_freqs : int, num_times : int,
                beamtype : int, precision = 'double'):
    """All ctype arrays in `components` should have been allocated memory,
    so we just have to copy across the values from `components_py`.
    
    You can't just assign an array"""
    
    if comp_type == CompTypes.POINT:
        
        n_comps = source_py.n_points
        components_py = source_py.point_components
        components = source.point_components
        
        n_comp_powers = source_py.n_point_powers
        n_comp_curves = source_py.n_point_curves
        n_comp_lists = source_py.n_point_lists
        n_lists = source_py.point_components.total_num_flux_entires
        
    elif comp_type == CompTypes.GAUSSIAN:
        
        n_comps = source_py.n_gauss
        components_py = source_py.gauss_components
        components = source.gauss_components
        
        n_comp_powers = source_py.n_gauss_powers
        n_comp_curves = source_py.n_gauss_curves
        n_comp_lists = source_py.n_gauss_lists
        n_lists = source_py.gauss_components.total_num_flux_entires
        
    elif comp_type == CompTypes.SHAPELET:
        
        n_comps = source_py.n_shapes
        components_py = source_py.shape_components
        components = source.shape_components
        
        n_comp_powers = source_py.n_shape_powers
        n_comp_curves = source_py.n_shape_curves
        n_comp_lists = source_py.n_shape_lists
        n_lists = source_py.shape_components.total_num_flux_entires
        
    ##Numbers we need down the line
    components.num_primarybeam_values = components_py.num_primarybeam_values
    components.total_num_flux_entires = n_lists
    
    ##is this gonna be super slow? can do something like
    # carray = narray.data_as(ctypes.POINTER(ctypes.c_double*num_comps)).contents
    ##to convert a numpy to a carray
    
    for comp_ind in range(n_comps):
        components.ras[comp_ind] = components_py.ras[comp_ind]
        components.decs[comp_ind] = components_py.decs[comp_ind]
        
    ##We have power-law type components, so fill them in
    for pow_in in range(n_comp_powers):
        components.power_ref_freqs[pow_in] = components_py.power_ref_freqs[pow_in]
        components.power_ref_stokesI[pow_in] = components_py.power_ref_stokesI[pow_in]
        components.power_ref_stokesQ[pow_in] = components_py.power_ref_stokesQ[pow_in]
        components.power_ref_stokesU[pow_in] = components_py.power_ref_stokesU[pow_in]
        components.power_ref_stokesV[pow_in] = components_py.power_ref_stokesV[pow_in]
        components.power_SIs[pow_in] = components_py.power_SIs[pow_in]
        
        components.power_comp_inds[pow_in] = components_py.power_comp_inds[pow_in]
        
    ##We have curved power-law type components, so fill them in
    for cur_in in range(n_comp_curves):
        components.curve_ref_freqs[cur_in] = components_py.curve_ref_freqs[cur_in]
        components.curve_ref_stokesI[cur_in] = components_py.curve_ref_stokesI[cur_in]
        components.curve_ref_stokesQ[cur_in] = components_py.curve_ref_stokesQ[cur_in]
        components.curve_ref_stokesU[cur_in] = components_py.curve_ref_stokesU[cur_in]
        components.curve_ref_stokesV[cur_in] = components_py.curve_ref_stokesV[cur_in]
        components.curve_SIs[cur_in] = components_py.curve_SIs[cur_in]
        components.curve_qs[cur_in] = components_py.curve_qs[cur_in]
        
        components.curve_comp_inds[cur_in] = components_py.curve_comp_inds[cur_in]
        
    for lis_in in range(n_lists):
        components.list_freqs[lis_in] = components_py.list_freqs[lis_in]
        components.list_stokesI[lis_in] = components_py.list_stokesI[lis_in]
        components.list_stokesQ[lis_in] = components_py.list_stokesQ[lis_in]
        components.list_stokesU[lis_in] = components_py.list_stokesU[lis_in]
        components.list_stokesV[lis_in] = components_py.list_stokesV[lis_in]
        
    for lis_in in range(n_comp_lists):
        components.num_list_values[lis_in] = components_py.num_list_values[lis_in]
        components.list_start_indexes[lis_in] = components_py.list_start_indexes[lis_in]
        
        components.list_comp_inds[lis_in] = components_py.list_comp_inds[lis_in]
        
    if comp_type == CompTypes.GAUSSIAN or comp_type == CompTypes.SHAPELET:
        for comp_ind in range(n_comps):
            components.majors[comp_ind] = components_py.majors[comp_ind]
            components.minors[comp_ind] = components_py.minors[comp_ind]
            components.pas[comp_ind] = components_py.pas[comp_ind]
    
    if comp_type == CompTypes.SHAPELET:
        for coeff in range(source_py.n_shape_coeffs):
            components.shape_coeffs[coeff] = components_py.shape_coeffs[coeff]
            components.n1s[coeff] = components_py.n1s[coeff]
            components.n2s[coeff] = components_py.n2s[coeff]
            components.param_indexes[coeff] = components_py.param_indexes[coeff]
            
            
    if beamtype == BeamTypes.GAUSS_BEAM.value or beamtype == BeamTypes.MWA_ANALY.value:
        
        for beam_in in range(n_comps*num_times):
        
            components.beam_has[beam_in] = components_py.beam_has[beam_in]
            components.beam_decs[beam_in] = components_py.beam_decs[beam_in]
        
    ##only the NO_BEAM and GAUSS_BEAM options don't need az,za coords
    if beamtype == BeamTypes.GAUSS_BEAM.value or beamtype == BeamTypes.NO_BEAM.value:
        pass
    else:
        for beam_in in range(n_comps*num_times):
        
            components.azs[beam_in] = components_py.azs[beam_in]
            components.zas[beam_in] = components_py.zas[beam_in]
    
    return

def _convert_source_to_ctype(source_py : Source_Noctype,
                source : Union[Source_Float, Source_Double],
                num_freqs : int, num_time_steps : int,
                beamtype : int):
    
    source.n_points = source_py.n_points
    source.n_point_lists = source_py.n_point_lists
    source.n_point_powers = source_py.n_point_powers
    source.n_point_curves = source_py.n_point_curves
    source.n_gauss = source_py.n_gauss
    source.n_gauss_lists = source_py.n_gauss_lists
    source.n_gauss_powers = source_py.n_gauss_powers
    source.n_gauss_curves = source_py.n_gauss_curves
    source.n_shapes = source_py.n_shapes
    source.n_shape_lists = source_py.n_shape_lists
    source.n_shape_powers = source_py.n_shape_powers
    source.n_shape_curves = source_py.n_shape_curves
    source.n_shape_coeffs = source_py.n_shape_coeffs
    source.n_comps = source_py.n_comps
    
    if source_py.n_points:
        _convert_components_to_ctype(source_py, source,
                                     CompTypes.POINT,
                                     num_freqs, num_time_steps,
                                     beamtype)
        
    if source_py.n_gauss:
        _convert_components_to_ctype(source_py, source,
                                     CompTypes.GAUSSIAN,
                                     num_freqs, num_time_steps,
                                     beamtype)
        
    if source_py.n_shapes:
        _convert_components_to_ctype(source_py, source,
                                     CompTypes.SHAPELET,
                                     num_freqs, num_time_steps,
                                     beamtype)
    
    return 

        
    

def convert_source_catalogue_to_ctypes(
                                source_catalogue_py : Source_Catalogue_Noctype,
                                chunked_skymodel_maps : list,
                                num_freqs : int, num_time_steps : int,
                                beamtype : int, precision='double'):
    
    # print(source_catalogue_py.num_sources)
    # print(source_catalogue_py.num_shapelets)

    source_catalogue = setup_source_catalogue(source_catalogue_py.num_sources,
                                        source_catalogue_py.num_shapelets,
                                        precision)
    
    for chunk_ind, chunk_map in enumerate(chunked_skymodel_maps):
        chunked_source = setup_chunked_source(chunk_map, num_freqs, num_time_steps,
                                              beamtype, precision)
        
        # chunked_sources.append(chunked_source)
        source_catalogue.sources[chunk_ind] = chunked_source
        
        _convert_source_to_ctype(source_catalogue_py.sources[chunk_ind],
                                 source_catalogue.sources[chunk_ind],
                                 num_freqs, num_time_steps, beamtype)

    return source_catalogue