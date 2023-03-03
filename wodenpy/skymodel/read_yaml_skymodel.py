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
    
    from woden_skymodel import Component_Type_Counter, Component_Info, CompTypes
    from chunk_sky_model import Skymodel_Chunk_Map
    from skymodel_structs import setup_chunked_source, setup_source_catalogue, Source_Catalogue_Float, Source_Catalogue_Double, add_info_to_source_catalogue
    
except KeyError:
    from wodenpy.skymodel.woden_skymodel import Component_Type_Counter, Component_Info, CompTypes
    from wodenpy.skymodel.chunk_sky_model import Skymodel_Chunk_Map
    from wodenpy.use_libwoden.skymodel_structs import setup_chunked_source, setup_source_catalogue, Source_Catalogue_Float, Source_Catalogue_Double, add_info_to_source_catalogue


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
                    # print("HAH", current_source)
                    source_name = line[:-1]
                    comp_counter.new_source()
                    # print("WHY", comp_counter.comp_index,
                    #              comp_counter.source_index)
                
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



def read_yaml_skymodel_chunks(yaml_path : str,
                              chunked_skymodel_maps : list,
                              num_freqs, num_time_steps,
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
                                                precision=precision)
        
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
        # print("orig_index_map", chunk_ind, chunk_map.all_orig_inds )
        
        all_chunk_comp_indexes[low_ind:low_ind+len(chunk_map.all_orig_inds)] = chunk_map.all_orig_inds
        
        map_comp_to_chunk[low_ind:low_ind+len(chunk_map.all_orig_inds)] = chunk_ind
        
        low_ind += chunk_map.n_points + chunk_map.n_gauss + chunk_map.n_shape_coeffs
        
        lowest_file_numbers.append(chunk_map.lowest_file_number)
        
    # print("all_chunk_comp_indexes", all_chunk_comp_indexes)
    # print("map_comp_to_chunk", map_comp_to_chunk)
        
    ##as we iterate through, as soon as we find a component, we add one to
    ##the component index. So start counter one below the
    comp_ind = all_chunk_comp_indexes.min() - 1
    
    with open(yaml_path) as file:
        
        current_source = -1
        
        use_component = False
        comp_info = False
        
        collected_comps = 0
        
        source_name = "first in read"
        
        for line_ind, line in enumerate(file):
            # comp_counter.new_file_line()
            
            ##Stop iterating if we've collected everything we want
            if collected_comps == num_comps_all_chunks:
                break
            
            if line_ind < min(lowest_file_numbers):
                pass
            else:
            
                if line != '---\n' and '#' not in line and line != ''  and line != ' ' and line != '\n':
                    
                    if line[0] != ' ':
                        current_source += 1
                        source_name = line[:-1]
                        # comp_counter.new_source()
                        # print("WHY", comp_counter.comp_index,
                        #              comp_counter.source_index)
                    
                    elif 'ra:' in line:
                        
                        ##ra is the first bit of information
                        if comp_info and use_component:
                            
                            empty_fluxes = comp_info.finalise_comp()
                            
                            if len(empty_fluxes) > 0:
                                print(f"WARNING: {len(empty_fluxes)} components in source '{source_name}' have no flux values. Setting to zero")
                            
                            collected_comps = add_info_to_source_catalogue(chunked_skymodel_maps,
                                     source_catalogue, comp_ind, comp_info,
                                     map_comp_to_chunk, all_chunk_comp_indexes,
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
                        
                    elif 'dec:' in line and use_component:
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
            collected_comps = add_info_to_source_catalogue(chunked_skymodel_maps,
                                     source_catalogue, comp_ind, comp_info,
                                     map_comp_to_chunk, all_chunk_comp_indexes,
                                     collected_comps)

    # print(collected_comps, num_comps_all_chunks)
    
    
    ##TODO some kind of consistency check between the chunk_maps and the
    ##sources in the catalogue - make sure we read in the correct information
    
    return source_catalogue


# def read_yaml_all_info(yaml_path : str,
#                        comp_counter : Component_Type_Counter,
#                        include_flags : np.array, 
#                        full_catalogue = False,
#                        chunk_size = 1e+9):
#     """Read all the information and dump into chunk ctype classes """

#     ##Read the whole thing including sources below the horizon
#     if full_catalogue:
#         include_flags = np.ones(comp_counter.total_comps)
        
#     with open(yaml_path) as file, open('outlines.txt', 'w') as outfile:
#     # lines = file.read().split('\n')
    
#         current_source = 0
        
#         ##These are all things we have to count to be able to malloc correct
#         ##amount of memory down the line

#         for line in file:
#             if line != '---\n' and '#' not in line and line != ''  and line != ' ' and line != '\n':
                
#                 if line[0] != ' ':
#                     # print(current_source)
#                     outfile.write(line[:-2] + '\n')
#                     comp_counter.source_names.append(line[:-2])
#                     current_source += 1
                
#                 elif 'ra:' in line:

#                     ##ra should be the first thing in a component, so we need
#                     ##to append all the previously found values and reset the
#                     ##counters
#                     comp_counter.add_component_reset_counters()

#                     comp_counter.comp_ras.append(float(line.split()[-1]))
#                     comp_counter.source_indexes.append(current_source)

#                 elif 'dec:' in line:
#                     comp_counter.comp_decs.append(float(line.split()[-1]))

#                 elif 'comp_type: point' in line:
#                     comp_counter.point = 1
#                 elif 'gaussian:' in line:
#                     comp_counter.gaussian = 1
#                 elif 'shapelet:' in line:
#                     comp_counter.shapelet = 1

#                 elif 'n1:' in line:

#                     comp_counter.num_shape_basis += 1

#                 elif 'freq:' in line:
#                     comp_counter.count_list_flux()
        

