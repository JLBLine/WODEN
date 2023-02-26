import numpy as np
import sys
import os

##If we are performing a ctest, this check means we use the code we are
##testing and NOT what has been pip or conda installed
try:
    testdir = os.environ['CMAKE_CURRENT_SOURCE_DIR']
    sys.path.append('{:s}/../../../wodenpy'.format(testdir))
    
    from skymodel.woden_skymodel import Component_Type_Counter
    from skymodel.chunk_sky_model import Skymodel_Chunk_Map
    
except KeyError:
    from wodenpy.skymodel.woden_skymodel import Component_Type_Counter
    from wodenpy.skymodel.chunk_sky_model import Skymodel_Chunk_Map


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
                    source_name = line[:-2]
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
    
    comp_counter.total_components()
        
    return comp_counter


def read_yaml_skymodel_chunks(yaml_path : str,
                              chunked_skymodel_maps : list,
                              chunk_indexes = False,
                              precision = "double"):

    if chunk_indexes:
        pass
    else:
        chunk_indexes = range(len(cchunked_skymodel_maps))

    populated_skymodel_chunks = []

    for chunk_ind in chunk_indexes:
        skymodel_chunk = setup_source(chunked_skymodel_maps[chunk_ind], precision)

        populated_skymodel_chunks.append(skymodel_chunk)

    ##TODO
    ## - find some kind of map (fucking again??) of all the components to
    ##   read in, and which chunk they belong to
    ## - find lowest filenumber from all the selected sky chunks
    ## - find what component number that corresponds to, so we can count
    ##   from that number onward
    ## - then we iterate over all lines in the file, grabbing components
    ##   that we want, and shoving them into the correct skymodel_chunk


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
        

