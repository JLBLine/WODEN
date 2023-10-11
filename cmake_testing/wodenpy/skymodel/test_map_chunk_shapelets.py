"""

"""

from sys import path
import os
import unittest
import numpy as np

from wodenpy.skymodel.woden_skymodel import Component_Type_Counter, CompTypes
from wodenpy.skymodel.chunk_sky_model import map_chunk_shapelets, Skymodel_Chunk_Map, create_shape_basis_maps

from common_skymodel_test import fill_comp_counter, Expec_Counter, BaseChunkTest

D2R = np.pi/180.0

##for now, WODEN has three flux types: power law, curved power law, and list
NUM_FLUX_TYPES = 3

##Vehicle for running tests
class Test(BaseChunkTest):
    """Test `wodenpy.skymodel.chunk_sky_model.map_chunk_shapelets`
    works correctly"""
    
    def run_test_map_chunk_shapelet(self, max_num_visibilities : int,
                                    num_shapes : int, num_coeff_per_shape : int,
                                    num_time_steps : int, num_list_values : int):
        """Makes a populted Component_Type_Counter based off the given
        inputs, runs it through
        `wodenpy.skymodel.chunk_sky_model.map_chunk_pointgauss` and
        checks that gives the correct answers"""
        
        ## Only testing SHAPELETs, set these to zero
        num_points = 0
        num_gauss = 0
        
        ##set these same for everything
        num_baselines = 5
        num_freqs = 2
                
        full_comp_counter = fill_comp_counter(num_points, num_gauss,
                                            num_shapes, num_coeff_per_shape,
                                            num_list_values, num_time_steps)
        
        full_comp_counter.total_components()
        
        coeffs_per_chunk = int(np.floor(max_num_visibilities / (num_baselines * num_freqs * num_time_steps)))
        num_coeffs_to_chunk = full_comp_counter.total_shape_basis
        num_chunks = int(np.ceil(num_coeffs_to_chunk / coeffs_per_chunk))
        
        # ##Counters to see how far many point/gauss have already been assigned
        # ##to previous chunks
        # expec_counter = Expec_Counter()
        
        shape_basis_to_orig_comp_index_map, shape_basis_to_comp_type_map, shape_basis_param_index = create_shape_basis_maps(full_comp_counter)
        
        for chunk_ind in range(num_chunks):
            chunk_map = map_chunk_shapelets(full_comp_counter,
                                        shape_basis_to_orig_comp_index_map,
                                        shape_basis_to_comp_type_map,
                                        shape_basis_param_index,
                                        chunk_ind, coeffs_per_chunk,
                                        text_file=True)
            
            self.check_shapelet_chunking(chunk_ind, num_coeff_per_shape,
                                         coeffs_per_chunk,
                                         num_list_values,  num_shapes,
                                         full_comp_counter, chunk_map)
            
    ##Ok now run with many many combinations
    
    def test_Shapelet010_Coeff10_Chunk1000_Time004(self):
        max_num_visibilities = 1000
        num_shapes = 10
        num_coeff_per_shape = 10
        num_time_steps = 4

        num_list_values = 3
        
        self.run_test_map_chunk_shapelet(max_num_visibilities,
                                           num_shapes, num_coeff_per_shape,
                                           num_time_steps, num_list_values)
        
        
    def test_Shapelet010_Coeff13_Chunk1000_Time004(self):
        max_num_visibilities = 1000
        num_shapes = 10
        num_coeff_per_shape = 13
        num_time_steps = 4
        num_list_values = 4
        self.run_test_map_chunk_shapelet(max_num_visibilities,
                                           num_shapes, num_coeff_per_shape,
                                           num_time_steps, num_list_values)

    def test_Shapelet193_Coeff1266_Chunk3913_Time007(self):
        max_num_visibilities = 3913
        num_shapes = 193
        num_coeff_per_shape = 1266
        num_time_steps = 7

        num_list_values = 5
        self.run_test_map_chunk_shapelet(max_num_visibilities,
                                           num_shapes, num_coeff_per_shape,
                                           num_time_steps, num_list_values)
        
##Run the test
if __name__ == '__main__':
    unittest.main()
    