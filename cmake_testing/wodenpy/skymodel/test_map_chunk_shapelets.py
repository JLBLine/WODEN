"""

"""

from sys import path
import os
import unittest
import numpy as np

from wodenpy.skymodel.woden_skymodel import Component_Type_Counter, CompTypes
from wodenpy.skymodel.chunk_sky_model import map_chunk_shapelets, Skymodel_Chunk_Map, create_shape_basis_maps, find_num_dirs_per_chunk

from common_skymodel_test import fill_comp_counter_for_chunking, Expec_Counter, BaseChunkTest, Skymodel_Settings

D2R = np.pi/180.0

##for now, WODEN has three flux types: power law, curved power law, and list
NUM_FLUX_TYPES = 3

##Vehicle for running tests
class Test(BaseChunkTest):
    """Test `wodenpy.skymodel.chunk_sky_model.map_chunk_shapelets`
    works correctly"""
    
    def run_test_map_chunk_shapelet(self, settings : Skymodel_Settings):
        """Makes a populted Component_Type_Counter based off the given
        inputs, runs it through
        `wodenpy.skymodel.chunk_sky_model.map_chunk_pointgauss` and
        checks that gives the correct answers"""
        
        max_num_visibilities = settings.max_num_visibilities
        num_shapes = settings.num_shapes
        num_coeff_per_shape = settings.num_coeff_per_shape
        num_time_steps = settings.num_time_steps
        num_list_values = settings.num_list_values
        stokesV_cadence = settings.stokesV_cadence
        linpol_cadence = settings.linpol_cadence
        ## Only testing SHAPELETs, set these to zero
        settings.num_points = 0
        settings.num_gauss = 0
        
        ##set these same for everything
        num_baselines = 5
        num_freqs = 2
        settings.num_baselines = num_baselines
        settings.num_freqs = num_freqs
        
        full_comp_counter = fill_comp_counter_for_chunking(settings)
        
        full_comp_counter.total_components()
        
        coeffs_per_chunk = int(np.floor(max_num_visibilities / (num_baselines * num_freqs * num_time_steps)))
        num_coeffs_to_chunk = full_comp_counter.total_shape_basis
        num_chunks = int(np.ceil(num_coeffs_to_chunk / coeffs_per_chunk))
        
        # ##Counters to see how far many point/gauss have already been assigned
        # ##to previous chunks
        # expec_counter = Expec_Counter()
        
        shape_basis_to_orig_comp_index_map, shape_basis_to_comp_type_map, shape_basis_param_index = create_shape_basis_maps(full_comp_counter)
        
        num_threads=1
        num_shape_dirs = find_num_dirs_per_chunk(full_comp_counter.total_shape_comps, coeffs_per_chunk,
                                                 num_threads)
        
        
        chunk_maps = map_chunk_shapelets(full_comp_counter,
                                        shape_basis_to_orig_comp_index_map,
                                        shape_basis_to_comp_type_map,
                                        shape_basis_param_index,
                                        coeffs_per_chunk, num_shape_dirs)
        
        for chunk_ind, chunk_map in enumerate(chunk_maps):
            
            
            self.check_shapelet_chunking(chunk_ind, num_coeff_per_shape,
                                         coeffs_per_chunk, num_shapes,
                                         settings,
                                         full_comp_counter, chunk_map)
            
    ##Ok now run with many many combinations
    
    def test_Shapelet010_Coeff10_Chunk1000_Time004(self):
        settings = Skymodel_Settings(max_num_visibilities = 1000,
                                     num_shapes = 10,
                                     num_coeff_per_shape = 10,
                                     num_time_steps = 4,
                                     num_list_values = 3)
        
        self.run_test_map_chunk_shapelet(settings)
        
        
    def test_Shapelet010_Coeff13_Chunk1000_Time004(self):
        settings = Skymodel_Settings(max_num_visibilities = 1000,
                                     num_shapes = 10,
                                     num_coeff_per_shape = 13,
                                     num_time_steps = 4,
                                     num_list_values = 4)

        self.run_test_map_chunk_shapelet(settings)

    def test_Shapelet193_Coeff1266_Chunk3913_Time007(self):
        settings = Skymodel_Settings(max_num_visibilities = 3913,
                                     num_shapes = 193,
                                     num_coeff_per_shape = 1266,
                                     num_time_steps = 7,
                                     num_list_values = 5)

        self.run_test_map_chunk_shapelet(settings)
        
    def test_Shapelet193_Coeff1266_Chunk3913_Time007_Lincadence4_Vcadence3(self):
        settings = Skymodel_Settings(max_num_visibilities = 3913,
                                     num_shapes = 193,
                                     num_coeff_per_shape = 1266,
                                     num_time_steps = 7,
                                     num_list_values = 5,
                                     linpol_cadence = 4,
                                     linpol_num_list=3,
                                     linpol_num_p_list=4,
                                     stokesV_cadence = 3,
                                     stokesV_num_list=2)
        
        self.run_test_map_chunk_shapelet(settings)
        
##Run the test
if __name__ == '__main__':
    unittest.main()
    