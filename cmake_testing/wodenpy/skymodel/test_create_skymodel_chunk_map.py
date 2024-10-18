"""

"""

from sys import path
import os
import unittest
import numpy as np

##This is where our code lives
# code_dir = os.environ['CMAKE_CURRENT_SOURCE_DIR']

# ##Code we are testing
from wodenpy.skymodel import read_yaml_skymodel
# import wodenpy
from wodenpy.skymodel.woden_skymodel import Component_Type_Counter, CompTypes
from wodenpy.skymodel.chunk_sky_model import create_skymodel_chunk_map, Skymodel_Chunk_Map

from common_skymodel_test import fill_comp_counter_for_chunking, Expec_Counter, BaseChunkTest, Skymodel_Settings

import binpacking

D2R = np.pi/180.0

##for now, WODEN has three flux types: power law, curved power law, and list
NUM_FLUX_TYPES = 3

##based class on an existing class test that includes methods for
##for 
class Test(BaseChunkTest):
    
    def run_test_create_skymodel_chunk_map(self, settings : Skymodel_Settings):
        """Makes a populted Component_Type_Counter based off the given
        inputs, runs it through
        `wodenpy.skymodel.chunk_sky_model.map_chunk_pointgauss` and
        checks that gives the correct answers"""
        
        
        max_num_visibilities = settings.max_num_visibilities
        num_time_steps = settings.num_time_steps
        num_baselines = settings.num_baselines
        num_freqs = settings.num_freqs
        
        num_points = settings.num_points
        num_gauss = settings.num_gauss
        num_shapes = settings.num_shapes
        num_coeff_per_shape = settings.num_coeff_per_shape
        
        ##make the fake sky model        
        comp_counter = fill_comp_counter_for_chunking(settings)
        
        comp_counter.total_components()
        
        comp_counter.print_info()
        
        
        ##Run the code we are testing!
        chunked_skymodel_counters = create_skymodel_chunk_map(comp_counter,
                                        max_num_visibilities, num_baselines,
                                        num_freqs, num_time_steps,
                                        max_chunks_per_set=1e5)
        
        # print("DIS", chunked_skymodel_counters.shape)
        # print("DIS", len(chunked_skymodel_counters[0,0]), chunked_skymodel_counters[0,0])
        
        comps_per_chunk = int(np.floor(max_num_visibilities / (num_baselines * num_freqs * num_time_steps)))
        
        ##pray this never happens, probably means we're going to run out of
        ##GPU memory TODO don't pray, submit a warning?
        if comps_per_chunk < 1: comps_per_chunk = 1
        
        num_point_chunks = int(np.ceil(NUM_FLUX_TYPES*num_points / comps_per_chunk))
        num_gauss_chunks = int(np.ceil(NUM_FLUX_TYPES*num_gauss / comps_per_chunk))
        
        ##we split SHAPELET by the basis components (number of coeffs)
        num_coeff_chunks = int(np.ceil(comp_counter.total_shape_basis / comps_per_chunk))
        
        ##Counters to see how far many point/gauss have already been assigned
        ##to previous chunks
        expec_counter = Expec_Counter()
        
        for point_ind in range(num_point_chunks):
            
            # chunk_map = chunked_skymodel_counters[point_ind][0][0]
            chunk_map = chunked_skymodel_counters[0,0][point_ind]
            
            expec_counter = self.check_pointgauss_chunking(point_ind,
                                            comps_per_chunk,
                                            num_point_chunks, num_points,
                                            settings,
                                            CompTypes.POINT,
                                            comp_counter, chunk_map,
                                            expec_counter)
            
        for gauss_ind in range(num_gauss_chunks):
            
            # chunk_map = chunked_skymodel_counters[num_point_chunks + gauss_ind][0][0]
            chunk_map = chunked_skymodel_counters[0,0][num_point_chunks + gauss_ind]
            
            expec_counter = self.check_pointgauss_chunking(gauss_ind,
                                            comps_per_chunk,
                                            num_gauss_chunks, num_gauss,
                                            settings,
                                            CompTypes.GAUSSIAN,
                                            comp_counter, chunk_map,
                                            expec_counter)
        
        if num_coeff_chunks > 0:
            num_shapes_per_comp = []
                
            for coeff_ind in range(num_coeff_chunks):
                chunk_map = 'meh'
                
                num_shapes_comp = self.check_shapelet_chunking(coeff_ind, num_coeff_per_shape,
                                            comps_per_chunk, num_shapes,
                                            settings,
                                            comp_counter, chunk_map,
                                            total_point_comps=NUM_FLUX_TYPES*num_points,
                                            total_gauss_comps=NUM_FLUX_TYPES*num_gauss,
                                            do_check=False)
                num_shapes_per_comp.append(num_shapes_comp)
                
            ##We will have some unedfined number of chunks, so we want to split
            ##things as evenly as possible in the available number of threads
            indexed_shape_chunk_sizes = [(i, n_shape) for i,n_shape in enumerate(num_shapes_per_comp)]  # List of (index, value) tuples
            target_volume = comps_per_chunk  # Set the target volume for each bin

            # Step 2: Partition the numbers while keeping track of indices using the `to_constant_volume` function
            binned_shape_chunk_sizes = binpacking.to_constant_volume(indexed_shape_chunk_sizes, target_volume, weight_pos=1)
            
            # print(len(binned_shape_chunk_sizes), binned_shape_chunk_sizes)
            num_threads = 1
            if len(binned_shape_chunk_sizes) > num_threads:
                while len(binned_shape_chunk_sizes) > num_threads:
                    # Find the two smallest binned_shape_chunk_sizes and merge them
                    binned_shape_chunk_sizes = sorted(binned_shape_chunk_sizes, key=lambda bin: sum(item[1] for item in bin))  # Sort binned_shape_chunk_sizes by their total sum
                    binned_shape_chunk_sizes[0].extend(binned_shape_chunk_sizes[1])  # Merge the two smallest binned_shape_chunk_sizes
                    binned_shape_chunk_sizes.pop(1)  # Remove the now-empty bin

            shape_comp_chunk_order = []

            for bin_index_size in binned_shape_chunk_sizes:
                for index, value in bin_index_size:
                    shape_comp_chunk_order.append(index)
                    
            for new_ind, coeff_ind in enumerate(shape_comp_chunk_order):
                # print(num_coeff_chunks, num_point_chunks, num_gauss_chunks)
                # print(len(chunked_skymodel_counters[0,0])-num_point_chunks-num_gauss_chunks)
                
                # chunk_map = chunked_skymodel_counters[-1][0][coeff_ind]
                chunk_map = chunked_skymodel_counters[0,0][num_point_chunks + num_gauss_chunks + new_ind]
                
                self.check_shapelet_chunking(coeff_ind, num_coeff_per_shape,
                                            comps_per_chunk, num_shapes,
                                            settings,
                                            comp_counter, chunk_map,
                                            total_point_comps=NUM_FLUX_TYPES*num_points,
                                            total_gauss_comps=NUM_FLUX_TYPES*num_gauss)
            
    ##Ok now run with many many combinations
    
    def test_P3_G0000_S00_00_C1_Time001(self):
        max_num_visibilities = 2
        
        settings = Skymodel_Settings(num_points = 3,
                                     num_gauss = 0,
                                     num_shapes = 0,
                                     num_coeff_per_shape = 0,
                                     num_time_steps = 1,
                                     num_freqs = 1,
                                     num_baselines = 8128,
                                     num_list_values = 3)
        self.run_test_create_skymodel_chunk_map(settings),
        


    def test_P4544_G1736_S20_14_C1e8_Time014(self):

        settings = Skymodel_Settings(num_points = 4544,
                                    num_gauss = 1736,
                                    num_shapes = 20,
                                    num_coeff_per_shape = 14,
                                    max_num_visibilities = 1e8,
                                    num_time_steps = 14,
                                    num_baselines = 8128,
                                    num_freqs = 16,
                                    num_list_values = 4)
        self.run_test_create_skymodel_chunk_map(settings)


    def test_P1000_G0000_S00_00_C1e8_Time014(self):

        settings = Skymodel_Settings(num_points = 1000,
                                    num_gauss = 0,
                                    num_shapes = 0,
                                    num_coeff_per_shape = 0,
                                    max_num_visibilities = 1e9,
                                    num_time_steps = 14,
                                    num_baselines = 8128,
                                    num_freqs = 16,
                                    num_list_values = 6)
        self.run_test_create_skymodel_chunk_map(settings)


    def test_P0000_G1000_S00_00_C1e8_Time014(self):

        settings = Skymodel_Settings(num_points = 0,
                                    num_gauss = 1000,
                                    num_shapes = 0,
                                    num_coeff_per_shape = 0,
                                    max_num_visibilities = 1e9,
                                    num_time_steps = 14,
                                    num_baselines = 8128,
                                    num_freqs = 16,
                                    num_list_values = 5)
        self.run_test_create_skymodel_chunk_map(settings)


    def test_P0000_G0000_S67_78_C1e8_Time014(self):

        settings = Skymodel_Settings(num_points = 0,
                                    num_gauss = 0,
                                    num_shapes = 67,
                                    num_coeff_per_shape = 78,
                                    max_num_visibilities = 1e9,
                                    num_time_steps = 14,
                                    num_baselines = 8128,
                                    num_freqs = 16,
                                    num_list_values = 2)
        self.run_test_create_skymodel_chunk_map(settings)


    def test_P0050_G0050_S10_25_C4576_Time014(self):

        settings = Skymodel_Settings(num_points = 50,
                                    num_gauss = 50,
                                    num_shapes = 10,
                                    num_coeff_per_shape = 25,
                                    max_num_visibilities = 4576,
                                    num_time_steps = 14,
                                    num_baselines = 8128,
                                    num_freqs = 16,
                                    num_list_values = 8)
        self.run_test_create_skymodel_chunk_map(settings)


    def test_P87760_G12207_S121_678_C1e10_Time056(self):

        settings = Skymodel_Settings(num_points = 87760,
                                    num_gauss = 12207,
                                    num_shapes = 121,
                                    num_coeff_per_shape = 678,
                                    max_num_visibilities = 1e10,
                                    num_time_steps = 56,
                                    num_baselines = 8128,
                                    num_freqs = 32,
                                    num_list_values = 16)
        self.run_test_create_skymodel_chunk_map(settings)
        
    def test_P87760_G12207_S121_678_C1e10_Time056_Lin70_V100(self):

        settings = Skymodel_Settings(num_points = 87760,
                                    num_gauss = 12207,
                                    num_shapes = 121,
                                    num_coeff_per_shape = 678,
                                    max_num_visibilities = 1e10,
                                    num_time_steps = 56,
                                    num_baselines = 8128,
                                    num_freqs = 32,
                                    num_list_values = 16,
                                    linpol_cadence = 70,
                                    linpol_num_list=2,
                                    linpol_num_p_list=3,
                                    stokesV_cadence = 80,
                                    stokesV_num_list=4)
        self.run_test_create_skymodel_chunk_map(settings)


        
##Run the test
if __name__ == '__main__':
    unittest.main()
    