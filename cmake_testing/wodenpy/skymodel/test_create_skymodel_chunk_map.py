"""

"""

from sys import path
import os
import unittest
import numpy as np

##This is where our code lives
code_dir = os.environ['CMAKE_CURRENT_SOURCE_DIR']

##Append the location of run_woden.py to the sys.path to import it
path.append('{:s}/../../../wodenpy/skymodel'.format(code_dir))
path.append('{:s}/../../../wodenpy'.format(code_dir))

# ##Code we are testing
import read_yaml_skymodel
# import wodenpy
from woden_skymodel import Component_Type_Counter, CompTypes
from chunk_sky_model import create_skymodel_chunk_map, Skymodel_Chunk_Map

from common_skymodel_test import fill_comp_counter, Expec_Counter, BaseChunkTest

D2R = np.pi/180.0

##for now, WODEN has three flux types: power law, curved power law, and list
NUM_FLUX_TYPES = 3

##based class on an existing class test that includes methods for
##for 
class Test(BaseChunkTest):
    
    def run_test_create_skymodel_chunk_map(self, max_num_visibilities : int,
                            num_points : int, num_gauss : int,
                            num_shapes : int, num_coeff_per_shape : int,
                            num_list_values : int, num_time_steps : int,
                            num_baselines : int, num_freqs : int):
        """Makes a populted Component_Type_Counter based off the given
        inputs, runs it through
        `wodenpy.skymodel.chunk_sky_model.map_chunk_pointgauss` and
        checks that gives the correct answers"""
        
        ##make the fake sky model        
        comp_counter = fill_comp_counter(num_points, num_gauss,
                                         num_shapes, num_coeff_per_shape,
                                         num_list_values, num_time_steps)
        
        
        ##Run the code we are testing!
        chunked_skymodel_counters = create_skymodel_chunk_map(comp_counter,
                                        max_num_visibilities, num_baselines,
                                        num_freqs, num_time_steps)
        
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
            
            chunk_map = chunked_skymodel_counters[point_ind]
            
            expec_counter = self.check_pointgauss_chunking(point_ind,
                                            comps_per_chunk,
                                            num_point_chunks, num_list_values,
                                            num_points,
                                            CompTypes.POINT,
                                            comp_counter, chunk_map,
                                            expec_counter)
            
        for gauss_ind in range(num_gauss_chunks):
            
            chunk_map = chunked_skymodel_counters[num_point_chunks + gauss_ind]
            
            expec_counter = self.check_pointgauss_chunking(gauss_ind,
                                            comps_per_chunk,
                                            num_gauss_chunks, num_list_values,
                                            num_gauss,
                                            CompTypes.GAUSSIAN,
                                            comp_counter, chunk_map,
                                            expec_counter)
            
            
        for coeff_ind in range(num_coeff_chunks):
            
            chunk_map = chunked_skymodel_counters[num_point_chunks + num_gauss_chunks + coeff_ind]
            
            self.check_shapelet_chunking(coeff_ind, num_coeff_per_shape,
                                         comps_per_chunk,
                                         num_list_values,  num_shapes,
                                         comp_counter, chunk_map,
                                         total_point_comps=NUM_FLUX_TYPES*num_points,
                                         total_gauss_comps=NUM_FLUX_TYPES*num_gauss)
            
    ##Ok now run with many many combinations
    
    def test_P3_G0000_S00_00_C1_Time001(self):
        max_num_visibilities = 2
        
        num_points = 3
        num_gauss = 0
        num_shapes = 0
        num_coeff_per_shape = 0

        num_time_steps = 1
        num_baselines = 1
        num_freqs = 1

        num_list_values = 3
        
        self.run_test_create_skymodel_chunk_map(max_num_visibilities,
                                            num_points, num_gauss,
                                            num_shapes, num_coeff_per_shape,
                                            num_list_values, num_time_steps,
                                            num_baselines, num_freqs)
        


    def test_P4544_G1736_S20_14_C1e8_Time014(self):

        num_points = 4544
        num_gauss = 1736
        num_shapes = 20
        num_coeff_per_shape = 14

        max_num_visibilities = 1e8

        num_time_steps = 14
        num_baselines = 8128
        num_freqs = 16

        num_list_values = 4
        self.run_test_create_skymodel_chunk_map(max_num_visibilities,
                                            num_points, num_gauss,
                                            num_shapes, num_coeff_per_shape,
                                            num_list_values, num_time_steps,
                                            num_baselines, num_freqs)


    def test_P1000_G0000_S00_00_C1e8_Time014(self):

        num_points = 1000
        num_gauss = 0
        num_shapes = 0
        num_coeff_per_shape = 0

        max_num_visibilities = 1e9

        num_time_steps = 14
        num_baselines = 8128
        num_freqs = 16

        num_list_values = 6
        self.run_test_create_skymodel_chunk_map(max_num_visibilities,
                                            num_points, num_gauss,
                                            num_shapes, num_coeff_per_shape,
                                            num_list_values, num_time_steps,
                                            num_baselines, num_freqs)


    def test_P0000_G1000_S00_00_C1e8_Time014(self):

        num_points = 0
        num_gauss = 1000
        num_shapes = 0
        num_coeff_per_shape = 0

        max_num_visibilities = 1e9

        num_time_steps = 14
        num_baselines = 8128
        num_freqs = 16

        num_list_values = 5
        self.run_test_create_skymodel_chunk_map(max_num_visibilities,
                                            num_points, num_gauss,
                                            num_shapes, num_coeff_per_shape,
                                            num_list_values, num_time_steps,
                                            num_baselines, num_freqs)


    def test_P0000_G0000_S67_78_C1e8_Time014(self):

        num_points = 0
        num_gauss = 0
        num_shapes = 67
        num_coeff_per_shape = 78

        max_num_visibilities = 1e9

        num_time_steps = 14
        num_baselines = 8128
        num_freqs = 16

        num_list_values = 2
        self.run_test_create_skymodel_chunk_map(max_num_visibilities,
                                            num_points, num_gauss,
                                            num_shapes, num_coeff_per_shape,
                                            num_list_values, num_time_steps,
                                            num_baselines, num_freqs)


    def test_P0050_G0050_S10_25_C4576_Time014(self):

        num_points = 50
        num_gauss = 50
        num_shapes = 10
        num_coeff_per_shape = 25

        max_num_visibilities = 4576

        num_time_steps = 14
        num_baselines = 8128
        num_freqs = 16

        num_list_values = 8
        self.run_test_create_skymodel_chunk_map(max_num_visibilities,
                                            num_points, num_gauss,
                                            num_shapes, num_coeff_per_shape,
                                            num_list_values, num_time_steps,
                                            num_baselines, num_freqs)


    def test_P87760_G12207_S121_678_C1e10_Time056(self):

        num_points = 87760
        num_gauss = 12207
        num_shapes = 121
        num_coeff_per_shape = 678

        max_num_visibilities = 1e10

        num_time_steps = 56
        num_baselines = 8128
        num_freqs = 32

        num_list_values = 16
        self.run_test_create_skymodel_chunk_map(max_num_visibilities,
                                            num_points, num_gauss,
                                            num_shapes, num_coeff_per_shape,
                                            num_list_values, num_time_steps,
                                            num_baselines, num_freqs)


        
##Run the test
if __name__ == '__main__':
    unittest.main()
    