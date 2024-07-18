"""

"""

from sys import path
import os
import unittest
import numpy as np

from wodenpy.skymodel.woden_skymodel import Component_Type_Counter, CompTypes
from wodenpy.skymodel.chunk_sky_model import map_chunk_pointgauss, Skymodel_Chunk_Map

from common_skymodel_test import fill_comp_counter_for_chunking, Expec_Counter, BaseChunkTest

D2R = np.pi/180.0

##for now, WODEN has three flux types: power law, curved power law, and list
NUM_FLUX_TYPES = 3

class Test(BaseChunkTest):
    
    def run_test_map_chunk_pointgauss(self, max_num_visibilities : int,
                            num_points : int, num_gauss : int,
                            num_list_values : int, num_time_steps : int,
                            num_baselines : int, num_freqs : int,
                            stokesV_cadence : int = 0,
                            linpol_cadence : int = 0):
        """Makes a populted Component_Type_Counter based off the given
        inputs, runs it through
        `wodenpy.skymodel.chunk_sky_model.map_chunk_pointgauss` and
        checks that gives the correct answers"""
        
        ## Not testing SHAPELETs, set these to zero
        num_shapes = 0
        num_coeff_per_shape = 0
                
        full_comp_counter = fill_comp_counter_for_chunking(num_points, num_gauss,
                                            num_shapes, num_coeff_per_shape,
                                            num_list_values, num_time_steps,
                                            stokesV_cadence=stokesV_cadence,
                                            linpol_cadence=linpol_cadence)
        
        full_comp_counter.total_components()
        
        comps_per_chunk = int(np.floor(max_num_visibilities / (num_baselines * num_freqs * num_time_steps)))
        
        num_point_chunks = int(np.ceil(NUM_FLUX_TYPES*num_points / comps_per_chunk))
        num_gauss_chunks = int(np.ceil(NUM_FLUX_TYPES*num_gauss / comps_per_chunk))
        
        ##Counters to see how far many point/gauss have already been assigned
        ##to previous chunks
        expec_counter = Expec_Counter()
        
        for chunk_ind in range(num_point_chunks):
            chunk_map = map_chunk_pointgauss(full_comp_counter,
                                             chunk_ind, comps_per_chunk,
                                             point_source = True)
            
            expec_counter = self.check_pointgauss_chunking(chunk_ind,
                                            comps_per_chunk,
                                            num_point_chunks, num_list_values,
                                            num_points,
                                            linpol_cadence, stokesV_cadence,
                                            CompTypes.POINT,
                                            full_comp_counter, chunk_map,
                                            expec_counter)
            
        for chunk_ind in range(num_gauss_chunks):
            chunk_map = map_chunk_pointgauss(full_comp_counter,
                                                chunk_ind, comps_per_chunk,
                                                gaussian_source = True)
            
            expec_counter = self.check_pointgauss_chunking(chunk_ind,
                                            comps_per_chunk,
                                            num_point_chunks, num_list_values,
                                            num_gauss,
                                            linpol_cadence, stokesV_cadence,
                                            CompTypes.GAUSSIAN,
                                            full_comp_counter, chunk_map,
                                            expec_counter)
            
    ##Ok now run with many many combinations
    
    def test_Point100_Gauss000_Chunk1000_Time004(self):
        max_num_visibilities = 1000
        num_points = 100
        num_gauss = 0
        num_time_steps = 4
        num_baselines = 5
        num_freqs = 2
        num_list_values = 2
        
        self.run_test_map_chunk_pointgauss(max_num_visibilities,
                                           num_points,
                                           num_gauss, num_list_values,
                                           num_time_steps, num_baselines,
                                           num_freqs)
        
    def test_Point000_Gauss100_Chunk1000_Time004(self):
        max_num_visibilities = 1000
        num_points = 0
        num_gauss = 100
        num_time_steps = 4
        num_baselines = 5
        num_freqs = 2

        num_list_values = 3
        self.run_test_map_chunk_pointgauss(max_num_visibilities,
                                           num_points,
                                           num_gauss, num_list_values,
                                           num_time_steps, num_baselines,
                                           num_freqs)

    def test_Point100_Gauss100_Chunk1000_Time004(self):
        max_num_visibilities = 1000    
        num_points = 100    
        num_gauss = 100    
        num_time_steps = 4    
        num_baselines = 5    
        num_freqs = 2    

        num_list_values = 4
    
        self.run_test_map_chunk_pointgauss(max_num_visibilities,
                                           num_points,
                                           num_gauss, num_list_values,
                                           num_time_steps, num_baselines,
                                           num_freqs)
    

    def test_Point100_Gauss000_Chunk1000_Time003(self):
        max_num_visibilities = 1000    
        num_points = 100    
        num_gauss = 0    
        num_time_steps = 3    
        num_baselines = 5    
        num_freqs = 2    

        num_list_values = 5    
        self.run_test_map_chunk_pointgauss(max_num_visibilities,
                                           num_points,
                                           num_gauss, num_list_values,
                                           num_time_steps, num_baselines,
                                           num_freqs)   
    

    def test_Point000_Gauss100_Chunk1000_Time003(self):
        max_num_visibilities = 1000    
        num_points = 0    
        num_gauss = 100    
        num_time_steps = 3    
        num_baselines = 5    
        num_freqs = 2    

        num_list_values = 6    
        self.run_test_map_chunk_pointgauss(max_num_visibilities,
                                           num_points,
                                           num_gauss, num_list_values,
                                           num_time_steps, num_baselines,
                                           num_freqs)   
    

    def test_Point100_Gauss100_Chunk1000_Time003(self):
        max_num_visibilities = 1000    
        num_points = 100    
        num_gauss = 100    
        num_time_steps = 3    
        num_baselines = 5    
        num_freqs = 2    

        num_list_values = 5    
        self.run_test_map_chunk_pointgauss(max_num_visibilities,
                                           num_points,
                                           num_gauss, num_list_values,
                                           num_time_steps, num_baselines,
                                           num_freqs)   
    

    def test_Point100_Gauss000_Chunk0173_Time005(self):
        max_num_visibilities = 173    
        num_points = 100    
        num_gauss = 0    
        num_time_steps = 5    
        num_baselines = 5    
        num_freqs = 2    

        num_list_values = 4    
        self.run_test_map_chunk_pointgauss(max_num_visibilities,
                                           num_points,
                                           num_gauss, num_list_values,
                                           num_time_steps, num_baselines,
                                           num_freqs)   
    

    def test_Point000_Gauss100_Chunk0173_Time005(self):
        max_num_visibilities = 173    
        num_points = 0    
        num_gauss = 100    
        num_time_steps = 5    
        num_baselines = 5    
        num_freqs = 2    

        num_list_values = 3    
        self.run_test_map_chunk_pointgauss(max_num_visibilities,
                                           num_points,
                                           num_gauss, num_list_values,
                                           num_time_steps, num_baselines,
                                           num_freqs)   
    

    def test_Point100_Gauss100_Chunk0173_Time005(self):
        max_num_visibilities = 173    
        num_points = 100    
        num_gauss = 100    
        num_time_steps = 5    
        num_baselines = 5    
        num_freqs = 2    

        num_list_values = 2    
        self.run_test_map_chunk_pointgauss(max_num_visibilities,
                                           num_points,
                                           num_gauss, num_list_values,
                                           num_time_steps, num_baselines,
                                           num_freqs)   
    

    def test_Point10743_Gauss00000_Chunk4345_Time012(self):
        max_num_visibilities = 4345    
        num_points = 10743    
        num_gauss = 0    
        num_time_steps = 12    
        num_baselines = 5    
        num_freqs = 2    

        num_list_values = 16    
        self.run_test_map_chunk_pointgauss(max_num_visibilities,
                                           num_points,
                                           num_gauss, num_list_values,
                                           num_time_steps, num_baselines,
                                           num_freqs)   
    

    def test_Point00000_Gauss16789_Chunk4345_Time012(self):
        max_num_visibilities = 4345    
        num_points = 0    
        num_gauss = 16789    
        num_time_steps = 12    
        num_baselines = 5    
        num_freqs = 2    

        num_list_values = 20    
        self.run_test_map_chunk_pointgauss(max_num_visibilities,
                                           num_points,
                                           num_gauss, num_list_values,
                                           num_time_steps, num_baselines,
                                           num_freqs)   
    

    def test_Point10743_Gauss16789_Chunk4345_Time012(self):
        max_num_visibilities = 4345    
        num_points = 10743    
        num_gauss = 16789    
        num_time_steps = 12    
        num_baselines = 5    
        num_freqs = 2    

        num_list_values = 10    
        self.run_test_map_chunk_pointgauss(max_num_visibilities,
                                           num_points,
                                           num_gauss, num_list_values,
                                           num_time_steps, num_baselines,
                                           num_freqs)   
    

    def test_Point4544_Gauss1736_Chunk1e8_Time014(self):
        max_num_visibilities = 1e8    
        num_points = 4544    
        num_gauss = 1736    
        num_time_steps = 14    
        num_baselines = 8128    
        num_freqs = 16    

        num_list_values = 12    
    
        self.run_test_map_chunk_pointgauss(max_num_visibilities,
                                           num_points,
                                           num_gauss, num_list_values,
                                           num_time_steps, num_baselines,
                                           num_freqs) 
        
    def test_Point100_Gauss100_Chunk0173_Time005_Vcadence15(self):
        max_num_visibilities = 173    
        num_points = 100    
        num_gauss = 100    
        num_time_steps = 5    
        num_baselines = 5    
        num_freqs = 2
        num_list_values = 2
        stokesV_cadence = 15
        self.run_test_map_chunk_pointgauss(max_num_visibilities,
                                           num_points,
                                           num_gauss, num_list_values,
                                           num_time_steps, num_baselines,
                                           num_freqs,
                                           stokesV_cadence=stokesV_cadence)
        
    def test_Point100_Gauss100_Chunk0173_Time005_Lincadence20(self):
        max_num_visibilities = 173    
        num_points = 100    
        num_gauss = 100    
        num_time_steps = 5    
        num_baselines = 5    
        num_freqs = 2
        num_list_values = 2
        linpol_cadence = 20
        self.run_test_map_chunk_pointgauss(max_num_visibilities,
                                           num_points,
                                           num_gauss, num_list_values,
                                           num_time_steps, num_baselines,
                                           num_freqs,
                                           linpol_cadence=linpol_cadence)
        
    def test_Point100_Gauss100_Chunk0173_Time005_Lincadence10_Vcadence7(self):
        max_num_visibilities = 173    
        num_points = 100    
        num_gauss = 100    
        num_time_steps = 5    
        num_baselines = 5    
        num_freqs = 2
        num_list_values = 2
        linpol_cadence = 10
        stokesV_cadence = 7
        self.run_test_map_chunk_pointgauss(max_num_visibilities,
                                           num_points,
                                           num_gauss, num_list_values,
                                           num_time_steps, num_baselines,
                                           num_freqs,
                                           linpol_cadence=linpol_cadence,
                                           stokesV_cadence=stokesV_cadence)
        
    def test_Point4544_Gauss1736_Chunk1e8_Time014_Lincadence4_Vcadence3(self):
        max_num_visibilities = 1e8    
        num_points = 4544    
        num_gauss = 1736    
        num_time_steps = 14    
        num_baselines = 8128    
        num_freqs = 16    
        num_list_values = 12    
        linpol_cadence = 4
        stokesV_cadence = 3
    
        self.run_test_map_chunk_pointgauss(max_num_visibilities,
                                           num_points,
                                           num_gauss, num_list_values,
                                           num_time_steps, num_baselines,
                                           num_freqs,
                                           linpol_cadence=linpol_cadence,
                                           stokesV_cadence=stokesV_cadence) 
        
##Run the test
if __name__ == '__main__':
    unittest.main()
    