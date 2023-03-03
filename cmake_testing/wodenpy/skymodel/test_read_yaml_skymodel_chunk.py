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
path.append('{:s}/../../../wodenpy/use_libwoden'.format(code_dir))
path.append('{:s}/../../../wodenpy'.format(code_dir))

# ##Code we are testing
import read_yaml_skymodel
# import wodenpy
from woden_skymodel import Component_Type_Counter, CompTypes, crop_below_horizon
from chunk_sky_model import create_skymodel_chunk_map, Skymodel_Chunk_Map

from skymodel_structs import setup_chunked_source, _Ctype_Source_Into_Python

from common_skymodel_test import fill_comp_counter, Expec_Counter, BaseChunkTest, write_full_test_skymodel

D2R = np.pi/180.0
MWA_LATITUDE = -26.7*D2R

##for now, WODEN has three flux types: power law, curved power law, and list
NUM_FLUX_TYPES = 3

##based class on an existing class test that includes methods for
##for 
class Test(BaseChunkTest):
    
    def run_test_read_yaml_skymodel_chunk(self, deg_between_comps : int,
                       num_coeff_per_shape : int, num_list_values : int,
                       comps_per_source : int, lst : float):
        """Let's go do a bit"""
        
        num_freqs = 16
        num_baselines = 8128
        num_time_steps = 14
        max_num_visibilities = int(1e7)
        
        write_full_test_skymodel(deg_between_comps,
                                 num_coeff_per_shape,
                                 num_list_values,
                                 comps_per_source)
        
        # filename = f"{code_dir}/test_full_skymodel.yaml"
        filename = f"test_full_skymodel.yaml"
        
        comp_counter = read_yaml_skymodel.read_yaml_radec_count_components(filename)
        
        
        comp_counter = crop_below_horizon(lst, MWA_LATITUDE,
                                          comp_counter, 
                                          crop_by_component=True)
        
        ##Create a chunking map
        chunked_skymodel_maps = create_skymodel_chunk_map(comp_counter,
                                        max_num_visibilities, num_baselines,
                                        num_freqs, num_time_steps)
        
        
        source_catalogue = read_yaml_skymodel.read_yaml_skymodel_chunks(filename, chunked_skymodel_maps,
                                                     num_freqs, num_time_steps)
        
        for source_ind in range(source_catalogue.num_sources):
            source = source_catalogue.sources[source_ind]
            
            print("new source who dis--------------------------")
            python_source = _Ctype_Source_Into_Python(source)
            
 
            if python_source.n_points:
                # print(python_source.point_components)
                
                print("point_ras", python_source.point_components.ras)
                print("point_decs", python_source.point_components.decs)
                
                print("point power_ref_freqs", python_source.point_components.power_ref_freqs)
                print("point power_ref_stokesI", python_source.point_components.power_ref_stokesI)
                print("point power_ref_stokesQ", python_source.point_components.power_ref_stokesQ)
                print("point power_ref_stokesU", python_source.point_components.power_ref_stokesU)
                print("point power_ref_stokesV", python_source.point_components.power_ref_stokesV)
                print("point power_SIs", python_source.point_components.power_SIs)
                
                print("point curve_ref_freqs", python_source.point_components.curve_ref_freqs)
                print("point curve_ref_stokesI", python_source.point_components.curve_ref_stokesI)
                print("point curve_ref_stokesQ", python_source.point_components.curve_ref_stokesQ)
                print("point curve_ref_stokesU", python_source.point_components.curve_ref_stokesU)
                print("point curve_ref_stokesV", python_source.point_components.curve_ref_stokesV)
                print("point curve_SIs", python_source.point_components.curve_SIs)
                print("point curve_qs", python_source.point_components.curve_qs)
                
                print("point list_freqs", python_source.point_components.list_freqs)
                print("point list_stokesI", python_source.point_components.list_stokesI)
                print("point list_stokesQ", python_source.point_components.list_stokesQ)
                print("point list_stokesU", python_source.point_components.list_stokesU)
                print("point list_stokesV", python_source.point_components.list_stokesV)
            
            
            if python_source.n_gauss:
                # print(python_source.shape_components)
                
                print("gauss_ras", python_source.gauss_components.ras)
                print("gauss_decs", python_source.gauss_components.decs)
                
                print("gauss power_ref_freqs", python_source.gauss_components.power_ref_freqs)
                print("gauss power_ref_stokesI", python_source.gauss_components.power_ref_stokesI)
                print("gauss power_ref_stokesQ", python_source.gauss_components.power_ref_stokesQ)
                print("gauss power_ref_stokesU", python_source.gauss_components.power_ref_stokesU)
                print("gauss power_ref_stokesV", python_source.gauss_components.power_ref_stokesV)
                print("gauss power_SIs", python_source.gauss_components.power_SIs)
                
                print("gauss curve_ref_freqs", python_source.gauss_components.curve_ref_freqs)
                print("gauss curve_ref_stokesI", python_source.gauss_components.curve_ref_stokesI)
                print("gauss curve_ref_stokesQ", python_source.gauss_components.curve_ref_stokesQ)
                print("gauss curve_ref_stokesU", python_source.gauss_components.curve_ref_stokesU)
                print("gauss curve_ref_stokesV", python_source.gauss_components.curve_ref_stokesV)
                print("gauss curve_SIs", python_source.gauss_components.curve_SIs)
                print("gauss curve_qs", python_source.gauss_components.curve_qs)
                
                print("gauss list_freqs", python_source.gauss_components.list_freqs)
                print("gauss list_stokesI", python_source.gauss_components.list_stokesI)
                print("gauss list_stokesQ", python_source.gauss_components.list_stokesQ)
                print("gauss list_stokesU", python_source.gauss_components.list_stokesU)
                print("gauss list_stokesV", python_source.gauss_components.list_stokesV)
                
                print("gauss majors", python_source.gauss_components.majors * (3600 / D2R))
                print("gauss minors", python_source.gauss_components.minors * (3600 / D2R))
                print("gauss pas", python_source.gauss_components.pas / D2R)
                
            if python_source.n_shapes:
                # print(python_source.shape_components)
                
                print("shape_ras", python_source.shape_components.ras)
                print("shape_decs", python_source.shape_components.decs)
                
                print("shape power_ref_freqs", python_source.shape_components.power_ref_freqs)
                print("shape power_ref_stokesI", python_source.shape_components.power_ref_stokesI)
                print("shape power_ref_stokesQ", python_source.shape_components.power_ref_stokesQ)
                print("shape power_ref_stokesU", python_source.shape_components.power_ref_stokesU)
                print("shape power_ref_stokesV", python_source.shape_components.power_ref_stokesV)
                print("shape power_SIs", python_source.shape_components.power_SIs)
                
                print("shape curve_ref_freqs", python_source.shape_components.curve_ref_freqs)
                print("shape curve_ref_stokesI", python_source.shape_components.curve_ref_stokesI)
                print("shape curve_ref_stokesQ", python_source.shape_components.curve_ref_stokesQ)
                print("shape curve_ref_stokesU", python_source.shape_components.curve_ref_stokesU)
                print("shape curve_ref_stokesV", python_source.shape_components.curve_ref_stokesV)
                print("shape curve_SIs", python_source.shape_components.curve_SIs)
                print("shape curve_qs", python_source.shape_components.curve_qs)
                
                print("shape list_freqs", python_source.shape_components.list_freqs)
                print("shape list_stokesI", python_source.shape_components.list_stokesI)
                print("shape list_stokesQ", python_source.shape_components.list_stokesQ)
                print("shape list_stokesU", python_source.shape_components.list_stokesU)
                print("shape list_stokesV", python_source.shape_components.list_stokesV)
                
                print("shape list_freqs", python_source.shape_components.list_freqs)
                print("shape list_stokesI", python_source.shape_components.list_stokesI)
                print("shape list_stokesQ", python_source.shape_components.list_stokesQ)
                print("shape list_stokesU", python_source.shape_components.list_stokesU)
                print("shape list_stokesV", python_source.shape_components.list_stokesV)
                
                print("shape majors", python_source.shape_components.majors * (3600 / D2R))
                print("shape minors", python_source.shape_components.minors * (3600 / D2R))
                print("shape pas", python_source.shape_components.pas / D2R)
                
                print("shape param_indexes", python_source.shape_components.param_indexes)
                print("shape n1s", python_source.shape_components.n1s)
                print("shape n2s", python_source.shape_components.n2s)
                print("shape shape_coeffs", python_source.shape_components.shape_coeffs)
                # print("shape minors", python_source.shape_components.minors * (3600 / D2R))
                # print("shape pas", python_source.shape_components.pas / D2R)
                
        # for precision in ['double']:
            
            
                
            
            
            # print(chunked_source.point_components.ras[:1])
            # print(chunked_source.point_components.decs[:1])
            
            # self.assertTrue(np.array_equal(np.array([20.0]), chunked_source.point_components.ras[:1]))
            
            
            
            # np.ctypeslib.as_array(visibility_set.us_metres, shape=(num_visi,))
            
            
            
        
        # comps_per_chunk = int(np.floor(max_num_visibilities / (num_baselines * num_freqs * num_time_steps)))
        
        # ##pray this never happens, probably means we're going to run out of
        # ##GPU memory TODO don't pray, submit a warning?
        # if comps_per_chunk < 1: comps_per_chunk = 1
        
        # num_point_chunks = int(np.ceil(NUM_FLUX_TYPES*num_points / comps_per_chunk))
        # num_gauss_chunks = int(np.ceil(NUM_FLUX_TYPES*num_gauss / comps_per_chunk))
        
        # ##we split SHAPELET by the basis components (number of coeffs)
        # num_coeff_chunks = int(np.ceil(comp_counter.total_shape_basis / comps_per_chunk))
        
        # ##Counters to see how far many point/gauss have already been assigned
        # ##to previous chunks
        # expec_counter = Expec_Counter()
        
        # for point_ind in range(num_point_chunks):
            
        #     chunk_map = chunked_skymodel_counters[point_ind]
            
        #     expec_counter = self.check_pointgauss_chunking(point_ind,
        #                                     comps_per_chunk,
        #                                     num_point_chunks, num_list_values,
        #                                     num_points,
        #                                     CompTypes.POINT,
        #                                     comp_counter, chunk_map,
        #                                     expec_counter)
            
        # for gauss_ind in range(num_gauss_chunks):
            
        #     chunk_map = chunked_skymodel_counters[num_point_chunks + gauss_ind]
            
        #     expec_counter = self.check_pointgauss_chunking(gauss_ind,
        #                                     comps_per_chunk,
        #                                     num_gauss_chunks, num_list_values,
        #                                     num_gauss,
        #                                     CompTypes.GAUSSIAN,
        #                                     comp_counter, chunk_map,
        #                                     expec_counter)
            
            
        # for coeff_ind in range(num_coeff_chunks):
            
        #     chunk_map = chunked_skymodel_counters[num_point_chunks + num_gauss_chunks + coeff_ind]
            
        #     self.check_shapelet_chunking(coeff_ind, num_coeff_per_shape,
        #                                  comps_per_chunk,
        #                                  num_list_values,  num_shapes,
        #                                  comp_counter, chunk_map,
        #                                  total_point_comps=NUM_FLUX_TYPES*num_points,
        #                                  total_gauss_comps=NUM_FLUX_TYPES*num_gauss)
            
    ##Ok now run with many many combinations
    
    def test_the_first(self):
        deg_between_comps = 90
        num_coeff_per_shape = 4
        num_list_values = 4
        comps_per_source = 10
        lst = 0.0
        
        self.run_test_read_yaml_skymodel_chunk(deg_between_comps,
                                               num_coeff_per_shape,
                                               num_list_values,
                                               comps_per_source, lst)
        


    # def test_P4544_G1736_S20_14_C1e8_Time014(self):

    #     num_points = 4544
    #     num_gauss = 1736
    #     num_shapes = 20
    #     num_coeff_per_shape = 14

    #     max_num_visibilities = 1e8

    #     num_time_steps = 14
    #     num_baselines = 8128
    #     num_freqs = 16

    #     num_list_values = 4
    #     self.run_test_read_yaml_skymodel_chunk(max_num_visibilities,
    #                                         num_points, num_gauss,
    #                                         num_shapes, num_coeff_per_shape,
    #                                         num_list_values, num_time_steps,
    #                                         num_baselines, num_freqs)


    # def test_P1000_G0000_S00_00_C1e8_Time014(self):

    #     num_points = 1000
    #     num_gauss = 0
    #     num_shapes = 0
    #     num_coeff_per_shape = 0

    #     max_num_visibilities = 1e9

    #     num_time_steps = 14
    #     num_baselines = 8128
    #     num_freqs = 16

    #     num_list_values = 6
    #     self.run_test_read_yaml_skymodel_chunk(max_num_visibilities,
    #                                         num_points, num_gauss,
    #                                         num_shapes, num_coeff_per_shape,
    #                                         num_list_values, num_time_steps,
    #                                         num_baselines, num_freqs)


    # def test_P0000_G1000_S00_00_C1e8_Time014(self):

    #     num_points = 0
    #     num_gauss = 1000
    #     num_shapes = 0
    #     num_coeff_per_shape = 0

    #     max_num_visibilities = 1e9

    #     num_time_steps = 14
    #     num_baselines = 8128
    #     num_freqs = 16

    #     num_list_values = 5
    #     self.run_test_read_yaml_skymodel_chunk(max_num_visibilities,
    #                                         num_points, num_gauss,
    #                                         num_shapes, num_coeff_per_shape,
    #                                         num_list_values, num_time_steps,
    #                                         num_baselines, num_freqs)


    # def test_P0000_G0000_S67_78_C1e8_Time014(self):

    #     num_points = 0
    #     num_gauss = 0
    #     num_shapes = 67
    #     num_coeff_per_shape = 78

    #     max_num_visibilities = 1e9

    #     num_time_steps = 14
    #     num_baselines = 8128
    #     num_freqs = 16

    #     num_list_values = 2
    #     self.run_test_read_yaml_skymodel_chunk(max_num_visibilities,
    #                                         num_points, num_gauss,
    #                                         num_shapes, num_coeff_per_shape,
    #                                         num_list_values, num_time_steps,
    #                                         num_baselines, num_freqs)


    # def test_P0050_G0050_S10_25_C4576_Time014(self):

    #     num_points = 50
    #     num_gauss = 50
    #     num_shapes = 10
    #     num_coeff_per_shape = 25

    #     max_num_visibilities = 4576

    #     num_time_steps = 14
    #     num_baselines = 8128
    #     num_freqs = 16

    #     num_list_values = 8
    #     self.run_test_read_yaml_skymodel_chunk(max_num_visibilities,
    #                                         num_points, num_gauss,
    #                                         num_shapes, num_coeff_per_shape,
    #                                         num_list_values, num_time_steps,
    #                                         num_baselines, num_freqs)


    # def test_P87760_G12207_S121_678_C1e10_Time056(self):

    #     num_points = 87760
    #     num_gauss = 12207
    #     num_shapes = 121
    #     num_coeff_per_shape = 678

    #     max_num_visibilities = 1e10

    #     num_time_steps = 56
    #     num_baselines = 8128
    #     num_freqs = 32

    #     num_list_values = 16
    #     self.run_test_read_yaml_skymodel_chunk(max_num_visibilities,
    #                                         num_points, num_gauss,
    #                                         num_shapes, num_coeff_per_shape,
    #                                         num_list_values, num_time_steps,
    #                                         num_baselines, num_freqs)


        
##Run the test
if __name__ == '__main__':
    unittest.main()
    