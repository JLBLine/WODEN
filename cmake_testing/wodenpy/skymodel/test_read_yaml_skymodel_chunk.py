"""

"""

from sys import path
import os
import unittest
import numpy as np
import erfa

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
from chunk_sky_model import create_skymodel_chunk_map, Skymodel_Chunk_Map, increment_flux_type_counters
from beam_settings import BeamTypes

from skymodel_structs import setup_chunked_source, _Ctype_Source_Into_Python

from common_skymodel_test import fill_comp_counter, Expec_Counter, BaseChunkTest, Expected_Sky_Chunk, Expected_Components

import woden_settings as ws

from read_skymodel_common import check_components, check_all_sources, populate_pointgauss_chunk, populate_shapelet_chunk, make_expected_chunks


D2R = np.pi/180.0
# MWA_LATITUDE = -26.7*D2R

MWA_LAT = -26.703319405555554*D2R

##for now, WODEN has three flux types: power law, curved power law, and list
NUM_FLUX_TYPES = 3

D2R = np.pi/180.0

##for now, WODEN has three flux types: power law, curved power law, and list
NUM_FLUX_TYPES = 3

##limits for "all sky" sky model
LOW_DEC = -90.0*D2R
HIGH_DEC = 30.0*D2R

RTOL=1e-10


def add_power_law_yaml(outfile, value):
    """add power law infor for a component"""
    outfile.write(f"    flux_type:\n")
    outfile.write(f"      power_law:\n")
    outfile.write(f"        si: {float(value)/100.0}\n")
    outfile.write(f"        fd:\n")
    freq = 100e+6 + (value + 1)*1e+4
    outfile.write(f"          freq: {freq:.1f}\n")
    outfile.write(f"          i: {float(value)}\n")
    outfile.write(f"          q: {float(value)}\n")
    outfile.write(f"          u: {float(value)}\n")
    outfile.write(f"          v: {float(value)}\n")
    
def add_curved_power_law_yaml(outfile, value):
    """add curved power law infor for a component"""
    outfile.write(f"    flux_type:\n")
    outfile.write(f"      curved_power_law:\n")
    outfile.write(f"        si: {float(value)/100.0}\n")
    outfile.write(f"        fd:\n")
    freq = 100e+6 + (value + 1)*1e+4
    outfile.write(f"          freq: {freq:.1f}\n")
    outfile.write(f"          i: {float(value)}\n")
    outfile.write(f"          q: {float(value)}\n")
    outfile.write(f"          u: {float(value)}\n")
    outfile.write(f"          v: {float(value)}\n")
    outfile.write(f"        q: {float(value)/100.0}\n")
    
def add_list_flux_yaml(outfile, flux_index, num_list_values):
    """add curved power law infor for a component"""
    outfile.write(f"    flux_type:\n")
    outfile.write(f"      list:\n")
    
    for freq in np.arange(num_list_values)*1e+6:
        # outfile.write(f"         - freq: {float(flux_index)}\n")
        outfile.write(f"         - freq: {float(freq)}\n")
        outfile.write(f"           i: {float(flux_index)}\n")
        outfile.write(f"           q: {float(flux_index)}\n")
        outfile.write(f"           u: {float(flux_index)}\n")
        outfile.write(f"           v: {float(flux_index)}\n")
        
        flux_index += 1
        
    return flux_index
    

def add_point_yaml(outfile, ra : float, dec: float,
              point_index : int, gauss_index : int, shape_index : int,
              comp_type : CompTypes,
              flux_index : int, num_list_values : int,
              source_index : int, comps_per_source : int):
    
    comp_index = point_index + gauss_index + shape_index
    if comp_index % comps_per_source == 0:
        new_source = True
    else:
        new_source = False
    
    if new_source:
        outfile.write(f"source{source_index:08d}\n")
        
    outfile.write(f"  - ra: {ra/D2R:.5f}\n")
    outfile.write(f"    dec: {dec/D2R:.5f}\n")
    outfile.write(f"    comp_type: point\n")
    
    if comp_type == CompTypes.POINT_POWER:
        add_power_law_yaml(outfile, point_index)
        
    elif comp_type == CompTypes.POINT_CURVE:
        add_curved_power_law_yaml(outfile, point_index)
        
    elif comp_type == CompTypes.POINT_LIST:
        flux_index = add_list_flux_yaml(outfile, flux_index, num_list_values)
        
    if new_source: source_index += 1
    point_index += 1
        
    return source_index, point_index, flux_index

def add_gauss_yaml(outfile, ra : float, dec: float,
                    point_index : int, gauss_index : int, shape_index : int,
                    comp_type : CompTypes,
                    flux_index : int, num_list_values : int,
                    source_index : int, comps_per_source : int):
    
    comp_index = point_index + gauss_index + shape_index
    if comp_index % comps_per_source == 0:
        new_source = True
    else:
        new_source = False
    
    if new_source:
        outfile.write(f"source{source_index:08d}\n")
        
    outfile.write(f"  - ra: {ra/D2R:.12f}\n")
    outfile.write(f"    dec: {dec/D2R:.12f}\n")
    outfile.write(f"    comp_type:\n")
    outfile.write(f"      gaussian:\n")
    outfile.write(f"        maj: {gauss_index}.\n")
    outfile.write(f"        min: {gauss_index}.\n")
    outfile.write(f"        pa: {gauss_index}.\n")
    
    if comp_type == CompTypes.GAUSS_POWER:
        add_power_law_yaml(outfile, gauss_index)
        
    elif comp_type == CompTypes.GAUSS_CURVE:
        add_curved_power_law_yaml(outfile, gauss_index)
        
    elif comp_type == CompTypes.GAUSS_LIST:
        flux_index = add_list_flux_yaml(outfile, flux_index, num_list_values)
        
    if new_source: source_index += 1
    gauss_index += 1
        
    return source_index, gauss_index, flux_index
    

def add_shapelet_yaml(outfile, ra : float, dec: float,
                    point_index : int, gauss_index : int, shape_index : int,
                    comp_type : CompTypes,
                    flux_index : int, num_list_values : int,
                    source_index : int, comps_per_source : int,
                    num_coeff_per_shape : int, basis_index : int):
    
    comp_index = point_index + gauss_index + shape_index
    if comp_index % comps_per_source == 0:
        new_source = True
    else:
        new_source = False
    
    if new_source:
        outfile.write(f"source{source_index:08d}\n")
        
    outfile.write(f"  - ra: {ra/D2R:.12f}\n")
    outfile.write(f"    dec: {dec/D2R:.12f}\n")
    outfile.write(f"    comp_type:\n")
    outfile.write(f"      shapelet:\n")
    outfile.write(f"        maj: {shape_index}.\n")
    outfile.write(f"        min: {shape_index}.\n")
    outfile.write(f"        pa: {shape_index}.\n")
    outfile.write(f"        coeffs:\n")
    
    for basis in range(basis_index, basis_index + num_coeff_per_shape):
    
        outfile.write(f"          - n1: {basis}\n")
        outfile.write(f"            n2: {basis}\n")
        outfile.write(f"            value: {basis}\n")
    
    basis_index += num_coeff_per_shape
    
    if comp_type == CompTypes.SHAPE_POWER:
        add_power_law_yaml(outfile, shape_index)
        
    elif comp_type == CompTypes.SHAPE_CURVE:
        add_curved_power_law_yaml(outfile, shape_index)
        
    elif comp_type == CompTypes.SHAPE_LIST:
        flux_index = add_list_flux_yaml(outfile, flux_index, num_list_values)
        
    if new_source: source_index += 1
    shape_index += 1
        
    return source_index, shape_index, flux_index, basis_index


def write_full_test_skymodel_yaml(deg_between_comps : float,
                             num_coeff_per_shape : int,
                             num_list_values : int,
                             comps_per_source : int):
    """Write a sky model covering the whole sky"""
    # POINT_POWER
    # POINT_CURVE
    # POINT_LIST
    # GAUSS_POWER
    # GAUSS_CURVE
    # GAUSS_LIST
    # SHAPE_POWER
    # SHAPE_CURVE
    # SHAPE_LIST
    
    ra_range = np.arange(0, 360.0*D2R, deg_between_comps*D2R)
    dec_range = np.arange(LOW_DEC, HIGH_DEC, deg_between_comps*D2R)
    
    
    
    ra_range, dec_range = np.meshgrid(ra_range, dec_range)
    ra_range, dec_range = ra_range.flatten(), dec_range.flatten()
    
    dec_range[np.abs(dec_range) < 1e-10] = 0.0
    # print('SKY MODEL MAKING', len(ra_range))
    
    source_index = 0
    
    point_index = 0
    gauss_index = 0
    shape_index = 0
    
    flux_index = 0
    basis_index = 0
    
    new_source = False
    
    point = False
    gaussian = False
    shapelet = False
    
    point = True
    shapelet = True
    gaussian = True
    
    with open("test_full_skymodel.yaml", 'w') as outfile:
        
        for ra, dec in zip(ra_range, dec_range):
            
            if point:
                source_index, point_index, flux_index = add_point_yaml(outfile, ra, dec,
                                                point_index, gauss_index, shape_index,
                                                CompTypes.POINT_POWER,
                                                flux_index, num_list_values,
                                                source_index, comps_per_source)
                
                source_index, point_index, flux_index = add_point_yaml(outfile, ra, dec,
                                                point_index, gauss_index, shape_index,
                                                CompTypes.POINT_CURVE,
                                                flux_index, num_list_values,
                                                source_index, comps_per_source)
                
                source_index, point_index, flux_index = add_point_yaml(outfile, ra, dec,
                                                point_index, gauss_index, shape_index,
                                                CompTypes.POINT_LIST,
                                                flux_index, num_list_values,
                                                source_index, comps_per_source)
                
            if gaussian:
                source_index, gauss_index, flux_index = add_gauss_yaml(outfile, ra, dec,
                                                gauss_index, gauss_index, shape_index,
                                                CompTypes.GAUSS_POWER,
                                                flux_index, num_list_values,
                                                source_index, comps_per_source)
                
                source_index, gauss_index, flux_index = add_gauss_yaml(outfile, ra, dec,
                                                gauss_index, gauss_index, shape_index,
                                                CompTypes.GAUSS_CURVE,
                                                flux_index, num_list_values,
                                                source_index, comps_per_source)
                
                source_index, gauss_index, flux_index = add_gauss_yaml(outfile, ra, dec,
                                                gauss_index, gauss_index, shape_index,
                                                CompTypes.GAUSS_LIST,
                                                flux_index, num_list_values,
                                                source_index, comps_per_source)
           
            if shapelet: 
                source_index, shape_index, flux_index, basis_index = add_shapelet_yaml(outfile, ra, dec,
                                                point_index, gauss_index, shape_index,
                                                CompTypes.SHAPE_POWER,
                                                flux_index, num_list_values,
                                                source_index, comps_per_source,
                                                num_coeff_per_shape, basis_index)
                
                source_index, shape_index, flux_index, basis_index = add_shapelet_yaml(outfile, ra, dec,
                                                point_index, gauss_index, shape_index,
                                                CompTypes.SHAPE_CURVE,
                                                flux_index, num_list_values,
                                                source_index, comps_per_source,
                                                num_coeff_per_shape, basis_index)
                
                source_index, shape_index, flux_index, basis_index = add_shapelet_yaml(outfile, ra, dec,
                                                point_index, gauss_index, shape_index,
                                                CompTypes.SHAPE_LIST,
                                                flux_index, num_list_values,
                                                source_index, comps_per_source,
                                                num_coeff_per_shape, basis_index)
                
    return ra_range, dec_range

##based class on an existing class test that includes methods for
##for 
class Test(BaseChunkTest):
    
    def run_test_read_yaml_skymodel_chunk(self, skymodel_filename, expected_chunks,
                                          max_num_visibilities, num_baselines,
                                          num_freqs, num_time_steps,
                                          lst,
                                          crop_by_component=True):
        
        lst = 0.0
    
        woden_settings = ws.Woden_Settings_Double()
        
        woden_settings.time_res = 1.0
        woden_settings.latitude = -0.46606083776035967
        woden_settings.latitude_obs_epoch_base = -0.46606083776035967
        
        woden_settings.lst_base = 0.0
        woden_settings.lst_obs_epoch_base = 0.0
        
        woden_settings.jd_date = 2457278.201145833
        woden_settings.num_time_steps = num_time_steps
        
        woden_settings.do_precession = 1
        lsts = ws.setup_lsts_and_phase_centre(woden_settings)
        
        # print(lsts)
        
        comp_counter = read_yaml_skymodel.read_yaml_radec_count_components(skymodel_filename)
        
        
        comp_counter = crop_below_horizon(lst, MWA_LAT,
                                          comp_counter, 
                                          crop_by_component=crop_by_component)
        
        ##Create a chunking map
        chunked_skymodel_maps = create_skymodel_chunk_map(comp_counter,
                                        max_num_visibilities, num_baselines,
                                        num_freqs, num_time_steps)
        
        beamtype = BeamTypes.FEE_BEAM.value

        source_catalogue = read_yaml_skymodel.read_yaml_skymodel_chunks(
                                              skymodel_filename, chunked_skymodel_maps,
                                              num_freqs, num_time_steps,
                                              beamtype, lsts, MWA_LAT)
        
        check_all_sources(expected_chunks, source_catalogue)
        
    def run_write_model_test_read_yaml_skymodel_chunk(self, deg_between_comps : int,
                       num_coeff_per_shape : int, num_list_values : int,
                       comps_per_source : int, lst : float,
                       max_num_visibilities : float):
        """Let's go do a bit"""
        
        num_freqs = 16
        num_baselines = 8128
        num_time_steps = 14
        
        ra_range, dec_range = write_full_test_skymodel_yaml(deg_between_comps,
                                 num_coeff_per_shape,
                                 num_list_values,
                                 comps_per_source)
        
        # filename = f"{code_dir}/test_full_skymodel.yaml"
        filename = "test_full_skymodel.yaml"
        
        comps_per_chunk = int(np.floor(max_num_visibilities / (num_baselines * num_freqs * num_time_steps)))
        
        
        expec_skymodel_chunks = make_expected_chunks(ra_range, dec_range,
                             num_coeff_per_shape, num_list_values,
                             comps_per_source, comps_per_chunk)
        
        
        self.run_test_read_yaml_skymodel_chunk(filename, expec_skymodel_chunks,
                                          max_num_visibilities, num_baselines,
                                          num_freqs, num_time_steps, lst)
    
    # def test_the_first(self):
    #     deg_between_comps = 240
    #     num_coeff_per_shape = 4
    #     num_list_values = 4
    #     comps_per_source = 10
    #     lst = 0.0
    #     max_num_visibilities = 1e7
        
    #     self.run_write_model_test_read_yaml_skymodel_chunk(deg_between_comps,
    #                                            num_coeff_per_shape,
    #                                            num_list_values,
    #                                            comps_per_source, lst,
    #                                            max_num_visibilities)
    #     print('-----------------------------')
        
    def test_the_second(self):
        deg_between_comps = 10
        num_coeff_per_shape = 4
        num_list_values = 4
        comps_per_source = 11
        lst = 0.0
        max_num_visibilities = 1e8
        
        self.run_write_model_test_read_yaml_skymodel_chunk(deg_between_comps,
                                               num_coeff_per_shape,
                                               num_list_values,
                                               comps_per_source, lst,
                                               max_num_visibilities)
        print('-----------------------------')
        
    # def test_the_third(self):
    #     deg_between_comps = 7
    #     num_coeff_per_shape = 6
    #     num_list_values = 4
    #     comps_per_source = 20
    #     lst = np.pi
    #     max_num_visibilities = 1e10
        
    #     self.run_write_model_test_read_yaml_skymodel_chunk(deg_between_comps,
    #                                            num_coeff_per_shape,
    #                                            num_list_values,
    #                                            comps_per_source, lst,
    #                                            max_num_visibilities)
    #     print('-----------------------------')



        
##Run the test
if __name__ == '__main__':
    unittest.main()
    