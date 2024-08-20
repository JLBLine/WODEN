"""

"""

from sys import path
import os
import unittest
import numpy as np
import erfa
from astropy.table import Column, Table
from astropy.io import fits

# ##Code we are testing
from wodenpy.skymodel import read_fits_skymodel
# import wodenpy
from wodenpy.skymodel.woden_skymodel import Component_Type_Counter, CompTypes, crop_below_horizon
from wodenpy.skymodel.chunk_sky_model import create_skymodel_chunk_map, Skymodel_Chunk_Map, increment_flux_type_counters
from wodenpy.use_libwoden.beam_settings import BeamTypes
from wodenpy.use_libwoden.skymodel_structs import setup_source_catalogue, setup_chunked_source
import wodenpy.use_libwoden.woden_settings as ws

from common_skymodel_test import fill_comp_counter_for_chunking, Expec_Counter, BaseChunkTest, Expected_Sky_Chunk, Expected_Components, Skymodel_Settings
from read_skymodel_common import check_components, check_all_sources, populate_pointgauss_chunk, populate_shapelet_chunk, make_expected_chunks
from fits_skymodel_common import write_full_test_skymodel_fits
from wodenpy.use_libwoden.create_woden_struct_classes import Woden_Struct_Classes
import wodenpy.use_libwoden.skymodel_structs


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

##based class on an existing class test that includes methods for
##for 
class Test(BaseChunkTest):
    
    def run_test_read_fits_skymodel_chunk(self, skymodel_filename, expected_chunks,
                                          max_num_visibilities, num_baselines,
                                          num_freqs, num_time_steps,
                                          lst,
                                          crop_by_component=True):
        
        precision = "double"
        woden_struct_classes = Woden_Struct_Classes(precision)
        woden_settings = woden_struct_classes.Woden_Settings()
        
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
        
        comp_counter = read_fits_skymodel.read_fits_radec_count_components(skymodel_filename)
        
        # comp_counter.print_info()
        
        
        comp_counter = crop_below_horizon(lst, MWA_LAT,
                                          comp_counter, 
                                          crop_by_component=crop_by_component)
        
        comp_counter.print_info()
        
        ##Create a chunking map
        chunked_skymodel_maps = create_skymodel_chunk_map(comp_counter,
                                        max_num_visibilities, num_baselines,
                                        num_freqs, num_time_steps)
        
        beamtype = BeamTypes.FEE_BEAM.value
        
        main_table = Table.read(skymodel_filename, hdu=1)
        shape_table = Table.read(skymodel_filename, hdu=2)
        
        with fits.open(skymodel_filename) as hdus:
            num_hdus = len(hdus)
            hdu_names = [hdu.name for hdu in hdus]
            
        if 'V_LIST_FLUXES' in hdu_names:
                v_table = Table.read(skymodel_filename, hdu='V_LIST_FLUXES')
        else:
            v_table = False
            
        if 'Q_LIST_FLUXES' in hdu_names:
                q_table = Table.read(skymodel_filename, hdu='Q_LIST_FLUXES')
        else:
            q_table = False    
            
        if 'U_LIST_FLUXES' in hdu_names:
                u_table = Table.read(skymodel_filename, hdu='U_LIST_FLUXES')
        else:
            u_table = False
            
        if 'P_LIST_FLUXES' in hdu_names:
                p_table = Table.read(skymodel_filename, hdu='P_LIST_FLUXES')
        else:
            p_table = False
            
        # print("TABLES", type(v_table), type(q_table), type(u_table), type(p_table))

        source_catalogue = read_fits_skymodel.read_fits_skymodel_chunks(woden_struct_classes,
                                              main_table, shape_table, chunked_skymodel_maps,
                                              num_freqs, num_time_steps,
                                              beamtype, lsts, MWA_LAT,
                                              v_table, q_table, u_table, p_table,)
        
        check_all_sources(expected_chunks, source_catalogue,
                           fits_skymodel=True)
        
    def run_write_model_test_read_fits_skymodel_chunk(self, 
                                                      skymodel_settings : Skymodel_Settings,
                                                      lst : float,
                       max_num_visibilities : float):
        """Let's go do a bit"""
        
        num_freqs = 16
        num_baselines = 8128
        num_time_steps = 14
        
        ra_range, dec_range = write_full_test_skymodel_fits(skymodel_settings)
        
        # filename = f"{code_dir}/test_full_skymodel.fits"
        filename = "test_full_skymodel.fits"
        
        comps_per_chunk = int(np.floor(max_num_visibilities / (num_baselines * num_freqs * num_time_steps)))
        
        
        expec_skymodel_chunks = make_expected_chunks(ra_range, dec_range,
                             skymodel_settings, comps_per_chunk, lst=lst,
                             fits_skymodel=True)
        
        
        self.run_test_read_fits_skymodel_chunk(filename, expec_skymodel_chunks,
                                          max_num_visibilities, num_baselines,
                                          num_freqs, num_time_steps, lst)
    
    def test_the_first(self):
        deg_between_comps = 240
        num_coeff_per_shape = 4
        num_list_values = 4
        comps_per_source = 10
        lst = 0.0
        max_num_visibilities = 1e7
        
        settings = Skymodel_Settings(deg_between_comps,
                                          num_coeff_per_shape,
                                          num_list_values,
                                          comps_per_source)
        
        self.run_write_model_test_read_fits_skymodel_chunk(settings, lst,
                                               max_num_visibilities)
        print('-----------------------------')
        
    def test_the_second(self):
        deg_between_comps = 10
        num_coeff_per_shape = 4
        num_list_values = 4
        comps_per_source = 11
        lst = 0.0
        max_num_visibilities = 1e8
        
        settings = Skymodel_Settings(deg_between_comps,
                                          num_coeff_per_shape,
                                          num_list_values,
                                          comps_per_source)
        
        self.run_write_model_test_read_fits_skymodel_chunk(settings, lst,
                                               max_num_visibilities)
        print('-----------------------------')
        
    def test_the_third(self):
        deg_between_comps = 7
        num_coeff_per_shape = 6
        num_list_values = 4
        comps_per_source = 20
        lst = 0.0
        max_num_visibilities = 1e10
        
        settings = Skymodel_Settings(deg_between_comps,
                                          num_coeff_per_shape,
                                          num_list_values,
                                          comps_per_source)
        
        self.run_write_model_test_read_fits_skymodel_chunk(settings, lst,
                                               max_num_visibilities)
        print('-----------------------------')
        
    def test_the_fourth(self):
        deg_between_comps = 60
        num_coeff_per_shape = 6
        num_list_values = 4
        comps_per_source = 20
        lst = 0.0
        max_num_visibilities = 1e10
        
        settings = Skymodel_Settings(deg_between_comps,
                                     num_coeff_per_shape,
                                     num_list_values,
                                     comps_per_source,
                                     stokesV_frac_cadence = 5)
        
        self.run_write_model_test_read_fits_skymodel_chunk(settings, lst,
                                               max_num_visibilities)
        print('-----------------------------')
        
    def test_the_fifth(self):
        deg_between_comps = 60
        num_coeff_per_shape = 6
        num_list_values = 4
        comps_per_source = 20
        lst = 0.0
        max_num_visibilities = 1e10
        
        settings = Skymodel_Settings(deg_between_comps,
                                     num_coeff_per_shape,
                                     num_list_values,
                                     comps_per_source,
                                     stokesV_pl_cadence = 4)
        
        self.run_write_model_test_read_fits_skymodel_chunk(settings, lst,
                                               max_num_visibilities)
        print('-----------------------------')
        
    def test_the_sixth(self):
        deg_between_comps = 60
        num_coeff_per_shape = 6
        num_list_values = 4
        comps_per_source = 20
        lst = 0.0
        max_num_visibilities = 1e10
        
        settings = Skymodel_Settings(deg_between_comps,
                                     num_coeff_per_shape,
                                     num_list_values,
                                     comps_per_source,
                                     stokesV_cpl_cadence = 3)
        
        self.run_write_model_test_read_fits_skymodel_chunk(settings, lst,
                                               max_num_visibilities)
        print('-----------------------------')
        
    def test_the_seventh(self):
        deg_between_comps = 60
        num_coeff_per_shape = 6
        num_list_values = 4
        comps_per_source = 20
        lst = 0.0
        max_num_visibilities = 1e10
        
        settings = Skymodel_Settings(deg_between_comps,
                                     num_coeff_per_shape,
                                     num_list_values,
                                     comps_per_source,
                                     linpol_frac_cadence = 3)
        
        self.run_write_model_test_read_fits_skymodel_chunk(settings, lst,
                                               max_num_visibilities)
        print('-----------------------------')
        
    def test_the_eighth(self):
        deg_between_comps = 60
        num_coeff_per_shape = 6
        num_list_values = 4
        comps_per_source = 20
        lst = 0.0
        max_num_visibilities = 1e10
        
        settings = Skymodel_Settings(deg_between_comps,
                                     num_coeff_per_shape,
                                     num_list_values,
                                     comps_per_source,
                                     linpol_pl_cadence = 5)
        
        self.run_write_model_test_read_fits_skymodel_chunk(settings, lst,
                                               max_num_visibilities)
        print('-----------------------------')
        
    def test_the_ninth(self):
        deg_between_comps = 60
        num_coeff_per_shape = 6
        num_list_values = 4
        comps_per_source = 20
        lst = 0.0
        max_num_visibilities = 1e10
        
        settings = Skymodel_Settings(deg_between_comps,
                                     num_coeff_per_shape,
                                     num_list_values,
                                     comps_per_source,
                                     linpol_cpl_cadence = 4)
        
        self.run_write_model_test_read_fits_skymodel_chunk(settings, lst,
                                               max_num_visibilities)
        print('-----------------------------')
        
    def test_the_tenth(self):
        deg_between_comps = 60
        num_coeff_per_shape = 6
        num_list_values = 4
        comps_per_source = 20
        lst = 0.0
        max_num_visibilities = 1e10
        
        settings = Skymodel_Settings(deg_between_comps,
                                     num_coeff_per_shape,
                                     num_list_values,
                                     comps_per_source,
                                     stokesV_list_cadence=5,
                                     stokesV_num_list=3)
        
        self.run_write_model_test_read_fits_skymodel_chunk(settings, lst,
                                               max_num_visibilities)
        print('-----------------------------')
        
    def test_the_eleventh(self):
        deg_between_comps = 60
        num_coeff_per_shape = 6
        num_list_values = 4
        comps_per_source = 20
        lst = 0.0
        max_num_visibilities = 1e10
        
        settings = Skymodel_Settings(deg_between_comps,
                                     num_coeff_per_shape,
                                     num_list_values,
                                     comps_per_source,
                                     linpol_list_cadence=6,
                                     linpol_num_list=4)
        
        self.run_write_model_test_read_fits_skymodel_chunk(settings, lst,
                                               max_num_visibilities)
        print('-----------------------------')
        
    def test_the_twelfth(self):
        deg_between_comps = 60
        num_coeff_per_shape = 6
        num_list_values = 4
        comps_per_source = 20
        lst = 0.0
        max_num_visibilities = 1e10
        
        settings = Skymodel_Settings(deg_between_comps,
                                     num_coeff_per_shape,
                                     num_list_values,
                                     comps_per_source,
                                     linpol_p_list_cadence=3,
                                     linpol_num_p_list=5)
        
        self.run_write_model_test_read_fits_skymodel_chunk(settings, lst,
                                               max_num_visibilities)
        print('-----------------------------')

        
    def test_the_big_one(self):
        deg_between_comps = 7
        num_coeff_per_shape = 6
        num_list_values = 4
        comps_per_source = 20
        lst = 0.0
        max_num_visibilities = 1e10
        
        settings = Skymodel_Settings(deg_between_comps,
                                     num_coeff_per_shape,
                                     num_list_values,
                                     comps_per_source,
                                     stokesV_frac_cadence = 5,
                                     stokesV_pl_cadence = 4,
                                     stokesV_cpl_cadence = 3,
                                     linpol_frac_cadence = 3,
                                     linpol_pl_cadence = 5,
                                     linpol_cpl_cadence = 4,
                                     stokesV_list_cadence=7,
                                     stokesV_num_list=3,
                                     linpol_list_cadence=7,
                                     linpol_num_list=4,
                                     linpol_p_list_cadence=11,
                                     linpol_num_p_list=5)
        
        self.run_write_model_test_read_fits_skymodel_chunk(settings, lst,
                                               max_num_visibilities)
        print('-----------------------------')

##Run the test
if __name__ == '__main__':
    unittest.main()
    