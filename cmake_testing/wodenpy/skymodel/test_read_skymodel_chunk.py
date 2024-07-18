"""

"""

from sys import path
import os
import unittest
import numpy as np
import erfa

# ##Code we are testing
from wodenpy.skymodel.read_skymodel import read_skymodel_chunks, read_radec_count_components
from wodenpy.skymodel import read_fits_skymodel
from wodenpy.skymodel import read_yaml_skymodel
from wodenpy.skymodel.woden_skymodel import Component_Type_Counter, CompTypes, crop_below_horizon
from wodenpy.skymodel.chunk_sky_model import create_skymodel_chunk_map, Skymodel_Chunk_Map, increment_flux_type_counters
from wodenpy.use_libwoden.beam_settings import BeamTypes

from wodenpy.use_libwoden.skymodel_structs import setup_chunked_source, _Ctype_Source_Into_Python

from common_skymodel_test import fill_comp_counter_for_chunking, Expec_Counter, BaseChunkTest, Expected_Sky_Chunk, Expected_Components

from read_skymodel_common import check_components, check_all_sources, populate_pointgauss_chunk, populate_shapelet_chunk, make_expected_chunks, Skymodel_Settings

import wodenpy.use_libwoden.woden_settings as ws

from test_read_FITS_skymodel_chunk import write_full_test_skymodel_fits
from test_read_yaml_skymodel_chunk import write_full_test_skymodel_yaml
from test_read_text_skymodel_chunk import write_full_test_skymodel_text, make_expected_chunks_text

D2R = np.pi/180.0
MWA_LAT = -26.703319405555554*D2R

##for now, WODEN has three flux types: power law, curved power law, and list
NUM_FLUX_TYPES = 3

RTOL=1e-10


##test both FITS and yaml versions with real simple models, as other
##tests of relevant functions do more
num_freqs = 16
num_baselines = 8128
num_time_steps = 14
deg_between_comps = 240
num_coeff_per_shape = 4
num_list_values = 4
comps_per_source = 10
lst = 0.0
max_num_visibilities = 1e7

##based class on an existing class test that includes methods for
##for 
class Test(BaseChunkTest):
    
    def test_read_skymodel_chunks_with_fits(self):
        
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
        
        settings = Skymodel_Settings(deg_between_comps, num_coeff_per_shape,
                                     num_list_values, comps_per_source)
        
        
        ##create a test FITS skymodel to read
        ra_range, dec_range = write_full_test_skymodel_fits(settings)
        
        skymodel_filename = "test_full_skymodel.fits"
        
        ##come up with expected values
        comps_per_chunk = int(np.floor(max_num_visibilities / (num_baselines * num_freqs * num_time_steps)))
        
        expec_skymodel_chunks = make_expected_chunks(ra_range, dec_range,
                             num_coeff_per_shape, num_list_values,
                             comps_per_source, comps_per_chunk,
                             fits_skymodel=True)
        comp_counter = read_radec_count_components(skymodel_filename)
        
        
        comp_counter = crop_below_horizon(lst, MWA_LAT,
                                          comp_counter, 
                                          crop_by_component=True)
        
        ##Create a chunking map
        chunked_skymodel_maps = create_skymodel_chunk_map(comp_counter,
                                        max_num_visibilities, num_baselines,
                                        num_freqs, num_time_steps)
        
        
        beamtype = BeamTypes.FEE_BEAM.value

        source_catalogue = read_skymodel_chunks(skymodel_filename,
                                                chunked_skymodel_maps,
                                                num_freqs, num_time_steps,
                                                beamtype, lsts, MWA_LAT)
        
        check_all_sources(expec_skymodel_chunks, source_catalogue,
                           fits_skymodel=True)
        
        
    def test_read_skymodel_chunks_with_yaml(self):
        
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
        
        
        ##create a test FITS skymodel to read
        ra_range, dec_range = write_full_test_skymodel_yaml(deg_between_comps,
                                 num_coeff_per_shape,
                                 num_list_values,
                                 comps_per_source)
        
        skymodel_filename = "test_full_skymodel.yaml"
        
        ##come up with expected values
        comps_per_chunk = int(np.floor(max_num_visibilities / (num_baselines * num_freqs * num_time_steps)))
        
        expec_skymodel_chunks = make_expected_chunks(ra_range, dec_range,
                             num_coeff_per_shape, num_list_values,
                             comps_per_source, comps_per_chunk)
        
        comp_counter = read_radec_count_components(skymodel_filename)
        
        
        comp_counter = crop_below_horizon(lst, MWA_LAT,
                                          comp_counter, 
                                          crop_by_component=True)
        
        ##Create a chunking map
        chunked_skymodel_maps = create_skymodel_chunk_map(comp_counter,
                                        max_num_visibilities, num_baselines,
                                        num_freqs, num_time_steps)
        
        
        beamtype = BeamTypes.FEE_BEAM.value

        source_catalogue = read_skymodel_chunks(skymodel_filename,
                                                chunked_skymodel_maps,
                                                num_freqs, num_time_steps,
                                                beamtype, lsts, MWA_LAT)
        
        check_all_sources(expec_skymodel_chunks, source_catalogue)

    def test_read_skymodel_chunks_with_text(self):
        
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
        
        
        ##create a test FITS skymodel to read
        ra_range, dec_range = write_full_test_skymodel_text(deg_between_comps,
                                 num_coeff_per_shape,
                                 comps_per_source)
        
        skymodel_filename = "test_full_skymodel.txt"
        
        ##come up with expected values
        comps_per_chunk = int(np.floor(max_num_visibilities / (num_baselines * num_freqs * num_time_steps)))
        
        expec_skymodel_chunks = make_expected_chunks_text(ra_range,
                                                          dec_range,
                             num_coeff_per_shape, comps_per_source,
                             comps_per_chunk)
        
        comp_counter = read_radec_count_components(skymodel_filename)
        
        
        comp_counter = crop_below_horizon(lst, MWA_LAT,
                                          comp_counter, 
                                          crop_by_component=True)
        
        ##Create a chunking map
        chunked_skymodel_maps = create_skymodel_chunk_map(comp_counter,
                                        max_num_visibilities, num_baselines,
                                        num_freqs, num_time_steps)
        
        
        beamtype = BeamTypes.FEE_BEAM.value

        source_catalogue = read_skymodel_chunks(skymodel_filename,
                                                chunked_skymodel_maps,
                                                num_freqs, num_time_steps,
                                                beamtype, lsts, MWA_LAT)
        
        check_all_sources(expec_skymodel_chunks, source_catalogue)

##Run the test
if __name__ == '__main__':
    unittest.main()
    