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

from common_skymodel_test import fill_comp_counter, Expec_Counter, BaseChunkTest, Expected_Sky_Chunk, Expected_Components
from read_skymodel_common import check_components, check_all_sources, populate_pointgauss_chunk, populate_shapelet_chunk, make_expected_chunks

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


def add_power_law_fits(table_dict, row_ind, value):
    """add power law info for a component"""
    
    table_dict['mod_type'].append('pl')
    table_dict['norm_comp_pl'][row_ind] = float(value)
    table_dict['alpha_pl'][row_ind] = float(value)
    
def add_curved_power_law_fits(table_dict, row_ind, value):
    """add curved power law infor for a component"""
    table_dict['mod_type'].append('cpl')
    table_dict['norm_comp_cpl'][row_ind] = float(value)
    table_dict['alpha_cpl'][row_ind] = float(value)
    table_dict['curve_cpl'][row_ind] = float(value)
    
def add_list_flux_fits(table_dict, all_flux_cols, row_ind, flux_index, num_list_values):
    """add list law infor for a component"""
    table_dict['mod_type'].append('nan')
    
    for flux_col in range(num_list_values):
        
        all_flux_cols[flux_col][row_ind] = float(flux_index)
        
        flux_index += 1
        
    return flux_index
    

def add_point_fits(table_dict : dict, all_flux_cols : list,
              ra : float, dec: float,
              point_index : int, gauss_index : int, shape_index : int,
              comp_type : CompTypes,
              flux_index : int, num_list_values : int,
              source_index : int, comps_per_source : int):
    
    comp_index = point_index + gauss_index + shape_index
    if comp_index % comps_per_source == 0:
        source_index += 1
    
    table_dict['unq_source_id'].append(f"{source_index:07d}")
    table_dict['name'].append(f"{source_index:07d}_C{comp_index:05d}")
    table_dict['comp_type'].append('P')
    
    # print("DIS", source_index, comps_per_source, comp_index)
    
    # row_ind = source_index*comps_per_source + comp_index % comps_per_source
    row_ind = comp_index
    table_dict['ra'][row_ind] = ra
    table_dict['dec'][row_ind] = dec
    
    if comp_type == CompTypes.POINT_POWER:
        add_power_law_fits(table_dict, row_ind, point_index)
        
    elif comp_type == CompTypes.POINT_CURVE:
        add_curved_power_law_fits(table_dict, row_ind, point_index)
        
    elif comp_type == CompTypes.POINT_LIST:
        flux_index = add_list_flux_fits(table_dict, all_flux_cols, row_ind, flux_index, num_list_values)
        
    point_index += 1
        
    return source_index, point_index, flux_index

def add_gauss_fits(table_dict : dict, all_flux_cols : list,
                    ra : float, dec: float,
                    point_index : int, gauss_index : int, shape_index : int,
                    comp_type : CompTypes,
                    flux_index : int, num_list_values : int,
                    source_index : int, comps_per_source : int):
    
    comp_index = point_index + gauss_index + shape_index
    if comp_index % comps_per_source == 0:
        source_index += 1
    
    table_dict['unq_source_id'].append(f"{source_index:07d}")
    table_dict['name'].append(f"{source_index:07d}_C{comp_index:05d}")
    table_dict['comp_type'].append('G')
    
    row_ind = comp_index
    table_dict['ra'][row_ind] = ra
    table_dict['dec'][row_ind] = dec
    
    table_dict['major_dc'][row_ind] = float(gauss_index) / 3600.0
    table_dict['minor_dc'][row_ind] = float(gauss_index) / 3600.0
    table_dict['pa_dc'][row_ind] = float(gauss_index)
    
    if comp_type == CompTypes.GAUSS_POWER:
        add_power_law_fits(table_dict, row_ind, gauss_index)
        
    elif comp_type == CompTypes.GAUSS_CURVE:
        add_curved_power_law_fits(table_dict, row_ind, gauss_index)
        
    elif comp_type == CompTypes.GAUSS_LIST:
        flux_index = add_list_flux_fits(table_dict, all_flux_cols, row_ind, flux_index, num_list_values)

    gauss_index += 1
        
    return source_index, gauss_index, flux_index
    

def add_shapelet_fits(table_dict : dict, all_flux_cols : list,
                    ra : float, dec: float,
                    point_index : int, gauss_index : int, shape_index : int,
                    comp_type : CompTypes,
                    flux_index : int, num_list_values : int,
                    source_index : int, comps_per_source : int,
                    num_coeff_per_shape : int, basis_index : int):
    
    comp_index = point_index + gauss_index + shape_index
    if comp_index % comps_per_source == 0:
        source_index += 1
    
    table_dict['unq_source_id'].append(f"{source_index:07d}")
    table_dict['name'].append(f"{source_index:07d}_C{comp_index:05d}")
    table_dict['comp_type'].append('S')
        
    row_ind = comp_index
    table_dict['ra'][row_ind] = ra
    table_dict['dec'][row_ind] = dec
    table_dict['major_dc'][row_ind] = float(shape_index) / 3600.0
    table_dict['minor_dc'][row_ind] = float(shape_index) / 3600.0
    table_dict['pa_dc'][row_ind] = float(shape_index)
    
    for basis in range(basis_index, basis_index + num_coeff_per_shape):
    
        table_dict['s_names'].append(f"{source_index:07d}_C{comp_index:05d}")
        table_dict['s_n1s'][basis] = float(basis)
        table_dict['s_n2s'][basis] = float(basis)
        table_dict['s_coeffs'][basis] = float(basis)
    
    basis_index += num_coeff_per_shape
    
    if comp_type == CompTypes.SHAPE_POWER:
        add_power_law_fits(table_dict, row_ind, shape_index)
        
    elif comp_type == CompTypes.SHAPE_CURVE:
        add_curved_power_law_fits(table_dict, row_ind, shape_index)
        
    elif comp_type == CompTypes.SHAPE_LIST:
        flux_index = add_list_flux_fits(table_dict, all_flux_cols, row_ind, flux_index, num_list_values)
        
    shape_index += 1
        
    return source_index, shape_index, flux_index, basis_index


def write_full_test_skymodel_fits(deg_between_comps : float,
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
    print('There are', len(ra_range), 'coord locations')
    
    source_index = -1
    
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
    gaussian = True
    shapelet = True
    
    
    # with open("test_full_skymodel.fits", 'w') as outfile:
    
    num_comp_types = point + gaussian + shapelet
    num_list_freqs = num_comp_types * len(ra_range) * num_list_values
    
    total_num_comps = NUM_FLUX_TYPES*num_comp_types*len(ra_range)
    
    
    ##Start by making empty fits Columns to populate
    # unq_source_id = Column(data=np.full(total_num_comps, np.nan, dtype='str'), name="UNQ_SOURCE_ID")
    # name = Column(data=np.full(total_num_comps, np.nan, dtype='str'), name="NAME")
    
    ##Anything that is a string column is easier to make as a list and the
    ##convert to column
    unq_source_id = []
    name = []
    
    ra_col = Column(data=np.full(total_num_comps, np.nan) ,name="RA")
    dec_col = Column(data=np.full(total_num_comps, np.nan), name="DEC")
    major_dc = Column(data=np.full(total_num_comps, np.nan), name="MAJOR_DC")
    minor_dc = Column(data=np.full(total_num_comps, np.nan), name="MINOR_DC")
    pa_dc = Column(data=np.full(total_num_comps, np.nan), name="PA_DC")
    # mod_type = Column(data=np.full(total_num_comps, np.nan, dtype='str'), name="MOD_TYPE")
    # comp_type = Column(data=np.full(total_num_comps, np.nan, dtype='str'), name="COMP_TYPE")
    mod_type = []
    comp_type = []
    norm_comp_pl = Column(data=np.full(total_num_comps, np.nan), name="NORM_COMP_PL")
    alpha_pl = Column(data=np.full(total_num_comps, np.nan), name="ALPHA_PL")
    norm_comp_cpl = Column(data=np.full(total_num_comps, np.nan), name="NORM_COMP_CPL")
    alpha_cpl = Column(data=np.full(total_num_comps, np.nan), name="ALPHA_CPL")
    curve_cpl = Column(data=np.full(total_num_comps, np.nan), name="CURVE_CPL")
    
    ##What these columns into a table so we can easily feed them into functions below
    table_dict = {}
    table_dict['unq_source_id'] = unq_source_id
    table_dict['name'] = name
    table_dict['ra'] = ra_col
    table_dict['dec'] = dec_col
    table_dict['major_dc'] = major_dc
    table_dict['minor_dc'] = minor_dc
    table_dict['pa_dc'] = pa_dc
    table_dict['mod_type'] = mod_type
    table_dict['comp_type'] = comp_type
    table_dict['norm_comp_pl'] = norm_comp_pl
    table_dict['alpha_pl'] = alpha_pl
    table_dict['norm_comp_cpl'] = norm_comp_cpl
    table_dict['alpha_cpl'] = alpha_cpl
    table_dict['curve_cpl'] = curve_cpl
    
    ##Let's populate the source and component names
    
    ##OK in the FITS format, every different frequency has a different column
    ##For the test I setup, where every list entry has a different flux, this
    ##is a little mental, but easier to use the testing logic made for the yaml
    ##So just cop it
    all_flux_cols = [Column(data=np.full(total_num_comps, np.nan), name=f"INT_FLX{flux_col_ind:07.3f}") for flux_col_ind in range(num_list_values)]
    
    num_shape_coeffs = shapelet*num_coeff_per_shape*NUM_FLUX_TYPES*len(ra_range)
    
    # s_names = Column(data=np.full(num_shape_coeffs, np.nan, dtype='str'), name="NAME")
    s_names = []
    s_n1s = Column(data=np.full(num_shape_coeffs, np.nan, dtype=int), name="N1")
    s_n2s = Column(data=np.full(num_shape_coeffs, np.nan, dtype=int), name="N2")
    s_coeffs = Column(data=np.full(num_shape_coeffs, np.nan), name="COEFF")
    
    table_dict['s_names'] = s_names
    table_dict['s_n1s'] = s_n1s
    table_dict['s_n2s'] = s_n2s
    table_dict['s_coeffs'] = s_coeffs
    
    for ra, dec in zip(ra_range/D2R, dec_range/D2R):
        
        if point:
            source_index, point_index, flux_index = add_point_fits(table_dict, all_flux_cols,
                                            ra, dec,
                                            point_index, gauss_index, shape_index,
                                            CompTypes.POINT_POWER,
                                            flux_index, num_list_values,
                                            source_index, comps_per_source)
            
            source_index, point_index, flux_index = add_point_fits(table_dict, all_flux_cols,
                                            ra, dec,
                                            point_index, gauss_index, shape_index,
                                            CompTypes.POINT_CURVE,
                                            flux_index, num_list_values,
                                            source_index, comps_per_source)
            
            source_index, point_index, flux_index = add_point_fits(table_dict, all_flux_cols,
                                            ra, dec,
                                            point_index, gauss_index, shape_index,
                                            CompTypes.POINT_LIST,
                                            flux_index, num_list_values,
                                            source_index, comps_per_source)
            
        if gaussian:
            source_index, gauss_index, flux_index = add_gauss_fits(table_dict, all_flux_cols,
                                            ra, dec,
                                            point_index, gauss_index, shape_index,
                                            CompTypes.GAUSS_POWER,
                                            flux_index, num_list_values,
                                            source_index, comps_per_source)
            
            source_index, gauss_index, flux_index = add_gauss_fits(table_dict, all_flux_cols,
                                            ra, dec,
                                            point_index, gauss_index, shape_index,
                                            CompTypes.GAUSS_CURVE,
                                            flux_index, num_list_values,
                                            source_index, comps_per_source)
            
            source_index, gauss_index, flux_index = add_gauss_fits(table_dict, all_flux_cols,
                                            ra, dec,
                                            point_index, gauss_index, shape_index,
                                            CompTypes.GAUSS_LIST,
                                            flux_index, num_list_values,
                                            source_index, comps_per_source)
        
        if shapelet: 
            source_index, shape_index, flux_index, basis_index = add_shapelet_fits(table_dict, all_flux_cols,
                                            ra, dec,
                                            point_index, gauss_index, shape_index,
                                            CompTypes.SHAPE_POWER,
                                            flux_index, num_list_values,
                                            source_index, comps_per_source,
                                            num_coeff_per_shape, basis_index)
            
            source_index, shape_index, flux_index, basis_index = add_shapelet_fits(table_dict, all_flux_cols,
                                            ra, dec,
                                            point_index, gauss_index, shape_index,
                                            CompTypes.SHAPE_CURVE,
                                            flux_index, num_list_values,
                                            source_index, comps_per_source,
                                            num_coeff_per_shape, basis_index)
            
            source_index, shape_index, flux_index, basis_index = add_shapelet_fits(table_dict, all_flux_cols,
                                            ra, dec,
                                            point_index, gauss_index, shape_index,
                                            CompTypes.SHAPE_LIST,
                                            flux_index, num_list_values,
                                            source_index, comps_per_source,
                                            num_coeff_per_shape, basis_index)
        
    unq_source_id = Column(data=unq_source_id, name="UNQ_SOURCE_ID")
    name = Column(data=name, name="NAME")
    mod_type = Column(data=mod_type, name="MOD_TYPE")
    comp_type = Column(data=comp_type, name="COMP_TYPE")
    s_names = Column(data=s_names, name="NAME")
    
    out_columns = [unq_source_id, name, ra_col, dec_col, major_dc, minor_dc, pa_dc, mod_type, comp_type,
                   norm_comp_pl, alpha_pl, norm_comp_cpl, alpha_cpl, curve_cpl]

    for flux_col in all_flux_cols:
        out_columns.append(flux_col)

    main_table = Table()
    main_table.add_columns(out_columns)
    
    shape_table = Table()
    shape_table.add_columns([s_names, s_n1s, s_n2s, s_coeffs])

    hdu_list = fits.HDUList([
        fits.PrimaryHDU(),
        fits.table_to_hdu(main_table),
        fits.table_to_hdu(shape_table),
    ])
    
    hdu_list.writeto("test_full_skymodel.fits", overwrite=True)
                
    return ra_range, dec_range

##based class on an existing class test that includes methods for
##for 
class Test(BaseChunkTest):
    
    def run_test_read_fits_skymodel_chunk(self, skymodel_filename, expected_chunks,
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
        
        comp_counter = read_fits_skymodel.read_fits_radec_count_components(skymodel_filename)
        
        
        comp_counter = crop_below_horizon(lst, MWA_LAT,
                                          comp_counter, 
                                          crop_by_component=crop_by_component)
        
        ##Create a chunking map
        chunked_skymodel_maps = create_skymodel_chunk_map(comp_counter,
                                        max_num_visibilities, num_baselines,
                                        num_freqs, num_time_steps)
        
        beamtype = BeamTypes.FEE_BEAM.value

        source_catalogue = read_fits_skymodel.read_fits_skymodel_chunks(
                                              skymodel_filename, chunked_skymodel_maps,
                                              num_freqs, num_time_steps,
                                              beamtype, lsts, MWA_LAT)
        
        check_all_sources(expected_chunks, source_catalogue,
                           fits_skymodel=True)
        
    def run_write_model_test_read_fits_skymodel_chunk(self, deg_between_comps : int,
                       num_coeff_per_shape : int, num_list_values : int,
                       comps_per_source : int, lst : float,
                       max_num_visibilities : float):
        """Let's go do a bit"""
        
        num_freqs = 16
        num_baselines = 8128
        num_time_steps = 14
        
        ra_range, dec_range = write_full_test_skymodel_fits(deg_between_comps,
                                 num_coeff_per_shape,
                                 num_list_values,
                                 comps_per_source)
        
        # filename = f"{code_dir}/test_full_skymodel.fits"
        filename = "test_full_skymodel.fits"
        
        comps_per_chunk = int(np.floor(max_num_visibilities / (num_baselines * num_freqs * num_time_steps)))
        
        
        expec_skymodel_chunks = make_expected_chunks(ra_range, dec_range,
                             num_coeff_per_shape, num_list_values,
                             comps_per_source, comps_per_chunk,
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
        
        self.run_write_model_test_read_fits_skymodel_chunk(deg_between_comps,
                                               num_coeff_per_shape,
                                               num_list_values,
                                               comps_per_source, lst,
                                               max_num_visibilities)
        print('-----------------------------')
        
    def test_the_second(self):
        deg_between_comps = 10
        num_coeff_per_shape = 4
        num_list_values = 4
        comps_per_source = 11
        lst = 0.0
        max_num_visibilities = 1e8
        
        self.run_write_model_test_read_fits_skymodel_chunk(deg_between_comps,
                                               num_coeff_per_shape,
                                               num_list_values,
                                               comps_per_source, lst,
                                               max_num_visibilities)
        print('-----------------------------')
        
    def test_the_third(self):
        deg_between_comps = 7
        num_coeff_per_shape = 6
        num_list_values = 4
        comps_per_source = 20
        lst = np.pi
        max_num_visibilities = 1e10
        
        self.run_write_model_test_read_fits_skymodel_chunk(deg_between_comps,
                                               num_coeff_per_shape,
                                               num_list_values,
                                               comps_per_source, lst,
                                               max_num_visibilities)
        print('-----------------------------')



        
##Run the test
if __name__ == '__main__':
    unittest.main()
    