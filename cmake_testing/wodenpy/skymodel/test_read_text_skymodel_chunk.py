"""

"""

from sys import path
import os
import unittest
import numpy as np
import erfa

# ##Code we are testing
from wodenpy.skymodel import read_text_skymodel
from wodenpy.skymodel import read_skymodel
from wodenpy.skymodel.woden_skymodel import Component_Type_Counter, CompTypes, crop_below_horizon
from wodenpy.skymodel.chunk_sky_model import create_skymodel_chunk_map, Skymodel_Chunk_Map, increment_flux_type_counters
from wodenpy.use_libwoden.beam_settings import BeamTypes
from wodenpy.skymodel.read_skymodel import create_source_catalogue_from_python_sources
from wodenpy.use_libwoden.skymodel_structs import setup_chunked_source, _Ctype_Source_Into_Python
from wodenpy.skymodel.read_fits_skymodel import read_fits_skymodel_chunks
from common_skymodel_test import fill_comp_counter_for_chunking, Expec_Counter, BaseChunkTest, Expected_Sky_Chunk, Expected_Components, Skymodel_Settings, Args

from wodenpy.use_libwoden.create_woden_struct_classes import Woden_Struct_Classes
import wodenpy.use_libwoden.woden_settings as ws

from read_skymodel_common import check_components, check_all_sources, populate_pointgauss_chunk, populate_shapelet_chunk, make_expected_chunks
import binpacking

D2R = np.pi/180.0
# MWA_LATITUDE = -26.7*D2R

MWA_LAT = -26.703319405555554*D2R

D2R = np.pi/180.0

##the WODEN skymodel format is ooooold, so can only handle power law models
NUM_FLUX_TYPES = 1

##we have point, gaussian, and shapelet types
NUM_COMP_TYPES = 3

##limits for "all sky" sky model
LOW_DEC = -90.0*D2R
HIGH_DEC = 30.0*D2R

RTOL=1e-10


def make_expected_chunks_text(ra_range, dec_range,
                         num_coeff_per_shape,
                         comps_per_source, comps_per_chunk, lst = 0.0,
                         fits_skymodel = False):
    
    num_list_values = 0
    num_of_each_comp = len(ra_range)
    
    settings = Skymodel_Settings(120,
                                num_coeff_per_shape,
                                num_list_values,
                                comps_per_source)
    
    # comp_index_range = np.arange(num_of_each_comp)
    coeff_range = np.arange(NUM_FLUX_TYPES*num_of_each_comp*num_coeff_per_shape)
    flux_list_range = np.arange(NUM_FLUX_TYPES*num_of_each_comp*num_list_values)
    
    has = lst - ra_range
    azs, els = erfa.hd2ae(has, dec_range, MWA_LAT)
    
    ##Crop below the horizon based on COMPONENT/SOURCE
    # include_flags = np.zeros(len(all_comp_types))
    above_horizon = np.where(els >= 0)[0]
    
    expec_ra = ra_range[above_horizon]
    expec_dec = dec_range[above_horizon]
    
    print(f"There are {len(expec_ra)} coords above horizon")
    
    expec_pow_fluxes = np.arange(0, num_of_each_comp*NUM_FLUX_TYPES, NUM_FLUX_TYPES)[above_horizon]
    expec_cur_fluxes = False
    num_crop_comp = len(above_horizon)
    
    num_point_chunks = int(np.ceil((NUM_FLUX_TYPES*num_crop_comp) / comps_per_chunk))
    num_gauss_chunks = int(np.ceil((NUM_FLUX_TYPES*num_crop_comp) / comps_per_chunk))
    num_coeff_chunks = int(np.ceil((NUM_FLUX_TYPES*num_coeff_per_shape*num_crop_comp) / comps_per_chunk))
    
    expec_skymodel_chunks = np.empty(num_point_chunks + num_gauss_chunks + num_coeff_chunks, dtype=Expected_Sky_Chunk)
    
    n_powers = num_crop_comp
    n_curves = 0
    n_lists = 0
    
    polvalues = False
    
    for chunk_ind in range(num_point_chunks):
        expec_chunk = populate_pointgauss_chunk(CompTypes.POINT, chunk_ind,
                            comps_per_chunk, n_powers,
                            n_curves, n_lists, settings,
                            num_of_each_comp, above_horizon,
                            expec_ra, expec_dec,
                            expec_pow_fluxes, expec_cur_fluxes,
                            polvalues,
                            fits_skymodel=fits_skymodel)
        
        expec_skymodel_chunks[chunk_ind] = expec_chunk
        
    for chunk_ind in range(num_gauss_chunks):
        expec_chunk = populate_pointgauss_chunk(CompTypes.GAUSSIAN, chunk_ind,
                            comps_per_chunk, n_powers,
                            n_curves, n_lists, settings,
                            num_of_each_comp, above_horizon,
                            expec_ra, expec_dec,
                            expec_pow_fluxes, expec_cur_fluxes,
                            polvalues,
                            fits_skymodel=fits_skymodel)
        
        expec_skymodel_chunks[num_point_chunks + chunk_ind] = expec_chunk
        
    total_shape_basis = num_coeff_per_shape*num_crop_comp*NUM_FLUX_TYPES
    
    
    ##OK OK OK so when chunking, the code resets the order to be
    ##all power first, then all curve second, then all list third.
    ##So even though it's written in the skymodel as power, curve, list,
    ##power, curve, list etc write a different expected order here
    
    shape_comp_ind_to_comp_type = np.empty(num_crop_comp*NUM_FLUX_TYPES, dtype=CompTypes)
    
    shape_comp_ind_to_comp_type[:num_crop_comp] = CompTypes.SHAPE_POWER
    # shape_comp_ind_to_comp_type[num_crop_comp:2*num_crop_comp] = CompTypes.SHAPE_CURVE
    # shape_comp_ind_to_comp_type[-num_crop_comp:] = CompTypes.SHAPE_LIST
    
    shape_basis_to_comp_ind = np.repeat(np.arange(num_crop_comp*NUM_FLUX_TYPES), num_coeff_per_shape)
    
    shape_expec_pow_fluxes = np.repeat(expec_pow_fluxes, num_coeff_per_shape)
    shape_expec_cur_fluxes = np.repeat(expec_cur_fluxes, num_coeff_per_shape)

    orig_comp_inds = np.empty(num_crop_comp*NUM_FLUX_TYPES)

    orig_comp_inds[:num_crop_comp] = above_horizon*NUM_FLUX_TYPES
    # orig_comp_inds[num_crop_comp:2*num_crop_comp] = above_horizon*NUM_FLUX_TYPES + 1
    # orig_comp_inds[-num_crop_comp:] = above_horizon*NUM_FLUX_TYPES + 2
    
    ##Ok, so we originally iterate through power, curved, list types
    ##of shapelets as we change the ra,dec. During this, we iterate the
    ##basis function values. The chunking however changes this to be all
    ##power first, then curved, then list. So write a reordered version here
    ##so we can make predictions easily for each chunk
    shape_basis_values = np.empty(total_shape_basis)
    
    for comp_ind, orig_ind in enumerate(above_horizon):
        
        ##this slots in the power law components
        low_coeff = orig_ind*num_coeff_per_shape*NUM_FLUX_TYPES
        low_ind = comp_ind*num_coeff_per_shape
        high_ind = low_ind + num_coeff_per_shape
        shape_basis_values[low_ind:high_ind] = range(low_coeff, low_coeff + num_coeff_per_shape)
        
    ##Need to account for the thread-based chunking reordering the chunk
    ##to better distribute the work
    
    
    num_shapes_per_comp = []
    
    for chunk_ind in range(num_coeff_chunks):
        
        coeff_lower = chunk_ind*comps_per_chunk
        coeff_higher = (chunk_ind + 1)*comps_per_chunk
        
        ##Are there are enough coeffs to fill the chunk?
        if (total_shape_basis >= coeff_higher):
            n_shape_coeffs = comps_per_chunk;
            
        else:
            n_shape_coeffs = total_shape_basis % comps_per_chunk
            
        basis_inds = np.arange(coeff_lower, coeff_lower+n_shape_coeffs)
        comp_inds = np.array(np.unique(shape_basis_to_comp_ind[basis_inds]), dtype=int)

        power_inds = np.unique(np.where(shape_comp_ind_to_comp_type[comp_inds] == CompTypes.SHAPE_POWER)[0])
        # curve_inds = np.unique(np.where(shape_comp_ind_to_comp_type[comp_inds] == CompTypes.SHAPE_CURVE)[0])
        # list_inds = np.unique(np.where(shape_comp_ind_to_comp_type[comp_inds] == CompTypes.SHAPE_LIST)[0])
        
        num_chunk_power = len(power_inds)
        num_chunk_curve = 0
        num_chunk_list = 0
        num_shapes_per_comp .append(num_chunk_power)
        
    ##We will have some unedfined number of chunks, so we want to split
    ##things as evenly as possible in the available number of threads
    indexed_shape_chunk_sizes = [(i, n_shape*num_coeff_per_shape) for i,n_shape in enumerate(num_shapes_per_comp)]  # List of (index, value) tuples
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
        
    for chunk_ind in range(num_coeff_chunks):
        
        coeff_lower = chunk_ind*comps_per_chunk
        coeff_higher = (chunk_ind + 1)*comps_per_chunk
        
        ##Are there are enough coeffs to fill the chunk?
        if (total_shape_basis >= coeff_higher):
            n_shape_coeffs = comps_per_chunk;
            
        else:
            n_shape_coeffs = total_shape_basis % comps_per_chunk
            
        basis_inds = np.arange(coeff_lower, coeff_lower+n_shape_coeffs)
        comp_inds = np.array(np.unique(shape_basis_to_comp_ind[basis_inds]), dtype=int)

        power_inds = np.unique(np.where(shape_comp_ind_to_comp_type[comp_inds] == CompTypes.SHAPE_POWER)[0])
        # curve_inds = np.unique(np.where(shape_comp_ind_to_comp_type[comp_inds] == CompTypes.SHAPE_CURVE)[0])
        # list_inds = np.unique(np.where(shape_comp_ind_to_comp_type[comp_inds] == CompTypes.SHAPE_LIST)[0])
        
        num_chunk_power = len(power_inds)
        num_chunk_curve = 0
        num_chunk_list = 0
        
        expec_chunk = Expected_Sky_Chunk()
        
        expec_chunk.init_shape_components(num_chunk_power, num_chunk_curve,
                                            num_chunk_list,
                                            num_list_values, comps_per_chunk)
        
        
        chunk_basis_values = shape_basis_values[basis_inds]
        
        chunk_basis_param_indexes = shape_basis_to_comp_ind[basis_inds]
        chunk_basis_param_indexes -= chunk_basis_param_indexes[0]
        
        expec_chunk = populate_shapelet_chunk(expec_chunk,
                            num_list_values, num_coeff_per_shape,
                            above_horizon, expec_ra, expec_dec,
                            comp_inds, orig_comp_inds[comp_inds],
                            chunk_basis_param_indexes,
                            chunk_basis_values,
                            settings,
                            fits_skymodel=fits_skymodel)
        
        new_ind = np.argsort(shape_comp_chunk_order)[chunk_ind]
        
        expec_skymodel_chunks[num_point_chunks + num_gauss_chunks + new_ind] = expec_chunk
        
    ##not doing polarisation in the FITS skymodel format, so go away
    
    return expec_skymodel_chunks


def add_power_law_text(outfile, value):
    """add power law infor for a component"""
    
    freq = 100e+6 + (value + 1)*1e+4
    outfile.write(f"LINEAR {freq:.1f} {float(value)} 0.0 0.0 0.0 {float(value)/100.0}\n")
    
def add_point_text(outfile, ra : float, dec: float,
                   point_index : int):
    
    outfile.write(f"COMPONENT POINT {ra/(D2R*15.0):.12f} {dec/D2R:.12f}\n")
        
    add_power_law_text(outfile, point_index)
    
    outfile.write(f"ENDCOMPONENT\n")
    
    point_index += 1
        
    return point_index

def add_gauss_text(outfile, ra : float, dec: float,
                    gauss_index : int):
    
    outfile.write(f"COMPONENT GAUSSIAN {ra/(D2R*15.0):.12f} {dec/D2R:.12f}\n")
    outfile.write(f"GPARAMS {gauss_index:.5f} {gauss_index/60.0:.12f} {gauss_index/60.0:.12f}\n")
    
    add_power_law_text(outfile, gauss_index)
    
    outfile.write(f"ENDCOMPONENT\n")
    
    gauss_index += 1
        
    return gauss_index
    

def add_shapelet_text(outfile, ra : float, dec: float,
                    shape_index : int,
                    num_coeff_per_shape : int, basis_index : int):
    
    outfile.write(f"COMPONENT SHAPELET {ra/(D2R*15.0):.12f} {dec/D2R:.12f}\n")
    outfile.write(f"SPARAMS {shape_index:.5f} {shape_index/60.0:.12f} {shape_index/60.0:.12f}\n")
    
    add_power_law_text(outfile, shape_index)
    
    for basis in range(basis_index, basis_index + num_coeff_per_shape):
    
        outfile.write(f"SCOEFF {basis} {basis} {basis}\n")
        
    basis_index += num_coeff_per_shape
    
    outfile.write(f"ENDCOMPONENT\n")
    
    shape_index += 1
    
    return shape_index, basis_index


def check_and_write_source(outfile, point_index, gauss_index, shape_index,
                               comps_per_source, source_index, num_coeff_per_shape,
                               num_radec):
        
    if (point_index + gauss_index + shape_index) % comps_per_source == 0:
    
        p_in_source = 0
        g_in_source = 0
        s_in_source = 0
        sum_in_new_source = 0
        total_comps = num_radec*NUM_COMP_TYPES
        
        sum_until_now = point_index + gauss_index + shape_index
        
        overall_sum = sum_in_new_source + sum_until_now
        
        if point_index == gauss_index == shape_index:
            
            for comp in range(sum_until_now, total_comps):
                
                p_in_source += 1
                sum_in_new_source = p_in_source + g_in_source + s_in_source
                overall_sum = sum_in_new_source + sum_until_now
                if sum_in_new_source == comps_per_source or overall_sum == total_comps:
                    break
                
                g_in_source += 1
                sum_in_new_source = p_in_source + g_in_source + s_in_source
                overall_sum = sum_in_new_source + sum_until_now
                if sum_in_new_source == comps_per_source or overall_sum == total_comps:
                    break
                
                s_in_source += 1
                sum_in_new_source = p_in_source + g_in_source + s_in_source
                overall_sum = sum_in_new_source + sum_until_now
                if sum_in_new_source == comps_per_source or overall_sum == total_comps:
                    break
                
                # print(p_in_source, g_in_source, s_in_source, sum_in_new_source, overall_sum)
        
        elif point_index > gauss_index and point_index > shape_index:
            
            for comp in range(sum_until_now, total_comps):
                
                g_in_source += 1
                sum_in_new_source = p_in_source + g_in_source + s_in_source
                overall_sum = sum_in_new_source + sum_until_now
                if sum_in_new_source == comps_per_source or overall_sum == total_comps:
                    break
                
                s_in_source += 1
                sum_in_new_source = p_in_source + g_in_source + s_in_source
                overall_sum = sum_in_new_source + sum_until_now
                if sum_in_new_source == comps_per_source or overall_sum == total_comps:
                    break
                
                p_in_source += 1
                sum_in_new_source = p_in_source + g_in_source + s_in_source
                overall_sum = sum_in_new_source + sum_until_now
                if sum_in_new_source == comps_per_source or overall_sum == total_comps:
                    break
                
        else:
            
            for comp in range(sum_until_now, total_comps):
                
                s_in_source += 1
                sum_in_new_source = p_in_source + g_in_source + s_in_source
                overall_sum = sum_in_new_source + sum_until_now
                if sum_in_new_source == comps_per_source or overall_sum == total_comps:
                    break
                
                p_in_source += 1
                sum_in_new_source = p_in_source + g_in_source + s_in_source
                overall_sum = sum_in_new_source + sum_until_now
                if sum_in_new_source == comps_per_source or overall_sum == total_comps:
                    break
                
                g_in_source += 1
                sum_in_new_source = p_in_source + g_in_source + s_in_source
                overall_sum = sum_in_new_source + sum_until_now
                if sum_in_new_source == comps_per_source or overall_sum == total_comps:
                    break
        
        
        if not point_index == gauss_index == shape_index == 0:
            outfile.write('ENDSOURCE\n')
        
        outfile.write(f"SOURCE {source_index:05d} P {p_in_source} G {g_in_source} S {s_in_source} {s_in_source*num_coeff_per_shape}\n")
    
        source_index += 1
                
                    
    return source_index


def write_full_test_skymodel_text(deg_between_comps : float,
                             num_coeff_per_shape : int,
                             comps_per_source : int):
    """Write a sky model covering the whole sky"""
    
    
    
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
    
    basis_index = 0
    
    # new_source = False
    
    point = False
    gaussian = False
    shapelet = False
    
    point = True
    gaussian = True
    shapelet = True
    
    with open("test_full_skymodel.txt", 'w') as outfile:
        
        num_radec = len(ra_range)
        for ra, dec in zip(ra_range, dec_range):
            
            source_index = check_and_write_source(outfile, point_index, gauss_index, shape_index,
                               comps_per_source, source_index, num_coeff_per_shape,
                                num_radec)
            if point:
                point_index = add_point_text(outfile, ra, dec,
                                             point_index)
            
            source_index = check_and_write_source(outfile, point_index, gauss_index, shape_index,
                               comps_per_source, source_index, num_coeff_per_shape,
                                num_radec)
            if gaussian:
                gauss_index = add_gauss_text(outfile, ra, dec,
                                             gauss_index)
            
            source_index = check_and_write_source(outfile, point_index, gauss_index, shape_index,
                               comps_per_source, source_index, num_coeff_per_shape,
                                num_radec)
            if shapelet: 
                shape_index, basis_index = add_shapelet_text(outfile,
                                            ra, dec, shape_index,
                                            num_coeff_per_shape, basis_index)
                
        outfile.write('ENDSOURCE')
                
    return ra_range, dec_range

##based class on an existing class test that includes methods for
##for 
class Test(BaseChunkTest):
    
    def run_test_read_text_skymodel_chunk(self, skymodel_filename, expected_chunks,
                                          max_num_visibilities, num_baselines,
                                          num_freqs, num_time_steps,
                                          lst,
                                          crop_by_component=True):
        
        precision = "double"
        woden_struct_classes = Woden_Struct_Classes(precision)
        woden_settings = ws.Woden_Settings_Python()
        woden_settings.ra0 = 0.0
        woden_settings.dec0 = MWA_LAT
        
        woden_settings.time_res = 1.0
        woden_settings.latitude = -0.46606083776035967
        woden_settings.latitude_obs_epoch_base = -0.46606083776035967
        
        woden_settings.lst_base = lst
        woden_settings.lst_obs_epoch_base = lst
        
        woden_settings.jd_date = 2457278.201145833
        woden_settings.num_time_steps = num_time_steps
        
        woden_settings.do_precession = 1
        lsts, latitudes = ws.setup_lsts_and_phase_centre(woden_settings)
        
        comp_counter = read_text_skymodel.read_text_radec_count_components(skymodel_filename)
        
        comp_counter.print_info()
        
        comp_counter = crop_below_horizon(lst, MWA_LAT,
                                          comp_counter, 
                                          crop_by_component=crop_by_component)
        
        ##Create a chunking map
        chunked_skymodel_maps = create_skymodel_chunk_map(comp_counter,
                                        max_num_visibilities, num_baselines,
                                        num_freqs, num_time_steps)
        chunked_skymodel_maps = chunked_skymodel_maps[0,0]
        
        beamtype = BeamTypes.FEE_BEAM.value
        
        ##TODO if using everybeam, need args to have correct values
        args = Args()
        args.precision = precision

        num_beams = 1
        
        main_table, shape_table, v_table, q_table, u_table, p_table = read_skymodel.get_skymodel_tables(skymodel_filename)

        python_sources = read_fits_skymodel_chunks(args, main_table, shape_table,
                                                    chunked_skymodel_maps,
                                                    num_freqs, num_time_steps,
                                                    beamtype,
                                                    lsts, latitudes,
                                                    v_table, q_table,
                                                    u_table, p_table,
                                                    args.precision)
        
        source_catalogue = create_source_catalogue_from_python_sources(python_sources,
                                                                       woden_struct_classes,
                                                                       beamtype, precision)
        
        check_all_sources(expected_chunks, source_catalogue)
        
    def run_write_model_test_read_text_skymodel_chunk(self, deg_between_comps : int,
                       num_coeff_per_shape : int,
                       comps_per_source : int, lst : float,
                       max_num_visibilities : float):
        """Let's go do a bit"""
        
        num_freqs = 16
        num_baselines = 8128
        num_time_steps = 14
        
        ra_range, dec_range = write_full_test_skymodel_text(deg_between_comps,
                                 num_coeff_per_shape,
                                 comps_per_source)
        
        filename = "test_full_skymodel.txt"
        
        comps_per_chunk = int(np.floor(max_num_visibilities / (num_baselines * num_freqs * num_time_steps)))
        
        
        expec_skymodel_chunks = make_expected_chunks_text(ra_range,
                                                          dec_range,
                             num_coeff_per_shape, comps_per_source,
                             comps_per_chunk, lst=lst)
        
        # for expec_chunk in expec_skymodel_chunks:
        #     expec_chunk._print_things()
        #     print('-----------------')
        
        
        self.run_test_read_text_skymodel_chunk(filename, expec_skymodel_chunks,
                                          max_num_visibilities, num_baselines,
                                          num_freqs, num_time_steps, lst)
    
    def test_the_first(self):
        deg_between_comps = 240
        num_coeff_per_shape = 4
        comps_per_source = 10
        lst = 0.0
        max_num_visibilities = 1e7
        
        self.run_write_model_test_read_text_skymodel_chunk(deg_between_comps,
                                               num_coeff_per_shape,
                                               comps_per_source, lst,
                                               max_num_visibilities)
        print('-----------------------------')
        
    def test_the_second(self):
        deg_between_comps = 10
        num_coeff_per_shape = 4
        num_list_values = 4
        comps_per_source = 11
        lst = 1e-6
        max_num_visibilities = 1e8
        
        self.run_write_model_test_read_text_skymodel_chunk(deg_between_comps,
                                               num_coeff_per_shape,
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
        
        self.run_write_model_test_read_text_skymodel_chunk(deg_between_comps,
                                               num_coeff_per_shape,
                                               comps_per_source, lst,
                                               max_num_visibilities)
        print('-----------------------------')

##Run the test
if __name__ == '__main__':
    unittest.main()
    