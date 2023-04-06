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

from common_skymodel_test import fill_comp_counter, Expec_Counter, BaseChunkTest, write_full_test_skymodel, Expected_Sky_Chunk, Expected_Components

import woden_settings as ws



D2R = np.pi/180.0
# MWA_LATITUDE = -26.7*D2R

MWA_LAT = -26.703319405555554*D2R

##for now, WODEN has three flux types: power law, curved power law, and list
NUM_FLUX_TYPES = 3

RTOL=1e-10

##based class on an existing class test that includes methods for
##for 
class Test(BaseChunkTest):
    
    def run_test_components(self, found_comps, expec_comps,
                            n_powers, n_curves, n_lists,
                            rtol=RTOL):
        
        # print(n_powers, n_curves, n_lists)
        # print("found", found_comps.ras)
        # print("expec", expec_comps.ras)
        
        self.assertTrue(np.allclose(found_comps.ras, expec_comps.ras,
                                    rtol=rtol))
                
        self.assertTrue(np.allclose(found_comps.decs, expec_comps.decs,
                                    rtol=rtol))
        
        if n_powers > 0:
            self.assertTrue(np.allclose(found_comps.power_ref_freqs,
                                    expec_comps.power_ref_freqs, rtol=rtol))
            self.assertTrue(np.allclose(found_comps.power_ref_stokesI,
                                    expec_comps.power_ref_stokesI, rtol=rtol))
            self.assertTrue(np.allclose(found_comps.power_ref_stokesQ,
                                    expec_comps.power_ref_stokesQ, rtol=rtol))
            self.assertTrue(np.allclose(found_comps.power_ref_stokesU,
                                    expec_comps.power_ref_stokesU, rtol=rtol))
            self.assertTrue(np.allclose(found_comps.power_ref_stokesV,
                                    expec_comps.power_ref_stokesV, rtol=rtol))
            self.assertTrue(np.allclose(found_comps.power_SIs,
                                    expec_comps.power_SIs, rtol=rtol))
            
            self.assertTrue(np.allclose(found_comps.power_comp_inds,
                                        expec_comps.power_comp_inds, rtol=rtol))
            
        if n_curves > 0:
            self.assertTrue(np.allclose(found_comps.curve_ref_freqs,
                                    expec_comps.curve_ref_freqs, rtol=rtol))
            self.assertTrue(np.allclose(found_comps.curve_ref_stokesI,
                                    expec_comps.curve_ref_stokesI, rtol=rtol))
            self.assertTrue(np.allclose(found_comps.curve_ref_stokesQ,
                                    expec_comps.curve_ref_stokesQ, rtol=rtol))
            self.assertTrue(np.allclose(found_comps.curve_ref_stokesU,
                                    expec_comps.curve_ref_stokesU, rtol=rtol))
            self.assertTrue(np.allclose(found_comps.curve_ref_stokesV,
                                    expec_comps.curve_ref_stokesV, rtol=rtol))
            self.assertTrue(np.allclose(found_comps.curve_SIs,
                                    expec_comps.curve_SIs, rtol=rtol))
            self.assertTrue(np.allclose(found_comps.curve_qs,
                                    expec_comps.curve_qs, rtol=rtol))
            
            self.assertTrue(np.allclose(found_comps.curve_comp_inds,
                                        expec_comps.curve_comp_inds, rtol=rtol))
            
        if n_lists > 0:
            self.assertTrue(np.allclose(found_comps.list_freqs,
                                        expec_comps.list_freqs, rtol=rtol))
            self.assertTrue(np.allclose(found_comps.list_stokesI,
                                        expec_comps.list_stokesI, rtol=rtol))
            self.assertTrue(np.allclose(found_comps.list_stokesQ,
                                        expec_comps.list_stokesQ, rtol=rtol))
            self.assertTrue(np.allclose(found_comps.list_stokesU,
                                        expec_comps.list_stokesU, rtol=rtol))
            self.assertTrue(np.allclose(found_comps.list_stokesV,
                                        expec_comps.list_stokesV, rtol=rtol))
            
            self.assertTrue(np.allclose(found_comps.list_comp_inds,
                                        expec_comps.list_comp_inds, rtol=rtol))
            self.assertTrue(np.allclose(found_comps.num_list_values,
                                        expec_comps.num_list_values, rtol=rtol))
            self.assertTrue(np.allclose(found_comps.list_start_indexes,
                                        expec_comps.list_start_indexes, rtol=rtol))
    
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
        
        print(lsts)
        
        comp_counter = read_yaml_skymodel.read_yaml_radec_count_components(skymodel_filename)
        
        
        comp_counter = crop_below_horizon(lst, MWA_LAT,
                                          comp_counter, 
                                          crop_by_component=crop_by_component)
        
        ##Create a chunking map
        chunked_skymodel_maps = create_skymodel_chunk_map(comp_counter,
                                        max_num_visibilities, num_baselines,
                                        num_freqs, num_time_steps)
        
        # beamtype = BeamTypes.NO_BEAM.value
        
        beamtype = BeamTypes.FEE_BEAM.value

        source_catalogue = read_yaml_skymodel.read_yaml_skymodel_chunks(
                                              skymodel_filename, chunked_skymodel_maps,
                                              num_freqs, num_time_steps,
                                              beamtype, lsts, MWA_LAT)
        
        for source_ind, expec_chunk in enumerate(expected_chunks):
            # print("DO LE TESTING")
        # for source_ind in range(source_catalogue.num_sources):
            source = source_catalogue.sources[source_ind]
            
            # print("new source who dis--------------------------")
            python_source = _Ctype_Source_Into_Python(source)
            
            # print(source_ind, python_source.n_points, python_source.n_gauss,
            #           python_source.n_shapes)
            
            if python_source.n_points:
                
                n_powers = expec_chunk.n_point_powers
                n_curves = expec_chunk.n_point_curves
                n_lists = expec_chunk.n_point_lists
                
                self.assertEqual(python_source.n_point_powers, n_powers)
                self.assertEqual(python_source.n_point_curves, n_curves)
                self.assertEqual(python_source.n_point_lists, n_lists)
                
                found_comps = python_source.point_components
                expec_comps = expec_chunk.point_components
                
                self.run_test_components(found_comps, expec_comps,
                            n_powers, n_curves, n_lists)
                
            if python_source.n_gauss:
                
                n_powers = expec_chunk.n_gauss_powers
                n_curves = expec_chunk.n_gauss_curves
                n_lists = expec_chunk.n_gauss_lists
                
                self.assertEqual(python_source.n_gauss_powers, n_powers)
                self.assertEqual(python_source.n_gauss_curves, n_curves)
                self.assertEqual(python_source.n_gauss_lists, n_lists)
                
                found_comps = python_source.gauss_components
                expec_comps = expec_chunk.gauss_components
                
                self.run_test_components(found_comps, expec_comps,
                            n_powers, n_curves, n_lists)
                
                self.assertTrue(np.allclose(found_comps.majors,
                                            expec_comps.majors, rtol=RTOL))
                self.assertTrue(np.allclose(found_comps.minors,
                                            expec_comps.minors, rtol=RTOL))
                self.assertTrue(np.allclose(found_comps.pas,
                                            expec_comps.pas, rtol=RTOL))
                
            if python_source.n_shapes:
                
                n_powers = expec_chunk.n_shape_powers
                n_curves = expec_chunk.n_shape_curves
                n_lists = expec_chunk.n_shape_lists
                
                self.assertEqual(python_source.n_shape_powers, n_powers)
                self.assertEqual(python_source.n_shape_curves, n_curves)
                self.assertEqual(python_source.n_shape_lists, n_lists)
                
                found_comps = python_source.shape_components
                expec_comps = expec_chunk.shape_components

                self.run_test_components(found_comps, expec_comps,
                            n_powers, n_curves, n_lists)
                
                self.assertTrue(np.allclose(found_comps.majors,
                                            expec_comps.majors, rtol=RTOL))
                self.assertTrue(np.allclose(found_comps.minors,
                                            expec_comps.minors, rtol=RTOL))
                self.assertTrue(np.allclose(found_comps.pas,
                                            expec_comps.pas, rtol=RTOL))
                
                self.assertTrue(np.allclose(found_comps.n1s,
                                            expec_comps.n1s, rtol=RTOL))
                self.assertTrue(np.allclose(found_comps.n2s,
                                            expec_comps.n2s, rtol=RTOL))
                self.assertTrue(np.allclose(found_comps.shape_coeffs,
                                            expec_comps.shape_coeffs, rtol=RTOL))
                self.assertTrue(np.allclose(found_comps.param_indexes,
                                            expec_comps.param_indexes, rtol=RTOL))

                
    def populate_pointgauss_chunk(self, comp_type : CompTypes, chunk_ind : int,
                                comps_per_chunk : int, n_powers : int,
                                n_curves : int, n_lists : int,
                                num_list_values : int,
                                num_of_each_comp : int, above_horizon : np.ndarray,
                                expec_ra : np.ndarray, expec_dec : np.ndarray,
                                expec_pow_fluxes : np.ndarray,
                                expec_cur_fluxes : np.ndarray) -> Expected_Sky_Chunk:
        
        power_iter = 0
        curve_iter = 0
        list_iter = 0
        
        num_chunk_power = 0
        num_chunk_curve = 0
        num_chunk_list = 0
        
        ##Lower and upper indexes of components covered in this chunk
        lower_comp_ind = chunk_ind * comps_per_chunk
        upper_comp_ind = (chunk_ind + 1) * comps_per_chunk
        
        power_iter, curve_iter, list_iter, num_chunk_power, num_chunk_curve, num_chunk_list = increment_flux_type_counters(power_iter, curve_iter, list_iter, num_chunk_power, num_chunk_curve, num_chunk_list, n_powers, n_curves, n_lists, comps_per_chunk, lower_comp_ind, upper_comp_ind)
        
        expec_chunk = Expected_Sky_Chunk()
        
        if comp_type == CompTypes.POINT:
        
            expec_chunk.init_point_components(num_chunk_power, num_chunk_curve,
                                              num_chunk_list,
                                              num_list_values, comps_per_chunk)
            components = expec_chunk.point_components
            
            expec_orig_pow_inds = np.arange(0, num_of_each_comp*NUM_FLUX_TYPES*3, NUM_FLUX_TYPES*3)[above_horizon]
            
        elif comp_type == CompTypes.GAUSSIAN:
        
            expec_chunk.init_gauss_components(num_chunk_power, num_chunk_curve,
                                              num_chunk_list,
                                              num_list_values, comps_per_chunk)
            components = expec_chunk.gauss_components
            
            
        if num_chunk_power:
            
            low_pow_chunk = 0
            high_pow_chunk = num_chunk_power
            low_pow_coord = power_iter
            high_pow_coord = power_iter+num_chunk_power
            components.ras[low_pow_chunk:high_pow_chunk] = expec_ra[low_pow_coord:high_pow_coord]
            components.decs[low_pow_chunk:high_pow_chunk] = expec_dec[low_pow_coord:high_pow_coord]
            
            components.power_ref_freqs[:num_chunk_power] = expec_pow_fluxes[low_pow_coord:high_pow_coord]
            components.power_ref_stokesI[:num_chunk_power] = expec_pow_fluxes[low_pow_coord:high_pow_coord]
            components.power_ref_stokesQ[:num_chunk_power] = expec_pow_fluxes[low_pow_coord:high_pow_coord]
            components.power_ref_stokesU[:num_chunk_power] = expec_pow_fluxes[low_pow_coord:high_pow_coord]
            components.power_ref_stokesV[:num_chunk_power] = expec_pow_fluxes[low_pow_coord:high_pow_coord]
            components.power_SIs[:num_chunk_power] = expec_pow_fluxes[low_pow_coord:high_pow_coord]
            
            components.power_comp_inds = np.arange(num_chunk_power)
            
            if comp_type == CompTypes.GAUSSIAN:
                components.majors[low_pow_chunk:high_pow_chunk] = expec_pow_fluxes[low_pow_coord:high_pow_coord]*(D2R/3600.0)
                components.minors[low_pow_chunk:high_pow_chunk] = expec_pow_fluxes[low_pow_coord:high_pow_coord]*(D2R/3600.0)
                components.pas[low_pow_chunk:high_pow_chunk] = expec_pow_fluxes[low_pow_coord:high_pow_coord]*D2R
            
        if num_chunk_curve:
            low_cur_chunk = num_chunk_power
            high_cur_chunk = num_chunk_power+num_chunk_curve
            low_cur_coord = curve_iter
            high_cur_coord = curve_iter+num_chunk_curve
            
            components.ras[low_cur_chunk:high_cur_chunk] = expec_ra[low_cur_coord:high_cur_coord]
            components.decs[low_cur_chunk:high_cur_chunk] = expec_dec[low_cur_coord:high_cur_coord]
            
            components.curve_ref_freqs[:num_chunk_curve] = expec_cur_fluxes[low_cur_coord:high_cur_coord]
            components.curve_ref_stokesI[:num_chunk_curve] = expec_cur_fluxes[low_cur_coord:high_cur_coord]
            components.curve_ref_stokesQ[:num_chunk_curve] = expec_cur_fluxes[low_cur_coord:high_cur_coord]
            components.curve_ref_stokesU[:num_chunk_curve] = expec_cur_fluxes[low_cur_coord:high_cur_coord]
            components.curve_ref_stokesV[:num_chunk_curve] = expec_cur_fluxes[low_cur_coord:high_cur_coord]
            components.curve_SIs[:num_chunk_curve] = expec_cur_fluxes[low_cur_coord:high_cur_coord]
            components.curve_qs[:num_chunk_curve] = expec_cur_fluxes[low_cur_coord:high_cur_coord]
            
            components.curve_comp_inds = np.arange(num_chunk_curve) + num_chunk_power
            
            if comp_type == CompTypes.GAUSSIAN:
                components.majors[low_cur_chunk:high_cur_chunk] = expec_cur_fluxes[low_cur_coord:high_cur_coord]*(D2R/3600.0)
                components.minors[low_cur_chunk:high_cur_chunk] = expec_cur_fluxes[low_cur_coord:high_cur_coord]*(D2R/3600.0)
                components.pas[low_cur_chunk:high_cur_chunk] = expec_cur_fluxes[low_cur_coord:high_cur_coord]*D2R
        
        if num_chunk_list:
            low_lis_chunk = num_chunk_power + num_chunk_curve
            high_lis_chunk = num_chunk_power + num_chunk_curve + num_chunk_list
            low_lis_coord = list_iter
            high_lis_coord = list_iter+num_chunk_list
            
            components.ras[low_lis_chunk:high_lis_chunk] = expec_ra[low_lis_coord:high_lis_coord]
            components.decs[low_lis_chunk:high_lis_chunk] = expec_dec[low_lis_coord:high_lis_coord]
            
            expec_flux = np.empty(num_list_values*num_chunk_list)
            ##OK, so when writing the test sky model, we constantly iterate
            ##the flux values for each new flux entry and each new component
            ##we iterate over ALL component types, and add components in
            ##order of POINT, GAUSSIAN, SHAPELET. So to get expected flux
            ##entries, we have to offset via what is above horizon, how
            ##many flux entries per component, and what comp_type we have
            
            for flux_ind, flux_init in enumerate(above_horizon[low_lis_coord:high_lis_coord]):
                
                low = flux_ind*num_list_values
                high = (flux_ind + 1)*num_list_values
                
                base_flux_vals = np.arange(num_list_values) + flux_init*num_list_values*NUM_FLUX_TYPES
            
                if comp_type == CompTypes.POINT:
                    expec_flux[low:high] = base_flux_vals
                    
                elif comp_type == CompTypes.GAUSSIAN:
                    expec_flux[low:high] = base_flux_vals + num_list_values
                    
            components.list_freqs = expec_flux
            components.list_stokesI = expec_flux
            components.list_stokesQ = expec_flux
            components.list_stokesU = expec_flux
            components.list_stokesV = expec_flux
            
            components.list_comp_inds = np.arange(num_chunk_list) + num_chunk_power + num_chunk_curve
            
            components.num_list_values = np.full(num_chunk_list, num_list_values)
            components.list_start_indexes = np.arange(0, num_list_values*num_chunk_list, num_list_values)
            
            if comp_type == CompTypes.GAUSSIAN:
                
                expec_params = expec_pow_fluxes + 2
                
                components.majors[low_lis_chunk:high_lis_chunk] = expec_params[low_lis_coord:high_lis_coord]*(D2R/3600.0)
                components.minors[low_lis_chunk:high_lis_chunk] = expec_params[low_lis_coord:high_lis_coord]*(D2R/3600.0)
                components.pas[low_lis_chunk:high_lis_chunk] = expec_params[low_lis_coord:high_lis_coord]*D2R
        
        return expec_chunk
    
    def populate_shapelet_chunk(self,
                                expec_chunk : Expected_Sky_Chunk,
                                num_list_values : int, num_coeff_per_shape : int,
                                above_horizon : np.ndarray,
                                expec_ra : np.ndarray,
                                expec_dec : np.ndarray,
                                comp_inds : np.ndarray,
                                orig_comp_inds : np.ndarray,
                                chunk_basis_param_indexes : np.ndarray,
                                chunk_basis_values : np.ndarray):
        
        # print(orig_comp_inds)
        num_crop_comp = len(expec_ra)

        components = expec_chunk.shape_components
        
        n_powers = expec_chunk.n_shape_powers
        n_curves = expec_chunk.n_shape_curves
        n_lists = expec_chunk.n_shape_lists
        
        power_comp_ind = 0
        curve_comp_ind = 0
        list_comp_ind = 0

        for new_comp_ind, old_comp_ind in enumerate(comp_inds):
            comp_type_ind = int(old_comp_ind // num_crop_comp)
            coord_expec_ind = int(old_comp_ind % num_crop_comp)
            orig_ind = orig_comp_inds[new_comp_ind]

            components.ras[new_comp_ind] = expec_ra[coord_expec_ind]
            components.decs[new_comp_ind] = expec_dec[coord_expec_ind]
            
            ##this number is equal to the index of the shapelet within
            ##all shapelets in the original catalogue
            components.majors[new_comp_ind] = orig_ind*(D2R/3600.0)
            components.minors[new_comp_ind] = orig_ind*(D2R/3600.0)
            components.pas[new_comp_ind] = orig_ind*D2R
            
            
            ##these are the basis function values for this particular chunk
            ##HOWEVER in the code that reads things in, it doesn't care
            ##about the power,curve,list ordering because there is an
            ##array to do that ordering (param_indexes). So these things
            ##are read in as they appear in the srclist, meaning we
            ##need to reorded our prediction by argsorting via the basis
            ##values (as they increment with appearence.) ffffffuuuuuuuu
            
            reorder = np.argsort(chunk_basis_values)
            components.n1s = chunk_basis_values[reorder]
            components.n2s = chunk_basis_values[reorder]
            components.shape_coeffs = chunk_basis_values[reorder]
            components.param_indexes = chunk_basis_param_indexes[reorder]
            

            ##This means we have a power law source
            if old_comp_ind < num_crop_comp:
                components.power_ref_freqs[power_comp_ind] = orig_ind
                components.power_ref_stokesI[power_comp_ind] = orig_ind
                components.power_ref_stokesQ[power_comp_ind] = orig_ind
                components.power_ref_stokesU[power_comp_ind] = orig_ind
                components.power_ref_stokesV[power_comp_ind] = orig_ind
                components.power_SIs[power_comp_ind] = orig_ind
                
                components.power_comp_inds[power_comp_ind] = new_comp_ind
                power_comp_ind += 1
            
            ##This means we have a curved power law source
            elif old_comp_ind < 2*num_crop_comp:
                components.curve_ref_freqs[curve_comp_ind] = orig_ind
                components.curve_ref_stokesI[curve_comp_ind] = orig_ind
                components.curve_ref_stokesQ[curve_comp_ind] = orig_ind
                components.curve_ref_stokesU[curve_comp_ind] = orig_ind
                components.curve_ref_stokesV[curve_comp_ind] = orig_ind
                components.curve_SIs[curve_comp_ind] = orig_ind
                components.curve_qs[curve_comp_ind] = orig_ind
                
                components.curve_comp_inds[curve_comp_ind] = new_comp_ind
                curve_comp_ind += 1
                
            else:
                ##RIGHT in the original sky model, we keep iterating 
                ##the flux values of all point, gaussian, and shapelets,
                ##so we can do some maths here to work out what the flux
                ##values should be
                flux_start = num_list_values*orig_ind
                
                for flux_ind in range(num_list_values):
                    
                    list_ind = list_comp_ind*num_list_values

                    components.list_freqs[list_ind + flux_ind] = flux_start + flux_ind
                    components.list_stokesI[list_ind + flux_ind] = flux_start + flux_ind
                    components.list_stokesQ[list_ind + flux_ind] = flux_start + flux_ind
                    components.list_stokesU[list_ind + flux_ind] = flux_start + flux_ind
                    components.list_stokesV[list_ind + flux_ind] = flux_start + flux_ind
                    
                components.list_comp_inds[list_comp_ind] = new_comp_ind
                
                components.num_list_values[list_comp_ind] = num_list_values
                components.list_start_indexes[list_comp_ind] = list_comp_ind*num_list_values
                
                
                list_comp_ind += 1
                
        return expec_chunk
        
        
    def make_expected_chunks(self, ra_range, dec_range,
                             num_coeff_per_shape, num_list_values,
                             comps_per_source, comps_per_chunk, lst = 0.0):
        
        num_of_each_comp = len(ra_range)
        
        # comp_index_range = np.arange(num_of_each_comp)
        coeff_range = np.arange(NUM_FLUX_TYPES*num_of_each_comp*num_coeff_per_shape)
        flux_list_range = np.arange(NUM_FLUX_TYPES*num_of_each_comp*num_list_values)
        
        has = lst - ra_range
        azs, els = erfa.hd2ae(has, dec_range, MWA_LAT)
        
        np.savez('/home/jline/software/WODEN_dev/cmake_testing/wodenpy/skymodel/test_sky_crop.npz', comp_has=has, comp_decs=dec_range,
             comp_azs=azs, comp_els=els)
        
        ##Crop below the horizon based on COMPONENT/SOURCE
        # include_flags = np.zeros(len(all_comp_types))
        above_horizon = np.where(els >= 0)[0]
        
        expec_ra = ra_range[above_horizon]
        expec_dec = dec_range[above_horizon]
        
        print(f"There are {len(expec_ra)} coords above horizon")
        
        expec_pow_fluxes = np.arange(0, num_of_each_comp*NUM_FLUX_TYPES, NUM_FLUX_TYPES)[above_horizon]
        expec_cur_fluxes = np.arange(1, num_of_each_comp*NUM_FLUX_TYPES, NUM_FLUX_TYPES)[above_horizon]
        
        num_crop_comp = len(above_horizon)
        
        num_point_chunks = int(np.ceil((NUM_FLUX_TYPES*num_crop_comp) / comps_per_chunk))
        num_gauss_chunks = int(np.ceil((NUM_FLUX_TYPES*num_crop_comp) / comps_per_chunk))
        num_coeff_chunks = int(np.ceil((NUM_FLUX_TYPES*num_coeff_per_shape*num_crop_comp) / comps_per_chunk))
        
        expec_skymodel_chunks = []
        
        ##start with the POINTS
        
        n_powers = num_crop_comp
        n_curves = num_crop_comp
        n_lists = num_crop_comp
        
        for chunk_ind in range(num_point_chunks):
            expec_chunk = self.populate_pointgauss_chunk(CompTypes.POINT, chunk_ind,
                                comps_per_chunk, n_powers,
                                n_curves, n_lists, num_list_values,
                                num_of_each_comp, above_horizon,
                                expec_ra, expec_dec,
                                expec_pow_fluxes, expec_cur_fluxes)
            
            expec_skymodel_chunks.append(expec_chunk)
            
        for chunk_ind in range(num_gauss_chunks):
            expec_chunk = self.populate_pointgauss_chunk(CompTypes.GAUSSIAN, chunk_ind,
                                comps_per_chunk, n_powers,
                                n_curves, n_lists, num_list_values,
                                num_of_each_comp, above_horizon,
                                expec_ra, expec_dec,
                                expec_pow_fluxes, expec_cur_fluxes)
            
            expec_skymodel_chunks.append(expec_chunk)
            
        total_shape_basis = num_coeff_per_shape*num_crop_comp*NUM_FLUX_TYPES
        
        
        ##OK OK OK so when chunking, the code resets the order to be
        ##all power first, then all curve second, then all list third.
        ##So even though it's written in the skymodel as power, curve, list,
        ##power, curve, list etc write a different expected order here
        
        shape_comp_ind_to_comp_type = np.empty(num_crop_comp*NUM_FLUX_TYPES, dtype=CompTypes)
        
        shape_comp_ind_to_comp_type[:num_crop_comp] = CompTypes.SHAPE_POWER
        shape_comp_ind_to_comp_type[num_crop_comp:2*num_crop_comp] = CompTypes.SHAPE_CURVE
        shape_comp_ind_to_comp_type[-num_crop_comp:] = CompTypes.SHAPE_LIST
        
        shape_basis_to_comp_ind = np.repeat(np.arange(num_crop_comp*NUM_FLUX_TYPES), num_coeff_per_shape)
        
        shape_expec_pow_fluxes = np.repeat(expec_pow_fluxes, num_coeff_per_shape)
        shape_expec_cur_fluxes = np.repeat(expec_cur_fluxes, num_coeff_per_shape)

        orig_comp_inds = np.empty(num_crop_comp*NUM_FLUX_TYPES)

        orig_comp_inds[:num_crop_comp] = above_horizon*NUM_FLUX_TYPES
        orig_comp_inds[num_crop_comp:2*num_crop_comp] = above_horizon*NUM_FLUX_TYPES + 1
        orig_comp_inds[-num_crop_comp:] = above_horizon*NUM_FLUX_TYPES + 2
        
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
            
            ##this slots in the curved power law components
            low_coeff = orig_ind*num_coeff_per_shape*NUM_FLUX_TYPES + num_coeff_per_shape
            low_ind = (comp_ind + num_crop_comp)*num_coeff_per_shape
            high_ind = low_ind + num_coeff_per_shape
            shape_basis_values[low_ind:high_ind] = range(low_coeff, low_coeff + num_coeff_per_shape)
            
            ##this slots in the curved list components
            low_coeff = orig_ind*num_coeff_per_shape*NUM_FLUX_TYPES + 2*num_coeff_per_shape
            low_ind = (comp_ind + 2*num_crop_comp)*num_coeff_per_shape
            high_ind = low_ind + num_coeff_per_shape
            shape_basis_values[low_ind:high_ind] = range(low_coeff, low_coeff + num_coeff_per_shape)
            
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
            curve_inds = np.unique(np.where(shape_comp_ind_to_comp_type[comp_inds] == CompTypes.SHAPE_CURVE)[0])
            list_inds = np.unique(np.where(shape_comp_ind_to_comp_type[comp_inds] == CompTypes.SHAPE_LIST)[0])
            
            num_chunk_power = len(power_inds)
            num_chunk_curve = len(curve_inds)
            num_chunk_list = len(list_inds)
            
            expec_chunk = Expected_Sky_Chunk()
            
            expec_chunk.init_shape_components(num_chunk_power, num_chunk_curve,
                                              num_chunk_list,
                                              num_list_values, comps_per_chunk)
            
            
            chunk_basis_values = shape_basis_values[basis_inds]
            
            chunk_basis_param_indexes = shape_basis_to_comp_ind[basis_inds]
            chunk_basis_param_indexes -= chunk_basis_param_indexes[0]
            
            expec_chunk = self.populate_shapelet_chunk(expec_chunk,
                                num_list_values, num_coeff_per_shape,
                                above_horizon, expec_ra, expec_dec,
                                comp_inds, orig_comp_inds[comp_inds],
                                chunk_basis_param_indexes,
                                chunk_basis_values)
            
            expec_skymodel_chunks.append(expec_chunk)
        
        return expec_skymodel_chunks
        
        
    
    def run_write_model_test_read_yaml_skymodel_chunk(self, deg_between_comps : int,
                       num_coeff_per_shape : int, num_list_values : int,
                       comps_per_source : int, lst : float,
                       max_num_visibilities : float):
        """Let's go do a bit"""
        
        num_freqs = 16
        num_baselines = 8128
        num_time_steps = 14
        
        ra_range, dec_range = write_full_test_skymodel(deg_between_comps,
                                 num_coeff_per_shape,
                                 num_list_values,
                                 comps_per_source)
        
        # filename = f"{code_dir}/test_full_skymodel.yaml"
        filename = "test_full_skymodel.yaml"
        
        comps_per_chunk = int(np.floor(max_num_visibilities / (num_baselines * num_freqs * num_time_steps)))
        
        
        expec_skymodel_chunks = self.make_expected_chunks(ra_range, dec_range,
                             num_coeff_per_shape, num_list_values,
                             comps_per_source, comps_per_chunk)
        
        
        self.run_test_read_yaml_skymodel_chunk(filename, expec_skymodel_chunks,
                                          max_num_visibilities, num_baselines,
                                          num_freqs, num_time_steps, lst)
    
    def test_the_first(self):
        deg_between_comps = 240
        num_coeff_per_shape = 4
        num_list_values = 4
        comps_per_source = 10
        lst = 0.0
        max_num_visibilities = 1e7
        
        self.run_write_model_test_read_yaml_skymodel_chunk(deg_between_comps,
                                               num_coeff_per_shape,
                                               num_list_values,
                                               comps_per_source, lst,
                                               max_num_visibilities)
        print('-----------------------------')
        
    # def test_the_second(self):
    #     deg_between_comps = 10
    #     num_coeff_per_shape = 4
    #     num_list_values = 4
    #     comps_per_source = 11
    #     lst = 0.0
    #     max_num_visibilities = 1e8
        
    #     self.run_write_model_test_read_yaml_skymodel_chunk(deg_between_comps,
    #                                            num_coeff_per_shape,
    #                                            num_list_values,
    #                                            comps_per_source, lst,
    #                                            max_num_visibilities)
    #     print('-----------------------------')
        
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
    