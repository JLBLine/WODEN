from sys import path
import os
import unittest
import numpy as np
import erfa
import numpy.testing as npt

# ##Code we are testing
from wodenpy.skymodel import read_yaml_skymodel
# import wodenpy
from wodenpy.skymodel.woden_skymodel import Component_Type_Counter, CompTypes, crop_below_horizon
from wodenpy.skymodel.chunk_sky_model import create_skymodel_chunk_map, map_chunk_pointgauss, Skymodel_Chunk_Map, increment_flux_type_counters

from wodenpy.use_libwoden.beam_settings import BeamTypes

from wodenpy.use_libwoden.skymodel_structs import setup_chunked_source, _Ctype_Source_Into_Python

from common_skymodel_test import fill_comp_counter, Expec_Counter, BaseChunkTest, Expected_Sky_Chunk, Expected_Components

import wodenpy.use_libwoden.woden_settings as ws



D2R = np.pi/180.0
# MWA_LATITUDE = -26.7*D2R

MWA_LAT = -26.703319405555554*D2R

##for now, WODEN has three flux types: power law, curved power law, and list
NUM_FLUX_TYPES = 3

RTOL=1e-10

##limits for "all sky" sky model
LOW_DEC = -90.0*D2R
HIGH_DEC = 30.0*D2R

def check_components(found_comps, expec_comps,
                     n_powers, n_curves, n_lists,
                     rtol=RTOL,fits_skymodel=True):
        
    # print(n_powers, n_curves, n_lists)
    # print("found", found_comps.ras)
    # print("expec", expec_comps.ras)
    
    
    
    npt.assert_allclose(found_comps.ras, expec_comps.ras,
                                rtol=rtol)
            
    npt.assert_allclose(found_comps.decs, expec_comps.decs,
                                rtol=rtol)
    
    if n_powers > 0:
        
        npt.assert_allclose(found_comps.power_ref_freqs,
                                expec_comps.power_ref_freqs, rtol=rtol)
        npt.assert_allclose(found_comps.power_ref_stokesI,
                                expec_comps.power_ref_stokesI, rtol=rtol)
        
        # if not fits_skymodel:
        #     npt.assert_allclose(found_comps.power_ref_stokesQ,
        #                             expec_comps.power_ref_stokesQ, rtol=rtol)
        #     npt.assert_allclose(found_comps.power_ref_stokesU,
        #                             expec_comps.power_ref_stokesU, rtol=rtol)
        #     npt.assert_allclose(found_comps.power_ref_stokesV,
        #                             expec_comps.power_ref_stokesV, rtol=rtol)
        npt.assert_allclose(found_comps.power_SIs,
                                expec_comps.power_SIs, rtol=rtol)
        
        npt.assert_allclose(found_comps.power_comp_inds,
                                    expec_comps.power_comp_inds, rtol=rtol)
        
    if n_curves > 0:
        npt.assert_allclose(found_comps.curve_ref_freqs,
                                expec_comps.curve_ref_freqs, rtol=rtol)
        npt.assert_allclose(found_comps.curve_ref_stokesI,
                                expec_comps.curve_ref_stokesI, rtol=rtol)
        # if not fits_skymodel:
        #     npt.assert_allclose(found_comps.curve_ref_stokesQ,
        #                             expec_comps.curve_ref_stokesQ, rtol=rtol)
        #     npt.assert_allclose(found_comps.curve_ref_stokesU,
        #                             expec_comps.curve_ref_stokesU, rtol=rtol)
        #     npt.assert_allclose(found_comps.curve_ref_stokesV,
        #                             expec_comps.curve_ref_stokesV, rtol=rtol)
        npt.assert_allclose(found_comps.curve_SIs,
                                expec_comps.curve_SIs, rtol=rtol)
        npt.assert_allclose(found_comps.curve_qs,
                                expec_comps.curve_qs, rtol=rtol)
        
        npt.assert_allclose(found_comps.curve_comp_inds,
                                    expec_comps.curve_comp_inds, rtol=rtol)
        
    if n_lists > 0:
        npt.assert_allclose(found_comps.list_freqs,
                                    expec_comps.list_freqs, rtol=rtol)
        npt.assert_allclose(found_comps.list_stokesI,
                                    expec_comps.list_stokesI, rtol=rtol)
        
        # if not fits_skymodel:
        #     npt.assert_allclose(found_comps.list_stokesQ,
        #                                 expec_comps.list_stokesQ, rtol=rtol)
        #     npt.assert_allclose(found_comps.list_stokesU,
        #                                 expec_comps.list_stokesU, rtol=rtol)
        #     npt.assert_allclose(found_comps.list_stokesV,
        #                                 expec_comps.list_stokesV, rtol=rtol)
        
        npt.assert_allclose(found_comps.list_comp_inds,
                                    expec_comps.list_comp_inds, rtol=rtol)
        npt.assert_allclose(found_comps.num_list_values,
                                    expec_comps.num_list_values, rtol=rtol)
        npt.assert_allclose(found_comps.list_start_indexes,
                                    expec_comps.list_start_indexes, rtol=rtol)


def check_all_sources(expected_chunks, source_catalogue,
                      fits_skymodel=True):
    
    rtol = RTOL

    for source_ind, expec_chunk in enumerate(expected_chunks):
        source = source_catalogue.sources[source_ind]
        
        python_source = _Ctype_Source_Into_Python(source)
        
        # print(source_ind, python_source.n_points, python_source.n_gauss,
        #           python_source.n_shapes)
        
        if python_source.n_points:
            
            n_powers = expec_chunk.n_point_powers
            n_curves = expec_chunk.n_point_curves
            n_lists = expec_chunk.n_point_lists
            
            npt.assert_equal(python_source.n_point_powers, n_powers)
            npt.assert_equal(python_source.n_point_curves, n_curves)
            npt.assert_equal(python_source.n_point_lists, n_lists)
            
            found_comps = python_source.point_components
            expec_comps = expec_chunk.point_components
            
            check_components(found_comps, expec_comps,
                        n_powers, n_curves, n_lists,
                        fits_skymodel=fits_skymodel)
            
        if python_source.n_gauss:
            
            n_powers = expec_chunk.n_gauss_powers
            n_curves = expec_chunk.n_gauss_curves
            n_lists = expec_chunk.n_gauss_lists
            
            npt.assert_equal(python_source.n_gauss_powers, n_powers)
            npt.assert_equal(python_source.n_gauss_curves, n_curves)
            npt.assert_equal(python_source.n_gauss_lists, n_lists)
            
            found_comps = python_source.gauss_components
            expec_comps = expec_chunk.gauss_components
            
            check_components(found_comps, expec_comps,
                        n_powers, n_curves, n_lists,
                        fits_skymodel=fits_skymodel)
            
            npt.assert_allclose(found_comps.majors,
                                        expec_comps.majors, rtol=rtol)
            npt.assert_allclose(found_comps.minors,
                                        expec_comps.minors, rtol=rtol)
            npt.assert_allclose(found_comps.pas,
                                        expec_comps.pas, rtol=rtol)
            
        if python_source.n_shapes:
            
            n_powers = expec_chunk.n_shape_powers
            n_curves = expec_chunk.n_shape_curves
            n_lists = expec_chunk.n_shape_lists
            
            npt.assert_equal(python_source.n_shape_powers, n_powers)
            npt.assert_equal(python_source.n_shape_curves, n_curves)
            npt.assert_equal(python_source.n_shape_lists, n_lists)
            
            found_comps = python_source.shape_components
            expec_comps = expec_chunk.shape_components

            check_components(found_comps, expec_comps,
                        n_powers, n_curves, n_lists,
                        fits_skymodel=fits_skymodel)
            
            npt.assert_allclose(found_comps.majors,
                                        expec_comps.majors, rtol=rtol)
            npt.assert_allclose(found_comps.minors,
                                        expec_comps.minors, rtol=rtol)
            npt.assert_allclose(found_comps.pas,
                                        expec_comps.pas, rtol=rtol)
            
            npt.assert_allclose(found_comps.n1s,
                                        expec_comps.n1s, rtol=rtol)
            npt.assert_allclose(found_comps.n2s,
                                        expec_comps.n2s, rtol=rtol)
            npt.assert_allclose(found_comps.shape_coeffs,
                                        expec_comps.shape_coeffs, rtol=rtol)
            npt.assert_allclose(found_comps.param_indexes,
                                        expec_comps.param_indexes, rtol=rtol)

            
def populate_pointgauss_chunk(comp_type : CompTypes, chunk_ind : int,
                            comps_per_chunk : int, n_powers : int,
                            n_curves : int, n_lists : int,
                            num_list_values : int,
                            num_of_each_comp : int, above_horizon : np.ndarray,
                            expec_ra : np.ndarray, expec_dec : np.ndarray,
                            expec_pow_fluxes : np.ndarray,
                            expec_cur_fluxes : np.ndarray,
                            fits_skymodel = False) -> Expected_Sky_Chunk:
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
        
        # ##power law for FITS sky model is locked to a reference freq of 200MHz
        components.power_ref_freqs[:num_chunk_power] = 200e+6
            
        ##Actual FITS model is read in at 200MHz
        if fits_skymodel:
            components.power_ref_stokesI[:num_chunk_power] = expec_pow_fluxes[low_pow_coord:high_pow_coord]
            components.power_SIs[:num_chunk_power] = expec_pow_fluxes[low_pow_coord:high_pow_coord]
        ##Anything can be written a some other reference frequency, which
        ##we should be extrapolating to 200MHz
        else:
            ref_freqs = 100e+6 + (expec_pow_fluxes[low_pow_coord:high_pow_coord] + 1)*1e+4
            components.power_ref_stokesI[:num_chunk_power] = expec_pow_fluxes[low_pow_coord:high_pow_coord]*(200e+6 / ref_freqs)**(expec_pow_fluxes[low_pow_coord:high_pow_coord]/100.0)
            components.power_SIs[:num_chunk_power] = expec_pow_fluxes[low_pow_coord:high_pow_coord]/100.0
        
        components.power_ref_stokesQ[:num_chunk_power] = expec_pow_fluxes[low_pow_coord:high_pow_coord]
        components.power_ref_stokesU[:num_chunk_power] = expec_pow_fluxes[low_pow_coord:high_pow_coord]
        components.power_ref_stokesV[:num_chunk_power] = expec_pow_fluxes[low_pow_coord:high_pow_coord]
        
        
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
        
        ##curved power law for FITS sky model is locked to a reference freq of 200MHz
        if fits_skymodel:
            components.curve_ref_stokesI[:num_chunk_curve] = expec_cur_fluxes[low_cur_coord:high_cur_coord]
            components.curve_SIs[:num_chunk_curve] = expec_cur_fluxes[low_cur_coord:high_cur_coord]
            components.curve_qs[:num_chunk_curve] = expec_cur_fluxes[low_cur_coord:high_cur_coord]
        else:
            
            ref_freqs = 100e+6 + (expec_cur_fluxes[low_cur_coord:high_cur_coord] + 1)*1e+4
            sis =  expec_cur_fluxes[low_cur_coord:high_cur_coord] / 100.0
            curve_qs =  expec_cur_fluxes[low_cur_coord:high_cur_coord] / 100.0
            
            
            si_ratio = (200e+6 / ref_freqs)**sis
            exp_bit = np.exp(curve_qs*np.log(200e+6 / ref_freqs)**2)
            expec_Is = expec_cur_fluxes[low_cur_coord:high_cur_coord]*si_ratio*exp_bit
            
            
            logratio = np.log(ref_freqs / 200e+6)
            new_sis = (np.log(expec_cur_fluxes[low_cur_coord:high_cur_coord] / expec_Is) - curve_qs*logratio**2) / logratio
            
            
            components.curve_ref_stokesI[:num_chunk_curve] = expec_Is
            components.curve_SIs[:num_chunk_curve] = new_sis
            components.curve_qs[:num_chunk_curve] = curve_qs
            
        components.curve_ref_freqs[:num_chunk_curve] = 200e+6
        components.curve_ref_stokesQ[:num_chunk_curve] = expec_cur_fluxes[low_cur_coord:high_cur_coord]
        components.curve_ref_stokesU[:num_chunk_curve] = expec_cur_fluxes[low_cur_coord:high_cur_coord]
        components.curve_ref_stokesV[:num_chunk_curve] = expec_cur_fluxes[low_cur_coord:high_cur_coord]
        
        
        
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
            
        ##OK, so for the FITS skymodel format, we can't have a new freq for
        ##every single list entry, as each different freq has a column. So
        ##we have a different behaviour where we just have as many columns as
        ##`num_list_values`. Title the column by index, which is read in as MHz,
        ##so need to mulitply expected vaulues by 1e+6
        components.list_freqs = np.tile(np.arange(num_list_values)*1e+6, num_chunk_list)
            
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

def populate_shapelet_chunk(expec_chunk : Expected_Sky_Chunk,
                            num_list_values : int, num_coeff_per_shape : int,
                            above_horizon : np.ndarray,
                            expec_ra : np.ndarray,
                            expec_dec : np.ndarray,
                            comp_inds : np.ndarray,
                            orig_comp_inds : np.ndarray,
                            chunk_basis_param_indexes : np.ndarray,
                            chunk_basis_values : np.ndarray,
                            fits_skymodel=False):
    
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
        
        # ##these are the basis function values for this particular chunk
        components.n1s = chunk_basis_values
        components.n2s = chunk_basis_values
        components.shape_coeffs = chunk_basis_values
        components.param_indexes = chunk_basis_param_indexes
        # else:
            

        ##This means we have a power law source
        if old_comp_ind < num_crop_comp:
            ##Actual FITS model is read in at 200MHz
            if fits_skymodel:
                components.power_ref_stokesI[power_comp_ind] = orig_ind
                components.power_SIs[power_comp_ind] = orig_ind
            ##Anything can be written a some other reference frequency, which
            ##we should be extrapolating to 200MHz
            else:
                ref_freq = 100e+6 + (orig_ind + 1)*1e+4
                components.power_ref_stokesI[power_comp_ind] = orig_ind*(200e+6 / ref_freq)**(orig_ind/100.0)
                components.power_SIs[power_comp_ind] = orig_ind / 100.0
                
                
            components.power_ref_freqs[power_comp_ind] = 200e+6
            components.power_ref_stokesQ[power_comp_ind] = orig_ind
            components.power_ref_stokesU[power_comp_ind] = orig_ind
            components.power_ref_stokesV[power_comp_ind] = orig_ind
            
            components.power_comp_inds[power_comp_ind] = new_comp_ind
            power_comp_ind += 1
        
        ##This means we have a curved power law source
        elif old_comp_ind < 2*num_crop_comp:
            if fits_skymodel:
                components.curve_ref_stokesI[curve_comp_ind] = orig_ind
                components.curve_SIs[curve_comp_ind] = orig_ind
                components.curve_qs[curve_comp_ind] = orig_ind
            else:
                ref_freq = 100e+6 + (orig_ind + 1)*1e+4
                si =  orig_ind / 100.0
                curve_q =  orig_ind / 100.0
                
                si_ratio = (200e+6 / ref_freq)**si
                exp_bit = np.exp(curve_q*np.log(200e+6 / ref_freq)**2)
                expec_I = orig_ind*si_ratio*exp_bit
                
                logratio = np.log(ref_freq / 200e+6)
                new_si = (np.log(orig_ind / expec_I) - curve_q*logratio**2) / logratio
                
                components.curve_ref_stokesI[curve_comp_ind] = expec_I
                components.curve_SIs[curve_comp_ind] = new_si
                components.curve_qs[curve_comp_ind] = curve_q
                
            components.curve_ref_freqs[curve_comp_ind] = 200e+6    
            components.curve_ref_stokesQ[curve_comp_ind] = orig_ind
            components.curve_ref_stokesU[curve_comp_ind] = orig_ind
            components.curve_ref_stokesV[curve_comp_ind] = orig_ind
            
            components.curve_comp_inds[curve_comp_ind] = new_comp_ind
            curve_comp_ind += 1
            
        else:
            ##RIGHT in the original sky model, we keep iterating 
            ##the flux values of all point, gaussian, and shapelets,
            ##so we can do some maths here to work out what the flux
            ##values should be
            flux_start = num_list_values*orig_ind
            
            for flux_ind in range(num_list_values):
                
                ##OK, so for the FITS skymodel format, we can't have a new freq for
                ##every single list entry, as each different freq has a column. So
                ##we have a different behaviour where we just have as many columns as
                ##`num_list_values`. Title the column by index, which is read in as MHz,
                ##so need to mulitply expected vaulues by 1e+6
                freq_value = flux_ind*1e+6
                list_ind = list_comp_ind*num_list_values

                components.list_freqs[list_ind + flux_ind] = freq_value
                components.list_stokesI[list_ind + flux_ind] = flux_start + flux_ind
                components.list_stokesQ[list_ind + flux_ind] = flux_start + flux_ind
                components.list_stokesU[list_ind + flux_ind] = flux_start + flux_ind
                components.list_stokesV[list_ind + flux_ind] = flux_start + flux_ind
                
        
            components.list_comp_inds[list_comp_ind] = new_comp_ind
            components.num_list_values[list_comp_ind] = num_list_values
            components.list_start_indexes[list_comp_ind] = list_comp_ind*num_list_values
            
            list_comp_ind += 1
            
    return expec_chunk
    
    
def make_expected_chunks(ra_range, dec_range,
                         num_coeff_per_shape, num_list_values,
                         comps_per_source, comps_per_chunk, lst = 0.0,
                         fits_skymodel = False):
    
    num_of_each_comp = len(ra_range)
    
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
        expec_chunk = populate_pointgauss_chunk(CompTypes.POINT, chunk_ind,
                            comps_per_chunk, n_powers,
                            n_curves, n_lists, num_list_values,
                            num_of_each_comp, above_horizon,
                            expec_ra, expec_dec,
                            expec_pow_fluxes, expec_cur_fluxes,
                            fits_skymodel=fits_skymodel)
        
        expec_skymodel_chunks.append(expec_chunk)
        
    for chunk_ind in range(num_gauss_chunks):
        expec_chunk = populate_pointgauss_chunk(CompTypes.GAUSSIAN, chunk_ind,
                            comps_per_chunk, n_powers,
                            n_curves, n_lists, num_list_values,
                            num_of_each_comp, above_horizon,
                            expec_ra, expec_dec,
                            expec_pow_fluxes, expec_cur_fluxes,
                            fits_skymodel=fits_skymodel)
        
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
        
        expec_chunk = populate_shapelet_chunk(expec_chunk,
                            num_list_values, num_coeff_per_shape,
                            above_horizon, expec_ra, expec_dec,
                            comp_inds, orig_comp_inds[comp_inds],
                            chunk_basis_param_indexes,
                            chunk_basis_values,
                            fits_skymodel=fits_skymodel)
        
        expec_skymodel_chunks.append(expec_chunk)
    
    return expec_skymodel_chunks


class Skymodel_Settings:
    """Something to hold all the various settings and pass around between
    functions"""
    def __init__(self, deg_between_comps : float,
                 num_coeff_per_shape : int,
                 num_list_values : int,
                 comps_per_source : int,
                 stokesV_frac_cadence : int = 0,
                 stokesV_pl_cadence : int = 0,
                 stokesV_cpl_cadence : int = 0,
                 linpol_frac_cadence : int = 0,
                 linpol_pl_cadence : int = 0,
                 linpol_cpl_cadence : int = 0):
        
        self.deg_between_comps = deg_between_comps
        self.num_coeff_per_shape = num_coeff_per_shape
        self.num_list_values = num_list_values
        self.comps_per_source = comps_per_source
        self.stokesV_frac_cadence = stokesV_frac_cadence
        self.stokesV_pl_cadence = stokesV_pl_cadence
        self.stokesV_cpl_cadence = stokesV_cpl_cadence
        self.linpol_frac_cadence = linpol_frac_cadence
        self.linpol_pl_cadence = linpol_pl_cadence
        self.linpol_cpl_cadence = linpol_cpl_cadence
        
        
def make_expected_comp_counter(settings : Skymodel_Settings):
    """This doesn't 100% create everything that `read_fits_radec_count_components`
    would fill, but it makes all the arrays that other values are derived from.
    So if these are good, I would just be copying code from the function
    which I think is a waste o time"""
    
    ra_range = np.arange(0, 360.0*D2R, settings.deg_between_comps*D2R)
    dec_range = np.arange(LOW_DEC, HIGH_DEC, settings.deg_between_comps*D2R)
    ra_range, dec_range = np.meshgrid(ra_range, dec_range)
    ra_range, dec_range = ra_range.flatten(), dec_range.flatten()
    
    num_radec = len(ra_range)
    num_comp_types = 3 #point, gauss, shape
    total_num_comps = num_radec*NUM_FLUX_TYPES*num_comp_types
    
    expec_comp_counter = Component_Type_Counter(initial_size=total_num_comps)
    
    ##The index of each component to parent source, set via settings.comps_per_source
    expec_comp_counter.source_indexes = np.arange(total_num_comps) // settings.comps_per_source
    
    ##This is the order things are written when making example sky model
    expec_comp_counter.comp_types = np.empty(total_num_comps)
    
    expec_comp_counter.point_power_inds = []
    expec_comp_counter.point_curve_inds = []
    expec_comp_counter.point_list_inds = []
    expec_comp_counter.gauss_power_inds = []
    expec_comp_counter.gauss_curve_inds = []
    expec_comp_counter.gauss_list_inds = []
    expec_comp_counter.shape_power_inds = []
    expec_comp_counter.shape_curve_inds = []
    expec_comp_counter.shape_list_inds = []
    
    for coord in range(num_radec):
        base = coord*NUM_FLUX_TYPES*num_comp_types
        expec_comp_counter.comp_types[base + 0] = CompTypes.POINT_POWER.value
        expec_comp_counter.point_power_inds.append(base + 0)
        expec_comp_counter.comp_types[base + 1] = CompTypes.POINT_CURVE.value
        expec_comp_counter.point_curve_inds.append(base + 1)
        expec_comp_counter.comp_types[base + 2] = CompTypes.POINT_LIST.value
        expec_comp_counter.point_list_inds.append(base + 2)
        expec_comp_counter.comp_types[base + 3] = CompTypes.GAUSS_POWER.value
        expec_comp_counter.gauss_power_inds.append(base + 3)
        expec_comp_counter.comp_types[base + 4] = CompTypes.GAUSS_CURVE.value
        expec_comp_counter.gauss_curve_inds.append(base + 4)
        expec_comp_counter.comp_types[base + 5] = CompTypes.GAUSS_LIST.value
        expec_comp_counter.gauss_list_inds.append(base + 5)
        expec_comp_counter.comp_types[base + 6] = CompTypes.SHAPE_POWER.value
        expec_comp_counter.shape_power_inds.append(base + 6)
        expec_comp_counter.comp_types[base + 7] = CompTypes.SHAPE_CURVE.value
        expec_comp_counter.shape_curve_inds.append(base + 7)
        expec_comp_counter.comp_types[base + 8] = CompTypes.SHAPE_LIST.value
        expec_comp_counter.shape_list_inds.append(base + 8)
        
    ##We stick a set number of flux entries for the list types, and they
    ##appear every third component
    expec_comp_counter.num_list_fluxes = np.zeros(total_num_comps)
    expec_comp_counter.num_list_fluxes[np.arange(2, total_num_comps, 3)] = settings.num_list_values
    
    
    ##We stick a set number of shapelet coeffs for the shapelet type
    expec_comp_counter.num_shape_coeffs = np.zeros(total_num_comps)
    
    all_shapes = np.where((expec_comp_counter.comp_types == CompTypes.SHAPE_POWER.value) |
                          (expec_comp_counter.comp_types == CompTypes.SHAPE_CURVE.value) |
                          (expec_comp_counter.comp_types == CompTypes.SHAPE_LIST.value))[0]
    
    expec_comp_counter.num_shape_coeffs[all_shapes] = settings.num_coeff_per_shape
    
    stokesV_cpl_cadence = settings.stokesV_cpl_cadence
    stokesV_pl_cadence = settings.stokesV_pl_cadence
    stokesV_frac_cadence = settings.stokesV_frac_cadence
    
    expec_comp_counter.orig_v_point_power_inds = []
    expec_comp_counter.orig_v_point_curve_inds = []
    expec_comp_counter.orig_v_point_pol_frac_inds = []
    expec_comp_counter.orig_v_gauss_power_inds = []
    expec_comp_counter.orig_v_gauss_curve_inds = []
    expec_comp_counter.orig_v_gauss_pol_frac_inds = []
    expec_comp_counter.orig_v_shape_power_inds = []
    expec_comp_counter.orig_v_shape_curve_inds = []
    expec_comp_counter.orig_v_shape_pol_frac_inds = []
    
    ##Copy the logic from fits_skymodel_common.add_stokesV_fits to work this out
    ##Also need to factor in that point, gauss, shape have a PPPGGGSSS order
    
    if stokesV_cpl_cadence or stokesV_pl_cadence or stokesV_frac_cadence:
        expec_comp_counter.v_comp_types = np.full(total_num_comps, np.nan, dtype=np.float64)
        
        for comp_index in range(total_num_comps):
        
            while True:
            
                if stokesV_frac_cadence:
                    if comp_index % stokesV_frac_cadence == 0:
                        if comp_index % 9 < 3:
                            expec_comp_counter.v_comp_types[comp_index] = CompTypes.V_POINT_POL_FRAC.value
                            expec_comp_counter.orig_v_point_pol_frac_inds.append(comp_index)
                        elif comp_index % 9 < 6:
                            expec_comp_counter.v_comp_types[comp_index] = CompTypes.V_GAUSS_POL_FRAC.value
                            expec_comp_counter.orig_v_gauss_pol_frac_inds.append(comp_index)
                        else:
                            expec_comp_counter.v_comp_types[comp_index] = CompTypes.V_SHAPE_POL_FRAC.value
                            expec_comp_counter.orig_v_shape_pol_frac_inds.append(comp_index)
                        break
                        
                if stokesV_pl_cadence:
                    if comp_index % stokesV_pl_cadence == 0:
                        if comp_index % 9 < 3:
                            expec_comp_counter.v_comp_types[comp_index] = CompTypes.V_POINT_POWER.value
                            expec_comp_counter.orig_v_point_power_inds.append(comp_index)
                        elif comp_index % 9 < 6:
                            expec_comp_counter.v_comp_types[comp_index] = CompTypes.V_GAUSS_POWER.value
                            expec_comp_counter.orig_v_gauss_power_inds.append(comp_index)
                        else:
                            expec_comp_counter.v_comp_types[comp_index] = CompTypes.V_SHAPE_POWER.value
                            expec_comp_counter.orig_v_shape_power_inds.append(comp_index)
                        break
                        
                if stokesV_cpl_cadence:
                    if comp_index % stokesV_cpl_cadence == 0:
                        if comp_index % 9 < 3:
                            expec_comp_counter.v_comp_types[comp_index] = CompTypes.V_POINT_CURVE.value
                            expec_comp_counter.orig_v_point_curve_inds.append(comp_index)
                        elif comp_index % 9 < 6:
                            expec_comp_counter.v_comp_types[comp_index] = CompTypes.V_GAUSS_CURVE.value
                            expec_comp_counter.orig_v_gauss_curve_inds.append(comp_index)
                        else:
                            expec_comp_counter.v_comp_types[comp_index] = CompTypes.V_SHAPE_CURVE.value
                            expec_comp_counter.orig_v_shape_curve_inds.append(comp_index)
                        break
                break
            
            
    linpol_cpl_cadence = settings.linpol_cpl_cadence
    linpol_pl_cadence = settings.linpol_pl_cadence
    linpol_frac_cadence = settings.linpol_frac_cadence
    
    expec_comp_counter.orig_lin_point_power_inds = []
    expec_comp_counter.orig_lin_point_curve_inds = []
    expec_comp_counter.orig_lin_point_pol_frac_inds = []
    expec_comp_counter.orig_lin_gauss_power_inds = []
    expec_comp_counter.orig_lin_gauss_curve_inds = []
    expec_comp_counter.orig_lin_gauss_pol_frac_inds = []
    expec_comp_counter.orig_lin_shape_power_inds = []
    expec_comp_counter.orig_lin_shape_curve_inds = []
    expec_comp_counter.orig_lin_shape_pol_frac_inds = []
    
    ##Copy the logic from fits_skymodel_common.add_linpol_fits to work this out
    ##Also need to factor in that point, gauss, shape have a PPPGGGSSS order
    
    if linpol_cpl_cadence or linpol_pl_cadence or linpol_frac_cadence:
        expec_comp_counter.lin_comp_types = np.full(total_num_comps, np.nan, dtype=np.float64)
        
        for comp_index in range(total_num_comps):
        
            while True:
            
                if linpol_frac_cadence:
                    if comp_index % linpol_frac_cadence == 0:
                        if comp_index % 9 < 3:
                            expec_comp_counter.lin_comp_types[comp_index] = CompTypes.LIN_POINT_POL_FRAC.value
                            expec_comp_counter.orig_lin_point_pol_frac_inds.append(comp_index)
                        elif comp_index % 9 < 6:
                            expec_comp_counter.lin_comp_types[comp_index] = CompTypes.LIN_GAUSS_POL_FRAC.value
                            expec_comp_counter.orig_lin_gauss_pol_frac_inds.append(comp_index)
                        else:
                            expec_comp_counter.lin_comp_types[comp_index] = CompTypes.LIN_SHAPE_POL_FRAC.value
                            expec_comp_counter.orig_lin_shape_pol_frac_inds.append(comp_index)
                        break
                        
                if linpol_pl_cadence:
                    if comp_index % linpol_pl_cadence == 0:
                        if comp_index % 9 < 3:
                            expec_comp_counter.lin_comp_types[comp_index] = CompTypes.LIN_POINT_POWER.value
                            expec_comp_counter.orig_lin_point_power_inds.append(comp_index)
                        elif comp_index % 9 < 6:
                            expec_comp_counter.lin_comp_types[comp_index] = CompTypes.LIN_GAUSS_POWER.value
                            expec_comp_counter.orig_lin_gauss_power_inds.append(comp_index)
                        else:
                            expec_comp_counter.lin_comp_types[comp_index] = CompTypes.LIN_SHAPE_POWER.value
                            expec_comp_counter.orig_lin_shape_power_inds.append(comp_index)
                        break
                        
                if linpol_cpl_cadence:
                    if comp_index % linpol_cpl_cadence == 0:
                        if comp_index % 9 < 3:
                            expec_comp_counter.lin_comp_types[comp_index] = CompTypes.LIN_POINT_CURVE.value
                            expec_comp_counter.orig_lin_point_curve_inds.append(comp_index)
                        elif comp_index % 9 < 6:
                            expec_comp_counter.lin_comp_types[comp_index] = CompTypes.LIN_GAUSS_CURVE.value
                            expec_comp_counter.orig_lin_gauss_curve_inds.append(comp_index)
                        else:
                            expec_comp_counter.lin_comp_types[comp_index] = CompTypes.LIN_SHAPE_CURVE.value
                            expec_comp_counter.orig_lin_shape_curve_inds.append(comp_index)
                        break
                break

    return expec_comp_counter


def check_comp_counter(comp_counter : Component_Type_Counter,
                       settings : Skymodel_Settings):
    """Check the component counter matches expectations"""
    
    ##Make what we should have found for testing
    expec_comp_counter = make_expected_comp_counter(settings)
    
    npt.assert_array_equal(expec_comp_counter.source_indexes,
                           comp_counter.source_indexes)
    npt.assert_array_equal(expec_comp_counter.comp_types,
                           comp_counter.comp_types)
    npt.assert_array_equal(expec_comp_counter.num_list_fluxes,
                           comp_counter.num_list_fluxes)
    npt.assert_array_equal(expec_comp_counter.num_shape_coeffs,
                           comp_counter.num_shape_coeffs)
    
    npt.assert_array_equal(expec_comp_counter.point_power_inds,
                           comp_counter.point_power_inds)
    npt.assert_array_equal(expec_comp_counter.point_curve_inds,
                           comp_counter.point_curve_inds)
    npt.assert_array_equal(expec_comp_counter.point_list_inds,
                           comp_counter.point_list_inds)
    npt.assert_array_equal(expec_comp_counter.gauss_power_inds,
                           comp_counter.gauss_power_inds)
    npt.assert_array_equal(expec_comp_counter.gauss_curve_inds,
                           comp_counter.gauss_curve_inds)
    npt.assert_array_equal(expec_comp_counter.gauss_list_inds,
                           comp_counter.gauss_list_inds)
    npt.assert_array_equal(expec_comp_counter.shape_power_inds,
                           comp_counter.shape_power_inds)
    npt.assert_array_equal(expec_comp_counter.shape_curve_inds,
                           comp_counter.shape_curve_inds)
    npt.assert_array_equal(expec_comp_counter.shape_list_inds,
                           comp_counter.shape_list_inds)
    
    if settings.stokesV_cpl_cadence or settings.stokesV_pl_cadence or settings.stokesV_frac_cadence:
        npt.assert_array_equal(expec_comp_counter.v_comp_types,
                               comp_counter.v_comp_types)
        
        npt.assert_array_equal(expec_comp_counter.orig_v_point_power_inds,
                               comp_counter.orig_v_point_power_inds)
        npt.assert_array_equal(expec_comp_counter.orig_v_point_curve_inds,
                               comp_counter.orig_v_point_curve_inds)
        npt.assert_array_equal(expec_comp_counter.orig_v_point_pol_frac_inds,
                               comp_counter.orig_v_point_pol_frac_inds)
        npt.assert_array_equal(expec_comp_counter.orig_v_gauss_power_inds,
                               comp_counter.orig_v_gauss_power_inds)
        npt.assert_array_equal(expec_comp_counter.orig_v_gauss_curve_inds,
                               comp_counter.orig_v_gauss_curve_inds)
        npt.assert_array_equal(expec_comp_counter.orig_v_gauss_pol_frac_inds,
                               comp_counter.orig_v_gauss_pol_frac_inds)
        npt.assert_array_equal(expec_comp_counter.orig_v_shape_power_inds,
                               comp_counter.orig_v_shape_power_inds)
        npt.assert_array_equal(expec_comp_counter.orig_v_shape_curve_inds,
                               comp_counter.orig_v_shape_curve_inds)
        npt.assert_array_equal(expec_comp_counter.orig_v_shape_pol_frac_inds,
                               comp_counter.orig_v_shape_pol_frac_inds)
        
    if settings.linpol_cpl_cadence or settings.linpol_pl_cadence or settings.linpol_frac_cadence:
        npt.assert_array_equal(expec_comp_counter.lin_comp_types,
                               comp_counter.lin_comp_types)
        
        npt.assert_array_equal(expec_comp_counter.orig_lin_point_power_inds,
                               comp_counter.orig_lin_point_power_inds)
        npt.assert_array_equal(expec_comp_counter.orig_lin_point_curve_inds,
                               comp_counter.orig_lin_point_curve_inds)
        npt.assert_array_equal(expec_comp_counter.orig_lin_point_pol_frac_inds,
                               comp_counter.orig_lin_point_pol_frac_inds)
        npt.assert_array_equal(expec_comp_counter.orig_lin_gauss_power_inds,
                               comp_counter.orig_lin_gauss_power_inds)
        npt.assert_array_equal(expec_comp_counter.orig_lin_gauss_curve_inds,
                               comp_counter.orig_lin_gauss_curve_inds)
        npt.assert_array_equal(expec_comp_counter.orig_lin_gauss_pol_frac_inds,
                               comp_counter.orig_lin_gauss_pol_frac_inds)
        npt.assert_array_equal(expec_comp_counter.orig_lin_shape_power_inds,
                               comp_counter.orig_lin_shape_power_inds)
        npt.assert_array_equal(expec_comp_counter.orig_lin_shape_curve_inds,
                               comp_counter.orig_lin_shape_curve_inds)
        npt.assert_array_equal(expec_comp_counter.orig_lin_shape_pol_frac_inds,
                               comp_counter.orig_lin_shape_pol_frac_inds)