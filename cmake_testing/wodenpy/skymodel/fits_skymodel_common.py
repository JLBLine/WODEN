from sys import path
import os
import unittest
import numpy as np
import erfa
from astropy.table import Column, Table
from astropy.io import fits
from wodenpy.skymodel.woden_skymodel import Component_Type_Counter, Component_Info, CompTypes, calc_pl_norm_at_200MHz, calc_cpl_norm_at_200MHz

from read_skymodel_common import Skymodel_Settings

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

def add_power_law_fits(table_dict, row_ind, value):
    """add power law info for a component"""
    
    table_dict['mod_type'].append('pl')
    table_dict['norm_comp_pl'][row_ind] = float(value)
    table_dict['alpha_pl'][row_ind] = float(value)
    
def add_curved_power_law_fits(table_dict, row_ind, value):
    """add curved power law info for a component"""
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

def add_stokesV_fits(table_dict, comp_index, stokesV_frac_cadence,
                     stokesV_pl_cadence, stokesV_cpl_cadence):
    """Add a Stokes V polarisation fraction type to a component"""
    
    ##OK we only want to add the Stokes V polarisation info for one specific
    ##flux model, and we're doing things on a cadence where we might get overlaps,
    ##so do a logic test to work out who to add.
    
    ##Important for many choices is that any combination of the different types
    ##may be in play; if they are zero, don't use them at all
    
    while True:
        
        if stokesV_frac_cadence:
            if comp_index % stokesV_frac_cadence == 0:
                table_dict['V_MOD_TYPE'][comp_index] = 'pf'
                table_dict['V_POL_FRAC'][comp_index] = float(comp_index)
                break
                
        if stokesV_pl_cadence:
            if comp_index % stokesV_pl_cadence == 0:
                table_dict['V_MOD_TYPE'][comp_index] = 'pl'
                table_dict['V_NORM_COMP_PL'][comp_index] = comp_index
                table_dict['V_ALPHA_PL'][comp_index] = comp_index
                break
                
        if stokesV_cpl_cadence:
            if comp_index % stokesV_cpl_cadence == 0:
                table_dict['V_MOD_TYPE'][comp_index] = 'cpl'
                table_dict['V_NORM_COMP_CPL'][comp_index] = comp_index
                table_dict['V_ALPHA_CPL'][comp_index] = comp_index
                table_dict['V_CURVE_CPL'][comp_index] = comp_index
                break
        break
    
def add_linpol_fits(table_dict, comp_index, linpol_frac_cadence,
                     linpol_pl_cadence, linpol_cpl_cadence):
    """Add a Stokes V polarisation fraction type to a component"""
    
    ##OK we only want to add the Stokes V polarisation info for one specific
    ##flux model, and we're doing things on a cadence where we might get overlaps,
    ##so do a logic test to work out who to add.
    
    ##Important for many choices is that any combination of the different types
    ##may be in play; if they are zero, don't use them at all
    
    while True:
        
        if linpol_frac_cadence:
            if comp_index % linpol_frac_cadence == 0:
                table_dict['LIN_MOD_TYPE'][comp_index] = 'pf'
                # table_dict['LIN_POL_FRAC'][comp_index] = 0.01*float(comp_index/linpol_frac_cadence)
                table_dict['LIN_POL_FRAC'][comp_index] = float(comp_index)
                table_dict['RM'][comp_index] = float(comp_index*((2*np.pi)/360))
                table_dict['INTR_POL_ANGLE'][comp_index] = 0.1*float(comp_index*((2*np.pi)/360))
                
                break
                
        if linpol_pl_cadence:
            if comp_index % linpol_pl_cadence == 0:
                table_dict['LIN_MOD_TYPE'][comp_index] = 'pl'
                table_dict['LIN_NORM_COMP_PL'][comp_index] = comp_index
                table_dict['LIN_ALPHA_PL'][comp_index] = comp_index
                table_dict['RM'][comp_index] = float(comp_index*((2*np.pi)/360))
                table_dict['INTR_POL_ANGLE'][comp_index] = 0.1*float(comp_index*((2*np.pi)/360))
                break
                
        if linpol_cpl_cadence:
            if comp_index % linpol_cpl_cadence == 0:
                table_dict['LIN_MOD_TYPE'][comp_index] = 'cpl'
                table_dict['LIN_NORM_COMP_CPL'][comp_index] = comp_index
                table_dict['LIN_ALPHA_CPL'][comp_index] = comp_index
                table_dict['LIN_CURVE_CPL'][comp_index] = comp_index
                table_dict['RM'][comp_index] = float(comp_index*((2*np.pi)/360))
                table_dict['INTR_POL_ANGLE'][comp_index] = 0.1*float(comp_index*((2*np.pi)/360))
                break
        break
                
                
def add_point_fits(table_dict : dict, all_flux_cols : list,
              ra : float, dec: float,
              point_index : int, gauss_index : int, shape_index : int,
              comp_type : CompTypes,
              flux_index : int, source_index : int,
              settings : Skymodel_Settings):
    
    comp_index = point_index + gauss_index + shape_index
    if comp_index % settings.comps_per_source == 0:
        source_index += 1
    
    table_dict['unq_source_id'].append(f"{source_index:07d}")
    table_dict['name'].append(f"{source_index:07d}_C{comp_index:05d}")
    table_dict['comp_type'].append('P')
    
    row_ind = comp_index
    table_dict['ra'][row_ind] = ra
    table_dict['dec'][row_ind] = dec
    
    if comp_type == CompTypes.POINT_POWER:
        add_power_law_fits(table_dict, row_ind, point_index)
        
    elif comp_type == CompTypes.POINT_CURVE:
        add_curved_power_law_fits(table_dict, row_ind, point_index)
        
    elif comp_type == CompTypes.POINT_LIST:
        flux_index = add_list_flux_fits(table_dict, all_flux_cols, row_ind, flux_index, settings.num_list_values)
        
    ##This will add the Stokes V polarisation fraction style if required
    add_stokesV_fits(table_dict, comp_index, settings.stokesV_frac_cadence,
                     settings.stokesV_pl_cadence, settings.stokesV_cpl_cadence)
    
    add_linpol_fits(table_dict, comp_index, settings.linpol_frac_cadence,
                     settings.linpol_pl_cadence, settings.linpol_cpl_cadence)
        
        
    point_index += 1
        
    return source_index, point_index, flux_index

def add_gauss_fits(table_dict : dict, all_flux_cols : list,
                    ra : float, dec: float,
                    point_index : int, gauss_index : int, shape_index : int,
                    comp_type : CompTypes,
                    flux_index : int, source_index : int,
                    settings : Skymodel_Settings):
    
    comp_index = point_index + gauss_index + shape_index
    if comp_index % settings.comps_per_source == 0:
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
        flux_index = add_list_flux_fits(table_dict, all_flux_cols, row_ind, flux_index, settings.num_list_values)
        
    ##This will add the Stokes V polarisation fraction style if required
    add_stokesV_fits(table_dict, comp_index, settings.stokesV_frac_cadence,
                     settings.stokesV_pl_cadence, settings.stokesV_cpl_cadence)
    
    add_linpol_fits(table_dict, comp_index, settings.linpol_frac_cadence,
                     settings.linpol_pl_cadence, settings.linpol_cpl_cadence)

    gauss_index += 1
        
    return source_index, gauss_index, flux_index
    

def add_shapelet_fits(table_dict : dict, all_flux_cols : list,
                    ra : float, dec: float,
                    point_index : int, gauss_index : int, shape_index : int,
                    comp_type : CompTypes,
                    flux_index : int, source_index : int, basis_index : int,
                    settings : Skymodel_Settings):
    
    comp_index = point_index + gauss_index + shape_index
    if comp_index % settings.comps_per_source == 0:
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
    
    for basis in range(basis_index, basis_index + settings.num_coeff_per_shape):
    
        table_dict['s_names'].append(f"{source_index:07d}_C{comp_index:05d}")
        table_dict['s_n1s'][basis] = float(basis)
        table_dict['s_n2s'][basis] = float(basis)
        table_dict['s_coeffs'][basis] = float(basis)
    
    basis_index += settings.num_coeff_per_shape
    
    if comp_type == CompTypes.SHAPE_POWER:
        add_power_law_fits(table_dict, row_ind, shape_index)
        
    elif comp_type == CompTypes.SHAPE_CURVE:
        add_curved_power_law_fits(table_dict, row_ind, shape_index)
        
    elif comp_type == CompTypes.SHAPE_LIST:
        flux_index = add_list_flux_fits(table_dict, all_flux_cols, row_ind, flux_index, settings.num_list_values)
        
    ##This will add the Stokes V polarisation fraction style if required
    add_stokesV_fits(table_dict, comp_index, settings.stokesV_frac_cadence,
                     settings.stokesV_pl_cadence, settings.stokesV_cpl_cadence)
    
    add_linpol_fits(table_dict, comp_index, settings.linpol_frac_cadence,
                     settings.linpol_pl_cadence, settings.linpol_cpl_cadence)
        
    shape_index += 1
        
    return source_index, shape_index, flux_index, basis_index


# def write_full_test_skymodel_fits(deg_between_comps : float,
#                              num_coeff_per_shape : int,
#                              num_list_values : int,
#                              comps_per_source : int,
#                              stokesV_frac_cadence : int = 0,
#                              stokesV_pl_cadence : int = 0,
#                              stokesV_cpl_cadence : int = 0):
def write_full_test_skymodel_fits(settings : Skymodel_Settings):
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
    
    ra_range = np.arange(0, 360.0*D2R, settings.deg_between_comps*D2R)
    dec_range = np.arange(LOW_DEC, HIGH_DEC, settings.deg_between_comps*D2R)
    
    
    
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
    
    # point = False
    # gaussian = False
    # shapelet = False
    
    point = True
    gaussian = True
    shapelet = True
    
    
    num_comp_types = point + gaussian + shapelet
    num_list_freqs = num_comp_types * len(ra_range) * settings.num_list_values
    
    total_num_comps = NUM_FLUX_TYPES*num_comp_types*len(ra_range)
    
    ##Anything that is a string column is easier to make as a list and then
    ##convert to column
    unq_source_id = []
    name = []
    
    ra_col = Column(data=np.full(total_num_comps, np.nan) ,name="RA")
    dec_col = Column(data=np.full(total_num_comps, np.nan), name="DEC")
    major_dc = Column(data=np.full(total_num_comps, np.nan), name="MAJOR_DC")
    minor_dc = Column(data=np.full(total_num_comps, np.nan), name="MINOR_DC")
    pa_dc = Column(data=np.full(total_num_comps, np.nan), name="PA_DC")
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
    
    ##Now we add the optional circular polarisation columns
    if settings.stokesV_frac_cadence or settings.stokesV_pl_cadence or settings.stokesV_cpl_cadence:
        vpol_type = Column(data=np.full(total_num_comps, '', dtype='|S3'), name="V_MOD_TYPE", dtype='|S3')
        table_dict['V_MOD_TYPE'] = vpol_type
    
    if settings.stokesV_frac_cadence:
        vpol_frac = Column(data=np.full(total_num_comps, np.nan), name="V_POL_FRAC")
        table_dict['V_POL_FRAC'] = vpol_frac
        
    if settings.stokesV_pl_cadence:
        v_norm_comp_pl = Column(data=np.full(total_num_comps, np.nan), name="V_NORM_COMP_PL")
        table_dict['V_NORM_COMP_PL'] = v_norm_comp_pl
        v_alpha_pl = Column(data=np.full(total_num_comps, np.nan), name="V_ALPHA_PL")
        table_dict['V_ALPHA_PL'] = v_alpha_pl
        
    if settings.stokesV_cpl_cadence:
        v_norm_comp_cpl = Column(data=np.full(total_num_comps, np.nan), name="V_NORM_COMP_CPL")
        table_dict['V_NORM_COMP_CPL'] = v_norm_comp_cpl
        v_alpha_cpl = Column(data=np.full(total_num_comps, np.nan), name="V_ALPHA_CPL")
        table_dict['V_ALPHA_CPL'] = v_alpha_cpl
        v_curve_cpl = Column(data=np.full(total_num_comps, np.nan), name="V_CURVE_CPL")
        table_dict['V_CURVE_CPL'] = v_curve_cpl
    
    ##Now we add the optional linear polarisation columns    
    if settings.linpol_frac_cadence or settings.linpol_pl_cadence or settings.linpol_cpl_cadence:
        linpol_type = Column(data=np.full(total_num_comps, '', dtype='|S3'), name="LIN_MOD_TYPE", dtype='|S3')
        table_dict['LIN_MOD_TYPE'] = linpol_type
        rm = Column(data=np.full(total_num_comps, np.nan), name='RM')
        table_dict['RM'] = rm
        intr_pol_angle = Column(data=np.full(total_num_comps, np.nan), name='INTR_POL_ANGLE')
        table_dict['INTR_POL_ANGLE'] = intr_pol_angle
    
    if settings.linpol_frac_cadence:
        linpol_frac = Column(data=np.full(total_num_comps, np.nan), name="LIN_POL_FRAC")
        table_dict['LIN_POL_FRAC'] = linpol_frac
        
    if settings.linpol_pl_cadence:
        lin_norm_comp_pl = Column(data=np.full(total_num_comps, np.nan), name="LIN_NORM_COMP_PL")
        table_dict['LIN_NORM_COMP_PL'] = lin_norm_comp_pl
        lin_alpha_pl = Column(data=np.full(total_num_comps, np.nan), name="LIN_ALPHA_PL")
        table_dict['LIN_ALPHA_PL'] = lin_alpha_pl
        
    if settings.linpol_cpl_cadence:
        lin_norm_comp_cpl = Column(data=np.full(total_num_comps, np.nan), name="LIN_NORM_COMP_CPL")
        table_dict['LIN_NORM_COMP_CPL'] = lin_norm_comp_cpl
        lin_alpha_cpl = Column(data=np.full(total_num_comps, np.nan), name="LIN_ALPHA_CPL")
        table_dict['LIN_ALPHA_CPL'] = lin_alpha_cpl
        lin_curve_cpl = Column(data=np.full(total_num_comps, np.nan), name="LIN_CURVE_CPL")
        table_dict['LIN_CURVE_CPL'] = lin_curve_cpl
    
    ##Let's populate the source and component names
    
    ##OK in the FITS format, every different frequency has a different column
    ##For the test I setup, where every list entry has a different flux, this
    ##is a little mental, but easier to use the testing logic made for the yaml
    ##So just cop it
    all_flux_cols = [Column(data=np.full(total_num_comps, np.nan), name=f"INT_FLX{flux_col_ind:07.3f}") for flux_col_ind in range(settings.num_list_values)]
    
    num_shape_coeffs = shapelet*settings.num_coeff_per_shape*NUM_FLUX_TYPES*len(ra_range)
    
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
                                            flux_index, source_index, settings)
            
            source_index, point_index, flux_index = add_point_fits(table_dict, all_flux_cols,
                                            ra, dec,
                                            point_index, gauss_index, shape_index,
                                            CompTypes.POINT_CURVE,
                                            flux_index, source_index, settings)
            
            source_index, point_index, flux_index = add_point_fits(table_dict, all_flux_cols,
                                            ra, dec,
                                            point_index, gauss_index, shape_index,
                                            CompTypes.POINT_LIST,
                                            flux_index, source_index, settings)
            
        if gaussian:
            source_index, gauss_index, flux_index = add_gauss_fits(table_dict, all_flux_cols,
                                            ra, dec,
                                            point_index, gauss_index, shape_index,
                                            CompTypes.GAUSS_POWER,
                                            flux_index, source_index, settings)
            
            source_index, gauss_index, flux_index = add_gauss_fits(table_dict, all_flux_cols,
                                            ra, dec,
                                            point_index, gauss_index, shape_index,
                                            CompTypes.GAUSS_CURVE,
                                            flux_index, source_index, settings)
            
            source_index, gauss_index, flux_index = add_gauss_fits(table_dict, all_flux_cols,
                                            ra, dec,
                                            point_index, gauss_index, shape_index,
                                            CompTypes.GAUSS_LIST,
                                            flux_index, source_index, settings)
        
        if shapelet: 
            source_index, shape_index, flux_index, basis_index = add_shapelet_fits(table_dict, all_flux_cols,
                                            ra, dec,
                                            point_index, gauss_index, shape_index,
                                            CompTypes.SHAPE_POWER,
                                            flux_index, source_index,
                                            basis_index, settings)
            
            source_index, shape_index, flux_index, basis_index = add_shapelet_fits(table_dict, all_flux_cols,
                                            ra, dec,
                                            point_index, gauss_index, shape_index,
                                            CompTypes.SHAPE_CURVE,
                                            flux_index, source_index,
                                            basis_index, settings)
            
            source_index, shape_index, flux_index, basis_index = add_shapelet_fits(table_dict, all_flux_cols,
                                            ra, dec,
                                            point_index, gauss_index, shape_index,
                                            CompTypes.SHAPE_LIST,
                                            flux_index, source_index,
                                            basis_index, settings)
        
    unq_source_id = Column(data=unq_source_id, name="UNQ_SOURCE_ID")
    name = Column(data=name, name="NAME")
    mod_type = Column(data=mod_type, name="MOD_TYPE")
    comp_type = Column(data=comp_type, name="COMP_TYPE")
    s_names = Column(data=s_names, name="NAME")
    
    out_columns = [unq_source_id, name, ra_col, dec_col, major_dc, minor_dc, pa_dc, mod_type, comp_type,
                   norm_comp_pl, alpha_pl, norm_comp_cpl, alpha_cpl, curve_cpl]
    
    ##Add circ pol columns if required
    if settings.stokesV_frac_cadence or settings.stokesV_pl_cadence or settings.stokesV_cpl_cadence:
        out_columns.append(vpol_type)
    
    if settings.stokesV_frac_cadence:
        out_columns.append(vpol_frac)
        
    if settings.stokesV_pl_cadence:
        out_columns.append(v_norm_comp_pl)
        out_columns.append(v_alpha_pl)
        
    if settings.stokesV_cpl_cadence:
        out_columns.append(v_norm_comp_cpl)
        out_columns.append(v_alpha_cpl)
        out_columns.append(v_curve_cpl)
        
    ##Add circ pol columns if required
    if settings.linpol_frac_cadence or settings.linpol_pl_cadence or settings.linpol_cpl_cadence:
        out_columns.append(linpol_type)
        out_columns.append(rm)
        out_columns.append(intr_pol_angle)
        
    if settings.linpol_frac_cadence:
        out_columns.append(linpol_frac)
        
    if settings.linpol_pl_cadence:
        out_columns.append(lin_norm_comp_pl)
        out_columns.append(lin_alpha_pl)
        
    if settings.linpol_cpl_cadence:
        out_columns.append(lin_norm_comp_cpl)
        out_columns.append(lin_alpha_cpl)
        out_columns.append(lin_curve_cpl)

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