import numpy as np
import sys
import os
from typing import Union

from wodenpy.skymodel.woden_skymodel import Component_Type_Counter, Component_Info, CompTypes
from wodenpy.skymodel.chunk_sky_model import Skymodel_Chunk_Map, Components_Map
from wodenpy.use_libwoden.skymodel_structs import setup_chunked_source, setup_source_catalogue, _Ctype_Source_Into_Python
from wodenpy.use_libwoden.beam_settings import BeamTypes
from wodenpy.use_libwoden.create_woden_struct_classes import Woden_Struct_Classes
from astropy.table import Table, Column
import erfa
from astropy.io import fits
from wodenpy.primary_beam.everybeam import run_everybeam, load_OSKAR_telescope, get_everybeam_norm
from sys import exit
from wodenpy.use_libwoden.create_woden_struct_classes import Woden_Struct_Classes
import argparse
from astropy.time import Time, TimeDelta

##This call is so we can use it as a type annotation
woden_struct_classes = Woden_Struct_Classes()
Woden_Settings = woden_struct_classes.Woden_Settings


REF_FREQ = 200e+6

D2R = np.pi/180.0

def _error_for_missing_columns(message_start : str, fits_path : str, missing_cols : list,
                                col_present_dict : dict):
    """Given a list of missing columns to test for, error out with a helpful message
    detailing which columns are missing from the FITS file at `fits_path`.

    Parameters
    ----------
    message_start : str
        Initial part of the error message.
    fits_path : str
        Path to FITS sky model file.
    missing_cols : list
        List of missing columns.
    col_present_dict : dict
        Dictionary of all columns in the FITS file, and whether they are present or not.
    """
    
    present = 0
    
    for col in missing_cols:
        present += col_present_dict[col]
        
    if present < len(missing_cols):
        missing = ""
        for col in missing_cols:
            if not col_present_dict[col]:
                missing += f"{col} "
        
        sys.exit(f"Found {message_start} in the FITS file {fits_path}, but missing {missing}columns. Exiting now.")
        
def _error_for_missing_int_flx(table, col_entry, col_prepend, mod_col_name, fits_path, hdu_name=False):
    """
    Check if the table has any columns starting with `col_prepend` and raise an error if not found.
    Parameters
    ----------
    table: 
        The table to check for columns.
    col_entry: 
        The name of the SED model type column entry (e.g. "nan")
    col_prepend: 
        The prefix of the columns to search for (e.g. "INT_FLX")
    mod_col_name: 
        The name of the model column in the FITS file.
    fits_path: 
        The path to the FITS file.
    hdu_name: 
        The name of the HDU in the FITS file (default is False).
    Raises:
    - SystemExit: If no columns starting with `col_prepend` are found in the table.
                 If `hdu_name` is provided, the error message includes the HDU name.
                 If no `hdu_name` is provided, the error message does not include the HDU name.
                 If `col_prepend` is not 'INT_FLX' and 'NAME' column is not found in the table.
    """
    
    
    have_any_int_flx = False
    for key in table.columns:
        if key[:len(col_prepend)] == col_prepend:
            have_any_int_flx = True
    
    if not have_any_int_flx:
        
        if hdu_name:
            sys.exit(f"Found {col_entry} (list flux) in {mod_col_name} column in the FITS file {fits_path}, but no columns start with {col_prepend} in the {hdu_name} HDU, which are how list-type fluxes are added. Exiting now.")
        else:
            sys.exit(f"Found {col_entry} (list flux) in {mod_col_name} column in the FITS file {fits_path}, but no columns start with {col_prepend}, which are how list-type fluxes are added. Exiting now.")
    
    if col_prepend != 'INT_FLX':
    
        if 'NAME' not in table.columns:
            sys.exit(f"Found {col_entry} (list flux) in {mod_col_name} column in the FITS file {fits_path}, but no NAME column in the {hdu_name} HDU. The NAME links the flux entries to the ra/dec in the main table so is necessary. Exiting now.")
                
def _error_for_missing_flx_table_cols(hdu_names, hdu_name,  col_entry, col_prepend, mod_col_name, fits_path):
    """
    Check if the specified HDU name exists in the given list of HDU names.
    If the HDU name exists, read the FITS file using the specified HDU name and perform additional checks.
    If the HDU name does not exist, exit the program with an error message.
    Parameters:
    - hdu_names (list): A list of HDU names.
    - hdu_name (str): The name of the HDU to check for existence.
    - col_entry (str): The name of the column entry.
    - col_prepend (str): The prepend value for the column.
    - mod_col_name (str): The name of the modified column.
    - fits_path (str): The path to the FITS file.
    Returns:
    None
    """
    
    
    if hdu_name in hdu_names:
        list_table = Table.read(fits_path, hdu=hdu_name )
        _error_for_missing_int_flx(list_table, col_entry, col_prepend, mod_col_name, fits_path, hdu_name)
    else:
        exit(f"Found {col_entry} (list flux) in {mod_col_name} column in the FITS file {fits_path}, but no HDU named {hdu_name} was present. You need an HDU containing the fluxes, and to name that HDU {hdu_name}. Exiting now.")

def check_columns_fits(fits_path : str):
    """Checks that the FITS file at `fits_path` contains all the necessary
    information to be read in correctly. Also checks for optional columns, and
    certain combinations that must be present e.g. if there are shapelets
    specified in `COMP_TYPE`, then there must be a second shapelet HDU.
    
    Any problems and the function errors out with a helpful error message.

    Parameters
    ----------
    fits_path : str
        Path to FITS sky model file.
    """
    
    with fits.open(fits_path) as hdus:
        num_hdus = len(hdus)
        hdu_names = [hdu.name for hdu in hdus]
    
    main_table = Table.read(fits_path, hdu=1)
    
    all_col_names = ["UNQ_SOURCE_ID", "NAME", "RA", "DEC", "COMP_TYPE",
                     "MAJOR_DC", "MINOR_DC", "PA_DC", "MOD_TYPE",
                     "NORM_COMP_PL", "ALPHA_PL", "NORM_COMP_CPL", "ALPHA_CPL",
                     "CURVE_CPL", "V_MOD_TYPE", "V_POL_FRAC", "V_NORM_COMP_PL",
                     "V_ALPHA_PL", "V_NORM_COMP_CPL", "V_ALPHA_CPL",
                     "V_CURVE_CPL", "LIN_MOD_TYPE", "RM", "INTR_POL_ANGLE",
                     "LIN_POL_FRAC", "LIN_NORM_COMP_PL", "LIN_ALPHA_PL",
                     "LIN_NORM_COMP_CPL", "LIN_ALPHA_CPL", "LIN_CURVE_CPL"]
    
    ##set up a dictionary, set all columns as unfound
    col_present_dict = {}
    for col_name in all_col_names:
        col_present_dict[col_name] = False
    
    ##loop through all the columns in the FITS file, and if they are present,
    ##mark them as found
    for col_name in all_col_names:
        if col_name in main_table.columns:
            col_present_dict[col_name] = True
            
    ##Every catalogue must contain these columns, so error if missing
    essentials = ["UNQ_SOURCE_ID", "NAME", "RA", "DEC", "COMP_TYPE", "MOD_TYPE"]
    for essential in essentials:
        if col_present_dict[essential] == False:
            sys.exit(f"Essential column {essential} is missing from the FITS file {fits_path}. Exiting now.")
            
    ##Check for Gaussians, and verify that the major, minor, PA columns are present
    if 'G' in main_table['COMP_TYPE']:
        _error_for_missing_columns("Gaussian components", fits_path,
                                   ["MAJOR_DC", "MINOR_DC", "PA_DC"],
                                   col_present_dict)
        
    ##Check for shapelets, and verify that the major, minor, PA columns are present
    ##Also check the second shapelet table is present, and contains columns it needs
    if 'S' in main_table['COMP_TYPE']:
        _error_for_missing_columns("Shapelet components", fits_path,
                                   ["MAJOR_DC", "MINOR_DC", "PA_DC"],
                                   col_present_dict)
            
        ##Minimum number of HDUs is 2, as header is first, main table second
        if num_hdus < 3:
            sys.exit(f"Found Shapelet components in the FITS file {fits_path}, but missing the shapelet table. Exiting now.")
            
        ##check if there are names on the HDUs. If not, assume shapelets is in the second table HDU
        ##we didn't name HDUs before we started with the whole polarisation gamut
        if 'SHAPELET' in hdu_names:
            shape_table = Table.read(fits_path, hdu='SHAPELET')
        else:
            shape_table = Table.read(fits_path, hdu=2)
            
        shape_table_cols = ["NAME", "N1", "N2", "COEFF"]
        for shape_col in shape_table_cols:
            if shape_col not in shape_table.columns:
                sys.exit(f"Found Shapelet components in the FITS file {fits_path}, but missing the {shape_col} column in the shapelet table. Exiting now.")
                
        if np.max(shape_table['N1']) > 100 or np.max(shape_table['N2']) > 100:
            sys.exit(f"Found Shapelet components in the FITS file {fits_path}, but the maximum N1 or N2 value is greater than 100. Basis functions are only stored up to 100, so this isn't possible. Exiting now.")
                
    
    ##Check for Stokes I flux models making sense-------------------------------
    if 'pl' in main_table['MOD_TYPE']:
        _error_for_missing_columns("pl (power-law) in MOD_TYPE column", fits_path,
                                   ["NORM_COMP_PL", "ALPHA_PL"],
                                   col_present_dict)
        
    if 'cpl' in main_table['MOD_TYPE']:
        _error_for_missing_columns("cpl (curved power-law) in MOD_TYPE column", fits_path,
                                   ["NORM_COMP_CPL", "ALPHA_CPL", "CURVE_CPL"],
                                   col_present_dict)
    
    ##This one is trickier because the name starts with 'INT_FLX' but can
    ##end with any number (which is the frequency in MHz)
    if 'nan' in main_table['MOD_TYPE']:
        _error_for_missing_int_flx(main_table, 'nan', 'INT_FLX', 'MOD_TYPE', fits_path)
                
    
    ##Check for Stokes V information, and verify that V_MOD_TYPE is present-------
    if col_present_dict["V_POL_FRAC"] or col_present_dict["V_NORM_COMP_PL"] or col_present_dict["V_ALPHA_PL"] or col_present_dict["V_NORM_COMP_CPL"] or col_present_dict["V_ALPHA_CPL"] or col_present_dict["V_CURVE_CPL"]:
        
        if not col_present_dict["V_MOD_TYPE"]:
            sys.exit(f"Found polarisation information in the FITS file {fits_path}, but no V_MOD_TYPE column. Needed to distinguish between Stokes V model types. Exiting now.")
            
    if col_present_dict["V_MOD_TYPE"]:
            
        ##Check for Stokes V flux models making sense-------------------------------
        if 'pl' in main_table['V_MOD_TYPE']:
            _error_for_missing_columns("pl (power-law) in V_MOD_TYPE column", fits_path,
                                    ["V_NORM_COMP_PL", "V_ALPHA_PL"],
                                    col_present_dict)
            
        if 'cpl' in main_table['V_MOD_TYPE']:
            _error_for_missing_columns("cpl (curved power-law) in V_MOD_TYPE column", fits_path,
                                    ["V_NORM_COMP_CPL", "V_ALPHA_CPL", "V_CURVE_CPL"],
                                    col_present_dict)
            
        if 'pf' in main_table['V_MOD_TYPE']:
            _error_for_missing_columns("pf (polarisation fraction) in V_MOD_TYPE column", fits_path,
                                    ["V_POL_FRAC"],
                                    col_present_dict)
        
        if 'nan' in main_table['V_MOD_TYPE']:
            _error_for_missing_flx_table_cols(hdu_names, 'V_LIST_FLUXES',  'nan',
                                              'V_INT_FLX', 'V_MOD_TYPE', fits_path)
            
    ##Check for linear polarisation information, and verify that LIN_MOD_TYPE is present-------
    if col_present_dict["LIN_POL_FRAC"] or col_present_dict["LIN_NORM_COMP_PL"] or col_present_dict["LIN_ALPHA_PL"] or col_present_dict["LIN_NORM_COMP_CPL"] or col_present_dict["LIN_ALPHA_CPL"] or col_present_dict["LIN_CURVE_CPL"] or col_present_dict["RM"] or col_present_dict["INTR_POL_ANGLE"]:
        
        if not col_present_dict["LIN_MOD_TYPE"]:
            sys.exit(f"Found polarisation information in the FITS file {fits_path}, but no LIN_MOD_TYPE column. Needed to distinguish between linear polarisation model types. Exiting now.")
            
    if col_present_dict["LIN_MOD_TYPE"]:
            
        ##Check for Stokes V flux models making sense-------------------------------
        if 'pl' in main_table['LIN_MOD_TYPE']:
            _error_for_missing_columns("pl (power-law) in LIN_MOD_TYPE column", fits_path,
                                    ["LIN_NORM_COMP_PL", "LIN_ALPHA_PL", "RM"],
                                    col_present_dict)
            
        if 'cpl' in main_table['LIN_MOD_TYPE']:
            _error_for_missing_columns("cpl (curved power-law) in LIN_MOD_TYPE column", fits_path,
                                    ["LIN_NORM_COMP_CPL", "LIN_ALPHA_CPL", "LIN_CURVE_CPL", "RM"],
                                    col_present_dict)
            
        if 'pf' in main_table['LIN_MOD_TYPE']:
            _error_for_missing_columns("pf (polarisation fraction) in LIN_MOD_TYPE column", fits_path,
                                    ["LIN_POL_FRAC", "RM"],
                                    col_present_dict)
        
        ##'nan' means list types for both Q and U. So test that both the
        ##HDUs exist
        if 'nan' in main_table['LIN_MOD_TYPE']:
            _error_for_missing_flx_table_cols(hdu_names, 'Q_LIST_FLUXES',  'nan',
                                              'Q_INT_FLX', 'LIN_MOD_TYPE', fits_path)
            _error_for_missing_flx_table_cols(hdu_names, 'U_LIST_FLUXES',  'nan',
                                              'U_INT_FLX', 'LIN_MOD_TYPE', fits_path)
        
        ##'p_nan' means list type for linear polarised flux. So test that the polarised
        ##HDU exists. And we need at least rotation measure, as we'll be doing
        ##the whole RM rotation thing
        if 'p_nan' in main_table['LIN_MOD_TYPE']:
            _error_for_missing_flx_table_cols(hdu_names, 'P_LIST_FLUXES',  'p_nan',
                                              'P_INT_FLX', 'LIN_MOD_TYPE', fits_path)
            _error_for_missing_columns("p_nan (list flux) in LIN_MOD_TYPE column",
                                       fits_path, ["RM"], col_present_dict)
    

def count_num_list_fluxes(flux_mod_col_name : str, flux_col_prepend : str,
                          num_list_fluxes : np.ndarray,
                          main_table : Table, comp_counter : Component_Type_Counter,
                          model_type='nan',
                          optional_table : Table = False):
    """Counts the number of list-type fluxes for a given model type. Default
    is to look for Stokes I fluxes in the main table, but can be one of the 
    optional HDUs like `V_LIST_FLUXES`; point `optional_table` to correct
    table if needed. Fills `num_list_fluxes` with the number of list-type fluxes
    for each component of the given model type.
    
    Also checks that the number of references to that model type (e.g. 'nan'
    in the `MOD_TYPE` column) matches the number of rows in the table that
    have list fluxes. If they don't match, the function exits.

    Parameters
    ----------
    flux_mod_col_name : str
        The model column name, e.g. 'MOD_TYPE' or `V_MOD_TYPE`.
    flux_col_prepend : str
        Prepend that flux columns start with, e.g. 'INT_FLX' or 'V_INT_FLX'.
    num_list_fluxes : np.ndarray
        Array to fill with the number of list-type fluxes for each component (e.g. `comp_counter.num_list_fluxes`).
    main_table : Table
        The main table (e.g. `Table.read(fits_path, hdu=1)`)
    comp_counter : Component_Type_Counter
        _description_
    model_type : str, optional
        _description_, by default 'nan'
    optional_table : Table, optional
        _description_, by default False
    """
    
    flux_col_names = []
    list_type_inds = np.where(main_table[flux_mod_col_name] == model_type)[0]
            
    if optional_table:
        ##Do an extra check here that the number of model types in the main table
        ##match 
        use_table = optional_table
        num_flux_lists = len(use_table['NAME'])
        
        if num_flux_lists != len(list_type_inds):
            exit(f"The number of {model_type} entries in main table column {flux_mod_col_name} does not match the number of rows in the table {[flux_col_prepend[0]]}_LIST_FLUXES. They should match otherwise things will go wrong. Exiting now.")
    else:
        use_table = main_table
        
    for key in use_table.columns:
        if key[:len(flux_col_prepend)] == flux_col_prepend:
            flux_col_names.append(key)
            
    # print('HALLO', flux_col_names)
        
    for flux_col in flux_col_names:
        
        if type(use_table[flux_col]) == Column:
            present_fluxes = np.where(np.isnan(use_table[flux_col]) == False)[0]
            num_list_fluxes[list_type_inds] += 1
        ##OK, astropy.Table can turns things into masked arrays. Here we are
        ##accessing a column `main_table[flux_col]`, reading it's `mask` which
        ##returns a bunch of booleans (True == masked). We want everything that
        ##isn't masked, so use the ~ to swap. Finally, we want ints, not booleans
        else:
            present_fluxes = (~use_table[flux_col].mask).astype(int)
            num_list_fluxes[list_type_inds] += present_fluxes[list_type_inds]

# @profile
def read_fits_radec_count_components(fits_path : str):
    """
    Reads a FITS file sky model. Reads just the ra, dec, and counts how many POINT/GAUSS/SHAPE and POWER/CURVE/LIST, consolidating the
    information into a Component_Type_Counter object.
    
    Parameters
    -----------
    fits_path : str
        The path to the FITS file containing the sky model information.
    
    Returns
    --------
    comp_counter : Component_Type_Counter
        A Component_Type_Counter object that counts the number of components of each type and their properties.
    """
    if not os.path.isfile(fits_path):
        sys.exit(f"Cannot read sky model from {fits_path}. Please check your paths, exiting now.")
        
        
    ##grab all the relevant information out of the tables
    
    main_table = Table.read(fits_path, hdu=1)
    
    ras = np.array(main_table['RA'], dtype=np.float64)
    decs = np.array(main_table['DEC'], dtype=np.float64)
    comp_types = np.array(main_table['COMP_TYPE'], dtype=str)
    flux_types = np.array(main_table['MOD_TYPE'], dtype=str)
    unq_source_ID = np.array(main_table['UNQ_SOURCE_ID'], dtype=str)
    comp_names = np.array(main_table['NAME'], dtype=str)
    
    num_comps = len(ras)
    
    comp_counter = Component_Type_Counter(initial_size=num_comps)
    
    comp_counter.comp_ras = ras*D2R
    comp_counter.comp_decs = decs*D2R
    
    ##Right, we want an index of every component to a parent source, which
    ##is named via `source_names`.
    ##we can get this using np.unique
    
    source_names, comp_source_indexes = np.unique(unq_source_ID, return_inverse=True)
    
    comp_counter.source_indexes = comp_source_indexes
    
    ##point component stuff
    point_power = np.where((comp_types == 'P') & (flux_types == 'pl'))
    point_curve = np.where((comp_types == 'P') & (flux_types == 'cpl'))
    point_list = np.where((comp_types == 'P') & (flux_types == 'nan'))
    
    comp_counter.comp_types[point_power] = CompTypes.POINT_POWER.value
    comp_counter.comp_types[point_curve] = CompTypes.POINT_CURVE.value
    comp_counter.comp_types[point_list] = CompTypes.POINT_LIST.value
    
    ##gaussian component stuff
    gauss_power = np.where((comp_types == 'G') & (flux_types == 'pl'))
    gauss_curve = np.where((comp_types == 'G') & (flux_types == 'cpl'))
    gauss_list = np.where((comp_types == 'G') & (flux_types == 'nan'))
    
    comp_counter.comp_types[gauss_power] = CompTypes.GAUSS_POWER.value
    comp_counter.comp_types[gauss_curve] = CompTypes.GAUSS_CURVE.value
    comp_counter.comp_types[gauss_list] = CompTypes.GAUSS_LIST.value
    
    ##shape component stuff
    shape_power = np.where((comp_types == 'S') & (flux_types == 'pl'))
    shape_curve = np.where((comp_types == 'S') & (flux_types == 'cpl'))
    shape_list = np.where((comp_types == 'S') & (flux_types == 'nan'))
    
    comp_counter.comp_types[shape_power] = CompTypes.SHAPE_POWER.value
    comp_counter.comp_types[shape_curve] = CompTypes.SHAPE_CURVE.value
    comp_counter.comp_types[shape_list] = CompTypes.SHAPE_LIST.value
    
    ##for the list flux types, we want to count how many freqs we have.
    ##need them to malloc the ctype sky model later on
    count_num_list_fluxes('MOD_TYPE', 'INT_FLX', comp_counter.num_list_fluxes,
                    main_table, comp_counter)
    
    ##Count how many hdus there are; first table should always be the second hdu?
    with fits.open(fits_path) as hdus:
        num_hdus = len(hdus)
    
    if num_hdus > 2:
        shape_table = Table.read(fits_path, hdu=2)
        shape_comp_names = np.array(shape_table['NAME'], dtype=str)
        
        ##now count up how many shapelet basis functions each component has
        for comp_ind in np.where(comp_types == 'S')[0]:
            comp_name = comp_names[comp_ind]
            
            basis_func_inds = np.where(shape_comp_names == comp_name)[0]
            comp_counter.num_shape_coeffs[comp_ind] = len(basis_func_inds)
    else:
        print("INFO: couldn't find second table containing shapelet information, so not attempting to load any shapelets.")
        
    ##OK OK so the following indexing seems like overkill, but it's gathering
    ##indexes that we will use in the lazy-loading of the FITS file later on
    ##Given we chunk by point, gauss, shapelet, need to keep things separated
    ##by component type
    ##Doing it once here means we won't be doing it multiple times later on
    ##Check for (optional) polarisation information and store if needed
    if 'V_MOD_TYPE' in main_table.columns:
        ##Get column, find index of different types
        v_mod_types = np.array(main_table['V_MOD_TYPE'], dtype=str)
        
        v_point_power = np.where((comp_types == 'P') & (v_mod_types == 'pl'))
        v_point_curve = np.where((comp_types == 'P') & (v_mod_types == 'cpl'))
        v_point_pol_frac = np.where((comp_types == 'P') & (v_mod_types == 'pf'))
        v_point_list = np.where((comp_types == 'P') & (v_mod_types == 'nan'))
        v_gauss_power = np.where((comp_types == 'G') & (v_mod_types == 'pl'))
        v_gauss_curve = np.where((comp_types == 'G') & (v_mod_types == 'cpl'))
        v_gauss_pol_frac = np.where((comp_types == 'G') & (v_mod_types == 'pf'))
        v_gauss_list = np.where((comp_types == 'G') & (v_mod_types == 'nan'))
        v_shape_power = np.where((comp_types == 'S') & (v_mod_types == 'pl'))
        v_shape_curve = np.where((comp_types == 'S') & (v_mod_types == 'cpl'))
        v_shape_pol_frac = np.where((comp_types == 'S') & (v_mod_types == 'pf'))
        v_shape_list = np.where((comp_types == 'S') & (v_mod_types == 'nan'))
        
        ##Stick em in the comp counter with specific values
        comp_counter.v_comp_types = np.full(num_comps, np.nan, dtype=np.float64)
        comp_counter.v_comp_types[v_point_power] = CompTypes.V_POINT_POWER.value
        comp_counter.v_comp_types[v_point_curve] = CompTypes.V_POINT_CURVE.value
        comp_counter.v_comp_types[v_point_pol_frac] = CompTypes.V_POINT_POL_FRAC.value
        comp_counter.v_comp_types[v_point_list] = CompTypes.V_POINT_LIST.value
        comp_counter.v_comp_types[v_gauss_power] = CompTypes.V_GAUSS_POWER.value
        comp_counter.v_comp_types[v_gauss_curve] = CompTypes.V_GAUSS_CURVE.value
        comp_counter.v_comp_types[v_gauss_pol_frac] = CompTypes.V_GAUSS_POL_FRAC.value
        comp_counter.v_comp_types[v_gauss_list] = CompTypes.V_GAUSS_LIST.value
        comp_counter.v_comp_types[v_shape_power] = CompTypes.V_SHAPE_POWER.value
        comp_counter.v_comp_types[v_shape_curve] = CompTypes.V_SHAPE_CURVE.value
        comp_counter.v_comp_types[v_shape_pol_frac] = CompTypes.V_SHAPE_POL_FRAC.value
        comp_counter.v_comp_types[v_shape_list] = CompTypes.V_SHAPE_LIST.value
        
    ##Check for (optional) polarisation information and store if needed
    if 'LIN_MOD_TYPE' in main_table.columns:
        ##Get column, find index of different types
        lin_mod_types = np.array(main_table['LIN_MOD_TYPE'], dtype=str)
        lin_point_power = np.where((comp_types == 'P') & (lin_mod_types == 'pl'))
        lin_point_curve = np.where((comp_types == 'P') & (lin_mod_types == 'cpl'))
        lin_point_pol_frac = np.where((comp_types == 'P') & (lin_mod_types == 'pf'))
        lin_point_list = np.where((comp_types == 'P') & (lin_mod_types == 'nan'))
        lin_point_p_list = np.where((comp_types == 'P') & (lin_mod_types == 'p_nan'))
        lin_gauss_power = np.where((comp_types == 'G') & (lin_mod_types == 'pl'))
        lin_gauss_curve = np.where((comp_types == 'G') & (lin_mod_types == 'cpl'))
        lin_gauss_pol_frac = np.where((comp_types == 'G') & (lin_mod_types == 'pf'))
        lin_gauss_list = np.where((comp_types == 'G') & (lin_mod_types == 'nan'))
        lin_gauss_p_list = np.where((comp_types == 'G') & (lin_mod_types == 'p_nan'))
        lin_shape_power = np.where((comp_types == 'S') & (lin_mod_types == 'pl'))
        lin_shape_curve = np.where((comp_types == 'S') & (lin_mod_types == 'cpl'))
        lin_shape_pol_frac = np.where((comp_types == 'S') & (lin_mod_types == 'pf'))
        lin_shape_list = np.where((comp_types == 'S') & (lin_mod_types == 'nan'))
        lin_shape_p_list = np.where((comp_types == 'S') & (lin_mod_types == 'p_nan'))
        
        ##Stick em in the comp counter with specific values
        comp_counter.lin_comp_types = np.full(num_comps, np.nan, dtype=np.float64)
        comp_counter.lin_comp_types[lin_point_power] = CompTypes.LIN_POINT_POWER.value
        comp_counter.lin_comp_types[lin_point_curve] = CompTypes.LIN_POINT_CURVE.value
        comp_counter.lin_comp_types[lin_point_list] = CompTypes.LIN_POINT_LIST.value
        comp_counter.lin_comp_types[lin_point_p_list] = CompTypes.LIN_POINT_P_LIST.value
        comp_counter.lin_comp_types[lin_point_pol_frac] = CompTypes.LIN_POINT_POL_FRAC.value
        comp_counter.lin_comp_types[lin_gauss_power] = CompTypes.LIN_GAUSS_POWER.value
        comp_counter.lin_comp_types[lin_gauss_curve] = CompTypes.LIN_GAUSS_CURVE.value
        comp_counter.lin_comp_types[lin_gauss_list] = CompTypes.LIN_GAUSS_LIST.value
        comp_counter.lin_comp_types[lin_gauss_p_list] = CompTypes.LIN_GAUSS_P_LIST.value
        comp_counter.lin_comp_types[lin_gauss_pol_frac] = CompTypes.LIN_GAUSS_POL_FRAC.value
        comp_counter.lin_comp_types[lin_shape_power] = CompTypes.LIN_SHAPE_POWER.value
        comp_counter.lin_comp_types[lin_shape_curve] = CompTypes.LIN_SHAPE_CURVE.value
        comp_counter.lin_comp_types[lin_shape_list] = CompTypes.LIN_SHAPE_LIST.value
        comp_counter.lin_comp_types[lin_shape_p_list] = CompTypes.LIN_SHAPE_P_LIST.value
        comp_counter.lin_comp_types[lin_shape_pol_frac] = CompTypes.LIN_SHAPE_POL_FRAC.value
        
        ##check to see if they've included a INTR_POL_ANGLE angle
        if 'INTR_POL_ANGLE' in main_table.columns:
            comp_counter.has_intr_pol_angle = True
    
    comp_counter.total_components()
    ##If we have list-style polarisation information, we need to count how many
    ##flux entries there are for each component (needed for mallocing down the line)
    
    if comp_counter.num_v_point_lists + comp_counter.num_v_gauss_lists + comp_counter.num_v_shape_lists > 0:
        comp_counter.num_v_list_fluxes = np.zeros(num_comps, dtype=np.int32)
    
        count_num_list_fluxes('V_MOD_TYPE', 'V_INT_FLX', comp_counter.num_v_list_fluxes,
                            main_table, comp_counter,
                            optional_table=Table.read(fits_path, hdu='V_LIST_FLUXES'))
        
    if comp_counter.num_lin_point_lists + comp_counter.num_lin_gauss_lists + comp_counter.num_lin_shape_lists > 0:
        comp_counter.num_q_list_fluxes = np.zeros(num_comps, dtype=np.int32)
        comp_counter.num_u_list_fluxes = np.zeros(num_comps, dtype=np.int32)
    
        count_num_list_fluxes('LIN_MOD_TYPE', 'Q_INT_FLX', comp_counter.num_q_list_fluxes,
                            main_table, comp_counter,
                            optional_table=Table.read(fits_path, hdu='Q_LIST_FLUXES'))
        count_num_list_fluxes('LIN_MOD_TYPE', 'U_INT_FLX', comp_counter.num_u_list_fluxes,
                            main_table, comp_counter,
                            optional_table=Table.read(fits_path, hdu='U_LIST_FLUXES'))
        
    if comp_counter.num_lin_point_p_lists + comp_counter.num_lin_gauss_p_lists + comp_counter.num_lin_shape_p_lists > 0:
        comp_counter.num_p_list_fluxes = np.zeros(num_comps, dtype=np.int32)
    
        count_num_list_fluxes('LIN_MOD_TYPE', 'P_INT_FLX', comp_counter.num_p_list_fluxes,
                            main_table, comp_counter,
                            optional_table=Table.read(fits_path, hdu='P_LIST_FLUXES'),
                            model_type='p_nan')
        
    return comp_counter


def find_all_indexes_of_x_in_y(x : np.ndarray, y : np.ndarray):
    """
    Given two arrays `x` and `y`, find the indexes of matching elements from `y`
    in `x`.
    
    See here https://stackoverflow.com/questions/8251541/numpy-for-every-element-in-one-array-find-the-index-in-another-array
    for the origins of this function (github copypasta)
        
    Parameters
    ----------
    x : np.ndarray
        Array of elements to find in `y`.
    y : np.ndarray
        Array to search for elements in `x`.

    Returns
    -------
    np.ndarray
        Array of indexes of `y` that match the elements in `x`.
    """
    xsorted = np.argsort(y)
    ypos = np.searchsorted(y[xsorted], x)
    
    return xsorted[ypos]

woden_struct_classes = Woden_Struct_Classes()
Source = woden_struct_classes.Source
Components = woden_struct_classes.Components

def add_list_flux_info(mod_type : CompTypes, n_lists : int,
                       key_prepend : str, list_table : Table, map_components : Components_Map,
                       source_components : Components, main_table : Table = False,
                       comp_orig_inds : np.ndarray = False,
                       list_comp_ind_offset : int = 0):
    
    if mod_type == CompTypes.I_LIST:
        list_comp_inds = source_components.list_comp_inds
        num_list_values = source_components.num_list_values
        list_start_indexes = source_components.list_start_indexes
        list_freqs = source_components.list_freqs
        list_ref_flux = source_components.list_stokesI
        list_orig_inds = map_components.list_orig_inds
        new_list_indexes = np.arange(n_lists)
        
    elif mod_type == CompTypes.V_LIST:
        list_comp_inds = source_components.stokesV_list_comp_inds
        num_list_values = source_components.stokesV_num_list_values
        list_start_indexes = source_components.stokesV_list_start_indexes
        list_freqs = source_components.stokesV_list_ref_freqs
        list_ref_flux = source_components.stokesV_list_ref_flux
        main_orig_inds = map_components.v_list_orig_inds
        
    elif mod_type == CompTypes.Q_LIST:
        list_comp_inds = source_components.stokesQ_list_comp_inds
        num_list_values = source_components.stokesQ_num_list_values
        list_start_indexes = source_components.stokesQ_list_start_indexes
        list_freqs = source_components.stokesQ_list_ref_freqs
        list_ref_flux = source_components.stokesQ_list_ref_flux
        main_orig_inds = map_components.lin_list_orig_inds
        
    elif mod_type == CompTypes.U_LIST:
        list_comp_inds = source_components.stokesU_list_comp_inds
        num_list_values = source_components.stokesU_num_list_values
        list_start_indexes = source_components.stokesU_list_start_indexes
        list_freqs = source_components.stokesU_list_ref_freqs
        list_ref_flux = source_components.stokesU_list_ref_flux
        main_orig_inds = map_components.lin_list_orig_inds
        
    elif mod_type == CompTypes.LIN_LIST:
        list_comp_inds = source_components.linpol_p_list_comp_inds
        num_list_values = source_components.linpol_p_num_list_values
        list_start_indexes = source_components.linpol_p_list_start_indexes
        list_freqs = source_components.linpol_p_list_ref_freqs
        list_ref_flux = source_components.linpol_p_list_ref_flux
        main_orig_inds = map_components.lin_p_list_orig_inds
        
        # print(list_comp_inds)
        # print(num_list_values)
        # print(list_start_indexes)
        # print(list_freqs)
        # print(list_ref_flux)
        # print("main_orig_inds",main_orig_inds)
    
    ##If we're not doing Stokes I, we're using list information from a different
    ##table. We have the original indexes w.r.t the main_table, but the 
    ##list information is stored in the list table. So we need to match the 
    ##list information, which is done via the NAME column
    if mod_type != CompTypes.I_LIST:
        
        list_name_subset = np.array(main_table['NAME'])[main_orig_inds]
        all_list_names = np.array(list_table['NAME'], dtype=str)
        list_orig_inds = find_all_indexes_of_x_in_y(list_name_subset, all_list_names)
        ##also have to find the index of these list-type components w.r.t the
        ##chunk subset we're going into
        new_list_indexes = find_all_indexes_of_x_in_y(main_orig_inds, comp_orig_inds)
        
    ##Get all the flux column names/freqs and put em in a list
    flux_col_freqs = []
    flux_col_names = []
    for key in list_table.columns:
        if key[:len(key_prepend)] == key_prepend:
            flux_col_names.append(key)
            flux_col_freqs.append(float(key[len(key_prepend):])*1e+6)
    flux_col_freqs = np.array(flux_col_freqs)
    
    ##Read all flux info for the current chunk into an array so we can do
    ##faster array manipulation on it
    
    all_list_fluxes = np.zeros((n_lists, len(flux_col_names)), dtype=np.float64)
    
    for flux_ind, flux_col in enumerate(flux_col_names):
        all_list_fluxes[:, flux_ind] = list_table[flux_col][list_orig_inds]
        
    # ##TODO vectorise this please, is probs sloooow
    list_start_index = 0
    for list_ind, new_list_ind in enumerate(new_list_indexes):
        
        list_comp_inds[list_ind] = list_comp_ind_offset + new_list_ind
        
        ##Gather all frequency column information for this component
        these_fluxes = all_list_fluxes[list_ind, :]
        
        ##Figure out what isn't a nan to grab actual values
        use_fluxes = np.where(np.isnan(these_fluxes) == False)
        
        these_fluxes = these_fluxes[use_fluxes]
        these_freqs = flux_col_freqs[use_fluxes]
        
        num_list_values[list_ind] = int(len(these_freqs))
        list_start_indexes[list_ind] = list_start_index
        
        find = 0
        for freq, flux in zip(these_freqs, these_fluxes):
            list_freqs[list_start_index + find] = freq
            list_ref_flux[list_start_index + find] = flux
            find += 1
        
        list_start_index += int(len(these_freqs))

def add_fits_info_to_source_catalogue(comp_type : CompTypes,
                        main_table : Table, shape_table : Table,
                        chunk_source : Source,
                        chunk_map : Skymodel_Chunk_Map,
                        beamtype : int, lsts : np.ndarray, latitude : float,
                        v_table : Table = False, q_table : Table = False,
                        u_table : Table = False, p_table : Table = False, 
                        ra0 = False, dec0 = False, telescope = False,
                        beam_norms = False,
                        all_times : np.ndarray = False,
                        all_freqs : np.ndarray = False):
    """Given the desired components as detailed in the `chunk_map`, add
    the relevant information from the FITS file `main_table`, `shape_table` objects to the `chunk_source` object. As well as the skymodel information, this function adds either
    az/za or ha/dec information, depending on the `beamtype`.

    Parameters
    ----------
    comp_type : CompTypes
        The type of component we are adding information for.
    main_table : Table
        The main Table (with RA,Dec etc) from the FITS file.
    shape_table : Table
        The shapelet Table from the FITS file.
    chunk_source : Source
        The ctypes Source object to add information to.
    chunk_map : Skymodel_Chunk_Map
        The map object containing information about components for this chunk.
    beamtype : int
        The type of beam used (BeamTypes)
    lsts : np.ndarray
        LSTs for all time steps for this simulation.
    latitude : float
        Latitude of the array
    """
    
    num_time_steps = len(lsts)
    
    if comp_type == CompTypes.POINT:
        
        n_powers = chunk_map.n_point_powers
        n_curves = chunk_map.n_point_curves
        n_lists = chunk_map.n_point_lists
        
        power_inds = chunk_map.point_components.power_orig_inds
        curve_inds = chunk_map.point_components.curve_orig_inds
        list_inds = chunk_map.point_components.list_orig_inds
        
        source_components = chunk_source.point_components
        map_components = chunk_map.point_components
        
    elif comp_type == CompTypes.GAUSSIAN:
        
        n_powers = chunk_map.n_gauss_powers
        n_curves = chunk_map.n_gauss_curves
        n_lists = chunk_map.n_gauss_lists
        
        power_inds = chunk_map.gauss_components.power_orig_inds
        curve_inds = chunk_map.gauss_components.curve_orig_inds
        list_inds = chunk_map.gauss_components.list_orig_inds
        
        source_components = chunk_source.gauss_components
        map_components = chunk_map.gauss_components
        
    elif comp_type == CompTypes.SHAPELET:
        n_powers = chunk_map.n_shape_powers
        n_curves = chunk_map.n_shape_curves
        n_lists = chunk_map.n_shape_lists
        
        power_inds = chunk_map.shape_components.power_orig_inds
        curve_inds = chunk_map.shape_components.curve_orig_inds
        list_inds = chunk_map.shape_components.list_orig_inds
        
        source_components = chunk_source.shape_components
        map_components = chunk_map.shape_components
        
    num_components = n_powers + n_curves + n_lists
        
    ##chunk_map.all_orig_inds contains indexes of all comp type, i.e.
    ##possibly POINT and GAUSSIAN, so find all indexes for this component
    ##type to iterate through, in order of power,curve,list flux type
    comp_orig_inds = np.empty(n_powers + n_curves + n_lists, dtype=int)
    comp_orig_inds[:n_powers] = power_inds
    comp_orig_inds[n_powers:n_powers+n_curves] = curve_inds
    comp_orig_inds[n_powers+n_curves:] = list_inds
        
    ##TODO get some kind of numpy -> ctype conversion, rather than iterate?
    ##fill in positional information - all compoment types / flux models need this
    ##this should grab all RA/Dec for all components regardless of flux model type
    for new_comp_ind, old_comp_ind in enumerate(comp_orig_inds):
        source_components.ras[new_comp_ind] = main_table['RA'][old_comp_ind]*D2R
        source_components.decs[new_comp_ind] = main_table['DEC'][old_comp_ind]*D2R
        
        # print('ra, dec', main_table['RA'][old_comp_ind], main_table['DEC'][old_comp_ind])
        # print('ra, dec', source_components.ras[new_comp_ind], source_components.decs[new_comp_ind])
        
        if beamtype == BeamTypes.GAUSS_BEAM.value or beamtype == BeamTypes.MWA_ANALY.value:
                comp_has = lsts - source_components.ras[new_comp_ind]
                
                ##OK these ctype arrays cannot be sliced, so let's increment
                ##over them at a snail's pace
                for time_ind in range(num_time_steps):
                    hadec_low = new_comp_ind*num_time_steps
                    source_components.beam_has[hadec_low + time_ind] = comp_has[time_ind]
                    source_components.beam_decs[hadec_low + time_ind] = source_components.decs[new_comp_ind]
        
        ##Some beam models don't need az/za coords
        skip_beams = [BeamTypes.NO_BEAM.value, BeamTypes.GAUSS_BEAM.value,
                     BeamTypes.EB_OSKAR.value]
        ##only the NO_BEAM and GAUSS_BEAM options don't need az,za coords
        if beamtype in skip_beams:
            pass
        else:
            ##Calculate lst, and then azimuth/elevation
            comp_has = lsts - source_components.ras[new_comp_ind]
            comp_azs, comp_els = erfa.hd2ae(comp_has, source_components.decs[new_comp_ind],
                                            latitude)
            
            ##OK these ctype arrays cannot be sliced, so let's increment
            ##over them at a snail's pace
            for time_ind in range(num_time_steps):
                azza_low = new_comp_ind*num_time_steps
                source_components.azs[azza_low + time_ind] = comp_azs[time_ind]
                source_components.zas[azza_low + time_ind] = np.pi/2 - comp_els[time_ind]
                
        eb_beams = [BeamTypes.EB_OSKAR.value]
        
        if beamtype in eb_beams:
            num_freqs = len(all_freqs)
            ##Everybeam models are calculated on the CPU. This function will only
        ##run if we have set args.primary_beam to an everybeam model
        
            for time_ind, time in enumerate(all_times):
                for freq_ind, freq in enumerate(all_freqs):
                    
                    jones = run_everybeam(source_components.ras[new_comp_ind],
                                          source_components.decs[new_comp_ind],
                                          ra0, dec0, time, freq, telescope,
                                          beam_norms=beam_norms[time_ind, freq_ind, :])
                    
                    beam_ind = num_freqs*time_ind*num_components + (num_components*freq_ind) + new_comp_ind
                    
                    ##TODO likely need some kind of reordering to do IAU
                    ##convention
                    source_components.gxs[beam_ind].real = jones[0,0].real
                    source_components.gxs[beam_ind].imag = jones[0,0].imag
                    source_components.Dxs[beam_ind].real = jones[0,1].real
                    source_components.Dxs[beam_ind].imag = jones[0,1].imag
                    source_components.Dys[beam_ind].real = jones[1,0].real
                    source_components.Dys[beam_ind].imag = jones[1,0].imag
                    source_components.gys[beam_ind].real = jones[1,1].real
                    source_components.gys[beam_ind].imag = jones[1,1].imag
    
    ##now handle flux related things    
    ##always shove things into the source as power, curve, list
        ## - chunk_flux_type_index is the index to access with etc
        ##   source_components.power_ref_freqs
        ##   source_components.curve_ref_stokesQ
        ## - chunk_comp_index is the index to access within e.g
        ##   source_components.ras
        
    for pow_ind, old_comp_ind in enumerate(map_components.power_orig_inds):
        
        source_components.power_comp_inds[pow_ind] = pow_ind
        
        source_components.power_ref_freqs[pow_ind] = REF_FREQ
        source_components.power_ref_stokesI[pow_ind] = main_table['NORM_COMP_PL'][old_comp_ind]
        source_components.power_SIs[pow_ind] = main_table['ALPHA_PL'][old_comp_ind]
        
        
    for cur_ind, old_comp_ind in enumerate(map_components.curve_orig_inds):
        
        source_components.curve_comp_inds[cur_ind] = n_powers + cur_ind
        
        source_components.curve_ref_freqs[cur_ind] = REF_FREQ
        source_components.curve_ref_stokesI[cur_ind] = main_table['NORM_COMP_CPL'][old_comp_ind]
        source_components.curve_SIs[cur_ind] = main_table['ALPHA_CPL'][old_comp_ind]
        source_components.curve_qs[cur_ind] = main_table['CURVE_CPL'][old_comp_ind]
    
    ##Add the Stokes I list flux information
    add_list_flux_info(CompTypes.I_LIST, n_lists,
                       'INT_FLX', main_table, map_components,
                       source_components,
                       list_comp_ind_offset=n_powers + n_curves)
        
    # ##only some people need major, minor, pas
    if comp_type == CompTypes.GAUSSIAN or comp_type == CompTypes.SHAPELET:
        
        for new_comp_ind, old_comp_ind in enumerate(comp_orig_inds):
        
            source_components.majors[new_comp_ind] = main_table['MAJOR_DC'][old_comp_ind]*D2R
            source_components.minors[new_comp_ind] = main_table['MINOR_DC'][old_comp_ind]*D2R
            source_components.pas[new_comp_ind] = main_table['PA_DC'][old_comp_ind]*D2R
        
    # ##now for the shaepelet only stuff
    # ##need to cycle through all the shapelet SOURCEs, find all that match
    # ##in the second `shape_table` hdu, and add their basis function info
    if comp_type == CompTypes.SHAPELET:
        ##We can consolidate over flux model type here given we have a table
        ##(which we can't do using a text file) so consolodate some arrays
        n_s_powers = len(map_components.power_shape_orig_inds)
        n_s_curves = len(map_components.curve_shape_orig_inds)
        n_s_lists = len(map_components.list_shape_orig_inds)
        
        s_comp_orig_inds = np.empty(n_s_powers + n_s_curves + n_s_lists, dtype=int)
        s_comp_orig_inds[:n_s_powers] = map_components.power_shape_orig_inds
        s_comp_orig_inds[n_s_powers:n_s_powers+n_s_curves] = map_components.curve_shape_orig_inds
        s_comp_orig_inds[n_s_powers+n_s_curves:] = map_components.list_shape_orig_inds
        
        n_b_powers = len(map_components.power_shape_basis_inds)
        n_b_curves = len(map_components.curve_shape_basis_inds)
        n_b_lists = len(map_components.list_shape_basis_inds)
        
        b_comp_orig_inds = np.empty(n_b_powers + n_b_curves + n_b_lists, dtype=int)
        b_comp_orig_inds[:n_b_powers] = map_components.power_shape_basis_inds
        b_comp_orig_inds[n_b_powers:n_b_powers+n_b_curves] = map_components.curve_shape_basis_inds
        b_comp_orig_inds[n_b_powers+n_b_curves:] = map_components.list_shape_basis_inds
        
        for new_b_ind, old_b_ind in enumerate(s_comp_orig_inds):
            
            shape_name = str(main_table['NAME'][old_b_ind])
            
            basis_indexes = np.where(shape_table['NAME'] == shape_name)[0]
            this_basis = basis_indexes[b_comp_orig_inds[new_b_ind]]
            
            source_components.n1s[new_b_ind] = shape_table['N1'][this_basis]
            source_components.n2s[new_b_ind] = shape_table['N2'][this_basis]
            source_components.shape_coeffs[new_b_ind] = shape_table['COEFF'][this_basis]

            
            comp_ind = np.where(np.array(main_table['NAME'][comp_orig_inds], dtype=str) == shape_name)[0][0]
            # print(comp_ind)
            source_components.param_indexes[new_b_ind] = comp_ind
            
            
    ##polarisation times--------------------------------------------------------
    
    n_stokesV_pol_frac = map_components.num_v_pol_fracs
    n_stokesV_power = map_components.num_v_powers
    n_stokesV_curve = map_components.num_v_curves
    n_stokesV_list = map_components.num_v_lists
    
    n_linpol_pol_frac = map_components.num_lin_pol_fracs
    n_linpol_power = map_components.num_lin_powers
    n_linpol_curve = map_components.num_lin_curves
    n_linpol_list = map_components.num_lin_lists
    n_linpol_p_list = map_components.num_lin_p_lists
    
    ##TODO maybe one day in a far off future, if someone is doing a big
    ##linear diffuse sky, it might be worth adding the option to only so
    ##stokesI, stokesQ, stokesU, and not stokesV. This would involve updating
    ##some GPU code to not bother using V when calculating XX, YY, XY, YX
    ##If we have polarisation information
    if n_stokesV_pol_frac + n_stokesV_power + n_stokesV_curve + n_stokesV_list + \
        n_linpol_pol_frac + n_linpol_power + n_linpol_curve + n_linpol_list + n_linpol_p_list > 0:
        source_components.do_QUV = 1
    else:
        source_components.do_QUV = 0
        
    v_pol_frac_inds = map_components.v_pol_frac_orig_inds
    v_power_inds = map_components.v_power_orig_inds
    v_curve_inds = map_components.v_curve_orig_inds
    lin_pol_frac_inds = map_components.lin_pol_frac_orig_inds
    lin_power_inds = map_components.lin_power_orig_inds
    lin_curve_inds = map_components.lin_curve_orig_inds
    lin_p_list_inds = map_components.lin_p_list_orig_inds
    
    if n_stokesV_pol_frac:
        ##find the indexes of all the polarisation fraction relative to the
        ##components in this chunk
        chunk_inds = find_all_indexes_of_x_in_y(v_pol_frac_inds, comp_orig_inds)
        ##iterate and fill in information
        for this_ind, orig_ind, chunk_ind in zip(np.arange(n_stokesV_pol_frac), v_pol_frac_inds, chunk_inds):
            source_components.stokesV_pol_frac_comp_inds[this_ind] = chunk_ind
            source_components.stokesV_pol_fracs[this_ind] = main_table['V_POL_FRAC'][orig_ind]
        
    if n_stokesV_power:
        chunk_inds = find_all_indexes_of_x_in_y(v_power_inds, comp_orig_inds)
        for this_ind, orig_ind, chunk_ind in zip(np.arange(n_stokesV_power), v_power_inds, chunk_inds):
            source_components.stokesV_power_comp_inds[this_ind] = chunk_ind
            source_components.stokesV_power_ref_flux[this_ind] = main_table['V_NORM_COMP_PL'][orig_ind]
            source_components.stokesV_power_SIs[this_ind] = main_table['V_ALPHA_PL'][orig_ind]
            
    if n_stokesV_curve:
        chunk_inds = find_all_indexes_of_x_in_y(v_curve_inds, comp_orig_inds)
        for this_ind, orig_ind, chunk_ind in zip(np.arange(n_stokesV_curve), v_curve_inds, chunk_inds):
            source_components.stokesV_curve_comp_inds[this_ind] = chunk_ind
            source_components.stokesV_curve_ref_flux[this_ind] = main_table['V_NORM_COMP_CPL'][orig_ind]
            source_components.stokesV_curve_SIs[this_ind] = main_table['V_ALPHA_CPL'][orig_ind]
            source_components.stokesV_curve_qs[this_ind] = main_table['V_CURVE_CPL'][orig_ind]
            
    if n_stokesV_list:
        ##Add the Stokes V list flux information
        add_list_flux_info(CompTypes.V_LIST, n_stokesV_list,
                        'V_INT_FLX', v_table, map_components,
                        source_components, main_table=main_table,
                        comp_orig_inds=comp_orig_inds)
    
    ##The RM and instrinsic polarisation angle are need for all types of
    ##linear polarisation model, so keep track of them using `linpol_ind`
    linpol_ind = 0
    if n_linpol_pol_frac:
        ##find the indexes of all the polarisation fraction relative to the
        ##components in this chunk
        chunk_inds = find_all_indexes_of_x_in_y(lin_pol_frac_inds, comp_orig_inds)
        ##iterate and fill in information
        for this_ind, orig_ind, chunk_ind in zip(np.arange(n_linpol_pol_frac), lin_pol_frac_inds, chunk_inds):
            source_components.linpol_pol_frac_comp_inds[this_ind] = chunk_ind
            source_components.linpol_pol_fracs[this_ind] = main_table['LIN_POL_FRAC'][orig_ind]
            source_components.linpol_angle_inds[linpol_ind] = chunk_ind
            source_components.rm_values[linpol_ind] = main_table['RM'][orig_ind]
            if chunk_map.has_intr_pol_angle:
                source_components.intr_pol_angle[linpol_ind] = main_table['INTR_POL_ANGLE'][orig_ind]
                ##TODO might need to catch empty entries (NaNs) are convert to 0?
            else:
                source_components.intr_pol_angle[linpol_ind] = 0.0
            linpol_ind += 1
        
    if n_linpol_power:
        chunk_inds = find_all_indexes_of_x_in_y(lin_power_inds, comp_orig_inds)
        for this_ind, orig_ind, chunk_ind in zip(np.arange(n_linpol_power), lin_power_inds, chunk_inds):
            source_components.linpol_power_comp_inds[this_ind] = chunk_ind
            source_components.linpol_power_ref_flux[this_ind] = main_table['LIN_NORM_COMP_PL'][orig_ind]
            source_components.linpol_power_SIs[this_ind] = main_table['LIN_ALPHA_PL'][orig_ind]
            source_components.linpol_angle_inds[linpol_ind] = chunk_ind
            source_components.rm_values[linpol_ind] = main_table['RM'][orig_ind]
            if chunk_map.has_intr_pol_angle:
                source_components.intr_pol_angle[linpol_ind] = main_table['INTR_POL_ANGLE'][orig_ind]
            else:
                source_components.intr_pol_angle[linpol_ind] = 0.0
            linpol_ind += 1
            
    if n_linpol_curve:
        chunk_inds = find_all_indexes_of_x_in_y(lin_curve_inds, comp_orig_inds)
        for this_ind, orig_ind, chunk_ind in zip(np.arange(n_linpol_curve), lin_curve_inds, chunk_inds):
            source_components.linpol_curve_comp_inds[this_ind] = chunk_ind
            source_components.linpol_curve_ref_flux[this_ind] = main_table['LIN_NORM_COMP_CPL'][orig_ind]
            source_components.linpol_curve_SIs[this_ind] = main_table['LIN_ALPHA_CPL'][orig_ind]
            source_components.linpol_curve_qs[this_ind] = main_table['LIN_CURVE_CPL'][orig_ind]
            source_components.linpol_angle_inds[linpol_ind] = chunk_ind
            source_components.rm_values[linpol_ind] = main_table['RM'][orig_ind]
            if chunk_map.has_intr_pol_angle:
                source_components.intr_pol_angle[linpol_ind] = main_table['INTR_POL_ANGLE'][orig_ind]
            else:
                source_components.intr_pol_angle[linpol_ind] = 0.0
            linpol_ind += 1
            
    if n_linpol_list:
        add_list_flux_info(CompTypes.Q_LIST, n_linpol_list,
                        'Q_INT_FLX', q_table, map_components,
                        source_components, main_table=main_table,
                        comp_orig_inds=comp_orig_inds)
        
        add_list_flux_info(CompTypes.U_LIST, n_linpol_list,
                        'U_INT_FLX', u_table, map_components,
                        source_components, main_table=main_table,
                        comp_orig_inds=comp_orig_inds)
        
    if n_linpol_p_list:
        add_list_flux_info(CompTypes.LIN_LIST, n_linpol_p_list,
                        'P_INT_FLX', p_table, map_components,
                        source_components, main_table=main_table,
                        comp_orig_inds=comp_orig_inds)
        
        chunk_inds = find_all_indexes_of_x_in_y(lin_p_list_inds, comp_orig_inds)
        for this_ind, orig_ind, chunk_ind in zip(np.arange(n_linpol_p_list), lin_p_list_inds, chunk_inds):
            source_components.linpol_angle_inds[linpol_ind] = chunk_ind
            source_components.rm_values[linpol_ind] = main_table['RM'][orig_ind]
            if chunk_map.has_intr_pol_angle:
                source_components.intr_pol_angle[linpol_ind] = main_table['INTR_POL_ANGLE'][orig_ind]
            else:
                source_components.intr_pol_angle[linpol_ind] = 0.0
            linpol_ind += 1
        
    return


# @profile
def read_fits_skymodel_chunks(woden_struct_classes : Woden_Struct_Classes,
                              woden_settings : Woden_Settings,
                              args : argparse.Namespace,
                              main_table : Table, shape_table : Table,
                              chunked_skymodel_maps : list,
                              num_freqs : int, num_time_steps : int,
                              beamtype : int,
                              lsts : np.ndarray, latitude : float,
                              v_table : Table = False, q_table : Table = False,
                              u_table : Table = False, p_table : Table = False,
                              precision = "double"):
    """
    Uses Tables read from a FITS file and returns a source catalogue
    that can be used by C/CUDA code to calculate visibilities. Uses the
    maps in `chunked_skymodel_maps` to determine which components to read in
    and add to the `source_catalogue`

    Parameters
    ----------
    woden_struct_classes : Woden_Struct_Classes
        A class containing all the ctypes structures with the correct precision, needed for the C/CUDA code.
    main_table : Table
        The main Table (with RA,Dec etc) from the FITS file.
    shape_table : Table
        The shapelet Table from the FITS file.
    chunked_skymodel_maps : list
        List of ChunkedSkyModelMap objects, each representing a chunk of the sky model.
    num_freqs : int
        Number of frequency channels in the sky model.
    num_time_steps : int
        Number of time steps in the sky model.
    beamtype : int
        Type of beam used in the sky model.
    lsts : np.ndarray
        Array of LST values for each time step in the sky model.
    latitude : float
        Latitude of the observation site.
    precision : str, optional
        Precision of the source catalogue (either "float" or "double"), by default "double".

    Returns
    -------
    source_catalogue : Union[Source_Catalogue_Float, Source_Catalogue_Double]
        A source catalogue that can be used by C/CUDA code to calculate visibilities.
    """
    
    ##want to know how many shapelets are in all the chunks (used later
    # by "calculate_visiblities.cu")
    num_shapelets = 0
    for chunk_map in chunked_skymodel_maps:
        num_shapelets += chunk_map.n_shapes
        
    ##setup the source catalogue, which is going to store all of the information
    ##of each source and be fed straight into C/CUDA
    source_catalogue = setup_source_catalogue(woden_struct_classes.Source,
                                              woden_struct_classes.Source_Catalogue,
                                              len(chunked_skymodel_maps), num_shapelets,
                                              precision = precision)
    
    ##if we have an everybeam primary beam, we will be calculating it
    ##as we load in the sky model, as it happens on the CPU. So need to set
    ##some extra arguments here
    
    ra0 = woden_settings.ra0
    dec0 = woden_settings.dec0
    eb_beams = [BeamTypes.EB_OSKAR.value]
    
    if beamtype in eb_beams:
        all_times = []
        obs_time = Time(args.date, scale='utc')
        for time_step in range(woden_settings.num_time_steps):
            time_current = obs_time + TimeDelta((time_step + 0.5)*woden_settings.time_res, format='sec')
            all_times.append(time_current)
        
        band_num = woden_settings.band_nums[0]
        
        base_band_freq = ((band_num - 1)*woden_settings.coarse_band_width) + woden_settings.base_low_freq
        all_freqs = base_band_freq + np.arange(woden_settings.num_freqs)*woden_settings.frequency_resolution
        
        if beamtype == BeamTypes.EB_OSKAR.value:
            telescope = load_OSKAR_telescope(args.beam_ms_path)
            
        beam_norms = np.empty((len(all_times), len(all_freqs), 2))
        
        for time_ind, time in enumerate(all_times):
            for freq_ind, freq in enumerate(all_freqs):
                norm_x, norm_y = get_everybeam_norm(ra0, dec0, time, freq, telescope)
                beam_norms[time_ind, freq_ind, 0] = 1 / norm_x
                beam_norms[time_ind, freq_ind, 1] = 1 / norm_y
                # print(time_ind, freq_ind, beam_norms[time_ind, freq_ind])
        
    else:
        beam_norms = False
        telescope = False
        all_times = False
        all_freqs = False
        
    ##for each chunk map, create a Source_Float or Source_Double ctype
    ##struct, and "malloc" the right amount of arrays to store required infor
    for chunk_ind, chunk_map in enumerate(chunked_skymodel_maps):
        
        setup_chunked_source(source_catalogue.sources[chunk_ind], chunk_map,
                             num_freqs, num_time_steps, beamtype, precision=precision)
        
        ##count up the total number of components across all chunks
        ##annoyingly, beacuse Jack sucks, we split shapelet us by basis 
        ##and not component, so this number is actually a combination of
        ##component and basis numbers
        # num_comps_all_chunks += chunk_map.n_points + chunk_map.n_gauss + chunk_map.n_shape_coeffs
        
        if chunk_map.n_points > 0:
            add_fits_info_to_source_catalogue(CompTypes.POINT,
                                      main_table, shape_table,
                                      source_catalogue.sources[chunk_ind], chunk_map,
                                      beamtype, lsts, latitude,
                                      v_table, q_table, u_table, p_table,
                                      ra0, dec0, telescope, beam_norms,
                                      all_times, all_freqs)
            
        if chunk_map.n_gauss > 0:
            add_fits_info_to_source_catalogue(CompTypes.GAUSSIAN,
                                      main_table, shape_table,
                                      source_catalogue.sources[chunk_ind], chunk_map,
                                      beamtype, lsts, latitude,
                                      v_table, q_table, u_table, p_table,
                                      ra0, dec0, telescope, beam_norms,
                                      all_times, all_freqs)
            
        if chunk_map.n_shapes > 0:
            add_fits_info_to_source_catalogue(CompTypes.SHAPELET,
                                      main_table, shape_table,
                                      source_catalogue.sources[chunk_ind], chunk_map,
                                      beamtype, lsts, latitude,
                                      v_table, q_table, u_table, p_table,
                                      ra0, dec0, telescope, beam_norms,
                                      all_times, all_freqs)
            
    ##TODO some kind of consistency check between the chunk_maps and the
    ##sources in the catalogue - make sure we read in the correct information
    
    return source_catalogue

