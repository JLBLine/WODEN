from wodenpy.skymodel.read_yaml_skymodel import read_yaml_radec_count_components
from wodenpy.skymodel.read_fits_skymodel import read_fits_skymodel_chunks, read_fits_radec_count_components, check_columns_fits
from wodenpy.skymodel.read_text_skymodel import read_text_radec_count_components
# from wodenpy.use_libwoden.skymodel_structs import Source_Catalogue
from wodenpy.skymodel.woden_skymodel import Component_Type_Counter
import numpy as np
from typing import Union, List
import sys
from astropy.table import Table
from wodenpy.skymodel.read_text_skymodel import read_full_text_into_fitstable
from wodenpy.skymodel.read_yaml_skymodel import read_full_yaml_into_fitstable
from wodenpy.use_libwoden.create_woden_struct_classes import Woden_Struct_Classes
from wodenpy.use_libwoden.skymodel_structs import Source_Python, Skymodel_Chunk_Map, setup_source_catalogue, copy_python_source_to_ctypes
from astropy.io import fits
import argparse


##This call is so we can use it as a type annotation
woden_struct_classes = Woden_Struct_Classes()
Woden_Settings = woden_struct_classes.Woden_Settings

def read_radec_count_components(skymodel_path : str) -> Component_Type_Counter:
    """
    Reads basic info from a sky model. Checks whether the sky model is a .fits,
    .yaml, or .txt file, and errors if not. Then calls the appropriate function
    if so. Reads just the ra, dec, and counts how many POINT/GAUSS/SHAPE and
    POWER/CURVE/LIST, consolidating the information into a
    Component_Type_Counter object.
    
    Parameters
    -----------
    skymodel_path : str
        The path to the sky model file.
    
    Returns
    --------
    comp_counter : Component_Type_Counter
        A Component_Type_Counter object that counts the number of components of each type and their properties.
    """
    
    ##Figure out if our skymodel is supported or not
    if skymodel_path[-5:] == '.fits':
        check_columns_fits(skymodel_path)
        comp_counter = read_fits_radec_count_components(skymodel_path)
    
    elif skymodel_path[-5:] == '.yaml':
    
        comp_counter = read_yaml_radec_count_components(skymodel_path)

    elif skymodel_path[-4:] == '.txt':
    
        comp_counter = read_text_radec_count_components(skymodel_path)
        
    else:
        sys.exit('The filename fed into `wodenpy/read_skymodel/read_radec_count_components` was not of a supported file type. Currently supported formats are: .fits, .yaml, .txt')
        
    return comp_counter

def get_skymodel_tables(skymodel_path : str) -> tuple:
    """
    Read the contents of a sky model file into a set of astropy Tables. If the
    sky model is a .fits file, the tables are read in directly. If the sky model
    is a .yaml or .txt file, the file is read in an converted into Tables. As
    the yaml and text formats are deprecated, they never return any polarisation
    infomration.
    
    Parameters
    ----------
    skymodel_path : str
        The path to the sky model file. Must end in .fits, .yaml, or .txt.
    Returns
    -------
    tuple: 
        A tuple containing the following elements:
        - main_table (Table): The main table read from the first HDU.
        - shape_table (Table or str): The shapelet table if present, otherwise 'no_shape_table'.
        - v_table (Table or bool): The V_LIST_FLUXES table if present, otherwise False.
        - q_table (Table or bool): The Q_LIST_FLUXES table if present, otherwise False.
        - u_table (Table or bool): The U_LIST_FLUXES table if present, otherwise False.
        - p_table (Table or bool): The P_LIST_FLUXES table if present, otherwise False.
    """
    
    if skymodel_path[-5:] == '.fits':
        main_table = Table.read(skymodel_path, hdu=1)
        
        has_shapelets = False
        if 'S' in main_table['COMP_TYPE']:
            has_shapelets = True
        
        with fits.open(skymodel_path) as hdus:
            num_hdus = len(hdus)
            hdu_names = [hdu.name for hdu in hdus]
            
        if 'SHAPELET' in hdu_names:
            shape_table = Table.read(skymodel_path, hdu='SHAPELET')
        elif has_shapelets:
            shape_table = Table.read(skymodel_path, hdu=2)
        else:
            shape_table = 'no_shape_table'
            
        if 'V_LIST_FLUXES' in hdu_names:
                v_table = Table.read(skymodel_path, hdu='V_LIST_FLUXES')
        else:
            v_table = False
            
        if 'Q_LIST_FLUXES' in hdu_names:
                q_table = Table.read(skymodel_path, hdu='Q_LIST_FLUXES')
        else:
            q_table = False    
            
        if 'U_LIST_FLUXES' in hdu_names:
                u_table = Table.read(skymodel_path, hdu='U_LIST_FLUXES')
        else:
            u_table = False
            
        if 'P_LIST_FLUXES' in hdu_names:
                p_table = Table.read(skymodel_path, hdu='P_LIST_FLUXES')
        else:
            p_table = False
            
    elif skymodel_path[-5:] == '.yaml':
        main_table, shape_table = read_full_yaml_into_fitstable(skymodel_path)
        v_table = False
        q_table = False    
        u_table = False
        p_table = False
        
    elif skymodel_path[-4:] == '.txt':
        main_table, shape_table = read_full_text_into_fitstable(skymodel_path)
        v_table = False
        q_table = False    
        u_table = False
        p_table = False
        
    else:
        sys.exit('The filename fed into `wodenpy/read_skymodel/read_skymodel_chunks` was not of a supported file type. Currently supported formats are: .fits, .yaml, .txt')
        
    return main_table, shape_table, v_table, q_table, u_table, p_table

##just call this here so we can type annotate it
woden_struct_classes = Woden_Struct_Classes()
Source_Catalogue = woden_struct_classes.Source_Catalogue

def create_source_catalogue_from_python_sources(python_sources : List[Source_Python],
                                                chunked_skymodel_maps : List[Skymodel_Chunk_Map],
                                                woden_struct_classes : Woden_Struct_Classes,
                                                beamtype : int, precision : str = 'double') -> Source_Catalogue: # type: ignore
    
    num_shapelets = 0

    for chunk_map in chunked_skymodel_maps:
        num_shapelets += chunk_map.n_shapes
        
    ##setup the source catalogue, which is going to store all of the information
    ##of each source and be fed straight into C/CUDA
    source_catalogue = setup_source_catalogue(woden_struct_classes.Source_Ctypes,
                                              woden_struct_classes.Source_Catalogue,
                                              len(chunked_skymodel_maps), num_shapelets,
                                              precision = precision)
    
    
    for chunk_ind, chunk_map, python_source in zip(range(len(chunked_skymodel_maps)),
                                                    chunked_skymodel_maps,
                                                    python_sources):
    
        copy_python_source_to_ctypes(python_source, source_catalogue.sources[chunk_ind],
                                     beamtype, precision=precision)
        
    return source_catalogue