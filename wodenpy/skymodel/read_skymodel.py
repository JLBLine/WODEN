from wodenpy.skymodel.read_yaml_skymodel import read_yaml_radec_count_components
from wodenpy.skymodel.read_fits_skymodel import read_fits_skymodel_chunks, read_fits_radec_count_components, check_columns_fits
from wodenpy.skymodel.read_text_skymodel import read_text_radec_count_components
# from wodenpy.use_libwoden.skymodel_structs import Source_Catalogue
from wodenpy.skymodel.woden_skymodel import Component_Type_Counter
import numpy as np
from typing import Union
import sys
from astropy.table import Table
from wodenpy.skymodel.read_text_skymodel import read_full_text_into_fitstable
from wodenpy.skymodel.read_yaml_skymodel import read_full_yaml_into_fitstable
from wodenpy.use_libwoden.create_woden_struct_classes import Woden_Struct_Classes
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


def read_skymodel_chunks(woden_struct_classes : Woden_Struct_Classes,
                         woden_settings : Woden_Settings,
                         args : argparse.Namespace,
                         skymodel_path : str, chunked_skymodel_maps : list,
                         num_freqs : int, num_time_steps : int,
                         beamtype : int,
                         lsts : np.ndarray, latitude : float,
                         precision = "double"):
    """Lazy loads chunks of a sky model at `skymodel_path` into a
    Source_Catalogue object, as mapped by `chunked_skymodel_maps`. If the
    sky model isn't already a FITS file, it is converted to one. The
    resultant can be passed into C/CUDA code to calculate visibilities.

    Parameters
    ----------
    woden_struct_classes : Woden_Struct_Classes
        This holds all the various ctype structure classes that are equivalent
        to the C/CUDA structs.
    skymodel_path : str
        The path to the sky model file.
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
    
    ##TODO reading all these tables should happen outside of this function
    ##for small yaml/txt skymodels it's not much overhead, but for large ones
    ## I/O is gonna be a limiting factor
    # ##Figure out if our skymodel is supported or not
    if skymodel_path[-5:] == '.fits':
        main_table = Table.read(skymodel_path, hdu=1)
        
        num_shapelets = 0
        for chunk_map in chunked_skymodel_maps:
            num_shapelets += chunk_map.n_shapes
            
        with fits.open(skymodel_path) as hdus:
            num_hdus = len(hdus)
            hdu_names = [hdu.name for hdu in hdus]
            
        ##Only read in a shape table if there are shapelets
        if num_shapelets:
            if 'SHAPELET' in hdu_names:
                shape_table = Table.read(skymodel_path, hdu='SHAPELET')
            else:
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
        
    source_catalogue = read_fits_skymodel_chunks(woden_struct_classes,
                              woden_settings, args,
                              main_table, shape_table,
                              chunked_skymodel_maps,
                              num_freqs, num_time_steps, beamtype,
                              lsts, latitude,
                              v_table, q_table, u_table, p_table,
                              precision = precision)
    
    # print("HERE", type(source_catalogue))
        
    return source_catalogue