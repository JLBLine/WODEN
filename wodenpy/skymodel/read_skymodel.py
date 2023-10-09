from wodenpy.skymodel.read_yaml_skymodel import read_yaml_skymodel_chunks, read_yaml_radec_count_components
from wodenpy.skymodel.read_fits_skymodel import read_fits_skymodel_chunks, read_fits_radec_count_components
from wodenpy.skymodel.read_text_skymodel import read_text_skymodel_chunks, read_text_radec_count_components
from wodenpy.use_libwoden.skymodel_structs import Source_Catalogue_Float, Source_Catalogue_Double
from wodenpy.skymodel.woden_skymodel import Component_Type_Counter
import numpy as np
from typing import Union
import sys


def read_radec_count_components(skymodel_path : str) -> Component_Type_Counter:
    
    ##Figure out if our skymodel is supported or not
    if skymodel_path[-5:] == '.fits':
        comp_counter = read_fits_radec_count_components(skymodel_path)
    
    elif skymodel_path[-5:] == '.yaml':
    
        comp_counter = read_yaml_radec_count_components(skymodel_path)

    elif skymodel_path[-4:] == '.txt':
    
        comp_counter = read_text_radec_count_components(skymodel_path)
        
    else:
        sys.exit('The filename fed into `wodenpy/read_skymodel/read_radec_count_components` was not of a supported file type. Currently supported formats are: .fits, .yaml, .txt')
        
    return comp_counter


def read_skymodel_chunks(skymodel_path : str, chunked_skymodel_maps : list,
                         num_freqs : int, num_time_steps : int,
                         beamtype : int,
                         lsts : np.ndarray, latitude : float,
                         precision = "double") -> Union[Source_Catalogue_Float, Source_Catalogue_Double]:
    
    # print("HALLO", skymodel_path[-5:])

    # ##Figure out if our skymodel is supported or not
    if skymodel_path[-5:] == '.fits':
        source_catalogue = read_fits_skymodel_chunks(skymodel_path, chunked_skymodel_maps,
                              num_freqs, num_time_steps, beamtype,
                              lsts, latitude, precision = precision)
    
    elif skymodel_path[-5:] == '.yaml':
    
        source_catalogue = read_yaml_skymodel_chunks(skymodel_path, chunked_skymodel_maps,
                              num_freqs, num_time_steps, beamtype,
                              lsts, latitude, precision = precision)
        
    elif skymodel_path[-4:] == '.txt':
    
        source_catalogue = read_text_skymodel_chunks(skymodel_path, chunked_skymodel_maps,
                              num_freqs, num_time_steps, beamtype,
                              lsts, latitude, precision = precision)
        
    else:
        sys.exit('The filename fed into `wodenpy/read_skymodel/read_skymodel_chunks` was not of a supported file type. Currently supported formats are: .fits, .yaml, .txt')
        
    return source_catalogue