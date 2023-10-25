import numpy as np
import sys
import os
from typing import Union, Tuple

from wodenpy.skymodel.woden_skymodel import Component_Type_Counter, Component_Info, CompTypes, calc_pl_norm_at_200MHz, calc_cpl_norm_at_200MHz
from wodenpy.skymodel.chunk_sky_model import Skymodel_Chunk_Map
from wodenpy.use_libwoden.skymodel_structs import setup_chunked_source, setup_source_catalogue, Source_Catalogue_Float, Source_Catalogue_Double, _Ctype_Source_Into_Python, Components_Float, Components_Double, Source_Float, Source_Double
from wodenpy.skymodel.read_fits_skymodel import add_fits_info_to_source_catalogue

from wodenpy.use_libwoden.beam_settings import BeamTypes

from astropy.io import fits
from astropy.table import Table, Column

import erfa

D2R = np.pi/180.0

def read_text_radec_count_components(text_path : str):
    """
    Reads a text file sky model. Reads just the ra, dec, and counts how many POINT/GAUSS/SHAPE and POWER/CURVE/LIST, consolidating the
    information into a Component_Type_Counter object.
    
    Parameters
    -----------
    text_path : str
        The path to the text file containing the sky model information.
    
    Returns
    --------
    comp_counter : Component_Type_Counter
        A Component_Type_Counter object that counts the number of components of each type and their properties.
    """
    
    if not os.path.isfile(text_path):
        sys.exit(f"Cannot read sky model from {text_path}. Please check your paths, exiting now.")

    comp_counter = Component_Type_Counter()

    with open(text_path) as file:
        
        current_source = -1
        first_component = True
        
        ##These are all things we have to count to be able to malloc correct
        ##amount of memory down the line
        for line in file:
            comp_counter.new_file_line()
            if line != '\n\n' and '#' not in line and line != ''  and line != ' ' and line != '\n':
                
                if 'SOURCE' in line and 'END' not in line:
                    current_source += 1
                    source_name = line.split()
                    comp_counter.new_source()
                
                elif 'COMPONENT' in line and 'END' not in line:
                    if first_component:
                        first_component = False
                    else:
                        
                        ##ra should be the first thing in a component, so we need
                        ##to append all the previously found values and reset the
                        ##counters
                        comp_counter.add_component_reset_counters()
                        
                    ##Start a new component
                    comp_counter.initiate_component()
                    
                    if len(line.split()) != 4:
                        sys.exit(f'error reading line {comp_counter.file_line_num} of skymodel {text_path}')
                        
                    _, comp_type, ra, dec = line.split()
                    
                    
                    comp_counter.ra = float(ra)*D2R*15.0
                    comp_counter.dec = float(dec)*D2R
                    
                    ##everything is a power law in the WODEN text sky model
                    comp_counter.flux_power = 1
                    
                    if comp_type == 'POINT':
                        comp_counter.point = 1
                    elif comp_type == 'GAUSSIAN':
                        comp_counter.gaussian = 1
                    elif comp_type == 'SHAPELET':
                        comp_counter.shapelet = 1
                        
                elif 'SCOEFF' in line:
                    comp_counter.num_shape_basis += 1

        ##The final component needs to be counted as we used a trigger of the
        ##next to count the pervious components
        comp_counter.add_component_reset_counters()
    
    ##shrink the array length of all info to the amount actually in the model
    comp_counter.remove_array_padding()
    
    ##total up some component counts
    comp_counter.total_components()
        
    return comp_counter


def read_full_text_into_fitstable(text_path : str) -> Tuple[Table, Table]:
    """Read in an old-school WODEN style text sky model and convert it into
    the LoBES-style FITS sky model format. Can then be read in by
    `read_fits_skymodel.read_fits_skymodel_chunks`.

    Parameters
    ----------
    text_path : str
        Path to the text sky model

    Returns
    -------
    Tuple(Table, Table)
        The `main_table` and `shape_table` to be used by during lazy loading
    """
    
    main_table = False
    shape_table = False
    
    all_freqs = []
    all_names = []

    with open(text_path) as file:

        components = []
        sources = []

        component = False
        source_name = False
        current_source = 0
        comp_count = 0

        source_indexes = []

        freq_count = 0
        freq_indent = 0
        
        for line in file:
        
            if line != '\n\n' and '#' not in line and line != ''  and line != ' ' and line != '\n':
                
                if 'SOURCE' in line and 'END' not in line:
                    current_source += 1
                    source_name = line.split()
                    all_names.append(source_name)
                    current_source += 1
                    comp_count = -1
                    
                elif 'COMPONENT' in line and 'END' not in line:

                    ##COMPONENT should be the first thing in a component, so we need
                    ##to append all the previously found values and reset the
                    ##counters

                    ##If a previous component exists, append in to the list
                    ##of all components, and then make a new one
                    if component:
                        ##Make some things into np.arrays so we can maths them
                        component.n1s = np.array(component.n1s)
                        component.n2s = np.array(component.n2s)
                        component.coeffs = np.array(component.coeffs)
                        component.fluxes = np.array(component.fluxes)
                        component.freqs = np.array(component.freqs)
                        component.comp_name = f"{component.source_name}_C{component.comp_count:02d}"
                        
                        components.append(component)

                    component = Component_Info()
                    freq_count = 0

                    component.source_name = source_name
                    comp_count += 1
                    component.comp_count = comp_count
                    component.ra = float(line.split()[-1])

                    source_indexes.append(current_source)

                    _, comp_type, ra, dec = line.split()
                    
                    component.ra = float(ra)*15.0
                    component.dec = float(dec)
                    
                    if comp_type == 'POINT':
                        component.comp_type = 'P'
                    elif comp_type == 'GAUSSIAN':
                        component.comp_type = 'G'
                    elif comp_type == 'SHAPELET':
                        component.comp_type = 'S'
                    
                    ##everything is a power law in the WODEN text sky model
                    component.flux_type = 'pl'

                elif 'SPARAMS' in line or 'GPARAMS' in line:
                    _, pa, major, minor = line.split()

                    component.major = float(major) / 60.0
                    component.minor = float(minor) / 60.0
                    component.pa = float(pa)
                    
                elif 'SCOEFF' in line:
                    _, n1, n2, coeff = line.split()

                    component.n1s.append(float(n1))
                    component.n2s.append(float(n2))
                    component.coeffs.append(float(coeff))

                ##Assume a powerlaw with SI = -0.8 for lines with FREQ
                ##Otherwise read in the SI
                elif 'FREQ' in line or 'LINEAR' in line:
                    freq_count += 1
                    
                    if 'FREQ' in line:
                    
                        _, freq, stokesI, stokesQ, stokesU, stokesV = line.split()
                        si = -0.8
                    else:
                        _, freq, stokesI, stokesQ, stokesU, stokesV, si = line.split()
                        si = float(si)
                        
                    component.freqs.append(float(freq))
                    
                    ##Stick in a Stokes I,Q,U,V
                    component.fluxes.append(np.array([float(stokesI),
                                                      float(stokesQ),
                                                      float(stokesU),
                                                      float(stokesV)]))
                    
                    component.si = si

    ##last one doesn't get added to list, so do that
    component.n1s = np.array(component.n1s)
    component.n2s = np.array(component.n2s)
    component.coeffs = np.array(component.coeffs)
    component.fluxes = np.array(component.fluxes)
    component.freqs = np.array(component.freqs)
    component.comp_name = f"{component.source_name}_C{component.comp_count:02d}"

    components.append(component)

    all_freqs = np.sort(np.unique(all_freqs))
    
    components = np.array(components)
    num_components = len(components)

    flux_types = np.array([component.flux_type for component in components])
    
    power_laws = np.where(flux_types == 'pl')[0]

    for comp_ind, component in enumerate(components[power_laws]):
        component = calc_pl_norm_at_200MHz(component)
        
    ##for all components, what SOURCE do they belong to?
    comp_source_names = np.array([component.source_name for component in components])
    comp_names = np.array([component.comp_name for component in components])

    unq_source_ID = Column(data=comp_source_names, name='UNQ_SOURCE_ID')
    name = Column(data=comp_names, name='NAME')


    ras = Column(data=np.array([component.ra for component in components]),
                           name='RA', unit='deg')

    decs = Column(data=np.array([component.dec for component in components]),
                           name='DEC', unit='deg')

    majors = Column(data=np.array([component.major for component in components]),
                           name='MAJOR_DC', unit='deg')

    minors = Column(data=np.array([component.minor for component in components]),
                           name='MINOR_DC', unit='deg')

    pas = Column(data=np.array([component.pa for component in components]),
                           name='PA_DC', unit='deg')

    # main_table.add_columns([unq_source_ID, name])

    mod_type = Column(data=np.array([component.flux_type for component in components]),
                           name='MOD_TYPE', unit='deg')

    norm_comp_pl = Column(data=np.array([component.norm_comp_pl for component in components]),
                          name="NORM_COMP_PL")
    alpha_pl = Column(data=np.array([component.si for component in components]),
                          name="ALPHA_PL")
    norm_comp_cpl = Column(data=np.array([component.norm_comp_cpl for component in components]),
                          name="NORM_COMP_CPL")
    alpha_cpl = Column(data=np.array([component.si for component in components]),
                          name="ALPHA_CPL")
    curve_cpl = Column(data=np.array([component.curve_q for component in components]),
                          name="CURVE_CPL")
    comp_type = Column(data=np.array([component.comp_type for component in components]),
                          name="COMP_TYPE")

    out_columns = [unq_source_ID, name, ras, decs, majors, minors, pas, mod_type, comp_type,
                   norm_comp_pl, alpha_pl, norm_comp_cpl, alpha_cpl, curve_cpl]

    for freq in all_freqs:
        flux_data = np.full(num_components, np.nan)

        for comp_ind, component in enumerate(components):
            # print(component.fluxes.shape)
            fluxes = component.fluxes[np.where(component.freqs == freq)[0]]
            if len(fluxes) == 1:
                stokesI = fluxes[0][0]
                flux_data[comp_ind] = stokesI

        flux_col = Column(data=flux_data, name=f"INT_FLX{freq/1e+6:.3f}", unit='Jy')
        out_columns.append(flux_col)

    main_table = Table()
    main_table.add_columns(out_columns)
    # main_table.write(, overwrite=True)

    ##gather the shapelet specific information

    shape_names = []
    shape_n1s = []
    shape_n2s = []
    shape_coeffs = []

    for component in components:
        if component.comp_type == 'S':
            for n1, n2, coeff in zip(component.n1s, component.n2s,
                                     component.coeffs):

                shape_names.append(component.comp_name)
                shape_n1s.append(n1)
                shape_n2s.append(n2)
                shape_coeffs.append(coeff)

    s_names = Column(data=np.array(shape_names), name="NAME")
    s_n1s = Column(data=np.array(shape_n1s, dtype=int), name="N1")
    s_n2s = Column(data=np.array(shape_n2s, dtype=int), name="N2")
    s_coeffs = Column(data=np.array(shape_coeffs), name="COEFF")

    shape_table = Table()
    shape_table.add_columns([s_names, s_n1s, s_n2s, s_coeffs])

    # hdu_list = fits.HDUList([
    #     fits.PrimaryHDU(),
    #     fits.table_to_hdu(main_table),
    #     fits.table_to_hdu(shape_table),
    # ])
    
    ##TODO - option to write out FITS version of input model?
    # hdu_list.writeto('converted_input.fits', overwrite=True)
    
    return main_table, shape_table