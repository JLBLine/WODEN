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

def read_yaml_radec_count_components(yaml_path : str):
    """
    Reads a yaml file sky model. Reads just the ra, dec, and counts how many POINT/GAUSS/SHAPE and POWER/CURVE/LIST, consolidating the
    information into a Component_Type_Counter object.
    
    Parameters
    -----------
    yaml_path : str
        The path to the yaml file containing the sky model information.
    
    Returns
    --------
    comp_counter : Component_Type_Counter
        A Component_Type_Counter object that counts the number of components of each type and their properties.
    """
    
    if not os.path.isfile(yaml_path):
        sys.exit(f"Cannot read sky model from {yaml_path}. Please check your paths, exiting now.")

    comp_counter = Component_Type_Counter()

    with open(yaml_path) as file:
        
        current_source = -1
        first_component = True
        
        ##These are all things we have to count to be able to malloc correct
        ##amount of memory down the line
        for line in file:
            comp_counter.new_file_line()
            if line != '---\n' and '#' not in line and line != ''  and line != ' ' and line != '\n':
                
                if line[0] != ' ' and line[0] != '-':
                    current_source += 1
                    source_name = line[:-1]
                    comp_counter.new_source()
                
                elif 'ra:' in line:
                    ##If this is the first component, no need to reset
                    ##and count things up
                    if first_component:
                        first_component = False
                    else:
                        
                        ##ra should be the first thing in a component, so we need
                        ##to append all the previously found values and reset the
                        ##counters
                        comp_counter.add_component_reset_counters()
                        
                    ##Start a new component
                    comp_counter.initiate_component()
                        
                    ra = float(line.split()[-1])*D2R
                    comp_counter.ra = ra
                    
                elif 'dec:' in line:
                    dec = float(line.split()[-1])*D2R
                    comp_counter.dec = dec

                ##component type related things
                elif 'point' in line:
                    comp_counter.point = 1
                elif 'gaussian:' in line:
                    comp_counter.gaussian = 1
                elif 'shapelet:' in line:
                    comp_counter.shapelet = 1
                elif 'n1:' in line:
                    comp_counter.num_shape_basis += 1
                
                ##flux behaviou related things
                elif 'power_law:' in line and 'curved' not in line:
                    comp_counter.flux_power = 1
                elif 'curved_power_law:' in line:
                    comp_counter.flux_curve = 1
                elif 'list:' in line:
                    comp_counter.flux_list = 1
                elif 'freq:' in line:
                    comp_counter.count_list_flux()

        ##The final component needs to be counted as we used a trigger of the
        ##next to count the pervious components
        comp_counter.add_component_reset_counters()
    
    ##shrink the array length of all info to the amount actually in the model
    comp_counter.remove_array_padding()
    
    ##total up some component counts
    comp_counter.total_components()
        
    return comp_counter

def read_full_yaml_into_fitstable(yaml_path : str) -> Tuple[Table, Table]:
    """Read in a hyperdrive style yaml sky model and convert it into
    the LoBES-style FITS sky model format. Can then be read in by
    `read_fits_skymodel.read_fits_skymodel_chunks`.

    Parameters
    ----------
    yaml_path : str
        Path to the yaml sky model

    Returns
    -------
    Tuple(Table, Table)
        The `main_table` and `shape_table` to be used by during lazy loading
    """
    
    main_table = False
    shape_table = False
    
    all_freqs = []
    all_names = []

    with open(yaml_path) as file:

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
            if line != '---\n' and '#' not in line and line != ''  and line != ' ' and line != '\n':

                if line[0] != ' ' and line[0] != '-':
                    # print(current_source)
                    source_name = line.split('\n')[0].strip(':')
                    all_names.append(source_name)
                    current_source += 1
                    comp_count = -1

                elif 'ra:' in line:

                    ##ra should be the first thing in a component, so we need
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

                elif 'dec:' in line:
                    component.dec = float(line.split()[-1])

                elif 'comp_type: point' in line:
                    component.comp_type = 'P'
                elif 'gaussian:' in line:
                    component.comp_type = 'G'
                elif 'shapelet:' in line:
                    component.comp_type = 'S'

                elif "maj:" in line:
                    component.major = float(line.split()[-1])*(1 / 3600.0)
                elif "min:" in line:
                    component.minor = float(line.split()[-1])*(1 / 3600.0)
                elif "pa:" in line:
                    component.pa = float(line.split()[-1])

                elif 'n1:' in line:
                    component.n1s.append(float(line.split()[-1]))
                elif 'n2:' in line:
                    component.n2s.append(float(line.split()[-1]))
                elif 'value:' in line:
                    component.coeffs.append(float(line.split()[-1]))

                elif 'power_law:' in line and 'curved' not in line:
                    component.flux_type = 'pl'
                elif 'curved_power_law:' in line:
                    component.flux_type = 'cpl'
                elif 'si:' in line:
                    component.si = float(line.split()[-1])

                elif 'list:' in line:
                    component.flux_type = 'nan'

                elif 'freq:' in line:
                    freq_count += 1
                    component.freqs.append(float(line.split()[-1]))
                    
                    ##OK, only want to write out flux column info if this
                    ##is a 'list' type source. So just append to all freqs if
                    ##correct type

                    if component.flux_type == 'nan':
                        all_freqs.append(float(line.split()[-1]))

                    ##Stick in an empty np.array for Stokes I,Q,U,V
                    component.fluxes.append(np.array([np.nan, np.nan, np.nan,np.nan]))

                    ##See what indent this freq entry starts at - used to
                    ##line up following freq entries, as `q` can either mean
                    ##stokes Q or q curvature param
                    freq_indent = line.index('f')


                elif ' i:' in line:
                    component.fluxes[freq_count - 1][0] = float(line.split()[-1])

                ##Gotta be fancy here to work out if this is a Stokes Q or a
                ##curved power law 'q' param
                elif ' q:' in line:
                    q = float(line.split()[-1])
                    if line.index('q') == freq_indent:
                        component.fluxes[freq_count - 1][1] = q
                    else:
                        if component.flux_type == 'cpl':
                            component.curve_q = q

                elif ' u:' in line:
                    component.fluxes[freq_count - 1][2] = float(line.split()[-1])

                elif ' v:' in line:
                    component.fluxes[freq_count - 1][3] = float(line.split()[-1])


    ##last one doesn't get added to list, so do that
    component.n1s = np.array(component.n1s)
    component.n2s = np.array(component.n2s)
    component.coeffs = np.array(component.coeffs)
    component.fluxes = np.array(component.fluxes)
    component.freqs = np.array(component.freqs)
    component.comp_name = f"{component.source_name}_C{component.comp_count:02d}"

    components.append(component)

    # all_freqs = np.unique(all_freqs).tolist()
    
    all_freqs = np.sort(np.unique(all_freqs))
    
    # print("HERE MAN", all_freqs)
    
    components = np.array(components)
    num_components = len(components)

    flux_types = np.array([component.flux_type for component in components])
    
    power_laws = np.where(flux_types == 'pl')[0]
    curve_laws = np.where(flux_types == 'cpl')[0]
    list_laws = np.where(flux_types == 'nan')[0]
    
    # print("Before fitting: num power, curved, list", len(power_laws), len(curve_laws), len(list_laws))

    for comp_ind, component in enumerate(components[power_laws]):
        component = calc_pl_norm_at_200MHz(component)
        
    for comp_ind, component in enumerate(components[curve_laws]):
        component = calc_cpl_norm_at_200MHz(component)
    
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