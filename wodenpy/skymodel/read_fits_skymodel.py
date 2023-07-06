import numpy as np
import sys
import os
from typing import Union

from wodenpy.skymodel.woden_skymodel import Component_Type_Counter, Component_Info, CompTypes
from wodenpy.skymodel.chunk_sky_model import Skymodel_Chunk_Map
from wodenpy.use_libwoden.skymodel_structs import setup_chunked_source, setup_source_catalogue, Source_Catalogue_Float, Source_Catalogue_Double, add_info_to_source_catalogue, _Ctype_Source_Into_Python, Components_Float, Components_Double, Source_Float, Source_Double
from wodenpy.use_libwoden.beam_settings import BeamTypes

from astropy.table import Table, Column
import erfa
from astropy.io import fits

REF_FREQ = 200e+6

D2R = np.pi/180.0

# @profile
def read_fits_radec_count_components(fits_path : str):
    """Read just the  ra, dec, and count how many POINT/GAUSS/SHAPE and
    POWER/CURVE/LIST entries there are"""
    
    if not os.path.isfile(fits_path):
        sys.exit(f"Cannot read sky model from {fits_path}. Please check your paths, exiting now.")
        
        
    ##grabd all the relevant information out of the tables
    
    main_table = Table.read(fits_path, hdu=1)
    
    
    ras = np.array(main_table['RA'], dtype=np.float64)
    decs = np.array(main_table['DEC'], dtype=np.float64)
    comp_types = np.array(main_table['COMP_TYPE'], dtype=str)
    flux_types = np.array(main_table['MOD_TYPE'], dtype=str)
    unq_source_ID = np.array(main_table['UNQ_SOURCE_ID'], dtype=str)
    comp_names = np.array(main_table['NAME'], dtype=str)
    
    flux_col_names = []
    for key in main_table.columns:
        if key[:7] == 'INT_FLX':
            flux_col_names.append(key)
            # flux_col_freqs.append(float(key[7:])*1e+6)
        
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
    
    list_type_inds = np.where(flux_types == 'nan')[0]
    for flux_col in flux_col_names:
        
        ##OK, astropy.Table can turns things into masked arrays. Here we are
        ##accessing a column `main_table[flux_col]`, reading it's `mask` which
        ##returns a bunch of booleans (True == masked). We want everything that
        ##isn't masked, so use the ~ to swap. Finally, we want ints, not booleans
        
        if type(main_table[flux_col]) == Column:
            present_fluxes = np.where(np.isnan(main_table[flux_col]) == False)[0]
            comp_counter.num_list_fluxes[list_type_inds] += 1
        else:
            present_fluxes = (~main_table[flux_col].mask).astype(int)
            comp_counter.num_list_fluxes[list_type_inds] += present_fluxes[list_type_inds]
    
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
        print("WARNING - couldn't find second table containing shapelet information, so not attempting to load any shapelets.")
    
    
    comp_counter.total_components()

    return comp_counter


def add_fits_info_to_source_catalogue(comp_type : CompTypes,
                        main_table : Table, shape_table : Table,
                        chunk_source : Skymodel_Chunk_Map,
                        chunk_map : Union[Source_Float, Source_Double],
                        num_freqs : int, num_time_steps : int,
                        beamtype : int,
                        lsts : np.ndarray, latitude : float,
                        precision='double'):
    
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
            
        ##only the NO_BEAM and GAUSS_BEAM options don't need az,za coords
        if beamtype == BeamTypes.GAUSS_BEAM.value or beamtype == BeamTypes.NO_BEAM.value:
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
        source_components.power_ref_stokesQ[pow_ind] = 0.0
        source_components.power_ref_stokesU[pow_ind] = 0.0
        source_components.power_ref_stokesV[pow_ind] = 0.0
        source_components.power_SIs[pow_ind] = main_table['ALPHA_PL'][old_comp_ind]
        
        
    for cur_ind, old_comp_ind in enumerate(map_components.curve_orig_inds):
        
        source_components.curve_comp_inds[cur_ind] = n_powers + cur_ind
        
        source_components.curve_ref_freqs[cur_ind] = REF_FREQ
        source_components.curve_ref_stokesI[cur_ind] = main_table['NORM_COMP_CPL'][old_comp_ind]
        source_components.curve_ref_stokesQ[cur_ind] = 0.0
        source_components.curve_ref_stokesU[cur_ind] = 0.0
        source_components.curve_ref_stokesV[cur_ind] = 0.0
        source_components.curve_SIs[cur_ind] = main_table['ALPHA_CPL'][old_comp_ind]
        source_components.curve_qs[cur_ind] = main_table['CURVE_CPL'][old_comp_ind]
        
        
    ##Get all the flux column names/freqs and put em in a list
    flux_col_freqs = []
    flux_col_names = []
    for key in main_table.columns:
        if key[:7] == 'INT_FLX':
            flux_col_names.append(key)
            flux_col_freqs.append(float(key[7:])*1e+6)
    flux_col_freqs = np.array(flux_col_freqs)
    
    ##Read all flux info for the current chunk into an array so we can do
    ##faster array manipulation on it
    
    all_list_fluxes = np.zeros((n_lists, len(flux_col_names)), dtype=np.float64)
    
    for flux_ind, flux_col in enumerate(flux_col_names):
        all_list_fluxes[:, flux_ind] = main_table[flux_col][map_components.list_orig_inds]
        
    # ##TODO vectorise this please, is probs sloooow
    list_start_index = 0
    for list_ind, old_comp_ind in enumerate(map_components.list_orig_inds):
        
        source_components.list_comp_inds[list_ind] = n_powers + n_curves + list_ind
        
        ##Gather all frequency column information for this component
        these_fluxes = all_list_fluxes[list_ind, :]
        
        ##Figure out what isn't a nan to grab actual values
        use_fluxes = np.where(np.isnan(these_fluxes) == False)
        
        these_fluxes = these_fluxes[use_fluxes]
        these_freqs = flux_col_freqs[use_fluxes]
        
        source_components.num_list_values[list_ind] = int(len(these_freqs))
        source_components.list_start_indexes[list_ind] = list_start_index
        
        find = 0
        for freq, flux in zip(these_freqs, these_fluxes):
            source_components.list_freqs[list_start_index + find] = freq
            source_components.list_stokesI[list_start_index + find] = flux
            source_components.list_stokesQ[list_start_index + find] = 0.0
            source_components.list_stokesU[list_start_index + find] = 0.0
            source_components.list_stokesV[list_start_index + find] = 0.0
            find += 1
        
        list_start_index += int(len(these_freqs))
        
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

            
            comp_ind = np.where(np.array(main_table['NAME'][comp_orig_inds], dtype=str) == shape_name)[0]
            source_components.param_indexes[new_b_ind] = comp_ind
    
    return


# @profile
def read_fits_skymodel_chunks(fits_path : str,
                              chunked_skymodel_maps : list,
                              num_freqs : int, num_time_steps : int,
                              beamtype : int,
                              lsts : np.ndarray, latitude : float,
                              precision = "double") -> Union[Source_Catalogue_Float, Source_Catalogue_Double]:
    
    ##want to know how many shapelets are in all the chunks (used later
    # by "calculate_visiblities.cu")
    num_shapelets = 0
    for chunk_map in chunked_skymodel_maps:
        num_shapelets += chunk_map.n_shapes
        
    ##setup the source catalogue, which is going to store all of the information
    ##of each source and be fed straight into C/CUDA
    source_catalogue = setup_source_catalogue(len(chunked_skymodel_maps), num_shapelets,
                                precision = precision)
    
    
    main_table = Table.read(fits_path, hdu=1)
    
    ##Only read in a shape table if there are shapelets
    if chunk_map.n_shapes:
        shape_table = Table.read(fits_path, hdu=2)
    else:
        shape_table = 'no_shape_table'
    
    ##for each chunk map, create a Source_Float or Source_Double ctype
    ##struct, and "malloc" the right amount of arrays to store required infor
    for chunk_ind, chunk_map in enumerate(chunked_skymodel_maps):
        source_catalogue.sources[chunk_ind] = setup_chunked_source(chunk_map,
                                                num_freqs, num_time_steps,
                                                beamtype, precision=precision)
        
        chunk_source = source_catalogue.sources[chunk_ind]
        
        # print('CHUNK IND', chunk_ind, "--------------------------")
        
        ##count up the total number of components across all chunks
        ##annoyingly, beacuse Jack sucks, we split shapelet us by basis 
        ##and not component, so this number is actually a combination of
        ##component and basis numbers
        # num_comps_all_chunks += chunk_map.n_points + chunk_map.n_gauss + chunk_map.n_shape_coeffs
        
        if chunk_map.n_points > 0:
            add_fits_info_to_source_catalogue(CompTypes.POINT,
                                      main_table, shape_table,
                                      chunk_source, chunk_map,
                                      num_freqs, num_time_steps,
                                      beamtype, lsts, latitude,
                                      precision=precision)
            
        if chunk_map.n_gauss > 0:
            add_fits_info_to_source_catalogue(CompTypes.GAUSSIAN,
                                      main_table, shape_table,
                                      chunk_source, chunk_map,
                                      num_freqs, num_time_steps,
                                      beamtype, lsts, latitude,
                                      precision=precision)
            
        if chunk_map.n_shapes > 0:
            add_fits_info_to_source_catalogue(CompTypes.SHAPELET,
                                      main_table, shape_table,
                                      chunk_source, chunk_map,
                                      num_freqs, num_time_steps,
                                      beamtype, lsts, latitude,
                                      precision=precision)
        

    ##TODO some kind of consistency check between the chunk_maps and the
    ##sources in the catalogue - make sure we read in the correct information
    
    return source_catalogue

