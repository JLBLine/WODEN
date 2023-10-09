##Delete any old reports hanging around
rm coverage.xml *.gcov

##Use the ctest C code outputs and run gcov over them for codecov
gcov ../build/cmake_testing/C_code/primary_beam/CMakeFiles/primary_beam_float.dir/__/__/__/src/primary_beam.c.gcda \
    ../build/cmake_testing/C_code/visibility_set/CMakeFiles/visibility_set_float.dir/__/__/__/src/visibility_set.c.gcda

##Things are coded up for ctest environment, so stick a relative path in here
export CMAKE_CURRENT_SOURCE_DIR=../cmake_testing/wodenpy/array_layout

##Run the python tests using python 'coverage'
##Can be grabbed with 'pip install coverage'
coverage run --source=create_array_layout ../cmake_testing/wodenpy/array_layout/test_calc_XYZ_diffs.py
coverage run --source=create_array_layout -a ../cmake_testing/wodenpy/array_layout/test_enh2xyz.py
coverage run --source=create_array_layout -a ../cmake_testing/wodenpy/array_layout/test_RTS_PrecessXYZtoJ2000.py
coverage run --source=precession ../cmake_testing/wodenpy/array_layout/test_RTS_precess.py

coverage run --source=calc_obs ../cmake_testing/wodenpy/observational/test_calc_jdcal.py
coverage run --source=calc_obs -a ../cmake_testing/wodenpy/observational/test_get_uvfits_date_and_position_constants.py

export CMAKE_CURRENT_SOURCE_DIR=../cmake_testing/wodenpy/skymodel

coverage run --source=chunk_sky_model ../cmake_testing/wodenpy/skymodel/test_create_skymodel_chunk_map.py
coverage run --source=chunk_sky_model -a ../cmake_testing/wodenpy/skymodel/test_map_chunk_pointgauss.py
coverage run --source=chunk_sky_model -a ../cmake_testing/wodenpy/skymodel/test_map_chunk_shapelets.py

coverage run --source=woden_skymodel ../cmake_testing/wodenpy/skymodel/test_crop_below_horizon.py

coverage run --source=read_skymodel ../cmake_testing/wodenpy/skymodel/test_read_skymodel_chunk.py
coverage run --source=read_fits_skymodel ../cmake_testing/wodenpy/skymodel/test_read_skymodel_chunk.py
coverage run --source=read_text_skymodel ../cmake_testing/wodenpy/skymodel/test_read_skymodel_chunk.py
coverage run --source=read_yaml_skymodel ../cmake_testing/wodenpy/skymodel/test_read_skymodel_chunk.py


coverage run --source=read_fits_skymodel -a ../cmake_testing/wodenpy/skymodel/test_read_FITS_skymodel_chunk.py
coverage run --source=read_text_skymodel -a ../cmake_testing/wodenpy/skymodel/test_read_text_skymodel_chunk.py
coverage run --source=read_yaml_skymodel -a ../cmake_testing/wodenpy/skymodel/test_read_yaml_skymodel_chunk.py


#coverage run --source=create_array_layout -a ../cmake_testing/wodenpy/array_layout/test_RTS_encoding.py
# coverage run --source=run_woden -a ../cmake_testing/run_woden/test_command.py
# coverage run --source=run_woden -a ../cmake_testing/run_woden/test_create_uvfits.py

# 
# coverage run --source=run_woden -a ../cmake_testing/run_woden/test_load_data.py
# coverage run --source=run_woden -a ../cmake_testing/run_woden/test_make_antenna_table.py
# coverage run --source=run_woden -a ../cmake_testing/run_woden/test_make_baseline_date_arrays.py
# coverage run --source=run_woden -a ../cmake_testing/run_woden/test_remove_phase_tracking.py
# coverage run --source=run_woden -a ../cmake_testing/run_woden/test_write_json.py

# coverage run --source=add_woden_uvfits -a ../cmake_testing/add_woden_uvfits/test_add_woden_uvfits.py
# coverage run --source=concat_woden_uvfits -a ../cmake_testing/concat_woden_uvfits/test_concat_woden_uvfits.py

##convert output to something that codecov accepts
coverage xml

##delete things that were written out by running tests
rm WODEN_array_layout.txt *.uvfits *.json \
    example.txt test_load_data.dat test_full_skymodel*
