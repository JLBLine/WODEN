##Use the ctest C code outputs and run gcov over them for codecov
gcov ../build/cmake_testing/array_layout/CMakeFiles/array_layout.dir/__/__/src/array_layout.c.gcda \
    ../build/cmake_testing/chunk_sky_model/CMakeFiles/chunk_sky_model.dir/__/__/src/chunk_sky_model.c.gcda \
    ../build/cmake_testing/create_sky_model/CMakeFiles/create_sky_model.dir/__/__/src/create_sky_model.c.gcda \
    ../build/cmake_testing/FEE_primary_beam/CMakeFiles/FEE_primary_beam.dir/__/__/src/FEE_primary_beam.c.gcda \
    ../build/cmake_testing/primary_beam/CMakeFiles/primary_beam.dir/__/__/src/primary_beam.c.gcda \
    ../build/cmake_testing/visibility_set/CMakeFiles/visibility_set.dir/__/__/src/visibility_set.c.gcda \
    ../build/cmake_testing/woden_settings/CMakeFiles/woden_settings.dir/__/__/src/woden_settings.c.gcda \
    ../build/cmake_testing/shapelet_basis/CMakeFiles/shapelet_basis.dir/__/__/src/shapelet_basis.c.gcda

##Things are coded up for ctest environment, so stick a relative path in here
export CMAKE_CURRENT_SOURCE_DIR=../cmake_testing/run_woden/

##Run the python tests using python 'coverage'
##Can be grabbed with 'pip install coverage'
coverage run --source=run_woden ../cmake_testing/run_woden/test_argument_inputs.py
coverage run --source=run_woden -a ../cmake_testing/run_woden/test_calc_jdcal.py
coverage run --source=run_woden -a ../cmake_testing/run_woden/test_command.py
coverage run --source=run_woden -a ../cmake_testing/run_woden/test_create_uvfits.py
coverage run --source=run_woden -a ../cmake_testing/run_woden/test_enh2xyz.py
coverage run --source=run_woden -a ../cmake_testing/run_woden/test_get_uvfits_date_and_position_constants.py
coverage run --source=run_woden -a ../cmake_testing/run_woden/test_load_data.py
coverage run --source=run_woden -a ../cmake_testing/run_woden/test_make_antenna_table.py
coverage run --source=run_woden -a ../cmake_testing/run_woden/test_make_baseline_date_arrays.py
coverage run --source=run_woden -a ../cmake_testing/run_woden/test_remove_phase_tracking.py
coverage run --source=run_woden -a ../cmake_testing/run_woden/test_RTS_encoding.py
coverage run --source=run_woden -a ../cmake_testing/run_woden/test_write_json.py

##convert output to something that codecov accepts
coverage xml

##delete things that were written out by running tests
rm WODEN_array_layout.txt unittest_example.uvfits test_write_gaussian_beam.json \
    test_write_minimum_json.json test_write_mwafee_beam.json example.txt test_load_data.dat
