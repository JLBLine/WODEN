
include(FindPythonInterp)

add_test(python_test_command ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_command.py)

add_test(python_test_calc_jdcal ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_calc_jdcal.py)

add_test(python_test_get_uvfits_date_and_position_constants ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_get_uvfits_date_and_position_constants.py)

add_test(python_test_RTS_encoding ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_RTS_encoding.py)

add_test(python_test_make_antenna_table ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_make_antenna_table.py)

add_test(python_test_create_uvfits ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_create_uvfits.py)

add_test(python_test_enh2xyz ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_enh2xyz.py)

add_test(python_test_load_data ${PYTHON_EXECUTABLE}
${CMAKE_CURRENT_SOURCE_DIR}/test_load_data.py)

add_test(python_test_write_json ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_write_json.py)

add_test(python_test_make_baseline_date_arrays ${PYTHON_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/test_make_baseline_date_arrays.py)

add_test(python_test_remove_phase_tracking ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_remove_phase_tracking.py)

##Set an environment variable for the ctest of this, so we can find a path
##to some test files that live in this dir
add_test(python_test_argument_inputs ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_argument_inputs.py)
set_tests_properties(python_test_argument_inputs PROPERTIES
         ENVIRONMENT CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR})

#
add_test(python_test_read_uvfits_into_pyuvdata ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_read_uvfits_into_pyuvdata.py)
