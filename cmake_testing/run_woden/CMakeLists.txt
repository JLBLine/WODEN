
include(FindPythonInterp)

add_test(python_test_command ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_command.py)

add_test(python_test_calc_jdcal ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_calc_jdcal.py)

add_test(python_test_get_LST ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_get_LST.py)

add_test(python_test_RTS_encoding ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_RTS_encoding.py)

add_test(python_test_make_antenna_table ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_make_antenna_table.py)

add_test(python_test_create_uvfits ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_create_uvfits.py)

add_test(python_test_write_json ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_write_json.py)

add_test(python_test_make_baseline_date_arrays ${PYTHON_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/test_make_baseline_date_arrays.py)

# add_test(python_test_remove_phase_tracking ${PYTHON_EXECUTABLE}
#          ${CMAKE_CURRENT_SOURCE_DIR}/test_remove_phase_tracking.py)

# add_test(python_test_calc_jdcal ${PYTHON_EXECUTABLE}
#          ${CMAKE_CURRENT_SOURCE_DIR}/test_calc_jdcal.py)
