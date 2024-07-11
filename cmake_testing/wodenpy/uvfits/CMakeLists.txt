add_test(python_test_RTS_encoding ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_RTS_encoding.py)

add_test(python_test_make_antenna_table ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_make_antenna_table.py)

add_test(python_test_create_uvfits ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_create_uvfits.py)

add_test(python_test_make_baseline_date_arrays ${PYTHON_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/test_make_baseline_date_arrays.py)

add_test(python_test_read_uvfits_into_pyuvdata ${PYTHON_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/test_read_uvfits_into_pyuvdata.py)