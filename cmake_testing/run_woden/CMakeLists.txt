
include(FindPythonInterp)

# add_test(test_command ${PYTHON_EXECUTABLE}
#          ${CMAKE_CURRENT_SOURCE_DIR}/test_command.py)

# add_test(test_calc_jdcal ${PYTHON_EXECUTABLE}
#          ${CMAKE_CURRENT_SOURCE_DIR}/test_calc_jdcal.py)

# add_test(test_get_LST ${PYTHON_EXECUTABLE}
#          ${CMAKE_CURRENT_SOURCE_DIR}/test_get_LST.py)

# add_test(test_RTS_encoding ${PYTHON_EXECUTABLE}
#          ${CMAKE_CURRENT_SOURCE_DIR}/test_RTS_encoding.py)

add_test(test_make_antenna_table ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_make_antenna_table.py)




# add_test(test_calc_jdcal ${PYTHON_EXECUTABLE}
#          ${CMAKE_CURRENT_SOURCE_DIR}/test_calc_jdcal.py)
