add_test(python_test_calc_jdcal ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_calc_jdcal.py)

add_test(python_test_get_uvfits_date_and_position_constants ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_get_uvfits_date_and_position_constants.py)