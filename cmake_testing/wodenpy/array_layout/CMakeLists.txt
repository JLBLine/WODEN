include(FindPythonInterp)

# add_test(python_test_enh2xyz ${PYTHON_EXECUTABLE}
#          ${CMAKE_CURRENT_SOURCE_DIR}/test_enh2xyz.py)

add_test(python_test_RTS_precess ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_RTS_precess.py)
set_tests_properties(python_test_RTS_precess PROPERTIES
            ENVIRONMENT CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR})

add_test(python_test_RTS_PrecessXYZtoJ2000 ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_RTS_PrecessXYZtoJ2000.py)
set_tests_properties(python_test_RTS_PrecessXYZtoJ2000 PROPERTIES
            ENVIRONMENT CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR})

# add_test(python_test_calc_XYZ_diffs ${PYTHON_EXECUTABLE}
#          ${CMAKE_CURRENT_SOURCE_DIR}/test_calc_XYZ_diffs.py)
# set_tests_properties(python_test_calc_XYZ_diffs PROPERTIES
#             ENVIRONMENT CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR})