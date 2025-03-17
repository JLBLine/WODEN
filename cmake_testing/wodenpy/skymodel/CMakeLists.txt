add_test(python_test_crop_below_horizon ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_crop_below_horizon.py)
set_tests_properties(python_test_crop_below_horizon PROPERTIES
         ENVIRONMENT CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR})

add_test(python_test_map_chunk_pointgauss ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_map_chunk_pointgauss.py)
set_tests_properties(python_test_map_chunk_pointgauss PROPERTIES
         ENVIRONMENT CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR})

add_test(python_test_map_chunk_shapelets ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_map_chunk_shapelets.py)
set_tests_properties(python_test_map_chunk_shapelets PROPERTIES
         ENVIRONMENT CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR})

add_test(python_test_create_skymodel_chunk_map ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_create_skymodel_chunk_map.py)
set_tests_properties(python_test_create_skymodel_chunk_map PROPERTIES
         ENVIRONMENT CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR})

add_test(python_test_read_yaml_skymodel_chunk ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_read_yaml_skymodel_chunk.py)
set_tests_properties(python_test_read_yaml_skymodel_chunk PROPERTIES
         ENVIRONMENT CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR})

add_test(python_test_read_FITS_skymodel_chunk ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_read_FITS_skymodel_chunk.py)
set_tests_properties(python_test_read_FITS_skymodel_chunk PROPERTIES
         ENVIRONMENT CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR})

add_test(python_test_read_text_skymodel_chunk ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_read_text_skymodel_chunk.py)
set_tests_properties(python_test_read_text_skymodel_chunk PROPERTIES
         ENVIRONMENT CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR})

add_test(python_test_read_skymodel_chunk ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_read_skymodel_chunk.py)
set_tests_properties(python_test_read_skymodel_chunk PROPERTIES
         ENVIRONMENT CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR})

add_test(python_test_check_columns_fits ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_check_columns_fits.py)