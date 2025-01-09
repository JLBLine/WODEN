## Make a library to read in the woden_settings inputs, depending on whether
## things are float or double precision
foreach(PRECISION IN LISTS FLOAT DOUBLE)

  ##Set the flags depending on the precision
  DEFINE_COMP_FLAGS(${PRECISION} C_FLAGS CUDA_FLAGS C_COVER_FLAGS)

  ##Create a combined sky model reading library that covers both the text and
  ##yaml formats, given the code expected_skymodel_outcomes.c used below
  ##relies on both
  add_library(read_woden_settings_${PRECISION} SHARED read_woden_settings.c)
  target_link_libraries(read_woden_settings_${PRECISION})

  target_include_directories(read_woden_settings_${PRECISION} PUBLIC
                "${CMAKE_SOURCE_DIR}/cmake_testing/wodenpy/use_libwoden/"
                "${CMAKE_SOURCE_DIR}/include")
  target_compile_options(read_woden_settings_${PRECISION} PRIVATE
                          ${C_FLAGS})


  add_library(read_source_catalogue_${PRECISION} SHARED read_source_catalogue.c)
  target_link_libraries(read_source_catalogue_${PRECISION})

  target_include_directories(read_source_catalogue_${PRECISION} PUBLIC
                "${CMAKE_SOURCE_DIR}/cmake_testing/wodenpy/use_libwoden/"
                "${CMAKE_SOURCE_DIR}/include")
  target_compile_options(read_source_catalogue_${PRECISION} PRIVATE
                          ${C_FLAGS})

endforeach()

##Set an environment variable for the ctest of this, so we can find a path
##to some test files that live in this dir
add_test(python_test_fill_woden_settings_python ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_fill_woden_settings_python.py)
set_tests_properties(python_test_fill_woden_settings_python PROPERTIES
         ENVIRONMENT CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR})

add_test(python_test_convert_woden_settings_to_ctypes ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_convert_woden_settings_to_ctypes.py)
set_tests_properties(python_test_convert_woden_settings_to_ctypes PROPERTIES
         ENVIRONMENT CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR})

##Set an environment variable for the ctest of this, so we can find a path
##to some test files that live in this dir
add_test(python_test_setup_lsts_and_phase_centre ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_setup_lsts_and_phase_centre.py)
set_tests_properties(python_test_setup_lsts_and_phase_centre PROPERTIES
         ENVIRONMENT CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR})

add_test(python_test_create_sbf ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_create_sbf.py)
set_tests_properties(python_test_create_sbf PROPERTIES
         ENVIRONMENT CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR})

##Set an environment variable for the ctest of this, so we can find a path
##to some test files that live in this dir
add_test(python_test_make_skymodel_structs ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_make_skymodel_structs.py)
set_tests_properties(python_test_make_skymodel_structs PROPERTIES
         ENVIRONMENT CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR})