DEFINE_COMP_FLAGS(DOUBLE C_FLAGS CUDA_FLAGS C_COVER_FLAGS)

##Create a combined sky model reading library that covers both the text and
##yaml formats, given the code expected_skymodel_outcomes.c used below
##relies on both
add_library(c_logger_funcs SHARED test_c_logging.c
                                  ${CMAKE_SOURCE_DIR}/src/logger.c)
target_link_libraries(c_logger_funcs)

target_include_directories(c_logger_funcs PUBLIC
            "${CMAKE_SOURCE_DIR}/include")
target_compile_options(c_logger_funcs PRIVATE
                        ${C_FLAGS})

##Set an environment variable for the ctest of this, so we can find a path
##to some test files that live in this dir
add_test(python_test_argument_inputs ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_argument_inputs.py)

add_test(python_test_get_code_version ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_get_code_version.py)

add_test(python_test_make_logger ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_make_logger.py)

add_test(python_test_log_chosen_beamtype ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_log_chosen_beamtype.py)