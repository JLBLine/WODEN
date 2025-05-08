add_test(python_test_run_woden ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_run_woden.py)

DEFINE_COMP_FLAGS(DOUBLE C_FLAGS CUDA_FLAGS C_COVER_FLAGS)

##Create a combined sky model reading library that covers both the text and
##yaml formats, given the code expected_skymodel_outcomes.c used below
##relies on both
add_library(do_simple_functions SHARED do_simple_functions.c)
target_link_libraries(do_simple_functions)

target_include_directories(do_simple_functions PUBLIC
            "${CMAKE_SOURCE_DIR}/include")
target_compile_options(do_simple_functions PRIVATE
                        ${C_FLAGS})

add_test(python_test_get_future_result ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_get_future_result.py)