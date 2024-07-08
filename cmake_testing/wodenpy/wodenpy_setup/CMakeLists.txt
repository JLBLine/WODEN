
include(FindPython)

##Set an environment variable for the ctest of this, so we can find a path
##to some test files that live in this dir
add_test(python_test_argument_inputs ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_argument_inputs.py)

add_test(python_test_get_code_version ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_get_code_version.py)