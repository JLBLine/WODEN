
include(FindPythonInterp)

##Set an environment variable for the ctest of this, so we can find a path
##to some test files that live in this dir
add_test(python_test_argument_inputs ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_argument_inputs.py)
set_tests_properties(python_test_argument_inputs PROPERTIES
         ENVIRONMENT CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR})