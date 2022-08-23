
include(FindPythonInterp)

##Set an environment variable for the ctest of this, so we can find a path
##to some test files that live in this dir
add_test(python_test_concat_woden_uvfits ${PYTHON_EXECUTABLE}
${CMAKE_CURRENT_SOURCE_DIR}/test_concat_woden_uvfits.py)
set_tests_properties(python_test_concat_woden_uvfits PROPERTIES
         ENVIRONMENT CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR})
