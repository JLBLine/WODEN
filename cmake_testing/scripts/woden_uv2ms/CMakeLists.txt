
include(FindPythonInterp)

##Set an environment variable for the ctest of this, so we can find a path
##to some test files that live in this dir
add_test(python_test_woden_uv2ms ${PYTHON_EXECUTABLE}
${CMAKE_CURRENT_SOURCE_DIR}/test_woden_uv2ms.py)