
include(FindPython)

##Set an environment variable for the ctest of this, so we can find a path
##to some test files that live in this dir
add_test(python_test_add_instrumental_effects_woden ${PYTHON_EXECUTABLE}
${CMAKE_CURRENT_SOURCE_DIR}/test_add_instrumental_effects_woden.py)