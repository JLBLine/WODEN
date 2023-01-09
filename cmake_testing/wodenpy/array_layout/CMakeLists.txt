include(FindPythonInterp)

add_test(python_test_enh2xyz ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_enh2xyz.py)