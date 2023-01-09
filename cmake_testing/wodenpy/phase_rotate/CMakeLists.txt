include(FindPythonInterp)

add_test(python_test_remove_phase_tracking ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_remove_phase_tracking.py)