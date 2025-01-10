add_test(python_test_run_everybeam_over_threads ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_run_everybeam_over_threads.py)

add_test(python_test_run_everybeam_over_threads_MWA ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_run_everybeam_over_threads_MWA.py)

add_test(python_test_run_woden_lofar_para_rotate ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_run_woden_lofar_para_rotate.py)

add_test(python_test_run_fast_everybeam ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/test_run_fast_everybeam.py)