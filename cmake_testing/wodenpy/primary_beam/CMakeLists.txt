if (USE_EVERYBEAM)
    add_test(python_test_run_everybeam_over_threads ${PYTHON_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_run_everybeam_over_threads.py)

    add_test(python_test_run_everybeam_over_threads_MWA ${PYTHON_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_run_everybeam_over_threads_MWA.py)

    add_test(python_test_check_ms_telescope_type_matches_element_response ${PYTHON_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_check_ms_telescope_type_matches_element_response.py)

    add_test(python_test_create_filtered_ms ${PYTHON_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_create_filtered_ms.py)

    add_test(python_test_run_everybeam_OSKAR ${PYTHON_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_run_everybeam_OSKAR.py)
            
endif()



add_test(python_test_run_uvbeam_MWA ${PYTHON_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/test_run_uvbeam_MWA.py)

add_test(python_test_calc_uvbeam_for_components ${PYTHON_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/test_calc_uvbeam_for_components.py)

add_test(python_test_run_uvbeam_HERA_CST ${PYTHON_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/test_run_uvbeam_HERA_CST.py)

add_test(python_test_run_uvbeam_HERA_FITS ${PYTHON_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/test_run_uvbeam_HERA_FITS.py)
