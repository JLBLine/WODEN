##Set some directories we need the include files from
include_directories("${CMAKE_SOURCE_DIR}/cmake_testing/calculate_visibilities"
                     ${UNITY_DIR})


## Create tests for both float and double precision
foreach(PRECISION IN LISTS FLOAT DOUBLE)
# foreach(PRECISION IN LISTS DOUBLE)
# foreach(PRECISION IN LISTS FLOAT)
    ##Set the flags depending on the precision
    DEFINE_COMP_FLAGS(${PRECISION} C_FLAGS CUDA_FLAGS C_COVER_FLAGS)

    if(USE_CUDA)
        set_source_files_properties( "${CMAKE_SOURCE_DIR}/src/calculate_visibilities.cpp"
              "${CMAKE_SOURCE_DIR}/src/source_components.cpp"
              "${CMAKE_SOURCE_DIR}/src/fundamental_coords.cpp"
              "${CMAKE_SOURCE_DIR}/src/primary_beam_gpu.cpp" PROPERTIES LANGUAGE CUDA)
              
    endif()

    ##First add necessary CUDA code to a library
    add_library(calculate_visibilitiesCUDA_${PRECISION} SHARED
                "${CMAKE_SOURCE_DIR}/src/calculate_visibilities.cpp"
                "${CMAKE_SOURCE_DIR}/src/source_components.cpp"
                "${CMAKE_SOURCE_DIR}/src/fundamental_coords.cpp"
                "${CMAKE_SOURCE_DIR}/src/primary_beam_gpu.cpp")

    target_compile_options(calculate_visibilitiesCUDA_${PRECISION} PRIVATE ${CUDA_FLAGS})

    ##Quite a large number of C depenecies here so stick em in a library
    add_library(calculate_visibilities_Cdepends_${PRECISION} SHARED
                ${CMAKE_SOURCE_DIR}/src/visibility_set.c
                ${CMAKE_SOURCE_DIR}/src/shapelet_basis.c)
    target_compile_options(calculate_visibilities_Cdepends_${PRECISION} PRIVATE ${C_FLAGS})

    ##----------------------------------------------------------------------------
    ##Now start adding the tests, linking the necessary libraries and compiling
    ##with the correct precision
    ##----------------------------------------------------------------------------

    # ##Test the main calculate_visibilities function using no primary beam
    add_executable(test_calculate_visibilities_nobeam_${PRECISION}_app
        test_calculate_visibilities_nobeam.c
        test_calculate_visibilities_common.c
    )
    target_link_libraries(test_calculate_visibilities_nobeam_${PRECISION}_app
        calculate_visibilities_Cdepends_${PRECISION}
        calculate_visibilitiesCUDA_${PRECISION}
        ${CC_LINKLIBS}
        Unity
    )
    target_compile_options(test_calculate_visibilities_nobeam_${PRECISION}_app
                                                                PRIVATE ${C_FLAGS})
    add_test(CUDA_test_calculate_visibilities_nobeam_${PRECISION}
            test_calculate_visibilities_nobeam_${PRECISION}_app)



    ##Test the main calculate_visibilities function using the gaussian primary beam
    add_executable(test_calculate_visibilities_gaussbeam_${PRECISION}_app
        test_calculate_visibilities_gaussbeam.c
        test_calculate_visibilities_common.c
    )
    target_link_libraries(test_calculate_visibilities_gaussbeam_${PRECISION}_app
        calculate_visibilities_Cdepends_${PRECISION}
        calculate_visibilitiesCUDA_${PRECISION}
        ${CC_LINKLIBS}
        Unity
    )
    target_compile_options(test_calculate_visibilities_gaussbeam_${PRECISION}_app
                                                            PRIVATE ${C_FLAGS})
    add_test(CUDA_test_calculate_visibilities_gaussbeam_${PRECISION}
            test_calculate_visibilities_gaussbeam_${PRECISION}_app)



    ##Test the main calculate_visibilities function using the EDA2 primary beam
    add_executable(test_calculate_visibilities_edabeam_${PRECISION}_app
        test_calculate_visibilities_edabeam.c
        test_calculate_visibilities_common.c
    )
    target_link_libraries(test_calculate_visibilities_edabeam_${PRECISION}_app
        calculate_visibilities_Cdepends_${PRECISION}
        calculate_visibilitiesCUDA_${PRECISION}
        ${CC_LINKLIBS}
        Unity
    )
    target_compile_options(test_calculate_visibilities_edabeam_${PRECISION}_app
                                                            PRIVATE ${C_FLAGS})
    add_test(CUDA_test_calculate_visibilities_edabeam_${PRECISION}
            test_calculate_visibilities_edabeam_${PRECISION}_app)



    ##Test the main calculate_visibilities function using the mwa fee primary beam
    add_executable(test_calculate_visibilities_mwafeebeam_${PRECISION}_app
        test_calculate_visibilities_mwafeebeam.c
        test_calculate_visibilities_common.c
    )
    target_link_libraries(test_calculate_visibilities_mwafeebeam_${PRECISION}_app
        calculate_visibilities_Cdepends_${PRECISION}
        calculate_visibilitiesCUDA_${PRECISION}
        ${CC_LINKLIBS}
        Unity
    )
    target_compile_options(test_calculate_visibilities_mwafeebeam_${PRECISION}_app
                                                            PRIVATE ${C_FLAGS})
    add_test(CUDA_test_calculate_visibilities_mwafeebeam_${PRECISION}
            test_calculate_visibilities_mwafeebeam_${PRECISION}_app)


    # ##Test the main calculate_visibilities function using the mwa fee primary beam
    add_executable(test_calculate_visibilities_mwafeebeaminterp_${PRECISION}_app
        test_calculate_visibilities_mwafeebeaminterp.c
        test_calculate_visibilities_common.c
    )
    target_link_libraries(test_calculate_visibilities_mwafeebeaminterp_${PRECISION}_app
        calculate_visibilities_Cdepends_${PRECISION}
        calculate_visibilitiesCUDA_${PRECISION}
        ${CC_LINKLIBS}
        Unity
    )
    target_compile_options(test_calculate_visibilities_mwafeebeaminterp_${PRECISION}_app
                                                            PRIVATE ${C_FLAGS})
    add_test(CUDA_test_calculate_visibilities_mwafeebeaminterp_${PRECISION}
    test_calculate_visibilities_mwafeebeaminterp_${PRECISION}_app)


    # ##Test the main calculate_visibilities function using the mwa anaytic beam
    add_executable(test_calculate_visibilities_mwaanalybeam_${PRECISION}_app
        test_calculate_visibilities_mwaanalybeam.c
        test_calculate_visibilities_common.c
    )
    target_link_libraries(test_calculate_visibilities_mwaanalybeam_${PRECISION}_app
        calculate_visibilities_Cdepends_${PRECISION}
        calculate_visibilitiesCUDA_${PRECISION}
        ${CC_LINKLIBS}
        Unity
    )
    target_compile_options(test_calculate_visibilities_mwaanalybeam_${PRECISION}_app
                                                            PRIVATE ${C_FLAGS})
    add_test(CUDA_test_calculate_visibilities_mwaanalybeam_${PRECISION}
    test_calculate_visibilities_mwaanalybeam_${PRECISION}_app)


    # ##Test the main calculate_visibilities function using the mwa anaytic beam
    add_executable(test_calculate_visibilities_multibeams_${PRECISION}_app
    test_calculate_visibilities_multibeams.c
    test_calculate_visibilities_common.c
    )
    target_link_libraries(test_calculate_visibilities_multibeams_${PRECISION}_app
    calculate_visibilities_Cdepends_${PRECISION}
    calculate_visibilitiesCUDA_${PRECISION}
    ${CC_LINKLIBS}
    Unity
    )
    target_compile_options(test_calculate_visibilities_multibeams_${PRECISION}_app
                                                        PRIVATE ${C_FLAGS})
    add_test(CUDA_test_calculate_visibilities_multibeams_${PRECISION}
    test_calculate_visibilities_multibeams_${PRECISION}_app)

endforeach()

# Emergency code for testing ctest environment
# add_test(run_xterm xterm)
