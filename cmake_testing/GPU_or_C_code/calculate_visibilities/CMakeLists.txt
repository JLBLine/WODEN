##Set some directories we need the include files from
include_directories("${CMAKE_SOURCE_DIR}/cmake_testing/calculate_visibilities"
                     ${UNITY_DIR})


## Create tests for both float and double precision
foreach(PRECISION IN LISTS FLOAT DOUBLE)
# foreach(PRECISION IN LISTS DOUBLE)
# foreach(PRECISION IN LISTS FLOAT)
    ##Set the flags depending on the precision
    DEFINE_COMP_FLAGS(${PRECISION} C_FLAGS GPU_FLAGS C_COVER_FLAGS)

    if(USE_CUDA)
        set_source_files_properties( "${CMAKE_SOURCE_DIR}/src/calculate_visibilities_gpu.cpp"
              "${CMAKE_SOURCE_DIR}/src/source_components_gpu.cpp"
              "${CMAKE_SOURCE_DIR}/src/fundamental_coords_gpu.cpp"
              "${CMAKE_SOURCE_DIR}/src/primary_beam_gpu.cpp" PROPERTIES LANGUAGE CUDA)
              
    endif()

    ##First add necessary CUDA code to a library
    add_library(calculate_visibilitiesGPU_${PRECISION} SHARED
                "${CMAKE_SOURCE_DIR}/src/calculate_visibilities_gpu.cpp"
                "${CMAKE_SOURCE_DIR}/src/source_components_gpu.cpp"
                "${CMAKE_SOURCE_DIR}/src/fundamental_coords_gpu.cpp"
                "${CMAKE_SOURCE_DIR}/src/primary_beam_gpu.cpp")

    target_compile_options(calculate_visibilitiesGPU_${PRECISION} PRIVATE ${GPU_FLAGS})

    ##Quite a large number of C depenecies here so stick em in a library
    add_library(calculate_visibilities_CPUdepends_${PRECISION} SHARED
                ${CMAKE_SOURCE_DIR}/src/calculate_visibilities_common.c
                ${CMAKE_SOURCE_DIR}/src/source_components_common.c
                ${CMAKE_SOURCE_DIR}/src/visibility_set.c
                ${CMAKE_SOURCE_DIR}/src/shapelet_basis.c)
    target_link_libraries(calculate_visibilities_CPUdepends_${PRECISION} gcov)
    target_compile_options(calculate_visibilities_CPUdepends_${PRECISION} PRIVATE ${C_COVER_FLAGS})

    # foreach(DEVICE IN LISTS CPU GPU)
    foreach(DEVICE IN LISTS GPU)
        ##Test the main calculate_visibilities function using no primary beam
        add_executable(test_calculate_visibilities_nobeam_${DEVICE}_${PRECISION}_app
            test_calculate_visibilities_nobeam_${DEVICE}.c
            calculate_visibilities_nobeam_common.c
            calculate_visibilities_common_common.c
        )
        target_link_libraries(test_calculate_visibilities_nobeam_${DEVICE}_${PRECISION}_app
            calculate_visibilities_CPUdepends_${PRECISION}
            calculate_visibilitiesGPU_${PRECISION}
            ${CC_LINKLIBS}
            Unity gcov
        )
        target_compile_options(test_calculate_visibilities_nobeam_${DEVICE}_${PRECISION}_app
                                                                    PRIVATE ${C_FLAGS})
        add_test(${DEVICE}_test_calculate_visibilities_nobeam_${PRECISION}
                test_calculate_visibilities_nobeam_${DEVICE}_${PRECISION}_app)

        # ##Test the main calculate_visibilities function using no primary beam
        add_executable(test_calculate_visibilities_gaussbeam_${DEVICE}_${PRECISION}_app
            test_calculate_visibilities_gaussbeam_${DEVICE}.c
            calculate_visibilities_gaussbeam_common.c
            calculate_visibilities_common_common.c
        )
        target_link_libraries(test_calculate_visibilities_gaussbeam_${DEVICE}_${PRECISION}_app
            calculate_visibilities_CPUdepends_${PRECISION}
            calculate_visibilitiesGPU_${PRECISION}
            ${CC_LINKLIBS}
            Unity gcov
        )
        target_compile_options(test_calculate_visibilities_gaussbeam_${DEVICE}_${PRECISION}_app
                                                                    PRIVATE ${C_FLAGS})
        add_test(${DEVICE}_test_calculate_visibilities_gaussbeam_${PRECISION}
                test_calculate_visibilities_gaussbeam_${DEVICE}_${PRECISION}_app)

        # ##Test the main calculate_visibilities function using no primary beam
        add_executable(test_calculate_visibilities_edabeam_${DEVICE}_${PRECISION}_app
            test_calculate_visibilities_edabeam_${DEVICE}.c
            calculate_visibilities_edabeam_common.c
            calculate_visibilities_common_common.c
        )
        target_link_libraries(test_calculate_visibilities_edabeam_${DEVICE}_${PRECISION}_app
            calculate_visibilities_CPUdepends_${PRECISION}
            calculate_visibilitiesGPU_${PRECISION}
            ${CC_LINKLIBS}
            Unity gcov
        )
        target_compile_options(test_calculate_visibilities_edabeam_${DEVICE}_${PRECISION}_app
                                                                    PRIVATE ${C_FLAGS})
        add_test(${DEVICE}_test_calculate_visibilities_edabeam_${PRECISION}
                test_calculate_visibilities_edabeam_${DEVICE}_${PRECISION}_app)

            ##Test the main calculate_visibilities function using no primary beam
        add_executable(test_calculate_visibilities_mwafeebeam_${DEVICE}_${PRECISION}_app
            test_calculate_visibilities_mwafeebeam_${DEVICE}.c
            calculate_visibilities_mwafeebeam_common.c
            calculate_visibilities_common_common.c
        )
        target_link_libraries(test_calculate_visibilities_mwafeebeam_${DEVICE}_${PRECISION}_app
            calculate_visibilities_CPUdepends_${PRECISION}
            calculate_visibilitiesGPU_${PRECISION}
            ${CC_LINKLIBS}
            Unity gcov
        )
        target_compile_options(test_calculate_visibilities_mwafeebeam_${DEVICE}_${PRECISION}_app
                                                                    PRIVATE ${C_FLAGS})
        add_test(${DEVICE}_test_calculate_visibilities_mwafeebeam_${PRECISION}
                test_calculate_visibilities_mwafeebeam_${DEVICE}_${PRECISION}_app)


        add_executable(test_calculate_visibilities_mwafeeinterpbeam_${DEVICE}_${PRECISION}_app
            test_calculate_visibilities_mwafeeinterpbeam_${DEVICE}.c
            calculate_visibilities_mwafeeinterpbeam_common.c
            calculate_visibilities_common_common.c
        )
        target_link_libraries(test_calculate_visibilities_mwafeeinterpbeam_${DEVICE}_${PRECISION}_app
            calculate_visibilities_CPUdepends_${PRECISION}
            calculate_visibilitiesGPU_${PRECISION}
            ${CC_LINKLIBS}
            Unity gcov
        )
        target_compile_options(test_calculate_visibilities_mwafeeinterpbeam_${DEVICE}_${PRECISION}_app
                                                                    PRIVATE ${C_FLAGS})
        add_test(${DEVICE}_test_calculate_visibilities_mwafeeinterpbeam_${PRECISION}
                test_calculate_visibilities_mwafeeinterpbeam_${DEVICE}_${PRECISION}_app)

        add_executable(test_calculate_visibilities_mwaanalybeam_${DEVICE}_${PRECISION}_app
            test_calculate_visibilities_mwaanalybeam_${DEVICE}.c
            calculate_visibilities_mwaanalybeam_common.c
            calculate_visibilities_common_common.c
        )
        target_link_libraries(test_calculate_visibilities_mwaanalybeam_${DEVICE}_${PRECISION}_app
            calculate_visibilities_CPUdepends_${PRECISION}
            calculate_visibilitiesGPU_${PRECISION}
            ${CC_LINKLIBS}
            Unity gcov
        )
        target_compile_options(test_calculate_visibilities_mwaanalybeam_${DEVICE}_${PRECISION}_app
                                                                    PRIVATE ${C_FLAGS})
        add_test(${DEVICE}_test_calculate_visibilities_mwaanalybeam_${PRECISION}
                test_calculate_visibilities_mwaanalybeam_${DEVICE}_${PRECISION}_app)

        add_executable(test_calculate_visibilities_multibeams_${DEVICE}_${PRECISION}_app
            test_calculate_visibilities_multibeams_${DEVICE}.c
            calculate_visibilities_multibeams_common.c
            calculate_visibilities_common_common.c
        )
        target_link_libraries(test_calculate_visibilities_multibeams_${DEVICE}_${PRECISION}_app
            calculate_visibilities_CPUdepends_${PRECISION}
            calculate_visibilitiesGPU_${PRECISION}
            ${CC_LINKLIBS}
            Unity gcov
        )
        target_compile_options(test_calculate_visibilities_multibeams_${DEVICE}_${PRECISION}_app
                                                                    PRIVATE ${C_FLAGS})
        add_test(${DEVICE}_test_calculate_visibilities_multibeams_${PRECISION}
                test_calculate_visibilities_multibeams_${DEVICE}_${PRECISION}_app)

    endforeach()

endforeach()

# Emergency code for testing ctest environment
# add_test(run_xterm xterm)
