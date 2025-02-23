##Set some directories we need the include files from
include_directories("${CMAKE_SOURCE_DIR}/cmake_testing/calculate_visibilities"
                     ${UNITY_DIR})


## Create tests for both float and double precision
foreach(PRECISION IN LISTS FLOAT DOUBLE)
# foreach(PRECISION IN LISTS DOUBLE)
# foreach(PRECISION IN LISTS FLOAT)
    ##Set the flags depending on the precision
    DEFINE_COMP_FLAGS(${PRECISION} C_FLAGS GPU_FLAGS C_COVER_FLAGS)

    foreach(DEVICE IN LISTS CPU GPU)
    # foreach(DEVICE IN LISTS GPU)
    # foreach(DEVICE IN LISTS CPU)
        # ##Test the main calculate_visibilities function using no primary beam
        add_executable(test_calculate_visibilities_nobeam_${DEVICE}_${PRECISION}_app
            test_calculate_visibilities_nobeam_${DEVICE}.c
            calculate_visibilities_nobeam_common.c
            calculate_visibilities_common_common.c
        )
        target_link_libraries(test_calculate_visibilities_nobeam_${DEVICE}_${PRECISION}_app
            calculate_visibilities_CPU_${PRECISION}
            calculate_visibilities_GPU_${PRECISION}
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
            calculate_visibilities_CPU_${PRECISION}
            calculate_visibilities_GPU_${PRECISION}
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
            calculate_visibilities_CPU_${PRECISION}
            calculate_visibilities_GPU_${PRECISION}
            ${CC_LINKLIBS}
            Unity gcov
        )
        target_compile_options(test_calculate_visibilities_edabeam_${DEVICE}_${PRECISION}_app
                                                                    PRIVATE ${C_FLAGS})
        add_test(${DEVICE}_test_calculate_visibilities_edabeam_${PRECISION}
                test_calculate_visibilities_edabeam_${DEVICE}_${PRECISION}_app)

        add_executable(test_calculate_visibilities_mwafeebeam_${DEVICE}_${PRECISION}_app
            test_calculate_visibilities_mwafeebeam_${DEVICE}.c
            calculate_visibilities_mwafeebeam_common.c
            calculate_visibilities_common_common.c
        )
        target_link_libraries(test_calculate_visibilities_mwafeebeam_${DEVICE}_${PRECISION}_app
            calculate_visibilities_CPU_${PRECISION}
            calculate_visibilities_GPU_${PRECISION}
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
            calculate_visibilities_CPU_${PRECISION}
            calculate_visibilities_GPU_${PRECISION}
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
            calculate_visibilities_CPU_${PRECISION}
            calculate_visibilities_GPU_${PRECISION}
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
            calculate_visibilities_CPU_${PRECISION}
            calculate_visibilities_GPU_${PRECISION}
            ${CC_LINKLIBS}
            Unity gcov
        )
        target_compile_options(test_calculate_visibilities_multibeams_${DEVICE}_${PRECISION}_app
                                                                    PRIVATE ${C_FLAGS})
        add_test(${DEVICE}_test_calculate_visibilities_multibeams_${PRECISION}
                test_calculate_visibilities_multibeams_${DEVICE}_${PRECISION}_app)

        add_executable(test_calculate_visibilities_everybeam_mwa_${DEVICE}_${PRECISION}_app
            test_calculate_visibilities_everybeam_mwa_${DEVICE}.c
            calculate_visibilities_everybeam_common.c
            calculate_visibilities_common_common.c
        )
        target_link_libraries(test_calculate_visibilities_everybeam_mwa_${DEVICE}_${PRECISION}_app
            calculate_visibilities_CPU_${PRECISION}
            calculate_visibilities_GPU_${PRECISION}
            ${CC_LINKLIBS}
            Unity gcov
        )
        target_compile_options(test_calculate_visibilities_everybeam_mwa_${DEVICE}_${PRECISION}_app
                                                                    PRIVATE ${C_FLAGS})
        add_test(${DEVICE}_test_calculate_visibilities_everybeam_mwa_${PRECISION}
                test_calculate_visibilities_everybeam_mwa_${DEVICE}_${PRECISION}_app)

        # add_executable(test_calculate_visibilities_everybeam_lofar_${DEVICE}_${PRECISION}_app
        #     test_calculate_visibilities_everybeam_lofar_${DEVICE}.c
        #     calculate_visibilities_everybeam_common.c
        #     calculate_visibilities_common_common.c
        # )
        # target_link_libraries(test_calculate_visibilities_everybeam_lofar_${DEVICE}_${PRECISION}_app
        #     calculate_visibilities_CPU_${PRECISION}
        #     calculate_visibilities_GPU_${PRECISION}
        #     ${CC_LINKLIBS}
        #     Unity gcov
        # )
        # target_compile_options(test_calculate_visibilities_everybeam_lofar_${DEVICE}_${PRECISION}_app
        #                                                             PRIVATE ${C_FLAGS})
        # add_test(${DEVICE}_test_calculate_visibilities_everybeam_lofar_${PRECISION}
        #         test_calculate_visibilities_everybeam_lofar_${DEVICE}_${PRECISION}_app)

        # add_executable(test_calculate_visibilities_everybeam_oskar_${DEVICE}_${PRECISION}_app
        #     test_calculate_visibilities_everybeam_oskar_${DEVICE}.c
        #     calculate_visibilities_everybeam_common.c
        #     calculate_visibilities_common_common.c
        # )
        # target_link_libraries(test_calculate_visibilities_everybeam_oskar_${DEVICE}_${PRECISION}_app
        #     calculate_visibilities_CPU_${PRECISION}
        #     calculate_visibilities_GPU_${PRECISION}
        #     ${CC_LINKLIBS}
        #     Unity gcov
        # )
        # target_compile_options(test_calculate_visibilities_everybeam_oskar_${DEVICE}_${PRECISION}_app
        #                                                             PRIVATE ${C_FLAGS})
        # add_test(${DEVICE}_test_calculate_visibilities_everybeam_oskar_${PRECISION}
        #         test_calculate_visibilities_everybeam_oskar_${DEVICE}_${PRECISION}_app)

    endforeach()
    ##Add an executable to profile the code when calling hyperbeam, but
    ##don't add the test as we'll design it to run for a while
    # Test hyperbeam with the coarse hdf5 file
    add_executable(profile_calc_visi_cpu_feebeam_${PRECISION}_app
        profile_calc_visi_cpu_feebeam.c
        # calculate_visibilities_mwafeebeam_common.c
        # calculate_visibilities_common_common.c
    )
    target_link_libraries(profile_calc_visi_cpu_feebeam_${PRECISION}_app
        calculate_visibilities_CPU_${PRECISION}
        calculate_visibilities_GPU_${PRECISION}
        ${CC_LINKLIBS}
        Unity gcov
    )
    target_compile_options(profile_calc_visi_cpu_feebeam_${PRECISION}_app
                                                    PRIVATE ${C_FLAGS})
endforeach()


        

# Emergency code for testing ctest environment
# add_test(run_xterm xterm)
