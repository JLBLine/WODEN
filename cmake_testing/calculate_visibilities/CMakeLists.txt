add_library(calculate_visibilitiesCUDA SHARED
            "${CMAKE_SOURCE_DIR}/src/calculate_visibilities.cu"
            "${CMAKE_SOURCE_DIR}/src/source_components.cu"
            "${CMAKE_SOURCE_DIR}/src/fundamental_coords.cu"
            "${CMAKE_SOURCE_DIR}/src/primary_beam_cuda.cu"
            "${CMAKE_SOURCE_DIR}/src/FEE_primary_beam_cuda.cu")

target_compile_options(calculate_visibilitiesCUDA PRIVATE -g -G --compiler-options -Wall)

## For CMake versions > 3.18, the arch flag settings are handled by
## CUDA_ARCHITECTURES. For now, just switch this warning off - in future,
## replace the CUDA_FLAG lines above with something like
set_property(TARGET calculate_visibilitiesCUDA PROPERTY CUDA_ARCHITECTURES OFF)

include_directories("${CMAKE_SOURCE_DIR}/cmake_testing/calculate_visibilities")

# ##Test the main calculate_visibilities function using no primary beam
add_executable(test_calculate_visibilities_nobeam_app
    test_calculate_visibilities_nobeam.c
    test_calculate_visibilities_common.c
    ../../src/FEE_primary_beam.c
    ../../src/visibility_set.c
    ../../src/woden_settings.c
    ../../src/shapelet_basis.c
)
target_link_libraries(test_calculate_visibilities_nobeam_app
    Unity
    calculate_visibilitiesCUDA
    ${CC_LINKLIBS}
)
add_test(CUDA_test_calculate_visibilities_nobeam test_calculate_visibilities_nobeam_app)


##Test the main calculate_visibilities function using the gaussian primary beam
add_executable(test_calculate_visibilities_gaussbeam_app
    test_calculate_visibilities_gaussbeam.c
    test_calculate_visibilities_common.c
    ../../src/FEE_primary_beam.c
    ../../src/visibility_set.c
    ../../src/woden_settings.c
    ../../src/shapelet_basis.c
)
target_link_libraries(test_calculate_visibilities_gaussbeam_app
    Unity
    calculate_visibilitiesCUDA
    ${CC_LINKLIBS}
)
add_test(CUDA_test_calculate_visibilities_gaussbeam test_calculate_visibilities_gaussbeam_app)

# ##Test the main calculate_visibilities function using the EDA2 primary beam
add_executable(test_calculate_visibilities_edabeam_app
    test_calculate_visibilities_edabeam.c
    test_calculate_visibilities_common.c
    ../../src/FEE_primary_beam.c
    ../../src/visibility_set.c
    ../../src/woden_settings.c
    ../../src/shapelet_basis.c
)
target_link_libraries(test_calculate_visibilities_edabeam_app
    Unity
    calculate_visibilitiesCUDA
    ${CC_LINKLIBS}
)
add_test(CUDA_test_calculate_visibilities_edabeam test_calculate_visibilities_edabeam_app)

# ##Test the main calculate_visibilities function using the mwa fee primary beam
add_executable(test_calculate_visibilities_mwafeebeam_app
    test_calculate_visibilities_mwafeebeam.c
    test_calculate_visibilities_common.c
    ${CMAKE_SOURCE_DIR}/src/FEE_primary_beam.c
    ${CMAKE_SOURCE_DIR}/src/visibility_set.c
    ${CMAKE_SOURCE_DIR}/src/woden_settings.c
    ${CMAKE_SOURCE_DIR}/src/shapelet_basis.c
)
target_link_libraries(test_calculate_visibilities_mwafeebeam_app
    Unity
    calculate_visibilitiesCUDA
    ${CC_LINKLIBS}
)
target_compile_options(test_calculate_visibilities_mwafeebeam_app PRIVATE -g -Wall)

add_test(CUDA_test_calculate_visibilities_mwafeebeam test_calculate_visibilities_mwafeebeam_app)

# Emergency code for testing ctest environment
# add_test(run_xterm xterm)
