add_library(calculate_visibilitiesCUDA SHARED
            "../../src/calculate_visibilities.cu"
            "../../src/source_components.cu"
            "../../src/fundamental_coords.cu"
            "../../src/primary_beam_cuda.cu"
            "../../src/FEE_primary_beam_cuda.cu")

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


# ##Test the main calculate_visibilities function using the gaussian primary beam
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
