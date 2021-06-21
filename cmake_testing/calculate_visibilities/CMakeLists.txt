add_library(calculate_visibilitiesCUDA SHARED
            "../../src/calculate_visibilities.cu")

## For CMake versions > 3.18, the arch flag settings are handled by
## CUDA_ARCHITECTURES. For now, just switch this warning off - in future,
## replace the CUDA_FLAG lines above with something like
set_property(TARGET calculate_visibilitiesCUDA PROPERTY CUDA_ARCHITECTURES OFF)

include_directories("${CMAKE_SOURCE_DIR}/cmake_testing/FEE_primary_beam_cuda")

# ##Test the point source visibility kernel
add_executable(test_calculate_visibilities_app
    test_calculate_visibilities.c
    ../../src/FEE_primary_beam.c
)
target_link_libraries(test_calculate_visibilities_app
    Unity
    calculate_visibilitiesCUDA
    ${CC_LINKLIBS}
)
add_test(CUDA_test_calculate_visibilities test_calculate_visibilities_app)
