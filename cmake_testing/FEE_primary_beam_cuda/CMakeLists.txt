add_library(FEE_primary_beam_cudaCUDA SHARED
            "../../src/FEE_primary_beam_cuda.cu")
#
## For CMake versions > 3.18, the arch flag settings are handled by
## CUDA_ARCHITECTURES. For now, just switch this warning off - in future,
## replace the CUDA_FLAG lines above with something like
set_property(TARGET FEE_primary_beam_cudaCUDA PROPERTY CUDA_ARCHITECTURES OFF)

include_directories("${CMAKE_SOURCE_DIR}/cmake_testing/FEE_primary_beam_cuda")

# ##Test the point source visibility kernel
add_executable(test_RTS_FEE_beam_app
    test_RTS_FEE_beam.c
    ../../src/FEE_primary_beam.c
)
target_link_libraries(test_RTS_FEE_beam_app
    Unity
    FEE_primary_beam_cudaCUDA
    ${CC_LINKLIBS}
)
add_test(CUDA_test_RTS_FEE_beam test_RTS_FEE_beam_app)
