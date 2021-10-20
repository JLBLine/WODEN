
add_library(primary_beam_cudaCUDA SHARED
            "${CMAKE_SOURCE_DIR}/src/primary_beam_cuda.cu"
            "${CMAKE_SOURCE_DIR}/src/fundamental_coords.cu")



## For CMake versions > 3.18, the arch flag settings are handled by
## CUDA_ARCHITECTURES. For now, just switch this warning off - in future,
## replace the CUDA_FLAG lines above with something like
set_property(TARGET primary_beam_cudaCUDA PROPERTY CUDA_ARCHITECTURES OFF)

##Test gaussian beam functions
add_executable(test_gaussian_beam_app
    test_gaussian_beam.c
)
target_link_libraries(test_gaussian_beam_app
    Unity
    primary_beam_cudaCUDA
    ${CC_LINKLIBS}
)
add_test(CUDA_test_gaussian_beam test_gaussian_beam_app)



##Test analytic beam functions
add_executable(test_analytic_dipole_beam_app
    test_analytic_dipole_beam.c
)
target_link_libraries(test_analytic_dipole_beam_app
    Unity
    primary_beam_cudaCUDA
    ${CC_LINKLIBS}
)
add_test(CUDA_test_analytic_dipole_beam test_analytic_dipole_beam_app)
