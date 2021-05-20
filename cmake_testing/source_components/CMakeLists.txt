
add_library(source_componentsCUDA SHARED
            "../../src/source_components.cu"
            "../../src/fundamental_coords.cu"
            "../../src/primary_beam_cuda.cu"
            "../../src/FEE_primary_beam_cuda.cu")

## For CMake versions > 3.18, the arch flag settings are handled by
## CUDA_ARCHITECTURES. For now, just switch this warning off - in future,
## replace the CUDA_FLAG lines above with something like
set_property(TARGET source_componentsCUDA PROPERTY CUDA_ARCHITECTURES OFF)

include_directories("${CMAKE_SOURCE_DIR}/cmake_testing/source_components")

##Test flux extrapolation
add_executable(test_extrap_stokes_app
    test_extrap_stokes.c
)
target_link_libraries(test_extrap_stokes_app
    Unity
    source_componentsCUDA
    ${CC_LINKLIBS}
)
add_test(CUDA_test_extrap_stokes test_extrap_stokes_app)

##Test measurement equation
add_executable(test_calc_measurement_equation_app
    test_calc_measurement_equation.c
)
target_link_libraries(test_calc_measurement_equation_app
    Unity
    source_componentsCUDA
    ${CC_LINKLIBS}
)
add_test(CUDA_test_calc_measurement_equation test_calc_measurement_equation_app)

##Test applying primary beam gains to visibilities
add_executable(test_apply_beam_gains_app
    test_apply_beam_gains.c
)
target_link_libraries(test_apply_beam_gains_app
    Unity
    source_componentsCUDA
    ${CC_LINKLIBS}
)
add_test(CUDA_test_apply_beam_gains test_apply_beam_gains_app)

##Test grabbing the correct primary beam gains
add_executable(test_get_beam_gains_app
    test_get_beam_gains.c
)
target_link_libraries(test_get_beam_gains_app
    Unity
    source_componentsCUDA
    ${CC_LINKLIBS}
)
add_test(CUDA_test_get_beam_gains test_get_beam_gains_app)

##Test applying primary beam gains and extrapolated fluxes to visibilities
add_executable(test_update_sum_visis_app
    test_update_sum_visis.c
)
target_link_libraries(test_update_sum_visis_app
    Unity
    source_componentsCUDA
    ${CC_LINKLIBS}
)
add_test(CUDA_test_update_sum_visis test_update_sum_visis_app)

##Test the common calls between all component types return expected values
add_executable(test_source_component_common_app
    test_source_component_common.c
    ../../src/FEE_primary_beam.c
)
target_link_libraries(test_source_component_common_app
    Unity
    source_componentsCUDA
    ${CC_LINKLIBS}
)
add_test(CUDA_test_source_component_common test_source_component_common_app)

# ##Test the point source visibility kernel
add_executable(test_kern_calc_visi_point_app
    test_kern_calc_visi_point.c
    test_kern_calc_visi_common.c
    ../../src/shapelet_basis.c
    ../../src/FEE_primary_beam.c
)
target_link_libraries(test_kern_calc_visi_point_app
    Unity
    source_componentsCUDA
    ${CC_LINKLIBS}
)
add_test(CUDA_test_kern_calc_visi_point test_kern_calc_visi_point_app)

# ##Test the Gaussian visibility kernel
add_executable(test_kern_calc_visi_gauss_app
    test_kern_calc_visi_gauss.c
    test_kern_calc_visi_common.c
    ../../src/shapelet_basis.c
    ../../src/FEE_primary_beam.c
)
target_link_libraries(test_kern_calc_visi_gauss_app
    Unity
    source_componentsCUDA
    ${CC_LINKLIBS}
)
add_test(CUDA_test_kern_calc_visi_gauss test_kern_calc_visi_gauss_app)

# ##Test the shapelet visibility kernel
add_executable(test_kern_calc_visi_shape_app
    test_kern_calc_visi_shape.c
    test_kern_calc_visi_common.c
    ../../src/shapelet_basis.c
    ../../src/FEE_primary_beam.c
)
target_link_libraries(test_kern_calc_visi_shape_app
    Unity
    source_componentsCUDA
    ${CC_LINKLIBS}
)
add_test(CUDA_test_kern_calc_visi_shape test_kern_calc_visi_shape_app)
