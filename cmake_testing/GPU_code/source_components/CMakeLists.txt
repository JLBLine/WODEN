##Set some directories we need the include files from
include_directories("${CMAKE_SOURCE_DIR}/cmake_testing/source_components")

## Create tests for both float and double precision
foreach(PRECISION IN LISTS FLOAT DOUBLE)
# foreach(PRECISION IN LISTS DOUBLE)
  ##Set the flags depending on the precision
  DEFINE_COMP_FLAGS(${PRECISION} C_FLAGS GPU_FLAGS C_COVER_FLAGS)

  if(USE_CUDA)
    set_source_files_properties("${CMAKE_SOURCE_DIR}/src/source_components.cpp"
              "${CMAKE_SOURCE_DIR}/src/fundamental_coords.cpp"
              "${CMAKE_SOURCE_DIR}/src/primary_beam_gpu.cpp" PROPERTIES LANGUAGE CUDA)
  endif()

  add_library(source_componentsGPU_${PRECISION} SHARED
              "${CMAKE_SOURCE_DIR}/src/source_components.cpp"
              "${CMAKE_SOURCE_DIR}/src/fundamental_coords.cpp"
              "${CMAKE_SOURCE_DIR}/src/primary_beam_gpu.cpp")

  target_compile_options(source_componentsGPU_${PRECISION} PRIVATE ${GPU_FLAGS})

  ##----------------------------------------------------------------------------
  ##Now to add the tests in
  ##----------------------------------------------------------------------------

  # Test flux extrapolation
  add_executable(test_extrap_stokes_${PRECISION}_app
      test_extrap_stokes.c
      common_testing_functions.c
  )
  target_link_libraries(test_extrap_stokes_${PRECISION}_app
      Unity
      source_componentsGPU_${PRECISION}
      ${CC_LINKLIBS}
  )
  target_compile_options(test_extrap_stokes_${PRECISION}_app PRIVATE ${C_FLAGS})
  add_test(GPU_test_extrap_stokes_${PRECISION} test_extrap_stokes_${PRECISION}_app)

  ##Test measurement equation
  add_executable(test_calc_measurement_equation_${PRECISION}_app
      test_calc_measurement_equation.c
  )
  target_link_libraries(test_calc_measurement_equation_${PRECISION}_app
      Unity
      source_componentsGPU_${PRECISION}
      ${CC_LINKLIBS}
  )
  target_compile_options(test_calc_measurement_equation_${PRECISION}_app
                         PRIVATE ${C_FLAGS})
  add_test(GPU_test_calc_measurement_equation_${PRECISION}
           test_calc_measurement_equation_${PRECISION}_app)

  ##Test ${PRECISION}_applying primary beam gains to visibilities
  add_executable(test_apply_beam_gains_${PRECISION}_app
      test_apply_beam_gains.c
  )
  target_link_libraries(test_apply_beam_gains_${PRECISION}_app
      Unity
      source_componentsGPU_${PRECISION}
      ${CC_LINKLIBS}
  )
  target_compile_options(test_apply_beam_gains_${PRECISION}_app
                         PRIVATE ${C_FLAGS})
  add_test(GPU_test_apply_beam_gains_${PRECISION}
           test_apply_beam_gains_${PRECISION}_app)

  ##Test grabbing the correct primary beam gains
  add_executable(test_get_beam_gains_${PRECISION}_app
      test_get_beam_gains.c
  )
  target_link_libraries(test_get_beam_gains_${PRECISION}_app
      Unity
      source_componentsGPU_${PRECISION}
      ${CC_LINKLIBS}
  )
  target_compile_options(test_get_beam_gains_${PRECISION}_app
                         PRIVATE ${C_FLAGS})
  add_test(GPU_test_get_beam_gains_${PRECISION}
           test_get_beam_gains_${PRECISION}_app)


  ##Test the common calls between all component types return expected values
  add_executable(test_source_component_common_${PRECISION}_app
      test_source_component_common.c
      common_testing_functions.c
  )
  target_link_libraries(test_source_component_common_${PRECISION}_app
      Unity
      source_componentsGPU_${PRECISION}
      ${CC_LINKLIBS}
  )
  target_compile_options(test_source_component_common_${PRECISION}_app
                         PRIVATE ${C_FLAGS})
  add_test(GPU_test_source_component_common_${PRECISION}
           test_source_component_common_${PRECISION}_app)

  ##Test the point source visibility kernel
  add_executable(test_kern_calc_visi_point_${PRECISION}_app
      test_kern_calc_visi_point.c
      test_kern_calc_visi_common.c
      common_testing_functions.c
      ${CMAKE_SOURCE_DIR}/src/shapelet_basis.c
  )
  target_link_libraries(test_kern_calc_visi_point_${PRECISION}_app
      Unity
      source_componentsGPU_${PRECISION}
      ${CC_LINKLIBS}
  )
  target_compile_options(test_kern_calc_visi_point_${PRECISION}_app
                         PRIVATE ${C_FLAGS})
  add_test(GPU_test_kern_calc_visi_point_${PRECISION}
           test_kern_calc_visi_point_${PRECISION}_app)

  ## Test the Gaussian visibility kernel
  add_executable(test_kern_calc_visi_gauss_${PRECISION}_app
      test_kern_calc_visi_gauss.c
      test_kern_calc_visi_common.c
      common_testing_functions.c
      ${CMAKE_SOURCE_DIR}/src/shapelet_basis.c
  )
  target_link_libraries(test_kern_calc_visi_gauss_${PRECISION}_app
      Unity
      source_componentsGPU_${PRECISION}
      ${CC_LINKLIBS}
  )
  target_compile_options(test_kern_calc_visi_gauss_${PRECISION}_app
                         PRIVATE ${C_FLAGS})
  add_test(GPU_test_kern_calc_visi_gauss_${PRECISION}
           test_kern_calc_visi_gauss_${PRECISION}_app)

  ##Test the shapelet visibility kernel
  add_executable(test_kern_calc_visi_shape_${PRECISION}_app
      test_kern_calc_visi_shape.c
      test_kern_calc_visi_common.c
      common_testing_functions.c
      ${CMAKE_SOURCE_DIR}/src/shapelet_basis.c
  )
  target_link_libraries(test_kern_calc_visi_shape_${PRECISION}_app
      Unity
      source_componentsGPU_${PRECISION}
      ${CC_LINKLIBS}
  )
  target_compile_options(test_kern_calc_visi_shape_${PRECISION}_app
                         PRIVATE ${C_FLAGS})
  add_test(GPU_test_kern_calc_visi_shape_${PRECISION}
           test_kern_calc_visi_shape_${PRECISION}_app)

  ##Test applying primary beam gains and extrapolated fluxes to visibilities
  add_executable(test_update_sum_visis_${PRECISION}_app
      test_update_sum_visis.c
  )
  target_link_libraries(test_update_sum_visis_${PRECISION}_app
      Unity
      source_componentsGPU_${PRECISION}
      ${CC_LINKLIBS}
  )
  target_compile_options(test_update_sum_visis_${PRECISION}_app
                         PRIVATE ${C_FLAGS})
  add_test(GPU_test_update_sum_visis_${PRECISION}
           test_update_sum_visis_${PRECISION}_app)


  # ##Test the function that calculates the autocorrelations
  add_executable(test_kern_calc_autos_${PRECISION}_app
      test_kern_calc_autos.c
      ${CMAKE_SOURCE_DIR}/src/visibility_set.c
  )
  target_link_libraries(test_kern_calc_autos_${PRECISION}_app
      source_componentsGPU_${PRECISION}
      ${CC_LINKLIBS}
      Unity
  )
  target_compile_options(test_kern_calc_autos_${PRECISION}_app
                         PRIVATE ${C_FLAGS})
  add_test(GPU_test_kern_calc_autos_${PRECISION}
           test_kern_calc_autos_${PRECISION}_app)


  ##Test grabbing the correct primary beam gains when using two per baseline
  add_executable(test_get_beam_gains_two_antennas_${PRECISION}_app
  test_get_beam_gains_two_antennas.c
  )
  target_link_libraries(test_get_beam_gains_two_antennas_${PRECISION}_app
      Unity
      source_componentsGPU_${PRECISION}
      ${CC_LINKLIBS}
  )
  target_compile_options(test_get_beam_gains_two_antennas_${PRECISION}_app
                      PRIVATE ${C_FLAGS})
  add_test(GPU_test_get_beam_gains_two_antennas_${PRECISION}
          test_get_beam_gains_two_antennas_${PRECISION}_app)

  add_executable(test_update_sum_visis_multiants_${PRECISION}_app
  test_update_sum_visis_multiants.c
  )
  target_link_libraries(test_update_sum_visis_multiants_${PRECISION}_app
  Unity
  source_componentsGPU_${PRECISION}
  ${CC_LINKLIBS}
  )
  target_compile_options(test_update_sum_visis_multiants_${PRECISION}_app
                    PRIVATE ${C_FLAGS})
  add_test(GPU_test_update_sum_visis_multiants_${PRECISION}
      test_update_sum_visis_multiants_${PRECISION}_app)

  ##Test the common calls between all component types return expected values
  add_executable(test_source_component_common_multiants_${PRECISION}_app
  test_source_component_common_multiants.c
  common_testing_functions.c
  )
  target_link_libraries(test_source_component_common_multiants_${PRECISION}_app
      Unity
      source_componentsGPU_${PRECISION}
      ${CC_LINKLIBS}
  )
  target_compile_options(test_source_component_common_multiants_${PRECISION}_app
                        PRIVATE ${C_FLAGS})
  add_test(GPU_test_source_component_common_multiants_${PRECISION}
          test_source_component_common_multiants_${PRECISION}_app)

   # ##Test the function that calculates the autocorrelations
  add_executable(test_kern_calc_autos_multiants_${PRECISION}_app
  test_kern_calc_autos_multiants.c
  ${CMAKE_SOURCE_DIR}/src/visibility_set.c
  )
  target_link_libraries(test_kern_calc_autos_multiants_${PRECISION}_app
  source_componentsGPU_${PRECISION}
  ${CC_LINKLIBS}
  Unity
  )
  target_compile_options(test_kern_calc_autos_multiants_${PRECISION}_app
                      PRIVATE ${C_FLAGS})
  add_test(GPU_test_kern_calc_autos_multiants_${PRECISION}
      test_kern_calc_autos_multiants_${PRECISION}_app)


endforeach()
