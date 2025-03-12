##Set some directories we need the include files from
include_directories("${CMAKE_SOURCE_DIR}/cmake_testing/source_components")

## Create tests for both float and double precision
foreach(PRECISION IN LISTS FLOAT DOUBLE)
# foreach(PRECISION IN LISTS DOUBLE)
  ##Set the flags depending on the precision
  DEFINE_COMP_FLAGS(${PRECISION} C_FLAGS GPU_FLAGS C_COVER_FLAGS)

  if(USE_CUDA)
    set_source_files_properties(${CMAKE_SOURCE_DIR}/src/fundamental_coords_gpu.cpp
                                ${CMAKE_SOURCE_DIR}/src/primary_beam_gpu.cpp
                                ${CMAKE_SOURCE_DIR}/src/source_components_gpu.cpp
                                extrap_stokes_common_gpu.cpp
                                calc_measurement_equation_common_gpu.cpp
                                apply_beam_gains_common_gpu.cpp
                                get_beam_gains_common_gpu.cpp
                                source_component_common_common_gpu.cpp
                                calc_visi_common_gpu.cpp
                                update_sum_visis_common_gpu.cpp
                                calc_autos_common_gpu.cpp
                                PROPERTIES LANGUAGE CUDA)
  endif()

  ##compile GPU versions of the code
  add_library(source_components_GPU_${PRECISION} SHARED
              ${CMAKE_SOURCE_DIR}/src/fundamental_coords_gpu.cpp
              ${CMAKE_SOURCE_DIR}/src/primary_beam_gpu.cpp
              ${CMAKE_SOURCE_DIR}/src/source_components_gpu.cpp
              extrap_stokes_common_gpu.cpp
              calc_measurement_equation_common_gpu.cpp
              apply_beam_gains_common_gpu.cpp
              get_beam_gains_common_gpu.cpp
              source_component_common_common_gpu.cpp
              calc_visi_common_gpu.cpp
              update_sum_visis_common_gpu.cpp
              calc_autos_common_gpu.cpp)
  target_compile_options(source_components_GPU_${PRECISION} PRIVATE ${GPU_FLAGS})

  ##----------------------------------------------------------------------------
  ##Now to add the tests in
  ##----------------------------------------------------------------------------

  foreach(DEVICE IN LISTS CPU GPU)
  # foreach(DEVICE IN LISTS GPU)
  # foreach(DEVICE IN LISTS CPU)
    ##Test measurement equation
    add_executable(test_calc_measurement_equation_${DEVICE}_${PRECISION}_app
      test_calc_measurement_equation_${DEVICE}.c
      calc_measurement_equation_common.c
    )
    target_link_libraries(test_calc_measurement_equation_${DEVICE}_${PRECISION}_app
      Unity gcov
      source_components_GPU_${PRECISION}
      calculate_visibilities_CPU_${PRECISION}
      calculate_visibilities_GPU_${PRECISION}
      ${CC_LINKLIBS}
    )
    target_compile_options(test_calc_measurement_equation_${DEVICE}_${PRECISION}_app
                        PRIVATE ${C_FLAGS})
    add_test(${DEVICE}_test_calc_measurement_equation_${PRECISION}
      test_calc_measurement_equation_${DEVICE}_${PRECISION}_app)

    add_executable(test_apply_beam_gains_${DEVICE}_${PRECISION}_app
      test_apply_beam_gains_${DEVICE}.c
      apply_beam_gains_common.c
    )
    target_link_libraries(test_apply_beam_gains_${DEVICE}_${PRECISION}_app
      Unity gcov
      source_components_GPU_${PRECISION}
      calculate_visibilities_CPU_${PRECISION}
        calculate_visibilities_GPU_${PRECISION}
      ${CC_LINKLIBS}
    )
    target_compile_options(test_apply_beam_gains_${DEVICE}_${PRECISION}_app
                        PRIVATE ${C_FLAGS})
    add_test(${DEVICE}_test_apply_beam_gains_${PRECISION}
          test_apply_beam_gains_${DEVICE}_${PRECISION}_app)

    # # test Stokes extrapolation
    add_executable(test_extrap_stokes_${DEVICE}_${PRECISION}_app
      test_extrap_stokes_${DEVICE}.c
      common_testing_functions.c
      extrap_stokes_common.c
    )
    target_link_libraries(test_extrap_stokes_${DEVICE}_${PRECISION}_app
      Unity gcov
      source_components_GPU_${PRECISION}
      calculate_visibilities_CPU_${PRECISION}
        calculate_visibilities_GPU_${PRECISION}
      ${CC_LINKLIBS}
    )
    target_compile_options(test_extrap_stokes_${DEVICE}_${PRECISION}_app
                        PRIVATE ${C_FLAGS})
    add_test(${DEVICE}_test_extrap_stokes_${PRECISION}
          test_extrap_stokes_${DEVICE}_${PRECISION}_app)



    add_executable(test_get_beam_gains_${DEVICE}_${PRECISION}_app
      test_get_beam_gains_${DEVICE}.c
      get_beam_gains_common.c
    )
    target_link_libraries(test_get_beam_gains_${DEVICE}_${PRECISION}_app
      Unity gcov
      source_components_GPU_${PRECISION}
      calculate_visibilities_CPU_${PRECISION}
        calculate_visibilities_GPU_${PRECISION}
      ${CC_LINKLIBS}
    )
    target_compile_options(test_get_beam_gains_${DEVICE}_${PRECISION}_app
                        PRIVATE ${C_FLAGS})
    add_test(${DEVICE}_test_get_beam_gains_${PRECISION}
          test_get_beam_gains_${DEVICE}_${PRECISION}_app)

    add_executable(test_get_beam_gains_multiants_${DEVICE}_${PRECISION}_app
      test_get_beam_gains_multiants_${DEVICE}.c
      get_beam_gains_multiants_common.c
      common_testing_functions.c
    )
    target_link_libraries(test_get_beam_gains_multiants_${DEVICE}_${PRECISION}_app
      Unity gcov
      source_components_GPU_${PRECISION}
      calculate_visibilities_CPU_${PRECISION}
        calculate_visibilities_GPU_${PRECISION}
      ${CC_LINKLIBS}
    )
    target_compile_options(test_get_beam_gains_multiants_${DEVICE}_${PRECISION}_app
                        PRIVATE ${C_FLAGS})
    add_test(${DEVICE}_test_get_beam_gains_multiants_${PRECISION}
          test_get_beam_gains_multiants_${DEVICE}_${PRECISION}_app)

    add_executable(test_update_sum_visis_${DEVICE}_${PRECISION}_app
      test_update_sum_visis_${DEVICE}.c
      update_sum_visis_common.c
      common_testing_functions.c
    )
    target_link_libraries(test_update_sum_visis_${DEVICE}_${PRECISION}_app
      Unity gcov
      source_components_GPU_${PRECISION}
      calculate_visibilities_CPU_${PRECISION}
        calculate_visibilities_GPU_${PRECISION}
      ${CC_LINKLIBS}
    )
    target_compile_options(test_update_sum_visis_${DEVICE}_${PRECISION}_app
                        PRIVATE ${C_FLAGS})
    add_test(${DEVICE}_test_update_sum_visis_${PRECISION}
          test_update_sum_visis_${DEVICE}_${PRECISION}_app)


    add_executable(test_update_sum_visis_multiants_${DEVICE}_${PRECISION}_app
      test_update_sum_visis_multiants_${DEVICE}.c
      update_sum_visis_multiants_common.c
      update_sum_visis_common.c
      common_testing_functions.c
    )
    target_link_libraries(test_update_sum_visis_multiants_${DEVICE}_${PRECISION}_app
      Unity gcov
      source_components_GPU_${PRECISION}
      calculate_visibilities_CPU_${PRECISION}
        calculate_visibilities_GPU_${PRECISION}
      ${CC_LINKLIBS}
    )
    target_compile_options(test_update_sum_visis_multiants_${DEVICE}_${PRECISION}_app
                        PRIVATE ${C_FLAGS})
    add_test(${DEVICE}_test_update_sum_visis_multiants_${PRECISION}
          test_update_sum_visis_multiants_${DEVICE}_${PRECISION}_app)

    add_executable(test_calc_autos_${DEVICE}_${PRECISION}_app
      test_calc_autos_${DEVICE}.c
      calc_autos_common.c
      common_testing_functions.c
      ${CMAKE_SOURCE_DIR}/src/visibility_set.c
    )
    target_link_libraries(test_calc_autos_${DEVICE}_${PRECISION}_app
      Unity gcov
      source_components_GPU_${PRECISION}
      calculate_visibilities_CPU_${PRECISION}
        calculate_visibilities_GPU_${PRECISION}
      ${CC_LINKLIBS}
    )
    target_compile_options(test_calc_autos_${DEVICE}_${PRECISION}_app
                        PRIVATE ${C_FLAGS})
    add_test(${DEVICE}_test_calc_autos_${PRECISION}
          test_calc_autos_${DEVICE}_${PRECISION}_app)

    add_executable(test_calc_autos_multiants_${DEVICE}_${PRECISION}_app
      test_calc_autos_multiants_${DEVICE}.c
      calc_autos_multiants_common.c
      calc_autos_common.c
      common_testing_functions.c
      ${CMAKE_SOURCE_DIR}/src/visibility_set.c
    )
    target_link_libraries(test_calc_autos_multiants_${DEVICE}_${PRECISION}_app
      Unity gcov
      source_components_GPU_${PRECISION}
      calculate_visibilities_CPU_${PRECISION}
        calculate_visibilities_GPU_${PRECISION}
      ${CC_LINKLIBS}
    )
    target_compile_options(test_calc_autos_multiants_${DEVICE}_${PRECISION}_app
                        PRIVATE ${C_FLAGS})
    add_test(${DEVICE}_test_calc_autos_multiants_${PRECISION}
          test_calc_autos_multiants_${DEVICE}_${PRECISION}_app)


    add_executable(test_source_component_common_${DEVICE}_${PRECISION}_app
      test_source_component_common_${DEVICE}.c
      source_component_common_common.c
      common_testing_functions.c
    )
    target_link_libraries(test_source_component_common_${DEVICE}_${PRECISION}_app
      Unity gcov
      source_components_GPU_${PRECISION}
      calculate_visibilities_CPU_${PRECISION}
        calculate_visibilities_GPU_${PRECISION}
      ${CC_LINKLIBS}
    )
    target_compile_options(test_source_component_common_${DEVICE}_${PRECISION}_app
                        PRIVATE ${C_FLAGS})
    add_test(${DEVICE}_test_source_component_common_${PRECISION}
          test_source_component_common_${DEVICE}_${PRECISION}_app)


    add_executable(test_source_component_common_multiants_${DEVICE}_${PRECISION}_app
      test_source_component_common_multiants_${DEVICE}.c
      source_component_common_multiants_common.c
      common_testing_functions.c
    )
    target_link_libraries(test_source_component_common_multiants_${DEVICE}_${PRECISION}_app
      Unity gcov
      source_components_GPU_${PRECISION}
      calculate_visibilities_CPU_${PRECISION}
        calculate_visibilities_GPU_${PRECISION}
      ${CC_LINKLIBS}
    )
    target_compile_options(test_source_component_common_multiants_${DEVICE}_${PRECISION}_app
                        PRIVATE ${C_FLAGS})
    add_test(${DEVICE}_test_source_component_common_multiants_${PRECISION}
          test_source_component_common_multiants_${DEVICE}_${PRECISION}_app)

    add_executable(test_calc_visi_point_${DEVICE}_${PRECISION}_app
      test_calc_visi_point_${DEVICE}.c
      calc_visi_common.c
      common_testing_functions.c
      ${CMAKE_SOURCE_DIR}/src/shapelet_basis.c
    )
    target_link_libraries(test_calc_visi_point_${DEVICE}_${PRECISION}_app
      Unity gcov
      source_components_GPU_${PRECISION}
      calculate_visibilities_CPU_${PRECISION}
        calculate_visibilities_GPU_${PRECISION}
      ${CC_LINKLIBS}
    )
    target_compile_options(test_calc_visi_point_${DEVICE}_${PRECISION}_app
                        PRIVATE ${C_FLAGS})
    add_test(${DEVICE}_test_calc_visi_point_${PRECISION}
          test_calc_visi_point_${DEVICE}_${PRECISION}_app)

    add_executable(test_calc_visi_gauss_${DEVICE}_${PRECISION}_app
      test_calc_visi_gauss_${DEVICE}.c
      calc_visi_common.c
      common_testing_functions.c
      ${CMAKE_SOURCE_DIR}/src/shapelet_basis.c
    )
    target_link_libraries(test_calc_visi_gauss_${DEVICE}_${PRECISION}_app
      Unity gcov
      source_components_GPU_${PRECISION}
      calculate_visibilities_CPU_${PRECISION}
      calculate_visibilities_GPU_${PRECISION}
      ${CC_LINKLIBS}
    )
    target_compile_options(test_calc_visi_gauss_${DEVICE}_${PRECISION}_app
                        PRIVATE ${C_FLAGS})
    add_test(${DEVICE}_test_calc_visi_gauss_${PRECISION}
          test_calc_visi_gauss_${DEVICE}_${PRECISION}_app)

    add_executable(test_calc_visi_shape_${DEVICE}_${PRECISION}_app
      test_calc_visi_shape_${DEVICE}.c
      calc_visi_common.c
      common_testing_functions.c
      ${CMAKE_SOURCE_DIR}/src/shapelet_basis.c
    )
    target_link_libraries(test_calc_visi_shape_${DEVICE}_${PRECISION}_app
      Unity gcov
      source_components_GPU_${PRECISION}
      calculate_visibilities_CPU_${PRECISION}
        calculate_visibilities_GPU_${PRECISION}
      ${CC_LINKLIBS}
    )
    target_compile_options(test_calc_visi_shape_${DEVICE}_${PRECISION}_app
                        PRIVATE ${C_FLAGS})
    add_test(${DEVICE}_test_calc_visi_shape_${PRECISION}
          test_calc_visi_shape_${DEVICE}_${PRECISION}_app)


  endforeach()
endforeach()
