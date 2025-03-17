

## Create tests for both float and double precision
foreach(PRECISION IN LISTS FLOAT DOUBLE)
# foreach(PRECISION IN LISTS DOUBLE)
  ##Set the flags depending on the precision
  DEFINE_COMP_FLAGS(${PRECISION} C_FLAGS GPU_FLAGS C_COVER_FLAGS)

  # ##Now loop over the cpu and gpu tests
  # foreach(DEVICE IN LISTS CPU GPU)
  set(DEVICE "GPU")

    ##Test l,m,n functions
    add_executable(test_run_woden_${DEVICE}_${PRECISION}_app
                   test_run_woden.c
                   ${CMAKE_SOURCE_DIR}/cmake_testing/GPU_or_C_code/calculate_visibilities/calculate_visibilities_common_common.c
    )
    target_link_libraries(test_run_woden_${DEVICE}_${PRECISION}_app
                          Unity gcov
                          calculate_visibilities_CPU_${PRECISION}
                          calculate_visibilities_GPU_${PRECISION}
                          ${CC_LINKLIBS}
    )
    target_include_directories(test_run_woden_${DEVICE}_${PRECISION}_app PRIVATE
                               ${CMAKE_SOURCE_DIR}/cmake_testing/GPU_or_C_code/calculate_visibilities
    )
    target_compile_options(test_run_woden_${DEVICE}_${PRECISION}_app PRIVATE ${C_FLAGS})
    add_test(${DEVICE}_test_run_woden_${PRECISION} test_run_woden_${DEVICE}_${PRECISION}_app)
  
    # endforeach()
endforeach()
