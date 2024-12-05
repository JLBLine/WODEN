

## Create tests for both float and double precision
foreach(PRECISION IN LISTS FLOAT DOUBLE)
# foreach(PRECISION IN LISTS DOUBLE)
  ##Set the flags depending on the precision
  DEFINE_COMP_FLAGS(${PRECISION} C_FLAGS GPU_FLAGS C_COVER_FLAGS)

  ##If using CUDA, not HIP, then we need to set the source files to be
  ##compiled with nvcc
  if(USE_CUDA)
    set_source_files_properties(${CMAKE_SOURCE_DIR}/src/fundamental_coords_gpu.cpp
              PROPERTIES LANGUAGE CUDA)
  endif()

  ##Compile a GPU code library
  add_library(fundamental_coordsGPU_${PRECISION} SHARED
              ${CMAKE_SOURCE_DIR}/src/fundamental_coords_gpu.cpp)
  target_compile_options(fundamental_coordsGPU_${PRECISION} PRIVATE ${GPU_FLAGS})

  ##complile the CPU version of fundamental coords, and link to gcov so we
  ##can get coverage information
  add_library(fundamental_coords_cpu_${PRECISION} SHARED ${CMAKE_SOURCE_DIR}/src/fundamental_coords_cpu.c)
  target_link_libraries(fundamental_coords_cpu_${PRECISION} gcov)
  target_compile_options(fundamental_coords_cpu_${PRECISION} PRIVATE ${C_COVER_FLAGS})

  ##Now loop over the cpu and gpu tests
  foreach(DEVICE IN LISTS CPU GPU)
  # foreach(DEVICE IN LISTS GPU)
  # foreach(DEVICE IN LISTS CPU)

    ##Test l,m,n functions
    add_executable(test_lmn_coords_${DEVICE}_${PRECISION}_app
                   test_lmn_coords_${DEVICE}.c
                   lmn_coords_common.c
    )
    target_link_libraries(test_lmn_coords_${DEVICE}_${PRECISION}_app
                          Unity gcov
                          fundamental_coordsGPU_${PRECISION}
                          fundamental_coords_cpu_${PRECISION}
                          ${CC_LINKLIBS}
    )
    target_compile_options(test_lmn_coords_${DEVICE}_${PRECISION}_app PRIVATE ${C_FLAGS})
    add_test(${DEVICE}_test_lmn_coords_${PRECISION} test_lmn_coords_${DEVICE}_${PRECISION}_app)
  
    ##Test u,v,w functions
    add_executable(test_uvw_coords_${DEVICE}_${PRECISION}_app
                   test_uvw_coords_${DEVICE}.c
                   uvw_coords_common.c
    )
    target_link_libraries(test_uvw_coords_${DEVICE}_${PRECISION}_app
                          Unity gcov
                          fundamental_coordsGPU_${PRECISION}
                          fundamental_coords_cpu_${PRECISION}
                          ${CC_LINKLIBS}
    )
    target_compile_options(test_uvw_coords_${DEVICE}_${PRECISION}_app PRIVATE ${C_FLAGS})
    add_test(${DEVICE}_test_uvw_coords_${PRECISION} test_uvw_coords_${DEVICE}_${PRECISION}_app)
  endforeach()

endforeach()
