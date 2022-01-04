include_directories("${CMAKE_SOURCE_DIR}/cmake_testing/FEE_primary_beam_cuda")

## Create tests for both float and double precision
foreach(PRECISION IN LISTS FLOAT DOUBLE)

  ##Set the flags depending on the precision
  DEFINE_COMP_FLAGS(${PRECISION} C_FLAGS CUDA_FLAGS C_COVER_FLAGS)

  add_library(FEE_primary_beam_cudaCUDA_${PRECISION} SHARED
              "${CMAKE_SOURCE_DIR}/src/FEE_primary_beam_cuda.cu")
  target_compile_options(FEE_primary_beam_cudaCUDA_${PRECISION} PRIVATE ${CUDA_FLAGS})

  ## For CMake versions > 3.18, the arch flag settings are handled by
  ## CUDA_ARCHITECTURES. For now, just switch this warning off - in future,
  ## replace the CUDA_FLAG lines above with something like
  set_property(TARGET FEE_primary_beam_cudaCUDA_${PRECISION} PROPERTY CUDA_ARCHITECTURES OFF)


  # ##Test the point source visibility kernel
  add_executable(test_RTS_FEE_beam_${PRECISION}_app
      test_RTS_FEE_beam.c
      ${CMAKE_SOURCE_DIR}/src/FEE_primary_beam.c
  )
  target_link_libraries(test_RTS_FEE_beam_${PRECISION}_app
      Unity
      FEE_primary_beam_cudaCUDA_${PRECISION}
      ${CC_LINKLIBS}
  )
  target_compile_options(test_RTS_FEE_beam_${PRECISION}_app PRIVATE ${C_FLAGS})
  add_test(CUDA_test_RTS_FEE_beam_${PRECISION} test_RTS_FEE_beam_${PRECISION}_app)

endforeach()
