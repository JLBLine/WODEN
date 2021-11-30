## Create tests for both float and double precision
foreach(PRECISION IN LISTS FLOAT DOUBLE)

  ##Set the flags depending on the precision
  DEFINE_COMP_FLAGS(${PRECISION} C_FLAGS CUDA_FLAGS C_COVER_FLAGS)

  add_library(FEE_primary_beam_${PRECISION} SHARED "${CMAKE_SOURCE_DIR}/src/FEE_primary_beam.c")
  target_link_libraries(FEE_primary_beam_${PRECISION} gcov)
  target_compile_options(FEE_primary_beam_${PRECISION} PRIVATE ${C_COVER_FLAGS})
  ## Accumulate the coverage from the tests
  # add_coverage(FEE_primary_beam)

  # ##Test Initialisation of the MWA FEE bem works
  add_executable(test_RTS_MWAFEEInit_${PRECISION}_app
      test_RTS_MWAFEEInit.c
  )
  target_link_libraries(test_RTS_MWAFEEInit_${PRECISION}_app
      Unity gcov
      FEE_primary_beam_${PRECISION}
      ${CC_LINKLIBS}
  )
  target_compile_options(test_RTS_MWAFEEInit_${PRECISION}_app PRIVATE ${C_FLAGS})
  add_test(C_test_RTS_MWAFEEInit_${PRECISION} test_RTS_MWAFEEInit_${PRECISION}_app)

  # ##Test Initialisation of the MWA FEE bem works
  add_executable(test_multifreq_RTS_MWAFEEInit_${PRECISION}_app
      test_multifreq_RTS_MWAFEEInit.c
      multifreq_common.c
  )
  target_link_libraries(test_multifreq_RTS_MWAFEEInit_${PRECISION}_app
      Unity gcov
      FEE_primary_beam_${PRECISION}
      ${CC_LINKLIBS}
  )
  target_compile_options(test_multifreq_RTS_MWAFEEInit_${PRECISION}_app PRIVATE ${C_FLAGS})
  add_test(C_test_multifreq_RTS_MWAFEEInit_${PRECISION} test_multifreq_RTS_MWAFEEInit_${PRECISION}_app)

endforeach()
