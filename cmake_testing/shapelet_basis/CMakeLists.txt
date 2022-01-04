## Create tests for both float and double precision
foreach(PRECISION IN LISTS FLOAT DOUBLE)

  ##Set the flags depending on the precision
  DEFINE_COMP_FLAGS(${PRECISION} C_FLAGS CUDA_FLAGS C_COVER_FLAGS)

  add_library(shapelet_basis_${PRECISION} SHARED "${CMAKE_SOURCE_DIR}/src/shapelet_basis.c")
  target_link_libraries(shapelet_basis_${PRECISION} gcov)
  target_compile_options(shapelet_basis_${PRECISION} PRIVATE ${C_COVER_FLAGS})
  ## Accumulate the coverage from the tests
  # add_coverage(shapelet_basis)

  # ##Test the shapelet basis array creation works
  add_executable(test_create_sbf_${PRECISION}_app
      test_create_sbf.c
  )
  target_link_libraries(test_create_sbf_${PRECISION}_app
      Unity gcov
      shapelet_basis_${PRECISION}
      ${CC_LINKLIBS}
  )
  target_compile_options(test_create_sbf_${PRECISION}_app PRIVATE ${C_FLAGS})
  add_test(C_test_create_sbf_${PRECISION} test_create_sbf_${PRECISION}_app)

endforeach()
