
## Create tests for both float and double precision
foreach(PRECISION IN LISTS FLOAT DOUBLE)
# foreach(PRECISION IN LISTS DOUBLE)

  ##Set the flags depending on the precision
  DEFINE_COMP_FLAGS(${PRECISION} C_FLAGS CUDA_FLAGS C_COVER_FLAGS)

  add_library(chunk_sky_model_${PRECISION} SHARED "${CMAKE_SOURCE_DIR}/src/chunk_sky_model.c")
  target_link_libraries(chunk_sky_model_${PRECISION} gcov)
  target_compile_options(chunk_sky_model_${PRECISION} PRIVATE ${C_COVER_FLAGS})
  ## Accumulate the coverage from the tests
  # add_coverage(chunk_sky_model_${PRECISION})

  ##Test the COMPONENT nulling function
  add_executable(test_null_comps_${PRECISION}_app
      test_null_comps.c
  )
  target_link_libraries(test_null_comps_${PRECISION}_app
      Unity gcov
      chunk_sky_model_${PRECISION}
      ${CC_LINKLIBS}
  )
  target_compile_options(test_null_comps_${PRECISION}_app PRIVATE ${C_FLAGS})
  add_test(C_test_null_comps_${PRECISION} test_null_comps_${PRECISION}_app)

  ##Test the POINT and GAUSSIAN chunk filling functions works
  add_executable(test_fill_chunk_src_with_pointgauss_${PRECISION}_app
      test_fill_chunk_src_with_pointgauss.c
      fill_chunks_common.c
  )
  target_link_libraries(test_fill_chunk_src_with_pointgauss_${PRECISION}_app
      Unity gcov
      chunk_sky_model_${PRECISION}
      ${CC_LINKLIBS}
  )
  target_compile_options(test_fill_chunk_src_with_pointgauss_${PRECISION}_app PRIVATE ${C_FLAGS})
  add_test(C_test_fill_chunk_src_with_pointgauss_${PRECISION}
           test_fill_chunk_src_with_pointgauss_${PRECISION}_app)

  ## Test the SHAPELET chunk filling functions works
  add_executable(test_fill_chunk_src_with_shapelets_${PRECISION}_app
      test_fill_chunk_src_with_shapelets.c
      fill_chunks_common.c
  )
  target_link_libraries(test_fill_chunk_src_with_shapelets_${PRECISION}_app
      Unity gcov
      chunk_sky_model_${PRECISION}
      ${CC_LINKLIBS}
  )
  target_compile_options(test_fill_chunk_src_with_shapelets_${PRECISION}_app PRIVATE ${C_FLAGS})
  add_test(C_test_fill_chunk_src_with_shapelets_${PRECISION}
           test_fill_chunk_src_with_shapelets_${PRECISION}_app)

  ## Test the overall chunking function that wraps all this functionality works
  add_executable(test_create_chunked_sky_models_${PRECISION}_app
      test_create_chunked_sky_models.c
      fill_chunks_common.c
  )
  target_link_libraries(test_create_chunked_sky_models_${PRECISION}_app
      Unity gcov
      chunk_sky_model_${PRECISION}
      ${CC_LINKLIBS}
  )
  target_compile_options(test_create_chunked_sky_models_${PRECISION}_app PRIVATE ${C_FLAGS})
  add_test(C_test_create_chunked_sky_models_${PRECISION}
           test_create_chunked_sky_models_${PRECISION}_app)

endforeach()
