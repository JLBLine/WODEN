

## Create tests for both float and double precision
foreach(PRECISION IN LISTS FLOAT DOUBLE)
# foreach(PRECISION IN LISTS DOUBLE)
  ##Set the flags depending on the precision
  DEFINE_COMP_FLAGS(${PRECISION} C_FLAGS CUDA_FLAGS C_COVER_FLAGS)

  if(USE_CUDA)
    set_source_files_properties(${CMAKE_SOURCE_DIR}/src/fundamental_coords.cpp
              PROPERTIES LANGUAGE CUDA)
  endif()

  add_library(fundamental_coordsCUDA_${PRECISION} SHARED
              ${CMAKE_SOURCE_DIR}/src/fundamental_coords.cpp)
  target_compile_options(fundamental_coordsCUDA_${PRECISION} PRIVATE ${CUDA_FLAGS})

  ##Test l,m,n functions
  add_executable(test_lmn_coords_${PRECISION}_app
      test_lmn_coords.c
  )
  target_link_libraries(test_lmn_coords_${PRECISION}_app
      Unity
      fundamental_coordsCUDA_${PRECISION}
  )
  target_compile_options(test_lmn_coords_${PRECISION}_app PRIVATE ${C_FLAGS})
  add_test(CUDA_test_lmn_coords_${PRECISION} test_lmn_coords_${PRECISION}_app)

  ##Test u,v,w functions
  add_executable(test_uvw_coords_${PRECISION}_app
      test_uvw_coords.c
  )
  target_link_libraries(test_uvw_coords_${PRECISION}_app
      Unity
      fundamental_coordsCUDA_${PRECISION}
      ${CC_LINKLIBS}
  )
  target_compile_options(test_uvw_coords_${PRECISION}_app PRIVATE ${C_FLAGS})
  add_test(CUDA_test_uvw_coords_${PRECISION} test_uvw_coords_${PRECISION}_app)

endforeach()
