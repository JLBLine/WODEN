## Will be doing a bunch of sym linking to files for testing. Make sure we
##get the absolute path OK here
find_path(LAYOUT_PATH "example_array_layout.txt" HINTS "${CMAKE_CURRENT_LIST_DIR}")

set(A1 "example_array_layout.txt")
## Loop through them and write a symlink so the testing executable can find them
foreach(X IN LISTS A1)
    ADD_CUSTOM_TARGET(target_${X} ALL
        COMMAND ${CMAKE_COMMAND} -E create_symlink "${LAYOUT_PATH}/${X}" "${X}")
endforeach()

foreach(PRECISION IN LISTS FLOAT DOUBLE)

  ##Set the flags depending on the precision
  DEFINE_COMP_FLAGS(${PRECISION} C_FLAGS CUDA_FLAGS C_COVER_FLAGS)

  add_library(array_layout_${PRECISION} SHARED "${CMAKE_SOURCE_DIR}/src/array_layout.c")
  target_link_libraries(array_layout_${PRECISION} gcov)
  target_compile_options(array_layout_${PRECISION} PRIVATE ${C_COVER_FLAGS})

  ## Accumulate the coverage from the tests
  # add_coverage(array_layout)

  # ##Test the east,north,height to X,Y,Z works
  add_executable(test_RTS_ENH2XYZ_local_${PRECISION}_app
      test_RTS_ENH2XYZ_local.c
  )
  target_link_libraries(test_RTS_ENH2XYZ_local_${PRECISION}_app
      Unity gcov
      array_layout_${PRECISION}
      ${CC_LINKLIBS}
  )
  target_compile_options(test_RTS_ENH2XYZ_local_${PRECISION}_app PRIVATE ${C_FLAGS})
  add_test(C_test_RTS_ENH2XYZ_local_${PRECISION} test_RTS_ENH2XYZ_local_${PRECISION}_app)


  # ##Test calculating baseline length in X,Y,Z works
  add_executable(test_calc_XYZ_diffs_${PRECISION}_app
      test_calc_XYZ_diffs.c
  )
  target_link_libraries(test_calc_XYZ_diffs_${PRECISION}_app
      Unity gcov
      array_layout_${PRECISION}
      ${CC_LINKLIBS}
  )
  target_compile_options(test_calc_XYZ_diffs_${PRECISION}_app PRIVATE ${C_FLAGS})
  add_test(C_test_calc_XYZ_diffs_${PRECISION} test_calc_XYZ_diffs_${PRECISION}_app)
  
  
  add_executable(test_RTS_PrecessXYZtoJ2000_${PRECISION}_app
      test_RTS_PrecessXYZtoJ2000.c
  )
  target_link_libraries(test_RTS_PrecessXYZtoJ2000_${PRECISION}_app
      Unity gcov
      array_layout_${PRECISION}
      ${CC_LINKLIBS}
  )
  target_compile_options(test_RTS_PrecessXYZtoJ2000_${PRECISION}_app PRIVATE ${C_FLAGS})
  add_test(C_test_RTS_PrecessXYZtoJ2000_${PRECISION} test_RTS_PrecessXYZtoJ2000_${PRECISION}_app)

endforeach()
