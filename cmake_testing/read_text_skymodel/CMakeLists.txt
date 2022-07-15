## Will be doing a bunch of sym linking to files for testing. Make sure we
##get the absolute path OK here
find_path(SRCLIST_PATH "srclist_singlepoint.txt" HINTS "${CMAKE_SOURCE_DIR}/cmake_testing/read_text_skymodel/")

##Make a bunch of symlinks to sky models for testing
if(SRCLIST_PATH)
  message(STATUS "SRCLIST_PATH path: ${SRCLIST_PATH}")
endif()

#
# ## Name a srclists used for testing
set(S1 "srclist_singlepoint.txt")
set(S2 "srclist_singlegauss.txt")
set(S3 "srclist_singleshape.txt")
set(S4 "srclist_threesources.txt")
set(S5 "srclist_threecomponents.txt")
set(S6 "srclist_mulitple_source-components.txt")
set(S7 "srclist_empty_line.txt")
set(S8 "srclist_comment.txt")
set(S9 "srclist_badspell.txt")
set(S10 "srclist_badcoeff.txt")
set(S11 "srclist_no-comp_numbers.txt")
set(S12 "srclist_linear.txt")


## Loop through them and write a symlink so the testing executable can find them
foreach(X IN LISTS S1 S2 S3 S4 S5 S6 S7 S8 S9 S10 S11 S12)
    ADD_CUSTOM_TARGET(target_${X} ALL
        COMMAND ${CMAKE_COMMAND} -E create_symlink "${SRCLIST_PATH}/${X}" "${X}")
endforeach()

## Create tests for both float and double precision
foreach(PRECISION IN LISTS FLOAT DOUBLE)

  ##Set the flags depending on the precision
  DEFINE_COMP_FLAGS(${PRECISION} C_FLAGS CUDA_FLAGS C_COVER_FLAGS)

  ##Create a combined sky model reading library that covers both the text and
  ##yaml formats, given the code expected_skymodel_outcomes.c used below
  ##relies on both
  add_library(read_skymodel_text_${PRECISION} SHARED "${CMAKE_SOURCE_DIR}/src/read_text_skymodel.c"
  "${CMAKE_SOURCE_DIR}/src/read_yaml_skymodel.c"
  "${CMAKE_SOURCE_DIR}/src/create_sky_model.c")
  target_link_libraries(read_skymodel_text_${PRECISION} gcov)

  target_include_directories(read_skymodel_text_${PRECISION} PUBLIC "${CMAKE_SOURCE_DIR}/cmake_testing/create_sky_model/"
  "${CMAKE_SOURCE_DIR}/include")

  target_compile_options(read_skymodel_text_${PRECISION} PRIVATE ${C_COVER_FLAGS})
  ## Accumulate the coverage from the tests
  # add_coverage(read_text_skymodel_${PRECISION})




  ##Test reading of WODEN catalogue works
  add_executable(test_read_text_skymodel_${PRECISION}_app
      test_read_text_skymodel.c
      "${CMAKE_SOURCE_DIR}/cmake_testing/create_sky_model/expected_skymodel_outcomes.c"
  )
  target_link_libraries(test_read_text_skymodel_${PRECISION}_app
      Unity gcov
      read_skymodel_text_${PRECISION}
      ${CC_LINKLIBS}
  )
  target_compile_options(test_read_text_skymodel_${PRECISION}_app PRIVATE ${C_FLAGS})
  add_test(C_test_read_text_skymodel_${PRECISION} test_read_text_skymodel_${PRECISION}_app)

endforeach()
