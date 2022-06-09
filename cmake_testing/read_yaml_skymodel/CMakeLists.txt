## Will be doing a bunch of sym linking to files for testing. Make sure we
##get the absolute path OK here
find_path(YAML_PATH "srclist_comment.yaml" HINTS "${CMAKE_CURRENT_LIST_DIR}")

##Make a bunch of symlinks to sky models for testing
if(YAML_PATH)
  message(STATUS "YAML_PATH path: ${YAML_PATH}")
endif()

## Name of json inputs used for testing
set(Y1 "srclist_comment.yaml")
set(Y2 "srclist_empty_line.yaml")
set(Y3 "srclist_singlepoint_power.yaml")
set(Y4 "srclist_singlegauss_power.yaml")
set(Y5 "srclist_singleshape_power.yaml")
set(Y6 "srclist_threesources_power.yaml")
set(Y7 "srclist_threecomponents_power.yaml")

set(Y8 "srclist_singlepoint_list.yaml")
set(Y9 "srclist_singlegauss_list.yaml")
set(Y10 "srclist_singleshape_list.yaml")
set(Y11 "srclist_threesources_list.yaml")
set(Y12 "srclist_threecomponents_list.yaml")

set(Y13 "srclist_singlepoint_curve.yaml")
set(Y14 "srclist_singlegauss_curve.yaml")
set(Y15 "srclist_singleshape_curve.yaml")
set(Y16 "srclist_threesources_curve.yaml")
set(Y17 "srclist_threecomponents_curve.yaml")

set(Y18 "srclist_mulitple_source-components.yaml")

set(Y19 "srclist_badspell.yaml")
set(Y20 "srclist_badcoeff.yaml")

## Loop through them and write a symlink so the testing executable can find them
foreach(X IN LISTS Y1 Y2 Y3 Y4 Y5 Y6 Y7 Y8 Y9 Y10
Y11 Y12 Y13 Y14 Y15 Y16 Y17 Y18 Y19 Y20 Y21)
    ADD_CUSTOM_TARGET(target_${X} ALL
        COMMAND ${CMAKE_COMMAND} -E create_symlink "${YAML_PATH}/${X}" "${X}")
endforeach()

## Create tests for both float and double precision
foreach(PRECISION IN LISTS FLOAT DOUBLE)
# foreach(PRECISION IN LISTS DOUBLE)

  ##Set the flags depending on the precision
  DEFINE_COMP_FLAGS(${PRECISION} C_FLAGS CUDA_FLAGS C_COVER_FLAGS)

  ##Create a combined sky model reading library that covers both the text and
  ##yaml formats, given the code expected_skymodel_outcomes.c used below
  ##relies on both
  add_library(read_skymodel_yaml_${PRECISION} STATIC "${CMAKE_SOURCE_DIR}/src/read_text_skymodel.c"
  "${CMAKE_SOURCE_DIR}/src/read_yaml_skymodel.c"
  "${CMAKE_SOURCE_DIR}/src/create_sky_model.c")
  target_link_libraries(read_skymodel_yaml_${PRECISION} gcov)

  target_include_directories(read_skymodel_yaml_${PRECISION} PUBLIC "${CMAKE_SOURCE_DIR}/cmake_testing/create_sky_model/"
  "${CMAKE_SOURCE_DIR}/include")

  target_compile_options(read_skymodel_yaml_${PRECISION} PRIVATE ${C_COVER_FLAGS})
  ## Accumulate the coverage from the tests
  # add_coverage(read_text_skymodel_${PRECISION})


  ##Test reading of WODEN catalogue works
  add_executable(test_read_yaml_skymodel_${PRECISION}_app
      test_read_yaml_skymodel.c
      "${CMAKE_SOURCE_DIR}/cmake_testing/create_sky_model/expected_skymodel_outcomes.c"
  )
  target_link_libraries(test_read_yaml_skymodel_${PRECISION}_app
      Unity gcov
      read_skymodel_yaml_${PRECISION}
      ${CC_LINKLIBS}
  )
  target_compile_options(test_read_yaml_skymodel_${PRECISION}_app PRIVATE ${C_FLAGS})
  add_test(C_test_read_yaml_skymodel_${PRECISION} test_read_yaml_skymodel_${PRECISION}_app)

endforeach()
