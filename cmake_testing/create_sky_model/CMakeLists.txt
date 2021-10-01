## Will be doing a bunch of sym linking to files for testing. Make sure we
##get the absolute path OK here
find_path(SRCLIST_PATH "srclist_singlepoint.txt" HINTS "${CMAKE_CURRENT_LIST_DIR}")

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


## Loop through them and write a symlink so the testing executable can find them
foreach(X IN LISTS S1 S2 S3 S4 S5 S6 S7 S8 S9 S10 S11)
    ADD_CUSTOM_TARGET(target_${X} ALL
        COMMAND ${CMAKE_COMMAND} -E create_symlink "${SRCLIST_PATH}/${X}" "${X}")
endforeach()

add_library(create_sky_model SHARED "${CMAKE_SOURCE_DIR}/src/create_sky_model.c")
## Accumulate the coverage from the tests
add_coverage(create_sky_model)

# ##Test reading of WODEN catalogue works
add_executable(test_read_source_catalogue_app
    test_read_source_catalogue.c
)
target_link_libraries(test_read_source_catalogue_app
    Unity
    create_sky_model
    ${CC_LINKLIBS}
)
add_test(C_test_read_source_catalogue test_read_source_catalogue_app)

# ##Test the horizon test function works
add_executable(test_horizon_test_app
    test_horizon_test.c
)
target_link_libraries(test_horizon_test_app
    Unity
    create_sky_model
    ${CC_LINKLIBS}
)
add_test(C_test_horizon_test test_horizon_test_app)

# ##Test the horizon cropping function creates sky models correctly
add_executable(test_crop_sky_model_app
    test_crop_sky_model.c
)
target_link_libraries(test_crop_sky_model_app
    Unity
    create_sky_model
    ${CC_LINKLIBS}
)
add_test(C_test_crop_sky_model test_crop_sky_model_app)
