# # ##Test the east,north,height to X,Y,Z works
add_executable(test_RTS_ENH2XYZ_local_app
    test_RTS_ENH2XYZ_local.c
    ../../src/read_and_write.c
)
target_link_libraries(test_RTS_ENH2XYZ_local_app
    Unity
    ${CC_LINKLIBS}
)
add_test(C_test_RTS_ENH2XYZ_local test_RTS_ENH2XYZ_local_app)


# ##Make a bunch of symlinks to sky models for testing
find_path(SRCLIST_PATH "srclist_singlepoint.txt" HINTS "${CMAKE_CURRENT_LIST_DIR}")
if(SRCLIST_PATH)
  message(STATUS "SRCLIST_PATH path: ${SRCLIST_PATH}")
endif()

## Name a srclists used for testing
set(S1 "srclist_singlepoint.txt")
set(S2 "srclist_singlegauss.txt")
set(S3 "srclist_singleshape.txt")
set(S4 "srclist_threesources.txt")
set(S5 "srclist_threecomponents.txt")

## Loop through them and write a symlink so the testing executable can find them
foreach(X IN LISTS S1 S2 S3 S4 S5)
    # set(link_singlegauss ${X})
    ADD_CUSTOM_TARGET(target_${X} ALL
        COMMAND ${CMAKE_COMMAND} -E create_symlink "${SRCLIST_PATH}/${X}" "${X}")
endforeach()



# ##Test reading of WODEN catalogue works
add_executable(test_read_source_catalogue_app
    test_read_source_catalogue.c
    ../../src/read_and_write.c
)
target_link_libraries(test_read_source_catalogue_app
    Unity
    ${CC_LINKLIBS}
)
add_test(C_test_read_source_catalogue test_read_source_catalogue_app)
