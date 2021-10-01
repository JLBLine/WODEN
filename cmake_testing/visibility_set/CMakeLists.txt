add_library(visibility_set SHARED "${CMAKE_SOURCE_DIR}/src/visibility_set.c")
## Accumulate the coverage from the tests
add_coverage(visibility_set)


##Test settings for visibility_set are added correctly
add_executable(test_fill_timefreq_visibility_set_app
    test_fill_timefreq_visibility_set.c
)
target_link_libraries(test_fill_timefreq_visibility_set_app
    Unity
    visibility_set
    ${CC_LINKLIBS}
)
add_test(C_test_fill_timefreq_visibility_set test_fill_timefreq_visibility_set_app)

##Test some malloc and freeing functions
add_executable(test_malloc_and_free_app
    test_malloc_and_free.c
)
target_link_libraries(test_malloc_and_free_app
    Unity
    visibility_set
    ${CC_LINKLIBS}
)
add_test(C_test_malloc_and_free test_malloc_and_free_app)

##Test writing the binary file
add_executable(test_write_visi_set_binary_app
    test_write_visi_set_binary.c
)
target_link_libraries(test_write_visi_set_binary_app
    Unity
    visibility_set
    ${CC_LINKLIBS}
)
add_test(C_test_write_visi_set_binary test_write_visi_set_binary_app)
