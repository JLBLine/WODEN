add_library(chunk_sky_model SHARED "${CMAKE_SOURCE_DIR}/src/chunk_sky_model.c")
## Accumulate the coverage from the tests
add_coverage(chunk_sky_model)

# ##Test the COMPONENT nulling function
add_executable(test_null_comps_app
    test_null_comps.c
)
target_link_libraries(test_null_comps_app
    Unity
    chunk_sky_model
    ${CC_LINKLIBS}
)
add_test(C_test_null_comps test_null_comps_app)

##Test the POINT and GAUSSIAN chunk filling functions works
add_executable(test_fill_chunk_src_with_pointgauss_app
    test_fill_chunk_src_with_pointgauss.c
    fill_chunks_common.c
)
target_link_libraries(test_fill_chunk_src_with_pointgauss_app
    Unity
    chunk_sky_model
    ${CC_LINKLIBS}
)
add_test(C_test_fill_chunk_src_with_pointgauss test_fill_chunk_src_with_pointgauss_app)

## Test the SHAPELET chunk filling functions works
add_executable(test_fill_chunk_src_with_shapelets_app
    test_fill_chunk_src_with_shapelets.c
    fill_chunks_common.c
)
target_link_libraries(test_fill_chunk_src_with_shapelets_app
    Unity
    chunk_sky_model
    ${CC_LINKLIBS}
)
add_test(C_test_fill_chunk_src_with_shapelets test_fill_chunk_src_with_shapelets_app)

# ## Test the overall chunking function that wraps all this functionality works
add_executable(test_create_chunked_sky_models_app
    test_create_chunked_sky_models.c
    fill_chunks_common.c
)
target_link_libraries(test_create_chunked_sky_models_app
    Unity
    chunk_sky_model
    ${CC_LINKLIBS}
)
add_test(C_test_create_chunked_sky_models test_create_chunked_sky_models_app)
