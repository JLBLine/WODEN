# ##Test the COMPONENT nulling function
# add_executable(test_null_comps_app
#     test_null_comps.c
#     ../../src/chunk_sky_model.c
# )
# target_link_libraries(test_null_comps_app
#     Unity
#     ${CC_LINKLIBS}
# )
# add_test(C_test_null_comps test_null_comps_app)

# ##Test the POINT and GAUSSIAN chunk filling functions works
add_executable(test_fill_chunk_src_with_pointgauss_app
    test_fill_chunk_src_with_pointgauss.c
    ../../src/chunk_sky_model.c
)
target_link_libraries(test_fill_chunk_src_with_pointgauss_app
    Unity
    ${CC_LINKLIBS}
)
add_test(C_test_fill_chunk_src_with_pointgauss test_fill_chunk_src_with_pointgauss_app)
