# ##Test the shapelet basis array creation works
add_executable(test_create_sbf_app
    test_create_sbf.c
    ../../src/shapelet_basis.c
)
target_link_libraries(test_create_sbf_app
    Unity
    ${CC_LINKLIBS}
)
add_test(C_test_create_sbf test_create_sbf_app)
