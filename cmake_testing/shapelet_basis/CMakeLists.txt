add_library(shapelet_basis SHARED "${CMAKE_SOURCE_DIR}/src/shapelet_basis.c")
## Accumulate the coverage from the tests
add_coverage(shapelet_basis)

# ##Test the shapelet basis array creation works
add_executable(test_create_sbf_app
    test_create_sbf.c
)
target_link_libraries(test_create_sbf_app
    Unity
    shapelet_basis
    ${CC_LINKLIBS}
)
add_test(C_test_create_sbf test_create_sbf_app)
