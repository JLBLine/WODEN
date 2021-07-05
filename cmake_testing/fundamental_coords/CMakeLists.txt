file(GLOB CUDASOURCES "../../src/fundamental_coords.cu")

add_library(fundamental_coordsCUDA SHARED ${CUDASOURCES})

## For CMake versions > 3.18, the arch flag settings are handled by
## CUDA_ARCHITECTURES. For now, just switch this warning off - in future,
## replace the CUDA_FLAG lines above with something like
set_property(TARGET fundamental_coordsCUDA PROPERTY CUDA_ARCHITECTURES OFF)

##Test l,m,n functions
add_executable(test_lmn_coords_app
    test_lmn_coords.c
)
target_link_libraries(test_lmn_coords_app
    Unity
    fundamental_coordsCUDA
)
add_test(CUDA_test_lmn_coords test_lmn_coords_app)



##Test u,v,w functions
add_executable(test_uvw_coords_app
    test_uvw_coords.c
)
target_link_libraries(test_uvw_coords_app
    Unity
    fundamental_coordsCUDA
    ${CC_LINKLIBS}
)
add_test(CUDA_test_uvw_coords test_uvw_coords_app)
