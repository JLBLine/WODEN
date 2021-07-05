##Test C code
add_subdirectory(array_layout)
add_subdirectory(chunk_sky_model)
add_subdirectory(create_sky_model)
add_subdirectory(FEE_primary_beam)
add_subdirectory(primary_beam)
add_subdirectory(visibility_set)
add_subdirectory(woden_settings)
add_subdirectory(shapelet_basis)

##Test CUDA code
add_subdirectory(calculate_visibilities)
add_subdirectory(fundamental_coords)
add_subdirectory(FEE_primary_beam_cuda)
add_subdirectory(primary_beam_cuda)
add_subdirectory(source_components)

##Test python code
add_subdirectory(run_woden)
