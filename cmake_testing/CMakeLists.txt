##Define some constants that are used in all of the testing
set(FLOAT "float")
set(DOUBLE "double")

##Turns on code coverage type outputs for gcc
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-arcs -ftest-coverage")

##Flags to make C code compile in WODEN float or double modes
set(FLOAT_C_FLAGS -Wall -g )
set(DOUBLE_C_FLAGS -Wall -g -DDOUBLE_PRECISION)

##Flags to make C code compile in WODEN float or double modes,
##and to create code coverage reports. Only need these flags on the things
##we want to cover with tracking
set(FLOAT_C_COVER_FLAGS -coverage -g -Wall)
set(DOUBLE_C_COVER_FLAGS -coverage -g -Wall -DDOUBLE_PRECISION)

##Flags to make CUDA code compile in WODEN float or double modes
set(FLOAT_CUDA_FLAGS -g -G --compiler-options -Wall )
set(DOUBLE_CUDA_FLAGS -DDOUBLE_PRECISION -g -G --compiler-options -Wall)

##define a function that switches between flags
function(DEFINE_COMP_FLAGS PRECISION C_FLAGS CUDA_FLAGS C_COVER_FLAGS)
  if(${PRECISION} MATCHES ${FLOAT})
    set(${C_FLAGS} ${FLOAT_C_FLAGS} PARENT_SCOPE)
    set(${CUDA_FLAGS} ${FLOAT_CUDA_FLAGS} PARENT_SCOPE)
    set(${C_COVER_FLAGS} ${FLOAT_C_COVER_FLAGS} PARENT_SCOPE)
  elseif(${PRECISION} MATCHES ${DOUBLE}) # optional block, can be repeated
    set(${C_FLAGS} ${DOUBLE_C_FLAGS} PARENT_SCOPE)
    set(${CUDA_FLAGS} ${DOUBLE_CUDA_FLAGS} PARENT_SCOPE)
    set(${C_COVER_FLAGS} ${DOUBLE_C_COVER_FLAGS} PARENT_SCOPE)
  endif()
endfunction()

##Test C code
add_subdirectory(array_layout)
add_subdirectory(chunk_sky_model)
add_subdirectory(create_sky_model)
add_subdirectory(primary_beam)
add_subdirectory(read_text_skymodel)
add_subdirectory(read_yaml_skymodel)
add_subdirectory(visibility_set)
add_subdirectory(shapelet_basis)
add_subdirectory(woden_settings)


###Test CUDA code
add_subdirectory(calculate_visibilities)
add_subdirectory(fundamental_coords)
add_subdirectory(primary_beam_cuda)
add_subdirectory(source_components)

####Test python code
add_subdirectory(run_woden)
add_subdirectory(add_woden_uvfits)
add_subdirectory(concat_woden_uvfits)
