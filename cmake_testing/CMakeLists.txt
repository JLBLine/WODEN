##Define some constants that are used in all of the testing
set(FLOAT "float")
set(DOUBLE "double")
set(GPU "gpu")
set(CPU "cpu")

##Turns on code coverage type outputs for gcc
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-arcs -ftest-coverage")

##Flags to make C code compile in WODEN float or double modes
set(FLOAT_C_FLAGS -Wall -g -lm )
# set(DOUBLE_C_FLAGS -O3 -march=native -Wall -g -lm -DDOUBLE_PRECISION)
set(DOUBLE_C_FLAGS -Wall -g -lm -DDOUBLE_PRECISION)

##Flags to make C code compile in WODEN float or double modes,
##and to create code coverage reports. Only need these flags on the things
##we want to cover with tracking
set(FLOAT_C_COVER_FLAGS -coverage -g -lm -Wall)
set(DOUBLE_C_COVER_FLAGS -coverage -g -lm -Wall -DDOUBLE_PRECISION)

if(USE_CUDA)
  ##Flags to make CUDA code compile in WODEN float or double modes
  set(FLOAT_GPU_FLAGS -g -G --compiler-options -Wall )
  set(DOUBLE_GPU_FLAGS -DDOUBLE_PRECISION -g -G --compiler-options -Wall -D__NVCC__ -rdc=true)
elseif(USE_HIP)
  ##Flags to make CUDA code compile in WODEN float or double modes
  set(FLOAT_GPU_FLAGS -g -Wall )
  set(DOUBLE_GPU_FLAGS -DDOUBLE_PRECISION -g -Wall -D__HIPCC__)
endif()



##define a function that switches between flags
function(DEFINE_COMP_FLAGS PRECISION C_FLAGS GPU_FLAGS C_COVER_FLAGS)
  if(${PRECISION} MATCHES ${FLOAT})
    set(${C_FLAGS} ${FLOAT_C_FLAGS} PARENT_SCOPE)
    set(${GPU_FLAGS} ${FLOAT_GPU_FLAGS} PARENT_SCOPE)
    set(${C_COVER_FLAGS} ${FLOAT_C_COVER_FLAGS} PARENT_SCOPE)
  elseif(${PRECISION} MATCHES ${DOUBLE}) # optional block, can be repeated
    set(${C_FLAGS} ${DOUBLE_C_FLAGS} PARENT_SCOPE)
    set(${GPU_FLAGS} ${DOUBLE_GPU_FLAGS} PARENT_SCOPE)
    set(${C_COVER_FLAGS} ${DOUBLE_C_COVER_FLAGS} PARENT_SCOPE)
  endif()
endfunction()

##make sure we have Python. Deal with differing versions of CMake causing
##slightly different definitions of PYTHON_EXECUTABLE
##Why did they change PYTHON_EXECUTABLE to Python_EXECUTABLE??

include(FindPython)

if (NOT DEFINED ${PYTHON_EXECUTABLE})
    set(PYTHON_EXECUTABLE ${Python_EXECUTABLE})
endif()

## Test GPU/C code
add_subdirectory(GPU_or_C_code)

# ## Test python library
add_subdirectory(wodenpy)

# ## Test python scripts
add_subdirectory(scripts)
