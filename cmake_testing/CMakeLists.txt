cmake_minimum_required(VERSION 3.10)

##Need to find CUDA compiler
include(CheckLanguage)
check_language(CUDA)

project(woden VERSION 0.1 LANGUAGES CUDA CXX C)

##Gotta find CUDA and get it to find libraries
find_package(CUDA)
if(CUDA_FOUND)
  message(STATUS "Found CUDA ${CUDA_VERSION_STRING} at ${CUDA_TOOLKIT_ROOT_DIR}")
  message(STATUS "CUDA library found; building GPU-enabled executables.")
  ##Sets the -std=c++11 nvcc compiler flag I think
  set(CMAKE_CUDA_STANDARD 11)
  #set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -std=c++11 ")
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)
endif()

# # file(GLOB CUDASOURCES "../src/*.cu")
# file(GLOB CUDASOURCES "../src/fundamental_coords.cu")
#
# add_library(wodenCUDA SHARED ${CUDASOURCES})

set(TARGET_GROUP production CACHE STRING "Group to build")

include_directories("${CMAKE_SOURCE_DIR}/../include/")


set(HINTS_ERFA_INC "/usr/local/include/;${ERFA_ROOT}/include/;${ERFA_ROOT}")
set(HINTS_ERFA_LIB "/usr/local/lib/;${ERFA_ROOT}/src/.libs/;${ERFA_ROOT}.libs/;${ERFA_ROOT};${ERFA_ROOT}/lib/")

find_path(ERFA_INC erfa.h HINTS ${HINTS_ERFA_INC})
find_library(ERFA_LIB erfa HINTS ${HINTS_ERFA_LIB} REQUIRED)
if(ERFA_LIB)
  message(STATUS "ERFA includes path: ${ERFA_INC}")
  message(STATUS "ERFA library path: ${ERFA_LIB}")
endif()

set(CC_LINKLIBS "-pthread -lm ${ERFA_LIB}")

if(TARGET_GROUP STREQUAL production)
    # add_subdirectory(main)
elseif(TARGET_GROUP STREQUAL test)

    include(CTest)

    set(HINTS_UNITY "${UNITY_ROOT};${UNITY_ROOT}/src/;/usr/local/Unity/src/")

    # set(HINTS_UNITY "Unity/src")

    find_path(UNITY_DIR unity.h HINTS ${HINTS_UNITY})

    if(UNITY_DIR)
      message(STATUS "UNITY src path: ${UNITY_DIR}")
    endif()

    add_library(Unity STATIC
        "${UNITY_DIR}/unity.c"
    )

    target_include_directories(Unity PUBLIC
        "${UNITY_DIR}"
    )

    add_subdirectory(unit_tests)
else()
    message(FATAL_ERROR "Given TARGET_GROUP unknown")
endif()
