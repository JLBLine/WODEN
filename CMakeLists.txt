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

##Pull out arch version if specified by user, and pass to nvcc
if (CUDA_ARCH EQUAL "1.3")
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_13")
elseif (CUDA_ARCH EQUAL "2.0")
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_20")
elseif (CUDA_ARCH EQUAL "2.1")
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_21")
elseif (CUDA_ARCH EQUAL "3.0")
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_30")
elseif (CUDA_ARCH EQUAL "3.2")
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_32")
elseif (CUDA_ARCH EQUAL "3.5")
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_35")
elseif (CUDA_ARCH EQUAL "3.7")
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_37")
elseif (CUDA_ARCH EQUAL "5.0")
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_50")
elseif (CUDA_ARCH EQUAL "5.2")
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_52")
elseif (CUDA_ARCH EQUAL "6.0")
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_60")
elseif (CUDA_ARCH EQUAL "6.1")
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_61")
elseif (CUDA_ARCH EQUAL "6.2")
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_62")
elseif (CUDA_ARCH EQUAL "7.0")
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_70")
elseif (CUDA_ARCH EQUAL "7.5")
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_75")
endif()

##Helps with version numbering
configure_file(wodenConfig.h.in wodenConfig.h)

##All the includes dirs are in include
#include_directories(include)

##Find the source codes
file(GLOB CSOURCES "src/*.c")
file(GLOB CUDASOURCES "src/*.cu")


##Find the following libraries
##Use can specify CFITSIO_ROOT and JSONC_LIB when running Cmake
set(HINTS_CFITSIO_LIB "${CFITSIO_ROOT}/lib;${CFITSIO_ROOT}")
set(HINTS_CFITSIO_INC "${CFITSIO_ROOT}/include;${CFITSIO_ROOT}")

find_path(CFITSIO_INC fitsio.h HINTS ${HINTS_CFITSIO_INC})
find_library(CFITSIO_LIB cfitsio HINTS ${HINTS_CFITSIO_LIB} REQUIRED)
if(CFITSIO_LIB)
  message(STATUS "CFITSIO include path: ${CFITSIO_INC}")
  message(STATUS "CFITSIO library path: ${CFITSIO_LIB}")
endif()


set(HINTS_JSONC_LIB "${JSONC_ROOT}/lib;${JSONC_ROOT}/lib64;${JSONC_ROOT}")
set(HINTS_JSONC_INC "${JSONC_ROOT}/include;${JSONC_ROOT}/include/json-c;${JSONC_ROOT};/usr/include/json-c/")

find_path(JSONC_INC json.h HINTS ${HINTS_JSONC_INC})
find_library(JSONC_LIB json-c HINTS ${HINTS_JSONC_LIB} REQUIRED)
if(JSONC_LIB)
  message(STATUS "JSONC includes path: ${JSONC_INC}")
  message(STATUS "JSONC library path: ${JSONC_LIB}")
endif()

##Set the linkers for the C code
set(CC_LINKLIBS "-lm -pthread ${CFITSIO_LIB} ${JSONC_LIB}")

##All the WODEN header dirs are in include dir
include_directories(include ${JSONC_INC} ${CFITSIO_INC})

##Generate libraries from the C and CUDA source code
add_library(wodenC SHARED ${CSOURCES})
add_library(wodenCUDA SHARED ${CUDASOURCES})

##make the final WODEN executable, linking the libraries
add_executable(woden src/woden.c)
target_link_libraries(woden wodenC wodenCUDA ${CC_LINKLIBS})
