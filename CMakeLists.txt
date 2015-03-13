cmake_minimum_required(VERSION 2.8)
project (SciPAL)
SET (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

find_package(OpenBlas REQUIRED)

find_package(CUDA QUIET REQUIRED)
set(CUDA_HOST_COMPILATION_CPP ON)
set(CUDA_PROPAGATE_HOST_FLAGS OFF)
set(CUDA_SEPARABLE_COMPILATION ON)
set(CUDA_NVCC_FLAGS --ptxas-options=-v -std=c++11 -g -Xcompiler -fexceptions -Xcompiler -fPIC -Xcompiler -fopenmp)
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-O3 -gencode arch=compute_20,code=sm_20 -gencode arch=compute_35,code=sm_35)
set(CMAKE_C_FLAGS "-g -fPIC -std=c99")
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-g -fPIC -std=c++11 -fopenmp")

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}" CACHE STRING "c++ flags")
SET (CMAKE_EXPORT_COMPILE_COMMANDS 1)
set (CMAKE_EXPORT_COMPILE_COMMANDS "${CMAKE_EXPORT_COMPILE_COMMANDS}" CACHE STRING "export compile_commands.json")

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)
set (CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING "build type")

FIND_PACKAGE(Qt4 COMPONENTS QtCore REQUIRED)
INCLUDE(${QT_USE_FILE})

include_directories(
  ${SciPAL_SOURCE_DIR}
  ${SciPAL_SOURCE_DIR}/include
  /usr/local/include
)

link_directories(
  /usr/local/lib
)

add_subdirectory(testsite/step-42)
