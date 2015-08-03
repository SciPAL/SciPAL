project(XRAY_DECONV C CXX)
# CMake
cmake_minimum_required(VERSION 2.8)

# support external builds
include_directories(${CMAKE_CURRENT_BINARY_DIR}) 
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/..)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../SciPAL/include)
#set mpic++ as compiler
#set(CMAKE_C_COMPILER mpicc)
#set(CMAKE_CXX_COMPILER mpicxx)
# Boost
find_package(Boost COMPONENTS serialization date_time REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})


find_library(FFTW_LIBRARY 
    NAMES fftw3 fftw)
set(FFTW_LIBRARIES "${FFTW_LIBRARY}")
find_library(FFTW_libm_LIBRARY
    NAMES m)
list(APPEND FFTW_LIBRARIES "${FFTW_libm_LIBRARY}")

#target_link_libraries(main ${FFTW_LIBRARIES})
#list(APPEND LocalisatorLibraries
#    ${FFTW_LIBRARIES}
#)

# Qt
#find_package(Qt4 REQUIRED)
#set(QT_USE_QTMAIN true)
#set(QT_USE_QTNETWORK true)
#include(${QT_USE_FILE})

# CUDA
find_package(CUDA)
include_directories(${CUDA_TOOLKIT_INCLUDE})
if (CUDA_FOUND)
	message(${CUDA_TOOLKIT_INCLUDE})
else()
	message("CUDA not found!")
endif()
set(CUDA_NVCC_FLAGS "-O3 --gpu-architecture sm_35 -Xcompiler -fopenmp -pthreads")

#cuda_add_library(loc cuda_kernel_step-35.cu)

cuda_add_executable(loc step-35.cpp)
#
# OpenMP
find_package(OpenMP)
#set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS};-ftz=true;-prec-div=false;-prec-sqrt=false;-gencode;arch=compute_30,code=sm_30;-gencode;arch=compute_20,code=sm_20;-gencode;arch=compute_13,code=sm_13;-gencode;arch=compute_11,code=sm_11")
#if(NOT THRUST_BASE)
#  set( THRUST_BASE ${CUDA_TOOLKIT_INCLUDE} CACHE PATH "Thrust library path" )
#endif()
#find_path(THRUST_INCLUDE_DIR thrust/version.h
#	PATHS ${THRUST_BASE}
#NO_DEFAULT_PATH
#)
#include_directories(${THRUST_INCLUDE_DIR})

#find_package(TIFF REQUIRED)
set(TIFF_INCLUDE_DIR /scratch/step-35/lib/include)
include_directories(${TIFF_INCLUDE_DIR})

list(APPEND LocalisatorSources
		step-35.cpp
)

list(APPEND LocalisatorHeaders
            cuda_driver_step-35.h 
            cuda_driver_step-35.hh 
            cuda_kernel_wrapper_step-35.cu.h 
            SimParams.h 
            step-35.hh
)

#list(APPEND LocalisatorMocs
#)

#qt4_wrap_cpp(LocalisatorMetaSource ${LocalisatorMocs})

#list(APPEND LocalisatorUIs
#)

#qt4_wrap_ui(LocalisatorUIHeaders ${LocalisatorUIs})
set(TIFF_LIBRARY /scratch/step-35/lib/lib/libtiff.so.5.2.0)
list(APPEND LocalisatorLibraries
	${TIFF_LIBRARY}
)
list(APPEND LocalisatorLibraries
    ${FFTW_LIBRARIES}
    )
#list(APPEND LocalisatorLibraries mpi)
add_executable(localisator
	${LocalisatorSources}
	${LocalisatorHeaders}
	${LocalisatorMetaSources}
	${LocalisatorUIHeaders}
)
target_link_libraries(localisator
	${LocalisatorLibraries}
)
##set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)
##'--std=c++11' is the better compiler flag, '--std=c++0x -pthread' is a fix for c4/5 (old gcc, ubuntu)
set(CMAKE_CXX_FLAGS_RELEASE "-O2 --std=c++0x -pthread ${OpenMP_CXX_FLAGS} -lgomp")
#set(CMAKE_CXX_FLAGS_RELEASE "-O2 Xcompiler -fopenmp -pthread")
