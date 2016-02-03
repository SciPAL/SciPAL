#This file is part of SciPAL.

#    SciPAL is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Lesser General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    SciPAL is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Lesser General Public License for more details.

#    You should have received a copy of the GNU Lesser General Public License
#    along with SciPAL.  If not, see <http://www.gnu.org/licenses/>.

#Copyright  S. C. Kramer , J. Hagemann  2010 - 2014



# rules for host part
# Remove any QT-dependency
# QT =
# Due to QThread we need at least
# QT = core
# thread-safe console application with rules for debug and release version



CONFIG += console \
    thread \
    debug_and_release# \
    #qt

# The following keeps the compiler from messing up the different signals&slots implementations from QT and boost.
CONFIG += no_keywords

# For the time being we have to disable the usage of the CUDATimer class on Mac OSX as it leads to linker errors.
# This also shows how to configure different compilation behaviors on different OSes.
    macx {
         DEFINES += DONT_USE_CUDATIMER
           }
    else {
        unix {

        }
    }
    win32 {

    }

# the following disables the complex number type provided by CUDA.
DEFINES += nUSE_CUDA_COMPLEX_VERSION
DEFINES += USE_CPP11




    HOME = $$(HOME) # your home directory
    PRAK = $$_PRO_FILE_PWD_/../../ #$$HOME/cuda-2014/Praktikum_2013 # path to your copy of the lab course folder. Typically, this 2 levels above the source folder of your step-
    SciPAL_DIR = $$PRAK
    STEP_PARENT_DIR = $$_PRO_FILE_PWD_/..
    STEP_DIR = $$_PRO_FILE_PWD_
    DEALHOME = /usr/local/dealii

message("SciPALs home :" $$SciPAL_DIR)
message("step home :" $$STEP_DIR)

# Qt considers OSX as a unix.
#    macx {
#         DEALHOME = /usr/local #
#           }
#    else {
#        unix {
#         DEALHOME = /usr/local/deal.II-8.1 # path to deal II in NAM
#        }
#    }


     #put here some non standard header includes for cuda
     # which are specific to your project.
     # This has to be done before the generic part of the configuration is added.
     # -I. is needed if you compiled within your source directory.
     # The argument -I../../.. is needed to go from the shadow-build directory
     # to the source directory of your project and assumes that the shadow build has three levels, e.g.
     #
     # build-step-2/Debug/Desktop_Qt_5_2_1_clang64/
     #
        CUDA_INCLUDES += -I$$STEP_PARENT_DIR \
                        -I$$STEP_DIR/cufftShift/Src/


        CUDA_INCLUDES +=-I$$SciPAL_DIR/include

#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scratch/step-35/lib/include
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scratch/step-35/lib/lib
#export BOOST_LIBRARYDIR=/scratch/step-35/lib/lib
#export TIFF_LIBRARY=/usr/lib/x86_64-linux-gnu
#export TIFF_INCLUDE_DIR=/usr/include/x86_64-linux-gnu
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/dealII-8.0/include


    # and here for the gcc
INCLUDEPATH += .. . \
                ./cufftShift/Src/
#INCLUDEPATH += /scratch/step-35/lib/include

#put here your non-standard libs

  macx {
INCLUDEPATH += /opt/local/include
INCLUDEPATH += /Library/deal.II-bundle/deal.II-8.1-serial/include
INCLUDEPATH += /Library/deal.II-bundle/deal.II-8.1-serial/include/deal.II/bundled
LIBS += /opt/local/lib/libtiff.dylib

} else {
DEFINES += USE_OMP
INCLUDEPATH += /usr/local/dealii/include
INCLUDEPATH += /usr/include/x86_64-linux-gnu


LIBS += /usr/lib/x86_64-linux-gnu/libtiff.so
LIBS += /usr/lib/x86_64-linux-gnu/libtiffxx.so
}

#include(qtlibs)
#LIBS += -lfftw
   # The PRAK variable is used inside scipal_conf.pro.
   # This file provides a basic configuration for compiling
   # a cuda project which depends on CUDA, deal.II and SciPAL.
   # It defines the rules for compiling the CUDA part. In this definition
   # the variable CUDA_INCLUDES is evaluated. Therefore, we have to add project-specific
   # includes BEFORE we include the general configuration.
include($$SciPAL_DIR/config/scipal_conf.pro)

# Depending an what is needed from deal.II we either have to include a full MPI-capable installation or the quick-and-easy default version
# which only provides mltithreading for shared memory machines.
    macx {
message("Load deal.II simple config")
         include($$SciPAL_DIR/config/dealii_simple_conf.pro) #
           }
    else {
        unix {
message("Load deal.II MPI config")
         include($$SciPAL_DIR/config/dealii_simple_conf.pro)
        }
    }



# Enter project specific source and header files here
SOURCES += \
    step-35.cpp \
    cuda_kernel_step-35.cu \
    ./cufftShift/Src/CUDA/Interfaces/in-place/cufftShift_1D_IP_impl.cu \
    ./cufftShift/Src/CUDA/Kernels/in-place/cufftShift_1D_IP.cu \
    ./cufftShift/Src/CUDA/Interfaces/in-place/cufftShift_2D_IP_impl.cu \
    ./cufftShift/Src/CUDA/Kernels/in-place/cufftShift_2D_IP.cu \
    ./cufftShift/Src/CUDA/Interfaces/in-place/cufftShift_3D_IP_impl.cu \
    ./cufftShift/Src/CUDA/Kernels/in-place/cufftShift_3D_IP.cu \
    ./cufftShift/Src/CUDA/Interfaces/out-of-place/cufftShift_1D_OP_impl.cu \
    ./cufftShift/Src/CUDA/Kernels/out-of-place/cufftShift_1D_OP.cu \
    ./cufftShift/Src/CUDA/Interfaces/out-of-place/cufftShift_2D_OP_impl.cu \
    ./cufftShift/Src/CUDA/Kernels/out-of-place/cufftShift_2D_OP.cu \
    ./cufftShift/Src/CUDA/Interfaces/out-of-place/cufftShift_3D_OP_impl.cu \
    ./cufftShift/Src/CUDA/Kernels/out-of-place/cufftShift_3D_OP.cu \
    ./cufftShift/Src/CXX/cufftShift_1D.cpp \
    ./cufftShift/Src/CXX/cufftShift_2D.cpp\
    ./cufftShift/Src/CXX/cufftShift_3D.cpp\
    ./cufftShift/Src/cuUtils/configGPU.cpp\
    ./cufftShift/Src/cxxUtils/PrintMemory.cpp
 #\
    #$$SciPAL_DIR/include/numerics/propagation_kernels.cu

HEADERS += \
    cuda_driver_step-35.h \
    cuda_driver_step-35.hh \
    cuda_kernel_wrapper_step-35.cu.h \
    cuda_kernel_step-35.cu.c \
    patch.h \
    cuda_helper.h \
    ADMMParams.h \
    preprocessor_directives.h \
    extremeValueStatisticsGenerator.h \
    step-35.hh \
    smre_problem.hh \
    *.h* \
    *.cu.c

   # the following variable contains files which should appear in the Projects view on the left of QtCreator
   # which are not subject to compilation.
OTHER_FILES = doxygen_filelist \
                doc/*.dox \
                prm/*.prm





















