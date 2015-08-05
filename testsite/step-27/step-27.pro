# rules for host part
# Remove any QT-dependency
# QT =
# Due to QThread we need at least
#QT = core

# maybe test this:
QT+= core #gui multimedia multimediawidgets
#greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

# thread-safe console application with rules for debug and release version
CONFIG += console \
    thread \
    debug_and_release #\
    #qt

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
DEFINES += "NUM_THREADS=512"
DEFINES += "BLOCK_DIM_Y=4"


    HOME = $$(HOME) # your home directory
    PRAK = $$_PRO_FILE_PWD_/../../ #$$HOME/cuda-2014/Praktikum_2013 # path to your copy of the lab course folder. Typically, this 2 levels above the source folder of your step-
    SciPAL_DIR = $$PRAK
    STEP_PARENT_DIR = $$_PRO_FILE_PWD_/..

message("SciPALs home :" $$SciPAL_DIR)
message("step home :" $$STEP_DIR)

# Qt considers OSX as a unix.
#    macx {
#         DEALHOME = /Library/deal.II-bundle/deal.II-8.1 # /usr/local #
#           }
#    else {
#        unix {
#         DEALHOME = /usr/local/deal.II-7.2.0 # path to deal II in NAM
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
        CUDA_INCLUDES = -I$$STEP_PARENT_DIR


        CUDA_INCLUDES +=-I$$SciPAL_DIR/include


    # and here for the gcc
INCLUDEPATH += ..

#put here your non-standard libs
LIBS +=

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
         include($$SciPAL_DIR/config/dealii_mpi_conf.pro)
        }
    }

CUDA_CFLAGS += -rdc=true

# Enter project specific source and header files here
SOURCES = \
    step-27.cpp \
    SimParams.cpp \
    NumericsParams.cpp \
    PhysicalParams.cpp



HEADERS = \
    SimParams.h \
    Architecture.h \
    drs_cell.hh \
    UnitTestProblemData.hh \
    BoundaryDoFTools.hh \
    drs_simulation.hh \
    assemble_system.hh \
    bemforms.h \
    simulation_handler.h \
    timing.h \
    geomgpc.h \
    bemforms.hh \
    coupling_type.h \
    NumericsParams.h \
    PhysicalParams.h \
    drs_lac.hh \
    LaplaceForm.hh


   # the following variable contains files which should appear in the Projects view on the left of QtCreator
   # which are not subject to compilation.
OTHER_FILES = doxygen_filelist\
                doc/*.dox \
                ../step-27-build/prm/step-27.prm \
                ../step-27-build/prm/step-27_drs_numerics.prm \
                prm/*.prm


