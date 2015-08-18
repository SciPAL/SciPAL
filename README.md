## This is the SciPAL readme file.

In order to use SciPAL you need the following dependencies:

* deal.II (www.dealii.org), see additional information in their tar-archive.
  Running deal.II's cmake-file should give you a working installation

* CUDA-Toolkit (https://developer.nvidia.com/cuda-downloads). For further
  information have a look in NVidia's getting started guide (also on link).

* CMake

* A compiler that is able to make use of OpenMP (gcc works fine, clang not)

* Python2 (necessary for the various scripts)

* QtCore (libqt4-core qt4-dev-tools)

* OpenBLAS (libopenblas-base libopenblas-dev)

* Boost (http://www.boost.org/)

* Intel TBB (https://www.threadingbuildingblocks.org/)

If you have these installed, put SciPAL in your preferred folder and run the
installation script install.py from /scripts/. This will set the paths of the
CUDA and deal.II installation.

The best way to get this running is probably to use Qt Creator as IDE, or
compile the project by hand with CMake.

Note: The installation script configures only for a basic deal.II installation.
      If you have a deal.II installation using MPI, petsc, etc. you can use the
      advanced configuration file from /config/dealii_mpi_conf.pro as template
      for your configuration.

Furhter notes on deal.II configuration:
* LAPACK, BLAS needed
* it shouldn't have MPI, Trilmos, etc. All these needs to be turned off.

deal.II can be found in /clusterfs/cuda-2015/deal.II on the computers C4 to C6 in the numerics dept. 
A working version of SciPAL will also be there soon.
