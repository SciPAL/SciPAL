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

# helper file for cleaning up the qmake config
# this files contains standard libs / includes and compiler config
# for the deal.II-MPI part of the example programs.

#DEFINES += __cplusplus

QMAKE_CXXFLAGS += -std=gnu++0x

# Qt considers OSX as a unix.
    macx {
DEALHOME =/clusterfs/cuda-2015/dealii-minimal-installation
           }
    else {
        unix {
DEALHOME =/clusterfs/cuda-2015/dealii-minimal-installation
        }
    }
#includes for deal.II
INCLUDEPATH += $$DEALHOME/include \
          /opt/local/ \ # this is boost's home on mac OSX its 1.52, check boost/version.hpp
          $$DEALHOME/include/deal.II/bundled/ \ # this is for the TBB headers
		  $$DEALHOME/bundled/boost-1.49.0/include/

#deal.II LIBS
LIBS += -L$$DEALHOME/lib

CONFIG(debug, debug|release) {
  LIBS += -ldeal_II.g
    # Compiler switches for debugging mode
    DEFINES += DEBUG
}

CONFIG(release, debug|release) {
 LIBS += -ldeal_II
}

# blas, lapack
LIBS += -lblas -llapack -lz -ldl

## tbb
#CONFIG(debug, debug|release) {
#LIBS += -ltbb_debug
#}

#CONFIG(release, debug|release) {
#LIBS +=  -ltbb
#}



