# helper file for cleaning up the qmake config
# this files contains standard libs / includes and compiler config
# for the deal.II-MPI part of the example programs.

DEFINES += __cplusplus

QMAKE_CXXFLAGS += -std=gnu++0x


#includes for deal.II
INCLUDEPATH += $$DEALHOME/include \
               $$DEALHOME/contrib/boost-1.46.1/include \
                /usr/lib/openmpi/include  \
                /usr/lib/openmpi/include/openmpi \
                /usr/local/share/petsc-3.2-p7/include  \
                /usr/local/share/petsc-3.2-p7/arch-linux2-c-debug/include  \
                /usr/local/share/trilinos-10.4.2/include \
                $$DEALHOME/contrib/tbb/tbb*/include \
                /usr/local/share/p4est/DEBUG/include/


# MPI libs
LIBS += -L/usr/lib \
        -lmpi_cxx \
        -lmpi
#deal.II LIBS
LIBS += -L$$DEALHOME/lib

CONFIG(debug, debug|release) {
    LIBS += -lpetscall.g  -ldeal_II.g
    # Compiler switches for debugging mode
    DEFINES += DEBUG
}

CONFIG(release, debug|release) {
    LIBS +=  -lpetscall -ldeal_II
}

# blas, lapack
LIBS += -lblas -llapack -lz -ldl

# trilinos
LIBS += -L/usr/local/share/trilinos-10.4.2/lib \
        -lstratimikos -lstratimikosbelos -lstratimikosaztecoo -lstratimikosamesos -lstratimikosml \
        -lstratimikosifpack \
        -lbelostpetra -lbelosepetra -lbelos -lml -lifpack -lamesos -lgaleri -laztecoo -lisorropia \
        -lthyratpetra -lthyraepetraext -lthyraepetra -lthyra -lepetraext -ltpetrainout -ltpetra \
        -ltriutils -lzoltan -lepetra -lkokkoslinalg \
        -lkokkosnodeapi -lkokkos -lrtop -lsacado -ltpi -lteuchos \
        /usr/lib/libnetcdf.so.6 /usr/lib/libnetcdf_c++.so.5

#petsc
LIBS += -L/usr/local/share/petsc-3.2-p7/arch-linux2-c-debug/lib -lpetsc

# metis
LIBS += -L/usr/local/share/metis-5.0.2/lib -lmetis

## tbb
#CONFIG(debug, debug|release) {
#LIBS += -ltbb_debug
#}

#CONFIG(release, debug|release) {
#LIBS +=  -ltbb
#}



