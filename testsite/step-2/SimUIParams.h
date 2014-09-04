/*This file is part of SciPAL.

    SciPAL is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SciPAL is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with SciPAL.  If not, see <http://www.gnu.org/licenses/>.

Copyright  S. C. Kramer , J. Hagemann  2010 - 2014
*/


#ifndef SimUIParams_H
#define SimUIParams_H

// Before we discuss the implementation of the benchmarking we briefly go through the
// parameters of the program which can be set at runtime.
//
// The runtime parameters are grouped into several
// simple structures. Their main purpose is to encapsulate the
// declare() and get() functions needed to load the parameters by means
// of the deal.II ParameterHandler class. These functions are all
// very similar to each other
// so that we do not present all of their definitions in this documentation.
// To understand the program it suffices to know which parameters
// have been declared. Hence, we merely include the definitions
// of the parameter structures.

#include <QString>
#include <deal.II/base/parameter_handler.h>

namespace step2 {

// @sect3{Enum: MVCase}
//
// These enumerated values will be used as tags for the different implementations of the matrix-vector product.
enum MVCase { Fujimoto_mv, atlas_mv, cublas_mv, none /* for future use : , openmp_mv, ...*/ };

// @sect3{struct: SimUIParams}
//
// This structure contains all parameters which are needed for each test.
// These are essentially some flags which indicate whether the MV-test
// has to be executed for a certain combination of numerical precision and
// BLAS-implementation.
struct SimUIParams {

    SimUIParams() {}

    static void declare(dealii::ParameterHandler & prm);

    void get(dealii::ParameterHandler & prm);

    // The attributes consist of a long list of flags
    // indicating which combinations of tests are to be run ...
    bool
    run_cpublas_vs_cublas_float,
    run_cpublas_vs_cublas_double,
    run_cpublas_vs_Fujimoto_float,
    run_cpublas_vs_Fujimoto_double,
    run_Fujimoto_vs_cublas_float,
    run_Fujimoto_vs_cublas_double;

    int fj_version;

    // and a variable for holding the run directory.
    QString run_dir;

    // The test runs are subdivided into groups for
    // single and double precision.
    // The list of matrix-vector product implementations which are to be executed are stored as vector.
    std::vector<MVCase> float_tests, double_tests;
    // For the comparisons of the runtimes we need the reverse, i.e. a mapping
    // from the name of the implementation to its index in the list of tests.
    std::map<MVCase, int> float_vs, double_vs;

private:
    SimUIParams (const SimUIParams & /*other*/) {}

    SimUIParams & operator = (const SimUIParams & /*other*/) { return *this; }
};

}
#endif // SimUIParams_H
