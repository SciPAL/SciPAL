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


#ifndef MVMultDriverInterface_H
#define MVMultDriverInterface_H

#include <vector>
#include <numeric>
#include <lac/blas++.h>

// To avoid boring code we declare the member functions
// of the classes defined in this file inline.
namespace step2 {

typedef    std::pair<size_t, size_t> matrixSize;

// @sect3{Class: MVMultDriverInterface}
//
// This is the common base class for all matrix-vector
// multiplications tested in this program. We define the
// type of dense matrix which is to be used on the host side in this
// class. As this is a template class w.r.t the data type of the matrix
// entries the MVMultDriverInterface class has to be a template class as well.
template<typename Number>
class MVMultDriverInterface {

public:

    typedef ::FullMatrixAccessor<Number> FullMatrixAccessor;

    MVMultDriverInterface()	{}

    virtual ~MVMultDriverInterface () {}

    // @sect4{Function: mvmult}
    //
    // This function has to implement the profiling of a
    // matrix-vector product in a derived class. There is
    // no reasonable default implementation thus we make this an
    // abstract virtual function.
    // @param y : The vector for the result of $Ax$.
    // @param A : The matrix to multiply with @p $x$.
    // @param x : The source vector.
    // @param n_repetitions : How many times the matrix-vector product is to be run.
    // @return : cumulative runtime.
    virtual double mvmult(std::vector<Number> & y,
                          const FullMatrixAccessor & A,
                          const std::vector<Number> & x,
                          int n_repetitions) = 0;
};

} // namespace step2 END
#endif // MVMultDriverInterface_H
