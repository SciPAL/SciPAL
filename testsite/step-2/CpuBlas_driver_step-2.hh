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


#ifndef CPUBLAS_DRIVER_STEP_2_HH
#define CPUBLAS_DRIVER_STEP_2_HH

#include <base/CUDATimer.h>

// @sect4{Function: mvmult}
//
// We overload the mvmult() function from the base class
// to delegate the execution of the matrix-vector product
// to the ATLAS library.
template<typename EntryType,typename blasType>
double step2::CpuBlasDriver<EntryType,blasType>::mvmult(std::vector<EntryType>& y,
                                                        const FullMatrixAccessor& A,
                                                        const std::vector<EntryType>& x,
                                                        int n_repetitions)
{
    // We introduce some convenience variables.
    int n_rows = A.n_rows();
    int n_cols = A.n_cols();

    // For better readability we
    // create local variables for
    // the arguments of the BLAS function.
    EntryType * d_y = &y[0];

    EntryType * d_A = const_cast<FullMatrixAccessor &>(A).val();

    EntryType * d_x = &const_cast<std::vector<EntryType> &>(x)[0];

    EntryType alpha = 1;
    EntryType beta  = 1;

    int incx = 1; int incy = 1;

    // Internally, BLAS assumes a column-major
    // storage of the matrix. Therefore, the leading
    // dimension is the number of columns.
    int lda = n_cols;

    // However, we use a row-major ordering and thus have to invoke
    // BLAS' gemv() function such that it performs the transposed matrix-vector product.
    // To do this, we have to swap the roles of the number of rows and columns
    // and pass a 't' as first argument to gemv().
    int m = n_cols;
    int n = n_rows;

    // We execute the matrix-vector product several times
    // for a given matrix size to get a reliable estimate of the runtime.
    double cumulative_runtime = 0.;
    for (int i=0; i<n_repetitions; i++)
    {
        // BLAS does not take care of properly initializing
        // the result vector to zero. Hence, we have to do it ourselves.
        for (unsigned int k = 0; k < n_rows; k++)
            d_y[k] = 0.;

        CUDATimer timer;

        blasType::gemv ('t', m, n, alpha,
                        d_A, lda,
                        d_x, incx,
                        beta,
                        d_y, incy);

        timer.stop();
        cumulative_runtime += timer.elapsed();
    }

    return cumulative_runtime;
}
#endif // CPUBLAS_DRIVER_STEP_2_HH
