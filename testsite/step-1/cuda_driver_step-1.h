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


#ifndef CUDADriver_STEP_1_H
#define CUDADriver_STEP_1_H

#include <step-1/SimParams.h>

#include <lac/FullMatrixAccessor.h>
#include <lac/blas++.h>

#include <step-1/cuda_kernel_wrapper_step-1.cu.h>



        // Each example program is put into a separate namespace so that it
        // is easier to combine classes from different examples.
namespace step1 {

        // @sect3{Class: CUDADriver}
        //
        // This class is responsible for managing host-device communication,
        // mainly data transfer and invoking kernels.
        // To simplify host-device data transfer, we use our SciPal-library
        // which encapsulates all the details. The implementation of the kernels
        // is inherited privately in order to express that
        // this class is <i>implemented with</i> the Kernels class.
template<typename T>
class CUDADriver : private Kernels<T> {

public:
        // Some typedefs to make life easier.
    typedef typename blas_pp<T, cublas>::blas_wrapper_type BW;
    typedef typename blas_pp<T, cublas>::FullMatrixAccessor FullMatrixAccessor;
    typedef typename blas_pp<T, cublas>::Matrix Matrix;
    typedef typename blas_pp<T, cublas>::SubMatrix SubMatrix;
    typedef typename blas_pp<T, cublas>::MatrixSubCol MatrixSubCol;
    typedef typename blas_pp<T, cublas>::Vector Vector;
    typedef typename blas_pp<T, cublas>::SubColVector SubColVector;
    typedef typename blas_pp<T, cublas>::SubVectorBase SubVectorBase;

    typedef std::map<std::string,double> TimerName2Value;



    CUDADriver(const SimParams &p);



    double factorize(FullMatrixAccessor& A);

    double factorize(dealii::FullMatrix<T> &A);

    void chol_fac(FullMatrixAccessor& A, TimerName2Value& times);

    void lu_fac(FullMatrixAccessor& A, TimerName2Value& times);

    void single_thread_cholesky(FullMatrixAccessor& A);


private:
    Matrix A_d;

    const SimParams * params;
};

} // namespace step1 END

#endif // CUDADriver_STEP_1_H
