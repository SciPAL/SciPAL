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


#ifndef CUDA_DRIVER_STEP_2_HH
#define CUDA_DRIVER_STEP_2_HH

#include <base/CUDATimer.h>

// @sect4{Constructor: CUBlasDriver}
//
// In case of CUBLAS the constructor and destructor
// take care of initializing the library and shutting it down again.
template<typename T,typename blasType>
step2::CUBlasDriver<T,blasType>::CUBlasDriver() : MVMultDriverInterface<T>()
{
    blasType::Init();
}


// @sect4{Destructor: CUBlasDriver}
template<typename T,typename blasType>
step2::CUBlasDriver<T,blasType>::~CUBlasDriver()
{
    blasType::Shutdown();
}


// @sect4{Function: mvmult}
//
// For details about the function arguments see the base class.
template<typename T,typename blasType>
double step2::CUBlasDriver<T,blasType>::mvmult(std::vector<T>& y,
                                               const FullMatrixAccessor& A,
                                               const std::vector<T>& x,
                                               int n_repetitions)
{   
    // As before, we introduce some convenience variables
    // and local variables for the arguments of the BLAS function.
    int n_rows = A.n_rows();
    int n_cols = A.n_cols();

    T * dst = &y[0];

    T * A_entries = const_cast<FullMatrixAccessor &>(A).val();

    T * src = &const_cast<std::vector<T> &>(x)[0];

    T alpha=1.;
    T beta=1.;
    int incx=1;
    int incy=1;

    // Since CUBLAS is also designed to work column-major everything we have said in the body of
    // CpuBlasDriver::mvmult() applies here as well.
    int lda=n_cols;
    int m = n_cols;
    int n = n_rows;

    // Up to this point there is no difference to CpuBlasDriver::mvmult().
    // In order to use CUBLAS we have to copy the matrix, the source and the destination vector
    // to the device.
    // To do this, we have to declare pointers in the device memory ...
    T *d_A, *d_x, *d_y;

    // allocate enough space ...
    cudaMalloc((void **) &d_A, n_rows * n_cols * sizeof(T));

    cudaMalloc((void **) &d_x, n_cols * sizeof(T));
    cudaMalloc((void **) &d_y, n_rows * sizeof(T));

    // and copy everything from host to device.
    cudaMemcpy(d_x, src, n_cols * sizeof(T), cudaMemcpyHostToDevice);

    cudaMemcpy(d_A, A_entries, n_rows * n_cols * sizeof(T), cudaMemcpyHostToDevice);

    // <p> It should be obvious that excessive use of @p cudaMalloc, @p cudaMemcpy and @p cudaFree quickly
    // becomes error-prone and leads to a lot of redundant code.
    // To remedy this, we have introduced our SciPal library which hides memory transfers behind assignment
    // operators and the bare blas functions
    // are wrapped up in an operator-based interface for linear algebra operations.</p>
    //
    // After these preparations we execute the CUBLAS version of gemv() via SciPal's CUBLAS wrapper
    // which allows us to use a precision-independent formulation. After each run we have to synchronize
    // in order to measure the actual runtime and not the time CUBLAS needs to start its kernels.
    double cumulative_runtime = 0.;
    for (int i=0; i<n_repetitions; i++)
    {
        cudaThreadSynchronize();
        CUDATimer timer;

        // Like BLAS, CUBLAS does not take care of properly initializing
        // the result vector to zero. Hence, for each run we copy the
        // vector from the host to the device.
        cudaMemcpy(d_y, dst, n_rows * sizeof(T), cudaMemcpyHostToDevice);

        blasType::gemv ('t', m, n, alpha,
                        d_A, lda,
                        d_x,  incx,
                        beta,
                        d_y,  incy);

        cudaThreadSynchronize();
        timer.stop();
        cumulative_runtime += timer.elapsed();
    }

    // Once we are done, we can copy the result back to the host
    // and return the cumulative runtime.
    cudaMemcpy(dst, d_y, n_rows * sizeof(T), cudaMemcpyDeviceToHost);

    // After the test is complete, we deallocate the arrays on the device.
    cudaFree(d_y);
    cudaFree(d_x);
    cudaFree(d_A);

    return cumulative_runtime;
}
#endif // CUDA_DRIVER_STEP_2_HH
