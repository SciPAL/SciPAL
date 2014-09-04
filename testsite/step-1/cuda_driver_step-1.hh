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


#ifndef CUDA_DRIVER_STEP_1_HH
#define CUDA_DRIVER_STEP_1_HH

#include <base/CUDATimer.h>

#include <step-1/cuda_driver_step-1.h>

#include <step-1/cuda_kernel_wrapper_step-1.cu.h>



#include <QTime>


        // @sect4{Constructor: CUDADriver}
        //
        // In this case the constructor of the driver class does not do much.
template<typename T>
step1::CUDADriver<T>::CUDADriver(const SimParams &p)
    :
      params(&p)
{}


        // @sect4{Function: factorize}
        //
        // This function copies the data from the host to the device, starts the CUDA-based
        // factorization and copies the result back to the host.
        // There are two version of this function. The first one executes
        // the Cholesky decomposition as a whole.
// @param A : Matrix to factorize
// @return Number of seconds spent in GPU version of Cholesky factorization. NOT milliseconds.
template<typename T>
double
step1::CUDADriver<T>::factorize(FullMatrixAccessor& A)
{
    // Copy from host to device
    this->A_d = A;

    QTime t;
    t.start();
    // Call the multi-threaded factorization.
    // To do this, we have to dereference the Matrix object and retrieve
    // the bare pointer to the array containing the matrix entries.
    this->cholesky.blackbox(this->A_d.array().val(), this->A_d.n_cols(), this->A_d.leading_dim );

    // Convert milliseconds into seconds. Cf. QT-doc.
    double kernel_time = t.elapsed()/1000.;

    // Finally, copy the result from device back to host.
    A = this->A_d;

    return kernel_time;
}

// @sect4{Function: chol_fac}
      //
// The second version allows to
// measure the performance of the different parts of the algorithm
template<typename T>
void
step1::CUDADriver<T>::chol_fac(FullMatrixAccessor& A, TimerName2Value& times)
{
    // Copy from host to device
    this->A_d = A;

    // For timing measurements we use
    QTime t;

    // We dereference the Matrix object and retrieve
    // the bare pointer to the linear array containing the matrix entries.
    T* a_d = this->A_d.array().val();

    int n_rows = this->A_d.n_rows();

    int leading_dim = this->A_d.leading_dim;


    // Compute the number of blocks needed to cover the matrix.
    int n_blocks = (A.n_rows()+int(DEFAULT_TILE_SIZE)-1)/int(DEFAULT_TILE_SIZE);

    times["factorize_diagonal_block"] = 0.;
    times["strip_update"] = 0.;
    times["diag_update"] = 0.;
    times["lo_update"] = 0.;

    // Loop over the virtual diagonal blocks of the matrix.
    for(int i = n_blocks; i > 2; --i)
    {
        t.restart();


        cudaError_t error = this->cholesky.factorize_diag_block(a_d, n_blocks-i, n_rows, leading_dim);

        times["factorize_diagonal_block"]+=t.elapsed();

        AssertThrow(error == cudaSuccess, dealii::ExcMessage( cudaGetErrorString(error) ) );


        t.restart();
        this->cholesky.strip_update(a_d, n_blocks-i, i, n_rows, leading_dim);

        times["strip_update"]+=t.elapsed();


        t.restart();
        this->cholesky.diag_update(a_d, n_blocks-i, i, n_rows, leading_dim);

        times["diag_update"]+=t.elapsed();


        t.restart();
        this->cholesky.lo_update(a_d, n_blocks-i, n_blocks, i, n_rows, leading_dim);

        times["lo_update"]+=t.elapsed();
    }

    // For the last 2x2-Block submatrix @p lo_update() is not needed anymore.
    if(n_blocks>1)
    {
        this->cholesky.factorize_diag_block(a_d, n_blocks-2, n_rows, leading_dim);

        this->cholesky.strip_update(a_d, n_blocks-2, 2, n_rows, leading_dim);

        this->cholesky.diag_update(a_d, n_blocks-2, 2,  n_rows, leading_dim);
    }

    // Factorize the last diagonal block.

    std::cout << "Cholesky decomposition..." << std::endl;
    this->cholesky.factorize_diag_block(a_d, n_blocks-1, n_rows, leading_dim);

    // Finally, copy the result from device back to host.
    A = this->A_d;
}




        // @sect4{Function: single_thread_cholesky}
        //
        // This factorization method uses only one thread on the GPU. This demonstrates that CUDA can
        // be basically used like C.
        // @param A : Matrix to factorize
template<typename T>
void
step1::CUDADriver<T>::single_thread_cholesky(FullMatrixAccessor& A)
{
    // Copy from host to device
    this->A_d = A;

    // Call the single-threaded factorization.
    // To this end, we have to dereference the Matrix object and retrieve
    // the bare pointer to the array containing the matrix entries.
    this->cholesky.single_thread(this->A_d.array().val(),
                                 this->A_d.n_rows(), this->A_d.leading_dim );

    // Finally, copy the result from device back to host.
    A = this->A_d;
}

#endif // CUDA_DRIVER_STEP_1_HH
