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


/*
   Copyright 2010,2011,2012,2013 Stephan Kramer, as of 2013: Dr. Stephan Kramer

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   This file includes

   - the CUDA kernels for computing a Cholesky factorization
     of a real, symmetric (and hopefully positive definite) matrix.

   - device functions needed by the kernels for mapping thread and block indices
     to row and column indices of the matrix which is to be factorized or
     to the positions in the linear array holding the matrix entries.

   - the definitions of the wrapper functions for the kernels. These wrapper functions
     are declared in a separate header file and either allow to call
     the kernels individually or to execute the complete factorization
     via the 'blackbox' function.

   All kernels and wrapper functions are enclosed in a namespace 'step1'
*/
#ifndef CUDA_KERNEL_STEP_1_CU_H
#define CUDA_KERNEL_STEP_1_CU_H


// The most interesting part of this program is the way parallelization
// is implemented with CUDA. Therefore, we walk through the source code
// in a bottom-up fashion. It begins with the kernels and ends with the main function.
// A general feature of matrix operations is that large matrices get tiled
// into smaller submatrices on which the work is done at the end.
// Typically, the choice of the tile size depends on the hardware.

// @sect3{Parallelization of Cholesky Factorization}
//
// Looking at the algorithm given in the introduction, we see that the main computational
// effort is caused by computing the auxiliary variable @p sum in the update of the
// off-diagonal elements.
//
// For an efficient parallelization we subdivide the matrix into blocks and start with
// factorizing the upper leftmost diagonal block.
// Afterwards we compute the blocks in the columns below the diagonal block.
// Finally, we update the lower right part of the matrix which is still to be factorized
// by subtracting the auxiliary variable  @p sum.
// This happens in diag_update() and lo_update().

// For CUDA, since the size of a warp is 32, the most efficient choice for the blocksize
// is 16 in most cases where double precision is the type of choice for floating point numbers
// and 32 should be the optimum for single precision.
//#define DEFAULT_TILE_SIZE 16


namespace step1 {

// As of CUDA 5.0 it is possible to pass the size of shared memory arrays via
// template parameters into the kernels. Therefore we use a static constant
// rather than some macro to set the value for the size of matrix tiles.
// For CUDA 2.x and 3.x it did not seem to work.
// For 4.x we did not try.
//
// The square of this number gives the number of threads in a thread block.
// 16 has turned out to be a fairly good compromise.
static const int DEFAULT_TILE_SIZE = 16;

// @sect3{Class: Kernels}
//
// The CUDA-based kernels for the parallel computation are encapsulated
// into wrapper functions which are all collected into one structure.
// From a performance point of view one of the more interesting things is the
// dependence on the precision. Therefore, the kernels are templatized with respect to the
// number type @p T.
// This minimizes the amount of work for creating all the implementations
// for the different number types and precisions.
// To do this,
// we provide explicit template-specializations of this class.
// The number of specializations
// is kept at a minimum by grouping the wrapper functions
// for the different factorization methods
// into different private internal classes. Access to the kernels is
// via the @p lu and @p cholesky attributes.
template<typename T>
class Kernels {

    // @sect4{Class: Cholesky}
    //
    // This internal class provides the interface to
    // the actual kernels of the Cholesky factorization.
    struct Cholesky {

        void single_thread(T * A, int n_cols, int leading_dim);

        cudaError_t factorize_diag_block(T *A,
                                         int n_blocks_done, int n_cols, int leading_dim );

        void strip_update(T *A,
                          int n_blocks_done, int n_remaining_blocks, int n_cols, int leading_dim);

        void diag_update(T *A,
                         int n_blocks_done, int n_remaining_blocks, int n_cols, int leading_dim);

        void lo_update(T *A,
                       int n_blocks_done, int n_blocks, int n_remaining_blocks, int n_cols, int leading_dim);

        void blackbox(T *A, int n_cols, int leading_dim);

    };

    // @sect4{Class: LU}
    //
    // A possible extension of this program would be to add an
    // internal class which provides the interface to
    // LU-specific kernels. This is the reason why the wrapper
    // functions for the Cholesky kernels have been encapsulated
    // in a private, internal class.


public:
    Cholesky cholesky;

};

} // namespace step1 END

#endif // CUDA_KERNEL_STEP_1_CU_H
