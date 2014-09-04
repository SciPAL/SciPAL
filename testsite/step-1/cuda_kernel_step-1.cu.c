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
#ifdef USE_CU_C

// In order to use printf() from within a kernel
// the good old C-header has to be included.
#include <stdio.h>

#include <step-1/cuda_kernel_wrapper_step-1.cu.h>

namespace step1 {

// @sect3{Device Functions}
//
// Before discussing the kernels we take a look at so-called device functions.
// They execute only on the GPU and up to CUDA 3.2 are automatically inline.
// Nowadays they may be not inline. To enforce inlining the keyword @p __forceinline__
// get introduced.
//
// @sect4{Device Function: lex_index_2D}
//
// Compute a position $r\cdot row\_length + c$ in a linear array
// from a row index $r$ and column index $c$
// of a matrix entry. It is assumed that the matrix is stored row-wise.
// @param r : row index
// @param c : column index
// @param leading_dim : number of matrix entries per row
__forceinline__
__device__ int lex_index_2D(int r, int c, int leading_dim)
{
    return c +  r*leading_dim;
}


// @sect4{Device Function: global_pos}
//
// Compute a global row or column index given the Block index and thread index.
// @param t_pos : local row or comlumn index within block
// @param n_blocks_done : Index of the block in C-counting
template<int TILE_SIZE>
__forceinline__
__device__ int global_pos(int t_pos, int n_blocks_done)
{
    return t_pos + TILE_SIZE*n_blocks_done;
}


// @sect4{Device Function: inv_sqrt}
//
// CUDA provides single and double-precision versions for the
// computation of the reciprocal of a square root.
// In order to use it in a template context we have to provide
// our own little wrapper function for unifying the name
// by employing C++'s polymorphism. As for the other device functions we would like
// to enforce inlining. However, this leads to compiler errors
// (at least when using the nvcc which ships with CUDA 5).
// Therefore, we do not inline.
// The function @p rsqrtf
// is documented in the CUDA Programming Guide.
// @param x : real to take the square root of.
__device__ float inv_sqrt(float x)
{
    return rsqrtf(x);
}


__device__ double inv_sqrt(double x)
{
    return rsqrt(x);
}

// @sect3{Cholesky Kernels}
//
// Kernels run on the GPU and are global functions. They cannot be members
// of a class but of a namespace. Therefore, in order to group them according
// to the factorization method they are used for we put them into different naemspaces.
//
namespace Chol {

// @sect4{Kernel: __single_thread}
//
// This kernel is only started once to
// compute the factorization of the whole matrix.
// This shows that the source of the CPU-based Cholesky factorization
// would also run unchanged on the GPU.
// In memory the matrix is stored as a linear array whose length is a
// multiple of @p leading_dim. This is for improved memory troughput. For instance,
// this allows to make rows to have a length which is a multiple of the cache line length
// which avoids misaligned accesses.
// Since for Cholesky factorizations a matrix must be square it does not matter
// whether it is stored in row- or column-major format.
// @param A : linear array containing the matrix entries in row-major order.
// @param n_rows : length of a column
// @param leading_dim :
template<typename T>
__global__
void
__single_thread(T *A, const int n_rows, const int leading_dim)
{
    // The outer loop runs over the rows <i>L</i>.
    for (unsigned int r = 0; r < n_rows; ++r)
    {
        // Compute diagonal entry of Cholesky factor.
        T sum = 0.;
        unsigned int idx;
        unsigned int idx_c;
        // Sum squares of entries computed so far in this row.
        for (unsigned int u = 0; u < r; ++u)
        {
            idx = lex_index_2D(r, u, leading_dim);
            sum += A[idx] * A[idx];
        }
        idx = lex_index_2D(r, r, leading_dim);
        A[idx] = sqrt(A[idx] - sum);

        // Off-diagonal entries. Here, we exploit the symmetry of <i>A</i>.
        // The auxiliary variable @p sum corresponds to step 2.1 of the algorithm
        // given in the introduction.
        for (unsigned int c = r+1; c < n_rows; ++c)
        {
            sum = 0.;

            for (unsigned int u = 0; u < r; ++u)
            {
                idx_c = lex_index_2D(c, u, leading_dim);
                idx   = lex_index_2D(r, u, leading_dim);
                sum += A[idx_c]*A[idx];
            }

            idx_c = lex_index_2D(c, r, leading_dim);
            idx   = lex_index_2D(r, c, leading_dim);
            A[idx_c]  = A[idx] - sum;

            idx   = lex_index_2D(r, r, leading_dim);
            A[idx_c] /= A[idx];
        }
    }
}


//
// @sect4{Kernel: factorize_diag_block}
//
// This kernel factorizes a diagonal block assuming that all previous
// diagonal blocks have already been factored.

// In contrast to a serial implementation we hide the summation
// of the off-diagonal elements from the factorized part
// in the usage of the  <i>thread</i> index and in the choice of
// synchronization points.
//
// Each instance of this kernel computes one matrix entry of the Cholesky factor $L$.
// @param A : Linear array containing the elements of matrix to factorize
// @param n_blocks_done : distance of the block which is to be factorized
// from the left uppermost diagonal block.
// @param n_cols : length of a row of $A$.
template<typename T, int TILE_SIZE>
__global__
void
__factorize_diag_block(T *A, int n_blocks_done,
                       int n_cols, int leading_dim)
{
    // In C, arrays are stored row-wise; thus the $x$ coordinate
    // of @p threadIdx indicates the column index.
    int col = threadIdx.x;

    // The $y$ coordinate
    // of @p threadIdx indicates the row.
    int row = threadIdx.y;

    // From the thread and block index we have to compute the index
    // of the matrix element this thread has to work on.
    // This is delegated to device functions.
    int global_row = global_pos<TILE_SIZE>(row, n_blocks_done);
    int global_col = global_pos<TILE_SIZE>(col, n_blocks_done);

    // For matrices whose number of rows is not
    // a multiple of @p TILE_SIZE we have to
    // take care that thread do not work on non-existing matrix entries.
    if ((global_row >= n_cols) || (global_col >= n_cols))
        return;

    int idx = lex_index_2D(global_row, global_col, leading_dim);

    // Simplify debugging especially of index problems we provide
    // possibility for some output.
#ifdef INDEX_BOUND_DEBUG
    if (row == 0 && col == 0)
    printf("%s:\n------------------------\n", __FUNCTION__);
    __syncthreads();

    printf("row, col : (%d, %d), g_row,g_col : (%d, %d), idx : %d\n ",
           row, col, global_row, global_col, idx);
#endif

    // Copy the diagonal block to <i>shared memory</i>
    // so that threads can exchange their results.
    // To avoid memory bank conflicts
    // when consecutively accessing the elements of a column we add one column.
    // This trick is discussed in the <i>CUDA Best Practices Guide</i> in the
    // chapter about multiplying a matrix with its transpose.
    __shared__ T L[TILE_SIZE][TILE_SIZE+1];

    // To minimize the number of accesses to global memory we copy
    // the matrix entries of the block to shared memory and synchronize
    // all threads within the block before we go on.
    L[row][col]= A[idx];
    __syncthreads();

    T fac;

    // Now, we can compute the entries of the Cholesky factors.
    // We have to distinguish between diagonal elements
    // i.e. $ r = c = k$,
    // elements of the uppermost row, and the rest.
    // For matrices whose number of rows is not
    // a multiple of @p TILE_SIZE we have to modify
    // the upper bound of the loop over the diagonal of
    // matrix block.
    int k_max = TILE_SIZE;
    // Next, figure out whether we are in the rightmost block column of the matrix.
    // Note that @p global_pos(0, bo) cannot exceed @p n_cols due to the way threads get started.
    if (n_cols - global_pos<TILE_SIZE>(0, n_blocks_done) < TILE_SIZE)
        k_max = n_cols%TILE_SIZE;

#ifdef INDEX_BOUND_DEBUG
      printf("k_max : %d\n", k_max);
#endif

    for(int k=0; k < k_max; k++)
    {
        __syncthreads();
        // Compute the inverse square root of diagonal element $1/\sqrt{A_{kk}}$
        // using the device function defined above,
        // and store the result in a register.
        // The device function is necessary to encapsulate the precision-dependent
        // function names of the reciprocal square root.
        // Precomputing $1/\sqrt{A_{kk}}$ allows to map the division by $\sqrt{A_{kk}}$
        // of the off-diagonal elements of $L$ to a multiplication which is computationally cheaper.
        fac = inv_sqrt(L[k][k]);
        __syncthreads();

        // We compute the diagonal element and the row to its right-hand side
        // \f{eqnarray}
        // L_{ck}  =  fac \cdot A_{ck} & = &
        // \left\lbrace
        // \begin{array}{ll}
        // \sqrt{A_{kk}} & c = k \\
        // \frac{A_{ck}}{\sqrt{A_{kk}}} & c > k
        // \end{array}
        // \right.
        // \f}
        if ((row==k)&&(col>=k)) L[col][row]=(L[col][row])*fac;

        // Again we have to synchronize.
        __syncthreads();

        // Next, compute the lower left triangle.

        // The off-diagonal entries follow from
        // \f{eqnarray}
        // s_k & = & L_{ck}L_{rk} \\
        // s_c & = & \sum_{k=c+1}^{TILE\_SIZE}(-s_k) \\
        // L_{rc} & = & \frac{1}{\sqrt{A_{kk}}} \left(A_{rc}
        // + s_c \right)\,.
        // \f}
        // We only have to perform the first of these steps, i.e. computing $s_k$,
        // and one addition in the sum.
        // Hence, $s_c$ is computed incrementally by the @p k loop.
        // Instead of storing  $s_c$, we subtract the individual terms  $s_k$ from $A_{rc}$
        // and store the modified  $A_{rc}$ in the Cholesky factor.
        //
        // When the condition
        // @p row == @p k is true, @p L[row][col] contains the final value
        // $A_{rc} + s_c$. The multiplication with @p fac, that is $1/\sqrt{A_{kk}}$,
        // is performed implicitly during this process.
        if ((row>=col)&&(col>k)) L[row][col] -= L[col][k]*L[row][k];
    }

    __syncthreads();

    // At the end, we copy the result back to global memory.
    if (row>=col) A[idx] = L[row][col];


#ifdef INDEX_BOUND_DEBUG
    printf("A_%d : %f\n", idx,  L[row][col]);
#endif
}


// @sect4{Kernel: strip_update}
//
// This kernel updates $L_{ij}$ in the columns below
// diagonal block <i>j</i>.
//
// @param A : global Matrix
// @param n_blocks_done : distance from the diagonal block just factorized.
// @param n_cols : length of a row of $A$.
template<typename T, int TILE_SIZE>
__global__
void
__strip_update(T *A, int n_blocks_done, int n_cols, int leading_dim)
{
    // @p boffy und @p boffx are the coordinates of the matrix block
    // this thread works on.
    int boffy=n_blocks_done;

    // The "+1" is needed since @p n_blocks_done denotes the position of the left uppermost
    // block.
    int boffx = blockIdx.x + boffy + 1;

    // As in the last kernel.
    int col = threadIdx.x;
    int row = threadIdx.y;

    // Again avoid bank conficts by adding a column
    __shared__ T topleft[TILE_SIZE][TILE_SIZE+1];
    __shared__ T workingmat[TILE_SIZE][TILE_SIZE+1];

    // Grab the data of the diagonal block just factorized ...
    int global_row = global_pos<TILE_SIZE>(row,n_blocks_done);
    int global_col = global_pos<TILE_SIZE>(col,n_blocks_done);

    // For matrices whose number of rows is not
    // a multiple of @p TILE_SIZE we have to
    // take care that thread do not work on non-existing matrix entries.
    if ((global_row >= n_cols) || (global_col >= n_cols))
        return;

    int idx = lex_index_2D(global_row, global_col, leading_dim);

    topleft[row][col]=A[idx];

    // and the transposed block which is to be processed
    global_row = global_pos<TILE_SIZE>(row,boffx);
    int idx_w = lex_index_2D(global_row, global_col, leading_dim);

    workingmat[col][row] = A[idx_w];

    __syncthreads();


    int k_max = TILE_SIZE;

    // Do step 2.2 of the algorithm given in the introduction.
    // Each thread works on one column.
    if(row==0)
        for (int k=0; k < k_max; k++)
        {
            T sum=0.;
            for (int m = 0; m < k; m++)
                sum += topleft[k][m]*workingmat[m][col];

            workingmat[k][col] = (workingmat[k][col] - sum)/topleft[k][k];
        }

    __syncthreads();

    A[idx_w] = workingmat[col][row];
}

// @sect4{Kernel: diag_update}
//
// This Kernel computes the contribution of the last factorized diagonal block
// to the auxiliary variable @p sum in step 2.1.
// @param A : global matrix
// @param n_blocks_done : block offset from upper left corner of the matrix.
// @param n_cols : length of a row of $A$.
template<typename T, int TILE_SIZE>
__global__
void
__diag_update(T *A, int n_blocks_done, int n_cols, int leading_dim)
{
    // Finding the global indices and setup of the shared memory is as above.
    int boffx = blockIdx.x + n_blocks_done + 1;

    int col = threadIdx.x;
    int row = threadIdx.y;

    int global_row = global_pos<TILE_SIZE>(row, boffx);
    int global_col = global_pos<TILE_SIZE>(col, n_blocks_done);

    // For matrices whose number of rows is not
    // a multiple of @p TILE_SIZE we have to
    // take care that threads do not work on non-existing matrix entries.
    if ((global_row >= n_cols) || (global_col >= n_cols))
        return;

    int idx = lex_index_2D(global_row, global_col, leading_dim);

    __shared__ T left[TILE_SIZE][TILE_SIZE+1];

    // Copy and synchronize.
    left[row][col]= A[idx];

    __syncthreads();

    // The thread with index (row, col) computes the corresponding
    // term from step 2.1.
    T sum = 0.f;


    int k_max = TILE_SIZE;

    if(row>=col)
    {
        for(int kk=0; kk<k_max; kk++) sum += left[row][kk]*left[col][kk];

        // Subtract the result from the global Matrix entry.
        global_col = global_pos<TILE_SIZE>(col, boffx);
        idx = lex_index_2D(global_row, global_col, leading_dim);

        A[idx] -= sum;
    }
}


// @sect4{Kernel: lo_update}
//
// This kernel applies the intermediate results produced by strip_update() and diag_update()
// to the rest of the matrix.
// @param A : global matrix
// @param n_blocks_done : Block-offset
// @param n_blocks : Number of blocks, in which a row of @p A is subdivided.
// @param n_cols : length of a row of @p A.
template<typename T, int TILE_SIZE>
__global__
void
__lo_update(T *A, int n_blocks_done, int n_blocks, int n_cols, int leading_dim)
{

    // Start with local and global Indices
    // of the entry this thread works on ...
    int col = threadIdx.x;
    int row = threadIdx.y;

    int boffy=blockIdx.y+n_blocks_done+1;
    int boffx=boffy+1;

    __shared__ T left[TILE_SIZE][TILE_SIZE];

    // The extra column is needed only for @p upt.
    __shared__ T upt[TILE_SIZE][TILE_SIZE+1];


    // Start reading data at the lower left.
    int global_row_src = global_pos<TILE_SIZE>(row, boffy);
    int global_col_src = global_pos<TILE_SIZE>(col, n_blocks_done);

    // For matrices whose number of rows is not
    // a multiple of @p TILE_SIZE we have to
    // take care that thread do not work on non-existing matrix entries.
    if ((global_row_src >= n_cols) || (global_col_src >= n_cols))
        return;

    int idx = lex_index_2D(global_row_src, global_col_src, leading_dim);

    upt[row][col]=A[idx];
    __syncthreads();


#ifdef nSCHUR_DEBUG
    if (row == 0 && col == 0)
    printf("%s block (%d,%d):\n------------------------\n",
           __FUNCTION__, blockIdx.x, blockIdx.y);
    __syncthreads();

#endif

    for (;boffx<n_blocks;boffx++)
    {
        int global_row = global_pos<TILE_SIZE>(row, boffx);
        idx = lex_index_2D(global_row, global_col_src, leading_dim);

        // Reset shared memory.
        left[row][col]= 0.;
        left[row][col]=A[idx];

#ifdef SCHUR_DEBUG
        printf("loading  left[%d][%d]=A[%d] == %f\n", row, col, idx, A[idx]);
#endif
        __syncthreads();

        if (global_row < n_cols)
        {
            T matrixprod=0.f;

            // The thread with index (row, col) computes the corresponding term from step 2.2.
            int k_max = TILE_SIZE;

            for (int kk=0;kk<k_max;kk++)
            {
                matrixprod+=left[row][kk]*upt[col][kk];

#ifdef SCHUR_DEBUG
                printf("row, col : (%d, %d), g_row,g_col : (%d, %d), idx : %d, mprod : %f, L_r%d : %f, U_1c : %f \n ",
                       row, col, global_row, global_col_src, idx, matrixprod, kk, left[row][kk], upt[col][kk]);
#endif
            }

            int global_col = global_pos<TILE_SIZE>(col, boffy);

            if (global_col < n_cols)
            {
                idx = lex_index_2D(global_row, global_col, leading_dim);
                A[idx] -= matrixprod;

#ifdef SCHUR_DEBUG
                if (row == 0 && col == 0)
                    printf("%s:\n------------------------\n", __FUNCTION__);
                __syncthreads();

                printf("row, col : (%d, %d), g_row,g_col : (%d, %d), idx : %d, mprod : %f\n ",
                       row, col, global_row, global_col, idx, matrixprod);
#endif
            }
        }
    }

}

} // namespace Chol END

} // namespace step1 END


// @sect3{Wrapper functions}
//
// The wrapper functions mainly repeat the arguments of the kernels. Besides that
// they manage the set up of the thread blocks and grids. Due to the tiling of
// matrix all kernels will be executed by grids of 2-dimensional thread blocks.
// The block size reflects the tiling. The grids are always 1-dimensional
// since the off-diagonal update operations are effectively independent
// updates of several rows or columns at once.
//
// @sect4{Function: blackbox}
//
// This function provides a complete GPU-based implementation
// of the Cholesky factorization and is an example of how to
// call several kernels from one wrapper function.
template<typename T>
void step1::Kernels<T>::Cholesky::blackbox(T * a_d, int n_cols, int leading_dim)
{
    cudaError_t error;

    // Compute the number of blocks needed to cover the matrix.
    int n_blocks = (n_cols+int(DEFAULT_TILE_SIZE)-1)/int(DEFAULT_TILE_SIZE);

    // A thread-block should be as large as a matrix block.
    dim3 threads(DEFAULT_TILE_SIZE, DEFAULT_TILE_SIZE);

    dim3 logrid;

    for(int i=n_blocks; i>2; --i)
    {
        logrid.x=1;
        logrid.y=i-2;

        dim3 stripgrid(i-1);

        // For the diagonal block we need only one block, thus the grid size is 1.
        Chol::__factorize_diag_block<T, DEFAULT_TILE_SIZE>
                <<<1, threads>>>(a_d, n_blocks-i, n_cols, leading_dim);
        cudaThreadSynchronize();

        Chol::__strip_update<T, DEFAULT_TILE_SIZE>
                <<<stripgrid, threads>>>(a_d, n_blocks-i, n_cols, leading_dim);
        cudaThreadSynchronize();

        Chol::__diag_update<T, DEFAULT_TILE_SIZE>
                <<<stripgrid, threads>>>(a_d, n_blocks-i, n_cols, leading_dim);
        cudaThreadSynchronize();
        Chol::__lo_update<T, DEFAULT_TILE_SIZE>
                <<< logrid, threads >>>(a_d, n_blocks-i, n_blocks, n_cols, leading_dim);
        cudaThreadSynchronize();
    }

    // For the last 2x2-Block submatrix @p lo_update() is not needed anymore.
    if(n_blocks>1)
    {
        Chol::__factorize_diag_block<T, DEFAULT_TILE_SIZE><<<1, threads>>>(a_d, n_blocks-2, n_cols,
                                                   leading_dim);
        cudaThreadSynchronize();

        Chol::__strip_update<T, DEFAULT_TILE_SIZE><<<1, threads>>>(a_d, n_blocks-2, n_cols, leading_dim);
        cudaThreadSynchronize();

        Chol::__diag_update<T, DEFAULT_TILE_SIZE><<<1, threads>>>(a_d, n_blocks-2, n_cols, leading_dim);
        cudaThreadSynchronize();

    }

    // Factorize the last diagonal block.
    Chol::__factorize_diag_block<T, DEFAULT_TILE_SIZE><<<1, threads>>>(a_d, n_blocks-1, n_cols, leading_dim);

    // Make all threads finish before the host is allowed to go on.
    // Remember, kernel starts are asynchronous.
    cudaThreadSynchronize();

    // ... check error state of GPU.
    error=cudaGetLastError();
    if (error != cudaSuccess)
    {
        printf("     Error code %d: %s.\n",error,cudaGetErrorString(error));
        exit(-1);
    }
}


// @sect4{Function: single_thread}
//
// This function wraps the kernel call.
//
// @param a_d : Pointer to device memory containing the matrix entries.
// @param n_rows : the number of rows of the matrix. This implies the number of columns and total number of elements as the matrix behind @p a_d is supposed to be square.
template<typename T>
void step1::Kernels<T>::Cholesky::single_thread(T * a_d, int n_cols, int leading_dim)
{
    // Start the kernel which is supposed to work only in one thread.
    Chol::__single_thread<<<1,1>>>(a_d, n_cols, leading_dim);
}



// @sect4{Function: factorize_diag_block}
//
template<typename T>
cudaError_t step1::Kernels<T>::Cholesky::factorize_diag_block(T * a_d,
                                                              int n_blocks_done,
                                                              int n_cols,
                                                              int leading_dim)
{
    dim3 threads(DEFAULT_TILE_SIZE, DEFAULT_TILE_SIZE);
    Chol::__factorize_diag_block<T, DEFAULT_TILE_SIZE><<<1, threads>>>(a_d, n_blocks_done,  n_cols, leading_dim);
    cudaThreadSynchronize();

    return cudaGetLastError();
}




// @sect4{Function: strip_update}
//
template<typename T>
void
step1::Kernels<T>::Cholesky::strip_update(T *a_d,
                                          int n_blocks_done,
                                          int n_remaining_blocks, int n_cols,
                                          int leading_dim)
{
    cudaError_t error;
    dim3 stripgrid(n_remaining_blocks-1);
    dim3 threads(DEFAULT_TILE_SIZE, DEFAULT_TILE_SIZE);

    Chol::__strip_update<T, DEFAULT_TILE_SIZE><<<stripgrid, threads>>>(a_d,
                                                 n_blocks_done, n_cols,
                                                 leading_dim);

    // Every update must be synchronized, because we can't continue the factorization before
    // having updated the submatrix of A.
    cudaThreadSynchronize();

    // The last task is to query the error state of the CUDA context.
    error=cudaGetLastError();
    if (error != cudaSuccess)
    {
        printf("     Error code %d: %s.\n",error,cudaGetErrorString(error));

        // In case of an error
        exit(-1);
    }
}





// @sect4{Function: diag_update}
//
template<typename T>
void step1::Kernels<T>::Cholesky::diag_update(T *a_d,
                                              int n_blocks_done,
                                              int n_remaining_blocks,
                                              int n_cols,
                                              int leading_dim)
{
    cudaError_t error;
    dim3 stripgrid(n_remaining_blocks-1);
    dim3 threads(DEFAULT_TILE_SIZE, DEFAULT_TILE_SIZE);

    Chol::__diag_update<T, DEFAULT_TILE_SIZE><<<stripgrid, threads>>>(a_d,
                                                n_blocks_done, n_cols,
                                                leading_dim);

    cudaThreadSynchronize();
    error=cudaGetLastError();
    if (error != cudaSuccess)
    {
        printf("     Error code %d: %s.\n",error,cudaGetErrorString(error));
        exit(-1);
    }
}





// @sect4{Function: lo_update}
//
template<typename T>
void step1::Kernels<T>::Cholesky::lo_update(T *a_d,
                                            int n_blocks_done,
                                            int n_blocks,
                                            int n_remaining_blocks ,
                                            int n_cols,
                                            int leading_dim)
{
    cudaError_t error;
    dim3 logrid;
    logrid.x=1;
    logrid.y=n_remaining_blocks-2;
    dim3 threads(DEFAULT_TILE_SIZE, DEFAULT_TILE_SIZE);

    Chol::__lo_update<T, DEFAULT_TILE_SIZE><<< logrid, threads >>>(a_d,
                                             n_blocks_done, n_blocks,  n_cols, leading_dim);
    cudaThreadSynchronize();
    error=cudaGetLastError();
    if (error != cudaSuccess)
    {
        printf("     Error code %d: %s.\n",error,cudaGetErrorString(error));
        exit(-1);
    }
}



// This 2-liner provides all possible template sepcializations
// for real-valued matrices.
template class step1::Kernels<float>;
template class step1::Kernels<double>;

#endif
