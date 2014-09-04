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


#include <step-2/cuda_kernel_wrapper_step-2.cu.h>

#include <base/CUDATimer.h>

// @sect3{Kernels}
//
// Collect all kernels in a separate namespace to avoid collisions with other steps.
//
namespace step2 {

// abbreviations for built-in variables
#define bx blockIdx.x
#define tx threadIdx.x
#define ty threadIdx.y

#ifndef nUSE_TEX
texture<float4, 2, cudaReadModeElementType> texRefA;



// @sect4{Kernel: _mv_tex}
//
// This is the CUDA kernel for the matrix-vector multiplication
// invented by Noriyuki Fujimoto.
// y = Ax
// A : m-by-n matrix, x : n elements vector, y : m elements vector
// m and n are arbitrary positive integers
//
// Copyright (C) 2008 Noriyuki Fujimoto, All Rights Reserved
// Fujimoto@mi.s.osakafu-u.ac.jp
//
// Please refer the paper below if you use Fujimoto's algorithm in your published work :
//
// Noriyuki Fujimoto, Faster Matrix-Vector Multiplication on GeForce 8800GTX,
// In the Proceedings of the 22nd IEEE International Parallel and
// Distributed Processing Symposium (IPDPS), LSPP-402, pp.1-8, April 2008
__global__
void
_mv_tex( float* y, cudaArray* A, float* x, int m, int n)
{
    // 16x16 threadblock deals with a 16x16 submatrix.
    __shared__ float xs[16][16];
    // The intermediate partial sums are bufferd in shared memory.
    __shared__ float Ps[16][16];
    float4 a;
    // The 2D matrix indices must be converted to a 1D index
    // (leftshift of 4 = multiplication with 16 (2^4) which is the blocksize)
    float *Psptr = (float *) Ps + (ty << 4) + tx;
    // $y$ coordinate of position in the texture
    int ay = (bx << 4) + ty;
    // The entries of the solution vector
    // and the partial sums which will be computed by this thread
    // are addressed via pointers.
    float *xptr = x + (ty << 4) + tx;
    float *xsptr = (float *) xs + (tx << 2);

    // Each thread intializes its partial sum to 0.
    *Psptr = 0.0f;
    int i;
    // Each row of a thread block deals with one row of the matrix
    // each thread walks through the matrix with a stride of 64, but partial loop unrolling leads to an outer stride of 256
    for (i = 0; i < (n & ~255); i += 256, xptr += 256) {
        // copy source vector entries into shared memory
        xs[ty][tx] = *xptr;
        __syncthreads();

        // The core of the algorithm is to compute the scalar product
        // of 4 consecutive entries of a matrix row
        // with the corresponding subvector of the source vector.
        // The matrix is treaed like a texture where each entry
        // represents 4 consecutive entries of a matrix row as float4, i.e. a vector of 4 floats.
        //
        // $x$ coordinate of position in the texture
        int ax = tx + (i >> 2);

        a = tex2D(texRefA, ax     , ay);
        *Psptr += a.x * *xsptr         + a.y * *(xsptr +   1) + a.z * *(xsptr +   2) + a.w * *(xsptr +   3);

        a = tex2D(texRefA, ax + 16, ay);
        *Psptr += a.x * *(xsptr +  64) + a.y * *(xsptr +  65) + a.z * *(xsptr +  66) + a.w * *(xsptr +  67);

        a = tex2D(texRefA, ax + 32, ay);
        *Psptr += a.x * *(xsptr + 128) + a.y * *(xsptr + 129) + a.z * *(xsptr + 130) + a.w * *(xsptr + 131);

        a = tex2D(texRefA, ax + 48, ay);
        *Psptr += a.x * *(xsptr + 192) + a.y * *(xsptr + 193) + a.z * *(xsptr + 194) + a.w * *(xsptr + 195);

        __syncthreads();

    }

    // The following lines deal with the remaining (n mod 256) columns
    if (i + (ty << 4) + tx < n) {
        xs[ty][tx] = *xptr;
    }
    __syncthreads();
    int j;
    for (j = 0; j < ((n - i) >> 6); j++, xsptr += 61) {
        a = tex2D(texRefA, tx + (i >> 2) + (j << 4), ay);
        *Psptr += a.x * *xsptr++ + a.y * *xsptr++ + a.z * *xsptr++ + a.w * *xsptr;
    }
    __syncthreads();
    int remain = (n - i) & 63;
    if ((tx << 2) < remain) {
        a = tex2D(texRefA, tx + (i >> 2) + (j << 4), ay);
        *Psptr += a.x * *xsptr++;
    }
    if ((tx << 2) + 1 < remain) *Psptr += a.y * *xsptr++;
    if ((tx << 2) + 2 < remain) *Psptr += a.z * *xsptr++;
    if ((tx << 2) + 3 < remain) *Psptr += a.w * *xsptr;
    __syncthreads();

    // The last step before we can copy the result back to global memory
    // is to reduce the partial sums into one.
    if (tx < 8) *Psptr += *(Psptr + 8);
    if (tx < 4) *Psptr += *(Psptr + 4);
    if (tx < 2) *Psptr += *(Psptr + 2);
    if (tx < 1) *Psptr += *(Psptr + 1);

    __syncthreads();
    if (ty == 0 && (bx << 4) + tx < m) y[(bx << 4) + tx] = Ps[tx][0];
}
#endif

// @sect4{Function: mv_tex}
//
// Wrapper function for the matrix-vector product kernel.
// In contrast to the usual way, this function manages not only
// launching threads but also memory transfer.
//
#ifndef nUSE_TEX

template<typename T>
struct IsFloat;

template<>
struct IsFloat<float>
{
    static const bool value = true;
};

template<typename T>
struct IsFloat
{
    static const bool value = false;
};

template<typename T>
void Kernels<T>::mv_tex(T *y, const T *A, const T *x, int m, int n, int n_repetitions, double&  elapsedTime)
{

    if (!IsFloat<T>::value == true)
    {
        printf("Fujimoto double version not implemented\n");
        return;
    }

     // To compute the number of blocks we have to divide the number of rows @p m
    // by 16 and add 1, if (m mod 16) <>0
    int blkNum = (m >> 4) + ((m & 15) ? 1 : 0);
    int height = blkNum << 4;

    // For the width we have to do the same but with 16 replaced by 256.
    int width = (n & 255) ? (((n >> 8) + 1) << 8) : n;

    // Each row of a thread block deals with one row of the matrix
    dim3 threads(16, 16); dim3 grid(blkNum, 1); cudaArray *d_A; float *d_x, *d_y;

    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float4>();

    cudaMallocArray(&d_A, &channelDesc, width >> 2, height);

    cudaMemcpy2DToArray(d_A, 0, 0, A,
                        n * sizeof(float),
                        n * sizeof(float),
                        m,
                        cudaMemcpyHostToDevice);

    cudaBindTextureToArray(texRefA, d_A);

    cudaMalloc((void **) &d_x, n * sizeof(float));
    cudaMalloc((void **) &d_y, m * sizeof(float));

    cudaMemcpy(d_x, x, n * sizeof(float), cudaMemcpyHostToDevice);

    CUDATimer timer;
    for (int i=1; i<=n_repetitions;i++)
        _mv_tex<<< grid, threads >>>(d_y, d_A, d_x, m, n);

    cudaThreadSynchronize();
    timer.stop();
    elapsedTime = timer.elapsed() ;

    #ifdef DEBUG
            timer.print_elapsed("Time spent in Fujimotos MV product:");
    #endif

    cudaMemcpy(y, d_y, m * sizeof(float), cudaMemcpyDeviceToHost);

    cudaFree(d_y);
    cudaFree(d_x);
    cudaUnbindTexture(texRefA);
    cudaFreeArray(d_A);
}

#endif


// Bit stuff:
// & = The bits in the result at set to 1 if the corresponding bits in the
// two operands are both 1.

#define TB_dim_r 16

// The following must be a multiple of 4!
#define TB_dim_c 16



// @sect4{Kernel: __mv}
//
template<typename T>
__global__
void
__mv( T* y, T* A, const T* const x, int m, int n)
{
    // All threads work on the same set of vector entries.
    // So we store them in shared memory.
    T x_shm[TB_dim_r*TB_dim_c];
    // Storing only the intermediate partial sums in shared memory
    // Allows to reach the maximum number of resident blocks
    // on a multiprocessor.
    __shared__ T ps_shm[TB_dim_r*TB_dim_c];

    // The unavoidable loop over the columns of a row will
    // be partially unrolled into instruction groups of 4.
    float a[4];

    // The access to the correct place in the partial sum
    // array is managed by a pointer. Declaring it volatile
    // forces the compiler to perform all load and store operations
    // that we want. On the Fermi architecture they are sometimes
    // a victim of over-optimization.
    volatile T *ps = ps_shm + threadIdx.x + blockDim.x * threadIdx.y;

    volatile T * x_el = x_shm + threadIdx.x + blockDim.x * threadIdx.y;

    // Initialize the partial sum this thread
    // has to keep track of with zero.
    *ps = 0.;

    int block_size = TB_dim_r*TB_dim_c;

#ifdef CUDA_DEBUG
    printf("tx : %d, ty : %d \n", threadIdx.x, threadIdx.y);
#endif

    int c_stride = 0;
    int n_el_per_thread = min(4, n/TB_dim_c);
    // Threads are arranged in a 2D array where @p threadIdx.x is running fastest.
    // Each thread processes a subset of the elements of one row. All threads with the
    // same @p threadIdx.x work on the same row.
    // Each one of them will read four matrix elements.
    // All threads with the same @p threadIdx.y belong to the same warp
    // (depending on the dimensions of the thread block).
    // Thus, at the beginning all threads of a warp load the same source vector element from
    // shared memory. On the Fermi architecture this is supported by the broadcast operation.
    //
    //
    c_stride = n_el_per_thread * TB_dim_c;

    // The outer loop loads entries from the source vector.
    for (int cb = 0; cb < n; cb += block_size) {

        // Use all threads to load elements of
        // the source vector into shared memory.
        int global_index = threadIdx.x + blockDim.x * threadIdx.y + cb;

        // Let the thread read whatever it finds in memory.
        *x_el = x[global_index];
        __syncthreads();

#ifdef CUDA_DEBUG
        if (global_index < n)
        {
            printf("x_%d : %f\n", global_index, *x_el );
        }
#endif

        int a_begin = 0;
        // Offset due to rows processed by previous thread blocks.
        int row = TB_dim_r*blockIdx.x;
        // Offset due to rows processed by previous threads.
        row +=  threadIdx.x;

        a_begin += n * row;

        a_begin += threadIdx.y *n_el_per_thread;

        // The inner loop computes the local partial sums.
        int c_end = min(cb + block_size, n);
        for (int c_begin = cb; c_begin < c_end; c_begin += c_stride)
        {
#ifdef CUDA_DEBUG
            printf("row : %d, bid : %d, tx : %d, ty : %d, c_begin : %d, a_begin : %d\n",
                   row, blockIdx.x, threadIdx.x, threadIdx.y,  c_begin, a_begin);
#endif

            for (int k = 0; k < n_el_per_thread; k++)
            {
                a[k] = A[a_begin+k];
#ifdef CUDA_DEBUG
                printf("A(%d, %d) : %f - %d = %f\n",
                       row, threadIdx.y *n_el_per_thread + c_begin,
                       a[k], n*row +  threadIdx.y *n_el_per_thread + c_begin +k + 1,
                       a[k] - (n*row +  threadIdx.y *n_el_per_thread + c_begin +k + 1));
#endif
            }

            T result = 0.;

            for (int k = 0; k < n_el_per_thread; k++)
            {
                result += a[k] * x[threadIdx.y *n_el_per_thread + c_begin +k];//x_el[n_el_per_thread*threadIdx.y + k + c_begin -cb];
#ifdef CUDA_DEBUG
                printf("result : %f, a_%d : %f, x_%d : %f\n", result , k, a[k],
                       n_el_per_thread*threadIdx.y + k + c_begin -cb,
                       x[threadIdx.y *n_el_per_thread + c_begin +k]
                       // x_el[n_el_per_thread*threadIdx.y + k + c_begin -cb]
                       );
#endif
            }
            *ps += result;
            __syncthreads();
        }
    }
    __syncthreads();

#ifdef CUDA_DEBUG
    printf("PS[%d][%d] = %f\n", threadIdx.x, threadIdx.y, *ps);
#endif

    //    if (threadIdx.y < 8) *ps += *(ps + 8);
    //     __syncthreads();
    //    if (threadIdx.y < 4) *ps += *(ps + 4);
    //     __syncthreads();
    //    if (threadIdx.y < 2) *ps += *(ps + 2);
    //     __syncthreads();
    //    if (threadIdx.y < 1) *ps += *(ps + 1);
    //
    //    __syncthreads();
    if (ty == 0 /*&& (bx << 4) + tx < m*/)
    {
        T y_el = *ps;
        for (int k = 1; k < blockDim.y; k++)
            y_el += *(ps + blockDim.x * k);
        y[blockIdx.x * blockDim.x + threadIdx.x] = y_el; //ps_shm[threadIdx.x];
    }
#ifdef ksdfsdf
    // Loop over columns
    for (int i = 0; i < (n & ~255); i += 256, xptr += 256) {
        // The following relies on having an efficient
        // load operation of broadcast type.
        x_el = x[threadIdx.y];


        ax = A[m * ax + ay];
        ay = A[m * ax + ay + 1];
        az = A[m * ax + ay + 2];
        aw = A[m * ax + ay + 3];
        // a = tex2D(texRefA, ax     , ay);
        *Psptr += ax * *xsptr         + ay * *(xsptr +   1) + az * *(xsptr +   2) + aw * *(xsptr +   3);
        // a = tex2D(texRefA, ax + 16, ay);
        *Psptr += ax * *(xsptr +  64) + ay * *(xsptr +  65) + az * *(xsptr +  66) + aw * *(xsptr +  67);
        // a = tex2D(texRefA, ax + 32, ay);
        *Psptr += ax * *(xsptr + 128) + ay * *(xsptr + 129) + az * *(xsptr + 130) + aw * *(xsptr + 131);
        // a = tex2D(texRefA, ax + 48, ay);
        *Psptr += ax * *(xsptr + 192) + ay * *(xsptr + 193) + az * *(xsptr + 194) + aw * *(xsptr + 195);
        __syncthreads();
    }

    if (i + (ty << 4) + tx < n) {
        xs[ty][tx] = *xptr;
    }
    __syncthreads();
    int j;
    for (j = 0; j < ((n - i) >> 6); j++, xsptr += 61) {
        // a = tex2D(texRefA, tx + (i >> 2) + (j << 4), ay);
        *Psptr += ax * *xsptr++ + ay * *xsptr++ + az * *xsptr++ + aw * *xsptr;
    }
    __syncthreads();
    int remain = (n - i) & 63;
    if ((tx << 2) < remain) {
        // a = tex2D(texRefA, tx + (i >> 2) + (j << 4), ay);
        *Psptr += ax * *xsptr++;
    }
    if ((tx << 2) + 1 < remain) *Psptr += ay * *xsptr++;
    if ((tx << 2) + 2 < remain) *Psptr += az * *xsptr++;
    if ((tx << 2) + 3 < remain) *Psptr += aw * *xsptr;
    __syncthreads();

    if (tx < 8) *Psptr += *(Psptr + 8);
    if (tx < 4) *Psptr += *(Psptr + 4);
    if (tx < 2) *Psptr += *(Psptr + 2);
    if (tx < 1) *Psptr += *(Psptr + 1);

    __syncthreads();
    if (ty == 0 && (bx << 4) + tx < m) y[(bx << 4) + tx] = Ps[tx][0];
#endif


}




// @sect4{Function: _mv}
//
template<typename T>
void
_mv(T * y, const T * const A, const T * x, const int m, const int n, double& time1)
{
    const int dim_r = TB_dim_r;
    const int dim_c = TB_dim_c;

    dim3 threads(dim_r, dim_c);

    int blkNum = (m + dim_r -1)/dim_r;   //(m >> 4) + ((m & 15) ? 1 : 0);


    dim3 grid(blkNum, 1);


    T *d_A, *d_x, *d_y;

    cudaMalloc((void **) &d_A, n * m * sizeof(T));
    //cudaMalloc((void **) &d_x, std::max(n,m) * sizeof(T)); //n
    //cudaMalloc((void **) &d_y, std::max(n,m) * sizeof(T)); //m

    cudaMalloc((void **) &d_x, n * sizeof(T)); //n
    cudaMalloc((void **) &d_y, m * sizeof(T)); //m

    cudaMemcpy(d_x, x, n * sizeof(T),     cudaMemcpyHostToDevice);
    cudaMemcpy(d_A, A, n * m * sizeof(T), cudaMemcpyHostToDevice);
    cudaMemcpy(d_y, y, m * sizeof(T),     cudaMemcpyHostToDevice);

    cudaThreadSynchronize();
    CUDATimer timer;

    //__mv<<< grid, threads >>>mv_tex(d_y, d_A, d_x, m, n);

    __mv<<< grid, threads >>>(d_y, d_A, d_x,  m,n);

    //  __mv<<< grid, threads >>>(d_x, d_A, d_y, n, m);

    cudaThreadSynchronize();
    timer.stop();
    time1 = timer.elapsed();

    // timer.print_elapsed("Time spent in Fujimotos MV product without textures:");
    // timer.print_elapsed(" ");

    cudaMemcpy(y, d_y, m * sizeof(T), cudaMemcpyDeviceToHost);

    cudaFree(d_y);
    cudaFree(d_x);
    cudaFree(d_A);
}

// @sect4{Function: mv}
template<typename T>
void
Kernels<T>::mv(T * y, const T * const A, const T * const x,
               const int m, const int n, double& time1)
{

    if (!IsFloat<T>::value == true)
    {
        printf("Fujimoto double version not implemented\n");
        return;
    }

    _mv(y, A, x, m, n, time1);

}

// This 2-liner provides all possible template specializations
// for real-valued matrices.
template class Kernels<float>;
template class Kernels<double>;

} // namespace step2 END


