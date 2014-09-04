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

// On OSX10.9 the CUDATimer class causes linker errors. Therefore, we disable it.
//#define DONT_USE_CUDATIMER

#ifndef DONT_USE_CUDATIMER
#include <base/CUDATimer.h>
#endif

// With the advent of the Fermi architecture
// it became possible to use the printf command from within a kernel.
// To enable this feature one has to include the corresponding header from the C standard library.
#include <stdio.h>

// @sect3{Kernels}
//
// The core of any CUDA-based application are its kernels. They are the place where the
// actual numerical labour is executed.
// Like on the host side we put all kernels in a separate namespace
// to avoid collisions with other steps.
//
namespace step2 {

// To shorten the kernel code we introduce some abbreviations for built-in variables.
#define bx blockIdx.x
#define tx threadIdx.x
#define ty threadIdx.y


// The crucial idea of Fujimoto's method is to use
// the texture cache(s) to speed up the reading of the matrix by tiling and accessing it as a texture,
// as we have already mentioned in the introduction.
// To manage textures CUDA allows to either define texture objects or texture references
// (for the details have a look at the programming guide, Sec. 3.2.10).
// We follow Fujimoto and use references.
// A subtlety of texture references is that they must exist as global variables in the same file as the kernel
// which is going to access them.
// In computer graphics it is often desirable to normalize the values read from a texture to the
// range [-1., +1] or [0, +1.] depending on whether they are signed or unsigned.
// Defining the access as <i>cudaReadModeElementType</i> assures that
// during the texture fetch the data is read "as is".
// The immutable parameters needed at compile time are the element type, the spatial dimension and the type of access.
texture<float4, 2, cudaReadModeElementType> fTexRefA;

// @sect4{Kernel: _mv_fujimoto}
//
// This is the CUDA kernel for the matrix-vector multiplication
//
// \f{equation}
// y = Ax
// \f}
// introduced by Noriyuki Fujimoto.
//
// Copyright (C) 2008 Noriyuki Fujimoto, All Rights Reserved
// Fujimoto@mi.s.osakafu-u.ac.jp
//
// Please refer the paper below if you use Fujimoto's algorithm in your published work :
//
// Noriyuki Fujimoto, Faster Matrix-Vector Multiplication on GeForce 8800GTX,
// In the Proceedings of the 22nd IEEE International Parallel and
// Distributed Processing Symposium (IPDPS), LSPP-402, pp.1-8, April 2008
//
// This CUDA kernel heavily relies on bitshift operations.
// A good review can be found <a href="http://www.cprogramming.com/tutorial/bitwise_operators.html">here</a>.
//
// We start with the original version to explain the details of the parallelization.
// Afterwards we generalize it to double precision, consider the dependence of performance on the bitshifting and discuss
// an improved way of reading a double-precision matrix.
//
// To parallelize the multiplication of a large, dense matrix with a dense vector the matrix gets subdivided
// into small tiles.
// One row of tiles is assigned to one thread block of 16 x 16 threads. These sizes are chosen such that no thread of a
// block is idle when the entries of the vector are loaded into shared memory.
// Each thread block loads one tile at a time and computes
// its product with the corresponding slice from the source vectors.
// @param y : pointer to array containing the entries of the results vector.
// @param A : pointer to a cudaArray which the matrix has been copied.
// @param x : pointer to array containing the entries of the source vector.
// @param m : Number of rows of the matrix. The internal dimensions of the cudaArray may be different.
// @param n : Number of columns of the matrix. The internal dimensions of the cudaArray may be different.
__global__
void
_mv_fujimoto(float* y, cudaArray* A, float* x, int m, int n)
{
    // In the original version each 16x16 threadblock computes the matrix-vector product
    // of a 16x16 submatrix of the texture
    // with 256 entries of the source vector $x$ at a time.
    // Due to the frequent reuse of the vector entries they are stored in shared memory.
    __shared__ float xs[16][16];

    // The intermediate partial sums are buffered in shared memory
    __shared__ float Ps[16][16];

    // and the matrix entries read from the texture in a small vector.
    float4 a;

    // The 2D matrix indices must be converted to a 1D index
    // (leftshift of 4 = multiplication with 16 (2^4) which is the blocksize)
    float *Psptr = (float *) Ps + (ty << 4) + tx;
    // Then we determine the $y$ coordinate of the position in the texture, i.e. the global row index
    int ay = (bx << 4) + ty;
    // The entries of the solution vector
    // and the partial sums which will be computed by this thread
    // are addressed via pointers to hide the index arithmetic.
    float *xptr = x + (ty << 4) + tx;
    float *xsptr = (float *) xs + (tx << 2);

    // Each thread initializes its partial sum to 0.
    *Psptr = 0.0f;

    // Since the loop over the column index gets partially unrolled we have to define
    // it outside the loop unlike standard C++ procedure.
    int i;
    // Each row of a thread block deals with one row of the matrix.
    // Each thread walks through the texture with a stride of 64 and computes
    // the scalar product of 4 consecutive matrix elements
    // with 4 consecutive elements of the source vector.
    // The loop for the scalar product is unrolled.
    // With respect to the individual matrix elements this leads to an outer stride of 256.
    // Hence, one row of threads in the thread block simultaneously works on 64 texture elements
    // or 256 matrix and source vector elements respectively.
    // Thus, there are no idle threads for most of the time especially when the next bunch of source vector entries is loaded.
    // Exceptions may occur when the threads have reached the right boundary of the matrix.
    for (i = 0; i < (n & ~255); i += 256, xptr += 256)
    {
        // Since the source vector entries are the only ones which can be reused
        // they are copied into shared memory.
        xs[ty][tx] = *xptr;
        __syncthreads();

        // The core of the algorithm is to compute the scalar product
        // of 4 consecutive entries of a matrix row
        // with the corresponding subvector of the source vector.
        // The matrix is read through a texture to take advantage of their optimization with respect
        // to 2D spatial locality. Each texture element
        // represents 4 consecutive entries of a matrix row as float4, i.e. a vector of 4 floats.
        //
        // The texture elements are formed by 4-component vectors.
        // The 16 threads of one row of the thread block are supposed to read 16 consecutive texture elements.
        // Therefore, the $x$ coordinate for the texture fetch is given by the $x$ component of the
        // thread index @p tx plus the column index @p i divided by 4.
        // The division is realized by shifting the bits of @p i to the right.
        // This algorithm is also illustrated in
        // <a href="http://ch.nvidia.com/docs/IO/47905/fujimoto_lspp2008.pdf">Fujimoto's paper</a>.
        int ax = tx + (i >> 2);

        a = tex2D(fTexRefA, ax     , ay);
        // After reading the texture element we can compute the inner product of @p a with the corresponding
        // 4 entries of @p x. This is done by all threads and thus simultaneously for 64 consecutive elements
        // of $x$.
        *Psptr += a.x * *xsptr         + a.y * *(xsptr +   1) + a.z * *(xsptr +   2) + a.w * *(xsptr +   3);

        // The next statements repeat this procedure for the remaining 3 subgroups comprising
        // 64 elements of $x$ each.
        a = tex2D(fTexRefA, ax + 16, ay);
        *Psptr += a.x * *(xsptr +  64) + a.y * *(xsptr +  65) + a.z * *(xsptr +  66) + a.w * *(xsptr +  67);

        a = tex2D(fTexRefA, ax + 32, ay);
        *Psptr += a.x * *(xsptr + 128) + a.y * *(xsptr + 129) + a.z * *(xsptr + 130) + a.w * *(xsptr + 131);

        a = tex2D(fTexRefA, ax + 48, ay);
        *Psptr += a.x * *(xsptr + 192) + a.y * *(xsptr + 193) + a.z * *(xsptr + 194) + a.w * *(xsptr + 195);

        __syncthreads();
    }

    // At this point the number of remaining columns does not suffice to keep all threads busy with the unrolling strategy.
    // Especially, there are less than 256 vector entries which still can be read.
    if (i + (ty << 4) + tx < n) {
        xs[ty][tx] = *xptr;
    }
    __syncthreads();

    // After reading what is left over the 16 threads of a row in the thread block
    // work on the remaining columns in chunks of 64. Increasing @p xsptr only by 61
    // in each iteration is due to the fact that it is already incremented in the body of the
    // loop three times.
    int j;
    for (j = 0; j < ((n - i) >> 6); j++, xsptr += 61) {
        a = tex2D(fTexRefA, tx + (i >> 2) + (j << 4), ay);
        *Psptr += a.x * *xsptr++ + a.y * *xsptr++ + a.z * *xsptr++ + a.w * *xsptr;
    }
    __syncthreads();

    // When there are less than 64 columns left some threads have to stay idle while the rest
    // finishes the job.
    int remain = (n - i) & 63;

    if ((tx << 2) < remain) {
        a = tex2D(fTexRefA, tx + (i >> 2) + (j << 4), ay);
        *Psptr += a.x * *xsptr++;
    }
    if ((tx << 2) + 1 < remain) *Psptr += a.y * *xsptr++;
    if ((tx << 2) + 2 < remain) *Psptr += a.z * *xsptr++;
    if ((tx << 2) + 3 < remain) *Psptr += a.w * *xsptr;
    __syncthreads();

    // The last step, before we can write the result to global memory,
    // is to reduce the partial sums into one.
    if (tx < 8) *Psptr += *(Psptr + 8);
    if (tx < 4) *Psptr += *(Psptr + 4);
    if (tx < 2) *Psptr += *(Psptr + 2);
    if (tx < 1) *Psptr += *(Psptr + 1);
    __syncthreads();

    if (ty == 0 && (bx << 4) + tx < m) y[(bx << 4) + tx] = Ps[tx][0];
}



// @sect4{Kernel: _mv_fujimoto_T}
//
// A limitation of Fujimoto's original kernel is that it is implemented
// only for single-precision, real-valued matrices.
// Here, we show how to extend it to double precision. Complex numbers should be similar.
//
// The first thing we need is a structure which allows us to distinguish floats from doubles at compile-time.
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

// Then we need a template structure which depending on the number type @p T of the matrix elements provides
// the correct data types for the texture elements and the buffer of the matrix elements.
// Further, it has to provide the stride for the texture access which is needed in the determination
// of the effective width of the cudaArray through which the texture is accessed, see below in the wrapper function.
template<typename T>
struct TexEl;

// For single-precision, real numbers we need the types and values from the original implementation.
template<>
struct TexEl<float>
{
    typedef float4 value_type;
    typedef float4 alt_value_type;
    typedef float4 texel_type;
    static const int tex_stride = 4;
};

// For double precision we can use @p double4 as buffer but there are no such texture fetches. Therefore,
// we have to compose the texture of elements of type @p int4 and each double has to be reconstructed from 2 ints
// using CUDA's @p __hiloint2double function.
//
// For later use we introduce @p double2 as alternative value type.
template<>
struct TexEl<double>
{
    typedef double4 value_type;
    typedef double2 alt_value_type;
    typedef int4 texel_type;
    static const int tex_stride = 2;
};

// Last but not least we need a separate texture reference for double precision.
texture<TexEl<double>::texel_type, 2, cudaReadModeElementType> dTexRefA;


// The problem with the texture references is that we have to define one for each data type.
// Yet, we want to abstract the fetching from the data type. To do this, we define a little structure
// for the texture access which only provides an overloaded operator().
template<typename T>
struct texAccessor {

    __device__
    typename TexEl<T>::value_type
    operator() (int ax, int ay);

};

// Depending on the template specialization we can access the different texture references.
// For single precision we return what we fetch right away.
template<>
TexEl<float>::value_type texAccessor<float>::operator() (int ax, int ay)
{
    return tex2D(fTexRefA, ax, ay);
}

// In case of double precision this is also place where we reconstruct the matrix entries.
template<>
TexEl<double>::value_type texAccessor<double>::operator() (int ax, int ay)
{
    TexEl<double>::texel_type tmp;
    TexEl<double>::value_type a;

    tmp = tex2D(dTexRefA, 2*ax, ay);
    a.x = __hiloint2double(tmp.y, tmp.x);
    a.y = __hiloint2double(tmp.w, tmp.z);

    tmp = tex2D(dTexRefA, 2*ax+1, ay);
    a.z = __hiloint2double(tmp.y, tmp.x);
    a.w = __hiloint2double(tmp.w, tmp.z);

#ifndef nDEBUG
    printf("a : %f, %f, %f, %f\n", a.x, a.y, a.z, a.w);
#endif
    return a;
}

// Since there are only so few changes we only comment the differences.
// First of all, all occurences of @p float are replaced by @p T.
template<typename T>
__global__
void
_mv_fujimoto_T( T* y, cudaArray* A, T* x, int m, int n)
{
    __shared__ T xs[16][16];

    __shared__ T Ps[16][16];

    // Instead of a hard-coded type for buffer for the matrix elements
    // we let the TexEl structure determine the correct type from
    // which particular @T is passed in as template argument.
    // This is typical template metaprogramming.
    typename TexEl<T>::value_type a;

    T *Psptr = (T *) Ps + (ty << 4) + tx;

    int ay = (bx << 4) + ty;

    T *xptr = x + (ty << 4) + tx;
    T *xsptr = (T *) xs + (tx << 2);

    *Psptr = 0.0f;
    int i;

    // In contrast to the original kernel we now need
    // an accessor object for the texture fetches which
    // eliminiates the texture reference from the list of arguments of the texture fetch.
    texAccessor<T> tex_2D;

    // Except for the changes in the texture fetches everything else remains unchanged.
    for (i = 0; i < (n & ~255); i += 256, xptr += 256)
    {
        xs[ty][tx] = *xptr;
        __syncthreads();

        int ax = tx + (i >> 2);

        a = tex_2D(ax, ay);
        *Psptr += a.x * *xsptr         + a.y * *(xsptr +   1) + a.z * *(xsptr +   2) + a.w * *(xsptr +   3);

        a = tex_2D(ax+16, ay);
        *Psptr += a.x * *(xsptr +  64) + a.y * *(xsptr +  65) + a.z * *(xsptr +  66) + a.w * *(xsptr +  67);

        a = tex_2D(ax+32, ay);
        *Psptr += a.x * *(xsptr + 128) + a.y * *(xsptr + 129) + a.z * *(xsptr + 130) + a.w * *(xsptr + 131);

        a = tex_2D(ax+48, ay);
        *Psptr += a.x * *(xsptr + 192) + a.y * *(xsptr + 193) + a.z * *(xsptr + 194) + a.w * *(xsptr + 195);

        __syncthreads();
    }

    if (i + (ty << 4) + tx < n) {
        xs[ty][tx] = *xptr;
    }
    __syncthreads();

    int j;
    for (j = 0; j < ((n - i) >> 6); j++, xsptr += 61) {
        a = tex_2D(tx + (i >> 2) + (j << 4), ay);
        *Psptr += a.x * *xsptr++ + a.y * *xsptr++    +     a.z * *xsptr++ + a.w * *xsptr;
    }
    __syncthreads();

    int remain = (n - i) & 63;
    if ((tx << 2) < remain) {
        a = tex_2D(tx + (i >> 2) + (j << 4), ay);

        *Psptr += a.x * *xsptr++;
    }
    if ((tx << 2) + 1 < remain) *Psptr += a.y * *xsptr++;
    if ((tx << 2) + 2 < remain) *Psptr += a.z * *xsptr++;
    if ((tx << 2) + 3 < remain) *Psptr += a.w * *xsptr;
    __syncthreads();


    if (tx < 8) *Psptr += *(Psptr + 8);
    if (tx < 4) *Psptr += *(Psptr + 4);
    if (tx < 2) *Psptr += *(Psptr + 2);
    if (tx < 1) *Psptr += *(Psptr + 1);

    __syncthreads();
    if (ty == 0 && (bx << 4) + tx < m) y[(bx << 4) + tx] = Ps[tx][0];
}


// @sect4{Kernel: _mv_fujimoto_T2}
//
// Bitshifting maybe efficient. Yet it limits the readability of the code.
// Therefore, we provide one version where the bitshift operations are replaced by
// multiplications and divisions. Running the different Fujimoto versions then gives an overview
// of whether that technique pays off. Except for the bitshifts everything else is a verbatim copy
// of the previous kernel.
template<typename T>
__global__
void
_mv_fujimoto_T2( T* y, cudaArray* A, T* x, int m, int n)
{
    __shared__ T xs[16][16];

    __shared__ T Ps[16][16];

    typename TexEl<T>::value_type a;

    T *Psptr = (T *) Ps + (16*ty) + tx;

    int ay = (16*bx) + ty;

    T *xptr = x + (16*ty) + tx;
    T *xsptr = (T *) xs + (4*tx);

    *Psptr = 0.0f;
    int i;

    texAccessor<T> tex_2D;


    for (i = 0; i < (n & ~255); i += 256, xptr += 256)
    {
        xs[ty][tx] = *xptr;
        __syncthreads();

        int ax = tx + (i/4);

        a = tex_2D(ax, ay);
        *Psptr += a.x * *xsptr         + a.y * *(xsptr +   1) + a.z * *(xsptr +   2) + a.w * *(xsptr +   3);

        a = tex_2D(ax+16, ay);
        *Psptr += a.x * *(xsptr +  64) + a.y * *(xsptr +  65) + a.z * *(xsptr +  66) + a.w * *(xsptr +  67);

        a = tex_2D(ax+32, ay);
        *Psptr += a.x * *(xsptr + 128) + a.y * *(xsptr + 129) + a.z * *(xsptr + 130) + a.w * *(xsptr + 131);

        a = tex_2D(ax+48, ay);
        *Psptr += a.x * *(xsptr + 192) + a.y * *(xsptr + 193) + a.z * *(xsptr + 194) + a.w * *(xsptr + 195);

        __syncthreads();
    }

    if (i + (16*ty) + tx < n) {
        xs[ty][tx] = *xptr;
    }
    __syncthreads();

    int j;
    for (j = 0; j < ((n - i)/64); j++, xsptr += 61) {
        a = tex_2D(tx + (i/4) + (16*j), ay);
        *Psptr += a.x * *xsptr++ + a.y * *xsptr++    +     a.z * *xsptr++ + a.w * *xsptr;
    }
    __syncthreads();

    int remain = (n - i) & 63;
    if ((4*tx) < remain) {
        a = tex_2D(tx + (i/4) + (16*j), ay);

        *Psptr += a.x * *xsptr++;
    }
    if ((4*tx) + 1 < remain) *Psptr += a.y * *xsptr++;
    if ((4*tx) + 2 < remain) *Psptr += a.z * *xsptr++;
    if ((4*tx) + 3 < remain) *Psptr += a.w * *xsptr;
    __syncthreads();


    if (tx < 8) *Psptr += *(Psptr + 8);
    if (tx < 4) *Psptr += *(Psptr + 4);
    if (tx < 2) *Psptr += *(Psptr + 2);
    if (tx < 1) *Psptr += *(Psptr + 1);

    __syncthreads();
    if (ty == 0 && (16*bx) + tx < m) y[(16*bx) + tx] = Ps[tx][0];
}


// @sect4{Kernel: _mv_fujimoto_T3}
//
// The next issue is the optimization of double precision performance. The way it is handled up to now should
// lead to warps running at only half of their possible memory bandwidth.
// To change that we have to modify the texAccessor class and the magic numbers in the unrolled innermost loop.
template<typename T>
struct texAccessorOpt {

    __device__
    typename TexEl<T>::alt_value_type
    operator() (int ax, int ay);

};

// Depending on the template specialization we can access the different texture references.
// For single precision nothing changes.
template<>
TexEl<float>::alt_value_type texAccessorOpt<float>::operator() (int ax, int ay)
{
    return tex2D(fTexRefA, ax, ay);
}

// In case of double precision we now read only one texture element.
// This has to be compensated in the magic numbers for the partial unrolling of the loop over the columns of a row.
template<>
TexEl<double>::alt_value_type texAccessorOpt<double>::operator() (int ax, int ay)
{
    TexEl<double>::texel_type tmp;
    TexEl<double>::alt_value_type a;

    tmp = tex2D(dTexRefA, ax, ay);
    a.x = __hiloint2double(tmp.y, tmp.x);
    a.y = __hiloint2double(tmp.w, tmp.z);

#ifdef DEBUG
    printf("a : %f, %f\n", a.x, a.y);
#endif
    return a;
}

// In the kernel we distinguish between floats and doubles using our IsFloat structure.
template<typename T>
__global__
void
_mv_fujimoto_T3( T* y, cudaArray* A, T* x, int m, int n)
{
    __shared__ T xs[16][16];

    __shared__ T Ps[16][16];


    T *Psptr = (T *) Ps + (16*ty) + tx;

    int ay = (16*bx) + ty;

    T *xptr = x + (16*ty) + tx;


    *Psptr = 0.0f;
    int i;

    if(IsFloat<T>::value==true)
    {
        T *xsptr = (T *) xs + (4*tx);
        typename TexEl<T>::value_type a;
        texAccessor<T> tex_2D;
        for (i = 0; i < (n & ~255); i += 256, xptr += 256)
        {
            xs[ty][tx] = *xptr;
            __syncthreads();

            int ax = tx + (i/4);

            a = tex_2D(ax, ay);
             printf("a : %f, %f, %f, %f\n", a.x, a.y, a.z, a.w);
            *Psptr += a.x * *xsptr         + a.y * *(xsptr +   1) + a.z * *(xsptr +   2) + a.w * *(xsptr +   3);

            a = tex_2D(ax+16, ay);
            *Psptr += a.x * *(xsptr +  64) + a.y * *(xsptr +  65) + a.z * *(xsptr +  66) + a.w * *(xsptr +  67);

            a = tex_2D(ax+32, ay);
            *Psptr += a.x * *(xsptr + 128) + a.y * *(xsptr + 129) + a.z * *(xsptr + 130) + a.w * *(xsptr + 131);

            a = tex_2D(ax+48, ay);
            *Psptr += a.x * *(xsptr + 192) + a.y * *(xsptr + 193) + a.z * *(xsptr + 194) + a.w * *(xsptr + 195);

            __syncthreads();
        }

        if (i + (16*ty) + tx < n) {
            xs[ty][tx] = *xptr;
        }
        __syncthreads();


        int j;
        for (j = 0; j < ((n - i)/64); j++, xsptr += 61) {
            a = tex_2D(tx + (i/4) + (16*j), ay);
            *Psptr += a.x * *xsptr++ + a.y * *xsptr++ + a.z * *xsptr++ + a.w * *xsptr;
        }
        __syncthreads();


        int remain = (n - i) & 63;
        if ((4*tx) < remain) {
            a = tex_2D(tx + (i/4) + (16*j), ay);

            *Psptr += a.x * *xsptr++;
        }
        if ((4*tx) + 1 < remain) *Psptr += a.y * *xsptr++;
        if ((4*tx) + 2 < remain) *Psptr += a.z * *xsptr++;
        if ((4*tx) + 3 < remain) *Psptr += a.w * *xsptr;
        __syncthreads();
    }
    // For double precision we still have 16 multiplications in the innermost unrolled loop but we have to read
    // twice as often from the texture. Yet, it is done such that threads read contiguous pieces of memory when
    // the instruction for a texture fetch is issued unlike before.
    else {
        T *xsptr = (T *) xs + (2*tx);

        typename TexEl<T>::alt_value_type a;

        texAccessorOpt<T> tex_2D;

        for (i = 0; i < (n & ~255); i += 256, xptr += 256)
        {
            xs[ty][tx] = *xptr;

            __syncthreads();

            int ax = tx + (i/2);


            a = tex_2D(ax, ay);
            *Psptr += a.x * *xsptr         + a.y * *(xsptr +   1);

            a = tex_2D(ax+16, ay);
            *Psptr += a.x * *(xsptr +  32) + a.y * *(xsptr +  33);

            a = tex_2D(ax+32, ay);
            *Psptr += a.x * *(xsptr + 64) + a.y * *(xsptr + 65);

            a = tex_2D(ax+48, ay);
            *Psptr += a.x * *(xsptr + 96) + a.y * *(xsptr + 97);

            a = tex_2D(ax+64, ay);
            *Psptr += a.x * *(xsptr + 128) + a.y * *(xsptr + 129);

            a = tex_2D(ax+80, ay);
            *Psptr += a.x * *(xsptr + 160) + a.y * *(xsptr + 161);

            a = tex_2D(ax+96, ay);
            *Psptr += a.x * *(xsptr + 192) + a.y * *(xsptr + 193);

            a = tex_2D(ax+112, ay);
            *Psptr += a.x * *(xsptr + 224) + a.y * *(xsptr + 225);

            __syncthreads();
        }

        if (i + (16*ty) + tx < n) xs[ty][tx] = *xptr;

        __syncthreads();


        int j;

        for (j = 0; j < ((n - i)/32); j++, xsptr += 31)
        {
            a = tex_2D(tx + (i/2) + (16*j), ay);
            *Psptr += a.x * *xsptr++ + a.y * *xsptr;
        }
        __syncthreads();


        int remain = (n - i) & 31;

        if ((2*tx) < remain) {
            a = tex_2D(tx + (i/2) + (16*j), ay);

            *Psptr += a.x * *xsptr++;
        }

        if ( (2*tx) + 1 < remain)
        {
            *Psptr += a.y * *xsptr;
        }

        __syncthreads();
    }

    if (tx < 8) *Psptr += *(Psptr + 8);
    if (tx < 4) *Psptr += *(Psptr + 4);
    if (tx < 2) *Psptr += *(Psptr + 2);
    if (tx < 1) *Psptr += *(Psptr + 1);

    __syncthreads();
    if (ty == 0 && (16*bx) + tx < m) y[(16*bx) + tx] = Ps[tx][0];
}


// @sect4{Function: mv_fujimoto}
//
// After all the kernels we finally have to define the wrapper function for the
// kernel call.
// Its task is to setup the size of the thread blocks and the size of the grid
// which depends on the dimensions of the matrix.
// Furthermore, it binds the texture, allocates the memory and copies the
// matrix and vectors to the device.
// after all this is done it starts the selected kernel and when it has finished copies the result
// back to the host.
// To accomodate for the original kernel we have to provide an auxiliary wrapper function
// which in the case float calls Fujimoto's kernel and in all other cases does nothing.
void __do_FJ_orig(float* d_y, cudaArray* d_A, float* d_x, int m, int n, dim3 grid, dim3 threads)
{
    _mv_fujimoto<<< grid, threads >>>(d_y, d_A, d_x, m, n);
}

void __do_FJ_orig(double* , cudaArray* , double* , int , int , dim3 , dim3 )
{
    printf("For double precision the original version of Fujimoto is not available!\n");
}

// After these preparations we can define the generic wrapper function.
template<typename T>
void Kernels<T>::mv_fujimoto(T *y, const T *A, const T *x, const int m, const int n,
                             const int n_repetitions,
                             const int fj_version,
                             double&  elapsed_time)
{

    // To compute the number of blocks we have to divide the number of rows @p m
    // by 16 and add 1, if (m mod 16) != 0.
    // We keep the original version by Fujimoto as comment.
    int blkNum = (m + 15)/16; // (m >> 4) + ((m & 15) ? 1 : 0);
    int height = blkNum*16; // blkNum << 4;

    // For the width we have to do the same but with 16 replaced by 256.
    int width = 256*((n+255)/256); // (n & 255) ? (((n >> 8) + 1) << 8) : n;

    // Each row of a thread block deals with one row of the matrix.
    dim3 threads(16, 16);
    // Therefore, we need a one-dimensional grid of thread blocks which is large enough
    // so that all rows of the matrix are covered.
    dim3 grid(blkNum, 1);

    // The crucial idea of Fujimoto is to read the matrix via the texture cache.
    // To do this, the matrix must be stored in a cudaArray.
    cudaArray *d_A;
    T *d_x, *d_y;

    // In contrast to Fujimoto we do not
    // pass the type of the texture element directly to the channel description
    // but ask the TexEl structure for it. The original version is kept in a comment.
    cudaChannelFormatDesc
            channelDesc = cudaCreateChannelDesc<typename TexEl<T>::texel_type/*float4*/>();

    // In case of float we can store the matrix entries as float4.
    // Therefore, the width of the cudaArray is only
    // @p width/4.
    // For double precision we store 2 doubles as one int4 which requires
    // twice as much texels per row. Thus, the width of the cudaArray is
    // @p width/2.
    // To select these numbers automatically we use the static constant
    // @p tex_stride from the TeXEl structure.
    cudaMallocArray(&d_A, &channelDesc, width/TexEl<T>::tex_stride, height);

    size_t size_of_T = sizeof(T);

    cudaMemcpy2DToArray(d_A, 0, 0, A,
                        n * size_of_T,
                        n * size_of_T,
                        m,
                        cudaMemcpyHostToDevice);


    // The reference to the texture is created at runtime.
    // Depending on the number type we either bind the float4 or int4 texture.
    if (IsFloat<T>::value)
        cudaBindTextureToArray(fTexRefA, d_A);
    else
        cudaBindTextureToArray(dTexRefA, d_A);

    cudaMalloc((void **) &d_x, n * size_of_T );
    cudaMalloc((void **) &d_y, m * size_of_T );

    cudaMemcpy(d_x, x, n * size_of_T, cudaMemcpyHostToDevice);

    // Although the final values of matrix-vector product are assigned to  @p d_y
    // the result is only correct if we initialize
    // @p d_y with zeros. To do this, we copy for each run the
    // vector from the host to the device.
#ifndef DONT_USE_CUDATIMER
    CUDATimer timer;
#endif
    for (int i=0; i<n_repetitions; i++)
    {
        cudaMemcpy(d_y, y, m * size_of_T, cudaMemcpyHostToDevice);
        switch (fj_version) {
        case 0:
            __do_FJ_orig(d_y, d_A, d_x, m, n, grid, threads);
            break;
        case 1:
            _mv_fujimoto_T<<< grid, threads >>>(d_y, d_A, d_x, m, n);
            break;
        case 2:
            _mv_fujimoto_T2<<< grid, threads >>>(d_y, d_A, d_x, m, n);
            break;
        case 3:
            _mv_fujimoto_T3<<< grid, threads >>>(d_y, d_A, d_x, m, n);
            break;
        default:
            break;
        }
    }
    cudaThreadSynchronize();
#ifndef DONT_USE_CUDATIMER
    timer.stop();
    elapsed_time = timer.elapsed() ;
#else
    elapsed_time = 3.1415926;
#endif
#ifdef DEBUG
#ifndef DONT_USE_CUDATIMER
    timer.print_elapsed("Time spent in Fujimotos MV product:");
#endif
#endif

    // After the matrix-vector products are done we copy the result back to the host and clean up.
    cudaMemcpy(y, d_y, m * size_of_T, cudaMemcpyDeviceToHost);

    cudaFree(d_y);
    cudaFree(d_x);

    // Depending on which texture we have bound we have to selectively unbind.
    if(IsFloat<T>::value)
        cudaUnbindTexture( fTexRefA);
    else
        cudaUnbindTexture( dTexRefA);
    cudaFreeArray(d_A);
}

// This 2-liner provides all possible template specializations
// for real-valued matrices and finally explains why we enclosed the wrapper functions in a structure.
template class Kernels<float>;
template class Kernels<double>;

} // namespace step2 END


