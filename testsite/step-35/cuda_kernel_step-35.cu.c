//@sect3{File: cuda_kernel_step-35.cu.c}
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

Copyright Stephan Kramer, Johannes Hagemann,  Lutz Künneke, Jan Lebert 2014-2016
*/

//Header containing the declarations of the wrapper functions - not the kernels.
//This header is the interface between the code compiled by nvcc and the one compiled by gcc.

//std
#include <stdio.h>
#include <omp.h>

//CUDA
#include <cuComplex.h>

//Thrust
#include <thrust/functional.h>

#ifndef nUSE_CPP11
#include <thrust/transform.h>
#include <thrust/sort.h>
#include <thrust/execution_policy.h>
#include <thrust/device_ptr.h>
#include <thrust/fill.h>
#endif


//SciPAL
#include <base/PrecisionTraits.h>

//Our stuff
#include "cuda_helper.h"

#include <step-35/cuda_kernel_wrapper_step-35.cu.h>

//@sect4{Constant CUDA memory}
//Maximum subset size is 512 (reasonable choices are more like 15-30)
#ifdef DOUBLE_PRECISION
__constant__ double cs[512];
#else
__constant__ float cs[512];
#endif

//@sect5{Function: set_cs}
//@brief Copy $c_s$ to constant memory
//@param cs_h pointer to $c_s$ on host
//@param maxnum number of elements in cs_h
template<typename T, ParallelArch arch>
void step35::Kernels<T, arch>::set_cs(T* cs_h, int maxnum) {
    checkCudaErrors(cudaMemcpyToSymbol(cs, cs_h, sizeof(T)*maxnum));
}

//@sect4{Device Functions}
//
//Prior to the kernels we have to define the device functions that we need.

//@sect5{Struct: SharedMemory}
//Workaround as CUDA doesn't immediately support shared memory arrays in templated functions,
//see <a href="http://docs.nvidia.com/cuda/cuda-samples/#simple-templates">CUDA Simple Templates example</a>.
template <typename T>
struct SharedMemory {
    // Ensure that we won't compile any un-specialized types
    __device__ T *getPointer() {
        extern __device__ void error(void);
        error();
        return NULL;
    }
};

//Float specialization
template <>
struct SharedMemory <float> {
    __device__ float *getPointer() {
        extern __shared__ float s_float[];
        return s_float;
    }
};

//Double specialization
template <>
struct SharedMemory <double> {
    __device__ double *getPointer() {
        extern __shared__ double s_double[];
        return s_double;
    }
};

//@sect5{Struct: __abs}
//@brief Unary functions to calculate the absolute value for an element
template <typename T>
struct __abs : public thrust::unary_function<T,T> {
    // Ensure that we won't compile any un-specialized types
    __device__
    T operator()(T x) {
        extern __device__ void error(void);
        error();
        return NULL;
    }
};

//Float specialization: use fabsf()
template <>
struct __abs <float> : public thrust::unary_function<float,float> {
    __device__
    float operator()(float x) {
        return fabsf(x);
    }
};

//Double specialization: use fabs()
template <>
struct __abs <double> : public thrust::unary_function<double,double> {
    __device__
    double operator()(double x) {
        return fabs(x);
    }
};

//@sect5{Struct: __element_norm_product}
//@brief Unary functions to calculate $\frac{x\cdot y}{n}$ where $x,y$ are complex numbers and $n$ a constant factor.
template <typename T>
struct __element_norm_product {
    __element_norm_product(T _n) {}
    // Ensure that we won't compile any un-specialized types
    __device__
    T operator()(T x) {
        extern __device__ void error(void);
        error();
        return NULL;
    }
};

//Float specialization
template <>
struct __element_norm_product <float> {
    cuFloatComplex normalization;

    __element_norm_product(float _n) {
        normalization = make_cuFloatComplex(1.0/_n, 0);
    }

    __host__ __device__
    cuFloatComplex operator()(const cuFloatComplex& x, const cuFloatComplex& y) const {
        return cuCmulf(cuCmulf(x, y), normalization);
    }
};

//Double specialization
template <>
struct __element_norm_product <double> {
    cuDoubleComplex normalization;

    __element_norm_product(double _n) {
        normalization = make_cuDoubleComplex(1.0/_n, 0);
    }

    __host__ __device__
    cuDoubleComplex operator()(const cuDoubleComplex& x, const cuDoubleComplex& y) const {
        return cuCmul(cuCmul(x, y), normalization);
    }
};


//@sect5{Function: power_of_two}
//@param x exponent
//@return $2^x$ for $x \ge 0$. For $x < 0$ returns 0.
__device__ int pow_of_two(int x) {
    if (x < 0) {
        return 0;
    }
    return 1 << x;
}

//@sect5{Function: mysign}
//Sign function
template<typename T>
__device__ T mysign(T in)
{
    if ( in < (T)0)
        return (T)(-1.0);
    else
        return (T)1;
}


// @sect5{Function: InverseErf}
//
// For the computation of the weights $c_S$ we need the inverse error function. The implementation is taken from the
// <a href="https://people.maths.ox.ac.uk/gilesm/files/gems_erfinv.pdf">paper</a> by Mike Giles.
__forceinline__ __device__ float MBG_erfinv(float x)
{
    float w, p;
    w = - __logf((1.0f-x)*(1.0f+x));

    if ( w < 5.000000f )
    {
        w = w - 2.500000f;
        p =  2.81022636e-08f;
        p =  3.43273939e-07f + p*w;
        p = -3.5233877e-06f  + p*w;
        p = -4.39150654e-06f + p*w;
        p =  0.00021858087f  + p*w;
        p = -0.00125372503f  + p*w;
        p=  -0.00417768164f  + p*w;
        p=   0.246640727f    + p*w;
        p=   1.50140941f     + p*w;
    }
    else
    {
        w = sqrtf(w) - 3.000000f;
        p = -0.000200214257f;
        p =  0.000100950558f + p*w;
        p =  0.00134934322f  + p*w;
        p = -0.00367342844f  + p*w;
        p =  0.00573950773f  + p*w;
        p = -0.0076224613f   + p*w;
        p =  0.00943887047f  + p*w;
        p =  1.00167406f     + p*w;
        p =  2.83297682f     + p*w;
    }
    return p*x;
}


// @sect4{Kernel}
//
// CUDA extends C by <i>Kernel</i>s, which are executed in parallel by
// <i>threads</i>. Furthermore, one can define <i>device functions</i>
// which can only be called from inside a kernel.
// <br>
// Note that unlike device functions kernels cannot be members of class.

//@sect5{Kernel: __dykstra}
//@brief This Kernel provides the Gaussian Dykstra projection.
//@param qg pointer to the q vectors in global memory
//@param A path to the image
//@param cinfo array of pointers where to find n,i,j and numpatch
//@param q_offset array of offsets where to find first $q$ value for a threadblock
//@param ni height of the image
//@param nj width of the image
//NOTICE: only 2D functionality
template<typename T>
__global__ void
__dykstra(T *qg, T *A, int **cinfo, const int * q_offset, const int ni, const int nj) {
    //Find out where I live

    int k=0;
    int border=0;
    int *n=cinfo[4*blockIdx.x];
    int *i=cinfo[4*blockIdx.x+1];
    int *j=cinfo[4*blockIdx.x+2];
    int numpatch=cinfo[4*blockIdx.x+3][0];

    while ( k < numpatch ) {
        if ( border+n[k]*n[k] > threadIdx.x )
            break;
        border+=n[k]*n[k];
        k++;
    }
    if ( threadIdx.x - border >= n[k]*n[k] || k >= numpatch) {
        //This is bad, exit
        return;
    }
    //Allocate shared memory
    SharedMemory<T> shared;
    T* s_mem = shared.getPointer();
    T* q=&(s_mem[0]); //takes first n2 sizeof(T) bytes
    T* f=&(s_mem[blockDim.x]); //takes second third
    T* e=&(s_mem[2*blockDim.x]); //takes last third

    int is = (threadIdx.x-border)/n[k] + i[k];
    int js = (threadIdx.x-border)%n[k] + j[k];

    //Pixel index in the current in the global image
    const int idx = is*nj + js;

    //Fill shared memory with variables we later use
    q[threadIdx.x] = qg[threadIdx.x+q_offset[blockIdx.x]];

    f[threadIdx.x] = A[idx] - q[threadIdx.x];
    e[threadIdx.x] = f[threadIdx.x]*f[threadIdx.x];
    __syncthreads();

    //m is the next power of 2 bigger than n divided by two
    int m=n[k]*n[k];
    int mmax=pow_of_two(floor(log2f((float)m)));
    //Sum over all pixels in one frame
    m=1;
    while (m <= mmax) {
        if (threadIdx.x - border + m < n[k]*n[k]) {
            e[threadIdx.x] += e[threadIdx.x + m];
        }
        m = m << 1;
        __syncthreads();
    }

    //$q = x_{r+1} - x_r$
    q[threadIdx.x] = -f[threadIdx.x];
    if (((T) cs[n[k]-1])*e[border] > 1.0) {
        f[threadIdx.x] = f[threadIdx.x]/(sqrt(((T) cs[n[k]-1])*e[border]));
    }

    //Update global memory (image and q)
    A[idx] = f[threadIdx.x];
    qg[threadIdx.x+q_offset[blockIdx.x]] = q[threadIdx.x] + f[threadIdx.x];
}



template<typename T>
__forceinline__ __device__ void L2_norm_on_sub_rows(int s, T * ps_sum)
{
    int row = threadIdx.y;
    int col = threadIdx.x;
    int tIx = row*blockDim.x + col;

    // Pointer to first element of a row which (i.e. the element) is processed by this thread.
    volatile T * m_ps = ps_sum + tIx;

    T x =  // 1.0; //
            *m_ps;
    __syncthreads();

    *m_ps = x * x;
    __syncthreads();

    int edge_x = pow_of_two(s);


    // ------------------
    //    for (int e = edge_x; e>0; e/= 2)
    //    {
    //        int lane_id = col % (edge_x);

    //        if (lane_id < e/2)
    //        {
    //            *m_ps += *(m_ps + e/2);
    //            if (blockIdx.x == 0)
    //                printf ("col : %d, e : %d, lane_id : %d \n", col, e, lane_id);
    //        }

    //__syncthreads();
    //    }
    // ------------------



    // Each line is processed by one warp
    for (int e = 2; e <= edge_x; e*= 2)
    {
        if ((col+e) % e == 0)
        {
            *m_ps += *(m_ps + e/2);
            // if (*m_ps != e)
            // printf("PS_%d(%d, %d) : %f, e/2 : %d, tIx : %d\n", e, row, col, *m_ps, e/2, tIx);
        }
        __syncthreads();
    }

}


template<typename T>
__forceinline__ __device__ void sum_of_squares_on_2D_subsets(int s, T * ps_sum)
{
    int row = threadIdx.y;
    int col = threadIdx.x;
    int tIx = row*blockDim.x + col;

    // Pointer to first element of a row which (i.e. the element) is processed by this thread.
    volatile T * m_ps = ps_sum + tIx;

    // For debugging the 2D summation set to 1. Then each scale results in a power of 4 as the partial sum over the pixels of the subset.
    T x = // 1.0; //
            *m_ps;
    __syncthreads();

    *m_ps = x * x;
    __syncthreads();

    int edge_x = pow_of_two(s);

    // Each line is processed by one warp
    for (int e = 2; e <= edge_x; e*= 2)
    {
        if ((col+e) % e == 0)
        {
            *m_ps += *(m_ps + e/2);
            //  printf("PS_%d(%d, %d) : %f, e/2 : %d, tIx : %d\n", e, row, col, *m_ps, e/2, tIx);
        }
        // __syncthreads();
    }

    __syncthreads();

    int edge_y = pow_of_two(s);

    for (int e = 2; e <= edge_y; e*= 2) // add rows
    {
        __syncthreads(); // Needed for some reason to ensure thread safety.
        int neighbor = blockDim.x*(row + e/2) + col;
        if ((row+e) % e == 0)
        {
            //            if (*m_ps != ps_sum[neighbor] /*&& row == 0 && col == 0*/)
            //            {
            //                printf("adding rows : ps_sum[%d] : %f, ps_sum[%d] : %f, e : %d\n", tIx, *m_ps, neighbor, ps_sum[neighbor], e);
            //            }
            *m_ps += ps_sum[neighbor];
        }

        __syncthreads();

    }

}


template<typename T>
__forceinline__ __device__ T project(int s_x, int s_y, T * ps_sum, T x, T data, T g_noise, T weight)
{
    int row = threadIdx.y;
    int col = threadIdx.x;
    int tIx = row*blockDim.x + col;

    // FIXME: move into kernel
    int edge_x = pow_of_two(s_x);
    int edge_y = pow_of_two(s_y);

    int ws_isx = blockDim.x*edge_y*(row/edge_y) +  edge_x*(col/edge_x);

    T ws = (ps_sum[ws_isx]/(g_noise*g_noise)) * weight;
    __syncthreads();

    if (ws > 1.)
        x *= 1.
                / sqrt(ws);

    // For debugging the 2D summation
    if (false) //  ( blockIdx.x == blockIdx.y) && (ws != T(edge_x*edge_y)) ) //  && blockIdx.y == gridDim.y/2)
        printf("ws_isx : %d, ws(%d, %d) : %f, x : %f, ps : %f\n ", ws_isx, col, row, ps_sum[ws_isx], x, T(edge_x*edge_y));

    if(false) // ( blockIdx.x == 0) && (blockIdx.y == 0) && row == 0 && col == 0)
        printf("weight[s=%d] : %f\n", s_x, weight);

    return x;
}


template<typename T>
__forceinline__ __device__ T L2_norm(T * ps_sum)
{
    int row = threadIdx.y;
    int col = threadIdx.x;
    int tIx = row*blockDim.x + col;


    ps_sum[tIx] *= ps_sum[tIx]; // Compute squares of entries.
    __syncthreads();

    int k_max = blockDim.x*blockDim.y/2;

    for (int k = k_max; k > 0; k/=2)
    {
        //  if (tIx == 0)  printf ("k in reduction : %d\n", k);
        if (tIx < k)
            ps_sum[tIx] += ps_sum[tIx+k];
        __syncthreads();
    }

    T result = ps_sum[0];
    __syncthreads();

    // Maybe we need the following for debugging again someday.
    if (false) // tIx == 0)
        printf ("ps_sum[%d] in reduction : %f, blockDim.x : %d\n", tIx, result, blockDim.x);

    return sqrt(result/(blockDim.x*blockDim.y)); // returns norm squared.
}

#define N_SCALES_2D 6
#define N_PX_X_2D 32
#define N_PX_Y_2D 32


#define N_SCALES_1D 10
#define N_PX_X_1D 512

/*
gnuplot script for computing weights for dyadic Dykstra without shifts

mu_S(S) = (S -.5)**.25
sigma_S(S) = sqrt(1./(8 * sqrt(S)))
q(a,S) = sqrt(2)*inverf(a**(1./S))
c_S(a,S) = (q(a,S) * sigma_S(S) + mu_S(S))**(-4)

do for [k=0:4] {; print c_S(4**k); }

using the image rows as largest subsets gives for a 512x512 image:

 do for [k=0:9] {; print c_S(.95, 512*2**k); }
0.00154391203759664
0.00081990498443775
0.000429039881840065
0.000221956966667552
0.000113828038818893
5.79917105865275e-05
2.93993246304473e-05
1.48495267507585e-05
7.48008622774254e-06
3.76036659965917e-06

These are the weights one has to use for the subsets within the rows. Therefore the factor of 512 in c_S.
*/


template<typename T, int offset_x, int offset_y>
struct DykstraStep {

    // In contrast to the formulation as pseudo code we only need a constant number of @p h variables because they are computed incrementally.

    // For the corrections we need the full history w.r.t. the number of levels in subset hierarchy.
    T Q[N_SCALES_2D];

    int n_scales;

    int  row, col, tIx, global_row, global_col, global_idx;

    __device__ DykstraStep (int n_s, const int width)
        :
          n_scales(n_s),

          row (threadIdx.y),
          col (threadIdx.x),
          tIx (row*blockDim.x + col), // This is for QTCreator's insufficient Intellisense.
          global_row (blockDim.y*blockIdx.y + row + offset_y),
          global_col (blockDim.x*blockIdx.x + col + offset_x),
          global_idx (global_row*width + global_col)
    {

        for (int j = 1; j <= n_scales; j++)
            Q[j-1] = 0;
    }

    __forceinline__
    __device__ void sum_of_squares_subsets(int s, T * ps_sum);


    __forceinline__ __device__ T project(int s_x, int s_y, T * ps_sum, T x,
                                         T g_noise, T weight)
    {
        // FIXME: move into kernel
        int edge_x = pow_of_two(s_x);
        int edge_y = pow_of_two(s_y);

        int ws_isx = blockDim.x*edge_y*(row/edge_y) +  edge_x*(col/edge_x);

        T ws = (ps_sum[ws_isx]/(g_noise*g_noise)) * weight;
        __syncthreads();

        if (ws > 1.)
            x *= 1.
                    / sqrt(ws);

        return x;
    }


    bool __device__ out_of_bound(const int width, const int height) const
    {
        // Let threads working outside the computational domain idle.
        if (
                (global_col < 0 || global_col >= width - offset_x)
                ||
                (global_row < 0 || global_row >= height - offset_y)
                )
            return true;
        else
            return false;
    }

    T __device__ sweep_fine_scales (
            T* h_old_in,  T* Q_full,
            const T* ICD_weights, T* ps_sum, const T g_noise)
    {
        volatile T * m_ps = ps_sum + tIx;

        Q[0] = Q_full[this->global_idx];

        T h_old = h_old_in[this->global_idx];
        T h_new = T(0);

        for (int j = 1; j <= this->n_scales; j++)         // loop over subsets
        {
            h_old -= this->Q[j-1];

            *m_ps  = h_old;

            __syncthreads();

            int s_x = n_scales - j; // we loop trough the scales from coarse to fine
            int s_y = s_x;

            this->sum_of_squares_subsets(s_x/*j-1*/, ps_sum);

            T weight =  // 0.25*
                    ICD_weights[ //
                    n_scales-j]; //   j-1]; // /sqrt(1.0*j); // n_scales - (j)]; // cs[j-1];

            h_new = this->project(s_x, //j-1 /* s_x*/,
                                  s_y, // j-1 /* s_y */,
                                  ps_sum, // this->h[j-1]
                                  h_old,
                                  g_noise, weight);

            __syncthreads();

            this->Q[j-1] = h_new - h_old;
            h_old = h_new;
        }
        return h_new; //  - h_init; // difference formed in on-chip registers
    }

};




//
template<typename T, int offset_x, int offset_y>
__forceinline__ __device__
void DykstraStep<T, offset_x, offset_y>::sum_of_squares_subsets(int s, T * ps_sum)
{
    // Pointer to first element of a row which (i.e. the element) is processed by this thread.
    volatile T * m_ps = ps_sum + tIx;

    // For debugging the 2D summation set to 1. Then each scale results in a power of 4 as the partial sum over the pixels of the subset.
    T x = // 1.0; //
            *m_ps;
    __syncthreads();

    *m_ps = x * x;
    __syncthreads();

    int edge_x = pow_of_two(s);

    // Each line is processed by one warp
    for (int e = 2; e <= edge_x; e*= 2)
    {
        if ((col+e) % e == 0)
        {
            *m_ps += *(m_ps + e/2);

        }
    }

    __syncthreads();

    int edge_y = pow_of_two(s);

    for (int e = 2; e <= edge_y; e*= 2) // add rows
    {
        __syncthreads(); // Needed for some reason to ensure thread safety.
        int neighbor = blockDim.x*(row + e/2) + col;
        if ((row+e) % e == 0)
        {
            *m_ps += ps_sum[neighbor];
        }

        __syncthreads();

    }

}





template<typename T, int offset_x, int offset_y>
__global__ void
__incomplete_dykstra_2D_fine_scales(
        T* h_iter, T* h_old, T* Q_full,
        const T g_noise,
        const int height, const int width, const int depth
        )
{
    const int n_scales = N_SCALES_2D;

    // 2D
    T ICD_weights[] = { 0.0181658,
                        0.0207498,
                        0.0155408,
                        0.00760004,
                        0.002742,
                        0.000825881};


    DykstraStep<T, offset_x, offset_y> dykstra( n_scales, width);

    T __shared__ ps_sum[N_PX_X_2D * N_PX_Y_2D];

    volatile T * m_ps = ps_sum + dykstra.tIx;

    *m_ps = T(0);

    __syncthreads();


    if (dykstra.out_of_bound(width, height))
        return;

    // Projection on scales suitable for intra-threadblock processing.
    h_iter[dykstra.global_idx] /* *m_ps_2*/ = dykstra.sweep_fine_scales(h_old, Q_full,
                                                                        ICD_weights, ps_sum, g_noise);


}


template<typename T, ParallelArch arch>
void step35::Kernels<T, arch>::dyadic_dykstra_fine_scale_part(
        T* h_iter, T* h_old, T* Q_full,
        const T g_noise,
        const int ni, const int nj, const int nk)
{
    int grid_2D_i= ni/N_PX_X_2D;
    int grid_2D_j= nj/N_PX_Y_2D;
    dim3 grid_2D(grid_2D_i, grid_2D_j);
    dim3 blocks_2D(N_PX_X_2D, N_PX_Y_2D);


    __incomplete_dykstra_2D_fine_scales<T, 0, 0><<<grid_2D, blocks_2D>>> (h_iter,
                                                                          h_old,  Q_full,
                                                                           g_noise,
                                                                           ni, nj, nk
                                                                           );

    getLastCudaError("__incomplete_dykstra_2D_fine_scales<T, 0, 0><<<>>> execution failed\n");
    cudaDeviceSynchronize();

}

template<typename T, ParallelArch arch>
void step35::Kernels<T, arch>::dyadic_dykstra_fine_scale_part_cpu(
        T* h_iter, T* h_old, T* Q_full,
        const T g_noise,
        const int ni, const int nj, const int nk)
{
    int grid_2D_i= ni/N_PX_X_2D;
    int grid_2D_j= nj/N_PX_Y_2D;
    dim3 grid_2D(grid_2D_i, grid_2D_j);
    dim3 blocks_2D(N_PX_X_2D, N_PX_Y_2D);

    T ps_sum[N_PX_X_2D][N_PX_Y_2D];
    T Q[N_SCALES_2D][N_PX_X_2D*N_PX_Y_2D];

    // TO DO: h_init = residual

    T ICD_weights[] = { 0.0181658,
                        0.0207498,
                        0.0155408,
                        0.00760004,
                        0.002742,
                        0.000825881};

#pragma omp parallel for
    for(uint bx = 0; bx  < grid_2D.x; bx++ )
    {
#pragma omp parallel for private(Q, ps_sum)
        for(uint by = 0; by < grid_2D.y; by++)
        {
            for(uint s = 1; s <= N_SCALES_2D; s++)
            {
//                printf("bx: %d, by: %d, s:%d \n", bx, by, s);
                uint s_x = N_SCALES_2D - s;

                uint s_y = s_x;

                uint edge_x = pow(2., s_x);
                uint edge_y = pow(2., s_y);

                T weight = ICD_weights[N_SCALES_2D - s];

                // loop over tiles in a block
                for(uint ti = 0; ti < blocks_2D.x; ti += edge_x)
                {
                    for(uint tj = 0; tj < blocks_2D.y; tj += edge_y)
                    {
                        for(uint ii = 0; ii < edge_x; ii++)
                            for(uint jj = 0; jj < edge_y; jj++)
                                ps_sum[ii][jj] = 0;

                        // loop over pixels in tile
                        for(uint ii = 0; ii < edge_x; ii++)
                        {
                            for(uint jj = 0; jj < edge_y; jj++)
                            {
                                // position of pixel wrt ps_sum
                                uint tile_idx = (tj + jj) * blocks_2D.x + ti + ii;
                                uint global_idx =
                                        by * blocks_2D.y * ni + ( tj + jj) * ni
                                        + bx * blocks_2D.x + ti + ii;

                                if(s == 1) // do that only on smallest level
                                    Q[0][tile_idx] = Q_full[global_idx];

                                T h_0 = h_old[global_idx] - Q[s - 1][tile_idx];

                                ps_sum[ti/edge_x][tj/edge_y] += h_0 * h_0;
                            }// jj end
                        }// ii end

                        // calculate weighted residua
                        T ws = sqrt(ps_sum[ti/edge_x][tj/edge_y] * weight) / g_noise;
                        if(ws > 1) // adapt pixels if statistic is violated
                        {
                            //adapt all pixels in the tile
                            for(uint ii = 0; ii < edge_x; ii++)
                                for(uint jj = 0; jj < edge_y; jj++)
                                {
                                    uint global_idx =
                                            by * blocks_2D.y * ni + ( tj + jj) * ni
                                            + bx * blocks_2D.x + ti + ii;
                                    h_iter[global_idx] = h_old[global_idx] / ws;
                                }
                        }

                        // update iteration variables for all pixels in a tile
                        for(uint ii = 0; ii < edge_x; ii++)
                            for(uint jj = 0; jj < edge_y; jj++)
                            {
                                uint tile_idx = (tj + jj) * blocks_2D.x + ti + ii;
                                uint global_idx =
                                        by * blocks_2D.y * ni + ( tj + jj) * ni
                                        + bx * blocks_2D.x + ti + ii;

                                Q[s - 1][tile_idx] = h_iter[global_idx] - h_old[global_idx];
                                h_old[global_idx] = h_iter[global_idx];
                            }

                    }//block x loop end
                }// block y loop end
            }// scale n end
        }//by end
    }// bx end
}



//@setc5{Kernel: __tv_derivative}
//@brief kernel for evaluation of the functional derivative of th TV regularization functional.
template<typename T>
__global__ void
__tv_derivative(T* dTV_du, const T* u, const T* f, T lambda, const int height, const int width, const int depth)
{
    T beta = 1e-4;

    // Arrays for the gradient at a given site.
    T grad_p[3/*dim*/];
    T grad_m[3];

    for (int d = 0; d < 3; d++)
    {
        grad_p[d] = 0.;
        grad_m[d] = 0.;
    }

    int row = threadIdx.y;
    int col = threadIdx.x;

    int global_row = blockDim.y*blockIdx.y + row;
    int global_col = blockDim.x*blockIdx.x + col;

    int site = global_row*width + global_col;

    int north = site + width;
    int south = site - width;

    int east = site + 1;
    int west = site - 1;

    // TO DO: 3D

    // The gradient is computed using forward differences.
    if (global_col < width-1)
        grad_p[0] = u[east]  - u[site];

    if (global_row < height-1)
        grad_p[1] = u[north] - u[site];
    // grad_p[2] = u[top] -u[site];

    if (global_col > 0)
        grad_m[0] = u[site] - u[west];

    if (global_row > 0)
        grad_m[1] = u[site] - u[south];

    // grad_m[2] = u[site] - u[bottom];

    // From the components of the gradient we can compute the local diffusion coefficient.
    T alpha_p = 0;
    T alpha_m = 0;

    for (int d = 0; d < 3; d++)
    {
        alpha_p += grad_p[d]*grad_p[d];
        alpha_m += grad_m[d]*grad_m[d];
    }

    alpha_p = 1./(sqrt(alpha_p + beta));
    alpha_m = 1./(sqrt(alpha_m + beta));

    // from the diffusion coefficient and the gradients
    // we finally compute the divergence weighted by
    // the regularization parameter $\lambda$.
    // To minimize multiplications we first add up the components
    // of the gradients and multiply with the diffusion coefficients only afterwards.
    T sum_grad_p = 0.;
    T sum_grad_m = 0.;

    for (int d = 0; d < 3; d++)
    {
        sum_grad_p += grad_p[d];
        sum_grad_m += grad_m[d];
    }

    dTV_du[site] = lambda * (
                (alpha_p *
                 sum_grad_p
                 -
                 alpha_m *
                 sum_grad_m)
                +
                2 * ( // f[site] -
                      u[site])
                );
}




//@setc5{Kernel: __tv_derivative}
//@brief kernel for evaluation of the functional derivative of th TV regularization functional.
template<typename T>
inline void
__tv_derivative_cpu(dim3 blockDim, dim3 blockId, T* dTV_du, const T* u,
                    const T* f, T lambda, const int height, const int width,
                    const int depth)
{
    T beta = 1e-4;

    // Arrays for the gradient at a given site.
    T grad_p[3/*dim*/];
    T grad_m[3];

    for (int d = 0; d < 3; d++)
    {
        grad_p[d] = 0.;
        grad_m[d] = 0.;
    }

    // CUDA
    // int row = threadIdx.y;
    // int col = threadIdx.x;

    // OpenMP
    for (int row = 0; row < blockDim.y; row++)
        for (int col = 0; col < blockDim.x; col++)
        {

            int global_row = blockDim.y*blockId.y + row;
            int global_col = blockDim.x*blockId.x + col;

            int site = global_row*width + global_col;

            int north = site + width;
            int south = site - width;

            int east = site + 1;
            int west = site - 1;

            // TO DO: 3D

            // The gradient is computed using forward differences.
            if (global_col < width-1)
                grad_p[0] = u[east]  - u[site];

            if (global_row < height-1)
                grad_p[1] = u[north] - u[site];
            // grad_p[2] = u[top] -u[site];

            if (global_col > 0)
                grad_m[0] = u[site] - u[west];

            if (global_row > 0)
                grad_m[1] = u[site] - u[south];

            // grad_m[2] = u[site] - u[bottom];

            // From the components of the gradient we can compute the local diffusion coefficient.
            T alpha_p = 0;
            T alpha_m = 0;

            for (int d = 0; d < 3; d++)
            {
                alpha_p += grad_p[d]*grad_p[d];
                alpha_m += grad_m[d]*grad_m[d];
            }

            alpha_p = 1./(sqrt(alpha_p + beta));
            alpha_m = 1./(sqrt(alpha_m + beta));

            // from the diffusion coefficient and the gradients
            // we finally compute the divergence weighted by
            // the regularization parameter $\lambda$.
            // To minimize multiplications we first add up the components
            // of the gradients and multiply with the diffusion coefficients only afterwards.
            T sum_grad_p = 0.;
            T sum_grad_m = 0.;

            for (int d = 0; d < 3; d++)
            {
                sum_grad_p += grad_p[d];
                sum_grad_m += grad_m[d];
            }

            dTV_du[site] = lambda * (
                        (alpha_p *
                         sum_grad_p
                         -
                         alpha_m *
                         sum_grad_m)
                        +
                        2 * ( // f[site] -
                              u[site])
                        );
        }
}


//@sect5{Function: tv_derivative}
//@brief This is a wrapper for the __tv_derivative Kernel.
//@param A pointer to the image
//@param ni height of the image
//@param nj width of the image
//@param offseti offset in vertical direction
//@param offsetj offset in horizontal direction
//@param smin minimum frame size is $2^{smin}$
template<typename T, ParallelArch arch>
void step35::Kernels<T, arch>::tv_derivative(T * dTV_du,
                                       const T *A_image,
                                       const T* f,
                                       const T lambda,
                                       const int ni, const int nj, const int nk)
{
    int gridsi= ni/N_PX_X_2D;
    int gridsj= nj/N_PX_Y_2D;
    dim3 grid(gridsi,gridsj);
    dim3 blocks(N_PX_X_2D, N_PX_Y_2D);

    // IF CUDA
    if (arch == gpu_cuda) {
    __tv_derivative<T><<<grid,blocks>>> (dTV_du,
                                         A_image,
                                         f,
                                         lambda,
                                         ni, nj, nk);
    getLastCudaError("__tv_derivative<<<>>> execution failed\n");
    cudaDeviceSynchronize();
    }
    else
    {
        // ELSE
#pragma omp parallel for
        for(uint bx = 0; bx  < grid.x; bx++ )
        {
#pragma omp parallel for
            for(uint by = 0; by < grid.y; by++)
            {
                dim3 blockId(bx,by);
                __tv_derivative_cpu<T> (blocks, blockId, dTV_du,
                                        A_image,
                                        f,
                                        lambda,
                                        ni, nj, nk);
            }
        }
    }
}



//@sect4{Kernel wrapper functions}
//The functions which call the kernels





//@sect5{Function: dyadic_dykstra}
//@brief This is a wrapper for the __dyadyc_dykstra Kernel function.
//@param A pointer to the image
//@param ni height of the image
//@param nj width of the image
//@param offseti offset in vertical direction
//@param offsetj offset in horizontal direction
//@param smin minimum frame size is $2^{smin}$



// Finally, we have to specialize the templates to force the compiler to actually compile something.
// This has to be at the end of file, because all functions have to be declared and their bodies defined before the
// class can be explictly instantiated by the compiler.
template class step35::Kernels<float, gpu_cuda>;
template class step35::Kernels<double, gpu_cuda>;
template class step35::Kernels<float, cpu>;
template class step35::Kernels<double, cpu>;

//void ImplCUDA<float>::apply<SDf, minus, SDf >(SDf&, SciPAL::DevBinaryExpr<SDf, minus, SDf > const&)", referenced from:

//      void SciPAL::LAOOperations::apply<float, cublas, Vc, BinX<Vcblas, minus, Vcblas > >(Vcblas&, SciPAL::Expr<BinX<Vcblas, minus, Vcblas > > const&)

