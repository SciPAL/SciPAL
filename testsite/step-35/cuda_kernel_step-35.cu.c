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

Copyright  Lutz Künneke, Jan Lebert 2014
*/

//Header containing the declarations of the wrapper functions - not the kernels.
//This header is the interface between the code compiled by nvcc and the one compiled by gcc.

//std
#include <stdio.h>

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
template<typename T>
void step35::Kernels<T>::set_cs(T* cs_h, int maxnum) {
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

#ifndef nFMM_VERSION
    if (ws > 1.)
      x *= 1.
              / sqrt(ws); // x /= sqrt(ws);
#else
  // Dantzig version:  x = data + (x - data)*g_noise/sqrt(ws);
  //  x = x * (1. + 1./sqrt(weight*ps_sum[ws_isx]));
// JVIM version
    if ( ws > 1)
    x = data // paper : +
            +
            (x - data) // already in ws : * g_noise
            / sqrt(ws);
#endif
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

#define N_SCALES_2D 5
#define N_PX_X_2D 16
#define N_PX_Y_2D 16


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


template<typename T, bool horizontal>
__global__ void
__incomplete_dykstra_1D(T *residual,
                     const T * data, // new argument w.r.t. LJ
                     const T g_noise,
                     const int height, const int width, const int depth, const int n_max_steps, const T Tol)
{
    const int n_scales = N_SCALES_1D;

    T h[N_SCALES_1D + 1];
    T Q[N_SCALES_1D];

    T ICD_weights[] = {
        0.0271138027865195,
        0.0267867372407835,
        0.0249654121336563,
        0.0212163669068055,
        0.016364546216846,
        0.0115308932622277,
        0.00751081908224469,
        0.00458500994380723,
        0.00265856153102295,
        0.00148168722893069,
        0.000801553557730089,
        //
       0.244266654223811,
        0.123736396966113,
        0.0729229695774703,
        0.0445585553148736,
        0.0270981150609612,
        0.0161298109104501,
        0.00934412994711567,
        0.00526764212249119,
        0.00289717856415689,
        0.00154391203759664,   // 512 * 512x1
                        0.00081990498443775,   // 256
                        0.000429039881840065,  // 128
                        0.000221956966667552,  // 64
                        0.000113828038818893,  // 32
                        5.79917105865275e-05,  // 16
                        2.93993246304473e-05,  // 8
                        1.48495267507585e-05,  // 4
                        7.48008622774254e-06,  // 2
                        3.76036659965917e-06   // 1
                      };


    T __shared__ ps_sum[N_PX_X_1D];
    T __shared__ ps_sum_2[N_PX_X_1D];

    int row = threadIdx.y;
    int col = threadIdx.x;
    int idx = row*blockDim.x + col; // This is for QTCreator's insufficient Intellisense.


    volatile T * m_ps = ps_sum + idx;
    volatile T * m_ps_2 = ps_sum_2 + idx;

    *m_ps = T(0);
    *m_ps_2 = T(0);
    __syncthreads();

    int global_row = blockDim.y*blockIdx.y + row;
    int global_col = blockDim.x*blockIdx.x + col;

    int global_idx = global_row*width + global_col;

    if (!horizontal)
    {
        // swap row and column-
        int tmp = global_row;
        global_row = global_col;
        global_col = tmp;

        // position in global memory this thread is working on.
        global_idx = global_row*width + global_col;
    }



    volatile T * m_px = residual + global_idx; // Ptr to pixel this thread is responsible for.

    const T m_data = data[global_idx];

    h[0] =  *m_px;
    for (int j = 1; j <= n_scales; j++)
         Q[j-1] = 0;

    int iter = 0;
    bool iterate = true;
    while (iterate)
    {
        for (int j = 1; j <= n_scales; j++)         // loop over subsets
        {
             h[j-1] -= Q[j-1];
#ifndef nFMM_VERSION
            *m_ps  = h[j-1];
#else
             *m_ps  = h[j-1] - m_data;
#endif
            __syncthreads();

            T weight = // 0.25*
                    ICD_weights[ //j-1]; //
            n_scales - (j)]; // cs[j-1];

            L2_norm_on_sub_rows(n_scales - j // j-1
                                , ps_sum);

            h[j] = project(n_scales-j, // j-1 /* s_x*/,
                           0 /* s_y */,
                           ps_sum, h[j-1], m_data, g_noise, weight);

            __syncthreads();

            Q[j-1] = h[j] - h[j-1];
        }

        // For debugging L2_norm:
//        h[0] = -.5*sqrt(T(idx));;
//        h[n_scales] = .5*sqrt(T(idx));

        *m_ps_2 = /* 1.0; */ h[n_scales] - h[0];
        __syncthreads();

        T norm = L2_norm(ps_sum_2);

         if ( false) // (blockIdx.x == blockIdx.y) && (60 < idx && idx < 90))
            printf("iter : %d, norm : %f, h[0] : %f, h[1] : %f, h[2] : %f, h[3] : %f, h[4] : %f, h[5] : %f, h[6] : %f\n", iter, norm, h[0], h[1], h[2], h[3], h[4], h[5], h[6]);

        iterate = ( norm > Tol &&
                   iter < n_max_steps);
        iterate ?  h[0] = h[n_scales] : *m_px = h[n_scales];
        iter++;
    }

  //  if (blockIdx.x == blockIdx.y) //  && blockIdx.y == gridDim.y/2)
    //    printf("e(%d, %d) : %f\n ", col, row, h[n_scales]);
}


template<typename T, int offset_x, int offset_y>
struct DykstraStep {

    // In contrast to the formulation as pseudo code we only need a constant number of @p h variables because they are computed incrementally.
    T h_init, h_old, h_new;

    // For the corrections we need the full history w.r.t. the number of levels in subset hierarchy.
    T Q[N_SCALES_2D];

    int min_scale;
    int n_scales;

    int  row, col, tIx, global_row, global_col, global_idx;
#ifndef nFMM_VERSION
#else
     const T m_data;
#endif
    __device__ DykstraStep (int m, int n_s, const T* residual,
                        #ifndef nFMM_VERSION
                        #else
                            const T* data,
                        #endif
                            const int width)
        :
          min_scale(m),
          n_scales(n_s),

          row (threadIdx.y),
          col (threadIdx.x),
          tIx (row*blockDim.x + col), // This is for QTCreator's insufficient Intellisense.
          global_row (blockDim.y*blockIdx.y + row + offset_y),
          global_col (blockDim.x*blockIdx.x + col + offset_x),
          global_idx (global_row*width + global_col)
    #ifndef nFMM_VERSION
    #else
        ,
          m_data (data[global_idx])
    #endif
    {
        h_old = h_init = residual[global_idx];

        for (int j = 1; j <= n_scales; j++)
             Q[j-1] = 0;
    }


#ifdef blurb
    __forceinline__
#endif
     __forceinline__
    __device__ void sum_of_squares_subsets(int s, T * ps_sum)
#ifdef blurb
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
#else
    ;
#endif


    __forceinline__ __device__ T project(int s_x, int s_y, T * ps_sum, T x,
                                     #ifndef nFMM_VERSION
                                     #else
                                         T data,
                                     #endif
                                         T g_noise, T weight)
    {
        // FIXME: move into kernel
        int edge_x = pow_of_two(s_x);
        int edge_y = pow_of_two(s_y);

        int ws_isx = blockDim.x*edge_y*(row/edge_y) +  edge_x*(col/edge_x);

        T ws = (ps_sum[ws_isx]/(g_noise*g_noise)) * weight;
        __syncthreads();

    #ifndef nFMM_VERSION
        if (ws > 1.)
          x *= 1.
                  / sqrt(ws); // x /= sqrt(ws);
    #else
      // Dantzig version:  x = data + (x - data)*g_noise/sqrt(ws);
      //  x = x * (1. + 1./sqrt(weight*ps_sum[ws_isx]));
    // JVIM version
        if ( ws > 1)
        x = data // paper : +
                +
                (x - data) // already in ws : * g_noise
                / sqrt(ws);
    #endif

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


    // shared memory arrays cannot be attributes, hence we pass it as argument.
    T __device__ operator () (const T* ICD_weights, T* ps_sum, const T g_noise)
    {
        volatile T * m_ps = ps_sum + tIx;

        // FIXME: increasing min_scale leads to a loss of the large scales!
        for (int j = 1; j <= this->n_scales /* play here for increasing min_scale */; j++)         // loop over subsets
        {
                    h_old -= this->Q[j-1];
                    // The following makes a subtle difference.
                    if (j == this->min_scale + 1)
                        h_init = h_old;
#ifndef nFMM_VERSION
            *m_ps  = h_old;
#else
            *m_ps  = this->h[j-1] - m_data;
#endif
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
                   h_old ,
                      #ifndef nFMM_VERSION
                      #else
                                  m_data,
                      #endif
                                  g_noise, weight);

            __syncthreads();

            this->Q[j-1] = h_new - h_old;
            h_old = h_new;
        }
        return h_new - h_init;
    }


    T __device__ sweep_fine_scales (
            T* h_iter, T* h_old_in,  T* Q_full,
            const T* ICD_weights, T* ps_sum, const T g_noise)
    {
        volatile T * m_ps = ps_sum + tIx;

        Q[0] = Q_full[this->global_idx];

        this->h_old = h_old_in[this->global_idx];


        for (int j = this->min_scale + 1; j <= this->n_scales; j++)         // loop over subsets
        {
                    h_old -= this->Q[j-1];
#ifndef nFMM_VERSION
            *m_ps  = h_old;
#else
            *m_ps  = this->h[j-1] - m_data;
#endif
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
                   h_old ,
                      #ifndef nFMM_VERSION
                      #else
                                  m_data,
                      #endif
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
        T *residual,
    #ifndef nFMM_VERSION
    #else
        const T * data, // new argument w.r.t. LJ
    #endif
        const T g_noise,
        const int height, const int width, const int depth
        )
{
    const int n_scales = N_SCALES_2D;

    // 2D
    T ICD_weights[] = {
//        // alpha = 0.5
//        0.0247175052482826,
//        0.0234101860882775,
//        0.0157064920498531,
//        0.00733315498850165,
//        0.00262402918350882,
//        0.000796077044123706,
//        0.000220114811039996,
        // alpha = 0.1
        0.0284849245029229,
        0.0258315596691237,
        0.0167229685717504,
        0.00760602018197696,
        0.00267687492501606,
//        // alpha = 0.9
//        0.0202187993308712,
//        0.0203274835770538,
//        0.0143412644022621,
//        0.00695182573457273,
//        0.00254827749287444,
        /* 1.,
                        0.5,
                       0.25,
                        0.125,
                        0.0625,
                        0.03125,
                        0.015625,
                        0.00390625,*/
//                        0.0009765625,
                        // From extreme value theory:
      0.28, //
        //            1, //
        0.180406, // 1 subset
     //      0.14, //
                     //
        //             0.25, //
        0.0636429, // 4 subsets
       //     0.06, //
          //           0.1,  //
        0.0253601, // 16 subsets
         //   0.02, //
                  //   0.02, //
        0.00904285, // 64
          //  0.005, //
                        0.00284996, // 256
                    0.000818855, // 1024 16x16 tiles in 512x512 image
                    0.000221662, // 4096  8x8
    5.79917105865275e-05,  // 4x4
    1.48495267507585e-05,   // 2x2
    3.76036659965917e-06,   // 1x1
    9.46523485047832e-07 };


    DykstraStep<T, offset_x, offset_y> dykstra(1 /*min_scale*/, n_scales, residual,
                                           #ifndef nFMM_VERSION
                                           #else
                                               data,
                                           #endif
                                               width);

    T __shared__ ps_sum[// 4 *
            N_PX_X_2D * N_PX_Y_2D];

    volatile T * m_ps = ps_sum + dykstra.tIx;

    *m_ps = T(0);

    __syncthreads();


    if (dykstra.out_of_bound(width, height))
        return;

        // Projection on scales suitable for intra-threadblock processing.
    h_iter[dykstra.global_idx] /* *m_ps_2*/ = dykstra.sweep_fine_scales(h_iter, h_old, Q_full,
                ICD_weights, ps_sum, g_noise);


}


template<typename T, int offset_x, int offset_y>
__global__ void
__incomplete_dykstra_2D(T *residual,
                        #ifndef nFMM_VERSION
                        #else
                     const T * data, // new argument w.r.t. LJ
                        #endif
                     const T g_noise,
                     const int height, const int width, const int depth, const int n_max_steps, const T Tol)
{
    const int n_scales = N_SCALES_2D;

    // 2D
    T ICD_weights[] = {
//        // alpha = 0.5
//        0.0247175052482826,
//        0.0234101860882775,
//        0.0157064920498531,
//        0.00733315498850165,
//        0.00262402918350882,
//        0.000796077044123706,
//        0.000220114811039996,
        // alpha = 0.1
        0.0284849245029229,
        0.0258315596691237,
        0.0167229685717504,
        0.00760602018197696,
        0.00267687492501606,
//        // alpha = 0.9
//        0.0202187993308712,
//        0.0203274835770538,
//        0.0143412644022621,
//        0.00695182573457273,
//        0.00254827749287444,
        /* 1.,
                        0.5,
                       0.25,
                        0.125,
                        0.0625,
                        0.03125,
                        0.015625,
                        0.00390625,*/
//                        0.0009765625,
                        // From extreme value theory:
      0.28, //
        //            1, //
        0.180406, // 1 subset
     //      0.14, //
                     //
        //             0.25, //
        0.0636429, // 4 subsets
       //     0.06, //
          //           0.1,  //
        0.0253601, // 16 subsets
         //   0.02, //
                  //   0.02, //
        0.00904285, // 64
          //  0.005, //
                        0.00284996, // 256
                    0.000818855, // 1024 16x16 tiles in 512x512 image
                    0.000221662, // 4096  8x8
    5.79917105865275e-05,  // 4x4
    1.48495267507585e-05,   // 2x2
    3.76036659965917e-06,   // 1x1
    9.46523485047832e-07 };


    DykstraStep<T, offset_x, offset_y> dykstra(1 /*min_scale*/, n_scales, residual, // data,
                                               width);

    T __shared__ ps_sum[// 4 *
            N_PX_X_2D * N_PX_Y_2D];
    T __shared__ ps_sum_2[// 4 *
            N_PX_X_2D * N_PX_Y_2D];

    volatile T * m_ps = ps_sum + dykstra.tIx;
    volatile T * m_ps_2 = ps_sum_2 + dykstra.tIx;

    *m_ps = T(0);
    *m_ps_2 = T(0);
    __syncthreads();


    if (dykstra.out_of_bound(width, height))
        return;

    volatile T * m_px = residual + dykstra.global_idx; // Ptr to pixel this thread is responsible for.


    int iter = 0;
    bool iterate = true;
    while (iterate)
    {
        // Projection on scales suitable for intra-threadblock processing.
        *m_ps_2 = dykstra(ICD_weights, ps_sum, g_noise);

        // Convergence control.
        {
           T norm = L2_norm(ps_sum_2);

            iterate = ( norm > Tol &&
                        iter < n_max_steps);

            iterate ? dykstra.h_init = dykstra.h_new : *m_px = dykstra.h_new;
            iter++;
        }
    }
}


template<typename T, int offset_x, int offset_y>
__global__ void
__incomplete_dykstra_2D_old(T *residual,
                     const T * data, // new argument w.r.t. LJ
                     const T g_noise,
                     const int height, const int width, const int depth, const int n_max_steps, const T Tol)
{
    const int n_scales = N_SCALES_2D;

    T h[N_SCALES_2D + 1];
    T Q[N_SCALES_2D];

    // 2D
    T ICD_weights[] = {
//        // alpha = 0.5
//        0.0247175052482826,
//        0.0234101860882775,
//        0.0157064920498531,
//        0.00733315498850165,
//        0.00262402918350882,
//        0.000796077044123706,
//        0.000220114811039996,
        // alpha = 0.1
        0.0284849245029229,
        0.0258315596691237,
        0.0167229685717504,
        0.00760602018197696,
        0.00267687492501606,
//        // alpha = 0.9
//        0.0202187993308712,
//        0.0203274835770538,
//        0.0143412644022621,
//        0.00695182573457273,
//        0.00254827749287444,
        /* 1.,
                        0.5,
                       0.25,
                        0.125,
                        0.0625,
                        0.03125,
                        0.015625,
                        0.00390625,*/
//                        0.0009765625,
                        // From extreme value theory:
      0.28, //
        //            1, //
        0.180406, // 1 subset
     //      0.14, //
                     //
        //             0.25, //
        0.0636429, // 4 subsets
       //     0.06, //
          //           0.1,  //
        0.0253601, // 16 subsets
         //   0.02, //
                  //   0.02, //
        0.00904285, // 64
          //  0.005, //
                        0.00284996, // 256
                    0.000818855, // 1024 16x16 tiles in 512x512 image
                    0.000221662, // 4096  8x8
    5.79917105865275e-05,  // 4x4
    1.48495267507585e-05,   // 2x2
    3.76036659965917e-06,   // 1x1
    9.46523485047832e-07 };


    T __shared__ ps_sum[// 4 *
            N_PX_X_2D * N_PX_Y_2D];
    T __shared__ ps_sum_2[// 4 *
            N_PX_X_2D * N_PX_Y_2D];

    int row = threadIdx.y;
    int col = threadIdx.x;
    int idx = row*blockDim.x + col; // This is for QTCreator's insufficient Intellisense.

    int min_scale = 0;


    volatile T * m_ps = ps_sum + idx;
    volatile T * m_ps_2 = ps_sum_2 + idx;

    *m_ps = T(0);
    *m_ps_2 = T(0);
    __syncthreads();

    int global_row = blockDim.y*blockIdx.y + row + offset_y;
    int global_col = blockDim.x*blockIdx.x + col + offset_x;

    int global_idx = global_row*width + global_col;


    {
        // Let threads working outside the computational domain idle.
        if (global_col < 0 || global_col >= width - offset_x)
            return;

        if (global_row < 0 || global_row >= height - offset_y)
            return;
    }



    volatile T * m_px = residual + global_idx; // Ptr to pixel this thread is responsible for.

    const T m_data = data[global_idx];

    h[min_scale] =  *m_px;
    for (int j = 1; j <= n_scales; j++)
         Q[j-1] = 0;

    int iter = 0;
    bool iterate = true;
    while (iterate)
    {
#ifndef nOLD_IMPL
        for (int j = min_scale+1; j <= n_scales; j++)         // loop over subsets
        {
             h[j-1] -= Q[j-1];
#ifndef nFMM_VERSION
            *m_ps  = h[j-1];
#else
             *m_ps  = h[j-1] - m_data;
#endif
            __syncthreads();

            int s_x = n_scales - j; // we loop trough the scales from coarse to fine
            int s_y = s_x;

            sum_of_squares_on_2D_subsets(s_x/*j-1*/, ps_sum);

            T weight = ICD_weights[n_scales-j];



            h[j]   = project(s_x, //j-1 /* s_x*/,
                             s_y, // j-1 /* s_y */,
                             ps_sum, h[j-1], m_data, g_noise, weight);

            __syncthreads();
            if ((threadIdx.x == 1) && blockIdx.x == 1)
               printf("n_scales= %i, j= %i, weight= %f\n", n_scales, j, weight);

            Q[j-1] = h[j] - h[j-1];
        }
#else

#endif
        // For debugging L2_norm:
//        h[0] = -.5*sqrt(T(idx));;
//        h[n_scales] = .5*sqrt(T(idx));

        *m_ps_2 = /* 1.0; */ h[n_scales] - h[min_scale];
        __syncthreads();

        T norm = L2_norm(ps_sum_2);

         if ( false) // (blockIdx.x == blockIdx.y) && (60 < idx && idx < 90))
            printf("iter : %d, norm : %f, h[0] : %f, h[1] : %f, h[2] : %f, h[3] : %f, h[4] : %f, h[5] : %f, h[6] : %f\n", iter, norm, h[0], h[1], h[2], h[3], h[4], h[5], h[6]);

        iterate = ( norm > Tol &&
                   iter < n_max_steps);
        iterate ?  h[min_scale] = h[n_scales] : *m_px = h[n_scales];
        iter++;
    }

  //  if (blockIdx.x == blockIdx.y) //  && blockIdx.y == gridDim.y/2)
    //    printf("e(%d, %d) : %f\n ", col, row, h[n_scales]);
}


template<typename T>
void step35::Kernels<T>::dyadic_dykstra_fine_scale_part(
        T* h_iter, T* h_old, T* Q_full,
        T *A_image,
        //    const T * data, // new argument w.r.t. LJ
        const T g_noise, // dto.
        const int ni, const int nj, const int nk,
        const int n_max_steps=200, const T Tol=1e-4)
{
    int grid_2D_i= ni/N_PX_X_2D;
    int grid_2D_j= nj/N_PX_Y_2D;
    dim3 grid_2D(grid_2D_i, grid_2D_j);
    dim3 blocks_2D(N_PX_X_2D, N_PX_Y_2D);


    __incomplete_dykstra_2D_fine_scales<T, 0, 0><<<grid_2D, blocks_2D>>> (
                                                                           h_iter,  h_old,  Q_full,
                                                                           A_image,
                                                                           // data,
                                                                           g_noise,
                                                                           ni, nj, nk
                                                                           );

    getLastCudaError("__incomplete_dykstra_2D_fine_scales<T, 0, 0><<<>>> execution failed\n");
    cudaDeviceSynchronize();

}

//@sect5{Function: dyadyc_dykstra}
//@brief This is a wrapper for the __dyadyc_dykstra Kernel function.
//@param A pointer to the image
//@param ni height of the image
//@param nj width of the image
//@param offseti offset in vertical direction
//@param offsetj offset in horizontal direction
//@param smin minimum frame size is $2^{smin}$
template<typename T>
void step35::Kernels<T>::dyadic_dykstra(T *A_image,
                                        const T * data, // new argument w.r.t. LJ
                                        const T g_noise, // dto.
                                        const int ni, const int nj, const int nk, const int n_max_steps=200, const T Tol=1e-4)
{
    int grid_2D_i= ni/N_PX_X_2D;
    int grid_2D_j= nj/N_PX_Y_2D;
    dim3 grid_2D(grid_2D_i, grid_2D_j);
    dim3 blocks_2D(N_PX_X_2D, N_PX_Y_2D);

    int grid_1D_i= ni/N_PX_X_1D;
    int grid_1D_j= nj;
    dim3 grid_1D(grid_1D_i, grid_1D_j);
    dim3 blocks_1D(N_PX_X_1D);



    static const bool horizontal = true;

    // We first project on the rows. Since the largest subsets span through the whole image
    // this gives a fairly good exchange of information.

//    __incomplete_dykstra_1D<T, horizontal><<<grid_1D, blocks_1D>>> (A_image,
//                                                          data, g_noise,
//                                                          ni, nj, nk, n_max_steps, Tol);
//    getLastCudaError("__incomplete_dykstra_1D<<<>>> execution failed\n");
//    cudaDeviceSynchronize();

// return;

    // Then we project onto squares for a more isotropic equilibration.

    __incomplete_dykstra_2D<T, 0, 0><<<grid_2D, blocks_2D>>> (A_image,
                                                          #ifndef nFMM_VERSION
                                                          #else
                                                           data,
                                                          #endif
                                                           g_noise,
                                                           ni, nj, nk, n_max_steps, Tol);

    getLastCudaError("__incomplete_dykstra_2D(0,0)<<<>>> execution failed\n");

    cudaDeviceSynchronize();

    return;




    // Then on the columns.
    __incomplete_dykstra_1D<T, !horizontal><<<grid_1D, blocks_1D>>> (A_image,
                                                           data,
                                                           g_noise,
                                                           ni, nj, nk, n_max_steps, Tol);
    getLastCudaError("__incomplete_dykstra_1D<<<>>> execution failed\n");
    cudaDeviceSynchronize();



    __incomplete_dykstra_2D<T, 0, 4><<<grid_2D, blocks_2D>>> (A_image,
                                                         //  data,
                                                           g_noise,
                                                           ni, nj, nk, n_max_steps, Tol);

    getLastCudaError("__incomplete_dykstra_2D(0,1)<<<>>> execution failed\n");
    cudaDeviceSynchronize();


//    __incomplete_dykstra_1D<T, horizontal><<<grid_1D, blocks_1D>>> (A_image,
//                                                           data,
//                                                           g_noise,
//                                                           ni, nj, nk, n_max_steps, Tol);
//    getLastCudaError("__incomplete_dykstra_1D<<<>>> execution failed\n");
//    cudaDeviceSynchronize();


    __incomplete_dykstra_2D<T, 4, 0><<<grid_2D, blocks_2D>>> (A_image,
                                                         //  data,
                                                           g_noise,
                                                           ni, nj, nk, n_max_steps, Tol);

    getLastCudaError("__incomplete_dykstra_2D(1,0)<<<>>> execution failed\n");
    cudaDeviceSynchronize();


//    __incomplete_dykstra_1D<T, !horizontal><<<grid_1D, blocks_1D>>> (A_image,
//                                                           data,
//                                                           g_noise,
//                                                           ni, nj, nk, n_max_steps, Tol);
//    getLastCudaError("__incomplete_dykstra_1D<<<>>> execution failed\n");
//    cudaDeviceSynchronize();

    __incomplete_dykstra_2D<T, 0, -4><<<grid_2D, blocks_2D>>> (A_image,
                                                        //   data,
                                                           g_noise,
                                                           ni, nj, nk, n_max_steps, Tol);

    getLastCudaError("__incomplete_dykstra_2D(0,-1)<<<>>> execution failed\n");
    cudaDeviceSynchronize();


//    __incomplete_dykstra_1D<T, horizontal><<<grid_1D, blocks_1D>>> (A_image,
//                                                           data,
//                                                           g_noise,
//                                                           ni, nj, nk, n_max_steps, Tol);
//    getLastCudaError("__incomplete_dykstra_1D<<<>>> execution failed\n");
//    cudaDeviceSynchronize();

    __incomplete_dykstra_2D<T, -4, 0><<<grid_2D, blocks_2D>>> (A_image,
                                                        //   data,
                                                           g_noise,
                                                           ni, nj, nk, n_max_steps, Tol);

    getLastCudaError("__incomplete_dykstra_2D(-1,0)<<<>>> execution failed\n");
    cudaDeviceSynchronize();

//    __incomplete_dykstra_1D<T, !horizontal><<<grid_1D, blocks_1D>>> (A_image,
//                                                           data,
//                                                           g_noise,
//                                                           ni, nj, nk, n_max_steps, Tol);
//    getLastCudaError("__incomplete_dykstra_1D<<<>>> execution failed\n");
//    cudaDeviceSynchronize();



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




//@sect5{Function: dyadyc_dykstra}
//@brief This is a wrapper for the __dyadyc_dykstra Kernel function.
//@param A pointer to the image
//@param ni height of the image
//@param nj width of the image
//@param offseti offset in vertical direction
//@param offsetj offset in horizontal direction
//@param smin minimum frame size is $2^{smin}$
template<typename T>
void step35::Kernels<T>::tv_derivative(T * dTV_du,
                                       const T *A_image,
                                       const T* f,
                                       const T lambda,
                                       const int ni, const int nj, const int nk)
{
    int gridsi= ni/N_PX_X_2D;
    int gridsj= nj/N_PX_Y_2D;
    dim3 grid(gridsi,gridsj);
    dim3 blocks(N_PX_X_2D, N_PX_Y_2D);

    __tv_derivative<T><<<grid,blocks>>> (dTV_du,
                                         A_image,
                                         f,
                                         lambda,
                                         ni, nj, nk);
    getLastCudaError("__tv_derivative<<<>>> execution failed\n");
    cudaDeviceSynchronize();
}




//@brief CUDA adapted approximate Dykstra
//@param A pointer to the image
//@param ni height of the image
//@param nj width of the image
//@param offseti offset in vertical direction
//@param offsetj offset in horizontal direction
//@param smin minimum frame size is $2^{smin}$
template<typename T>
__global__ void
__incomplete_dykstra_LJ(T *A_image, const int ni, const int nj, const int nk, const int offseti, const int offsetj, const int offsetk, const int smin) {
    //Shared memory
    __shared__ float q[6144];
    __shared__ float f[1024];
    __shared__ float e[1024];

    // Set q = 0
    for (int s=5; s>=0; s--) {
        q[s*1024+threadIdx.x]=0;
    }

    // Temporary variables
    int is, js, ks, m;
    // Pixel index this thread is processing
    int idx;

    // Iteration counter
    int it=0;

    // FIXME: make adjustable from outside
    // Maximum iteration count
    const int itmax=100;

    // FIXME: make adjustable from outside
    // Tolerance for convergence check
    const T tol = 1e-4;
    // Increment
    T delta=2.0*tol;

    while ( delta > tol && it < itmax ) {
        it++;
        delta=0;
        // In each threadblock we apply dykstra's algorithm to
        // subsets with edge lengths 32, 16, 8, 4, 2 and 1.
        //\image html s5.png
        //\image html s4.png
        //\image html s3.png
        //\image html s2.png
        //\image html s1.png
        //\image html s0.png
        // Wait for every step before starting the iteration
        __syncthreads();
        for (int s=5; s>=smin; s--) {
            //Edge length of one subset = $2^s$
            int subset_length = pow_of_two(s);
            //Number of subsets in one threadblock = $2^{5-s}$
            int n_subsets =pow_of_two(5-s);
            //Number of pixels in one threadblock = $2^{2\cdot s}$
            int n_px_per_subset = pow_of_two(2*s);

            //i = blockInd + row * MaxBlock + col * MaxRow * MaxBlock \n
            //row= ( i / MaxBlock ) % MaxRow                          \n
            //col = ( i / MaxRow ) / MaxBlock
int tci = threadIdx.x/n_px_per_subset;
int tci_s = tci*n_px_per_subset;
            //Line in global image
            is = subset_length*(tci%n_subsets) +
                 (threadIdx.x%n_px_per_subset)%subset_length + blockIdx.x*32 + offseti;
            //Column in global image
            js = subset_length*(tci/n_subsets) +
                 (threadIdx.x%n_px_per_subset)/subset_length + blockIdx.y*32 + offsetj;
            ks=offsetk;
            is=is%ni;
            js=js%nj;
            //Pixel index in the current global image
            idx = ks*ni*nj + is*nj + js;

#ifdef DY_DEBUG
            if (idx%891==0)
            printf("A_%d : %f, ", idx, A_image[idx] );
#endif
            //Fill shared memory with variables we later use
            f[threadIdx.x] = A_image[idx] - q[s*1024+threadIdx.x];
            e[threadIdx.x] = f[threadIdx.x]*f[threadIdx.x];
            __syncthreads();

            m=1;

            //Sum over all pixels of one chunck, write result to e[tci_s]
            while (m <= pow_of_two(2*s-1)) {
                if (threadIdx.x - tci_s + m < n_px_per_subset) {
                    e[threadIdx.x] += e[threadIdx.x + m];
                }
                m = m << 1; // m = m*2
                __syncthreads();
            }

            //$q = x_{r+1} - x_r$
            q[s*1024+threadIdx.x] = -f[threadIdx.x];

            if (cs[s]*e[tci_s] > 1.0) {
                f[threadIdx.x] = f[threadIdx.x]/(sqrt(cs[s]*e[tci_s]));
            }
            // Update $q$
            q[s*1024+threadIdx.x] += f[threadIdx.x];

            // Templatized abs function
            __abs<T> mabs;
            // Calculate increment
            delta+=mabs(A_image[idx]-f[threadIdx.x]);
            // Update image
            A_image[idx] = f[threadIdx.x];
            __syncthreads();
        }
    }
   // if (threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0)
     //   printf("ICD convergence of block (%d,%d) in step %d : res_init = %f, res = %f, tol = %f\n", blockIdx.x, blockIdx.y, it, delta_init, delta, tol);
}




//@setc5{Kernel: __tv_regularization}
//@brief kernel for tv regularization
template<typename T>
__global__ void
__tv_regularization(T* x,T* z,T* lag,T lambda,T rho,int nx,int ny,int nz)
{
    int i = threadIdx.x + blockIdx.x*blockDim.x;
    int j = threadIdx.y + blockIdx.y*blockDim.y;
    int k = threadIdx.z + blockIdx.z*blockDim.z;
    int maxit=100;
    T sw=0.3;
    int it=0;
    T tol=1e-5;
    T delta=2.0*tol;
    __abs<T> mabs;

    // FIXME: For convergence control we need global data exchange
    T __shared__ ps_sum[N_PX_X_2D * N_PX_Y_2D];

    int row = threadIdx.y;
    int col = threadIdx.x;
    int idx = row*blockDim.x + col; // This is for QTCreator's insufficient Intellisense.

    volatile T * m_ps = ps_sum + idx;

    *m_ps = T(0);

    __syncthreads();

    int k_i_j = k*nx*ny+i*ny+j;

    // We need pointers to the pixel (site) this thread is responsible for ...
    volatile T * z_kij = z + k_i_j; // Ptr to pixel this threads is responsible for.

    // and pointers to the neighbors. We assume that the x axis point from
    // left to right and the y axis upwards when drawn on paper.
    // This gives rise to the following notion of east, west, south and north.
    // For 3D we add bottom (z<0) and top (z>0).
    volatile T * z_east = z + k_i_j - ny;
    volatile T * z_west = z + k_i_j + ny;

    volatile T * z_south = z + k_i_j - 1;
    volatile T * z_north = z + k_i_j + 1;

    volatile T * z_bottom = z + k_i_j - nx*ny;
    volatile T * z_top = z + k_i_j + nx*ny;

    // The following values do not get changed and thus are loaded only once.
    T x_kij = x[k_i_j];
     T lag_kij = lag[k_i_j];

    while ( delta > tol && it < maxit )
    {
        delta=0;
        it++;
        T newval;

#ifndef USE_OLD_VERSION
        newval= -rho*(x[k_i_j]-z[k_i_j]) - lag[k_i_j];
        if ( nx > 1 && i < nx-1 )
            newval=newval-lambda*mysign( z[k*nx*ny+(i+1)*ny+j] - z[k_i_j] ); // west
        if ( nx > 1 && i > 0 )
            newval=newval-lambda*mysign( z[k*nx*ny+(i-1)*ny+j] - z[k_i_j]); // east
        if ( ny > 1 && j < ny-1 )
            newval=newval-lambda*mysign( z[k_i_j+1] - z[k_i_j] ); // north
        if ( ny > 1 && j > 0  )
            newval=newval-lambda*mysign( z[k_i_j-1] - z[k_i_j] ); // south
        if ( nz > 1 && k < nz-1 )
            newval=newval-lambda*mysign( z[(k+1)*nx*ny+i*ny+j] - z[k_i_j] ); // top
        if ( nz > 1 && k > 0 )
            newval=newval-lambda*mysign( z[(k-1)*nx*ny+i*ny+j] - z[k_i_j] ); // bottom


        // delta=delta+mabs(newval); // convergence control by LJ

        __syncthreads();
        z[k_i_j]=z[k_i_j]-sw*newval;
        if ( z[k_i_j] < 0 )
            z[k_i_j]=0;
        *m_ps = mabs(newval);
        __syncthreads();
         delta = L2_norm(ps_sum);
#else
        newval = -rho * (x_kij - *z_kij) - lag_kij;

        T nn_contrib = 0.;
        T grad_norm = 0.;
        // Pixels not on the boundary of the image.
        if ( nx > 1 && i < nx-1 )
        {
            T du =  (*z_west - *z_kij);
            nn_contrib -= du;
            grad_norm += du*du;
        }
        if ( nx > 1 && i > 0 )
        {
             T du =  (*z_east - *z_kij);
            nn_contrib -= du;
            grad_norm += du*du;
        }
        if ( ny > 1 && j < ny-1 ) {
            T du =  (*z_north - *z_kij);
            nn_contrib -= du;
            grad_norm += du*du;
        }
        if ( ny > 1 && j > 0  )
        {
          T du =  (*z_south - *z_kij);
          nn_contrib -= du;
          grad_norm += du*du;
      }
        if ( nz > 1 && k < nz-1 ) {
             T du =  (*z_top - *z_kij);
             nn_contrib -= du;
             grad_norm += du*du;
         }
        if ( nz > 1 && k > 0 )
        {
          T du = (*z_bottom - *z_kij);
          nn_contrib -= du;
          grad_norm += du*du;
      }
        grad_norm = sqrt(grad_norm+1e-14);

        newval += (nn_contrib/grad_norm)*lambda;

        *m_ps = mabs(newval);
        __syncthreads();

        *z_kij -= sw*newval;
        if (*z_kij < 0 )
            *z_kij=0;
        __syncthreads();

        // We only control whether the TV regulrization has
        // converged within the subset of the image covered by this thread block.
        delta = L2_norm(ps_sum);

#endif
    }
}




//@sect5{Kernel: __soft_thresholding}
//@brief Applies a soft threshold to the 1 Norm of a real vector
//@param arr target array
//@param lag lagrangian for soft threholding in ADMM
//@param xvec source array
//@param rho weight of quadratic penalty term in ADMM
//@param gamma penalty for non-vanishing elements in 1 Norm
//@param size number of elements in array
template<typename T>
__global__ void
__soft_thresholding(T* arr,T* lag,T* xvec,T rho, T gamma, int size) {
    int idx=threadIdx.x+blockIdx.x*blockDim.x;
    if ( idx < size ) {
        if ( arr[idx] > 0 ) {
            arr[idx]=xvec[idx] +(lag[idx] - gamma)/rho;///rho)/(rho+lagr[idx]);
            if (arr[idx] < 0 )
                arr[idx]=0;
        } else {
            arr[idx]=xvec[idx] +(lag[idx] + gamma)/rho;///rho)/(rho+lagr[idx]);
            if (arr[idx] > 0 )
                arr[idx]=0;
        }
    }
}

//@sect5{Kernel: __soft_thresholding_complex}
//@brief Kernel for comlex soft thresholding to the 0 Norm of a complex vector
template<typename T>
__global__ void
__soft_thresholding_complex(typename PrecisionTraits<T, gpu_cuda>::ComplexType* arr,
                            T gamma, int size) {
    int idx=threadIdx.x+blockIdx.x*blockDim.x;
    if ( idx < size ) {
        if ( arr[idx].x > 0 ) {
            arr[idx].x=arr[idx].x - gamma;///rho)/(rho+lagr[idx]);
            if (arr[idx].x < 0 )
                arr[idx].x=0;
        } else {
            arr[idx].x=arr[idx].x + gamma;///rho)/(rho+lagr[idx]);
            if (arr[idx].x > 0 )
                arr[idx].x=0;
        }
        if ( arr[idx].y > 0 ) {
            arr[idx].y=arr[idx].y - gamma;///rho)/(rho+lagr[idx]);
            if (arr[idx].y < 0 )
                arr[idx].y=0;
        } else {
            arr[idx].y=arr[idx].y + gamma;///rho)/(rho+lagr[idx]);
            if (arr[idx].y > 0 )
                arr[idx].y=0;
        }
    }
}


//@sect5{Kernel: __sparse_thresholding}
//@brief Kernel for real valued soft thresholding to enforce 0 Norm sparsity in ADMM
//@param arr target array
//@lag Lagrangian
//@param rho quadratic penalty term for ADMM
//@param gamma penalty for non-vanishing entries in 0-Norm
//@parma size size of the array
template<typename T>
__global__ void
__sparse_thresholding(T* arr, T *lag, T rho, T gamma, int size) {
    int idx=threadIdx.x+blockIdx.x*blockDim.x;
    if ( idx < size ) {
        if ( arr[idx] > 0 ) {
            arr[idx]=arr[idx]-(gamma/rho)/(rho+lag[idx]);
            if ( arr[idx] < 0 )
                arr[idx]=0;
        } else {
            arr[idx]=arr[idx]+(gamma/rho)/(rho+lag[idx]);
            if ( arr[idx] > 0 )
                arr[idx]=0;
        }
    }
}

//@sect5{Kernel: __pseudo_inverse}
//@brief Kernel to minimize the 2-Norm of arr with a constraint. The solution corresponds to the pseudo-inverse solution of the constraint
template<typename T>
__global__ void
__pseudo_inverse(T* arr, T *lag, T rho, T gamma, int size) {
    int idx=threadIdx.x+blockIdx.x*blockDim.x;
    if ( idx < size ) {
        arr[idx]=(rho*arr[idx]+lag[idx])/(rho+gamma);
    }
}
//@sect5{Haar Wavelet Transformation}
//@sect6{Kernel: __haar_horizontal}
//@brief Calculates 1D Haar Wavelet Transformation on all rows for one step with step width $w$.
//@see <a href="http://unix4lyfe.org/haar/">Here</a> and <a href="http://www.tilman.de/uni/ws05/scivis/2d-transformation.html">here</a> for a complete Haar Wavelet Transformation algorithm details.
//@param in input
//@param out output
//@param ny height of the image
//@param w step width for one step

//Treat the first $w$ elements of each row as $w/2$ pairs called $(a, b)$. \n
//Calculate $(a + b) / \sqrt{2}$ for each pair, these values will be the first half of the considered of elements in the output array. \n
//Calculate $(a - b) / \sqrt{2}$ for each pair, these values will be the second half of the considered of elements in the output array.\n
//For $w = 8$ and input array \n
//\htmlonly <pre>[  1.000   2.000   3.000   1.000   2.000   3.000   4.000   0.000  ]</pre>\endhtmlonly
//will result into \htmlonly <pre>[  2.121   2.828   3.536   2.828  -0.707   1.414  -0.707   2.828  ]</pre>\endhtmlonly
//Treating this output array as input array with $w=4$ will yield \n
//\htmlonly <pre>[  3.500   4.500  -0.500   0.500  -0.707   1.414  -0.707   2.828  ]</pre>\endhtmlonly
template<typename T>
__global__ void __haar_horizontal(T* in, T* out, const int ny, int w) {
    int x_idx = threadIdx.x + blockIdx.x*blockDim.x;
    int y_idx = threadIdx.y + blockIdx.y*blockDim.y;

    if (x_idx >= (w+1)/2 || y_idx >= w)
        return;

    int in_idx  = 2*x_idx + y_idx*ny;
    int sum_idx =   x_idx + y_idx*ny;

    const T one_over_sqrt_2 = 0.70710677; //= 1/sqrt(2)
    out[sum_idx]       = (in[in_idx] + in[in_idx+1])*one_over_sqrt_2;
    out[sum_idx + w/2] = (in[in_idx] - in[in_idx+1])*one_over_sqrt_2;
}

//@sect6{Kernel: __inverse_haar_horizontal}
//@brief Inverse 1D row Haar Wavelet Transformation step. Inverts __haar_horizontal
template<typename T>
__global__ void __inverse_haar_horizontal(T* in, T* out, const int ny, int w) {
    int x_idx = threadIdx.x + blockIdx.x*blockDim.x;
    int y_idx = threadIdx.y + blockIdx.y*blockDim.y;

    if (x_idx >= (w+1)/2 || y_idx >= w)
        return;

    int sum_idx      =   x_idx + y_idx*ny;
    int original_idx = 2*x_idx + y_idx*ny;

    const T one_over_sqrt_2 = 0.70710677; //= 1/sqrt(2)
    // $a = ((a+b)/\sqrt{2}+(a-b)/\sqrt{2})/\sqrt{2}$
    out[original_idx]   = (in[sum_idx] + in[sum_idx + w/2])*one_over_sqrt_2;
    // $b = ((a+b)/\sqrt{2}-(a-b)/\sqrt{2})/\sqrt{2}$
    out[original_idx+1] = (in[sum_idx] - in[sum_idx + w/2])*one_over_sqrt_2;
}


//@sect6{Kernel: __haar_vertical}
//@brief Calculates 1D Haar Transformation on all columns for one step with width $w$.
//@see __haar_horizontal
template<typename T>
__global__ void __haar_vertical(T* in, T* out, const int ny, int w) {
    int x_idx = threadIdx.x + blockIdx.x*blockDim.x;
    int y_idx = threadIdx.y + blockIdx.y*blockDim.y;

    if (y_idx >= (w+1)/2 || x_idx >= w)
        return;

    int in_idx1 = x_idx + 2*y_idx*ny;
    int in_idx2 = x_idx + (2*y_idx+1)*ny;
    int sum_idx = x_idx + y_idx*ny;

    const T one_over_sqrt_2 = 0.70710677; //= 1/sqrt(2)
    out[sum_idx]        = (in[in_idx1]+in[in_idx2])*one_over_sqrt_2;
    out[sum_idx+ny*w/2] = (in[in_idx1]-in[in_idx2])*one_over_sqrt_2;
}

//@sect6{Kernel: __inverse_haar_vertical}
//@brief Inverse 1D column Haar Wavelet Transformation step. Inverts __haar_vertical
template<typename T>
__global__ void __inverse_haar_vertical(T* in, T* out, const int ny, int w) {
    int x_idx = threadIdx.x + blockIdx.x*blockDim.x;
    int y_idx = threadIdx.y + blockIdx.y*blockDim.y;

    if (y_idx >= (w+1)/2 || x_idx >= w)
        return;

    int original_idx1 = x_idx + 2*y_idx*ny;
    int original_idx2 = x_idx + (2*y_idx+1)*ny;
    int sum_idx       = x_idx + y_idx*ny;

    const T one_over_sqrt_2 = 0.70710677; //= 1/sqrt(2)
    // $a = ((a+b)/\sqrt{2}+(a-b)/\sqrt{2})/\sqrt{2}$
    out[original_idx1] = (in[sum_idx]+in[sum_idx+ny*w/2])*one_over_sqrt_2;
    // $b = ((a+b)/\sqrt{2}-(a-b)/\sqrt{2})/\sqrt{2}$
    out[original_idx2] = (in[sum_idx]-in[sum_idx+ny*w/2])*one_over_sqrt_2;
}

//@sect4{Kernel wrapper functions}
//The functions which call the kernels

//@sect5{Haar Wavelet Transformation}

//@sect6{Function: haar}
//@brief Calculates full 2D nonstandard Haar Wavelet decomposition.
//@param in_d pointer to the array on device, will be transformed inplace. Has to be of quadratic size.
//@param tmp_d pointer to an array on device with same size as in_d. Will be filled with garbage.
//@param ny edge length of the matrix, has to be of the power of two

//Alternating row and column 1D Haar Wavelet steps
template<typename T>
void step35::Kernels<T>::haar(T* in_d, T* tmp_d, const int ny) {
    dim3 blocksize(32,32);
    dim3 gridsize;

    int w = ny;
    gridsize.x=(w+blocksize.x-1)/blocksize.x;
    gridsize.y=(w+blocksize.y-1)/blocksize.y;

    while(w>1) {
        __haar_horizontal<T><<<gridsize,blocksize>>>(in_d, tmp_d, ny, w);
        __haar_vertical<T><<<gridsize,blocksize>>>(tmp_d, in_d, ny, w);
        w /= 2;
    }
    cudaDeviceSynchronize();
    getLastCudaError("haar() execution failed\n");
}

//@sect6{Function: inverse_haar}
//@brief Inverse 2D nonstandard Haar Wavelet decomposition algorithm
//@param in_d pointer to the array on device, will be transformed inplace. Has to be of quadratic size.
//@param tmp_d pointer to an array on device with same size as in_d. Will be filled with garbage.
//@param ny edge length of the matrix, has to be of the power of two
//@see haar()
template<typename T>
void step35::Kernels<T>::inverse_haar(T* in_d, T* tmp_d, const int ny) {
    dim3 blocksize(32,32);
    dim3 gridsize;

    int w = 2;
    gridsize.x=(ny+blocksize.x-1)/blocksize.x;
    gridsize.y=(ny+blocksize.y-1)/blocksize.y;

    while(w <= ny) {
        __inverse_haar_horizontal<T><<<gridsize,blocksize>>>(in_d, tmp_d, ny, w);
        __inverse_haar_vertical<T><<<gridsize,blocksize>>>(tmp_d, in_d, ny, w);
        w *= 2;
    }
    cudaDeviceSynchronize();
    getLastCudaError("inverse_haar() execution failed\n");
}

//@sect5{Function: dykstra}
//@brief This is a wrapper for the __dykstra Kernel function.
//@param qg pointer to the q vectors on device
//@param A pointer to the image on device
//@param cinfo array of device pointers where to find n,i,j and numpatch
//@param q_offset array of offsets where to find first $q$ value for each cluster
//@param num_of_cluster numbers of clusters (threadblocks)
//@param ni height of the image
//@param nj width of the image
//@param maximum number of pixels per threadblock
//@param mystream a CUDA stream
template<typename T>
void step35::Kernels<T>::dykstra(T *qg, T *A,  int **cinfo, const int * q_offset, const int num_of_cluster,
                                 const int ni, const int nj, const int numpx, cudaStream_t *mystream) {
    dim3 grid(num_of_cluster,1);
    dim3 blocks(numpx,1);
    __dykstra<T><<<grid, blocks, 3*numpx*sizeof(T), *mystream>>> (qg, A, cinfo, q_offset, ni, nj);
}


//@sect5{Function: dyadyc_dykstra}
//@brief This is a wrapper for the __dyadyc_dykstra Kernel function.
//@param A pointer to the image
//@param ni height of the image
//@param nj width of the image
//@param offseti offset in vertical direction
//@param offsetj offset in horizontal direction
//@param smin minimum frame size is $2^{smin}$
template<typename T>
void step35::Kernels<T>::dyadic_dykstra_LJ(T *A_image, const int ni, const int nj, const int nk, const int offseti, const int offsetj, const int offsetk, const int so) {
    int gridsi=1; // ni/32;
    int gridsj=1; // nj/32;
    dim3 grid(gridsi,gridsj);
    dim3 blocks(128,1);
    __incomplete_dykstra_LJ<T><<<grid,blocks>>> (A_image, ni, nj, nk, offseti, offsetj, offsetk, so);
    getLastCudaError("__incomplete_dykstra<<<>>> execution failed\n");
    cudaDeviceSynchronize();
}


// @sect5{Function: element_norm_product}
//@brief Elementwise multiplication of two arrays with complex values and division by the number of elements in each array.
//@param arr1 first input array, output will be calculated inplace here
//@param arr2 second input array
//@param width width of the array
//@param height height of the array
template<typename T>
void step35::Kernels<T>::element_norm_product(typename PrecisionTraits<T, gpu_cuda>::ComplexType* arr1,
        typename PrecisionTraits<T, gpu_cuda>::ComplexType* arr2,
        int width, int height, int depth) {
    int size=depth*width*(height/2+1);
    T normalization=width*height*depth;
    __element_norm_product<T> op(normalization);
#ifndef nUSE_CPP11
    thrust::device_ptr<typename PrecisionTraits<T, gpu_cuda>::ComplexType> ptr1 = thrust::device_pointer_cast(arr1);
    thrust::device_ptr<typename PrecisionTraits<T, gpu_cuda>::ComplexType> ptr2 = thrust::device_pointer_cast(arr2);
    thrust::transform(thrust::device, ptr1, ptr1 + size, ptr2, ptr1, op);
#endif
    cudaDeviceSynchronize();
    getLastCudaError("element_norm_product<<<>>> execution failed\n");
}



//@sect5{Function: reset}
//@brief Fills input array with zeros
//@param arr input array
//@param size number of elements in the input array
template<typename T>
void step35::Kernels<T>::reset(T* arr, int size) {
#ifndef nUSE_CPP11
    thrust::device_ptr<T> ptr = thrust::device_pointer_cast(arr);
    thrust::fill(thrust::device, ptr, ptr + size, 0);
#endif
    cudaDeviceSynchronize();
    getLastCudaError("thrust::fill execution failed\n");
}

//@sect5{Function: abs}
//@brief Replaces every element in an array with it the absolute value of the element
//@param arr input array
//@param size number of elements in array
template<typename T>
void step35::Kernels<T>::abs(T* arr, int size) {
    __abs<T> op;
#ifndef nUSE_CPP11
    thrust::device_ptr<T> ptr = thrust::device_pointer_cast(arr);
    thrust::transform(thrust::device, ptr, ptr + size, ptr, op);
#endif
    cudaDeviceSynchronize();
    getLastCudaError("thrust::transform __abs execution failed\n");
}

#ifdef USE_TV_ITER
//@sect5{Function: tv_regularization}
//@brief TODO
template<typename T>
void step35::Kernels<T>::tv_regularization(T* x, T* z, T* lag, T lambda, T rho, int nx, int ny, int nz) {
    int gridsi=nx/N_PX_X_2D;
    int gridsj=ny/N_PX_Y_2D;
    int gridsk=nz;
    dim3 grid(gridsi,gridsj,gridsk);
    dim3 blocks(N_PX_X_2D, N_PX_Y_2D, 1);
    __tv_regularization<T><<<grid,blocks>>>(x,z,lag,lambda,rho,nx,ny,nz);
    cudaDeviceSynchronize();
    getLastCudaError("__tv_regularization<<<>>> execution failed\n");
}
#endif

//@sect5{Function: soft_threshold}
//@brief Applies a soft threshold. Absolute values smaller than a certain threshold will be to zero,
//       larger values will be shifted by the value of the threshold closer to zero.
//@param arr device array to sort, result will be inplace
//@param threshold threshold value
//@param size number of elements in array
template<typename T>
void step35::Kernels<T>::soft_threshold(T* arr, T* lag,T *xvec, T rho, T threshold, int size) {
    int grids=size/1024+1;
    dim3 grid(grids,1);
    dim3 blocks(1024,1);
    __soft_thresholding<T><<<grid,blocks>>>(arr, lag,xvec,rho, threshold, size);
    cudaDeviceSynchronize();
    getLastCudaError("__soft_threshold<<<>>> execution failed\n");
}

//@sect5{Function: soft_threshold_complex}
//@brief Soft thresholding for complex fields
//@param threshold threshold for both complex and real part
//@param arr inplace thresholding on arr
//@param size size of arr
template<typename T>
void step35::Kernels<T>::soft_threshold_complex(typename PrecisionTraits<T, gpu_cuda>::ComplexType* arr,
        T threshold, int size) {
    int grids=size/1024+1;
    dim3 grid(grids,1);
    dim3 blocks(1024,1);
    __soft_thresholding_complex<T><<<grid,blocks>>>(arr,threshold,size);
    cudaDeviceSynchronize();
    getLastCudaError("__soft_thresholding_complex<<<>>> execution failed\n");
}


//@sect5{Function: sparse}
//@brief Corresponds to soft thresholding in real fields
//@param gamma thresholding parameter
//@param arr array to perform inplace thresholding on
//@param lag Lagrangian
//@param rho parameter rho from ADMM
//@param size size of arr
template<typename T>
void step35::Kernels<T>::sparse(T* arr, T* lag, T rho, T gamma, int size) {
    int grids=size/1024+1;
    dim3 grid(grids,1);
    dim3 blocks(1024,1);
    __sparse_thresholding<T><<<grid,blocks>>>(arr,lag,rho,gamma,size);
    cudaDeviceSynchronize();
    getLastCudaError("sparse execution failed\n");
}
template<typename T>

//@sect5{Function: pseudo_inverse}
//@brief Pseudoinverse
void step35::Kernels<T>::pseudo_inverse(T* arr, T* lag, T rho, T gamma, int size) {
    int grids=size/1024+1;
    dim3 grid(grids,1);
    dim3 blocks(1024,1);
    __pseudo_inverse<T><<<grid,blocks>>>(arr,lag,rho,gamma,size);
    cudaDeviceSynchronize();
    getLastCudaError("__pseudo_inverse execution failed\n");
}

//@sect5{Function: sort}
//@brief Sorts array using thrust::sort
//@param arr device array to sort, result will be inplace
//@param size number of elements in array
//@brief Wrapper around thrust::sort
template<typename T>
void step35::Kernels<T>::sort(T* arr,int size) {
#ifndef nUSE_CPP11
    thrust::device_ptr<T> ptr = thrust::device_pointer_cast(arr);
    thrust::sort(thrust::device, ptr, ptr + size);
#endif
    cudaDeviceSynchronize();
    getLastCudaError("thrust::sort execution failed\n");
}

//@sect5{Function: sum_all}
//@return Sum over all elements
//@param arr device array
//@param size number of elements in array
//@brief Wrapper around thrust::reduce
template<typename T>
T step35::Kernels<T>::sum_all(T* arr,int size) {
#ifndef nUSE_CPP11
    thrust::device_ptr<T> ptr = thrust::device_pointer_cast(arr);
    T sum = thrust::reduce(thrust::device, ptr, ptr + size, 0, thrust::plus<T>());
#else
    T sum = -1;
#endif
    cudaDeviceSynchronize();
    getLastCudaError("thrust::reduce plus execution failed\n");
    return sum;
}

// Finally, we have to specialize the templates to force the compiler to actually compile something.
// This has to be at the end of file, because all functions have to be declared and their bodies defined before the
// class can be explictly instantiated by the compiler.
template class step35::Kernels<float>;
template class step35::Kernels<double>;

//void ImplCUDA<float>::apply<SDf, minus, SDf >(SDf&, SciPAL::DevBinaryExpr<SDf, minus, SDf > const&)", referenced from:

//      void SciPAL::LAOOperations::apply<float, cublas, Vc, BinX<Vcblas, minus, Vcblas > >(Vcblas&, SciPAL::Expr<BinX<Vcblas, minus, Vcblas > > const&)

