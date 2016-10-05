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

#ifdef HAS_OPENMP
#include <omp.h>
#endif

//CUDA
#include <cuComplex.h>




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
    {/*nothing to do*/}

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


    bool
    __device__
    out_of_bound(const int width, const int height) const
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

    T __device__
    sweep_fine_scales ( T* h_old_in,  T* Q_full,
                        const T* ICD_weights, T* ps_sum, const T g_noise, int image_size)
    {
        volatile T * m_ps = ps_sum + tIx;


        T h_old = h_old_in[this->global_idx];
        T h_new = T(0);

        for (int j = 1; j <= this->n_scales; j++)         // loop over subsets
        {
            int Q_idx = this->global_idx + j*image_size;
            h_old -= Q_full[Q_idx];

            *m_ps  = h_old;

            __syncthreads();

            int s_x = n_scales - j; // we loop trough the scales from coarse to fine
//            int s_x = j-1; // we loop trough the scales from coarse to fine
            int s_y = s_x;

            this->sum_of_squares_subsets(s_x/*j-1*/, ps_sum);

            T weight =  // 0.25*
                    ICD_weights[ //
                    n_scales-j];
//               j-1]; // /sqrt(1.0*j); // n_scales - (j)]; // cs[j-1];

            h_new = this->project(s_x, //j-1 /* s_x*/,
                                  s_y, // j-1 /* s_y */,
                                  ps_sum, // this->h[j-1]
                                  h_old,
                                  g_noise, weight);

            __syncthreads();

            Q_full[Q_idx] = h_new - h_old;
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
                                                                        ICD_weights, ps_sum, g_noise, width*height);


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
typedef unsigned int uint;

#pragma omp parallel for
    for(unsigned int bx = 0; bx  < grid_2D.x; bx++ )
    {
#pragma omp parallel for private(Q, ps_sum)
        for(unsigned int by = 0; by < grid_2D.y; by++)
        {
            for(uint s = 1; s <= N_SCALES_2D; s++)
            {
//                printf("bx: %d, by: %d, s:%d \n", bx, by, s);
                int s_x = N_SCALES_2D - s;

                int s_y = s_x;

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
//@brief kernel for evaluation of the functional derivative of the TV regularization functional.
// Any previous content of @p dTV_du gets overwritten.
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
        for(unsigned int bx = 0; bx  < grid.x; bx++ )
        {
#pragma omp parallel for
            for(unsigned int by = 0; by < grid.y; by++)
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

//@setc5{Kernel: __L1_derivative}
//@brief kernel for evaluation of the functional derivative of the L1
// regularization functional.
// This function is dimension-independent.
// @param lambda : strength of regularization
template<typename T>
__global__ void
__L1_derivative(T* dL1_du, const T* u, T lambda, const int height, const
int width, const int depth)
{
    T margin= 1e-8;

    int row = threadIdx.y;
    int col = threadIdx.x;

    int global_row = blockDim.y*blockIdx.y + row;
    int global_col = blockDim.x*blockIdx.x + col;

    int site = global_row*width + global_col;

    T x_k = u[site];

    dL1_du[site] = lambda * (
               x_k / sqrt(x_k * x_k + margin)
                );
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
void step35::Kernels<T, arch>::L1_derivative(T * dL1_du,
                                       const T *u,
                                       const T lambda,
                                       const int ni, const int nj, const int nk)
{
    int gridsi= ni/N_PX_X_2D;
    int gridsj= nj/N_PX_Y_2D;
    dim3 grid(gridsi,gridsj);
    dim3 blocks(N_PX_X_2D, N_PX_Y_2D);

    // IF CUDA
    if (arch == gpu_cuda) {
    __L1_derivative<T><<<grid,blocks>>> (dL1_du,
                                         u,
                                         lambda,
                                         ni, nj, nk);
    getLastCudaError("__tv_derivative<<<>>> execution failed\n");
    cudaDeviceSynchronize();
    }
//    else
//    {
//        // ELSE
//#pragma omp parallel for
//        for(uint bx = 0; bx  < grid.x; bx++ )
//        {
//#pragma omp parallel for
//            for(uint by = 0; by < grid.y; by++)
//            {
//                dim3 blockId(bx,by);
//                __tv_derivative_cpu<T> (blocks, blockId, dTV_du,
//                                        A_image,
//                                        f,
//                                        lambda,
//                                        ni, nj, nk);
//            }
//        }
//    }
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

