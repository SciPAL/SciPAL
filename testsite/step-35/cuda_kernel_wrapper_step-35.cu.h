//@sect3{File: cuda_kernel_wrapper_step-35.cu.h}
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

#ifndef CUDA_KERNEL_STEP_35_CU_H
#define CUDA_KERNEL_STEP_35_CU_H


#define USE_TV_ITER
#undef USE_TV_ITER

#define nFMM_VERSION
#undef nFMM_VERSION

//CUDA
#include <cuda.h>

//SciPAL
#include <base/PrecisionTraits.h>
namespace step35 {

// @sect4{Class: Kernels}
// Wrapper for our kernels
template<typename T>
struct Kernels {

    void set_cs(T* cs_h, int maxnum);

    void dykstra(T *qg, T *A, int **cinfo, const int * q_offset_d, const int num_of_cluster, const int ni, const int nj, const int numpx, cudaStream_t *mystream);

    void dyadic_dykstra(T *A, 
                        const T * data, // new argument w.r.t. LJ
                        const T g_noise, // dto.
                        const int ni, const int nj, const int nk, const int n_max_steps, const T Tol);

    void dyadic_dykstra_fine_scale_part(
            T* h_iter, T* h_old, T* Q_full,
            T *A_image,
            // const T * data, // new argument w.r.t. LJ
            const T g_noise, // dto.
            const int ni, const int nj, const int nk, const int n_max_steps, const T Tol);

    void dyadic_dykstra_fine_scale_part_cpu(
            T* h_iter, T* h_old, T* Q_full,
            T *A_image,
            // const T * data, // new argument w.r.t. LJ
            const T g_noise, // dto.
            const int ni, const int nj, const int nk, const int n_max_steps, const T Tol);



    void dyadic_dykstra_LJ(T *A, const int ni, const int nj, const int nk, const int offseti, const int offsetj, const int offsetk, const int so);

    void element_norm_product(typename PrecisionTraits<T, gpu_cuda>::ComplexType *arr1, typename PrecisionTraits<T, gpu_cuda>::ComplexType *arr2, int width, int height, int depth);

    void reset(T* arr, int size);

    void tv_regularization(T* x,T* z,T* lag,T lambda,T rho,int nx,int ny,int nz);

    void haar(T* in_h, T* tmp_d, const int ny);

    void inverse_haar(T* in_h, T* tmp_d, const int ny);

    void abs(T* arr, int size);

    void soft_threshold(T* arr,T* lag,T* xvec,T rho, T threshold, int size);

    void sort(T *arr, int size);

    T    sum_all(T* arr,int size);

    void sparse(T* arr, T* lag, T rho, T gamma, int size);

    void pseudo_inverse(T* arr,T* lag,T rho,T gamma,int size);

    void soft_threshold_complex(typename PrecisionTraits<T, gpu_cuda>::ComplexType* arr,T threshold,int size);

    void tv_derivative(T *dTV_du, const T *A_image, const T *f, const T lambda, const int ni, const int nj, const int nk);
};

}

#endif // CUDA_KERNEL_STEP_35_CU_H
