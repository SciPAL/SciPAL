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

Copyright  Stephan Kramer, Johannes Hagemann 2014 - 2016, Lutz Künneke, Jan Lebert 2014 - 2015
*/

#ifndef CUDA_KERNEL_STEP_35_CU_H
#define CUDA_KERNEL_STEP_35_CU_H




//CUDA
#include <cuda.h>

namespace step35 {

// @sect4{Class: Kernels}
// Wrapper for our kernels
template<typename T, ParallelArch arch>
struct Kernels {

    void set_cs(T* cs_h, int maxnum);

    void dyadic_dykstra_fine_scale_part(T* h_iter, T* h_old, T* Q_full,
                                        const T g_noise,
                                        const int ni, const int nj, const int nk);

    void dyadic_dykstra_fine_scale_part_cpu(T* h_iter, T* h_old, T* Q_full,
                                            const T g_noise,
                                            const int ni, const int nj, const int nk);

    void tv_derivative(T *dTV_du,
                       const T *A_image, const T *f, const T lambda,
                       const int ni, const int nj, const int nk);
};

}

#endif // CUDA_KERNEL_STEP_35_CU_H
