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


#ifndef CUDA_KERNEL_STEP_2_CU_H
#define CUDA_KERNEL_STEP_2_CU_H

// @sect3{Declaration of CUDA Interface}
//
namespace step2 {


// @sect3{Class: Kernels}
//
// In our example programs we try to encapsulate CUDA as much as possible such that
// CUDA-independent parts of the source code can be reused easily.
// Especially, we want to localize CUDA-specific extensions of the C/C++language in only a few files.
// The bridge between these two parts is formed by the Kernels class.
// Its member functions mirror the arguments of the CUDA kernels and their bodies
// provide the encapsulation of the non-standard syntax for the kernels calls.
// Currently, there is only one wrapper function for the kernel by Fujimoto.
// Often kernels can be implemented in a type-independent way. Therefore, the Kernels
// class is templated with respect to the number type.
template<typename T>
struct Kernels {

    void mv_fujimoto(T *y, const T *A, const T *x,
                     const int m, const int n,
                     const int n_repetitions, const int fj_version, double& elapsed_time);
};
}
#endif // CUDA_KERNEL_STEP_2_CU_H
