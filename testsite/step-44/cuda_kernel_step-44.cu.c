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


// Header containing the declarations of the wrapper functions - not the kernels.
// This header is the interface between the code compiled by nvcc and the one compiled by gcc.
#include <step-44/cuda_kernel_wrapper_step-44.cu.h>

// With the advent of the Fermi architecture
// it became possible to use the printf command from within a kernel.
// To enable this feature one has to include the corresponding header from the C standard library.
#include <stdio.h>

// @sect3{Device Functions}
//
// Prior to the kernels we have to define the device functions that we need.
//
// @sect4{Device Function: lex_index_2D}
//
// This function maps a row index $r$ and column index $c$
// of a matrix entry to the position $r\cdot row\_length + c$ in a linear array,
// which stores a matrix of row length
// $row\_length$ in row-major order.
// @param r : row index
// @param c : column index
// @param row_length : Zeilenlaenge
__device__ int lex_index_2D(int r, int c, int row_length)
{
    return c +  r*row_length;
}

// @sect3{Kernel}
//
// CUDA extends C by <i>Kernel</i>s, whicha re executed in parallel by
// <i>threads</i>. Furthermore, one can define <i>device functions</i>
// which can only be called from inside a kernel.
// <br>
// Note that unlike device functions kernels cannot be members of class.


// @sect4{Kernel: __dummy}
//
// This Kernel is just a mockup.
template<typename T>
__global__ void __step_44_dummy()
{
    int i = 2;
    int j = 3;
    int k = i+j;
    i += k;

    printf("Hi there from CUDA kernel. result : %d, thread Id : %d\n", k, threadIdx.x);
}


// @sect4{Function: step_44_dummy}
//
// Wrapper-Function for the dummy kernel.
// The naming convention is to prepend the name of the wrapper
// function by two underscores in order to get the name of the kernel.
template<typename T>
void step44::Kernels<T>::dummy()
{
    dim3 grid(2,1);
    dim3 blocks(16,1);

     printf("Hi from wrapper function\n");

    // Start the kernel. Note the template syntax.
    // Normally, it is not necessary to explicitly state the template parameters.
    // Here we have to because there are no function arguments from which the values
    // for the template paramter could be deduced.
    __step_44_dummy<T><<<grid, blocks>>>();

    cudaThreadSynchronize();
}


// Finally, we have to specialize the templates to force the compiler to actually compile something.
// This has to be at the end of file, because all functions have to be declared and their bodies defined before the
// class can be explictly instantiated by the compiler.
template class step44::Kernels<float>;
template class step44::Kernels<double>;
