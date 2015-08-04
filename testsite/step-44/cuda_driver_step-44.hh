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


#ifndef CUDA_DRIVER_STEP_44_HH
#define CUDA_DRIVER_STEP_44_HH


// The declaration of the interface to the CUDA-backend
// is contained in the following header.
#include <cuda_kernel_wrapper_step-44.cu.h>


// We have to include
#include <cuda_runtime_api.h>


        // @sect4{Constructor: CUDADriver}
        //
        // The constructor of the driver class allocates memory
        // on the GPU and copies data from host to device.
        // Furthermore, it keeps a pointer to the original host data for the
        // case that it has to be modified by the GPU results.
        // @param v_h : Pointer to a linear array in host-side memory that is to be copied to the GPU.
        // @param n : Number of entries of @p v_h.
step44::CUDADriver::CUDADriver(float * v_h, int n)
    :
    __n(n),
    __mem_size_v(sizeof(float)*n),
    __values_h(v_h)
{
            // The plain vanilla copy operation as it can be found in the <i>Programming Guide</i>.
            // You can do better than this by using <i>SciPal</i>. For an example of how to use SciPal
            // to hide details of the memory transfers have a look at the steps 4,5,6,7, and 12
            // from the lab course in 2011.
    cudaMalloc((void**)&__values_d, __mem_size_v);

    cudaMemcpy(__values_d, __values_h, __mem_size_v, cudaMemcpyHostToDevice) ;
}



        // @sect4{Destructor: CUDADriver}
        //
        // Free GPU memory and remove the pointer to the host data.
step44::CUDADriver::~CUDADriver()
{
   cudaFree(__values_d);
    __values_h = 0;
}



        // @sect4{Function: run}
        //
        // This function calls a kernel via the @p backend. The @p backend object encapsulates
        // the parallelization with CUDA.
void step44::CUDADriver::run()
{
           // Instantiate the @p backend. It could also be implemented as an attribute of this class.
           // In a worked out example the template parameter of @p Kernels will lead to CUDADriver being a template class as well.
    step44::Kernels<float> backend;

           // Delegate work to the @p backend and let the GPU do something ...
    backend.dummy();

           // For the sake of completeness we show how to get something back to the host.
    cudaMemcpy(__values_h, __values_d, __mem_size_v, cudaMemcpyDeviceToHost);
}


#endif // CUDA_DRIVER_STEP_44_HH
