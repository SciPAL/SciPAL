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


#ifndef FUJIMOTO_DRIVER_step2_H
#define FUJIMOTO_DRIVER_step2_H

#include <vector>
#include <lac/blas++.h>
#include <step-2/MVMultDriverInterface.h>
#include <step-2/cuda_kernel_wrapper_step-2.cu.h>

namespace step2 {

// @sect3{Class: FujimotoDriver}
//
// The most interesting variant is the matrix-vector product by Fujimoto.
// Although it is executed on the GPU its driver class is more similar to
// CpuBlasDriver because the management of the host-device memory transfers
// is moved into the wrapper function of the CUDA kernel.
// We use private inheritance for the Kernels<T> structure to express the "is implemented with" relationship.
template<typename T,typename blasType=cublas>
class FujimotoDriver : public MVMultDriverInterface<T>, private Kernels<T>
{
public:
    typedef
    typename MVMultDriverInterface<T>::FullMatrixAccessor FullMatrixAccessor;

    FujimotoDriver(const int v) : fj_version(v) {}

    virtual  ~FujimotoDriver () {}

    virtual double mvmult(std::vector<T>& y,
                          const FullMatrixAccessor& A,
                          const std::vector<T>& x,
                          int n_repetitions);

protected:
    const int fj_version;
};

} // namespace step2 END

#include <step-2/Fujimoto_driver_step-2.hh>
#endif // FUJIMOTO_DRIVER_step2_H
