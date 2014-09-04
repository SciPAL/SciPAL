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


#ifndef FUJIMOTO_DRIVER_step2_HH
#define FUJIMOTO_DRIVER_step2_HH

#include <step-2/cuda_kernel_wrapper_step-2.cu.h>
#include <base/CUDATimer.h>

// @sect4{Function: mvmult}
//
// For details about the function arguments see the base class.
// It is very similar to CpuBlasDriver::mvmult(). The major differences are
// that instead of a blas function a CUDA kernel is invoked and that the measurement of the runtime
// is moved into the kernel wrapper function which gives us a cumulative runtime and not an average.
template<typename Number,typename blasType>
double step2::FujimotoDriver<Number,blasType>::mvmult(std::vector<Number> & y,
                                                      const FullMatrixAccessor& A,
                                                      const std::vector<Number> & x,
                                                      int n_repetitions)
{
    int n_rows = A.n_rows();
    int n_cols = A.n_cols();

    Number * dst = &y[0];

    const Number * A_entries = A.val();

    const Number * src = &x[0];


    double cumulative_elapsed_time = 0;

    this->mv_fujimoto(dst, A_entries, src,
                      n_rows, n_cols,
                      n_repetitions,
                      this->fj_version,
                      cumulative_elapsed_time);

    return cumulative_elapsed_time;
}

#endif // FUJIMOTO_DRIVER_step2_HH
