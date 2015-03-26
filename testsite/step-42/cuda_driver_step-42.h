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


#ifndef CUDADriver_STEP_42_H
#define CUDADriver_STEP_42_H

#include <lac/release/blas_wrapper.hh>
#include <lac/release/cublas_wrapper.hh>

#include <numerics/SolverWrapper.h>

#include <base/CudaComplex.h>

// We encapsulate each project
// into a dedicated namespace
// in order to be able to re-use
// parts of a test program in others.
namespace step42 {

        // @sect3{Class: CUDADriver}
        //
        // This class manages the communication between host and device.
        // In particular the issue of memory transfers from and to the device.
        // The dummy implementation given here is supposed to give an
        // impression how this management could be done.
        // For worked out examples have a look at the other steps from
        // previous lab courses.
        // The documentation of the member functions is kept together
        // with their definitions.
class CUDADriver {

    typedef double Number;
    typedef SciPAL::CudaComplex<Number> cplxNumber;

    typedef cublas BW;

public:

    CUDADriver();

    ~CUDADriver() { BW::Shutdown(); }

    void gemm_tests();

    void gemv_tests();

    void complex_tests();

    void cusolver_demonstration();

    void feature_demonstration();

private:


};

} // namespace step42 END

#endif // CUDADriver_STEP_42_H
