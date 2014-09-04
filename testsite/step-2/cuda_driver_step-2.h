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


#ifndef CUBlasDriver_STEP_2_H
#define CUBlasDriver_STEP_2_H

#include <vector>
#include <lac/blas++.h>
#include <step-2/MVMultDriverInterface.h>

namespace step2 {

// @sect3{Class: CUBlasDriver}
//
// In contrast to the CpuBlasDriver class which has just wrapped
// the blas function into the mvmult() function for providing
// a more comfortable progammer's interface
// this class additionally has to manage the host-device-communication,
// i.e. the data transfer between the two different memories.
template<typename T,typename blasType=cublas>
class CUBlasDriver : public MVMultDriverInterface<T>
{
public:
    typedef
    typename MVMultDriverInterface<T>::FullMatrixAccessor
    FullMatrixAccessor;


    CUBlasDriver();

    virtual ~CUBlasDriver ();

    virtual double mvmult(std::vector<T>& y,
                          const FullMatrixAccessor& A,
                          const std::vector<T>& x,
                          int	n_repetitions);
};

} // namespace step2 END

#include <step-2/cuda_driver_step-2.hh>
#endif // CUBlasDriver_STEP_2_H

