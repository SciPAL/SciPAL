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


#ifndef CPUBLAS_DRIVER_STEP2_H
#define CPUBLAS_DRIVER_STEP2_H

#include <vector>
#include <lac/blas++.h>
#include <step-2/MVMultDriverInterface.h>

namespace step2 {

// @sect3{Class: CpuBlasDriver}
//
// The first of our driver classes implements the CPU-based
// matrix-vector multiplication.
template<typename EntryType, typename blasType>
class CpuBlasDriver : public MVMultDriverInterface<EntryType>
{
public:

    typedef
    typename MVMultDriverInterface< EntryType>::FullMatrixAccessor
    FullMatrixAccessor;

    CpuBlasDriver() {}

    virtual ~CpuBlasDriver () {}

    // The only function of interest of this class is the
    // one which profiles the matrix-vector multiplication.
    virtual double mvmult(std::vector<EntryType>& y,
                          const FullMatrixAccessor& A,
                          const std::vector<EntryType>& x,
                          int n_repetitions);
};

} // namespace step2 END

// In standard C++ it is a good habit to separate function declarations
// from their definitions. For template classes this does not work.
// Therefore, we put the definitions of the member functions into a
// separate header file to keep at least the improved overview provided
// by this separation.
#include <step-2/CpuBlas_driver_step-2.hh>
#endif // CPUBLAS_DRIVER_STEP2_H
