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


#ifndef CUDA_KERNEL_STEP_44_CU_H
#define CUDA_KERNEL_STEP_44_CU_H

#include <cuda.h>

namespace step44 {

    // @sect3{Class: Kernels}
    // Often the kernels need only one implementation yet it might be
    // interesting to study their behavior for different precisions.
    // This leads to the need to have a templatized implementation where
    // the precision or the particular number type is a template parameter.
    // This carries over to the class representing the interface to the
    // parallelized implementation. By adding a template parameter to this
    // class all its member functions are templatized as well.
    // Since we work with two compilers this requires us to provide full
    // specializations of this class at the end of the source file containing the kernels.
    // </br>
    // In large projects you will have several of these classes, probably in separate files.
template<typename T>
struct Kernels {

    void dummy();

    T a_dummy_attribute;
};

}

#endif // CUDA_KERNEL_STEP_44_CU_H
