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


#ifndef VTRAITS_H
#define VTRAITS_H
#include <base/ArchTraits.h>
#include <base/PrecisionTraits.h>
#include <lac/blas++.h>
#include <base/CudaComplex.h>

namespace SciPAL {

    // @sect3{struct: VTraits}
    // This struct holds some typedef for number and vector types which are used throughout
    //the program.
    template <typename T, ParallelArch arch>
    struct VTraits {

        typedef typename PrecisionTraits<T, arch>::NumberType NumberType;
        typedef typename PrecisionTraits<T, arch>::ComplexType ComplexType;
        typedef typename blas_pp<CudaComplex<NumberType>, typename archTraits<arch>::BlasType>::Vector cplxVector;
//        typedef typename blas_pp<ComplexType, typename archTraits<arch>::BlasType>::Vector cplxVector;
        typedef typename blas_pp<NumberType, typename archTraits<arch>::BlasType>::Vector Vector;
        typedef typename blas_pp<int, typename archTraits<arch>::BlasType>::Vector intVector;
        typedef typename blas_pp<uint, typename archTraits<arch>::BlasType>::Vector uintVector;

    };
}
#endif // VTRAITS_H
