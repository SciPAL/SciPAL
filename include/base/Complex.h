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


#ifndef COMPLEX_H
#define COMPLEX_H
#include <base/ParallelArch.h>
#include <complex>
#include "CudaComplex.h"

template<typename T, ParallelArch arch>
struct ComplexTraits;

template<typename T>
struct ComplexTraits<T, cpu>
{
//    typedef typename std::complex<T> ComplexType;
    typedef typename step16::CudaComplex<T> ComplexType;
};


template<typename T>
struct ComplexTraits<T, gpu_cuda>
{
    typedef typename step16::CudaComplex<T> ComplexType;
};


template<typename T, ParallelArch arch>
struct Complex
{
    typedef typename ComplexTraits<T, arch>::ComplexType ComplexType;
};


#endif // COMPLEX_H
