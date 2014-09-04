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


#ifndef PRECISIONTRAITS_H
#define PRECISIONTRAITS_H
#include <base/ParallelArch.h>
#include <complex>

typedef float2 cuComplex;
typedef double2 cuDoubleComplex;

namespace SciPAL{
    template <typename T> struct CudaComplex;
}

// @sect3{struct: PrecisionTraits}
template<typename T, ParallelArch arch> struct PrecisionTraits;

template <>
struct PrecisionTraits<cuComplex, gpu_cuda> {

    typedef float     NumberType;

    typedef cuComplex ComplexType;

};

template <>
struct PrecisionTraits<cuDoubleComplex, gpu_cuda> {

   typedef double       NumberType;

   typedef cuDoubleComplex ComplexType;

};

template <>
struct PrecisionTraits<float, gpu_cuda> {

    typedef float        NumberType;

    typedef cuComplex ComplexType;

};

template <>
struct PrecisionTraits<double, gpu_cuda> {

    typedef double          NumberType;

    typedef cuDoubleComplex ComplexType;

};


////////////////////////////////////////
//CPU specializations


template <>
struct PrecisionTraits<float2, cpu> {

    typedef float     NumberType;

    typedef float2    ComplexType;
};

template <>
struct PrecisionTraits<double2, cpu> {

   typedef double       NumberType;

   typedef double2 ComplexType;
};


template <>
struct PrecisionTraits<std::complex<float>, cpu> {

    typedef float     NumberType;

    typedef std::complex<float>    ComplexType;
};

template <>
struct PrecisionTraits<std::complex<double>, cpu> {

   typedef double       NumberType;

   typedef std::complex<double> ComplexType;
};


template <>
struct PrecisionTraits<float, cpu> {

    typedef float        NumberType;

    typedef float2 ComplexType;
};

template <>
struct PrecisionTraits<double, cpu> {

    typedef double          NumberType;

    typedef double2 ComplexType;
};

/////////////////////////////////////
// GPU / CPU specialization for int / uint

template <ParallelArch arch>
struct PrecisionTraits<int, arch> {

    typedef int         NumberType;

    typedef int2        ComplexType;

};

template <ParallelArch arch>
struct PrecisionTraits<int2, arch> {

    typedef int         NumberType;

    typedef int2        ComplexType;

};

template <ParallelArch arch>
struct PrecisionTraits<unsigned int, arch> {

    typedef unsigned int        NumberType;

    typedef uint2 ComplexType;

};

template <ParallelArch arch>
struct PrecisionTraits<uint2, arch> {

    typedef unsigned int    NumberType;

    typedef uint2           ComplexType;

};

//specialization for unsigned short
template <ParallelArch arch>
struct PrecisionTraits<unsigned short, arch> {

    typedef unsigned short        NumberType;

    typedef ushort2        ComplexType;

};

template <ParallelArch arch>
struct PrecisionTraits<ushort2, arch> {

    typedef unsigned short         NumberType;

    typedef ushort2        ComplexType;

};



//CudaComplex specialization

template <ParallelArch arch>
struct PrecisionTraits<SciPAL::CudaComplex<float>, arch> {

    typedef float  NumberType;

    typedef float2 ComplexType;

};

template < ParallelArch arch>
struct PrecisionTraits<SciPAL::CudaComplex<double>, arch> {

    typedef double  NumberType;

    typedef double2 ComplexType;

};


//pointer specialization
template <typename T, ParallelArch arch>
struct PrecisionTraits<T* , arch> {

    typedef typename PrecisionTraits<T, arch>::NumberType*  NumberType;

    typedef typename PrecisionTraits<T, arch>::ComplexType* ComplexType;
};



#endif // PRECISIONTRAITS_H
