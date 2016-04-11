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

#ifndef SCIPAL_KERNELS_WRAPPER_CU_H
#define SCIPAL_KERNELS_WRAPPER_CU_H
//OpenMP
#ifdef USE_OPENMP
#include <omp.h>
#endif
//SciPAL includes
#include <base/ForewardDeclarations.h>
#include <lac/expression_templates_device.h>
#include <lac/UnaryFunctions.h>
#include <base/PrecisionTraits.h>
#include <base/CudaComplex.h>
#include <lac/OperandInfo.h>



#include <iostream>

#include <base/ForewardDeclarations.h>
using namespace SciPAL;

template <typename T> class ImplCUDA {
public:

    typedef typename PrecisionTraits<T, gpu_cuda>::ComplexType NumberType2;
    typedef typename PrecisionTraits<T, gpu_cuda>::NumberType NumberType;

    template <typename X, typename op>
    static void apply(SciPAL::ShapeData<T> & d_dst,
                      const SciPAL::DevUnaryExpr<X, op> & Ax);

    template <typename L, typename op, typename R>
    static void apply(SciPAL::ShapeData<T>& d_dst,
                      const typename ::SciPAL::DevBinaryExpr<L, op, R> & Ax);

private:
};

template <typename T> class ImplOpenMP {
public:

    typedef typename PrecisionTraits<T, cpu>::ComplexType NumberType2;
    typedef typename PrecisionTraits<T, cpu>::NumberType  NumberType;

    template <typename X, typename op>
    static void apply(SciPAL::ShapeData<T> & d_dst,
                      const SciPAL::DevUnaryExpr<X, op> & Ax);

    template <typename L, typename op, typename R>
    static void apply(SciPAL::ShapeData<T> & d_dst,
                      const typename ::SciPAL::DevBinaryExpr<L, op, R> & Ax);





private:
};


namespace SciPAL {

// @sect3{Class: KernelsImpl}



// @sect3{Class: Kernels}
template <typename T, ParallelArch arch>
struct KernelBase;

template <typename T>
struct KernelBase<T, gpu_cuda>
{
    typedef  ImplCUDA<T> Type;
};

template <typename T>
struct KernelBase<T, cpu>
{
    typedef  ImplOpenMP<T> Type;
};

template<typename T, ParallelArch arch>
struct Kernels : KernelBase<T, arch>::Type
{
    Kernels(int /*num_omp_threads*/)
    {
        //omp_set_num_threads(num_omp_threads);
    }
};

} //end namespace SciPAL


#endif // SCIPAL_KERNELS_WRAPPER_CU_H
