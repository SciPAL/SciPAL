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


#ifndef VECTORCUSTOMOPERATIONS_H
#define VECTORCUSTOMOPERATIONS_H

#include <lac/expression_templates_host.h>
#include <lac/expression_templates_device.h>

#include <lac/scipal_kernels_wrapper.cu.h>
#include <lac/UnaryFunctions.h>

#include <base/ForewardDeclarations.h>
#include <lac/expression_template_helper.h>


// TODO: rename this file to CustomOperations.h. Currently, "custom" means generic array arithmetic applying both to matrices and vectors.

namespace SciPAL {

// @sect3{namespace: VectorCustomOperations}
//
namespace LAOOperations
{


//! @sect3{function: apply}
//! This is a partial specialization of the function below for setting a matrix or vector to a certain
//! value given by the Literal.

//! this function compile sometimes to some strange behavior in the development
//! process thus this directive to disable it
#ifndef DEV_DEBUGGING
template <typename T, //! numbertype
          typename BW, //! blas type
          template <typename, typename> class LAO //! template for result type
          >
static void apply(LAO<T, BW> &result,
                  const ::SciPAL::Expr<Literal<T, BW> > &lit)
{
   ::SciPAL::Kernels<T, BW::arch> bla(4);
   UnaryExpr<Literal<T, BW>, Setter<T> >parent(~lit) ;//Expr envelope for Literal
   DevUnaryExpr<Literal<T, BW>, Setter<T> >child(parent);

   bla.apply(result, child);
}
#endif

//! @sect3{function: apply}
//! This function evaluates unary expressions applied to matrices or vectors, for instance $ A = \sin(A)$ where $A$ is a matrix or vector.
//! Evaluation is element-wise and thus perfectly parallel. Mainly for CUDA the function internally has to mirror
//! the given expression @p parent and rebuild it as device expression @p child which does not contain any
//! high-level types like Matrix or Vector. This is necessary for passing @p child by-value to the kernel which
//! actually evaluates the expression. @p child only contains low-level types like ShapeData which can be copied
//! by-value from the host to the device without caution.
//!
//! Since these are all templates, especially for the (CUDA) kernel a programmer is required to provide the
//! necessary template specializations in the part of the code which not compiled by the g++. The simplest solution is to
//! let the linker generate its "unresolved reference" errors which in fact display the exact form of the missing
//! specializations and to copy the error messages (after cleaning) into the source for the CUDA part of the program.
//! ==> This is actually done by the python script which is invoked inside the building process.
//!
//! @param result: linear algebra object which will contain the result of the expression.
//! @param parent: expresion object built on the host-side of the program.

//! this function compile sometimes to some strange behavior in the development
//! process thus this directive to disable it
#ifndef DEV_DEBUGGING
//! Specialization for views
template <typename T, //! numbertype
          typename BW, //! blas type
          typename T_src, //! source type of view
          template <typename, typename, typename> class LAO, //! template for result type
          typename E/*short for Expression*/ >
static void apply(LAO<T, BW, T_src> &result,
                  const ::SciPAL::Expr<E> &parent)
{
   SciPAL::Kernels<T, BW::arch> bla(4);
//   result.reinit((~parent).get_l());
    typename ExprChooser<E::I_am, E>::DevEType child(~parent);
   bla.apply(result, child);
}

template <typename T, //! numbertype
          typename BW, //! blas type
          template <typename, typename> class LAO, //! template for result type
          typename E/*short for Expression*/ >
static void apply(LAO<T, BW> &result,
                  const ::SciPAL::Expr<E> &parent)
{
   SciPAL::Kernels<T, BW::arch> bla(4);
//   result.reinit((~parent).get_l());
    typename ExprChooser<E::I_am, E>::DevEType child(~parent);
   bla.apply(result, child);
}

#endif

} // END namespace LAOOperations

} // END namespace SciPAL
#endif // VECTORCUSTOMOPERATIONS_H
