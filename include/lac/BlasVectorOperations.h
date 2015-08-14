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


#ifndef VECTOROPERATIONS_H
#define VECTOROPERATIONS_H

#include <base/ForewardDeclarations.h>

#include <lac/expression_templates_host.h>
#include <lac/VectorCustomOperations.h>
#include <lac/scipal_kernels_wrapper.cu.h>

namespace SciPAL {


using namespace SciPAL;

template< typename T, typename BW>
struct BlasVecExp/*VectorExpressions*/
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;
    //  typedef SparseMatrix<T, BW> SpMtx;
    typedef ::SciPAL::Vector<T, BW> Vtr;

    typedef Literal<T> Lit;

    // $ y = \alpha x$
  //  typedef typename ::SciPAL::Mul<T, Vtr > scaledV;
    typedef typename  ::SciPAL::BinaryExpr<Lit, mult, Vtr> scaledV;

    // $ y = A \cdot x$
    typedef typename ::SciPAL::BinaryExpr<Mtx, mult, Vtr> // ::SciPAL::Mul<Mtx, Vtr >
    MV;

    // $ y = \alpha A \cdot x$, just a crutch
    typedef typename ::SciPAL::BinaryExpr<Lit, mult, Mtx> scaledM;

    typedef typename ::SciPAL::BinaryExpr<scaledM, mult, Vtr> // SciPAL::Mul<Mul<T, Mtx>, Vtr >
    sMV;

    // $ y = \alpha A \cdot x + \beta y$
    typedef typename SciPAL::BinaryExpr<sMV, plus, scaledV>
    sMVaV;

    // $ y = \alpha x + y$
    typedef typename SciPAL::BinaryExpr<scaledV, plus, Vtr>
    axpy;

    // $ z =  x^t \cdot y$
    typedef typename SciPAL::BinaryExpr< SciPAL::transpose<Vtr>, mult, Vtr>
    scalar_product;
};


// @sect3{namespace: LAOOperations}
//
// The different apply() functions are collected in
// a namespace. Using a namespace and not a structure or class
// allows to extend the list of operations without
// having to modify any class of the library.
// A new operation is introduced by opening the namespace,
// defining the apply function and closing the namespace before
// the header file of the matrix or vector class is included.
// The corresponding apply is chosen by type.
namespace LAOOperations
{

// @sect4{Function: apply}
//!
//! apply function for literals scalar product

template <typename T, typename BW>
static void apply(Literal<T> &result,
                  const typename BlasVecExp<T, BW>::scalar_product& expr)
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;
    typedef Vector<T, BW> Vtr;

    const Vtr x = expr.l.A;
    const Vtr & y = expr.r;

    result = BW::dot(x.size(), x.array().val(), 1, y.array().val(), 1);

}

template <typename T, typename BW>
static void apply(Vector<T, BW> &result,
                  const typename BlasVecExp<T, BW>::scaledV& expr)
{
    typedef Vector<T, BW> Vtr;

    T alpha = expr.l;
    const Vtr & A = expr.r;

    Vtr & C = result;

    C = A;

    int incx = 1;

    //int n = A.n_elements(); this is not accessiable
    int n = A.size();
    // BW::scal(n, alpha, &(C.array().val()[0]), incx);

    typedef SciPAL::ShapeData<T> LAOShape;

      SciPAL::Kernels<T, BW::arch> dev_exp(4);
    DevBinaryExpr<DevLiteral<T>, mult, LAOShape> dev_e(expr);
       dev_exp.apply(result, dev_e);
}

template <typename T, typename BW>
static void apply(Vector<T, BW> &result,
                  const typename BlasVecExp<T, BW>::axpy& expr)
{
    typedef Vector<T, BW> Vtr;

    T alpha = expr.l.l;
    const Vtr & A = expr.l.r;

    result = expr.r;

    int incx = 1;
    int incy = 1;

    int n = A.size();

    BW::axpy(n, alpha, &(A.array().val()[0]), incx,
            &(result.array().val()[0]), incy);
}

// @sect4{Function: apply}
//!
//! Matrix-vector multiplication.
//! $dst = A \cdot src$

template <typename T, typename BW>
static void apply(Vector<T, BW> &result,
                  const typename BlasVecExp<T, BW>::MV& expr)
{

    typedef ::SciPAL::Matrix<T, BW> Mtx;
    typedef Vector<T, BW> Vtr;

    T alpha = T(1);
    const Mtx & A = expr.l;
    const Vtr & x = expr.r;

    T beta = T(0);

    Vtr & dst = result;


    BW::gemv('n',
             A.n_rows(), // m
             A.n_cols(), // n
             alpha,
             A.data_ptr,
             A.n_rows(), // lda
             x.data_ptr,
             1, // incx
             beta,
             dst.data_ptr, // y
             1 // incy
             );
}

// @sect4{Function: scaled_vmult}
//!
//! Generic matrix-vector multiplication.
//! $dst = \alpha A \cdot src + \beta dst$

template <typename T, typename BW>
static void apply(Vector<T, BW> &result,
                  const typename BlasVecExp<T, BW>::sMVaV& expr)
{

    typedef ::SciPAL::Matrix<T, BW> Mtx;
    typedef Vector<T, BW> Vtr;

    T alpha = expr.l.l.l;
    const Mtx & A = expr.l.l.r;
    const Vtr & x = expr.l.r;

    T beta = expr.r.l;

    Vtr & dst = result;


    BW::gemv('n',
             A.n_rows(), // m
             A.n_cols(), // n
             alpha,
             A.data_ptr,
             A.n_rows(), // lda
             x.data_ptr,
             1, // incx
             beta,
             dst.data_ptr, // y
             1 // incy
             );
}



} // END namespace LAOOperations

} // END namespace SciPAL

#endif // VECTOROPERATIONS_H
