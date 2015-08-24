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
//#include <lac/VectorCustomOperations.h>
#include <lac/scipal_kernels_wrapper.cu.h>

namespace SciPAL {


using namespace SciPAL;

template< typename T, typename BW>
struct BlasVecExp/*VectorExpressions*/
{
//        typedef typename ::SciPAL::Matrix<T, BW>::MyShape Mtx;
//        typedef typename ::SciPAL::Vector<T, BW>::MyShape Vtr;
//        typedef ::SciPAL::Shape<T, matrix> Mtx;
//        typedef ::SciPAL::Shape<T, vector> Vtr;

    typedef ::SciPAL::Matrix<T, BW> Mtx;
    typedef ::SciPAL::Vector<T, BW> Vtr;

    typedef ::SciPAL::SubMatrixView<T, BW> SMtx;
    typedef ::SciPAL::SubVectorView<T, BW, Vtr> SVtr;

    typedef Literal<T, BW> Lit;

    // $ y = \alpha x$
  //  typedef typename ::SciPAL::Mul<T, Vtr > scaledV;
    typedef typename  ::SciPAL::BinaryExpr<Lit, mult, Vtr> scaledV;

    // $ y = A \cdot x$
    typedef typename ::SciPAL::BinaryExpr<Mtx, mult, Vtr> // ::SciPAL::Mul<Mtx, Vtr >
    MV;

    // $ y = A \cdot x$
    typedef typename ::SciPAL::BinaryExpr<UnaryExpr<Mtx, expr_transpose >, mult, Vtr> // ::SciPAL::Mul<Mtx, Vtr >
    MtV;

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
    typedef typename SciPAL::BinaryExpr< UnaryExpr<Vtr, expr_transpose>, mult, Vtr>
    scalar_product;

    // vector = (smatrix) * vector
    typedef typename SciPAL::BinaryExpr<
    SubMatrixView<T, BW>, mult, Vtr > SMmV;

    // vector = (smatrix) * svector
    typedef typename SciPAL::BinaryExpr<
    SubMatrixView<T, BW> , mult, SVtr > SMmSV;

    // vector = transpose(smatrix) * vector
    typedef typename SciPAL::BinaryExpr<
    UnaryExpr<SciPAL::SubMatrixView<T, BW>, expr_transpose>, mult, Vtr > SMtmV;

    // vector = transpose(smatrix) * svector
    typedef typename SciPAL::BinaryExpr<
    UnaryExpr<SciPAL::SubMatrixView<T, BW>, expr_transpose>, mult, SubVectorView<T, BW, Vtr> >
    SMtmSV;




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
static void apply(Literal<T, BW> &result,
                  const typename BlasVecExp<T, BW>::scalar_product& expr)
{
    typedef Vector<T, BW> Vtr;

    const Vtr& x = expr.l.l;
    const Vtr& y = expr.r;

    result = BW::dot(x.size(), x.data(), 1, y.data(), 1);
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
    // BW::scal(n, alpha, &(C.data()[0]), incx);

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

    BW::axpy(n, alpha, &(A.data()[0]), incx,
            &(result.data()[0]), incy);
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
             A.data(),
             A.n_rows(), // lda
             x.data(),
             1, // incx
             beta,
             dst.data(), // y
             1 // incy
             );
}

// @sect4{Function: apply}
//!
//! Matrix-vector multiplication.
//! $dst = A^t \cdot src$

template <typename T, typename BW>
static void apply(Vector<T, BW> &result,
                  const typename BlasVecExp<T, BW>::MtV& expr)
{

    typedef ::SciPAL::Matrix<T, BW> Mtx;
    typedef Vector<T, BW> Vtr;

    T alpha = T(1);
    const Mtx & A = expr.l;
    const Vtr & x = expr.r;

    T beta = T(0);

    Vtr & dst = result;


    BW::gemv('t',
             A.n_rows(), // m
             A.n_cols(), // n
             alpha,
             A.data(),
             A.n_rows(), // lda
             x.data(),
             1, // incx
             beta,
             dst.data(), // y
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
             A.data(),
             A.n_rows(), // lda
             x.data(),
             1, // incx
             beta,
             dst.data(), // y
             1 // incy
             );
}


// @sect4{Function: apply}
//!
//! Matrix-vector multiplication.
//! $dst = A^t \cdot x$

template <typename T, typename BW>
static void apply(Vector<T, BW> &result,
                  const typename BlasVecExp<T, BW>::SMtmV& expr)
{

    typedef ::SciPAL::SubMatrixView<T, BW> Mtx;
    typedef Vector<T, BW> Vtr;

    T alpha = T(1);
    const Mtx & A = expr.l.l;
    const Vtr & x = expr.r;

    T beta = T(0);

    Vtr & dst = result;


    BW::gemv('t',
             A.n_rows(), // m
             A.n_cols(), // n
             alpha,
             A.data(),
             A.leading_dim(), // lda
             x.data(),
             1, // incx
             beta,
             dst.data(), // y
             1 // incy
             );
}

// @sect4{Function: apply}
//!
//! subMatrix^T-subvector multiplication.
//! $dst = A^t \cdot src$

template <typename T, typename BW>
static void apply(Vector<T, BW> &result,
                  const typename BlasVecExp<T, BW>::SMtmSV& expr)
{

    typedef ::SciPAL::SubMatrixView<T, BW> Mtx;
    typedef SubVectorView<T, BW, typename BlasVecExp<T, BW>::Vtr> SVtr;
    typedef typename BlasVecExp<T, BW>::Vtr Vtr;

    T alpha = T(1);
    const Mtx & A = expr.l.l;
    const SVtr & x = expr.r;

    T beta = T(0);

    Vtr & dst = result;


    BW::gemv('t',
             A.n_rows(), // m
             A.n_cols(), // n
             alpha,
             A.data(),
             A.leading_dim(), // lda
             x.data(),
             1, // incx
             beta,
             dst.data(), // y
             1 // incy
             );
}

// @sect4{Function: apply}
//!
//! subMatrix^T-subvector multiplication.
//! $dst = A^t \cdot src$

template <typename T, typename BW>
static void apply(SubVectorView<T, BW, typename BlasVecExp<T, BW>::Vtr> &result,
                  const typename BlasVecExp<T, BW>::SMtmSV& expr)
{

    typedef ::SciPAL::SubMatrixView<T, BW> Mtx;
    typedef SubVectorView<T, BW,typename BlasVecExp<T, BW>::Vtr> SVtr;
//    typedef typename BlasVecExp<T, BW>::Vtr Vtr;

    T alpha = T(1);
    const Mtx & A = expr.l.l;
    const SVtr & x = expr.r;

    T beta = T(0);

    SVtr & dst = result;


    BW::gemv('t',
             A.n_rows(), // m
             A.n_cols(), // n
             alpha,
             A.data(),
             A.leading_dim(), // lda
             x.data(),
             1, // incx
             beta,
             dst.data(), // y
             1 // incy
             );
}


// @sect4{Function: apply}
//!
//! Matrix-vector multiplication.
//! $dst = A \cdot src$

template <typename T, typename BW>
static void apply(Vector<T, BW> &result,
                  const typename BlasVecExp<T, BW>::SMmV& expr)
{

    typedef ::SciPAL::SubMatrixView<T, BW> Mtx;
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
             A.data(),
             A.leading_dim(), // lda
             x.data(),
             1, // incx
             beta,
             dst.data(), // y
             1 // incy
             );
}

// @sect4{Function: apply}
//!
//! Matrix-vector multiplication.
//! $dst = A \cdot x$

template <typename T, typename BW>
static void apply(SubVectorView<T, BW, typename BlasVecExp<T, BW>::Vtr> &result,
                  const typename BlasVecExp<T, BW>::SMmSV& expr)
{

    typedef ::SciPAL::SubMatrixView<T, BW> Mtx;
    typedef typename BlasVecExp<T, BW>::SVtr SVtr;

    T alpha = T(1);
    const Mtx & A = expr.l;
    const SVtr & x = expr.r;

    T beta = T(0);

    SVtr & dst = result;


    BW::gemv('n',
             A.n_rows(), // m
             A.n_cols(), // n
             alpha,
             A.data(),
             A.leading_dim(), // lda
             x.data(),
             1, // incx
             beta,
             dst.data(), // y
             1 // incy
             );
}


} // END namespace LAOOperations

} // END namespace SciPAL

#endif // VECTOROPERATIONS_H
