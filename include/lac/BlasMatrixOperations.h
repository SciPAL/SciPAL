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


#ifndef MATRIXOPERATIONS_H
#define MATRIXOPERATIONS_H


#include <lac/expression_templates_host.h>

#include <base/ForewardDeclarations.h>

namespace SciPAL {

// @sect3{Struct: BlasMatExp}
//
// This structure summarizes all binary expressions
// which represent some operation, having a matrix as result, from the BLAS library.
template< typename T, typename BW> // IDEA:, MatrixStorage ms>
struct BlasMatExp/*MatrixExpressions*/
{
    typedef typename ::SciPAL::Matrix<T, BW>::MyShape Mtx;
    typedef typename ::SciPAL::Vector<T, BW>::MyShape Vtr;

//    typedef ::SciPAL::Matrix<T, BW> Mtx;
//    typedef ::SciPAL::Vector<T, BW> Vtr;
//    typedef ::SciPAL::Shape<T, matrix> Mtx;
//    typedef ::SciPAL::Shape<T, vector> Vtr;

//    typedef ::SciPAL::SubMatrixView<T, BW> SMtx;
//    typedef ::SciPAL::VectorView<T, BW> SVtr;

    typedef Literal<T> Lit;
    //scaled Matrix
    typedef typename SciPAL::BinaryExpr<Lit,  mult, Mtx> scaledM;

    //scaled matrix multiplied matrix
    typedef typename SciPAL::BinaryExpr<scaledM,  mult, Mtx> sMM;

    //matrix multiplied matrix
    typedef typename SciPAL::BinaryExpr<Mtx,  mult, Mtx> MM;

    // diagonal matrix multiplied matrix
    typedef typename SciPAL::BinaryExpr<SciPAL::diag<Vtr>,  mult, Mtx> dMM;

    // matrix multiplied diagonal matrix
    typedef typename SciPAL::BinaryExpr<Mtx,  mult, SciPAL::diag<Vtr> > MdM;


    typedef typename SciPAL::BinaryExpr< SciPAL::transpose<Mtx>, mult, Mtx> MtM;
    typedef typename SciPAL::BinaryExpr< SciPAL::adjoint<Mtx>, mult, Mtx> MaM;

    typedef typename SciPAL::BinaryExpr<Mtx, mult, SciPAL::transpose<Mtx> > MMt;
    typedef typename SciPAL::BinaryExpr<Mtx, mult, SciPAL::adjoint<Mtx> > MMa;

    //matrix times matrix plus matrix
    typedef typename SciPAL::BinaryExpr<MM, plus, Mtx>  MMaM;

    //scaled matrix times matrix plus matrix
    typedef typename SciPAL::BinaryExpr<sMM, plus, Mtx>  sMMaM;

    //generalized matrix matrix product
    typedef typename SciPAL::BinaryExpr<sMM, plus, scaledM> sMMasM;

    typedef typename SciPAL::adjoint<Mtx> adM;
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
// Similar to the pre-defined BLAS operations one can collect the specializations
// of BinaryExpr in a structure which has to be defined
// before the custom apply() function is defined.
namespace LAOOperations
{
template <typename T, typename BW> // IDEA:, MatrixStorage ms>
static void apply(::SciPAL::Matrix<T, BW> // IDEA:, ms>
                  &result,
                  const typename BlasMatExp<T, BW> // IDEA:, ms>
                  ::scaledM& expr)
{
    typedef ::SciPAL::Matrix<T, BW> // IDEA:, ms>
            Mtx;

    T alpha = expr.l;
    const Mtx & A = expr.r;

    Mtx & C = result;

    C = A; // FIXME: This could be optimized by swapping the data pointers.

    int incx = 1;

    int n = A.array().n_elements();

    BW::scal(n, alpha, &(C.array().val()[0]), incx);
}



template <typename T, typename BW> // IDEA:, MatrixStorage ms>
static void apply(::SciPAL::Matrix<T, BW> // IDEA:, ms>
                  &result,
                  const typename BlasMatExp<T, BW> // IDEA:, ms>
                  ::MM& expr)
{
    typedef ::SciPAL::Matrix<T, BW> // IDEA:, ms>
            Mtx;

    T alpha = 1.0;
    const Mtx & A = expr.l;
    const Mtx & B = expr.r;

    T beta = 0.0;

    Mtx & C = result;

    C.reinit(A.n_rows(), B.n_cols());

    A.scaled_mmult_add_scaled(C, B, 'n', 'n', alpha, beta);
}


template <typename T, typename BW>
static void apply(::SciPAL::Matrix<T, BW> &result,
                  const typename BlasMatExp<T, BW>::MtM& expr)
                  //SciPAL::BinaryExpr<SciPAL::transpose<SciPAL::Matrix<T, BW> >, SciPAL::mult, SciPAL::Matrix<T, BW> > & expr)

{
    typedef ::SciPAL::Matrix<T, BW> Mtx;

    T alpha = 1.0;
    const Mtx & A = expr.l.A;
    const Mtx & B = expr.r;

    T beta = 0.0;

    Mtx & C = result;

    C.reinit(A.n_cols(), B.n_cols());

    A.scaled_mmult_add_scaled(C, B, 't', 'n', alpha, beta);
}

template <typename T, typename BW>
static void apply(::SciPAL::Matrix<T, BW> &result,
                  const typename BlasMatExp<T, BW>::MaM& expr)
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;

    T alpha = 1.0;
    const Mtx & A = expr.l.A;
    const Mtx & B = expr.r;

    T beta = 0.0;

    Mtx & C = result;

    C.reinit(A.n_cols(), B.n_cols());

    A.scaled_mmult_add_scaled(C, B, 'c', 'n', alpha, beta);
}


template <typename T, typename BW>
static void apply(::SciPAL::Matrix<T, BW> &result,
                  const typename BlasMatExp<T, BW>::MMt& expr)
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;

    T alpha = 1.0;
    const Mtx & A = expr.l;
    const Mtx & B = expr.r.A;

    T beta = 0.0;

    Mtx & C = result;

    C.reinit(A.n_rows(), B.n_rows());

    A.scaled_mmult_add_scaled(C, B, 'n', 't', alpha, beta);
}


template <typename T, typename BW>
static void apply(::SciPAL::Matrix<T, BW> &result,
                  const typename BlasMatExp<T, BW>::MMa& expr)
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;

    T alpha = 1.0;
    const Mtx & A = expr.l;
    const Mtx & B = expr.r.A;

    T beta = 0.0;

    Mtx & C = result;

    C.reinit(A.n_rows(), B.n_rows());

    A.scaled_mmult_add_scaled(C, B, 'n', 'c', alpha, beta);
}


template <typename T, typename BW>
static void apply(::SciPAL::Matrix<T, BW> &result,
                  const typename BlasMatExp<T, BW>::dMM& expr)
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;
    typedef ::SciPAL::Vector<T, BW> Vtr;

    const Vtr & diag = expr.l.A;
    const Mtx & A = expr.r;


    Mtx & C = result;

    C.reinit(A.n_rows(), A.n_rows());

    BW::dgmm(CUBLAS_SIDE_LEFT,
             A.n_rows(), A.n_cols(),
             A.data_ptr, A.leading_dim,
             diag.data_ptr, 1,
             C.data_ptr, C.leading_dim);
}


template <typename T, typename BW>
static void apply(::SciPAL::Matrix<T, BW> &result,
                  const typename BlasMatExp<T, BW>::MdM& expr)
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;
    typedef ::SciPAL::Vector<T, BW> Vtr;

    const Vtr & diag = expr.r.A;
    const Mtx & A = expr.l;


    Mtx & C = result;

    C.reinit(A.n_rows(), A.n_rows());

    BW::dgmm(CUBLAS_SIDE_RIGHT,
             A.n_rows(), A.n_cols(),
             A.data_ptr, A.leading_dim,
             diag.data_ptr, 1,
             C.data_ptr, C.leading_dim);
}



template <typename T, typename BW>
static void apply(::SciPAL::Matrix<T, BW> &result,
                  const typename BlasMatExp<T, BW>::sMM& expr)
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;

    T alpha = expr.l.l;
    const Mtx & A = expr.l.r;
    const Mtx & B = expr.r;

    T beta = 0.0;

    Mtx & C = result;

    C.reinit(A.n_rows(), B.n_cols());

    A.scaled_mmult_add_scaled(C, B, 'n', 'n', alpha, beta);
}


template <typename T, typename BW>
static void apply(::SciPAL::Matrix<T, BW> & result,
                  const typename BlasMatExp<T, BW>::MMaM& expr)
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;

    T alpha = 1.0;
    const Mtx & A = expr.l.l;
    const Mtx & B = expr.l.r;

    T beta = 1.;

    Mtx & C = result; //this is equivalent to expr.r

    A.scaled_mmult_add_scaled(C, B, 'n', 'n', alpha, beta);
}


template <typename T, typename BW>
static void apply(::SciPAL::Matrix<T, BW> & result,
                  const typename BlasMatExp<T, BW>::sMMasM& expr)
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;

    T alpha = expr.l.l.l;
    const Mtx & A = expr.l.l.r;
    const Mtx & B = expr.l.r;

    T beta = expr.r.l;

    Mtx & C = result;

    A.scaled_mmult_add_scaled(C, B, 'n', 'n', alpha, beta);
}

template <typename T, typename BW>
static void apply(::SciPAL::Matrix<T, BW> & result,
                  const typename BlasMatExp<T, BW>::adM& expr)
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;

    T alpha = 1.0;
    const Mtx & A = expr.A;

    Mtx B(A.n_rows(), A.n_rows());
    for (unsigned int i = 0; i < B.n_rows(); i++)
        B(i, i, 1.0);

    T beta = 1.0;

    Mtx & C = result;

    A.scaled_mmult_add_scaled(C, B, 'c', 'n', alpha, beta);
}

} // END namespace LAOOperations

} // END namespace SciPAL


#endif // MATRIXOPERATIONS_H
