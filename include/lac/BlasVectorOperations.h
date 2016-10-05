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
#include <deal.II/base/exceptions.h>

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
    typedef ::SciPAL::Vector<T, BW> Vtr;

    typedef ::SciPAL::SubMatrixView<T, BW> SMtx;
    typedef ::SciPAL::SubVectorView<T, BW, Vtr> SVtr;



    typedef ::SciPAL::SubVectorView<T, BW, const Mtx> SVtr_mtx;
    typedef ::SciPAL::ColVectorView<T, BW, const Mtx> CVtr;
    typedef ::SciPAL::RowVectorView<T, BW, const Mtx> RVtr;

    typedef ::SciPAL::SubVectorView<T, BW, const Vtr> cSVtr;

    typedef Literal<T, BW> Lit;

    // $ y = \alpha x$
  //  typedef typename ::SciPAL::Mul<T, Vtr > scaledV;
    typedef typename  ::SciPAL::BinaryExpr<Lit, mult, Vtr> scaledV;
    typedef typename  ::SciPAL::BinaryExpr<Lit, mult, SVtr_mtx> scaledSV;
    typedef typename  ::SciPAL::BinaryExpr<Lit, mult, RVtr> scaledRV;
    // $ y = A \cdot x$
    typedef typename ::SciPAL::BinaryExpr<Mtx, mult, Vtr> // ::SciPAL::Mul<Mtx, Vtr >
    MV;

    // $ y = A \cdot x$
    typedef typename ::SciPAL::BinaryExpr<UnaryExpr<Mtx, expr_transpose >, mult, Vtr> // ::SciPAL::Mul<Mtx, Vtr >
    MtV;

    // $ y = A \cdot x$
    typedef typename ::SciPAL::BinaryExpr<UnaryExpr<SMtx, expr_transpose >, mult, UnaryExpr<typename ::SciPAL::SubVectorView<T, BW, Mtx>, expr_transpose> > // ::SciPAL::Mul<Mtx, Vtr >
    SMtrVt_mtx;

    // $ y = \alpha A \cdot x$, just a crutch
    typedef typename ::SciPAL::BinaryExpr<Lit, mult, Mtx> scaledM;

    typedef typename ::SciPAL::BinaryExpr<scaledM, mult, Vtr> // SciPAL::Mul<Mul<T, Mtx>, Vtr >
    sMV;

    // $ y = \alpha A \cdot x + \beta y$
    typedef typename SciPAL::BinaryExpr<sMV, plus, scaledV>
    sMVasV;

    // $ y = \alpha x + y$
    typedef typename SciPAL::BinaryExpr<scaledV, plus, Vtr>
    axpy;

    // $ z =  x^t \cdot y$
    typedef typename SciPAL::BinaryExpr< UnaryExpr<Vtr, expr_transpose>, mult, Vtr>
    scalar_product;

    // $ z =  x^h \cdot y$
    typedef typename SciPAL::BinaryExpr< UnaryExpr<Vtr, expr_adjoint>, mult, Vtr>
    complex_scalar_product;

    typedef typename SciPAL::BinaryExpr< SVtr_mtx, mult, SVtr_mtx> scalar_product_on_view;

    // vector = (smatrix) * vector
    typedef typename SciPAL::BinaryExpr<
    SubMatrixView<T, BW>, mult, Vtr > SMmV;

    // vector = (smatrix) * svector
    typedef typename SciPAL::BinaryExpr<
    SubMatrixView<T, BW> , mult, SVtr > SMmSV;


    // vector = (smatrix) * svector
    typedef typename SciPAL::BinaryExpr<
    SubMatrixView<T, BW> , mult, SVtr_mtx > SMmSV_mtx;

    // vector = RowVector * SubMatrix
    typedef typename SciPAL::BinaryExpr<
    RVtr , mult, SMtx > RVmSM;

    // vector = RowVector * SubMatrix
    typedef typename SciPAL::BinaryExpr<
    SVtr_mtx , mult, SMtx > SV_mtxmSM;


    // vector = SubMatrix * ColVector
    typedef typename SciPAL::BinaryExpr<
    SMtx , mult, CVtr > SMmCV;

    // vector = transpose(smatrix) * vector
    typedef typename SciPAL::BinaryExpr<
    UnaryExpr<SciPAL::SubMatrixView<T, BW>, expr_transpose>, mult, Vtr > SMtmV;

    // vector = transpose(smatrix) * svector
    typedef typename SciPAL::BinaryExpr<
    UnaryExpr<SciPAL::SubMatrixView<T, BW>, expr_transpose>, mult, SubVectorView<T, BW, Vtr> > SMtmSV;

    // vector = transpose(smatrix) * svector_mtx
    typedef typename SciPAL::BinaryExpr<
    UnaryExpr<SciPAL::SubMatrixView<T, BW>, expr_transpose>, mult, SubVectorView<T, BW, const Mtx> > SMtmSV_mtx;

    // vector = smatrix * transpose(svector)
    typedef typename SciPAL::BinaryExpr<
    SciPAL::SubMatrixView<T, BW>, mult, UnaryExpr<SubVectorView<T, BW, const Mtx> , expr_transpose> > SMmSVt_mtx;

    // vector = matrix * svector_matrix
    typedef typename SciPAL::BinaryExpr<Matrix<T, BW> , mult, SVtr_mtx> MtxmSVtr_mtx;

    // svector_matrix = matrix * const_svector_vector
    typedef typename SciPAL::BinaryExpr<Matrix<T, BW> , mult, cSVtr> Mtx_m_cSVtr;


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
    Assert(x.size()==y.size(), dealii::ExcMessage("Input vector dimensions do not match in scalar_product."));

    result = BW::dot(x.size(), x.data(), x.stride, y.data(), y.stride);
}

template <typename T, typename BW>
static void apply(Literal<T, BW> &result,
                  const typename BlasVecExp<T, BW>::complex_scalar_product& expr)
{
    typedef Vector<T, BW> Vtr;

    const Vtr& x = expr.l.l;
    const Vtr& y = expr.r;
    Assert(x.size()==y.size(), dealii::ExcMessage("Input vector dimensions do not match in scalar_product."));

    result = BW::dotc(x.size(), x.data(), x.stride, y.data(), y.stride);
}

template <typename T, typename BW>
static void apply(Literal<T, BW> &result,
                  const typename BlasVecExp<T, BW>::scalar_product_on_view& expr)
{

    const typename BlasVecExp<T, BW>::SVtr_mtx& x = expr.l;
    const typename BlasVecExp<T, BW>::SVtr_mtx& y = expr.r;
    Assert(x.size()==y.size(), dealii::ExcMessage("Input vector dimensions do not match in scalar_product_on_view."));

    result = BW::dot(x.size(), x.data(), x.stride, y.data(), y.stride);
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
 
   Assert(A.size()==C.size(), dealii::ExcMessage("Input vector dimensions do not match in scaledV."));
    typedef SciPAL::ShapeData<T> LAOShape;

      SciPAL::Kernels<T, BW::arch> dev_exp(4);
    DevBinaryExpr<DevLiteral<T>, mult, LAOShape> dev_e(expr);
       dev_exp.apply(result, dev_e);
}

template <typename T, typename BW>
static void apply(Vector<T, BW> &result,
                  const typename BlasVecExp<T, BW>::scaledSV& expr)
{
    typedef Vector<T, BW> Vtr;
    typedef SciPAL::Matrix<T, BW> Mtx;
    typedef SciPAL::SubVectorView<T, BW, const Mtx> SVtr;
    T alpha = expr.l;
    const SVtr & A = expr.r;

    Vtr & C = result;

    C = A;

    int incx = C.stride();

    //int n = A.n_elements(); this is not accessiable
    int n = A.size();

    Assert(A.size()==C.size(), dealii::ExcMessage("Input vector dimensions do not match in scaledSV."));

    BW::scal(n, alpha, &(C.data()[0]), incx);

}

template <typename T, typename BW>
static void apply(SciPAL::RowVectorView<T, BW, const SciPAL::Matrix<T, BW> > &result,
                  const typename BlasVecExp<T, BW>::scaledRV& expr)
{
   // typedef Vector<T, BW> Vtr;
    typedef SciPAL::Matrix<T, BW> Mtx;
   // typedef SciPAL::SubVectorView<T, BW, const Mtx> SVtr;
    typedef SciPAL::RowVectorView<T, BW, const Mtx > RVtr;

    T alpha = expr.l;
    const RVtr & A = expr.r;
    Assert(A.size()==result.size(), dealii::ExcMessage("Input vector dimensions do not match in scaledSV."));
    result = A;
    BW::scal(A.size(), alpha, &(result.data()[0]),  result.stride);
}

template <typename T, typename BW>
static void apply(SciPAL::RowVectorView<T, BW, const SciPAL::Matrix<T, BW> > &result,
                  const typename BlasVecExp<T, BW>::scaledSV& expr)
{
   // typedef Vector<T, BW> Vtr;
    typedef SciPAL::Matrix<T, BW> Mtx;
    typedef SciPAL::SubVectorView<T, BW, const Mtx> SVtr;
    //typedef SciPAL::RowVectorView<T, BW, const Mtx > RVtr;

    T alpha = expr.l;
    const SVtr & A = expr.r;
    Assert(A.size()==result.size(), dealii::ExcMessage("Input vector dimensions do not match in scaledSV."));
    result = A;
    BW::scal(A.size(), alpha, &(result.data()[0]),  result.stride);
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
    Assert(A.size()==result.size(), dealii::ExcMessage("Input vector dimensions do not match in axpy."));

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
    dst.reinit(A.n_rows());
    Assert(A.n_cols()==x.n_rows(), dealii::ExcMessage("Input matrix/vector dimensions do not match in MV."));
    Assert(A.n_rows()==dst.n_rows(),
           dealii::ExcMessage("Result vector does not match input matrix dimensions in MV."));

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
//! $dst = A \cdot src$

template <typename T, typename BW>
static void apply(Vector<T, BW> &result,
                  const typename BlasVecExp<T, BW>::MtxmSVtr_mtx& expr)
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;
    typedef SubVectorView<T, BW, const Mtx> SVtr_mtx;
    typedef Vector<T, BW> Vtr;

    T alpha = T(1);
    const Mtx & A = expr.l;
    const SVtr_mtx & x = expr.r;

    T beta = T(0);

    Vtr & dst = result;
    dst.reinit(A.n_rows());
    Assert(A.n_cols()==x.n_elements_active, dealii::ExcMessage("Input matrix/vector dimensions do not match in MtxmSV_mtx."));
    Assert(A.n_rows()==dst.n_rows(),
           dealii::ExcMessage("Result vector does not match input matrix dimensions in MtxmSV_mtx."));
    BW::gemv('n',
             A.n_rows(), // m
             A.n_cols(), // n
             alpha,
             A.data(),
             A.n_rows(), // lda
             x.data(),
             x.stride, // incx
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
static void apply(SubVectorView<T, BW, Matrix<T, BW> > &result,
                  const typename BlasVecExp<T, BW>::Mtx_m_cSVtr& expr)
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;
    typedef SubVectorView<T, BW, Mtx> SVtr_mtx;
    typedef Vector<T, BW> Vtr;
    typedef SubVectorView<T, BW, const Vtr> sVtr_const_mtx;

    T alpha = T(1);
    const Mtx & A = expr.l;
    const sVtr_const_mtx & x = expr.r;

    T beta = T(0);

    SVtr_mtx & dst = result;

    Assert(dst.n_elements_active == A.n_rows(), dealii::ExcMessage("The size of the result view is wrong"));
    Assert(A.n_cols()==x.n_elements_active, dealii::ExcMessage("Input matrix/vector dimensions do not match in MtxmSV_mtx."));
    Assert(A.n_rows()==dst.n_rows_active(),
           dealii::ExcMessage("Result vector does not match input matrix dimensions in MtxmSV_mtx."));
    BW::gemv('n',
             A.n_rows(), // m
             A.n_cols(), // n
             alpha,
             A.data(),
             A.n_rows(), // lda
             x.data(),
             x.stride, // incx
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
    const Mtx & A = expr.l.l;
    const Vtr & x = expr.r;

    T beta = T(0);

    Vtr & dst = result;

    dst.reinit(A.n_cols());
    Assert(A.n_rows()==x.n_rows(), dealii::ExcMessage("Input matrix/vector dimensions do not match in MtV."));
    Assert(A.n_cols()==dst.n_rows(),
           dealii::ExcMessage("Result vector does not match input matrix dimensions in MtV."));

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

// @sect4{Function: apply}
//!
//! Matrix-vector multiplication.
//! $dst = A^t \cdot src$
template <typename T, typename BW>
static void apply(Vector<T, BW> &result,
                  const typename BlasVecExp<T, BW>::SMtrVt_mtx& expr)
{

    typedef ::SciPAL::Matrix<T, BW> Mtx;
    typedef typename BlasVecExp<T, BW>::SMtx SMtx;
    typedef SubVectorView<T, BW, Mtx> SVtr;
    typedef Vector<T, BW> Vtr;

    T alpha = T(1);
    const SMtx& A = expr.l.l;
    const SVtr & x = expr.r.l;

    T beta = T(0);

    Vtr & dst = result;

    dst.reinit(A.n_cols());
    Assert(A.n_rows()==x.n_cols(), dealii::ExcMessage("Input matrix/vector dimensions do not match in SMtrVt_mtx."));
    Assert(A.n_cols()==dst.n_rows(),
           dealii::ExcMessage("Result vector does not match input matrix dimensions in SMtrVt_mtx."));

    BW::gemv('t',
             A.n_rows(), // m
             A.n_cols(), // n
             alpha,
             A.data(),
             A.leading_dim(), // lda
             x.data(),
             x.stride, // incx
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
                  const typename BlasVecExp<T, BW>::sMVasV& expr)
{

    typedef ::SciPAL::Matrix<T, BW> Mtx;
    typedef Vector<T, BW> Vtr;

    T alpha = expr.l.l.l;
    const Mtx & A = expr.l.l.r;
    const Vtr & x = expr.l.r;

    T beta = expr.r.l;

    Vtr & dst = result;

    Assert(A.n_cols()==x.n_rows(), dealii::ExcMessage("Input matrix/vector dimensions do not match in MV."));
    Assert(A.n_rows()==dst.n_rows(),
           dealii::ExcMessage("Result vector does not match input matrix dimensions in MV."));

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

    dst.reinit(A.n_cols());
    Assert(A.n_rows()==x.n_rows(), dealii::ExcMessage("Input matrix/vector dimensions do not match in SMtmV."));
    Assert(A.n_cols()==dst.n_rows(),
           dealii::ExcMessage("Result vector does not match input matrix dimensions in SMtmV."));

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

    dst.reinit(A.n_cols());
    Assert(A.n_rows()==x.n_rows(), dealii::ExcMessage("Input matrix/vector dimensions do not match in SMtmSV."));
    Assert(A.n_cols()==dst.n_rows(),
           dealii::ExcMessage("Result vector does not match input matrix dimensions in SMtmSV."));

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

    Assert(A.n_rows()==x.n_rows(), dealii::ExcMessage("Input matrix/vector dimensions do not match in SMtmSV."));
    Assert(A.n_cols()==dst.n_rows(),
           dealii::ExcMessage("Result vector does not match input matrix dimensions in SMtmSV."));
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
             dst.stride // incy
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

    dst.reinit(A.n_rows());
    Assert(A.n_cols()==x.n_rows(), dealii::ExcMessage("Input matrix/vector dimensions do not match in SMmV."));
    Assert(A.n_rows()==dst.n_rows(),
           dealii::ExcMessage("Result vector does not match input matrix dimensions in SMmV."));

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
//! subMatrix-vector multiplication.
//! $dst = A \cdot src$

template <typename T, typename BW>
static void apply(SubVectorView<T, BW, typename BlasVecExp<T, BW>::Vtr> &result,
                  const typename BlasVecExp<T, BW>::SMmV& expr)
{

    typedef ::SciPAL::SubMatrixView<T, BW> Mtx;
    typedef Vector<T, BW> Vtr;

    T alpha = T(1);
    const Mtx & A = expr.l;
    const Vtr & x = expr.r;

    T beta = T(0);

    typename BlasVecExp<T, BW>::SVtr & dst = result;

    Assert(A.n_cols()==x.n_rows(), dealii::ExcMessage("Input matrix/vector dimensions do not match in SMmV."));
    Assert(A.n_rows()==dst.n_rows(),
           dealii::ExcMessage("Result vector does not match input matrix dimensions in SMmV."));

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
             dst.stride // incy
             );
}


//! same as above, but different result type
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

    Assert(A.n_cols()==x.n_rows(), dealii::ExcMessage("Input matrix/vector dimensions do not match in SMmSV."));
    Assert(A.n_rows()==dst.n_rows(),
           dealii::ExcMessage("Result vector does not match input matrix dimensions in SMmSV."));

    BW::gemv('n',
             A.n_rows(), // m
             A.n_cols(), // n
             alpha,
             A.data(),
             A.leading_dim(), // lda
             x.data(),
             x.stride, // incx
             beta,
             dst.data(), // y
             dst.stride // incy
             );
}

//! same as above, but different result type
template <typename T, typename BW>
static void apply(ColVectorView<T, BW, const SciPAL::Matrix<T, BW> > &result,
                  const typename BlasVecExp<T, BW>::SMmSV_mtx& expr)
{

    //typedef Vector<T, BW> Vtr;
    typedef SciPAL::Matrix<T, BW> Mtx;
    typedef SciPAL::SubMatrixView<T, BW> SMtx;
    typedef SciPAL::SubVectorView<T, BW, const Mtx> SVtr;
    typedef SciPAL::ColVectorView<T, BW, const Mtx > CVtr;

    T alpha = T(1);
    const SMtx & A = expr.l;
    const SVtr & x = expr.r;

    T beta = T(0);

    CVtr & dst = result;

    Assert(A.n_cols()==x.n_rows, dealii::ExcMessage("Input matrix/vector dimensions do not match in SMmSV."));
    Assert(A.n_rows()==dst.n_rows,
           dealii::ExcMessage("Result vector does not match input matrix dimensions in SMmSV."));

    BW::gemv('n',
             A.n_rows(), // m
             A.n_cols(), // n
             alpha,
             A.data(),
             A.leading_dim(), // lda
             x.data(),
             x.stride, // incx
             beta,
             dst.data(), // y
             dst.stride // incy
             );
}

//! same as above, but different result type
template <typename T, typename BW>
static void apply(Vector<T, BW>& result,
                  const typename BlasVecExp<T, BW>::SMmSV_mtx& expr)
{

    typedef ::SciPAL::SubMatrixView<T, BW> Mtx;
    typedef typename BlasVecExp<T, BW>::SVtr_mtx SVtr;
    typedef Vector<T, BW> Vtr;

    T alpha = T(1);
    const Mtx & A = expr.l;
    const SVtr & x = expr.r;

    T beta = T(0);

    Vtr & dst = result;

    dst.reinit(A.n_rows());
    Assert(A.n_cols()==x.size()/*n_rows()*/, dealii::ExcMessage("Input matrix/vector dimensions do not match in SMmSV_mtx."));
    Assert(A.n_rows()==dst.n_rows(),
           dealii::ExcMessage("Result vector does not match input matrix dimensions in SMmSV_mtx."));

    BW::gemv('n',
             A.n_rows(), // m
             A.n_cols(), // n
             alpha,
             A.data(),
             A.leading_dim(), // lda
             x.data(),
             x.stride, // incx
             beta,
             dst.data(), // y
             dst.stride // incy
             );
}

//! Row = Row*A
template <typename T, typename BW>
static void apply(typename SciPAL::RowVectorView<T, BW, const typename SciPAL::Matrix<T, BW> > &result,
                  const typename BlasVecExp<T, BW>::SV_mtxmSM& expr)
{

    // typedef Vector<T, BW> Vtr;
    typedef SciPAL::Matrix<T, BW> Mtx;
    typedef SciPAL::SubMatrixView<T, BW> SMtx;
    typedef SciPAL::SubVectorView<T, BW, const Mtx> SVtr;
    typedef SciPAL::RowVectorView<T, BW, const Mtx > RVtr;

    T alpha = T(1);
    const SVtr & y = expr.l;
    const SMtx & A = expr.r;

    T beta = T(0);

    RVtr & dst = result;

    Assert(y.n_cols==A.n_rows(), dealii::ExcMessage("Input matrix/vector dimensions do not match in RVmSM."));
    Assert(A.n_cols()==dst.n_cols,
           dealii::ExcMessage("Result vector does not match input matrix dimensions in RVmSM."));

    BW::gemv('t',
             A.n_rows(), // m
             A.n_cols(), // n
             alpha,
             A.data(),
             A.leading_dim(), // lda
             y.data(),
             y.stride, // incx
             beta,
             dst.data(), // y
             dst.leading_dim // incy
             );
}

// @sect4{Function: apply}
//!
//! subMatrix-subvector^T multiplication.
//! $dst = A^t \cdot src$

template <typename T, typename BW>
static void apply(Vector<T, BW> &result,
                  const typename BlasVecExp<T, BW>::SMmSVt_mtx& expr)
{

    typedef ::SciPAL::SubMatrixView<T, BW> SMtx;
    typedef SubVectorView<T, BW, const typename BlasVecExp<T, BW>::Mtx> SVtr;
    typedef typename BlasVecExp<T, BW>::Vtr Vtr;

    T alpha = T(1);
    const SMtx & A = expr.l;
    const SVtr & x = expr.r.l;

    T beta = T(0);

    Vtr & dst = result;

    dst.reinit(A.n_rows());
    Assert(A.n_cols()==x.size()/*n_cols()*/, dealii::ExcMessage("Input matrix/vector dimensions do not match in SMmSVt_mtx."));
    Assert(A.n_rows()==dst.n_rows(),
           dealii::ExcMessage("Result vector does not match input matrix dimensions in SMmSVt_mtx."));

    BW::gemv('n',
             A.n_rows(), // m
             A.n_cols(), // n
             alpha,
             A.data(),
             A.leading_dim(), // lda
             x.data(),
             x.stride, // incx
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
                  const typename BlasVecExp<T, BW>::SMtmSV_mtx& expr)
{

    typedef ::SciPAL::SubMatrixView<T, BW> SMtx;
    typedef SubVectorView<T, BW, const typename BlasVecExp<T, BW>::Mtx> SVtr;
    typedef typename BlasVecExp<T, BW>::Vtr Vtr;

    T alpha = T(1);
    const SMtx & A = expr.l.l;
    const SVtr & x = expr.r;

    T beta = T(0);

    Vtr & dst = result;

    dst.reinit(A.n_cols());
    Assert(A.n_rows()==x.size()/*n_rows()*/, dealii::ExcMessage("Input matrix/vector dimensions do not match in SMtmSV_mtx."));
    Assert(A.n_cols()==dst.size(),
           dealii::ExcMessage("Result vector does not match input matrix dimensions in SMtmSV_mtx."));

    BW::gemv('t',
             A.n_rows(), // m
             A.n_cols(), // n
             alpha,
             A.data(),
             A.leading_dim(), // lda
             x.data(),
             x.stride, // incx
             beta,
             dst.data(), // y
             1 // incy
             );
}


} // END namespace LAOOperations

} // END namespace SciPAL

#endif // VECTOROPERATIONS_H
