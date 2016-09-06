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

#include <deal.II/base/exceptions.h>

#include <lac/expression_templates_host.h>
#include <lac/UnaryFunctions.h>
#include <base/ForewardDeclarations.h>
#include <lac/OperandInfo.h>

namespace SciPAL {

// @sect3{Struct: BlasMatExp}
//
// This structure summarizes all binary expressions
// which represent some operation, having a matrix as result, from the BLAS library.
template< typename T, typename BW> // IDEA:, MatrixStorage ms>
struct BlasMatExp/*MatrixExpressions*/
{
    //    typedef typename ::SciPAL::Matrix<T, BW>::MyShape Mtx;
    //    typedef typename ::SciPAL::Vector<T, BW>::MyShape Vtr;

    typedef ::SciPAL::Matrix<T, BW> Mtx;
    typedef ::SciPAL::Vector<T, BW> Vtr;
    //    typedef ::SciPAL::Shape<T, matrix> Mtx;
    //    typedef ::SciPAL::Shape<T, vector> Vtr;

    typedef ::SciPAL::SubMatrixView<T, BW> SMtx;
    typedef ::SciPAL::SubVectorView<T, BW, Vtr> SVtr;

    typedef Literal<T, BW> Lit;

//    template<EType ET, typename T2>
//    struct traverse_tree;

//    template< typename T2>
//    struct traverse_tree<leafE, T2>
//    {
//        typedef typename T2::Type Type;
//    };

//    template<EType ET, typename T2>
//    struct traverse_tree
//    {
//        typedef typename ExprChooser<ET, T2>::HostEType Type;
//    };

//    template <typename X>
//    struct generic_gemm
//    {
//        typedef typename
//        SciPAL::BinaryExpr<typename traverse_tree< X::L::I_am,
//                                         typename X::L::Type >::Type,
//        typename X::Operation,
//        typename traverse_tree< X::R::I_am,typename X::L>::Type> Type;
//    };

    //scaled Matrix
    typedef typename SciPAL::BinaryExpr<Lit,  mult, Mtx> scaledM;


    //scaled matrix multiplied matrix
    typedef typename SciPAL::BinaryExpr<scaledM,  mult, Mtx> sMM;

    //matrix multiplied matrix
    typedef typename SciPAL::BinaryExpr<Mtx,  mult, Mtx> MM;

    //submatrix multiplied matrix
    typedef typename SciPAL::BinaryExpr<SMtx,  mult, Mtx> SMM;

    //matrix multiplied submatrix
    typedef typename SciPAL::BinaryExpr<Mtx,  mult, SMtx> MSM;

    //submatrix multiplied submatrix
    typedef typename SciPAL::BinaryExpr<SMtx,  mult, SMtx> SMmSM;


    // diagonal matrix multiplied matrix
    typedef typename SciPAL::BinaryExpr<UnaryExpr<Vtr, expr_diag>,  mult, Mtx> dMM;

    // matrix multiplied diagonal matrix
    typedef typename SciPAL::BinaryExpr<Mtx,  mult, UnaryExpr<Vtr, expr_diag> > MdM;


    typedef typename SciPAL::BinaryExpr<UnaryExpr<Mtx, expr_transpose>, mult, Mtx> MtM;

    typedef typename SciPAL::BinaryExpr<Mtx, mult, UnaryExpr<Mtx, expr_transpose> > MMt;

    typedef typename SciPAL::BinaryExpr<UnaryExpr<Mtx, expr_transpose> , mult, UnaryExpr<Mtx, expr_transpose> > MtMt;

    typedef typename SciPAL::BinaryExpr<UnaryExpr<Mtx, expr_adjoint>, mult, Mtx> MaM;

    typedef typename SciPAL::BinaryExpr<Mtx, mult, UnaryExpr<Mtx, expr_adjoint> > MMa;

    typedef typename SciPAL::BinaryExpr<UnaryExpr<Mtx, expr_adjoint>, mult, UnaryExpr<Mtx, expr_adjoint> > MaMa;

    typedef typename SciPAL::BinaryExpr<UnaryExpr<Mtx, expr_transpose>, mult, UnaryExpr<Mtx, expr_adjoint> > MtMa;

    typedef typename SciPAL::BinaryExpr<UnaryExpr<Mtx, expr_adjoint>, mult, UnaryExpr<Mtx, expr_transpose> > MaMt;

    //matrix times matrix plus matrix
    typedef typename SciPAL::BinaryExpr<MM, plus, Mtx>  MMaM;

    //scaled matrix times matrix plus matrix
    typedef typename SciPAL::BinaryExpr<sMM, plus, Mtx>  sMMaM;

    //generalized matrix matrix product
    typedef typename SciPAL::BinaryExpr<sMM, plus, scaledM> sMMasM;

    typedef typename SciPAL::BinaryExpr<MtM, plus, Mtx> MtMpM;
    typedef typename SciPAL::BinaryExpr<MMt, plus, Mtx> MMtpM;

    typedef typename SciPAL::BinaryExpr<MaM, plus, Mtx> MaMpM;
    typedef typename SciPAL::BinaryExpr<MMa, plus, Mtx> MMapM;

    typedef typename SciPAL::UnaryExpr<Mtx, expr_adjoint> adM;

    typedef typename SciPAL::UnaryExpr<Mtx, expr_transpose> tM;


    /// (Should be a) Listing all kinds of 3 matrix multiplications (should be at least 432 (without submatrices...):
//! _ indicates bracketing

    /// MM * Something
    //matrix multiplied matrix multiplied matrix
    typedef typename SciPAL::BinaryExpr<MM,  mult, Mtx> MM_M;

    //matrix multiplied matrix multiplied transposed matrix
    typedef typename SciPAL::BinaryExpr<MM,  mult, tM> MM_Mt;

    //matrix multiplied matrix multiplied adjoint matrix
    typedef typename SciPAL::BinaryExpr<MM,  mult, adM> MM_Ma;

    /// MMt * Something
    //matrix multiplied transposed matrix multiplied matrix
    typedef typename SciPAL::BinaryExpr<MMt,  mult, Mtx> MMt_M;

    //matrix multiplied transposed matrix multiplied transposed matrix
    typedef typename SciPAL::BinaryExpr<MMt,  mult, tM> MMt_Mt;

    //matrix multiplied transposed matrix multiplied adjoint matrix
    typedef typename SciPAL::BinaryExpr<MMt,  mult, adM> MMt_Ma;

    /// MMa * Something
    //matrix multiplied adjoint matrix multiplied matrix
    typedef typename SciPAL::BinaryExpr<MMa,  mult, Mtx> MMa_M;

    //matrix multiplied adjoint matrix multiplied transposed matrix
    typedef typename SciPAL::BinaryExpr<MMa,  mult, tM> MMa_Mt;

    //matrix multiplied adjoint matrix multiplied adjoint matrix
    typedef typename SciPAL::BinaryExpr<MMa,  mult, adM> MMa_Ma;


    /// MtM * Something
    //transposed matrix multiplied matrix multiplied matrix
    typedef typename SciPAL::BinaryExpr<MtM,  mult, Mtx> MtM_M;

    //transposed matrix multiplied matrix multiplied transposed matrix
    typedef typename SciPAL::BinaryExpr<MtM,  mult, tM> MtM_Mt;

    //transposed matrix multiplied matrix multiplied adjoint matrix
    typedef typename SciPAL::BinaryExpr<MtM,  mult, adM> MtM_Ma;

    ///MtMt * Something
    //transposed matrix multiplied transposed matrix multiplied matrix
    typedef typename SciPAL::BinaryExpr<MtMt,  mult, Mtx> MtMt_M;

    //transposed matrix multiplied transposed matrix multiplied transposed matrix
    typedef typename SciPAL::BinaryExpr<MtMt,  mult, tM> MtMt_Mt;

    //transposed matrix multiplied transposed matrix multiplied adjoint matrix
    typedef typename SciPAL::BinaryExpr<MtMt,  mult, adM> MtMt_Ma;

    /// MtMa * Something
    //transposed matrix multiplied adjoint matrix multiplied matrix
    typedef typename SciPAL::BinaryExpr<MtMa,  mult, Mtx> MtMa_M;

    //transposed matrix multiplied adjoint matrix multiplied transposed matrix
    typedef typename SciPAL::BinaryExpr<MtMa,  mult, tM> MtMa_Mt;

    //transposed matrix multiplied adjoint matrix multiplied adjoint matrix
    typedef typename SciPAL::BinaryExpr<MtMa,  mult, adM> MtMa_Ma;


    ///MaM * Something
    //matrix multiplied matrix multiplied matrix
    typedef typename SciPAL::BinaryExpr<MaM,  mult, Mtx> MaM_M;

    //matrix multiplied matrix multiplied transposed matrix
    typedef typename SciPAL::BinaryExpr<MaM,  mult, tM> MaM_Mt;

    //matrix multiplied matrix multiplied adjoint matrix
    typedef typename SciPAL::BinaryExpr<MaM,  mult, adM> MaM_Ma;

    ///MaMt * Something
    //matrix multiplied transposed matrix multiplied matrix
    typedef typename SciPAL::BinaryExpr<MaMt,  mult, Mtx> MaMt_M;

    //matrix multiplied transposed matrix multiplied transposed matrix
    typedef typename SciPAL::BinaryExpr<MaMt,  mult, tM> MaMt_Mt;

    //matrix multiplied transposed matrix multiplied adjoint matrix
    typedef typename SciPAL::BinaryExpr<MaMt,  mult, adM> MaMt_Ma;

    ///MaMa * Something
    //matrix multiplied adjoint matrix multiplied matrix
    typedef typename SciPAL::BinaryExpr<MaMa,  mult, Mtx> MaMa_M;

    //matrix multiplied adjoint matrix multiplied transposed matrix
    typedef typename SciPAL::BinaryExpr<MaMa,  mult, tM> MaMa_Mt;

    //matrix multiplied adjoint matrix multiplied adjoint matrix
    typedef typename SciPAL::BinaryExpr<MaMa,  mult, adM> MaMa_Ma;



    /// Mtx * Something
    //matrix multiplied matrix multiplied matrix
    typedef typename SciPAL::BinaryExpr<Mtx,  mult, MM> M_MM;

    //matrix multiplied matrix multiplied transposed matrix
    typedef typename SciPAL::BinaryExpr<Mtx,  mult, MMt> M_MMt;

    //matrix multiplied matrix multiplied adjoint matrix
    typedef typename SciPAL::BinaryExpr<Mtx,  mult, MMa> M_MMa;

    //matrix multiplied transposed matrix multiplied matrix
    typedef typename SciPAL::BinaryExpr<Mtx,  mult, MtM> M_MtM;

    //matrix multiplied transposed matrix multiplied transposed matrix
    typedef typename SciPAL::BinaryExpr<Mtx,  mult, MtMt> M_MtMt;

    //matrix multiplied transposed matrix multiplied adjoint matrix
    typedef typename SciPAL::BinaryExpr<Mtx,  mult, MtMa> M_MtMa;

    //matrix multiplied adjoint matrix multiplied matrix
    typedef typename SciPAL::BinaryExpr<Mtx,  mult, MaM> M_MaM;

    //matrix multiplied adjoint matrix multiplied transposed matrix
    typedef typename SciPAL::BinaryExpr<Mtx,  mult, MaMt> M_MaMt;

    //matrix multiplied adjoint matrix multiplied adjoint matrix
    typedef typename SciPAL::BinaryExpr<Mtx,  mult, MaMa> M_MaMa;


    /// tM * Something
    //transposed matrix multiplied matrix multiplied matrix
    typedef typename SciPAL::BinaryExpr<tM,  mult, MM> Mt_MM;

    //transposed matrix multiplied matrix multiplied transposed matrix
    typedef typename SciPAL::BinaryExpr<tM,  mult, MMt> Mt_MMt;

    //transposed matrix multiplied matrix multiplied adjoint matrix
    typedef typename SciPAL::BinaryExpr<tM,  mult, MMa> Mt_MMa;

    //transposed matrix multiplied transposed matrix multiplied matrix
    typedef typename SciPAL::BinaryExpr<tM,  mult, MtM> Mt_MtM;

    //transposed matrix multiplied transposed matrix multiplied transposed matrix
    typedef typename SciPAL::BinaryExpr<tM,  mult, MtMt> Mt_MtMt;

    //transposed matrix multiplied transposed matrix multiplied adjoint matrix
    typedef typename SciPAL::BinaryExpr<tM,  mult, MtMa> Mt_MtMa;

    //transposed matrix multiplied adjoint matrix multiplied matrix
    typedef typename SciPAL::BinaryExpr<tM,  mult, MaM> Mt_MaM;

    //transposed matrix multiplied adjoint matrix multiplied transposed matrix
    typedef typename SciPAL::BinaryExpr<tM,  mult, MaMt> Mt_MaMt;

    //transposed matrix multiplied adjoint matrix multiplied adjoint matrix
    typedef typename SciPAL::BinaryExpr<tM,  mult, MaMa> Mt_MaMa;


    /// adM * Something
    //matrix multiplied matrix multiplied matrix
    typedef typename SciPAL::BinaryExpr<adM,  mult, MM> Ma_MM;

    //matrix multiplied matrix multiplied transposed matrix
    typedef typename SciPAL::BinaryExpr<adM,  mult, MMt> Ma_MMt;

    //matrix multiplied matrix multiplied adjoint matrix
    typedef typename SciPAL::BinaryExpr<adM,  mult, MMa> Ma_MMa;

    //matrix multiplied transposed matrix multiplied matrix
    typedef typename SciPAL::BinaryExpr<adM,  mult, MtM> Ma_MtM;

    //matrix multiplied transposed matrix multiplied transposed matrix
    typedef typename SciPAL::BinaryExpr<adM,  mult, MtMt> Ma_MtMt;

    //matrix multiplied transposed matrix multiplied adjoint matrix
    typedef typename SciPAL::BinaryExpr<adM,  mult, MtMa> Ma_MtMa;

    //matrix multiplied adjoint matrix multiplied matrix
    typedef typename SciPAL::BinaryExpr<adM,  mult, MaM> Ma_MaM;

    //matrix multiplied adjoint matrix multiplied transposed matrix
    typedef typename SciPAL::BinaryExpr<adM,  mult, MaMt> Ma_MaMt;

    //matrix multiplied adjoint matrix multiplied adjoint matrix
    typedef typename SciPAL::BinaryExpr<adM,  mult, MaMa> Ma_MaMa;

    /// TODO: ADD submatrixviews and scaled matrices...

    //matrix multiplied matrix multiplied matrix
    typedef typename SciPAL::BinaryExpr<sMM,  mult, Mtx> sMM_M;

    typedef typename SciPAL::BinaryExpr<Lit, mult, MM_M> s_MM_M;

    //matrix multiplied matrix multiplied matrix
    typedef typename SciPAL::BinaryExpr<Lit,  mult, M_MM> s_M_MM;

    //matrix multiplied matrix multiplied matrix
    typedef typename SciPAL::BinaryExpr<scaledM,  mult, MM> sM_MM;

    typedef typename SciPAL::BinaryExpr<SciPAL::BinaryExpr<scaledM,  mult, adM>,  mult, Mtx> sMaM_M;
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
//template <typename X, typename T, typename BW,
//          template<typename, typename> class LAO>
//static void apply(LAO<T, BW> & result,
//                  const typename BlasMatExp<T, BW>::template generic_gemm<X>::Type& expr)
//{
//    typedef ::SciPAL::Matrix<T, BW> Mtx;

//    T alpha = expr.l.l.l;
//    const Mtx & A = expr.l.l.r;
//    const Mtx & B = expr.l.r;

//    T beta = expr.r.l;

//    Mtx & C = result;

//    A.scaled_mmult_add_scaled(C, B, 'n', 'n', alpha, beta);
//}



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

    int n = A.size();
    Assert(A.n_rows()==C.n_rows() && A.n_cols()==C.n_cols(),
           dealii::ExcMessage("Result matrix does not match input matrix dimensions in scaledM."));
    BW::scal(n, alpha, &(C.data()[0]), incx);
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
    // A * B
    T beta = 0.0;

    Mtx & C = result;

    C.reinit(A.n_rows(), B.n_cols());

    Assert(A.n_cols()==B.n_rows(), dealii::ExcMessage("Input matrix dimensions do not match in MM."));
    Assert(A.n_rows()==C.n_rows() && B.n_cols()==C.n_cols(),
           dealii::ExcMessage("Result matrix does not match input matrix dimensions in MM."));
    A.scaled_mmult_add_scaled(C, B, 'n', 'n', alpha, beta);
}


template <typename T, typename BW>
static void apply(::SciPAL::SubMatrixView<T, BW>
                  &result,
                  const typename BlasMatExp<T, BW>
                  ::SMmSM& expr)
{
    typedef typename BlasMatExp<T, BW>::SMtx SMtx;
    const SMtx& A = expr.l;
    const SMtx& B = expr.r;

    int lda = A.leading_dim();
    int ldb = B.leading_dim(); //! src == B
    int ldc = result.leading_dim(); //! dst == C
    T alpha = (1.0);
    T beta = T(0.0);

    Assert(A.n_cols()==B.n_rows(), dealii::ExcMessage("Input matrix dimensions do not match in SMmSM."));
    Assert(A.n_rows()==result.n_rows() && B.n_cols()==result.n_cols(),
           dealii::ExcMessage("Result matrix does not match input matrix dimensions in SMmSM."));
    BW::gemm('n',
             'n',
             A.n_rows(), /*test for transpose...*/
             /* cublas doc : m == n_rows of op(A), i.e. n_cols for A^T*/
             B.n_cols(),
             /* cublas doc : n == n_cols of op(B), i.e. n_cols of C */
             A.n_cols(),
             /* cublas doc : k == n_cols of op(A), i.e. n_rows of op(B) or n_rows for A^T */
             alpha,
             A.data(), lda,
             B.data(), ldb,
             beta,
             result.data(), ldc);

}

template <typename T, typename BW>
static void apply(::SciPAL::Matrix<T, BW>
                  &result,
                  const typename BlasMatExp<T, BW>
                  ::MSM& expr)
{
    typedef typename BlasMatExp<T, BW>::SMtx SMtx;
    typedef ::SciPAL::Matrix<T, BW> Mtx;

    const Mtx& A = expr.l;
    const SMtx& B = expr.r;

    Mtx & C = result;

    C.reinit(A.n_rows(), B.n_cols());

    int lda = A.leading_dim;
    int ldb = B.leading_dim(); //! src == B
    int ldc = result.leading_dim; //! dst == C
    T alpha = (1.0);
    T beta = T(0.0);
    Assert(A.n_cols()==B.n_rows(), dealii::ExcMessage("Input matrix dimensions do not match in MSM."));
    Assert(A.n_rows()==C.n_rows() && B.n_cols()==C.n_cols(),
           dealii::ExcMessage("Result matrix does not match input matrix dimensions in MSM."));
    BW::gemm('n',
             'n',
             A.n_rows(), /*test for transpose...*/
             /* cublas doc : m == n_rows of op(A), i.e. n_cols for A^T*/
             B.n_cols(),
             /* cublas doc : n == n_cols of op(B), i.e. n_cols of C */
             A.n_cols(),
             /* cublas doc : k == n_cols of op(A), i.e. n_rows of op(B) or n_rows for A^T */
             alpha,
             A.data(), lda,
             B.data(), ldb,
             beta,
             result.data(), ldc);

}

template <typename T, typename BW>
static void apply(::SciPAL::SubMatrixView<T, BW>
                  &result,
                  const typename BlasMatExp<T, BW>
                  ::MSM& expr)
{
    typedef typename BlasMatExp<T, BW>::SMtx SMtx;
    typedef ::SciPAL::Matrix<T, BW> Mtx;

    const Mtx& A = expr.l;
    const SMtx& B = expr.r;

    SMtx & C = result;

    int lda = A.leading_dim;
    int ldb = B.leading_dim(); //! src == B
    int ldc = result.leading_dim(); //! dst == C
    T alpha = (1.0);
    T beta = T(0.0);
    Assert(A.n_cols()==B.n_rows(), dealii::ExcMessage("Input matrix dimensions do not match in MSM."));
    Assert(A.n_rows()==C.n_rows() && B.n_cols()==C.n_cols(),
           dealii::ExcMessage("Result matrix does not match input matrix dimensions in MSM."));
    BW::gemm('n',
             'n',
             A.n_rows(), /*test for transpose...*/
             /* cublas doc : m == n_rows of op(A), i.e. n_cols for A^T*/
             B.n_cols(),
             /* cublas doc : n == n_cols of op(B), i.e. n_cols of C */
             A.n_cols(),
             /* cublas doc : k == n_cols of op(A), i.e. n_rows of op(B) or n_rows for A^T */
             alpha,
             A.data(), lda,
             B.data(), ldb,
             beta,
             result.data(), ldc);

}

template <typename T, typename BW>
static void apply(::SciPAL::Matrix<T, BW>
                  &result,
                  const typename BlasMatExp<T, BW>
                  ::SMM& expr)
{
    typedef typename BlasMatExp<T, BW>::SMtx SMtx;
    typedef ::SciPAL::Matrix<T, BW> Mtx;

    const SMtx& A = expr.l;
    const Mtx& B = expr.r;

    Mtx & C = result;

    C.reinit(A.n_rows(), B.n_cols());

    int lda = A.leading_dim();
    int ldb = B.leading_dim; //! src == B
    int ldc = result.leading_dim; //! dst == C
    T alpha = (1.0);
    T beta = T(0.0);
    Assert(A.n_cols()==B.n_rows(), dealii::ExcMessage("Input matrix dimensions do not match in SMM."));
    Assert(A.n_rows()==C.n_rows() && B.n_cols()==C.n_cols(),
           dealii::ExcMessage("Result matrix does not match input matrix dimensions in SMM."));
    BW::gemm('n',
             'n',
             A.n_rows(), /*test for transpose...*/
             /* cublas doc : m == n_rows of op(A), i.e. n_cols for A^T*/
             B.n_cols(),
             /* cublas doc : n == n_cols of op(B), i.e. n_cols of C */
             A.n_cols(),
             /* cublas doc : k == n_cols of op(A), i.e. n_rows of op(B) or n_rows for A^T */
             alpha,
             A.data(), lda,
             B.data(), ldb,
             beta,
             result.data(), ldc);

}

template <typename T, typename BW>
static void apply(::SciPAL::Matrix<T, BW>
                  &result,
                  const typename BlasMatExp<T, BW>
                  ::SMmSM& expr)
{
    typedef typename BlasMatExp<T, BW>::SMtx SMtx;
    typedef ::SciPAL::Matrix<T, BW> Mtx;

    const SMtx& A = expr.l;
    const SMtx& B = expr.r;

    Mtx & C = result;

    C.reinit(A.n_rows(), B.n_cols());

    int lda = A.leading_dim();
    int ldb = B.leading_dim(); //! src == B
    int ldc = result.leading_dim; //! dst == C
    T alpha = (1.0);
    T beta = T(0.0);
    Assert(A.n_cols()==B.n_rows(), dealii::ExcMessage("Input matrix dimensions do not match in SMmSM."));
    Assert(A.n_rows()==C.n_rows() && B.n_cols()==C.n_cols(),
           dealii::ExcMessage("Result matrix does not match input matrix dimensions in SMmSM."));

    BW::gemm('n',
             'n',
             A.n_rows(), /*test for transpose...*/
             /* cublas doc : m == n_rows of op(A), i.e. n_cols for A^T*/
             B.n_cols(),
             /* cublas doc : n == n_cols of op(B), i.e. n_cols of C */
             A.n_cols(),
             /* cublas doc : k == n_cols of op(A), i.e. n_rows of op(B) or n_rows for A^T */
             alpha,
             A.data(), lda,
             B.data(), ldb,
             beta,
             result.data(), ldc);

}

template <typename T, typename BW>
static void apply(::SciPAL::Matrix<T, BW> &result,
                  const typename BlasMatExp<T, BW>::MtM& expr)

{
    typedef ::SciPAL::Matrix<T, BW> Mtx;

    T alpha = 1.0;
    const Mtx & A = expr.l.l;
    const Mtx & B = expr.r;

    T beta = 0.0;

    Mtx & C = result;

    C.reinit(A.n_cols(), B.n_cols());

    Assert(A.n_rows()==B.n_rows(), dealii::ExcMessage("Input matrix dimensions do not match in MtM."));
    Assert(A.n_cols()==C.n_rows() && B.n_cols()==C.n_cols(),
           dealii::ExcMessage("Result matrix does not match input matrix dimensions in MtM."));
    A.scaled_mmult_add_scaled(C, B, 't', 'n', alpha, beta);
}

template <typename T, typename BW>
static void apply(::SciPAL::Matrix<T, BW> &result,
                  const typename BlasMatExp<T, BW>::MaM& expr)
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;

    T alpha = 1.0;
    const Mtx & A = expr.l.l;
    const Mtx & B = expr.r;

    T beta = 0.0;

    Mtx & C = result;

    C.reinit(A.n_cols(), B.n_cols());

    Assert(A.n_rows()==B.n_rows(), dealii::ExcMessage("Input matrix dimensions do not match in MaM."));
    Assert(A.n_cols()==C.n_rows() && B.n_cols()==C.n_cols(),
           dealii::ExcMessage("Result matrix does not match input matrix dimensions in MaM."));

    A.scaled_mmult_add_scaled(C, B, 'c', 'n', alpha, beta);
}


template <typename T, typename BW>
static void apply(::SciPAL::Matrix<T, BW> &result,
                  const typename BlasMatExp<T, BW>::MMt& expr)
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;

    T alpha = 1.0;
    const Mtx & A = expr.l;
    const Mtx & B = expr.r.l;

    T beta = 0.0;

    Mtx & C = result;

    C.reinit(A.n_rows(), B.n_rows());

    Assert(A.n_cols()==B.n_cols(), dealii::ExcMessage("Input matrix dimensions do not match in MMt."));
    Assert(A.n_rows()==C.n_rows() && B.n_rows()==C.n_cols(),
           dealii::ExcMessage("Result matrix does not match input matrix dimensions in MMt."));
    A.scaled_mmult_add_scaled(C, B, 'n', 't', alpha, beta);
}

template <typename T, typename BW>
static void apply(::SciPAL::Matrix<T, BW> &result,
                  const typename BlasMatExp<T, BW>::MtMt& expr)
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;

    T alpha = 1.0;
    const Mtx & A = expr.l.l;
    const Mtx & B = expr.r.l;

    T beta = 0.0;

    Mtx & C = result;

    C.reinit(A.n_cols(), B.n_rows());
    Assert(A.n_rows()==B.n_cols(), dealii::ExcMessage("Input matrix dimensions do not match in MtMt."));
    Assert(A.n_cols()==C.n_rows() && B.n_rows()==C.n_cols(),
           dealii::ExcMessage("Result matrix does not match input matrix dimensions in MtMt."));

    A.scaled_mmult_add_scaled(C, B, 't', 't', alpha, beta);
}

template <typename T, typename BW>
static void apply(::SciPAL::Matrix<T, BW> &result,
                  const typename BlasMatExp<T, BW>::MMa& expr)
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;

    T alpha = 1.0;
    const Mtx & A = expr.l;
    const Mtx & B = expr.r.l;

    T beta = 0.0;

    Mtx & C = result;

    C.reinit(A.n_rows(), B.n_rows());
    Assert(A.n_cols()==B.n_cols(), dealii::ExcMessage("Input matrix dimensions do not match in MMa."));
    Assert(A.n_rows()==C.n_rows() && B.n_rows()==C.n_cols(),
           dealii::ExcMessage("Result matrix does not match input matrix dimensions in MMa."));

    A.scaled_mmult_add_scaled(C, B, 'n', 'c', alpha, beta);
}

template <typename T, typename BW>
static void apply(::SciPAL::Matrix<T, BW> &result,
                  const typename BlasMatExp<T, BW>::MaMa& expr)
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;

    T alpha = 1.0;
    const Mtx & A = expr.l.l;
    const Mtx & B = expr.r.l;

    T beta = 0.0;

    Mtx & C = result;

    C.reinit(A.n_cols(), B.n_rows());
    Assert(A.n_rows()==B.n_cols(), dealii::ExcMessage("Input matrix dimensions do not match in MaMa."));
    Assert(A.n_cols()==C.n_rows() && B.n_rows()==C.n_cols(),
           dealii::ExcMessage("Result matrix does not match input matrix dimensions in MaMa."));

    A.scaled_mmult_add_scaled(C, B, 'c', 'c', alpha, beta);
}

template <typename T, typename BW>
static void apply(::SciPAL::Matrix<T, BW> &result,
                  const typename BlasMatExp<T, BW>::MaMt& expr)
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;

    T alpha = 1.0;
    const Mtx & A = expr.l.l;
    const Mtx & B = expr.r.l;

    T beta = 0.0;

    Mtx & C = result;

    C.reinit(A.n_cols(), B.n_rows());
    Assert(A.n_rows()==B.n_cols(), dealii::ExcMessage("Input matrix dimensions do not match in MaMt."));
    Assert(A.n_cols()==C.n_rows() && B.n_rows()==C.n_cols(),
           dealii::ExcMessage("Result matrix does not match input matrix dimensions in MaMt."));

    A.scaled_mmult_add_scaled(C, B, 'c', 't', alpha, beta);
}

template <typename T, typename BW>
static void apply(::SciPAL::Matrix<T, BW> &result,
                  const typename BlasMatExp<T, BW>::MtMa& expr)
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;

    T alpha = 1.0;
    const Mtx & A = expr.l.l;
    const Mtx & B = expr.r.l;

    T beta = 0.0;

    Mtx & C = result;

    C.reinit(A.n_cols(), B.n_rows());
    Assert(A.n_rows()==B.n_cols(), dealii::ExcMessage("Input matrix dimensions do not match in MtMa."));
    Assert(A.n_cols()==C.n_rows() && B.n_rows()==C.n_cols(),
           dealii::ExcMessage("Result matrix does not match input matrix dimensions in MtMa."));

    A.scaled_mmult_add_scaled(C, B, 't', 'c', alpha, beta);
}

template <typename T, typename BW>
static void apply(::SciPAL::Matrix<T, BW> &result,
                  const typename BlasMatExp<T, BW>::dMM& expr)
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;
    typedef ::SciPAL::Vector<T, BW> Vtr;

    const Vtr & diag = expr.l.l;
    const Mtx & A = expr.r;


    Mtx & C = result;

    C.reinit(A.n_rows(), A.n_cols());

    Assert(diag.n_rows()==A.n_rows(), dealii::ExcMessage("Input matrix dimensions do not match in dMM."));

    BW::dgmm(CUBLAS_SIDE_LEFT,
             A.n_rows(), A.n_cols(),
             A.data(), A.leading_dim(),
             diag.data(), 1,
             C.data(), C.leading_dim());
}


template <typename T, typename BW>
static void apply(::SciPAL::Matrix<T, BW> &result,
                  const typename BlasMatExp<T, BW>::MdM& expr)
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;
    typedef ::SciPAL::Vector<T, BW> Vtr;

    const Vtr & diag = expr.r.l;
    const Mtx & A = expr.l;


    Mtx & C = result;

	 C.reinit(A.n_rows(), A.n_cols());

    Assert(diag.n_rows()==A.n_cols(), dealii::ExcMessage("Input matrix dimensions do not match in MdM."));

    BW::dgmm(CUBLAS_SIDE_RIGHT,
             A.n_rows(), A.n_cols(),
             A.data(), A.leading_dim(),
             diag.data(), 1,
             C.data(), C.leading_dim());
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
    Assert(A.n_cols()==B.n_rows(), dealii::ExcMessage("Input matrix dimensions do not match in sMM."));
    Assert(A.n_rows()==C.n_rows() && B.n_cols()==C.n_cols(),
           dealii::ExcMessage("Result matrix does not match input matrix dimensions in sMM."));
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

    Assert(A.n_cols()==B.n_rows(), dealii::ExcMessage("Input matrix dimensions do not match in MMaM."));
    Assert(A.n_rows()==C.n_rows() && B.n_cols()==C.n_cols(),
           dealii::ExcMessage("Result matrix does not match input matrix dimensions in MMaM."));

    A.scaled_mmult_add_scaled(C, B, 'n', 'n', alpha, beta);
}

template <typename T, typename BW>
static void apply(::SciPAL::Matrix<T, BW> & result,
                  const typename BlasMatExp<T, BW>::sMMaM& expr)
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;

    T alpha = expr.l.l.l;
    const Mtx & A = expr.l.l.r;
    const Mtx & B = expr.l.r;

    T beta = 1.0;

    Mtx & C = result;

    Assert(A.n_cols()==B.n_rows(), dealii::ExcMessage("Input matrix dimensions do not match in MMaM."));
    Assert(A.n_rows()==C.n_rows() && B.n_cols()==C.n_cols(),
           dealii::ExcMessage("Result matrix does not match input matrix dimensions in MMaM."));

    A.scaled_mmult_add_scaled(C, B, 'n', 'n', alpha, beta);
}
//typedef typename SciPAL::BinaryExpr<MtM, plus, Mtx> MtMpM;
template <typename T, typename BW>
static void apply(::SciPAL::Matrix<T, BW> & result,
                  const typename BlasMatExp<T, BW>::MtMpM& expr)
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;

    T alpha = 1.0, beta = 1.0;

    const Mtx & A = expr.l.l.l;
    const Mtx & B = expr.l.r;

#ifdef DEBUG
    const Mtx & C_in = expr.r;
#endif

    Mtx & C = result;

    Assert(&C_in == &C, dealii::ExcMessage("C_in not equal C in MtMpM."));

    Assert(A.n_rows()==B.n_rows(), dealii::ExcMessage("Input matrix dimensions do not match in MtMpM."));
    Assert(A.n_cols()==C_in.n_rows() && B.n_cols()==C_in.n_cols(),
           dealii::ExcMessage("Result matrix does not match input matrix dimensions in MtMpM."));

    A.scaled_mmult_add_scaled(C, B, 't', 'n', alpha, beta);
}


//typedef typename SciPAL::BinaryExpr<MMt, plus, Mtx> MMtpM;
template <typename T, typename BW>
static void apply(::SciPAL::Matrix<T, BW> & result,
                  const typename BlasMatExp<T, BW>::MMtpM& expr)
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;

    T alpha = 1.0, beta = 1.0;

    const Mtx & A = expr.l.l;
    const Mtx & B = expr.l.r.l;

    const Mtx & C_in = expr.r;

    Mtx & C = result;

    Assert(&C_in == &C, dealii::ExcMessage("C_in not equal C in MtMpM."));

    Assert(A.n_cols()==B.n_cols(), dealii::ExcMessage("Input matrix dimensions do not match in MtMpM."));
    Assert(A.n_rows()==C_in.n_rows() && B.n_rows()==C_in.n_cols(),
           dealii::ExcMessage("Result matrix does not match input matrix dimensions in MtMpM."));

    A.scaled_mmult_add_scaled(C, B, 'n', 't', alpha, beta);
}

//typedef typename SciPAL::BinaryExpr<MtM, plus, Mtx> MtMpM;
template <typename T, typename BW>
static void apply(::SciPAL::Matrix<T, BW> & result,
                  const typename BlasMatExp<T, BW>::MaMpM& expr)
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;

    T alpha = 1.0, beta = 1.0;

    const Mtx & A = expr.l.l.l;
    const Mtx & B = expr.l.r;
#ifdef DEBUG
    const Mtx & C_in = expr.r;
#endif

    Mtx & C = result;

    Assert(&C_in == &C, dealii::ExcMessage("C_in not equal C in MaMpM."));

    Assert(A.n_rows()==B.n_rows(), dealii::ExcMessage("Input matrix dimensions do not match in MtMpM."));
    Assert(A.n_cols()==C_in.n_rows() && B.n_cols()==C_in.n_cols(),
           dealii::ExcMessage("Result matrix does not match input matrix dimensions in MtMpM."));

    A.scaled_mmult_add_scaled(C, B, 'c', 'n', alpha, beta);
}


//typedef typename SciPAL::BinaryExpr<MMt, plus, Mtx> MMapM;
template <typename T, typename BW>
static void apply(::SciPAL::Matrix<T, BW> & result,
                  const typename BlasMatExp<T, BW>::MMapM& expr)
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;

    T alpha = 1.0, beta = 1.0;

    const Mtx & A = expr.l.l;
    const Mtx & B = expr.l.r.l;
#ifdef DEBUG
    const Mtx & C_in = expr.r;
#endif
    Mtx & C = result;

    Assert(&C_in == &C, dealii::ExcMessage("C_in not equal C in MMapM."));

    Assert(A.n_cols()==B.n_cols(), dealii::ExcMessage("Input matrix dimensions do not match in MMapM."));
    Assert(A.n_rows()==C_in.n_rows() && B.n_rows()==C_in.n_cols(),
           dealii::ExcMessage("Result matrix does not match input matrix dimensions in MMapM."));

    A.scaled_mmult_add_scaled(C, B, 'n', 'c', alpha, beta);
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

    Assert(A.n_cols()==B.n_rows(), dealii::ExcMessage("Input matrix dimensions do not match in MMaM."));
    Assert(A.n_rows()==C.n_rows() && B.n_cols()==C.n_cols(),
           dealii::ExcMessage("Result matrix does not match input matrix dimensions in MMaM."));

    A.scaled_mmult_add_scaled(C, B, 'n', 'n', alpha, beta);
}

template <typename T, typename BW>
static void apply(::SciPAL::Matrix<T, BW> & result,
                  const typename BlasMatExp<T, BW>::adM& expr)
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;

    T alpha = 1.0;
    const Mtx & A = expr.l;

    Mtx B(A.n_rows(), A.n_rows());
    for (unsigned int i = 0; i < B.n_rows(); i++)
        B(i, i, 1.0);

    T beta = 1.0;

    Mtx & C = result;
    C.reinit(A.n_cols(), A.n_rows());

    A.scaled_mmult_add_scaled(C, B, 'c', 'n', alpha, beta);
}

template <typename T, typename BW>
static void apply(::SciPAL::Matrix<T, BW> & result,
                  const typename BlasMatExp<T, BW>::tM& expr)
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;

    T alpha = 1.0;
    const Mtx & A = expr.l;

    Mtx B(A.n_rows(), A.n_rows());
    for (unsigned int i = 0; i < B.n_rows(); i++)
        B(i, i, 1.0);

    T beta = 1.0;

    Mtx & C = result;

    C.reinit(A.n_cols(), A.n_rows());
    A.scaled_mmult_add_scaled(C, B, 't', 'n', alpha, beta);
}

template <typename T, typename BW>
static void apply(::SciPAL::Matrix<T, BW> &result,
                  const typename BlasMatExp<T, BW>::s_MaM_M& expr)
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;

    T alpha = expr.l;
    const Mtx & A = expr.r.l.l.l;
    const Mtx & B = expr.r.l.r;
    const Mtx & C = expr.r.r;

    Mtx & D = result;

    D.reinit(A.n_cols(), C.n_cols());

    Assert(A.n_rows()==B.n_rows(), dealii::ExcMessage("Input matrix dimensions do not match in M_MM."));
    Assert(B.n_cols()==C.n_rows(), dealii::ExcMessage("Input matrix dimensions do not match in M_MM."));
    Assert(A.n_rows()==D.n_rows() && C.n_cols()==D.n_cols(),
           dealii::ExcMessage("Result matrix does not match input matrix dimensions in M_MM."));

    if (A.n_rows() > B.n_cols())
    {
        Mtx tmp(A.n_cols(), B.n_cols());
        tmp = alpha * adjoint(A) * B;
        D = tmp * C;
    }
    else
    {
        Mtx tmp(B.n_rows(), C.n_cols());
        tmp = B * C;
        D = alpha * adjoint(A) * tmp;
    }
}


template <typename T, typename BW>
static void apply(::SciPAL::Matrix<T, BW> &result,
                  const typename BlasMatExp<T, BW>::s_MM_M& expr)
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;

    T alpha = expr.l;
    const Mtx & A = expr.r.l.l;
    const Mtx & B = expr.r.l.r;
    const Mtx & C = expr.r.r;

    Mtx & D = result;

    D.reinit(A.n_rows(), C.n_cols());

    Assert(A.n_cols()==B.n_rows(), dealii::ExcMessage("Input matrix dimensions do not match in M_MM."));
    Assert(B.n_cols()==C.n_rows(), dealii::ExcMessage("Input matrix dimensions do not match in M_MM."));
    Assert(A.n_rows()==D.n_rows() && C.n_cols()==D.n_cols(),
           dealii::ExcMessage("Result matrix does not match input matrix dimensions in M_MM."));

    if (A.n_cols() > B.n_cols())
    {
        Mtx tmp(A.n_rows(), B.n_cols());
        tmp = alpha * A * B;
        D = tmp * C;
    }
    else
    {
        Mtx tmp(B.n_rows(), C.n_cols());
        tmp = B * C;
        D = alpha * A * tmp;
    }
}


template <typename T, typename BW>
static inline void apply(::SciPAL::Matrix<T, BW>& result,
                         const typename BlasMatExp<T, BW>::MM_M& expr)
{
    result = expr.l.l * (expr.l.r * expr.r);
}

template <typename T, typename BW>
static void apply(::SciPAL::Matrix<T, BW> &result,
                  const typename BlasMatExp<T, BW>::M_MM& expr)
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;

    const Mtx & A = expr.l;
    const Mtx & B = expr.r.l;
    const Mtx & C = expr.r.r;

    Mtx & D = result;

    D.reinit(A.n_rows(), C.n_cols());

    Assert(A.n_cols()==B.n_rows(), dealii::ExcMessage("Input matrix dimensions do not match in M_MM."));
    Assert(B.n_cols()==C.n_rows(), dealii::ExcMessage("Input matrix dimensions do not match in M_MM."));
    Assert(A.n_rows()==D.n_rows() && C.n_cols()==D.n_cols(),
           dealii::ExcMessage("Result matrix does not match input matrix dimensions in M_MM."));

    if (A.n_cols() > B.n_cols())
    {
        Mtx tmp(A.n_rows(), B.n_cols());
        tmp = A * B;
        D = tmp * C;
    }
    else
    {
        Mtx tmp(B.n_rows(), C.n_cols());
        tmp = B * C;
        D = A * tmp;
    }
}

template <typename T, typename BW>
static inline void apply(::SciPAL::Matrix<T, BW>& result,
                         const typename BlasMatExp<T, BW>::sMM_M& expr)
{
    result = (expr.l.l.l * expr.l.l.r) * (expr.l.r * expr.r);
}

template <typename T, typename BW>
static inline void apply(::SciPAL::Matrix<T, BW>& result,
                         const typename BlasMatExp<T, BW>::sMaM_M& expr)
{
    result = (expr.l.l.l * expr.l.l.r) * (expr.l.r * expr.r);
}


template <typename T, typename BW>
static void apply(::SciPAL::Matrix<T, BW> &result,
                  const typename BlasMatExp<T, BW>::sM_MM& expr)
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;

    T alpha = expr.l.l;
    const Mtx & A = expr.l.r;
    const Mtx & B = expr.r.l;
    const Mtx & C = expr.r.r;

    Mtx & D = result;

    D.reinit(A.n_rows(), C.n_cols());

    Assert(A.n_cols()==B.n_rows(), dealii::ExcMessage("Input matrix dimensions do not match in M_MM."));
    Assert(B.n_cols()==C.n_rows(), dealii::ExcMessage("Input matrix dimensions do not match in M_MM."));
    Assert(A.n_rows()==D.n_rows() && C.n_cols()==D.n_cols(),
           dealii::ExcMessage("Result matrix does not match input matrix dimensions in M_MM."));

    if (A.n_cols() > B.n_cols())
    {
        Mtx tmp(A.n_rows(), B.n_cols());
        tmp = alpha * A * B;
        D = tmp * C;
    }
    else
    {
        Mtx tmp(B.n_rows(), C.n_cols());
        tmp = B * C;
        D = alpha * A * tmp;
    }
}

template <typename T, typename BW>
static inline void apply(::SciPAL::Matrix<T, BW>& result,
                         const typename BlasMatExp<T, BW>::MM_Mt& expr)
{
    result = expr.l.l * (expr.l.r * expr.r);
}

template <typename T, typename BW>
static void apply(::SciPAL::Matrix<T, BW> &result,
                  const typename BlasMatExp<T, BW>::M_MMt& expr)
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;

    const Mtx & A = expr.l;
    const Mtx & B = expr.r.l;
    const Mtx & C = expr.r.r.l;

    Mtx & D = result;

    D.reinit(A.n_rows(), C.n_rows());

    Assert(A.n_cols()==B.n_rows(), dealii::ExcMessage("Input matrix dimensions do not match in M_MM."));
    Assert(B.n_cols()==C.n_cols(), dealii::ExcMessage("Input matrix dimensions do not match in M_MM."));
    Assert(A.n_rows()==D.n_rows() && C.n_rows()==D.n_cols(),
           dealii::ExcMessage("Result matrix does not match input matrix dimensions in M_MM."));

    if (A.n_cols() > B.n_cols())
    {
        Mtx tmp(A.n_rows(), B.n_cols());
        tmp = A * B;
        D = tmp * SciPAL::transpose(C);
    }
    else
    {
        Mtx tmp(B.n_rows(), C.n_rows());
        tmp = B * SciPAL::transpose(C);
        D = A * tmp;
    }
}


template <typename T, typename BW>
static inline void apply(::SciPAL::Matrix<T, BW>& result,
                         const typename BlasMatExp<T, BW>::MtM_M& expr)
{
    result = expr.l.l * (expr.l.r * expr.r);
}

template <typename T, typename BW>
static void apply(::SciPAL::Matrix<T, BW> &result,
                  const typename BlasMatExp<T, BW>::Mt_MM& expr)
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;

    const Mtx & A = expr.l.l;
    const Mtx & B = expr.r.l;
    const Mtx & C = expr.r.r;

    Mtx & D = result;

    D.reinit(A.n_cols(), C.n_cols());

    Assert(A.n_rows()==B.n_rows(), dealii::ExcMessage("Input matrix dimensions do not match in M_MM."));
    Assert(B.n_cols()==C.n_rows(), dealii::ExcMessage("Input matrix dimensions do not match in M_MM."));
    Assert(A.n_cols()==D.n_rows() && C.n_cols()==D.n_cols(),
           dealii::ExcMessage("Result matrix does not match input matrix dimensions in M_MM."));

    if (A.n_rows() > B.n_cols())
    {
        Mtx tmp(A.n_cols(), B.n_cols());
        tmp = SciPAL::transpose(A) * B;
        D = tmp * C;
    }
    else
    {
        Mtx tmp(B.n_rows(), C.n_cols());
        tmp = B * C;
        D = SciPAL::transpose(A) * tmp;
    }
}
template <typename T, typename BW>
static inline void apply(::SciPAL::Matrix<T, BW>& result,
                         const typename BlasMatExp<T, BW>::MM_Ma& expr)
{
    result = expr.l.l * (expr.l.r * expr.r);
}

template <typename T, typename BW>
static void apply(::SciPAL::Matrix<T, BW> &result,
                  const typename BlasMatExp<T, BW>::M_MMa& expr)
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;

    const Mtx & A = expr.l;
    const Mtx & B = expr.r.l;
    const Mtx & C = expr.r.r.l;

    Mtx & D = result;

    D.reinit(A.n_rows(), C.n_rows());

    Assert(A.n_cols()==B.n_rows(), dealii::ExcMessage("Input matrix dimensions do not match in M_MM."));
    Assert(B.n_cols()==C.n_cols(), dealii::ExcMessage("Input matrix dimensions do not match in M_MM."));
    Assert(A.n_rows()==D.n_rows() && C.n_rows()==D.n_cols(),
           dealii::ExcMessage("Result matrix does not match input matrix dimensions in M_MM."));

    if (A.n_cols() > B.n_cols())
    {
        Mtx tmp(A.n_rows(), B.n_cols());
        tmp = A * B;
        D = tmp * SciPAL::adjoint(C);
    }
    else
    {
        Mtx tmp(B.n_rows(), C.n_rows());
        tmp = B * SciPAL::adjoint(C);
        D = A * tmp;
    }
}


template <typename T, typename BW>
static inline void apply(::SciPAL::Matrix<T, BW>& result,
                         const typename BlasMatExp<T, BW>::MaM_M& expr)
{
    result = expr.l.l * (expr.l.r * expr.r);
}

template <typename T, typename BW>
static void apply(::SciPAL::Matrix<T, BW> &result,
                  const typename BlasMatExp<T, BW>::Ma_MM& expr)
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;

    const Mtx & A = expr.l.l;
    const Mtx & B = expr.r.l;
    const Mtx & C = expr.r.r;

    Mtx & D = result;

    D.reinit(A.n_cols(), C.n_cols());

    Assert(A.n_rows()==B.n_rows(), dealii::ExcMessage("Input matrix dimensions do not match in M_MM."));
    Assert(B.n_cols()==C.n_rows(), dealii::ExcMessage("Input matrix dimensions do not match in M_MM."));
    Assert(A.n_cols()==D.n_rows() && C.n_cols()==D.n_cols(),
           dealii::ExcMessage("Result matrix does not match input matrix dimensions in M_MM."));

    if (A.n_rows() > B.n_cols())
    {
        Mtx tmp(A.n_cols(), B.n_cols());
        tmp = SciPAL::adjoint(A) * B;
        D = tmp * C;
    }
    else
    {
        Mtx tmp(B.n_rows(), C.n_cols());
        tmp = B * C;
        D = SciPAL::adjoint(A) * tmp;
    }
}
} // END namespace LAOOperations

// trace(A^d * B)
// in this method the trace of a matrix product is calculated. In this case it is better to calculate the diagonal entries via vector-vector multiplication.
// ATTENTION: NOT YET WORKING FOR COMPLEX NUMBERS!!! (TODO; FIXME)
template<typename T1, typename T2>
 auto
trace(const SciPAL::BinaryExpr<SciPAL::UnaryExpr<T1, SciPAL::expr_adjoint>, SciPAL::mult, T2>& Ax)
-> decltype (Ax.l.l(0,0) * Ax.r(0,0))
{
    if (Ax.l.l.n_rows() != Ax.r.n_rows())
    {
        std::cerr << "Error in trace(adjoint(A)*B). Exiting" << std::endl;
        std::exit(1);
    }
    decltype (Ax.l.l(0,0) * Ax.r(0,0)) Output = 0;

    for (unsigned int i = 0; i < Ax.l.l.n_rows() && i < Ax.r.n_rows(); i++)
    {
        // We need the rows of the second matrix and
        const SciPAL::RowVectorView<typename T2::Number, typename T2::blas_wrapper_type, const T2> A_row(Ax.r, i, 0);
        // we need the rows of the first matrix
        const SciPAL::RowVectorView<typename T1::Number, typename T1::blas_wrapper_type, const T1> A_col(Ax.l.l, i, 0);

        Literal<decltype (Ax.l.l(0,0) * Ax.r(0,0)), typename T1::blas_wrapper_type> tmp = 0;
        tmp = A_row * A_col;
        Output += tmp;
    }
    return Output;
}

template<typename T1>
inline
const typename T1::NumberType
trace(const T1 & A)
{
    typename T1::NumberType Output = 0;
    for (unsigned int i = 0; i < A.n_cols() && A.n_rows(); i++)
    {
        SciPAL::ColVectorView<typename T1::NumberType, typename T1::BW_TYPE, T1> A_col(A, i);
        SciPAL::RowVectorView<typename T1::NumberType, typename T1::BW_TYPE, T1> A_row(A, i);
        Output += A_col * A_row;
    }
    return Output;
}

} // END namespace SciPAL


#endif // MATRIXOPERATIONS_H
