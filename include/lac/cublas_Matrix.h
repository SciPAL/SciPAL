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


#ifndef cublas_Matrix_H
#define cublas_Matrix_H


#include <iomanip>
#include <cmath>
#include <QtGlobal>

#include <deal.II/base/subscriptor.h>
#include <deal.II/lac/identity_matrix.h>

#include <lac/Array.h>
#include <lac/Expr.h>
#include <lac/BlasMatrixOperations.h>

// this is from SciPAL/include/
//
#include <lac/cublas_Vector.h>
#include <lac/Shape.h>
#include <base/ArchTraits.h>
#include <base/CUDA_error_check.h>
#include <base/Zero_One_Traits.h>
#include <base/ForewardDeclarations.h>



namespace SciPAL {

// IDEA: enum MatrixStorage { general, banded, symmetric, banded_packed, symmetric_packed };

// IDEA:
/*
template<typename T, typename BW>
class SymmetricMatrix : public Matrix<T, BW, symmetric> {


};
*/

// @sect3{Class: Matrix}
//! A class for full dense matrices which is compatible with the Krylov solvers of the deal.II library and which allows to be stored on a GPU.
template<typename T, typename BW> // IDEA:, MatrixStorage my_storage>
class Matrix
        :
        public SciPAL::Expr<Matrix<T, BW> >,
//        protected SciPAL::Array<T, BW>,
        public dealii::Subscriptor,
        public SciPAL::Shape<T, BW, matrix>
{

    friend class Vector<T, BW>;

    friend class SubVectorView<T, BW, Matrix<T, BW> >;

    friend class SubMatrixView<T, BW>;

    friend class FullMatrixAccessor<T>;

public:

    typedef T Number;

    typedef T value_type;

    typedef BW blas_wrapper_type;

    typedef Matrix<T, BW> Type;

    typedef const Type& ConstHandle;

    typedef SciPAL::ShapeData<T> DevType;

    typedef SciPAL::Shape<T, BW, matrix> MyShape;

    static const ParallelArch arch = BW::arch;

   // IDEA: MatrixStorage my_storage_scheme = ms;

    static const EType I_am = leafE;

    Matrix();

    Matrix(int n_rows, int n_cols);

    template<typename T2>
    Matrix(const unsigned int n_rows,
               const unsigned int n_cols,
               const std::vector<T2> &src);

    //FIX ME
//    template<typename X>
//    Matrix (const ::SciPAL::Expr<X> & e);

    Matrix(const Matrix<T, BW>& other);

//    Matrix(const Shape<T, BW, matrix>& other)
//        : Subscriptor(),
//          MyShape(other)
//    {}


// private:
//    Matrix(const Matrix & other) {
//        AssertThrow(false,
//                    dealii::ExcMessage("There should never be a need to invoke the copy ctor of a matrix!"));
//    }
public:
    template<typename T2>
    Matrix(const FullMatrixAccessor<T2> & matrix);

    Matrix(const dealii::IdentityMatrix & Id);

    void reinit(int n_rows, int n_cols);

    
    void swap(Matrix& other)
    {	
	Assert(false, dealii::ExcMessage("Seems not to be working, disabled... - Holger"));
        this->storage.swap(other);
        this->MyShape::swap(other);
    }

    template<typename T2>
    Matrix<T, BW> & operator= (const FullMatrixAccessor<T2> & matrix);

    Matrix<T, BW> & operator= (const dealii::IdentityMatrix & Id);

    //! Generate a deep copy of @p other
    Matrix<T, BW> & operator = (const Matrix<T, BW> & other)
    {
        this->reinit(other.n_rows(), other.n_cols());
        // element-wise copy of array.
        int inc_src  = 1;
        int inc_this = 1;
        BW::copy(this->size(), other.data(), inc_src,
                 this->data(), inc_this);
        return *this;
    }
	//! ToDo: rewrite with shapes
    //! Generate a deep copy of @p other
    Matrix<T, BW> & operator = (const SubMatrixView<T, BW> & other)
    {
        this->reinit(other.n_rows_active(), other.n_cols_active());
        int inc_src  = 1;
        int inc_this = 1;
        // TODO: i guess this could be optimized by checking the number of rows and cols
        // and therefore copying row or colwise to minimize the number of function calls

        // Matrix stored col major, so copy columnwise:
        for(unsigned int c = 0; c < other.n_cols_active(); c++)
        {
            BW::copy(other.n_rows_active(),
                     other.view_begin + c * other.leading_dim(),
                     inc_src,
                     this->data() + c * this->leading_dim,
                     inc_this);
        }
        return *this;
    }


    //! Generate a deep copy of @p other with different blas_wrapper_type
    template <typename BW2>
    Matrix<T, BW> & operator = (const Matrix<T, BW2> & other)
    {

        this->reinit(other.n_rows(), other.n_cols());

        //! copy from cublas matrix to blas matrix -> GetMatrix
        //! TODO: what is with asyn copy?
        if(typeid(BW) == typeid(blas) && typeid(BW2) == typeid(cublas) )
        {
            cublas::GetMatrix(other.n_rows(), other.n_cols(),
                              other.data(),
                              other.leading_dim,
                              this->data(),
                              this->leading_dim);
        }

        //! copy from cublas matrix to blas matrix -> SetMatrix
        //! TODO: what is with asyn copy?
        if(typeid(BW) == typeid(cublas) && typeid(BW2) == typeid(blas) )
        {
            cublas::SetMatrix(other.n_rows(), other.n_cols(),
                              other.data(),
                              other.leading_dim,
                              this->data(),
                              this->leading_dim);
        }
#ifdef DEBUG
        std::cout<<__FUNCTION__<<std::endl;
#endif
        return *this;
    }

    template<typename X>
    Matrix & operator= (const ::SciPAL::Expr<X> & e);

#ifndef nUSE_GENERIC_ADDITION
    // ToDO test if working correct
    template<typename X>
    Matrix & operator+= (const ::SciPAL::Expr<X> & e);
#else
    //! A += B
    Matrix<T, BW> & operator += (const Matrix<T, BW> & other);

    //! A += a * B * C * D
    Matrix<T, BW> & operator += (const SciPAL::BinaryExpr<Literal<T, BW>, mult, SciPAL::BinaryExpr<SciPAL::BinaryExpr<SciPAL::Matrix<T, BW>,  mult, SciPAL::Matrix<T, BW>>,  mult, SciPAL::Matrix<T, BW>>>& expr);

    //! A += a * adjoint(B) * C * D
    Matrix<T, BW> & operator += (const SciPAL::BinaryExpr<Literal<T, BW>, mult, SciPAL::BinaryExpr<SciPAL::BinaryExpr<SciPAL::UnaryExpr<SciPAL::Matrix<T, BW>, SciPAL::expr_adjoint>,  mult, SciPAL::Matrix<T, BW>>,  mult, SciPAL::Matrix<T, BW>>>& expr);

    //! A += a * B * C * adjoint(D)
    Matrix<T, BW> & operator += (const SciPAL::BinaryExpr<Literal<T, BW>, mult, SciPAL::BinaryExpr<SciPAL::BinaryExpr<SciPAL::Matrix<T, BW>,  mult, SciPAL::Matrix<T, BW>>,  mult, SciPAL::UnaryExpr<SciPAL::Matrix<T, BW>, SciPAL::expr_adjoint>>>& expr);

    //! A += a * B * C * transpose(D)
    Matrix<T, BW> & operator += (const SciPAL::BinaryExpr<Literal<T, BW>, mult, SciPAL::BinaryExpr<SciPAL::BinaryExpr<SciPAL::Matrix<T, BW>,  mult, SciPAL::Matrix<T, BW>>,  mult, SciPAL::UnaryExpr<SciPAL::Matrix<T, BW>, SciPAL::expr_transpose>>>& expr);

    //! A += trans(B) * C
    Matrix<T, BW> & operator += (const SciPAL::BinaryExpr<SciPAL::UnaryExpr<SciPAL::Matrix<T, BW>, SciPAL::expr_transpose>, mult, SciPAL::Matrix<T, BW>>& expr);

    //! A += B * trans(C)
    Matrix<T, BW> & operator += (const SciPAL::BinaryExpr<SciPAL::Matrix<T, BW>, mult, SciPAL::UnaryExpr<SciPAL::Matrix<T, BW>, SciPAL::expr_transpose>>& expr);
    //! A += B * adjoint(C)
    Matrix<T, BW> & operator += (const SciPAL::BinaryExpr<SciPAL::Matrix<T, BW>, mult, SciPAL::UnaryExpr<SciPAL::Matrix<T, BW>, SciPAL::expr_adjoint>>& expr);

    //! A += adjoint(B) * C
    Matrix<T, BW> & operator += (const SciPAL::BinaryExpr<SciPAL::UnaryExpr<SciPAL::Matrix<T, BW>, SciPAL::expr_adjoint>, mult, SciPAL::Matrix<T, BW>>& expr);


    //! A += B * C
    Matrix<T, BW> & operator += (const SciPAL::BinaryExpr<SciPAL::Matrix<T, BW>, mult, SciPAL::Matrix<T, BW>>& expr);
#endif

    Matrix<T, BW> & operator -= (const Matrix<T, BW> & other);

    template<typename T2>
    Matrix<T, BW> & operator *= (const T2 scale);

    template<typename VECTOR1, typename VECTOR2>
    void vmult(VECTOR1& dst, const VECTOR2& src) const;

    template<typename VECTOR1, typename VECTOR2>
    void Tvmult(VECTOR1& dst, const VECTOR2& src) const;

    void scaled_vmult_add_scaled(T alpha,
                                 bool A_is_transpose,
                                 Vector<T, BW>& dst /*y*/,
                                 T beta,
                                 const Vector<T, BW>& src /*x*/);


    //! $dst = this \cdot src$
    void mmult(Matrix<T, BW>& dst, const Matrix<T, BW>& src) const;

    void mmult(SubMatrixView<T, BW>& dst, const Matrix<T, BW>& src) const;

    //! $dst = this \cdot src^T$
    void mTmult(Matrix<T, BW>& dst, const Matrix<T, BW>& src) const;

    //! $dst = this^T \cdot src$
    void Tmmult(Matrix<T, BW>& dst, const Matrix<T, BW>& src) const;

    //! $dst = this^T \cdot src^T$
    void TmTmult(Matrix<T, BW>& dst, const Matrix<T, BW>& src) const;

    void scaled_mmult_add_scaled( Matrix<T, BW>& dst, const Matrix<T, BW>& src,
                                  char transpose_A='n', char transpose_B='n',
                                  T alpha = 1.0, T beta = 0) const;

    template<typename VECTOR1, typename VECTOR2>
    void add_scaled_outer_product(T alpha, const VECTOR1& col, const VECTOR2& row);



    inline size_t n_rows() const { return this->MyShape::n_rows_active(); }

    inline size_t n_cols() const { return this->MyShape::n_cols_active(); }

    void print(std::ostream &Output = std::cout) const;

    T operator () (const unsigned int i, const unsigned int j) const;

    void operator () (const unsigned int i, const unsigned int j, T data);

    typename PrecisionTraits<T, BW::arch>::NumberType l2_norm() const;

    T sum() const;

    inline MyShape & shape() { return *this; }
	inline const MyShape & shape() const { return *this; }

private:

    //Matrix<T, BW> & operator = (const Array<T, BW> & src);
};

}

// @sect3{SciPAL::Matrix methods}
//!
// @sect4{Constructors}
//! The default ctor creates a matrix object with zero entries.
//! You can pass it around but before actually using it you have to reinitialize it such that it contains a finite number of elements.
template<typename T, typename BW>
SciPAL::Matrix<T, BW>::Matrix()
    :
      MyShape()
{}

//! Allocate a matrix of @p n_rows and @p n_cols.
//! @param n_rows : Number of rows.
//! @param n_cols : Number of columns.
//! TO DO: Initialize the GPU-side memory of a Matrix with zeros.
template<typename T, typename BW>
SciPAL::Matrix<T, BW>::Matrix(int n_rows, int n_cols)
    :
      MyShape(0, n_rows, /*active rows*/
              0, n_cols/*active cols*/)
{}


//! Allocate a matrix of @p n_rows and @p n_cols and fill it with given values @p src.
//! @param n_rows : Number of rows.
//! @param n_cols : Number of columns.
//! @param src : values of the matrix elements in lexicographic order.
template<typename T, typename BW>
template<typename T2>
SciPAL::Matrix<T, BW>::Matrix(const unsigned int n_rows,
                                const unsigned int n_cols,
                                const std::vector<T2> & src)
    :
      MyShape(0, n_rows,
              0, n_cols)
{
    const T2 * const tmp_ptr = &src[0];

    T * this_data = this->data();

    BW::SetMatrix(n_rows, n_cols,
                  tmp_ptr,
                  n_rows,
                  this_data,
                  this->leading_dim);

}

//! Assign the result of a linear algebraic expression to a matrix.
//! @param e : Expression to evaluate.
template<typename T, typename BW>
template<typename X>
SciPAL::Matrix<T, BW> & SciPAL::Matrix<T, BW>::operator =
//(const typename ::SciPAL::BlasMatExp<T, BW>::sMMaM& e)
(const ::SciPAL::Expr<X> & e)
{
#ifdef DEBUG_MATRIX
    std::cout << "line :" << __LINE__ << ", Matrix<T,BW>" << std::endl;
    print_expr_info(__PRETTY_FUNCTION__);
#endif
    ::SciPAL::LAOOperations::apply(*this, ~e);

    return *this;
}

//! Copy ctor with constructs a matrix from an identity matrix.
//! @param Id : identity matrix which serves as source.
template<typename T, typename BW>
SciPAL::Matrix<T, BW>::Matrix(const dealii::IdentityMatrix & Id)
    :
      MyShape()
{
    *this = Id;
}



//! Copy ctor with constructs a matrix from an another matrix.
//! @param other : Matrix which serves as source.
template<typename T, typename BW>
SciPAL::Matrix<T, BW>::Matrix(const Matrix<T, BW> & other)
    :
      dealii::Subscriptor(),
      MyShape()
{
    *this = other;
}

//! Copy ctor with constructs a matrix from an another (host-side) matrix.
//! @param other : Matrix which serves as source.
template<typename T, typename BW>
template<typename T2>
SciPAL::Matrix<T, BW>::Matrix(const FullMatrixAccessor<T2> & other)
{
    (*this) = other;
}

//! Reinitialize the matrix to the given number of rows and columns.
template<typename T, typename BW>
void SciPAL::Matrix<T, BW>::reinit(int n_rows, int n_cols)
{
   this->MyShape::reinit(0, n_rows,
                         0, n_cols,
                         1 /*unit stride*/);
}

//! Initialize a matrix from an identity matrix.
//! @param Id : identity matrix which provides the information about the size of the matrix.
template<typename T, typename BW>
SciPAL::Matrix<T, BW> &
SciPAL::Matrix<T, BW>::operator = (const dealii::IdentityMatrix & Id)
{
    int __n_rows = Id.m();
    int __n_cols = Id.n();

    // Zeilen- und Spaltenzahl einer Identitaetsmatrix sind gleich
    // und entsprechen der Anzahl der Freiheitsgrade.
    int __n_dofs   = Id.m();

    this->MyShape::reinit(0, __n_rows,
                          0, __n_cols,
                          1 /*unit stride*/);


    FullMatrixAccessor<T> tmp_id(Id);

    // Laut CUBLAS Doku wird die 'leading dimension' immer durch die Anzahl
    // der Zeilen angegeben.
    // Ausserdem muss man sich bei Einheitsmatrizen nicht darum kuemmern,
    // ob sie column-major oder row-major sind.
    const T *  id_val = tmp_id.val();
    T * dst_val = this->data();

    BW::SetMatrix(__n_dofs, __n_dofs,
                  id_val, __n_dofs,
                  dst_val, this->leading_dim );

    return *this;
}

//! Copy the contents from a FullMatrixAccessor into a SciPAL::Matrix. FullMatrixAccessor is equivalent to a dealii::FullMatrix.
//! In case of cublas this assignment takes care of copying data from the host side to the device side.
//! @param src_matrix : full matrix in column-major format.
template<typename T, typename BW>
template<typename T2>
SciPAL::Matrix<T, BW> &
SciPAL::Matrix<T, BW>::operator = (const FullMatrixAccessor<T2> & src_matrix)
{
    Assert(src_matrix.is_column_major(),
           dealii::ExcMessage("SciPAL:Matrix expects a matrix in column major"
                              " format as source.") );


    int nr = src_matrix.n_rows();
    int nc = src_matrix.n_cols();

    this->MyShape::reinit(0, nr, /*active rows*/
                          0, nc, /*active cols*/
                          1 /*unit stride*/);



    const T * tmp_src = reinterpret_cast<const T*>(src_matrix.val());
    // Laut CUBLAS Doku wird die 'leading dimension' immer durch die Anzahl
    // der Zeilen angegeben.
    T * tmp_dst = this->data();
    BW::SetMatrix(nr, nc, tmp_src, nr, tmp_dst, nr);

    return *this;
}
#ifdef nUSE_GENERIC_ADDITION
// @sect4{Operator: +=}
//!
//! Incremental, element-wise addition of two matrices.
//! @param other : Matrix to add. Must have the same dimensions.
template<typename T, typename BW>
SciPAL::Matrix<T, BW> &
SciPAL::Matrix<T, BW>::operator += (const SciPAL::BinaryExpr<Literal<T, BW>, mult, SciPAL::BinaryExpr<SciPAL::BinaryExpr<SciPAL::Matrix<T, BW>,  mult, SciPAL::Matrix<T, BW>>,  mult, SciPAL::Matrix<T, BW>>>& expr)
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;

    T alpha = expr.l;
    const Mtx & A = expr.r.l.l;
    const Mtx & B = expr.r.l.r;
    const Mtx & C = expr.r.r;

    Assert((this->n_rows() == A.n_rows())
           && (this->n_cols() == C.n_cols()),
           dealii::ExcMessage("Dimension mismatch"));

    Assert(A.n_cols()==B.n_rows(), dealii::ExcMessage("Input matrix dimensions do not match in operator+=(s_MM_M)."));
    Assert(B.n_cols()==C.n_rows(), dealii::ExcMessage("Input matrix dimensions do not match in operator+=(s_MM_M)."));

    if (A.n_cols() > B.n_cols())
    {
        Mtx tmp(A.n_rows(), B.n_cols());
        tmp = alpha * A * B;
        *this = tmp * C + *this;
    }
    else
    {
        Mtx tmp(B.n_rows(), C.n_cols());
        tmp = B * C;
        *this = alpha * A * tmp + *this;
    }

    return *this;
}

// @sect4{Operator: +=}
//!
//! Incremental, element-wise addition of two matrices.
//! @param other : Matrix to add. Must have the same dimensions.
template<typename T, typename BW>
SciPAL::Matrix<T, BW> &
SciPAL::Matrix<T, BW>::operator += (const SciPAL::BinaryExpr<Literal<T, BW>, mult, SciPAL::BinaryExpr<SciPAL::BinaryExpr<SciPAL::UnaryExpr<SciPAL::Matrix<T, BW>, SciPAL::expr_adjoint>,  mult, SciPAL::Matrix<T, BW>>,  mult, SciPAL::Matrix<T, BW>>>& expr)
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;

    T alpha = expr.l;
    const Mtx & A = expr.r.l.l.l;
    const Mtx & B = expr.r.l.r;
    const Mtx & C = expr.r.r;

    Assert((this->n_rows() == A.n_cols())
           && (this->n_cols() == C.n_cols()),
           dealii::ExcMessage("Dimension mismatch"));

    Assert(A.n_rows()==B.n_rows(), dealii::ExcMessage("Input matrix dimensions do not match in operator+=(s_MaM_M)."));
    Assert(B.n_cols()==C.n_rows(), dealii::ExcMessage("Input matrix dimensions do not match in operator+=(s_MaM_M)."));

    Mtx tmp(A.n_cols(), B.n_cols());
    tmp = SciPAL::adjoint(A) * B;
    *this = alpha * tmp * C + *this;

    return *this;
}


// @sect4{Operator: +=}
//!
//! Incremental, element-wise addition of two matrices.
//! @param other : Matrix to add. Must have the same dimensions.
template<typename T, typename BW>
SciPAL::Matrix<T, BW> &
SciPAL::Matrix<T, BW>::operator += (const SciPAL::BinaryExpr<Literal<T, BW>, mult, SciPAL::BinaryExpr<SciPAL::BinaryExpr<SciPAL::Matrix<T, BW>,  mult, SciPAL::Matrix<T, BW>>,  mult, SciPAL::UnaryExpr<SciPAL::Matrix<T, BW>, SciPAL::expr_adjoint>>>& expr)
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;

    T alpha = expr.l;
    const Mtx & A = expr.r.l.l;
    const Mtx & B = expr.r.l.r;
    const Mtx & C = expr.r.r.l;

    Assert((this->n_rows() == A.n_rows())
           && (this->n_cols() == C.n_rows()),
           dealii::ExcMessage("Dimension mismatch"));

    Assert(A.n_cols()==B.n_rows(), dealii::ExcMessage("Input matrix dimensions do not match in operator+=(s_MaM_M)."));
    Assert(B.n_cols()==C.n_cols(), dealii::ExcMessage("Input matrix dimensions do not match in operator+=(s_MaM_M)."));

    Mtx tmp(B.n_rows(), C.n_rows());
    tmp = B * SciPAL::adjoint(C);
    *this = alpha * A * tmp + *this;

    return *this;
}

// @sect4{Operator: +=}
//!
//! Incremental, element-wise addition of two matrices.
//! @param other : Matrix to add. Must have the same dimensions.
template<typename T, typename BW>
SciPAL::Matrix<T, BW> &
SciPAL::Matrix<T, BW>::operator += (const SciPAL::BinaryExpr<Literal<T, BW>, mult, SciPAL::BinaryExpr<SciPAL::BinaryExpr<SciPAL::Matrix<T, BW>,  mult, SciPAL::Matrix<T, BW>>,  mult, SciPAL::UnaryExpr<SciPAL::Matrix<T, BW>, SciPAL::expr_transpose>>>& expr)
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;

    T alpha = expr.l;
    const Mtx & A = expr.r.l.l;
    const Mtx & B = expr.r.l.r;
    const Mtx & C = expr.r.r.l;

    Assert((this->n_rows() == A.n_rows())
           && (this->n_cols() == C.n_rows()),
           dealii::ExcMessage("Dimension mismatch"));

    Assert(A.n_cols()==B.n_rows(), dealii::ExcMessage("Input matrix dimensions do not match in operator+=(s_MaM_M)."));
    Assert(B.n_cols()==C.n_cols(), dealii::ExcMessage("Input matrix dimensions do not match in operator+=(s_MaM_M)."));

    Mtx tmp(B.n_rows(), C.n_rows());
    tmp = B * SciPAL::transpose(C);
    *this = alpha * A * tmp + *this;

    return *this;
}

// @sect4{Operator: +=}
//!
//! Incremental, element-wise addition of two matrices.
//! @param other : Matrix to add. Must have the same dimensions.
template<typename T, typename BW>
SciPAL::Matrix<T, BW> &
SciPAL::Matrix<T, BW>::operator += (const Matrix<T, BW> & other)
{
    Assert((this->n_rows() == other.n_rows())
           && (this->n_cols() == other.n_cols()),
           dealii::ExcMessage("Dimension mismatch"));

    int n = this->size();
    One<T> one;
    T alpha = one();

    const T * const x = other.data();
    int incx = 1;

    T * y = this->data();
    int incy = 1;

    BW::axpy(n, alpha, x, incx, y, incy);

    return *this;
}


// @sect4{Operator: +=}
//!
//! Incremental, element-wise addition of two matrices.
//! @param other : Matrix to add. Must have the same dimensions.
//! A += trans(B) * C
template<typename T, typename BW>
SciPAL::Matrix<T, BW> &
SciPAL::Matrix<T, BW>::operator += (const SciPAL::BinaryExpr<SciPAL::UnaryExpr<SciPAL::Matrix<T, BW>, SciPAL::expr_transpose>, mult, SciPAL::Matrix<T, BW>>& expr)
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;

    const Mtx & A = expr.l.l;
    const Mtx & B = expr.r;

    Assert((this->n_rows() == A.n_cols())
           && (this->n_cols() == B.n_cols()),
           dealii::ExcMessage("Dimension mismatch"));

    Assert(A.n_rows()==B.n_rows(), dealii::ExcMessage("Input matrix dimensions do not match in operator+=(MtM)."));

    *this = SciPAL::transpose(A) * B + *this;

    return *this;
}

// @sect4{Operator: +=}
//!
//! Incremental, element-wise addition of two matrices.
//! @param other : Matrix to add. Must have the same dimensions.
//! A += adjoint(B) * C
template<typename T, typename BW>
SciPAL::Matrix<T, BW> &
SciPAL::Matrix<T, BW>::operator += (const SciPAL::BinaryExpr<SciPAL::UnaryExpr<SciPAL::Matrix<T, BW>, SciPAL::expr_adjoint>, mult, SciPAL::Matrix<T, BW>>& expr)
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;

    const Mtx & A = expr.l.l;
    const Mtx & B = expr.r;

    Assert((this->n_rows() == A.n_cols())
           && (this->n_cols() == B.n_cols()),
           dealii::ExcMessage("Dimension mismatch"));

    Assert(A.n_rows()==B.n_rows(), dealii::ExcMessage("Input matrix dimensions do not match in operator+=(MaM)."));

    *this = SciPAL::adjoint(A) * B + *this;

    return *this;
}

// @sect4{Operator: +=}
//!
//! Incremental, element-wise addition of two matrices.
//! @param other : Matrix to add. Must have the same dimensions.
//! A += B * C
template<typename T, typename BW>
SciPAL::Matrix<T, BW> &
SciPAL::Matrix<T, BW>::operator += (const SciPAL::BinaryExpr<SciPAL::Matrix<T, BW>, mult, SciPAL::Matrix<T, BW>>& expr)
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;

    const Mtx & A = expr.l;
    const Mtx & B = expr.r;

    Assert((this->n_rows() == A.n_rows())
           && (this->n_cols() == B.n_cols()),
           dealii::ExcMessage("Dimension mismatch"));

    Assert(A.n_cols()==B.n_rows(), dealii::ExcMessage("Input matrix dimensions do not match in operator+=(MM)."));

    *this = A * B + *this;

    return *this;
}

// @sect4{Operator: +=}
//!
//! Incremental, element-wise addition of two matrices.
//! @param other : Matrix to add. Must have the same dimensions.
//! A += B * trans(C)
template<typename T, typename BW>
SciPAL::Matrix<T, BW> &
SciPAL::Matrix<T, BW>::operator += (const SciPAL::BinaryExpr<SciPAL::Matrix<T, BW>, mult, SciPAL::UnaryExpr<SciPAL::Matrix<T, BW>, SciPAL::expr_transpose>>& expr)
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;

    const Mtx & A = expr.l;
    const Mtx & B = expr.r.l;

    Assert((this->n_rows() == A.n_rows())
           && (this->n_cols() == B.n_rows()),
           dealii::ExcMessage("Dimension mismatch"));

    Assert(A.n_cols()==B.n_cols(), dealii::ExcMessage("Input matrix dimensions do not match in operator+=(MtM)."));

    *this = A * SciPAL::transpose(B) + *this;

    return *this;
}

// @sect4{Operator: +=}
//!
//! Incremental, element-wise addition of two matrices.
//! @param other : Matrix to add. Must have the same dimensions.
//! A += B * adjoint(C)
template<typename T, typename BW>
SciPAL::Matrix<T, BW> &
SciPAL::Matrix<T, BW>::operator += (const SciPAL::BinaryExpr<SciPAL::Matrix<T, BW>, mult, SciPAL::UnaryExpr<SciPAL::Matrix<T, BW>, SciPAL::expr_adjoint>>& expr)
{
    typedef ::SciPAL::Matrix<T, BW> Mtx;

    const Mtx & A = expr.l;
    const Mtx & B = expr.r.l;

    Assert((this->n_rows() == A.n_rows())
           && (this->n_cols() == B.n_rows()),
           dealii::ExcMessage("Dimension mismatch"));

    Assert(A.n_cols()==B.n_cols(), dealii::ExcMessage("Input matrix dimensions do not match in operator+=(MaM)."));

    *this = A * SciPAL::adjoint(B) + *this;

    return *this;
}
#endif 


// @sect4{Operator: -=}
//!
//! Incremental, element-wise difference of two matrices.
//! @param other : Matrix to substract. Must have the same dimensions.
template<typename T, typename BW>
SciPAL::Matrix<T, BW> &
SciPAL::Matrix<T, BW>::operator -= (const Matrix<T, BW> & other)
{

    Assert((this->n_rows() == other.n_rows())
           && (this->n_cols() == other.n_cols()),
           dealii::ExcMessage("Dimension mismatch"));

    int n = this->size();
    T alpha = -1.;

    const T * const x = other.data();
    int incx = 1;

    T * y = this->data();
    int incy = 1;

    BW::axpy(n, alpha, x, incx, y, incy);

    return *this;
}

// @sect4{Operator: *=}
//!
//! Rescaling by a constant factor. Works for complex matrices as well.
//! @param scale : Factor
template<typename T, typename BW>
template<typename T2>
SciPAL::Matrix<T, BW> &
SciPAL::Matrix<T, BW>::operator *= (const T2 scale)
{
    int elem_dist = 1;

    int n = this->size();

    BW::scal(n, scale, &(this->data()[0]), elem_dist);

    return *this;
}


// @sect4{Function: l2_norm}
//! Consider the contents of the matrix as a vector and compute its L2-norm.
template<typename T, typename BW>
typename PrecisionTraits<T, BW::arch>::NumberType SciPAL::Matrix<T, BW>::l2_norm() const
{
    typename PrecisionTraits<T, BW::arch>::NumberType result = BW::nrm2(this->size(), this->data(), 1/*incx*/);

    return result;
}


// @sect4{Function: sum}
//!
//! sum over all entries.

template<typename T, typename BW>
T SciPAL::Matrix<T, BW>::sum() const
{
    return BW::asum(this->__n, &(this->data()[0]), 1);
}

// @sect4{Function: vmult}
//!
//! Matrix-vector multiplication.
//! This is compatible with the Krylov solvers of deal.II, i.e. you can use SciPAL matrices in deal.II's Krylov solver suite.
//! $dst = A \cdot src$
//! @param dst : Vector which stores the result.
//! @param src : Vector which the matrix is multiplied with.
template<typename T, typename BW>
template<typename VECTOR1, typename VECTOR2>
inline
void SciPAL::Matrix<T, BW>::vmult(VECTOR1& dst, const VECTOR2& src) const
{
    //! y = alpha*A*x + beta*y
    One<T> one;
    Zero<T> zero;
    T alpha = one();
    T beta  = zero();
    int leading_dim_A = this->n_rows();

    T * dst_ptr = dst.data();
    const T * src_ptr = src.data();

    Assert(src.size() == this->n_cols(),
           dealii::ExcMessage("Dimension mismatch"));

    BW::gemv('n', this->n_rows(), this->n_cols(), alpha, this->data(),
             leading_dim_A, src_ptr, 1, beta, dst_ptr, 1);
}

// @sect4{Function: Tvmult}
//!
//!
//! Matrix-vector multiplication using the transpose of the matrix.
//! This is compatible with the Krylov solvers of deal.II, i.e. you can use SciPAL matrices in deal.II's Krylov solver suite.
//! $dst = A^T \cdot src$
//! @param dst : Vector which stores the result.
//! @param src : Vector which the matrix is multiplied with.
template<typename T, typename BW>
template<typename VECTOR1, typename VECTOR2>
inline
void SciPAL::Matrix<T, BW>::Tvmult(VECTOR1& dst, const VECTOR2& src) const
{
    // In general gemv does this: y = alpha*A*x + beta*y
    T alpha = 1.;
    T beta = 0.;
    int leading_dim_A = this->n_rows();

    BW::gemv('t', this->n_rows(), this->n_cols(), alpha, this->data(),
             leading_dim_A, src.data(), 1, beta, dst.data(), 1);
}



// @sect4{Function: scaled_vmult_add_scaled}
//!
//! General matrix-vector multiplication.
//! $dst = \alpha A \cdot src + \beta dst$
//! @param alpha : scalar coefficient of the MVP.
//! @param dst : Vector which stores the result.
//! @param beta : scalar coefficient for the vector which is added to the result of the MVP.
//! @param src : Vector which the matrix is multiplied with.
template<typename T, typename BW>
void
SciPAL::Matrix<T, BW>::scaled_vmult_add_scaled(T alpha,
                                                 bool A_is_transpose,
                                                 SciPAL::Vector<T, BW>& dst /*y*/,
                                                 T beta,
                                                 const SciPAL::Vector<T, BW>& src /*x*/)
{
    // y = alpha*A*x + beta*y

    int leading_dim_A = this->n_rows(); // TODO: is this always the case?

    T * dst_ptr = dst.data();
    const T * src_ptr = src.data();


    Assert(src.size() == this->n_cols(),
           dealii::ExcMessage("Dimension mismatch"));

    BW::gemv((A_is_transpose ? 't' : 'n'), this->n_rows(), this->n_cols(), alpha, this->data(),
             leading_dim_A, src_ptr, 1, beta, dst_ptr, 1);
}




// @sect4{Function: mmult}
//!
//! matrix-matrix multiplication
//! $dst = *this \cdot src$
//! @param dst : Matrix which stores the result.
//! @param src : Matrix which this matrix is multiplied with.
template<typename T, typename BW>
void
SciPAL::Matrix<T, BW>::mmult(Matrix<T, BW>& dst, const Matrix<T, BW>& src) const
{
    dst.reinit(this->n_rows(), src.n_cols());

    this->scaled_mmult_add_scaled(dst, src, false, false);
}

//! matrix-matrix multiplication to produce a Submatrix from the product of two matrices.
//! $dst = *this \cdot src$
//! @param dst : Submatrix which stores the result.
//! @param src : Matrix which this matrix is multiplied with.
template<typename T, typename BW>
void
SciPAL::Matrix<T, BW>::mmult(SubMatrixView<T, BW>& dst, const Matrix<T, BW>& src) const
{


    this->scaled_mmult_add_scaled(dst, src, false, false);
}



// @sect4{Function: mTmult}
//!
//! Matrix-Matrix multiplication with the second matrix transposed.
//! $dst = this \cdot src^T$
//! @param dst : Matrix which stores the result.
//! @param src : Matrix which this matrix is multiplied with.
template<typename T, typename BW>
void
SciPAL::Matrix<T, BW>::mTmult(Matrix<T, BW>& dst, const Matrix<T, BW>& src) const
{
    dst.reinit(this->n_rows(), src.n_rows());

    this->scaled_mmult_add_scaled(dst, src, false, true);

}


// @sect4{Function: Tmmult}
//!
//! Matrix-Matrix multiplication with the first matrix transposed.
//! $dst = this \cdot src^T$
//! @param dst : Matrix which stores the result.
//! @param src : Matrix which this matrix is multiplied with.
template<typename T, typename BW>
void
SciPAL::Matrix<T, BW>::Tmmult(Matrix<T, BW>& dst, const Matrix<T, BW>& src) const
{
    dst.reinit(this->n_cols(), src.n_cols());

    this->scaled_mmult_add_scaled(dst, src, true, false);

}

//! General matrix-matrix product.
//! $dst = \alpha \cdot this \cdot src + \beta dst$
//! @param dst : Matrix which stores the result.
//! @param src : Matrix which this matrix is multiplied with.
//! @param transpose_A : whether this matrix is transposed
//! @param transpose_B : whether @p src is transposed.
//! @param alpha : scalar coefficient of the MMP.
//! @param beta : coefficient to rescale @dst with before the result of the MMP is added.
template<typename T, typename BW>
void
SciPAL::Matrix<T, BW>::scaled_mmult_add_scaled( Matrix<T, BW>& dst,
                                                  const Matrix<T, BW>& src,
                                                  char transpose_A,
                                                  char transpose_B,
                                                  T alpha, T beta) const
{
    // T alpha = 1; T beta = 0;

    int lda = this->leading_dim; //! this == A
    int ldb = src.leading_dim; /* == this->n_cols() !!! */ //! src == B
    int ldc = dst.leading_dim; //! dst == C

    BW::gemm(transpose_A,
             transpose_B,
             transpose_A != 'n'? this->n_cols() : this->n_rows(),
             /* cublas doc : m == n_rows of op(A), i.e. n_cols for A^T*/
             dst.n_cols(),
             /* cublas doc : n == n_cols of op(B), i.e. n_cols of C */
             transpose_A != 'n'? this->n_rows() : this->n_cols(),
             /* cublas doc : k == n_cols of op(A), i.e. n_rows of op(B) or n_rows for A^T */
             alpha,
             this->data(), lda,
             src.data(), ldb,
             beta,
             dst.data(), ldc);
}


// @sect4{Function: TmTmult}
//!
//! Matrix-Matrix multiplication with the first matrix transposed.
//! $dst = this^T \cdot src^T$
//! @param dst : Matrix which stores the result.
//! @param src : Matrix which this matrix is multiplied with.
template<typename T, typename BW>
void
SciPAL::Matrix<T, BW>::TmTmult(Matrix<T, BW>& dst, const Matrix<T, BW>& src) const
{
    dst.reinit(this->n_cols(), src.n_rows());

    this->scaled_mmult_add_scaled(dst, src, true, true);
}



// @sect4{Function: add_scaled_outer_product}
//!
//! Needed for Householder transformations. Creats a matrix from multiplying a column by a row vector.
//! $ this += \alpha x \cdot y^T $
//! @param x : column vector.
//! @param y : row vector.
template<typename T, typename BW>
template<typename VECTOR1, typename VECTOR2>
void
SciPAL::Matrix<T, BW>::add_scaled_outer_product(T alpha,
                                                  const VECTOR1& x,
                                                  const VECTOR2& y)
{

    int m = this->n_rows();
    int n = this->n_cols();
    int lda = this->n_rows();

#ifdef DEBUG
    if (x.size() != m) std::cout << "Dimension mismatch. "
                                    "Vector x.size() should be " << m << " but is " << x.size() << std::endl;
    if (y.size() != n) std::cout << "Dimension mismatch. "
                                    "Vector y.size() should be " << n << " but is " << y.size() << std::endl;
#endif
    Assert(x.size() == m, dealii::ExcMessage("Dimension mismatch"));
    Assert(y.size() == n, dealii::ExcMessage("Dimension mismatch"));

    int incx = x.stride;
    int incy = y.stride;

    BW::ger(m, n, alpha,
            x.data(), incx,
            y.data(), incy,
            this->data(), lda);
}

// @sect4{Operator: ()}
//!
//! Read access to an individual matrix element.
//! @param i : row index in C-style, i.e. starting at 0.
//! @param j : column index in C-style, i.e. starting at 0.
template <typename T, typename BW>
inline
T
SciPAL::Matrix<T, BW>::operator () (const unsigned int r,
                                      const unsigned int c) const
{
#ifdef DEBUG_MATRIX
    if (r > this->n_rows())
    {
        std::cerr << "Out of range: Row " << r << " is wanted out of a " << this->n_rows() << "x" << this->n_cols() << " matrix." << std::endl;
        std::cerr << "line :" << __LINE__ << ", Matrix<T,BW>" << std::endl;
        //print_expr_info(__PRETTY_FUNCTION__);
        std::exit(-1);
    }
    else if (c > this->n_cols())
    {
        std::cerr << "Out of range: Column " << c << " is wanted out of a " << this->n_rows() << "x" << this->n_cols() << " matrix." << std::endl;
        std::cerr << "line :" << __LINE__ << ", Matrix<T,BW>" << std::endl;
        //print_expr_info(__PRETTY_FUNCTION__);
        std::exit(-1);
    }
#endif
    int lead_dim = this->leading_dim;
    const T * tmp_d =  & this->data()[c*lead_dim+r];
    T entry;
    T * p_e = &entry;
    BW::GetMatrix(1, 1,
                  tmp_d, lead_dim, p_e, 1);

    return entry;
}

// @sect4{Operator: ()}
//!
//! Write access to an individual matrix element.
//! @param i : row index in C-style, i.e. starting at 0.
//! @param j : column index in C-style, i.e. starting at 0.
//! @param data : Data that will be set.
template <typename T, typename BW>
inline
void
SciPAL::Matrix<T, BW>::operator () (const unsigned int r,
                                      const unsigned int c, T data)
{
#ifdef DEBUG_MATRIX
    if (r > this->n_rows())
    {
        std::cerr << "Out of range: Row " << r << " should be written in of a " << this->n_rows() << "x" << this->n_cols() << " matrix." << std::endl;
        std::cerr << "line :" << __LINE__ << ", Matrix<T,BW>" << std::endl;
        //print_expr_info(__PRETTY_FUNCTION__);
        Assert(false, dealii::ExcMessage("see above"));
    }
    else if (c > this->n_cols())
    {
        std::cerr << "Out of range: Column " << c << " should be written in of a " << this->n_rows() << "x" << this->n_cols() << " matrix." << std::endl;
        std::cerr << "line :" << __LINE__ << ", Matrix<T,BW>" << std::endl;
        //print_expr_info(__PRETTY_FUNCTION__);
        Assert(false, dealii::ExcMessage("see above"));
    }
#endif
    int lead_dim = this->leading_dim;
    T * tmp_d =  & this->data()[c*lead_dim+r];
    T * p_e = &data;
    BW::SetMatrix(1, 1,
                  p_e, lead_dim, tmp_d, 1);

    return;
}

// @sect4{Function: Matrix::print}
//!
//! Dump a matrix to a stream. In case of device matrices their content gets copied back to the host before sendign it to the stream.
template<typename T, typename BW>
void
SciPAL::Matrix<T, BW>::print(std::ostream& Output) const
{
    //! T numerical_zero = 2.5e-16;

    Output << "#Matrix dims : " << this->n_rows() << " "
              << this->n_cols() << std::endl;

    int n_el = this->n_rows() * this->n_cols();
    T * tmp = new T[n_el];

    BW::GetMatrix(this->n_rows(), this->n_cols(),
                  this->data(), this->n_rows(), tmp, this->n_rows());

    for (uint r = 0; r < this->n_rows(); ++r)
    {
        for (uint c = 0; c < this->n_cols(); ++c)
            Output << //std::setprecision(15) << std::fixed << std::setw(15) <<
                       std::setprecision (20) << std::scientific <<
                         //!(std::abs(tmp[c*this->n_rows() + r])> numerical_zero
                         //!           ?
                         tmp[c*this->n_rows() + r]
                         //!                 : 0.)
                      << " ";
        Output /*<<";"*/ << std::endl;

    }
    delete [] tmp;
}

#endif
