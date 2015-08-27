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


#ifndef cublas_SubMatrixView_H
#define cublas_SubMatrixView_H

#include <lac/cublas_Vector.h>

#include <base/ForewardDeclarations.h>

namespace SciPAL {

// @sect3{Klasse: SubMatrixView}
//!
//! Diese Klasse dient dazu, nicht interessierende Teile einer Matrix zu maskieren.
template<typename T, typename BW>
class SubMatrixView
        :
        public SciPAL::Expr<SubMatrixView<T, BW> >,
        public SciPAL::Shape<T, BW, matrix>
{

public:

    typedef T value_type;

    typedef BW blas_wrapper_type;

    typedef SubMatrixView<T, BW> Type;

    typedef const Type& ConstHandle;

    typedef  SciPAL::ShapeData<T> DevType;

    typedef SciPAL::Shape<T, BW, matrix> MyShape;

    static const EType I_am = leafE;

    SubMatrixView(Matrix<T, BW> & src, int r_begin, int c_begin);

    SubMatrixView(Matrix<T, BW> & src, int r_begin, int r_end,
                  int c_begin, int c_end);

    ~SubMatrixView()
    {std::cout<<"Destructor MatrixView"<<std::endl;}

    void print() const;

    SubMatrixView & operator = (const Matrix<T, BW>& col);

    //  SubMatrixView<T, BW> & operator = (const SMSMmult<T, BW> & AB);
    template <typename X>
    SubMatrixView<T, BW> & operator = (const SciPAL::Expr<X> & expr);

//    SubMatrixView<T, BW> & operator += (const SMSMmult<T, BW> & AB);

//    SubMatrixView<T, BW> & operator += (const SMSMTmult<T, BW> & AB);

//    SubMatrixView<T, BW> & operator += (const SMTSMmult<T, BW> & AB);

//    SubMatrixView<T, BW> & operator -= (const SMSMmult<T, BW> & AB);

    const T * val() const { return this->MyShape::view_begin; }

    T * val() { return  this->MyShape::view_begin; }

    template<typename VECTOR1, typename VECTOR2>
    void vmult(VECTOR1& dst, const VECTOR2& src) const;


    template<typename VECTOR1, typename VECTOR2>
    void Tvmult(VECTOR1& dst, const VECTOR2& src) const;

    void add_scaled_outer_product(T alpha,
                                  const Vector<T, BW>& x,
                                  const Vector<T, BW>& y);

    template<typename T_src>
    void add_scaled_outer_product(T alpha,
                                  const SubVectorView<T, BW, T_src>& x,
                                  const SubVectorView<T, BW, T_src>& y);

    SubMatrixView() {}

private:
    SubMatrixView(const SubMatrixView<T, BW>& other) {}


    SubMatrixView & operator = (const SubMatrixView<T, BW>& col) {}



    //! deal.II's SmartPointer sind keine eigenstaendigen Zeiger,
    //! sondern dazu gedacht, auf schon bestehende Objekte zu zeigen
    //! und eine Exception auszuloesen, wenn das Objekt auf
    //! das sie zeigen vorzeitig zerstoert wird.
    dealii::SmartPointer<Matrix<T, BW> >__src;

public:
    // @sect4{Inline Funktionen: SubMatrixView}
    //!
    //! Diese Funktionen erm&ouml;glichen den Zugriff auf die Objektattribute.

    inline int leading_dim() const{
        return this->MyShape::leading_dim;

    }


    inline int r_begin() const {
        return this->MyShape::r_begin_active;
    }

    inline int c_begin() const {
        return this->MyShape::c_begin_active;
    }

    inline int r_end() const {
        return this->MyShape::r_end_active;
    }

    inline int c_end() const {
        return this->MyShape::c_end_active;
    }

    size_t n_rows() const { return this->MyShape::n_rows_active(); }

    size_t n_cols() const { return this->MyShape::n_cols_active(); }


    inline const Matrix<T, BW>& origin() const {
        return *__src;
    }

    // @sect4{Funktion: data}
    //!
    //! Direkter Zugriff auf den den device-Zeiger.
    //! Zurueckgegeben wird der Zeiger auf das erste Element
    //! des betrachteten Teils eines Vektors.
    const T* data() const
    {
        //! Indizes in C-Z&auml;hlung!!!
        return this->view_begin;
    }


    T* data()
    {
        //! Indizes in C-Z&auml;hlung!!!
        return this->view_begin;
    }


    void shift(size_t m_r, size_t m_c);
    void reset(size_t new_r_begin, size_t new_r_end, size_t new_c_begin, size_t new_c_end);
    void reinit(size_t new_r_begin, size_t new_r_end, size_t new_c_begin, size_t new_c_end);
};
}



// @sect4{Funktion: SubMatrixView::shift}
//! @param m_r:        move row
//! @param m_c:        move column
//!
template<typename T, typename BW>
void SciPAL::SubMatrixView<T, BW>::shift(size_t m_r, size_t  m_c) {

    reset(this->MyShape::r_begin_active + m_r, this->MyShape::r_end_active + m_r,
          this->MyShape::c_begin_active + m_c, this->MyShape::c_end_active + m_c);

}

// @sect4{Funktion: SubMatrixView::reset}
//! @param new_r_begin:      Anfang der neuen Zeile
//! @param new_r_end:        Ende der neuen Zeile
//! @param new_c_begin:      Anfang der neuen Spalte
//! @param new_c_end:        Ende der neuen Spalte
//!
template<typename T, typename BW>
void SciPAL::SubMatrixView<T, BW>::reset(size_t new_r_begin, size_t new_r_end,
                                         size_t new_c_begin, size_t new_c_end) {
#ifdef DEBUG
    Assert(new_r_begin >= 0, dealii::ExcMessage("View out of matrix bounds."));
    Assert(new_r_begin < new_r_end, dealii::ExcMessage("View out of matrix bounds."));
    Assert(new_r_end <= __src->n_rows(), dealii::ExcMessage("View out of matrix bounds."));

    Assert(new_c_begin >= 0, dealii::ExcMessage("View out of matrix bounds."));
    Assert(new_c_begin < new_c_end, dealii::ExcMessage("View out of matrix bounds."));
    Assert(new_c_end <= __src->n_cols(), dealii::ExcMessage("View out of matrix bounds."));
#else
    AssertThrow(new_r_begin >= 0, dealii::ExcMessage("View out of matrix bounds."));
    AssertThrow(new_r_begin < new_r_end, dealii::ExcMessage("View out of matrix bounds."));
    AssertThrow(new_r_end <= __src->n_rows(), dealii::ExcMessage("View out of matrix bounds."));

    AssertThrow(new_c_begin >= 0, dealii::ExcMessage("View out of matrix bounds."));
    AssertThrow(new_c_begin < new_c_end, dealii::ExcMessage("View out of matrix bounds."));
    AssertThrow(new_c_end <= __src->n_cols(), dealii::ExcMessage("View out of matrix bounds."));
#endif

    this->MyShape::reinit_attr(new_r_begin, new_r_end,
                          new_c_begin, new_c_end, 1);
}


// @sect4{Funktion: SubMatrixView::reinit}
//! @param new_r_begin:      Anfang der neuen Zeile
//! @param new_r_end:        Ende der neuen Zeile
//! @param new_c_begin:      Anfang der neuen Spalte
//! @param new_c_end:        Ende der neuen Spalte
//!
template<typename T, typename BW>
void SciPAL::SubMatrixView<T, BW>::reinit(size_t new_r_begin, size_t new_r_end,
                                          size_t new_c_begin, size_t new_c_end) {
    reset(new_r_begin, new_r_end, new_c_begin, new_c_end);
}



//!
// @sect4{Konstruktoren}
//! @param src:       source Matrix
//! @param r_begin:   Anfang der Zeile
//! @param c_begin:   Anfang der Spalte
template<typename T, typename BW>
SciPAL::SubMatrixView<T, BW>::SubMatrixView(Matrix<T, BW> & src,
                                            int r_begin, int c_begin)
    :
      MyShape(src.shape(),
              r_begin, src.n_rows(),
              c_begin, src.n_cols()),
      __src(&src)
{
    Assert ((r_begin >= 0) && (r_begin < src.n_rows()),
            dealii::ExcIndexRange (r_begin, 0, src.n_rows()));
    Assert ((c_begin >= 0) && (c_begin < src.n_cols()),
            dealii::ExcIndexRange (c_begin, 0, src.n_cols()));

}


//! @param src:        source Matrix
//! @param r_begin:    Anfang der Zeile
//! @param r_end:      Ende der Zeile
//! @param c_begin:    Anfang der Spalte
//! @param c_end:      Ende der Spalte
template<typename T, typename BW>
SciPAL::SubMatrixView<T, BW>::SubMatrixView(Matrix<T, BW> & src,
                                            int r_begin, int r_end,
                                            int c_begin, int c_end)
    :
      MyShape(src.shape(),
              r_begin, r_end,
              c_begin, c_end),
      __src(&src)
{
    //! Pruefe im DEBUG-Modus Zulaessigkeit der Indexgrenzen.
    Assert ((r_begin >= 0) && (r_begin < src.n_rows()),
            dealii::ExcIndexRange (r_begin, 0, src.n_rows()));
    Assert ((c_begin >= 0) && (c_begin < src.n_cols()),
            dealii::ExcIndexRange (c_begin, 0, src.n_cols()));

    Assert ((r_end > r_begin) && (r_end <= src.n_rows()),
            dealii::ExcIndexRange (r_end, 1, src.n_rows()+1));
    Assert ((c_end > c_begin) && (c_end <= src.n_cols()),
            dealii::ExcIndexRange (c_end, 1, src.n_cols()+1));

}



// @sect4{Operator: =}
//! @param col:    Matrix
//!
template<typename T, typename BW>
SciPAL::SubMatrixView<T, BW> &
SciPAL::SubMatrixView<T, BW>::operator = (const Matrix<T, BW>& col)
{
    Assert(this->size() <= col.size(),
           dealii::ExcMessage("Dimension mismatch"));

    int incx = 1;
    int incy = 1;
    int n_rows_2_copy = this->n_rows_active();
    int n_cols_2_copy = this->n_cols_active();

    std::cout << "n_rows_2_copy : " << n_rows_2_copy << ", "
              << "n_cols_2_copy : " << n_cols_2_copy << std::endl;

    //! Kopiere spaltenweise.
    for (int c = this->c_begin_active; c < this->c_end_active; ++c)
        BW::copy(n_rows_2_copy, (col.data() + c*(col.n_rows()) + this->r_begin()),
                 col.stride,
                 this->__src->data() + this->r_begin() + c*this->leading_dim(),
                 this->stride);

    return *this;

}


//namespace SciPAL {

//namespace MatrixOperations {

//template <typename T, typename BW>
//static void apply(::SciPAL::SubMatrixView<T, BW> &result,
//                  const ::SciPAL::SMSMmult<T, BW> & // typename BlasMatExp<T, BW>::scaledM&
//                  AB)
//{
//#ifdef DEBUG
//    std::cout << " MY new and shiny SMV prod" << std::endl;
//#endif
//    const T * A = AB.l.data(); // SubMatrixView::val() automatically computes the begin of the view in memory  //matrix().array().data()[AB.l.__view_begin]) ;
//    const T * B = AB.r.data(); // &(AB.r.__src->val()[AB.r.__view_begin]) ;
//    T * C = result.data(); // &( result./*this->*/__src->val()[ result./*this->*/__view_begin]) ;

//    T alpha = +1;
//    T beta  = 0.;

//    int lda = AB.l.leading_dim();
//    int ldb = AB.r.leading_dim();
//    int ldc =  result./*this->*/leading_dim();

//    int m = AB.l.r_end() - AB.l.r_begin() ;
//    int n = AB.r.c_end() - AB.r.c_begin ();
//    int k = AB.l.c_end() - AB.l.c_begin ();

//    BW::gemm('n', 'n',
//             m, n, k,
//             alpha,
//             A, lda,
//             B, ldb,
//             beta,
//             C, ldc);

//}
//}
//}


// @sect4{Operator: =}
//! @param AB:        SMSMmult Objekt
//!
template<typename T, typename BW>
template<typename X>
SciPAL::SubMatrixView<T, BW> &
SciPAL::SubMatrixView<T, BW>::operator = (const SciPAL::Expr<X> &e)
{
#ifdef DEBUG
    std::cout << "line :" << __LINE__ << ", SubMatrixView<T,BW>" << std::endl;
    print_expr_info(__PRETTY_FUNCTION__);
#endif
    SciPAL::LAOOperations::apply(*this,  ~e);
    return *this;
}


// @sect4{Operator: +=}
//! @param AB:        SMSMmult Objekt
//!
/*template<typename T, typename BW>
SciPAL::SubMatrixView<T, BW> &
SciPAL::SubMatrixView<T, BW>::operator += (const SMSMmult<T, BW> & AB)
{
    const T * A = &(AB.l.__src->val()[AB.l.__view_begin]) ;
    const T * B = &(AB.r.__src->val()[AB.r.__view_begin]) ;
    T * C = &(this->__src->val()[this->__view_begin]) ;

    T alpha = +1;
    T beta = 1;

    int lda = AB.l.leading_dim;
    int ldb = AB.r.leading_dim;
    int ldc = this->leading_dim;

    int m = AB.l.__r_end - AB.l.__r_begin ;
    int n = AB.r.n_cols_active();
    int k = AB.l.n_cols_active();

    BW::gemm('n', 'n',
             m, n, k,
             alpha,
             A, lda,
             B, ldb,
             beta,
             C, ldc);

    return *this;
}*/



// @sect4{Operator: +=}
//! @param AB:        SMSMTmult Objekt
//!
/*template<typename T, typename BW>
SciPAL::SubMatrixView<T, BW> &
SciPAL::SubMatrixView<T, BW>::operator += (const SMSMTmult<T, BW> & AB)
{
    //!std::cout<<__PRETTY_FUNCTION__<<std::endl;

    const T * A = &(AB.l.__src->val()[AB.l.__view_begin]) ;
    const T * B = &(AB.r.A.__src->val()[AB.r.A.__view_begin]) ;
    T * C = &(this->__src->val()[this->__view_begin]) ;

    T alpha = +1;
    T beta = 1;

    int lda = AB.l.leading_dim;
    int ldb = AB.r.A.leading_dim;
    int ldc = this->leading_dim;

    int m = AB.l.__r_end - AB.l.__r_begin ;
    int n = AB.r.n_cols_active();
    int k = AB.l.n_cols_active();

    BW::gemm('n', 't',
             m, n, k,
             alpha,
             A, lda,
             B, ldb,
             beta,
             C, ldc);

    return *this;
}
*/

// @sect4{Operator: +=}
//! @param AB:        SMTSMmult Objekt
//!
/*template<typename T, typename BW>
SciPAL::SubMatrixView<T, BW> &
SciPAL::SubMatrixView<T, BW>::operator += (const SMTSMmult<T, BW> & AB)
{
    //! std::cout<<__PRETTY_FUNCTION__<<std::endl;
    const T * A = &(AB.l.A.__src->val()[AB.l.A.__view_begin]) ;
    const T * B = &(AB.r.__src->val()[AB.r.__view_begin]) ;
    T * C = &(this->__src->val()[this->__view_begin]) ;

    T alpha = +1;
    T beta = 1;

    int lda = AB.l.A.leading_dim;
    int ldb = AB.r.leading_dim;
    int ldc = this->leading_dim;

    int m = AB.l.A.n_cols_active();
    int n = AB.r.n_cols_active() ;
    int k = AB.l.A.n_cols_active();

    std::cout<<"m "<< m
            <<" n "<< n
            <<" k "<<k
            <<" lda "<<lda
            <<" ldb "<<ldb
            <<" ldc "<<ldc
            << std::endl;



    BW::gemm('t', 'n',
             m, n, k,
             alpha,
             A, lda,
             B, ldb,
             beta,
             C, ldc);

    return *this;
}
*/



// @sect4{Operator: -=}
//! @param AB:        SMSMult Objekt
//!
/*template<typename T, typename BW>
SciPAL::SubMatrixView<T, BW> &
SciPAL::SubMatrixView<T, BW>::operator -= (const SMSMmult<T, BW> & AB)
{
    const T * A = &(AB.l.__src->val()[AB.l.__view_begin]) ;
    const T * B = &(AB.r.__src->val()[AB.r.__view_begin]) ;
    T * C = &(this->__src->val()[this->__view_begin]) ;

    T alpha = -1;
    T beta = 1;

    int lda = AB.l.leading_dim;
    int ldb = AB.r.leading_dim;
    int ldc = this->leading_dim;

    int m = AB.l.n_rows_active() ;
    int n = AB.r.n_cols_active();
    int k = AB.l.n_cols_active();

    BW::gemm('n', 'n',
             m, n, k,
             alpha,
             A, lda,
             B, ldb,
             beta,
             C, ldc);

    return *this;
}
*/

// @sect4{Funktion: SubMatrixView::vmult}
//! @param src:        source Vektor
//! @param dst:        Zielvektor
//!
template<typename T, typename BW>
template<typename VECTOR1, typename VECTOR2>
void
SciPAL::SubMatrixView<T, BW>::vmult(VECTOR1& dst, const VECTOR2& src) const
{
    //! y = alpha*A*x + beta*y
    T alpha = 1.;
    T beta = 0.;
    int n_rows = this->n_rows_active();
    int n_cols = this->n_cols_active();

    Assert(src.size() >= n_cols, dealii::ExcMessage("Dimension mismatch"));
    Assert(dst.size() >= n_rows, dealii::ExcMessage("Dimension mismatch"));

    const int dst_val_begin = (VECTOR1::is_vector_view ? 0 : this->__r_begin );
    T *dst_val_ptr = dst.data() + dst_val_begin;


    const int src_val_begin = (VECTOR2::is_vector_view ? 0 : this->c_begin_active );
    const T * const src_val_ptr = src.data() + src_val_begin;

    BW::gemv('n', n_rows, n_cols, alpha, this->view_begin,
             this->leading_dim, src_val_ptr, 1, beta, dst_val_ptr, 1);
}



// @sect4{Funktion: SubMatrixView::Tvmult}
//! @param src:        source Vektor
//! @param dst:        Zielvektor
//!
template<typename T, typename BW>
template<typename VECTOR1, typename VECTOR2>
void
SciPAL::SubMatrixView<T, BW>::Tvmult(VECTOR1& dst, const VECTOR2& src) const
{
    //! y = alpha*A*x + beta*y
    T alpha = 1.;
    T beta = 0.;
    int n_rows = this->n_rows_active();
    int n_cols = this->n_cols_active();

    Assert(src.size() >= n_cols, dealii::ExcMessage("Dimension mismatch"));
    Assert(dst.size() >= n_rows, dealii::ExcMessage("Dimension mismatch"));

    const int dst_val_begin = (VECTOR1::is_vector_view ? 0 : this->c_begin_active );
    T *dst_val_ptr              = dst.data() + dst_val_begin;

    const int src_val_begin = (VECTOR2::is_vector_view ? 0 :this->__r_begin );
    const T * const src_val_ptr = src.data() + src_val_begin;

    //! AssertThrow(false, dealii::ExcNotImplemented() );

    BW::gemv('t', n_rows, n_cols, alpha, this->view_begin,
             this->leading_dim, src_val_ptr, 1, beta, dst_val_ptr, 1);
}



// @sect4{Funktion: add_scaled_outer_product}
//!
//! Wird fuer Householder-Transformationen gebraucht.
//! $ this += \alpha x \cdot y^T $
//! @param alpha:        source Vektor
//! @param x:            Vektor x
//! @param y:            Vektor y
//!
template<typename T, typename BW>
void
SciPAL::SubMatrixView<T, BW>::add_scaled_outer_product(T alpha,
                                                       const Vector<T, BW>& x,
                                                       const Vector<T, BW>& y)
{
    int m = this->n_rows();
    int n = this->n_cols();


    size_t lda = this->leading_dim();


    //! TO DO : correct ASSERTIONS
    //! Assert(x.size() == m, dealii::ExcMessage("Dimension mismatch"));
    //! Assert(y.size() == n, dealii::ExcMessage("Dimension mismatch"));

    int incx = 1;
    int incy = 1;

    //! @p x wird als Spaltenvektor betrachtet
//    const int x_val_begin = this->r_begin_active ;
    const T * const x_val_ptr = x.data();// + x_val_begin;

    //! @p y wird als Zeilenvektor betrachtet
//    const int y_val_begin = this->c_begin_active;
    const T * const y_val_ptr = y.data();// + y_val_begin;



    BW::ger(m, n, alpha,
            x_val_ptr, incx,
            y_val_ptr, incy,
            this->data(),
            lda);
}

template<typename T, typename BW>
template<typename T_src>
void
SciPAL::SubMatrixView<T, BW>::add_scaled_outer_product(T alpha,
                                                       const SubVectorView<T, BW, T_src>& x,
                                                       const SubVectorView<T, BW, T_src>& y)
{
    int m = this->n_rows();
    int n = this->n_cols();


    size_t lda = this->leading_dim();


    //! TO DO : correct ASSERTIONS
    //! Assert(x.size() == m, dealii::ExcMessage("Dimension mismatch"));
    //! Assert(y.size() == n, dealii::ExcMessage("Dimension mismatch"));

    int incx = 1;
    int incy = 1;

    //! @p x wird als Spaltenvektor betrachtet
//    const int x_val_begin = this->r_begin_active ;
    const T * const x_val_ptr = x.data();// + x_val_begin;

    //! @p y wird als Zeilenvektor betrachtet
//    const int y_val_begin = this->c_begin_active;
    const T * const y_val_ptr = y.data();// + y_val_begin;



    BW::ger(m, n, alpha,
            x_val_ptr, incx,
            y_val_ptr, incy,
            this->data(),
            lda);
}



// @sect4{Funktion: SubVectorView::print}
//!
template<typename T, typename BW>
void
SciPAL::SubMatrixView<T, BW>::print() const
{

    int n_rows_2_copy = this->n_rows_active();
    int n_cols_2_copy = this->n_cols_active();


    int n_el = this->size();
    T * tmp = new T[n_el];

    int lda = this->__src->leading_dim;
    int ldb = n_rows_2_copy;
    BW::GetMatrix(n_rows_2_copy, n_cols_2_copy,
                  (this->view_begin), lda,
                  tmp, ldb);

    for (int r = 0; r < n_rows_2_copy; ++r)
    {
        for (int c = 0; c < n_cols_2_copy; ++c)
            std::cout << std::setprecision(4)
                      << std::setw(15) << //! (std::abs(tmp[c*n_rows_2_copy + r])> numerical_zero
                         //! ?
                         tmp[c*n_rows_2_copy + r]
                         //!   : 0.)
                      << " ";
        std::cout << std::endl;
    }

    delete [] tmp;
}






#endif
