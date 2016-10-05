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


#ifndef cublas_SubVectorView_H
#define cublas_SubVectorView_H

#include <lac/Shape.h>
#include <lac/Expr.h>

#include <base/Zero_One_Traits.h>

#include <base/ForewardDeclarations.h>

namespace SciPAL {


// @sect3{Klasse: VectorView}
//!
//! Die View-Klassen dienen dazu, Teile von Vektoren oder Matrizen als
//! selbst&auml;ndige Objekte darzustellen, so dass man beispielsweise eine
//! Spalte einer Matrix als Quelle in einem Matrix-Vektorprodukt
//! nutzen kann.
template<typename T, typename BW, typename T_src>
class SubVectorView
        :
        public SciPAL::Expr<SubVectorView<T, BW, T_src> >,
        public SciPAL::Shape<T, BW, vector>
{

public:
    typedef BW blas_wrapper_type;

    typedef SciPAL::Shape<T, BW, vector> MyShape;

    friend class Vector<T, BW>;

    friend class Matrix<T, BW>;

    friend class SubMatrixView<T, BW>;

    template<typename, typename, typename> friend class  ColVectorView;
    template<typename, typename, typename> friend class  RowVectorView;

    typedef T value_type;

    typedef SubVectorView<T, BW, T_src> Type;

    typedef const Type& ConstHandle;

    typedef  SciPAL::ShapeData<T> DevType;


    //! Fuer Optimierungen der Matrix-Vektor-Operationen zwischen
    //! (Teil-)Matrizen und (Teil-)Vektoren muss bekannt sein, was fuer ein
    //! Typ die jeweiligen Objekte sind. Da statische Konstanten zur Compile-Zeit ausgewertet werden
    //! koennen dadurch bedingte Ausdruecke fuer die Komilierung geschrieben werden.
    //! Siehe beispielsweise SubMatrixView::vmult().
    static const bool is_vector_view = true;

    SubVectorView( T_src & src, int r_begin, int r_end, int c = 0);
    SubVectorView( T_src & src, int c_begin, int c_end, int r, bool iscol);

    ~SubVectorView()
    {
#ifdef DEBUG
        std::cout<<"Peng"<<std::endl;
#endif
    }

    typename PrecisionTraits<T, BW::arch>::NumberType l2_norm() const;

    template<typename VECTOR>
    T dot(const VECTOR & other) const;

    void print() const;


    template <typename T2>
    SubVectorView<T, BW, T_src> & operator = (const std::vector<std::complex<T2> > & other);

    unsigned int r_begin() const { return this->r_begin_active; }

    unsigned int c_begin() const { return this->c_begin_active; }

    unsigned int size() const { return this->n_elements_active; }

    const T* data() const;

    T* data();

    template<typename T2_src>
    SubVectorView & operator += (const SubVectorView <T, BW, T2_src> &other);
    SubVectorView & operator += (const SciPAL::BinaryExpr<SciPAL::Literal<T, BW>, SciPAL::mult, Vector<T, BW>>& expr);

    SubVectorView & operator -= (const SubVectorView <T, BW, T_src> &other);
    SubVectorView &operator -= (const Vector<T, BW>& col);

    SubVectorView & operator *= (const T alpha);
    SubVectorView & operator /= (const T alpha);
    SubVectorView & operator /=
    (const std::complex<typename PrecisionTraits<T, BW::arch>::NumberType> alpha);
    T operator () (int k) const;

    //! Copy the elements of other to this. This changes the elements of the
    //! source object of this.
    template <typename T_src2>
    SubVectorView & operator = (const SubVectorView<T, BW, T_src2> & other)
    {
        Assert(this->n_elements_active == other.n_elements_active,
               dealii::ExcMessage("Cannot copy subarrays of different lengths"));


        //! Elementweise Kopie des arrays.
        int inc_src = other.stride;
        int inc_dst = this->stride;

        BW::copy(this->n_elements_active, other.data(), inc_src,
                 this->data(), inc_dst);

        return *this;
    }

    SubVectorView & operator = (const SubVectorView& other)
    {
        Assert(this->n_elements_active == other.n_elements_active,
               dealii::ExcMessage("Cannot copy subarrays of different lengths"));


        //! Elementweise Kopie des arrays.
        int inc_src = other.stride;
        int inc_dst = this->stride;

        BW::copy(this->n_elements_active, other.data(), inc_src,
                 this->data(), inc_dst);

        return *this;
    }

    template<typename T_src2>
    void  copy(const SubVectorView<T, BW, T_src2> & other)
    {
        //        Base & self = *this;
        //        typedef typename ColVectorView<T, BW, T_src2>::Base Base2;
        //        const Base2 &src = other;
        //        self = src;

        //        return *this;

        Assert(this->n_elements_active == other.n_elements_active,
               dealii::ExcMessage("Cannot copy subarrays of different lengths"));

        int inc_src = other.stride;
        int inc_dst = this->stride;

        BW::copy(this->n_elements_active, other.data(), inc_src,
                 this->data(), inc_dst);
    }

    // @sect4{Operator: =}
    //! Applies expressions to view
    //!
    template<typename X>
    Type & operator = (const SciPAL::Expr<X> &e)
    {
#ifdef DEBUG_SUBVECTOR_VIEW
        std::cout << "line :" << __LINE__ << ", SubVectorView<T,BW>" << std::endl;
        print_expr_info(__PRETTY_FUNCTION__);
#endif

        SciPAL::LAOOperations::apply(*this,  ~e);
        return *this;
    }

    //

    // @sect4{Funktion: get_amax}
    //!
    //! returns the index of the absolutely largest value
    unsigned int
    get_amax()
    {
        Assert(false, dealii::ExcMessage("This function has been moved to 	algorithms.hh. and has been renamed to index_of_largest_element"))
    }

    // @sect4{Funktion: sadd}
    //!
    //! \f{eqnarray*} this & = \alpha \cdot other + this \f}
    //! @param alpha : Faktor
    //! @param other : zweiter Vector
    template<typename T_src2>
    void
    sadd (T alpha, const SubVectorView<T, BW,T_src2> & other)
    {
        Assert(this->n_elements_active == other.n_elements_active,
               dealii::ExcMessage("Cannot copy subarrays of different lengths in sadd"));

        int incx = 1;
        BW::scal(this->size(), alpha, this->data(), incx);

        //! cublasSaxpy
        int incy = 1;
        BW::axpy(this->size(), 1., other.data(), 1, this->data(), 1);
    }

private:


    //! deal.II's SmartPointer sind keine eigenst&auml;ndigen Zeiger,
    //! sondern dazu gedacht, auf schon bestehende Objekte zu zeigen
    //! und eine Exception auszul&ouml;sen, wenn das Objekt auf
    //! das sie zeigen vorzeitig zerst&ouml;rt wird.
    dealii::SmartPointer<T_src> __src;

protected:
    bool _is_col;
    SubVectorView() {}
    SubVectorView( const SubVectorView<T, BW, T_src> & other)
        :
          MyShape(other.__src->shape(),
                  other.r_begin_active, other.r_end_active,
                  other.c_begin_active, other.c_end_active),
          __src(other.__src),
          _is_col(true)
    {}
};



// @sect3{Klasse: ColVectorView}
//!
//! Ansicht auf einen Teil einer Spalte einer Matrix oder eines Spaltenvektors.
template<typename T, typename BW, typename T_src>
class ColVectorView
        :
        public SubVectorView<T, BW, T_src >
{

public:
    typedef BW blas_wrapper_type;
    typedef SubVectorView<T, BW, T_src> Base;
    typedef Shape<T, BW, vector> MyShape;
    ColVectorView(){};

    ColVectorView( T_src & src,
                   int c, int r_begin = 0 , int r_end = -1)
        : Base(src, r_begin, (r_end<0)?src.n_rows():r_end, c){}

    ColVectorView( const SubVectorView<T, BW, T_src>& other)
        : Base( other ) {}
    // @sect4{Operator: =}
    //! Applies expressions to view
    //!
    template<typename X>
    SciPAL::ColVectorView<T, BW, T_src> & operator = (const SciPAL::Expr<X> &e)
    {
#ifdef DEBUG_SUBVECTOR_VIEW
        std::cout << "line :" << __LINE__ << ", SubVectorView<T,BW>" << std::endl;
        print_expr_info(__PRETTY_FUNCTION__);
#endif

        SciPAL::LAOOperations::apply(*this,  ~e);
        return *this;
    }

    template<typename T2_src>
    ColVectorView<T, BW, T_src> & operator += (const SubVectorView<T, BW, T2_src> &other)
    {
        Base & self = *this;
        const typename ColVectorView<T, BW, T2_src>::Base & o = other;
        self += o;

        return *this;
    }

    template<typename T2_src>
    ColVectorView<T, BW, T_src> & operator += (const ColVectorView<T, BW, T2_src> &other)
    {
        Base & self = *this;
        const typename ColVectorView<T, BW, T2_src>::Base & o = other;
        self += o;

        return *this;
    }

    // @sect4{Operator: +=}
    //!
    //! Overloaded operator to add an expression "literal * vector" to a colvectorview
    //! col_vector_view_a += literal * vector_b
    //! @param other: expression "literal * vector"
    ColVectorView<T, BW, T_src> & operator +=(const SciPAL::BinaryExpr<SciPAL::Literal<T, BW>, SciPAL::mult, Vector<T, BW>>& expr)
    {
        Base & self = *this;
        const SciPAL::Vector<T, BW>& o = expr.r;
        T alpha = expr.l;

        self += alpha * o;

        return *this;
    }
     template<typename T2>
     SciPAL::ColVectorView<T, BW, T_src> & operator *= (const T2 scale);

    // Resets the column information of the vector view

    void reset(int c, int r_begin = 0 , int r_end = -1)
    {
        // Reset data_ptr in case the source has been reinitialized
        this->data_ptr = this->__src->data_ptr;
        this->leading_dim = this->__src->leading_dim;
        this->MyShape::reinit_attr(r_begin, (r_end<0)?this->r_end_active:r_end,
                                   c, c+1);
    }

    //somehow different behavior
    //    void reset(int r_begin, int c=0)
    //    {
    //        this->MyShape::reinit_attr(r_begin, this->r_end_active,
    //                              c, c+1);
    //    }



    template<typename T_src2>
    ColVectorView & operator = (const ColVectorView<T, BW, T_src2> & other)
    {
        Base & self = *this;
        typedef typename ColVectorView<T, BW, T_src2>::Base Base2;
        const Base2 &src = other;
        self = src;

        return *this;
    }
};//end class


// @sect3{Klasse: RowVectorView}
//!
//! Ansicht auf einen Teil einer Zeile einer Matrix oder eines Zeilenvektors.
template<typename T, typename BW, typename T_src>
class RowVectorView : public SubVectorView<T, BW, T_src> {

public:
    typedef SubVectorView<T, BW, T_src> Base;
    typedef BW blas_wrapper_type;

    typedef SciPAL::Shape<T, BW, vector> MyShape;
    RowVectorView(){};

    RowVectorView( T_src & src,
                   int r, int c_begin = 0, int c_end = -1)
        : Base(src, c_begin, (c_end<0)?src.n_cols():c_end, r, false)
    {}

    // Resets the row information of the vector view
    void reset(int r, int c_begin = 0, int c_end = -1)
    {
        // Reset data_ptr in case the source has been reinitialized
        this->data_ptr = this->__src->data_ptr;
        this->leading_dim = this->__src->leading_dim;
        this->MyShape::reinit_attr(r, r+1,
                                   c_begin, (c_end<0)?this->c_end_active:c_end,this->__src->leading_dim);
    }

    // @sect4{Operator: =}
    //! Applies expressions to view
    //!
    template<typename X>
    SciPAL::RowVectorView<T, BW, T_src> & operator = (const SciPAL::Expr<X> &e)
    {
#ifdef DEBUG_SUBVECTOR_VIEW
        std::cout << "line :" << __LINE__ << ", RowVectorView<T,BW>" << std::endl;
        print_expr_info(__PRETTY_FUNCTION__);
#endif
        SciPAL::LAOOperations::apply(*this,  ~e);
        return *this;
    }


    template<typename T_src2>
    RowVectorView & operator = (const RowVectorView<T, BW, T_src2> & other)
    {
        Base & self = *this;
        typedef typename RowVectorView<T, BW, T_src2>::Base Base2;
        const Base2 &src = other;
        self = src;

        return *this;
    }

    template <typename T_src2>
    RowVectorView & operator = (const SubVectorView<T, BW, T_src2> & other)
    {
        Assert(this->n_elements_active == other.n_elements_active,
               dealii::ExcMessage("Cannot copy subarrays of different lengths"));

        //! Elementweise Kopie des arrays.
        int inc_src = other.stride;
        int inc_dst = this->stride;

        BW::copy(this->n_elements_active, other.data(), inc_src,
                 this->data(), inc_dst);
        this->_is_col = false;
        return *this;
    }

};
}


// @sect3{SciPAL::VectorView Methoden}
//!

// @sect4{Konstruktor}
//!
//! Legt einen neuen VektorView an
//! @param src:      Urspr&uuml;ngliche Matrix oder Vektor
//! @param r_begin:  Zeilenindex bei dem der Vektor anf&auml;ngt
//! @param c:        kann sowohl Spalte als auch Zeile sein
template<typename T, typename BW, typename T_src>
SciPAL::SubVectorView<T, BW, T_src >::SubVectorView(T_src & src,
                                                    int r_begin, int r_end, int c)
    :
      MyShape(src.shape(),
              r_begin, r_end,
              c, c+1),
      __src(&src),
      _is_col(true)
{}
// @sect4{Konstruktor}
//!
//! Legt einen neuen VektorView an
//! @param src:      Urspr&uuml;ngliche Matrix oder Vektor
//! @param c_begin:  Spaltenindex bei dem der Vektor anf&auml;ngt
//! @param c_end:    Spaltenindex bei dem der Vektor anf&auml;ngt
//! @param r:        Zeilenindex
template<typename T, typename BW, typename T_src>
SciPAL::SubVectorView<T, BW, T_src >::SubVectorView(T_src & src,
                                                    int c_begin, int c_end, int r, bool iscol)
    :
      MyShape(src.shape(),
              r, r+1,
              c_begin, c_end),
      __src(&src),
      _is_col(iscol)
{
    this->n_elements_active  = c_end - c_begin;
    this->stride = src.n_rows();
}

//!@sect4{Funktion: l2_norm}
//!
//! L2-Norm des betrachteten Vektorteils.
template<typename T, typename BW, typename T_src>
typename PrecisionTraits<T, BW::arch>::NumberType
SciPAL::SubVectorView<T, BW, T_src>::
l2_norm() const
{
    typename PrecisionTraits<T, BW::arch>::NumberType result
            = BW::nrm2(this->size(), this->data(), this->stride);
    return result;
}

// @sect4{Funktion: data}
//!
//! Direkter Zugriff auf den den device-Zeiger.
//! Zurueckgegeben wird der Zeiger auf das erste Element
//! des betrachteten Teils eines Vektors.
template<typename T, typename BW, typename T_src>
const T *
SciPAL::SubVectorView<T, BW, T_src>::data() const
{
    //! Indizes in C-Z&auml;hlung!!!
    return this->view_begin;
}


template<typename T, typename BW, typename T_src>
T *
SciPAL::SubVectorView<T, BW, T_src>::data()
{
    //! Indizes in C-Z&auml;hlung!!!
    return this->view_begin;
}





// @sect4{Funktion: dot}
//!
//! Skalarprodukt von zwei Teilvektoren.
//! @param other: 2. Vektor
template<typename T, typename BW, typename T_src>
template<typename VECTOR>
T
SciPAL::SubVectorView<T, BW, T_src>::dot(const VECTOR & other) const
{

    Assert(this->n_elements_active == other.size(),
           dealii::ExcMessage("Dimension mismatch"));

    int incx = this->stride;
    int incy = other.stride;
    T result = BW::dot(this->n_elements_active,
                       this->data(), incx, other.data(), incy);

    return result;
}


// @sect4{Operator: +=}
//!
//! &Uuml;berladener Operator, mit dem man zwei Vektoren addieren kann
//! 1. Vektor += 2. Vektor
//! @param other: 2. Vektor
template<typename T, typename BW, typename T_src>
template<typename T2_src>
SciPAL::SubVectorView<T, BW, T_src> &
SciPAL::SubVectorView<T, BW, T_src> ::operator +=(const SciPAL::SubVectorView <T, BW, T2_src> &other)
{
    Assert(this->n_elements_active == other.size(),
           dealii::ExcMessage("Dimension mismatch"));

    One<T> one;

    int incx = other.stride;
    int incy = this->stride;

    BW::axpy(this->n_elements_active, one(), other.data(), incx, this->data(),
             incy);

    return *this;
}

// @sect4{Operator: +=}
//!
//! Overloaded operator to add an expression "literal * vector" to a subvectorview
//! sub_vector_view_a += literal * vector_b
//! @param other: expression "literal * vector"
template<typename T, typename BW, typename T_src>
SciPAL::SubVectorView<T, BW, T_src> &
SciPAL::SubVectorView<T, BW, T_src> ::operator +=(const SciPAL::BinaryExpr<SciPAL::Literal<T, BW>, SciPAL::mult, Vector<T, BW>>& expr)
{
    const Literal<T, BW>& a = expr.l;
    const Vector<T, BW>& A = expr.r;

    Assert(this->n_elements_active == A.n_elements_active,
           dealii::ExcMessage("Dimension mismatch"));

    int incx = A.stride;
    int incy = this->stride;

    BW::axpy(this->n_elements_active, a, A.data(), incx, this->data(),
             incy);

    return *this;
}

// @sect4{Operator: -=}
//!
//! &Uuml;berladener Operator, mit dem man zwei Vektoren subtrahieren kann
//! 1. Vektor -= 2. Vektor
//! @param other: 2. Vektor
template<typename T, typename BW, typename T_src>
SciPAL::SubVectorView<T, BW, T_src> &
SciPAL::SubVectorView<T, BW, T_src> ::operator -=(const SciPAL::SubVectorView <T, BW, T_src> &other)
{
    Assert(this->n_elements_active == other.size(),
           dealii::ExcMessage("Dimension mismatch"));

    One<T> one;

    int incx = other.stride;
    int incy = this->stride;

    BW::axpy(this->n_elements_active, one(false), other.data(), incx,this->data(),
             incy);

    return *this;
}

//! &Uuml;berladener Operator, bei dem ein Vektor mit einem zweiten Vektor gleichgesetzt wird
//! @param col: Vektor, dessen Werte gleichgesetzt werden
template<typename T, typename BW, typename T_src>
SciPAL::SubVectorView<T, BW, T_src> &
SciPAL::SubVectorView<T, BW, T_src>::operator -= (const Vector<T, BW>& other)
{
    Assert(this->n_elements_active <= other.size(),
           dealii::ExcMessage("Dimension mismatch"));

    int incx = 1; //! col._stride;
    int incy = this->stride;
    T alpha = -1.;
    //!BW::copy(this->n_elements_active, col.data(), incx,  this->data(), incy);
    BW::axpy(this->n_elements_active, alpha, other.data()+ this->__view_begin, incx,this->data(), incy);

    return *this;
}


// @sect4{Operator: *=}
//!
//! &Uuml;berladener Operator, bei dem man einen Vektor mit einem Skalar multipliziert
//! @param alpha: Skalar mit dem der Vektor multipliziert wird
template<typename T, typename BW, typename T_src>
SciPAL::SubVectorView<T, BW, T_src> &
SciPAL::SubVectorView<T, BW, T_src> ::operator *=(const T alpha)
{
    int incx = this->stride;

    BW::scal(this->n_elements_active, alpha, this->data(), incx);

    return *this;
}

template<typename T, typename BW, typename T_src>
template<typename T2>
SciPAL::ColVectorView<T, BW, T_src> & SciPAL::ColVectorView<T, BW, T_src>::operator *= (const T2 scale)
{
    int elem_dist = 1;

    int n = this->n_elements_active;

    Base::BW::scal(n, scale, &(this->data()[0]), elem_dist);

    return *this;
}


// @sect4{Operator: /=}
//!
//! &Uuml;berladener Operator, bei dem ein Vektor durch ein Skalar dividiert wird
//! @param alpha: Skalar mit dem der Vektor geteilt wird
template<typename T, typename BW, typename T_src>
SciPAL::SubVectorView<T, BW, T_src> &
SciPAL::SubVectorView<T, BW, T_src> ::operator /=(const T alpha)
{
#ifdef DEBUG_SUBVECTOR_VIEW
    Assert(alpha,
           dealii::ExcMessage("Div/0"));
#else
    AssertThrow(alpha,
                dealii::ExcMessage("Div/0"));
#endif

    int incx = this->stride;

    BW::scal(this->n_elements_active, (1/alpha), this->data(), incx);

    return *this;
}

template<typename T, typename BW, typename T_src>
SciPAL::SubVectorView<T, BW, T_src> &
SciPAL::SubVectorView<T, BW, T_src> ::operator /=(const std::complex<typename PrecisionTraits<T, BW::arch>::NumberType> alpha)
{
#ifdef DEBUG_SUBVECTOR_VIEW
    Assert(norm(alpha),
           dealii::ExcMessage("Div/0"));
#else
    AssertThrow(norm(alpha),
                dealii::ExcMessage("Div/0"));
#endif

    int incx = this->stride;

    BW::scal(this->n_elements_active, (1./alpha), this->data(), incx);

    return *this;
}


// @sect4{Operator: =}
//!
//! &Uuml;berladener Operator, bei dem ein Vektor mit einem zweiten Vektor gleichgesetzt wird
//! @param col: Vektor, dessen Werte gleichgesetzt werden
template<typename T, typename BW, typename T_src>
template <typename T2>
SciPAL::SubVectorView<T, BW, T_src> &
SciPAL::SubVectorView<T, BW, T_src>::operator = (const std::vector<std::complex<T2> > & dst)
{

    dst.resize(this->size());

    const T * const src_ptr = this->data();
    T * dst_ptr = &dst[0];

    static int inc_src = 1;

    static int inc_dst = 1;

    BW::GetVector(this->size(), src_ptr, inc_src, dst_ptr, inc_dst);


}

template<typename T, typename BW, typename T_src>
T
SciPAL::SubVectorView<T, BW,T_src>::operator () (int k) const
{
    //From cublas vector:
    //    std::vector<T> tmp(this->size());
    //    T * dst_ptr = &tmp[0];

    //    BW::GetVector(1, this->data()+k, 1, dst_ptr, 1);

    // TODO: check this and adapt to getvector:
    T tmp=0;
    tmp = this->data()[this->stride*k];
    return tmp;
}


// @sect4{Funktion: VectorView::print}
//!
//! Ausgabe des Vektors
template<typename T, typename BW, typename T_src>
void
SciPAL::SubVectorView<T, BW, T_src>::print() const
{
    std::vector<T> tmp(this->n_elements_active);
    int inc_src = this->stride;
    int inc_dst = 1;

    T * tmp_ptr = &tmp[0];
    BW::GetVector(this->n_elements_active, this->data(), inc_src, tmp_ptr, inc_dst);

    if (this->_is_col)
        for (int i = 0; i < this->n_elements_active; ++i)
            std::cout << tmp[i] << std::endl;
    else
        for (int i = 0; i < this->n_elements_active; ++i)
            std::cout << tmp[i] << " ";
    std::cout << std::endl;

}
#endif
