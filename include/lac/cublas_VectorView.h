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


#ifndef cublas_VectorView_H
#define cublas_VectorView_H

#include <lac/Shape.h>
#include <lac/Expr.h>

//! Vorwaertsdeklaration fuer expression templates
struct vmu_view;



namespace SciPAL {


    // @sect3{Klasse: VectorView}
    //!
    //! Die View-Klassen dienen dazu, Teile von Vektoren oder Matrizen als
    //! selbst&auml;ndige Objekte darzustellen, so dass man beispielsweise eine
    //! Spalte einer Matrix als Quelle in einem Matrix-Vektorprodukt
    //! nutzen kann.
    template<typename T, typename T_src>
    class VectorView
            :
            public SciPAL::Expr<VectorView<T, T_src> >,
            public SciPAL::Shape<T> {

   public:
        typedef typename T_src::blas_wrapper_type BW;

        typedef SciPAL::Shape<T> MyShape;

        friend class Vector<T, BW>;

        friend class Matrix<T, BW>;

        friend class SubMatrixView<T, BW>;

        template<typename, typename> friend class  ColVectorView;

        typedef T value_type;

        typedef VectorView<T, T_src> Type;

        typedef const Type& ConstHandle;

        typedef  SciPAL::ShapeData<T> DevType;


        //! Fuer Optimierungen der Matrix-Vektor-Operationen zwischen
        //! (Teil-)Matrizen und (Teil-)Vektoren muss bekannt sein, was fuer ein
        //! Typ die jeweiligen Objekte sind. Da statische Konstanten zur Compile-Zeit ausgewertet werden
        //! koennen dadurch bedingte Ausdruecke fuer die Komilierung geschrieben werden.
        //! Siehe beispielsweise SubMatrixView::vmult().
        static const bool is_vector_view = true;

        VectorView(T_src & src,
                      int r_begin, int r_end, int c);

        const T * val() const;

        typename PrecisionTraits<T, BW::arch>::NumberType l2_norm() const;

        template<typename VECTOR>
        T dot(const VECTOR & other) const;

        void print() const;

        VectorView<T, T_src> & operator = (const Vector<T, BW>& col);

        template <typename T2>
                VectorView<T, T_src> & operator = (const std::vector<std::complex<T2> > & other);

        int r_begin() const { return __r_begin; }

        int c_begin() const { return __col; }

        int size() const { return __n_el; }

        template<typename T2_src>
        VectorView & operator += (const VectorView <T, T2_src> &other);

        VectorView & operator -= (const VectorView <T, T_src> &other);
        VectorView &operator -= (const Vector<T, BW>& col);

        VectorView & operator *= (const T alpha);
        VectorView & operator /= (const T alpha);
        VectorView & operator /=(const std::complex<typename PrecisionTraits<T, BW::arch>::NumberType> alpha);


        //! Schreibzugriff f&uuml;r Matrix<T>::vmult.
        T * val();

        VectorView & operator = (const VectorView<T, T_src> & other)
        {
//!            this->__src = other.__src;
//!            this->__r_begin = other.__r_begin;
//!            this->__col     = other.__col;
            Assert(this->__n_el    == other.__n_el,
                   dealii::ExcMessage("Cannot copy subarrays of different lengths"));
//!            this->__view_begin = other.__view_begin;
//!            this->_is_col   = other._is_col;
//!            this->_stride   = other._stride;

            //! Elementweise Kopie des arrays.
            int inc_src = other._stride;
            int inc_dst = this->_stride;

            BW::copy(this->__n_el, other.val(), inc_src, this->val(), inc_dst);

            return *this;
        }
private:
        VectorView() {}

        VectorView(const VectorView<T, T_src> & ) {}

            //! deal.II's SmartPointer sind keine eigenst&auml;ndigen Zeiger,
            //! sondern dazu gedacht, auf schon bestehende Objekte zu zeigen
            //! und eine Exception auszul&ouml;sen, wenn das Objekt auf
            //! das sie zeigen vorzeitig zerst&ouml;rt wird.
        dealii::SmartPointer<T_src> __src;

        //! Zeilenindex
        int __r_begin;
        int __col;

    protected:
        int __n_el;

    private:
        int __view_begin;

    protected:
        bool _is_col;
    //! public:
        int _stride;
    };



    // @sect3{Klasse: ColVectorView}
    //!
    //! Ansicht auf einen Teil einer Spalte einer Matrix oder eines Spaltenvektors.
    template<typename T, typename T_src>
    class ColVectorView : public VectorView<T, T_src> {

    public:
        typedef VectorView<T, T_src> Base;

    ColVectorView(T_src & src,
                  int r_begin, int c=0) : Base(src, r_begin, c) { this->_stride = 1; }

#ifdef USE_OLD_ET
    template<typename M> //!, typename Op>
    ColVectorView(const //! X_read_read<M, Op, Vector<T, BW> >
                                    X_read_read<M, vmu_view, SciPAL::ColVectorView<T, T_src> >
                                    & Ax);


    template<typename M> //!, typename Op>
    ColVectorView<T, T_src> & operator = (const
                                    X_read_read<M, vmu_view, ColVectorView<T, T_src> >
                                    & Ax)
    {
        Ax.apply(*this);
        return *this;
    }
#endif

    template<typename T2>
    SciPAL::ColVectorView<T, T_src> & operator *= (const T2 scale)
                                                    {
        int elem_dist = 1;

        int n = this->__n_el;

        Base::BW::scal(n, scale, &(this->val()[0]), elem_dist);

        return *this;
    }

//!    template<typename T2>
//!    SciPAL::ColVectorView<T, T_src> & operator *= (const T2 scale)
//!                                                    {
//!        int elem_dist = 1;

//!        int n = this->__n_el;

//!        Base::BW::scal(n, scale, &(this->val()[0]), elem_dist);

//!        return *this;
//!    }

    ColVectorView<T, T_src> & operator = (const Vector<T, typename Base::BW>& col)
    {
        Assert(this->__n_el == col.size(),
               dealii::ExcMessage("Dimension mismatch"));

        int incx = 1; //! col._stride;
        int incy = this->_stride;
        Base::BW::copy(this->__n_el, col.val(), incx, this->val(), incy);

        return *this;
    }

    //! TO DO:
    /*
    template<typename X>
    ColVectorView<T, T_src> & operator = (const X & x)
    {

        x.apply(*this);

        return *this;
    }
*/
    template<typename T2_src>
    ColVectorView<T, T_src> & operator += (const ColVectorView<T, T2_src> &other)
    {
        Base & self = *this;
        const typename ColVectorView<T, T2_src>::Base & o = other;
        self += o;

        return *this;
    }

    void reset(int r_begin, int c=0)
    {
        this->__r_begin    = r_begin;
        this->__col        = c;
        this->__n_el       = this->__src->n_rows() - r_begin;
        this->__view_begin = c*this->__src->n_rows() + r_begin;
    }


    ColVectorView & operator = (const ColVectorView<T, T_src> & other)
    {
        Base & self = *this;
        const Base & src = other;
        self = src;

        return *this;
    }
    };



    // @sect3{Klasse: RowVectorView}
    //!
    //! Ansicht auf einen Teil einer Zeile einer Matrix oder eines Zeilenvektors.
    template<typename T, typename T_src>
    class RowVectorView : public VectorView<T, T_src> {

    public:
        typedef VectorView<T, T_src> Base;

        RowVectorView(T_src & src,
                      int r_begin, int c) : Base(src, r_begin, c)
        {
            this->__n_el  = src.n_cols() - c;
            this->_stride = src.n_rows(); this->_is_col = false;
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
template<typename T, typename T_src>
SciPAL::VectorView<T, T_src >::VectorView(T_src & src,
                                            int r_begin, int r_end, int c)
    :
      MyShape(reinterpret_cast<T*>(src.data_ptr + (c*src.n_rows() + r_begin)),
                r_end - r_begin, /*rows*/
                1/*cols*/,
                r_end - r_begin/*lda*/),
    __src(&src),
    __r_begin(r_begin),
    __col (c),
    __n_el(src.n_rows() - r_begin),
    __view_begin(c*src.n_rows() + r_begin),
    _is_col(true),
    _stride(1)
{}



// @sect4{Funktion: val}
//!
//! Direkter Zugriff auf den den device-Zeiger.
//! Zurueckgegeben wird der Zeiger auf das erste Element
//! des betrachteten Teils eines Vektors.
template<typename T, typename T_src>
const T *
SciPAL::VectorView<T, T_src>::val() const
{
    //! Indizes in C-Z&auml;hlung!!!
    return &(this->__src->val()[__view_begin]);
}


template<typename T, typename T_src>
T *
SciPAL::VectorView<T, T_src>::val()
{
    //! Indizes in C-Z&auml;hlung!!!
    return &(this->__src->val()[__view_begin]);
}


// @sect4{Funktion: l2_norm}
//!
//! L2-Norm des betrachteten Vektorteils.
template<typename T, typename T_src>
typename PrecisionTraits<T, T_src::blas_wrapper_type::arch>::NumberType
SciPAL::VectorView<T, T_src>::l2_norm() const
{
    typename PrecisionTraits<T, BW::arch>::NumberType result = BW::nrm2(this->__n_el, this->val(), this->_stride);
    return result;
}


// @sect4{Funktion: dot}
//!
//! Skalarprodukt von zwei Teilvektoren.
//! @param other: 2. Vektor
template<typename T, typename T_src>
template<typename VECTOR>
T
SciPAL::VectorView<T, T_src>::dot(const VECTOR & other) const
{

    Assert(this->__n_el == other.size(),
           dealii::ExcMessage("Dimension mismatch"));

    int incx = this->_stride;
    int incy = other._stride;
    T result = BW::dot(this->__n_el,
                       this->val(), incx, other.val(), incy);

    return result;
}


// @sect4{Operator: +=}
//!
//! &Uuml;berladener Operator, mit dem man zwei Vektoren addieren kann
//! 1. Vektor += 2. Vektor
//! @param other: 2. Vektor
template<typename T, typename T_src>
template<typename T2_src>
SciPAL::VectorView<T, T_src> &
SciPAL::VectorView<T, T_src> ::operator +=(const SciPAL::VectorView <T, T2_src> &other)
{
    Assert(this->__n_el == other.size(),
           dealii::ExcMessage("Dimension mismatch"));

    One<T> one;

    int incx = other._stride;
    int incy = this->_stride;

     BW::axpy(this->__n_el, one(), other.val(), incx,this->val(), incy);

    return *this;
}

// @sect4{Operator: -=}
//!
//! &Uuml;berladener Operator, mit dem man zwei Vektoren subtrahieren kann
//! 1. Vektor -= 2. Vektor
//! @param other: 2. Vektor
template<typename T, typename T_src>
SciPAL::VectorView<T, T_src> &
SciPAL::VectorView<T, T_src> ::operator -=(const SciPAL::VectorView <T, T_src> &other)
{
    Assert(this->__n_el == other.size(),
           dealii::ExcMessage("Dimension mismatch"));

    One<T> one;

    int incx = other._stride;
    int incy = this->_stride;

     BW::axpy(this->__n_el, one(false), other.val(), incx,this->val(), incy);

    return *this;
}

//! &Uuml;berladener Operator, bei dem ein Vektor mit einem zweiten Vektor gleichgesetzt wird
//! @param col: Vektor, dessen Werte gleichgesetzt werden
template<typename T, typename T_src>
SciPAL::VectorView<T, T_src> &
SciPAL::VectorView<T, T_src>::operator -= (const Vector<T, BW>& other)
{
    Assert(this->__n_el <= other.size(),
           dealii::ExcMessage("Dimension mismatch"));

    int incx = 1; //! col._stride;
    int incy = this->_stride;
    T alpha = -1.;
    //!BW::copy(this->__n_el, col.val(), incx,  this->val(), incy);
    BW::axpy(this->__n_el, alpha, other.val()+ this->__view_begin, incx,this->val(), incy);

    return *this;
}


// @sect4{Operator: *=}
//!
//! &Uuml;berladener Operator, bei dem man einen Vektor mit einem Skalar multipliziert
//! @param alpha: Skalar mit dem der Vektor multipliziert wird
template<typename T, typename T_src>
SciPAL::VectorView<T, T_src> &
SciPAL::VectorView<T, T_src> ::operator *=(const T alpha)
{
     int incx = this->_stride;

     BW::scal(this->__n_el, alpha, this->val(), incx);

    return *this;
}

// @sect4{Operator: /=}
//!
//! &Uuml;berladener Operator, bei dem ein Vektor durch ein Skalar dividiert wird
//! @param alpha: Skalar mit dem der Vektor geteilt wird
template<typename T, typename T_src>
SciPAL::VectorView<T, T_src> &
SciPAL::VectorView<T, T_src> ::operator /=(const T alpha)
{
#ifdef DEBUG
    Assert(alpha,
           dealii::ExcMessage("Div/0"));
#else
    AssertThrow(alpha,
           dealii::ExcMessage("Div/0"));
#endif

     int incx = this->_stride;

     BW::scal(this->__n_el, (1/alpha), this->val(), incx);

    return *this;
}

template<typename T, typename T_src>
SciPAL::VectorView<T, T_src> &
        SciPAL::VectorView<T, T_src> ::operator /=(const std::complex<typename PrecisionTraits<T, BW::arch>::NumberType> alpha)
{
#ifdef DEBUG
    Assert(norm(alpha),
           dealii::ExcMessage("Div/0"));
#else
    AssertThrow(norm(alpha),
           dealii::ExcMessage("Div/0"));
#endif

     int incx = this->_stride;

     BW::scal(this->__n_el, (1./alpha), this->val(), incx);

    return *this;
}

// @sect4{Operator: =}
//!
//! &Uuml;berladener Operator, bei dem ein Vektor mit einem zweiten Vektor gleichgesetzt wird
//! @param col: Vektor, dessen Werte gleichgesetzt werden
template<typename T, typename T_src>
SciPAL::VectorView<T, T_src> &
SciPAL::VectorView<T, T_src>::operator = (const Vector<T, BW>& col)
{
    Assert(this->__n_el == col.size(),
           dealii::ExcMessage("Dimension mismatch"));

    int incx = 1; //! col._stride;
    int incy = this->_stride;
    BW::copy(this->__n_el, col.val(), incx,  this->val(), incy);

    return *this;
}
template<typename T, typename T_src>
template <typename T2>
SciPAL::VectorView<T, T_src> &
SciPAL::VectorView<T, T_src>::operator = (const std::vector<std::complex<T2> > & dst)
{

dst.resize(this->size());

const T * const src_ptr = this->val();
T * dst_ptr = &dst[0];

static int inc_src = 1;

static int inc_dst = 1;

BW::GetVector(this->size(), src_ptr, inc_src, dst_ptr, inc_dst);


}


// @sect4{Funktion: VectorView::print}
//!
//! Ausgabe des Vektors
template<typename T, typename T_src>
void
SciPAL::VectorView<T, T_src>::print() const
{
    std::vector<T> tmp(this->__n_el);
    int inc_src = this->_stride;
    int inc_dst = 1;

    T * tmp_ptr = &tmp[0];
    BW::GetVector(this->__n_el, this->val(), inc_src, tmp_ptr, inc_dst);

    if (this->_is_col)
        for (int i = 0; i < this->__n_el; ++i)
            std::cout << tmp[i] << std::endl;
    else
        for (int i = 0; i < this->__n_el; ++i)
            std::cout << tmp[i] << " ";
    std::cout << std::endl;

}


//! --------------- EXPRESSION TEMPLATES ----------------------------


// @sect4{struct: vmu}
//!
//! Struktur zum Ausfuehren einer Matrix Vector Multiplikation
//! @param b : Zielvector
//! @param A : Matrix
//! @param x : Vector
struct vmu_view
{
    template<typename T, typename BW, typename T_src>
    static void apply(      SciPAL::ColVectorView<T, T_src> & b,
                      const SciPAL::Matrix<T, BW> & A,
                      const SciPAL::ColVectorView<T, T_src>  & x)
    {
        A.vmult(b,x);
    }


    template<typename T, typename BW, typename T_src>
    static void apply(      SciPAL::ColVectorView<T, T_src> & b,
                      const SciPAL::transpose<SciPAL::Matrix<T, BW> > & A_t,
                      const SciPAL::ColVectorView<T, T_src>  & x)
    {
         A_t.A.Tvmult(b,x);
    }


/*
    template<typename T, typename BW>
    static void apply(      SciPAL::Vector<T, BW> & b,
                      const SciPAL::SubMatrixView<T, BW> & A,
                      const SciPAL::Vector<T, BW> & x)
    {
        b.reinit(A.r_end() - A.r_begin());   A.vmult(b,x);
    }


    template<typename T, typename BW>
    static void apply(      SciPAL::Vector<T, BW> & b,
                      const transpose<SciPAL::Matrix<T, BW> > & A_t,
                      const SciPAL::Vector<T, BW> & x)
    {
        b.reinit(A_t.A.n_cols());    A_t.A.Tvmult(b,x);
    }



    template<typename T, typename BW>
    static void apply(      SciPAL::Vector<T, BW> & b,
                      const transpose<SciPAL::SubMatrixView<T, BW> > & A_t,
                      const SciPAL::Vector<T, BW> & x)
    {
        b.reinit(A_t.A.matrix().n_rows());   A_t.A.Tvmult(b,x);
    }
    */
};

#ifdef USE_OLD_ET
// @sect4{Operator *}
//!
//! Operator fuer Expression Template vmu
//! @param _l : Referenz auf rechte Seite (Matrix)
//! @param _r : Referenz auf linke Seite (Vector)
template<typename L, typename T, //! typename R>
typename BW //!, typename T_src
>
inline
X_read_read<L, vmu_view, SciPAL::ColVectorView<T, SciPAL::Matrix<T, BW>/*T_src*/ >
>
operator * (const L & _l, const SciPAL::ColVectorView<T, //!T_src
            SciPAL::Matrix<T, BW> >
            & _r)
{
    typedef SciPAL::Matrix<T, BW> T_src;
    return X_read_read<L, vmu_view, SciPAL::ColVectorView<T, T_src>
            > (_l,_r);
}
#endif


#ifdef USE_OLD_ET
template<typename T, typename T_src>
template<typename M>
SciPAL::ColVectorView<T, T_src>::ColVectorView(const
                                X_read_read<M, vmu_view, SciPAL::ColVectorView<T, T_src> >
                                & Ax)
{
    Ax.apply(*this);
}
#endif


/*
template<typename T, typename T_src>
template<typename M>
SciPAL::ColVectorView<T, T_src> & SciPAL::ColVectorView<T, T_src>::operator = (const //! X_read_read<M, Op, Vector<T, BW> >
                                X_read_read<M, vmu_view, SciPAL::ColVectorView<T, T_src> >
                                & Ax)
{
    Ax.apply(*this);
    return *this;
}
*/



/*
template<typename T, typename BW>
template<typename M, typename Op>
SciPAL::Vector<T, BW> &
SciPAL::Vector<T, BW>::operator = (const X_read_read<M, Op, Vector<T, BW> > & Ax)
{
    Ax.apply(*this);

    return *this;
}
*/







#endif
