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


#ifndef cublas_Vector_H
#define cublas_Vector_H


#include <lac/expression_template.h>
#include <lac/expression_templates_host.h>
#include <complex>
#include <float.h>
//! Vorwaertsdeklaration fuer expression templates
struct vmu;

//SciPAL implementation for expression templates

//#include <lac/SciPAL_kernels_wrapper.cu.h>
#include <deal.II/lac/vector.h>

#include <lac/Array.h>
#include <lac/Shape.h>

#include <lac/BlasVectorOperations.h>

#include <base/PrecisionTraits.h>
#include <base/Zero_One_Traits.h>

template<typename> class FullMatrixAccessor;


namespace SciPAL {

    template<typename, typename> class Array;
    template<typename, typename> class Matrix;

    template<typename, typename> class VectorView;
    template<typename, typename> class ColVectorView;

    template<typename, typename> class SubMatrixView;
}

namespace SciPAL {



// @sect3{Klasse: Vector}
//!
//! Mit irgendwas muss man ja rechnen.
template<typename T, typename BW>
class Vector
        :
        public SciPAL::Expr<Vector<T, BW> >,
        protected Array<T, BW>,
        public  dealii::Subscriptor,
        public SciPAL::Shape<T> {


    friend class SciPAL::Matrix<T, BW>;

    friend class SciPAL::SubMatrixView<T, BW>;

    template<typename, typename> friend class VectorView;


    template<typename, typename> friend class ColVectorView;

    typedef SciPAL::Shape<T> MyShape;

public:

    typedef T Number;

    typedef T value_type;

    typedef BW blas_wrapper_type;

    typedef Vector<T, BW> Type;

    typedef const Type& ConstHandle;

    typedef  SciPAL::ShapeData<T> DevType;

        //! Fuer Optimierungen der Matrix-Vektor-Operationen zwischen
        //! (Teil-)Matrizen und (Teil-)Vektoren muss bekannt sein, was fuer ein
        //! Typ die jeweiligen Objekte sind. Da statische Konstanten zur Compile-Zeit ausgewertet werden
        //! koennen dadurch bedingte Ausdruecke fuer die Komilierung geschrieben werden.
    //! Siehe beispielsweise SubMatrixView::vmult().
    static const bool is_vector_view = false;



    Vector();

    Vector(size_t n_elements);

    Vector(size_t n_elements, const Array<T, BW> & raw_data);

    Vector(const FullMatrixAccessor<T> & src,
                 int r_begin, int c);

    Vector(const Matrix<T, BW> & src,
                 int r_begin, int c);

    Vector(const Vector<T, BW> & other)
        :
          dealii::Subscriptor(),
          Array<T,BW>(),
        MyShape()
    {
        *this = other;
    }


    Vector (const std::vector<T> & src) { *this = src; }

        //! Lese- und Schreibzugriff auf das unterliegende array - falls man es doch mal braucht.
        //! Der Zwang explizit diese Funktion aufzurufen sollte Schutz genug sein gegen Missbrauch,
        //! d.h., ist man der Meinung diese Funktion nehmen zu muessen, hat man im Allgemeinen was
        //! falsch gemacht und sich nicht an die Philosophie der Objekt-Orientierung gehalten..
    inline Array<T, BW> & array() { return *this; }

    inline const Array<T, BW> & array() const { return *this; }

    template<typename BW2>
    Vector<T, BW> & operator = (const Vector<T, BW2> & other);

    Vector<T, BW> & operator = (const std::vector<T> & other);

    template<typename T2>
    void push_to(std::vector<T2> & dst) const;

    template<typename T2>
    void push_to(dealii::Vector<T2> & dst) const;

    Vector<T, BW> & operator = (const T value);

    template<typename T_src>
    Vector<T, BW> & operator = (const VectorView<T, T_src > & other);

    template<typename T_src>
    Vector<T, BW> & operator = (const dealii::Vector<T_src> & other);

    template<typename X>
     Vector<T, BW> & operator = (const SciPAL::Expr<X> & e);

    Vector<T, BW> & operator += (const Vector<T, BW> & other);

    template<typename T_src>
    Vector<T, BW> & operator += (const VectorView<T, T_src > & other);

    Vector<T, BW> & operator += (const T scalar);

    Vector<T, BW> & operator -= (const Vector<T, BW> & other);

    template<typename T_src>
    Vector<T, BW> & operator -= (const VectorView<T, T_src > & other);

    template<typename T2>
    Vector<T, BW> & operator *= (const T2 scale);

    Vector<T, BW> & operator *=(Vector<T, BW> &scale);

    T operator * (const Vector<T, BW> & other);

    Vector<T, BW> & operator /= (const T scale);

    Vector<T, BW> & operator /= (Vector<T, BW> &scale);

    template<typename VECTOR>
    T dot(const VECTOR & other) const;

    void sadd (T alpha, const Vector<T, BW> & other);

    void reinit(size_t new_size);

    void reinit(const Vector<T, BW> & other);


    void print() const;

    typename PrecisionTraits<T, BW::arch>::NumberType l2_norm() const;

    typename PrecisionTraits<T, BW::arch>::NumberType sum() const;

    int size()   const { return this->__n; }

    int n_rows() const { return this->__n; }

    int n_cols() const { return this->__n; }

    T operator () (int k) const;

    T& operator [](size_t el);

    void set(int k,const T value);

    void set(const int value);

    void add(int k,const T value);

    // additional functions deal.ii compatibility
    bool all_zero();

    void add(const Vector<T, BW> & other);
    void add(T alpha, const Vector<T, BW> & other);
    void add(T alpha, const Vector<T, BW> & other, T beta, const Vector<T, BW> & other2);

    void sadd(T alpha, T beta, const Vector<T, BW> & other);
    void sadd(T alpha, T beta, const Vector<T, BW> & other, T gamma, const Vector<T, BW> & other2);


    void equ(T alpha, const Vector<T, BW> & other);
    void equ(T alpha, const Vector<T, BW> & other, T beta, const Vector<T, BW> & other2);

protected:
    //! A Vector is just the hull for adding array math operations to
    //! a 1D array. Therefore, the stride is always 1.
   //static const int _stride = 1;

};
}

//namespace SciPAL {
//template<>
//template<>
//Vector<double, cublas> &
//Vector<double, cublas>::operator = (const Vector<double, blas> & other)
//{

//    if(this->size() != other.size())
//    {
//        size_t new_size = other.size();
//        this->reinit(new_size);
//    }
//    //! Elementweise Kopie des arrays.

//    const double * src_ptr = other.array().val();

//    size_t incx = 1;
//    size_t incy = 1;
//    cublas::SetVector(other.size(), src_ptr, incx, this->data(), incy);

//    return *this;
//}
//}

// @sect3{SciPAL::Vector Methoden}


// @sect4{Konstruktor: Vector}
//!
//! Standardkontruktor fuer Vector

template<typename T, typename BW>
SciPAL::Vector<T, BW>::Vector()
    :
      Array<T, BW>(),
      dealii::Subscriptor(),
      MyShape(this->data(), 0, 0, 0 /*TO DO: leading dim*/)
{}

// @sect4{Konstruktor: Vector(n)}
//!
//! Konstruktor fuer n-elementigen Vector
//! @param n_elements : Groesse des Vectors

template<typename T, typename BW>
SciPAL::Vector<T, BW>::Vector(size_t n_elements)
    :
    Array<T, BW>(n_elements),
    dealii::Subscriptor(),
    MyShape(this->array().val(),
              n_elements, /*n_rows*/
              1, /*n_cols*/
              n_elements, /*leading dim*/
              1) /*stride*/
{}

// @sect4{Konstruktor: Vector(n,raw_data)}
//!
//! Konstruktor fuer n-elementigen Vector, dessen Inhalt aus einem Array-Objekt kommen
//! @param n_elements : Groesse des Vectors
//! @param raw_data : Pointer auf ein Array Objekt

template<typename T, typename BW>
SciPAL::Vector<T, BW>::Vector(size_t n_elements,
                              const Array<T, BW> & raw_data)
    :
    Array<T, BW>(n_elements),
    dealii::Subscriptor(),
    MyShape(this->array().val(),
            n_elements, 1, n_elements/*TO DO: leading dim*/)
{
    *this = raw_data;
}


// @sect4{Konstruktor: Vector(FullMatrixAccessor,r_begin,c)}
//!
//! Konstruiert einen Spalten-Vector aus einer Matrix, indem das Startelement angegeben wird
//! @param src : Pointer auf einen FullMatrixAccessor
//! @param r_begin : Zeilenindex des der Matrix
//! @param c : Spaltenindex der Matrix

template<typename T, typename BW>
SciPAL::Vector<T, BW>::Vector(const FullMatrixAccessor<T> & src,
                              int r_begin, int c)
    :
    Array<T, BW>(src.n_rows() - r_begin),
    dealii::Subscriptor(),
    MyShape(this->array().val(), this->__n, 1, this->__n /*TODO: leading dim*/)
{
    int n_el = src.n_rows() - r_begin;

    int src_begin = c*src.n_rows() + r_begin;
    const T * src_val = src.val();

    int inc_src = 1;
    if (!src.is_column_major())
    {
        inc_src = src.n_cols();
        src_begin = r_begin*inc_src + c;
    }

    const T * src_ptr = &(src_val[src_begin]);

    BW::SetVector(n_el, src_ptr, inc_src, this->val(), 1);
}


// @sect4{Konstruktor: Vector(Matrix,r_begin,c)}
//!
//! Konstruiert einen Spalten-Vector aus einer Matrix, indem das Startelement angegeben wird
//! @param src : Pointer auf eine Matrix
//! @param r_begin : Zeilenindex des der Matrix
//! @param c : Spaltenindex der Matrix

template<typename T, typename BW>
SciPAL::Vector<T, BW>::Vector(const Matrix<T, BW> & src,
                              int r_begin, int c)
    :
    Array<T, BW>(src.n_rows() - r_begin),
    dealii::Subscriptor(),
    MyShape(this->array().val(), this->__n, 1, this->__n)
{
    int n_el = src.n_rows() - r_begin;

        //! Elementweise Kopie des arrays.
    int inc_src = 1;
    int inc_this = 1;

    const T * col = &(src.val()[c*src.n_rows() + r_begin]);

    //! cublasScopy
    BW::copy(this->__n, col, inc_src, this->val(), inc_this);
}



// @sect4{Funktion: reinit}
//!
//! Erzeugt Nullvektor der Lange @p new_size auf der GPU.
//! @param newsize : Groesse des Vectors

template<typename T, typename BW>
void SciPAL::Vector<T, BW>::reinit(size_t new_size)
{
    this->Array<T, BW>::reinit(new_size);

    size_t n_cols = 1;

    this->MyShape::reinit(this->array().val(), new_size, n_cols,
                          new_size /*leading_dim*/, 1/*stride*/);
}



template<typename T, typename BW>
void
SciPAL::Vector<T, BW>::reinit(const Vector<T, BW> & other)
{
    this->reinit( other.size() );
}



// @sect4{Funktion: l2_norm}
//!
//! Berechnet die L2-Norm eines Vectors

template<typename T, typename BW>
typename PrecisionTraits<T, BW::arch>::NumberType
SciPAL::Vector<T, BW>::l2_norm() const
{
    return BW::nrm2(this->__n, &(this->val()[0]), 1);
}

// @sect4{Funktion: sum}
//!
//! Berechnet die Summe über alle Vektorelemente

template<typename T, typename BW>
typename PrecisionTraits<T, BW::arch>::NumberType
SciPAL::Vector<T, BW>::sum() const
{
    return BW::asum(this->__n, &(this->val()[0]), 1);
}


// @sect4{Funktion: Vector::print}
//!
//! Kopiert die Daten vom device zurueck zum host
//! und gibt sie auf der Standardausgabe aus.

template <typename T>
std::ostream &operator << (std::ostream &ostr,SciPAL::CudaComplex<T> &a ){

    ostr << "(" << a.real() << "," << a.imag() << ")";
    return ostr;
}


template<typename T, typename BW>
void
SciPAL::Vector<T, BW>::print() const
{

    std::vector<T> tmp(this->__n);
    T * dst_ptr = &tmp[0];

    const T * src_ptr = this->val();

    BW::GetVector(this->__n, src_ptr, 1, dst_ptr, 1);

    for (size_t i = 0; i < this->__n; ++i)
         std::cout  << tmp[i] << std::endl;
}

// @sect4{Funktion: Vector::()}
//!
//! Elementzugriff auf Vektor - lesend
//!
//! @param k : Index des Vectors

template<typename T, typename BW>
T
SciPAL::Vector<T, BW>::operator () (int k) const
{
    std::vector<T> tmp(this->__n);
    T * dst_ptr = &tmp[0];

    BW::GetVector(1, this->val()+k, 1, dst_ptr, 1);
    return tmp[0];
}

// @sect4{Funktion: Vector::[]}
//!
//! Elementzugriff auf Vektor - lesend
//!
//! @param el : Index des Vectors

template<typename T, typename BW>
T&
SciPAL::Vector<T, BW>::operator [] (size_t el)
{
    return this->data_ptr[el];
}

// @sect4{Funktion: Vector::set}
//!
//! Elementzugriff auf Vektor - schreibend
//! @param k : Index des Vectors
//! @param value : zu schreibender Wert
template<typename T, typename BW>
void
SciPAL::Vector<T, BW>::set(int k, const T value)
{
    int inc_src = 1;
    int inc_dst = 1;

    BW::SetVector(1, &value, inc_src, &(this->val()[k]), inc_dst);
}

// @sect4{Function: Vector::set}
//!
//! Sets all elements of a vector to @p value.
//! @param value : zu schreibender Wert
template<typename T, typename BW>
void
SciPAL::Vector<T, BW>::set(const int value)
{

    BW::Memset(this->n_elements(), this->data_ptr, value);
}


// @sect4{Funktion: Vector::add}
//!
//! Elementzugriff auf Vektor - schreibend
//! @param k : Index des Vectors
//! @param value : zu addierender Wert
template<typename T, typename BW>
void
SciPAL::Vector<T, BW>::add(int k,const T value)
{
    Assert((k >= 0) && (k < this->size()),
           dealii::ExcMessage("Index out of range") );


    Vector<T, BW> tmp_d(1);
    tmp_d.set(0, value);


    BW::axpy(1/*other.size()*/, 1, tmp_d.val()/*other.val()*/, 1,
             &(this->val()[k]), 1);
}



// @sect4{Operator: =}
//!
//! Elementweise Kopie eines Vektors. Das Ziel wird bei Bedarf
//! in der Groesse der Quelle angepassst.
//! @param other : rechte Seite des = ist ein Vector
template<typename T, typename BW>
template<typename BW2>
SciPAL::Vector<T, BW> &
SciPAL::Vector<T, BW>::operator = (const Vector<T, BW2> & other)
{

    // element-wise copy of array.
    int inc_src  = 1;
    int inc_this = 1;
    //! same blas type no problem
    if(typeid(BW) == typeid(BW2) )
        BW::copy(this->n_elements(), other.val(), inc_src,
                 this->val(), inc_this);

    //! copy from cublas matrix to blas matrix -> GetMatrix
    //! TODO: what is with asyn copy?
    if(typeid(BW) == typeid(blas) && typeid(BW2) == typeid(cublas) )
    {
        cublas::GetMatrix(other.n_rows(), other.n_cols(), other.array().val(),
                          other.leading_dim, this->array().val(), this->leading_dim);
    }

    //! copy from cublas matrix to blas matrix -> SetMatrix
    //! TODO: what is with asyn copy?
    if(typeid(BW) == typeid(cublas) && typeid(BW2) == typeid(blas) )
    {
        cublas::SetMatrix(other.n_rows(), other.n_cols(),
                          other.array().val(),
                          other.leading_dim,
                          this->array().val(),
                          this->leading_dim);
    }

    std::cout<<__FUNCTION__<<std::endl;
    return *this;
}



// @sect4{Operator: =}
//!
//! Elementweise Kopie eines Vektors. Das Ziel wird bei Bedarf
//! in der Groesse der Quelle angepassst.
//! @param other : rechte Seite des = ist ein std::Vector
template<typename T, typename BW>
SciPAL::Vector<T, BW> &
SciPAL::Vector<T, BW>::operator = (const std::vector<T> & other)
{
    size_t new_size = other.size();
    this->reinit(new_size);

    const T * src_ptr = &other[0];

    size_t incx = 1;
    size_t incy = 1;
    BW::SetVector(other.size(), src_ptr, incx, this->data(), incy);

    return *this;
}


template <typename T, typename BW>
template <typename X>
SciPAL::Vector<T, BW> & SciPAL::Vector<T, BW>::operator =
(const SciPAL::Expr<X> & e)
{
#ifdef DEBUG
    std::cout << "line :" << __LINE__ << ", Vector<T,BW>  " << __FUNCTION__<< "\n"  << std::endl;
#endif

    SciPAL::LAOOperations::apply(*this, ~e);
    return *this;
}

template<typename T, typename BW>
template<typename T2>
void
SciPAL::Vector<T, BW>::push_to(std::vector<T2> & dst) const
{
    dst.resize(this->size());

    const T * const src_ptr = this->val();

    T * dst_ptr = reinterpret_cast<T*>(&dst[0]);

    static int inc_src = 1;

    static int inc_dst = 1;

    BW::GetVector(this->size(), src_ptr, inc_src, dst_ptr, inc_dst);
}

template<typename T, typename BW>
template<typename T2>
void
SciPAL::Vector<T, BW>::push_to(dealii::Vector<T2> & dst) const
{
    dst.reinit(this->size());

    const T * const src_ptr = this->val();

    T * dst_ptr = reinterpret_cast<T*>(dst.begin());

    static int inc_src = 1;

    static int inc_dst = 1;

    BW::GetVector(dst.size(), src_ptr, inc_src, dst_ptr, inc_dst);
}



// @sect4{Operator: =}
//!
//! Elementweise Kopie eines Vektors. Das Ziel wird bei Bedarf
//! in der Groesse der Quelle angepassst.
//! @param other : rechte Seite des = ist ein VectorView
//! Beim Kopieren wird das Koordinatensystem des Views uebernommen,
//! d.h. Beginn und Ende des Views bestimmen den Bereich in den kopiert wird.
//! Siehe Berechnung der Householder-Vektoren im QR-Beispiel.
template<typename T, typename BW>
template<typename T_src>
SciPAL::Vector<T, BW> &
SciPAL::Vector<T, BW>::operator = (const VectorView<T, T_src > & other)
{
    Assert(this->size() >= other.size(),
           dealii::ExcDimensionMismatch(this->size(), other.size()) );

    int incx = other._stride;
    int incy = 1;
    BW::copy(other.size(), other.val(), incx,
             &(this->val()[(other._is_col ?
                              other.r_begin() : other.c_begin())
                             ]), incy);

    return *this;
}


template<typename T, typename BW>
template<typename T_src>
SciPAL::Vector<T, BW> &
SciPAL::Vector<T, BW>::operator = (const dealii::Vector<T_src> & other)
{
    // int size = ;
    size_t new_size = 0;

    if(this->size() >= other.size() )
        new_size = this->size();
    else
        new_size = other.size();

    this->reinit(new_size);

//    Assert(this->size() >= other.size(),
//           dealii::ExcDimensionMismatch(this->size(), other.size()) );

    int incx = 1;
    int incy = 1;
//    BW::copy(other.size(), &(*other.begin()), incx,
//             this->data(), incy);

    const T * tmp_ptr = &(*other.begin());
    BW::SetVector(size(), tmp_ptr, incx, this->data(), incy);

    return *this;
}


// @sect4{Operator: =}
//!
template<typename T, typename BW>
SciPAL::Vector<T, BW> &
SciPAL::Vector<T, BW>::operator = (const T value)
{
    Vector<T, BW> tmp(1);
    tmp.set(0, value);
    int incx = 0; //! copy the same element all the time
    int incy = 1;
    BW::copy(tmp.size(), tmp.val(), incx, this->val(), incy);

    return *this;
}

// @sect4{Operator: +=}
//!
//! Addiert zwei Vectoren
//! @param other : rechte Seite des += ist ein Vector

template<typename T, typename BW>
SciPAL::Vector<T, BW> &
SciPAL::Vector<T, BW>::operator += (const Vector<T, BW> & other)
{
    One<T> one;
    BW::axpy(this->__n, one(), other.val(), 1, this->val(), 1);

    return *this;
}

// @sect4{Operator: +=}
//!
//! Addiert VectorView zu Vector
//! @param other : rechte Seite des += ist ein VectorView

template<typename T, typename BW>
template<typename T_src>
SciPAL::Vector<T, BW> &
SciPAL::Vector<T, BW>::operator += (const VectorView<T, T_src > & other)
{

    One<T> one;
    Assert(this->size() >= other.size(),
           dealii::ExcMessage("Dimension mismatch") );

    BW::axpy(other.size(), one(), other.val(), 1,
             &(this->val()[other.r_begin()]), 1);

    return *this;
}

// @sect4{Operator: +=}
//!
//! Addiert VectorView zu Vector
//! @param other : rechte Seite des += ist ein VectorView

template<typename T, typename BW>
SciPAL::Vector<T, BW> &
SciPAL::Vector<T, BW>::operator += (const T scalar)
{
    std::vector<T> scalar_vec(this->size(), scalar);
    Vector d_scalar;
    d_scalar = scalar_vec;
    One<T> one;

    BW::axpy(this->__n, one(), d_scalar.array().val(), 1,
             this->val(), 1);

    return *this;
}


// @sect4{Operator: *=}
//!
//! Skalare Vectormultiplikation
//! @param scale : Faktor

template<typename T, typename BW>
template<typename T2>
SciPAL::Vector<T, BW> &
SciPAL::Vector<T, BW>::operator *= (const T2 scale)
{
    int elem_dist = 1;

    BW::scal(this->__n, scale, &(this->val()[0]), elem_dist);

    return *this;
}


// @sect4{Operator: -=}
//!
//! Subtrahiert einen Vector von einem anderen
//! @param other : Vector, der subtrahiert werden soll


template<typename T, typename BW>
SciPAL::Vector<T, BW> &
SciPAL::Vector<T, BW>::operator -= (const Vector<T, BW> & other)
{
    One<T> one;
    bool plus = true;

    BW::axpy(this->__n, one(!plus), other.val(), 1, this->val(), 1 );

    return *this;
}


// @sect4{Operator: -=}
//!
//! Subtrahiert einen VectorView von einem Vector
//! @param other : VectorView, der subtrahiert werden soll

template<typename T, typename BW>
template<typename T_src>
SciPAL::Vector<T, BW> &
SciPAL::Vector<T, BW>::operator -= (const VectorView<T, T_src > & other)
{
    Assert(this->size() >= other.size(),
           dealii::ExcMessage("Dimension mismatch") );

    BW::axpy(other.size(), -1, other.val(), 1,
             &(this->val()[other.r_begin()]), 1);

    return *this;
}



 // @sect4{Operator: *}
 //!
//! Skalarprodukt zweier Vectoren
//! @param other : zweiter Vektor

template<typename T, typename BW>
T
SciPAL::Vector<T, BW>::operator * (const Vector<T, BW> & other)
{
     Assert(this->size() == other.size(),
           dealii::ExcMessage("Dimension mismatch") );
    return this->dot(other);
}



// @sect4{Operator: /=}
//!
//! Skalare Vectordivision
//! @param scale : Faktor

template<typename T, typename BW>
SciPAL::Vector<T, BW> &
SciPAL::Vector<T, BW>::operator /= (const T scale)
{
    int elem_dist = 1;

    Assert(scale,
           dealii::ExcMessage("Division by Zero") );

    BW::scal(this->__n, 1/scale, this->val(), elem_dist);

    return *this;
}


// @sect4{Funktion: dot}
//!
//! Berechnet das Skalarprodukt von zwei Vektoren.
//! @param other : zweiter Vector

template<typename T, typename BW>
template<typename VECTOR>
T
SciPAL::Vector<T, BW>::dot(const VECTOR & other) const
{

    //!AssertThrow(false,
     //!           dealii::ExcMessage("Strange implementation. "
     //!                              "Please check source code of the body of this function. SCK"));

     Assert(this->size() == other.size(),
     dealii::ExcMessage("Dimension mismatch"));

    int incx = 1; //! this->_stride;
    int incy = 1; //!other._stride;
    T result = BW::dot(this->__n,
                       this->val(), incx, other.val(), incy);

    return result;
}


// @sect4{Funktion: sadd}
//!
//! \f{eqnarray*} this & = & \alpha \cdot this \\ this & = & this + other \f}
//! @param alpha : Faktor
//! @param other : zweiter Vector

template<typename T, typename BW>
void
SciPAL::Vector<T, BW>::sadd (T alpha, const Vector<T, BW> & other)
{
    int incx = 1;
    BW::scal(this->__n, alpha, this->val(), incx);

    //! cublasSaxpy
    int incy = 1;
    BW::axpy(this->__n, 1., other.dev_ptr, 1, this->val(), 1);
}



// @sect4{struct: vmu}
//!
//! Struktur zum Ausfuehren einer Matrix Vector Multiplikation
//! @param b : Zielvector
//! @param A : Matrix
//! @param x : Vector
struct vmu
{
    template<typename T, typename BW>
    static void apply(      SciPAL::Vector<T, BW> & b,
                      const SciPAL::Matrix<T, BW> & A,
                      const SciPAL::Vector<T, BW> & x)
    {
        b.reinit(A.n_rows());   A.vmult(b,x);
    }

    template<typename T, typename BW>
    static void apply(      SciPAL::Vector<T, BW> & b,
                      const SciPAL::SubMatrixView<T, BW> & A,
                      const SciPAL::Vector<T, BW> & x)
    {
        b.reinit(A.r_end() - A.r_begin());   A.vmult(b,x);
    }


    template<typename T, typename BW>
    static void apply(      SciPAL::Vector<T, BW> & b,
                            const SciPAL::transpose<SciPAL::Matrix<T, BW> > & A_t,
                      const SciPAL::Vector<T, BW> & x)
    {
        b.reinit(A_t.A.n_cols());    A_t.A.Tvmult(b,x);
    }



    template<typename T, typename BW>
    static void apply(      SciPAL::Vector<T, BW> & b,
                      const SciPAL::transpose<SciPAL::SubMatrixView<T, BW> > & A_t,
                      const SciPAL::Vector<T, BW> & x)
    {
        b.reinit(A_t.A.matrix().n_rows());   A_t.A.Tvmult(b,x);
    }


    template<typename T, typename BW, typename T_src>
    static void apply(      SciPAL::Vector<T, BW> & b,
                      const SciPAL::Matrix<T, BW> & A,
                      const SciPAL::ColVectorView<T, T_src>  & x)
    {
        b.reinit(A.n_rows());   A.vmult(b,x);
    }
};

#ifdef USE_OLD_ET
// @sect4{Operator *}
//!
//! Operator fuer Expression Template vmu
//! @param _l : Referenz auf rechte Seite (Matrix)
//! @param _r : Referenz auf linke Seite (Vector)
template<typename L, typename T, typename BW>
inline
X_read_read<L, vmu, SciPAL::Vector<T, BW> >
operator * (const L & _l, const SciPAL::Vector<T, BW> & _r)
{
    return X_read_read<L, vmu, SciPAL::Vector<T, BW> > (_l,_r);
}

template<typename L, typename T, //!typename R> //!
typename BW, typename T_src>
inline
X_read_read<L, vmu, SciPAL::ColVectorView<T, T_src> >
operator * (const L & _l, const //! R
             SciPAL::ColVectorView<T, T_src>
            & _r)
{
    return X_read_read<L, vmu, //! R
            SciPAL::ColVectorView<T, T_src>
            > (_l,_r);
}
#endif

#ifdef USE_OLD_ET
template<typename T, typename BW>
template<typename M, typename Op>
SciPAL::Vector<T, BW>::Vector(const X_read_read<M, Op, Vector<T, BW> > & Ax)
{
    Ax.apply(*this);
}


template<typename T, typename BW>
template<typename M, typename Op
          , typename T_src
>
SciPAL::Vector<T, BW>::Vector(const //! X_read_read<M, Op, Vector<T, BW> >
                                X_read_read<M, Op, SciPAL::ColVectorView<T, T_src> >
                                & Ax)
{
    Ax.apply(*this);
}




template<typename T, typename BW>
template<typename M, typename Op>
SciPAL::Vector<T, BW> &
SciPAL::Vector<T, BW>::operator = (const X_read_read<M, Op, Vector<T, BW> > & Ax)
{
    Ax.apply(*this);

    return *this;
}


template<typename T, typename BW>
template<typename M, //! typename Op,
typename T_src>
SciPAL::Vector<T, BW> &
SciPAL::Vector<T, BW>::operator=(const
                                X_read_read<M, vmu, SciPAL::ColVectorView<T, T_src> >
                                & Ax)
{
    Ax.apply(*this);

     return *this;
}
#endif


// additional functions
template<typename T, typename BW>
bool
SciPAL::Vector<T, BW>::all_zero ()
{
    //return ( (x<10.0*DBL_MIN) && (x>-10.0*DBL_MIN) ? true : false)
    T nrm = BW::nrm2(this->__n, &(this->val()[0]), 1);

    return ( (nrm < 10.0*DBL_MIN) && (nrm >-10.0*DBL_MIN)? true:false);
}

template<typename T, typename BW>
void
SciPAL::Vector<T, BW>::add(const Vector<T, BW> & other)
{
    One<T> one;
    BW::axpy(this->__n, one(), other.val(), 1, this->val(), 1);
}



template<typename T, typename BW>
void
SciPAL::Vector<T, BW>::add(T alpha, const Vector<T, BW> & other)
{
    BW::axpy(this->__n, alpha, other.val(), 1, this->val(), 1);

}

template<typename T, typename BW>
void
SciPAL::Vector<T, BW>::add(T alpha, const Vector<T, BW> & other, T beta, const Vector<T, BW> & other2)
{
    BW::axpy(this->__n, alpha, other.val(), 1, this->val(), 1);
    BW::axpy(this->__n, beta, other2.val(), 1, this->val(), 1);

}

template<typename T, typename BW>
void
SciPAL::Vector<T, BW>::sadd(T alpha, T beta, const Vector<T, BW> & other)
{

    BW::scal(this->__n, alpha, &(this->val()[0]), 1);

    BW::axpy(this->__n, beta, other.val(), 1, this->val(), 1);
}

template<typename T, typename BW>
void
SciPAL::Vector<T, BW>::sadd(T alpha, T beta, const Vector<T, BW> & other, T gamma, const Vector<T, BW> & other2)
{

    BW::scal(this->__n, alpha, &(this->val()[0]), 1);

    BW::axpy(this->__n, beta, other.val(), 1, this->val(), 1);
    BW::axpy(this->__n, gamma, other2.val(), 1, this->val(), 1);
}

template<typename T, typename BW>
void
SciPAL::Vector<T, BW>::equ(T alpha, const Vector<T, BW> & other)
{

    // this = \alpha * other
    BW::scal(this->__n, 0.0, &(this->val()[0]), 1);

    BW::axpy(this->__n, alpha, other.val(), 1, this->val(), 1);

}

template<typename T, typename BW>
void
SciPAL::Vector<T, BW>::equ(T alpha, const Vector<T, BW> & other, T beta, const Vector<T, BW> & other2)
{

    // this = \alpha * other + \beta * other2

    dealii::Vector<T> tmp(other.size());
    tmp = 0.;

    int incx = 1;
    int incy = 1;

    const T * tmp_ptr = &(*tmp.begin());
    BW::SetVector(size(), tmp_ptr, incx, this->data(), incy);

    BW::axpy(this->__n, alpha, other.val(), incx, this->val(), incy);
    BW::axpy(this->__n, beta, other2.val(), incx, this->val(), incy);

}
#endif //! cublas_Vector_H

