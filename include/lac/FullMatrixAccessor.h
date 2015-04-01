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


#ifndef FullMatrixAccessor_H
#define FullMatrixAccessor_H

#include <deal.II/lac/full_matrix.h>


#include <lac/cublas_Matrix.h>

namespace SciPAL {

    template<typename, typename> class Array;
}

namespace SciPAL {

    template<typename, typename> class Matrix;

    template<typename, typename> class VectorView;
    template<typename, typename> class SubMatrixView;
}

template<typename T> struct ComplexAdapter;


template <>
        struct ComplexAdapter<float2>
{
    typedef std::complex<float> Type;
};

template <>
        struct ComplexAdapter<double2>
{
    typedef std::complex<double> Type;
};


template <typename T>
        struct ComplexAdapter {

    typedef T Type;
};

// @sect3{Klasse: FullMatrixDataAccessor}
//!
//! Fuer die Benutzung von CUBLAS wird der Zugriff auf die Rohdaten
//! einer Matrix benoetigt. Da deal.II's FullMatrix-Klasse dieses aber nicht
//! oeffentlich zur Verfuegung stellt, muessen wir eine abgeleitete Klasse
//! erstellen, die diesen Zugriff nachtraeglich implementiert.
//! Dazu leiten wir von dealii::FullMatrix ab und zwar nicht oeffentlich, um auszudruecken,
//! das diese Klasse mittels dealii::FullMatrix implementiert wird.
//! Da der Sinn dieser Klasse die Speicherung der Matrixeintraege im column-major-Format
//! ist, macht oeffentliche Vererbung keinen Sinn, da das ausdruecken wuerde, dass diese
//! Klasse wie eine dealii::FullMatrix waere. Durch die Aenderung der Speicheranordnung
//! ist diese "IS_A"-Beziehung aber nicht mehr gegeben.
template<typename _T>
class FullMatrixAccessor : protected dealii::FullMatrix<typename ComplexAdapter<_T>::Type>
{
    typedef typename ComplexAdapter<_T>::Type T;

    typedef FullMatrixAccessor<_T> Type;

public:
    //! Fuer die Kompatibilitat mit der dealii::FullMatrix-Klasse muss diese Klasse
    //! ein paar Datentypen definieren. Der Einfachheit halber nehmen wir diese direkt aus
    //! der Basisklasse
    typedef dealii::FullMatrix<T> Base;

    typedef typename Base::value_type value_type;

    // @sect4{Konstruktoren}
    //!
    //! Damit diese Klasse genauso benutzt werden kann wie
    //! deal.II's FullMatrix
    //! muessen alle Konstruktoren ueberladen werden.
    FullMatrixAccessor (const unsigned int n=0);


    FullMatrixAccessor (const unsigned int rows,
                        const unsigned int cols, bool is_column_major=false);

    FullMatrixAccessor (const dealii::FullMatrix<T> & other);

    FullMatrixAccessor (const Type & other);

    FullMatrixAccessor (const dealii::FullMatrix<T> & other,
                        bool transpose_copy);

    FullMatrixAccessor (const unsigned int rows,
                            const unsigned int cols,
                            const T * entries);


    FullMatrixAccessor (const dealii::IdentityMatrix &id, bool is_col_major=false);

    void reinit (size_t nr, size_t nc, bool use_col_major=false) {
        this->Base::reinit(nr,nc);
        this->__is_col_major = use_col_major;
    }

    // @sect4{Funktion: val}
    //!
    //! Das einzig neue ist der direkte Zugriff auf die Rohdaten.
    T * val();

    const T * val() const;

    template<typename BW>
    Type & operator = (const SciPAL::Matrix<_T, BW> & A_d);

    template<typename T2>
    void push_to(dealii::FullMatrix<T2> & dst) const;


    T & operator () (const unsigned int i, const unsigned int j);

    const T & operator () (const unsigned int i, const unsigned int j) const;


    Type & operator += (const Type & A);

    Type & operator -= (const Type & A);

    Type & operator += (const dealii::FullMatrix<T> & A_h);

    Type & operator -= (const dealii::FullMatrix<T> & A_h);

    Type & operator += (const dealii::IdentityMatrix & I_h);

    Type & operator -= (const dealii::IdentityMatrix & I_h);

    T frobenius_norm() const;

    bool is_column_major() const;

    inline unsigned int n_rows() const { return this->Base::n_rows(); }

    inline unsigned int n_cols() const { return this->Base::n_cols(); }

    inline unsigned int n_elements () const { return this->table_size[0]*this->table_size[1]; }

    void print() { this->dealii::FullMatrix<T>::print(std::cout, 11, 4); }

    // A quick hack for using this class in the dealii::Multigrid framework.
    void clear() { this->reinit(0,0); }

private:

    bool __is_col_major;

};











// @sect3{Klasse: FullMatrixDataAccessor}
//!
// @sect4{Konstruktoren}
//!
//! Damit diese Klasse genauso benutzt werden kann wie deal.II's FullMatrix
//! muessen alle Konstruktoren ueberladen werden.
//! Die ersten beiden allokieren lediglich den Speicher fuer
//! eine quadratische oder rechteckige matrix.
//! @param n : Anzahl der Zeilen und Spalten der Matrix.
template<typename _T>
FullMatrixAccessor<_T>::FullMatrixAccessor (const unsigned int n)
        :
        Base(n),
        __is_col_major(false)
{}

//! @param rows : Anzahl der Zeilen der anzulegenden Matrix.
//! @param cols : Anzahl der Spalten der anzulegenden Matrix.
template<typename _T>
FullMatrixAccessor<_T>::FullMatrixAccessor (const unsigned int rows,
                                           const unsigned int cols,bool is_column_major)
        :
        Base(rows, cols),
        __is_col_major(is_column_major)
{}


//! Die Copy-Konstruktoren kopieren eine Matrix entweder so, wie sie ist ...
//! @param other : zu kopierende Matrix.
template<typename _T>
FullMatrixAccessor<_T>::FullMatrixAccessor (const dealii::FullMatrix<T> & other)
        :
        Base(other),
        __is_col_major(false)
{}


//! Die Copy-Konstruktoren kopieren eine Matrix entweder so, wie sie ist ...
//! @param other : zu kopierende Matrix.
template<typename _T>
FullMatrixAccessor<_T>::FullMatrixAccessor (const Type & other)
        :
        Base(other),
        __is_col_major(other.is_column_major())
{}


//! oder sie erstellen eine Kopie der transponierten Matrix.
//! Da deal.II-Matrizen row-major sind, gibt @p transpose_copy den Wert
//! fuer das Attribut @p __is_col_major an.
//! @param other : zu kopierende Matrix.
//! @param transpose_copy : Gibt an, ob die Matrix spalten- anstatt zeilenweise
//! abgespeichert werden soll.
template<typename _T>
FullMatrixAccessor<_T>::FullMatrixAccessor (const dealii::FullMatrix<T> & other,
                                           bool transpose_copy)
    :
    __is_col_major(transpose_copy)
{
    if (!transpose_copy)
        this->copy_from(other);
    else {
         this->reinit( other.n_rows(), other.n_cols(), transpose_copy);

         const unsigned int this_n_rows = other.n_rows();
         const unsigned int this_n_cols = other.n_cols();

         T * entries = this->val();
         for (unsigned int c=0;c<this_n_cols;++c)
            for (unsigned int r=0;r<this_n_rows;++r)
                 entries[c*this_n_rows + r] = other(r,c);
    }
}


//! Ansonsten kann man noch ein lineares array in eine Matrix verwandeln ...
//! @param rows : Anzahl der Zeilen der anzulegenden Matrix.
//! @param cols : Anzahl der Spalten der anzulegenden Matrix.
//! @param entries : Zeiger auf Array in dem die Matrixeintraege
//! zeilenwise hintereinander abgespeichert sind.
template<typename _T>
FullMatrixAccessor<_T>::FullMatrixAccessor (const unsigned int rows,
                                           const unsigned int cols,
                                           const T * entries)
        :
        Base(rows, cols, entries),
        __is_col_major(false)
{}


//! oder eine Kopie der Identitaetsmatrix erstellen.
//! @param id : zu kopierende Identitaetsmatrix.
template<typename _T>
FullMatrixAccessor<_T>::FullMatrixAccessor (const dealii::IdentityMatrix &id, bool is_col_major/*=false*/)
        :
        Base(id),
        __is_col_major(is_col_major)
{}

// @sect4{Operator: =}
//!
//! Die Zuweisungsoperatoren funktionieren mit cublas-Matrizen und -Arrays.
//! @param A_d : zu kopierende Matrix, die im GPU-Speicher liegt.
template<typename _T>
template<typename BW>
FullMatrixAccessor<_T> &
FullMatrixAccessor<_T>::operator = (const SciPAL::Matrix<_T, BW> & A_d)
{
       //! Setze die Groesse entsprechend der der Quelle.

    int nr =  A_d.n_rows();
    int nc =  A_d.n_cols();
    this->reinit(nr, nc);
    this->__is_col_major = true;

    const _T * src = A_d.val();

    _T * dst = reinterpret_cast<_T*>(this->val());

    BW::GetMatrix(nr, nc,
                  src, nr,
                  dst, nr);

    src = 0;
    dst = 0;

    return *this;
}

// @sect4{Function: push_to}
//!
//! In order to copy the content of an object of this type into a dealii::FullMatrix
//! without modifying the dealii::FullMatrix class we provide this member function
//! which internally copies the contents element-wise, so that the entries of
//! the target matrix do not have to be of the same type as the entries of
//! the source matrix.
template <typename _T>
template<typename T2>
void FullMatrixAccessor<_T>::push_to(dealii::FullMatrix<T2> & dst) const
{

    if (this->is_column_major())
    {
        dst.reinit(this->n_rows(), this->n_cols());
        for (unsigned int r = 0; r < this->n_rows(); r++)
            for(unsigned int c = 0; c < this->n_cols(); c++)
                dst(r,c) = (*this)(r,c);
    }
    else {
        const Base & my_self = *this;

        dst = my_self;
    }

}



// @sect4{Operator: ()}
//!
//! Die elementweisen Zugriffe muessen im Falle nichtoeffentlicher
//! Vererbung auch ueberladen werden.
//! @param r : Zeilenindex in C-Zaehlung, d.h. bei 0 beginnend.
//! @param c : Spaltenindex in C-Zaehlung, d.h. bei 0 beginnend.
template <typename _T>
inline
typename FullMatrixAccessor<_T>::T &
FullMatrixAccessor<_T>::operator () (const unsigned int r,
                                    const unsigned int c)
{
  Assert (r < this->table_size[0],
          dealii::ExcIndexRange (r, 0, this->table_size[0]));
  Assert (c < this->table_size[1],
          dealii::ExcIndexRange (c, 0, this->table_size[1]));
  if (this->__is_col_major)
    return this->val()[c*this->table_size[0]+r];
  else
    return this->val()[r*this->table_size[1]+c];
}


//! @param r : Zeilenindex in C-Zaehlung, d.h. bei 0 beginnend.
//! @param c : Spaltenindex in C-Zaehlung, d.h. bei 0 beginnend.
template <typename _T>
inline
const typename FullMatrixAccessor<_T>::T &
FullMatrixAccessor<_T>::operator () (const unsigned int r,
                                    const unsigned int c) const
{
    Assert (r < this->table_size[0],
            dealii::ExcIndexRange (r, 0, this->table_size[0]));
    Assert (c < this->table_size[1],
            dealii::ExcIndexRange (c, 0, this->table_size[1]));
    if (this->__is_col_major)
      return this->val()[c*this->table_size[0]+r];
    else
      return this->val()[r*this->table_size[1]+c];
}


// @sect4{Funktion: val}
//!
//! Die wesentliche Existenzberechtigung dieser Klasse,
//! neben der Moeglichkeit
//! die Matrixeintraege spaltenweise abzupseichern,
//! ist der direkte Zugriff auf die Rohdaten.
template<typename _T>
typename FullMatrixAccessor<_T>::T * FullMatrixAccessor<_T>::val()
{
    //help
    return &this->Base::values[0];
}

template<typename _T>
const typename FullMatrixAccessor<_T>::T * FullMatrixAccessor<_T>::val() const
{
    //help
    return &this->Base::values[0];
}


// @sect4{Funktion: is_column_major()}
//!
//! Diese Funktion gibt Auskunft ueber die interne Datenabspeicherung.
template<typename _T>
bool
FullMatrixAccessor<_T>::is_column_major() const
{
    return this->__is_col_major;
}

// @sect4{Funktion: operator +=}
//!
//! Diese Funktion ueberlaedt den Additionoperator
//! @param A : zu addierende Matrix.
template<typename _T>
FullMatrixAccessor<_T> &
FullMatrixAccessor<_T>::operator += (const Type & A)
{
    Assert(this->n_rows()==A.n_rows(), dealii::ExcMessage("Dimension mismatch"));
    Assert(this->n_cols()==A.n_cols(), dealii::ExcMessage("Dimension mismatch"));
    for(int i=0; i<this->n_rows();i++){
        for(int j = 0; j < this->n_cols() ; j++){
            (*this)(i,j) += A(i,j);
        }
    }
    return *this;
}

// @sect4{Funktion: operator -=}
//!
//! Diese Funktion ueberlaedt den Subtraktionsoperator
//! @param A : zu subtrahiende Matrix.
template<typename _T>
FullMatrixAccessor<_T> &
FullMatrixAccessor<_T>::operator -= (const Type & A)
{
    Assert(this->n_rows()==A.n_rows(), dealii::ExcMessage("Dimension mismatch"));
    Assert(this->n_cols()==A.n_cols(), dealii::ExcMessage("Dimension mismatch"));
    for(int i=0; i<this->n_rows();i++){
        for(int j = 0; j < this->n_cols() ; j++){
            (*this)(i,j) -= A(i,j);
        }
    }
    return *this;
}



// @sect4{Funktion: operator +=}
//!
//! Diese Funktion ueberlaedt den Additionoperator
//! @param A_h : zu addierende Matrix.
template<typename _T>
FullMatrixAccessor<_T> &
FullMatrixAccessor<_T>::operator += (const dealii::FullMatrix<T> & A_h)
{
    Assert(this->n_rows()==A_h.n_rows(), dealii::ExcMessage("Dimension mismatch"));
    Assert(this->n_cols()==A_h.n_cols(), dealii::ExcMessage("Dimension mismatch"));
    for(int i=0; i<this->n_rows();i++){
        for(int j = 0; j < this->n_cols() ; j++){
            (*this)(i,j) += A_h(i,j);
        }
    }
    return *this;
}

// @sect4{Funktion: operator -=}
//!
//! Diese Funktion ueberlaedt den Subtraktionsoperator
//! @param A_h : zu subtrahiende Matrix.
template<typename _T>
FullMatrixAccessor<_T> &
FullMatrixAccessor<_T>::operator -= (const dealii::FullMatrix<T> & A_h)
{
    Assert(this->n_rows()==A_h.n_rows(), dealii::ExcMessage("Dimension mismatch"));
    Assert(this->n_cols()==A_h.n_cols(), dealii::ExcMessage("Dimension mismatch"));
    for(int i=0; i<this->n_rows();i++){
        for(int j = 0; j < this->n_cols() ; j++){
            (*this)(i,j) -= A_h(i,j);
        }
    }
    return *this;
}


// @sect4{Funktion: operator += IdentityMatrix}
//!
//! Diese Funktion ueberlaedt den Additionsoperator
//! @param I_h : zu addierende Identitaetsmatrix.
template<typename _T>
FullMatrixAccessor<_T> &
FullMatrixAccessor<_T>::operator += (const dealii::IdentityMatrix & I_h)
{
	 Assert(this->n_rows()==I_h.n(), dealii::ExcMessage("Dimension mismatch"));
	 Assert(this->n_cols()==I_h.n(), dealii::ExcMessage("Dimension mismatch"));
    for(int i=0; i<this->n_rows();i++){
				(*this)(i,i) += 1;
    }
    return *this;
}

// @sect4{Funktion: operator -= IdentityMatrix}
//!
//! Diese Funktion ueberlaedt den Subtraktionsoperator
//! @param I_h : zu subtrahiende Identitaetsmatrix.
template<typename _T>
FullMatrixAccessor<_T> &
FullMatrixAccessor<_T>::operator -= (const dealii::IdentityMatrix & I_h)
{
	 Assert(this->n_rows()==I_h.n(), dealii::ExcMessage("Dimension mismatch"));
	 Assert(this->n_cols()==I_h.n(), dealii::ExcMessage("Dimension mismatch"));
    for(int i=0; i<this->n_rows();i++){
				(*this)(i,i) -= 1;
    }
    return *this;
}



// @sect4{Funktion: frobenius_norm()}
//!
//! Diese Funktion gibt die frobeniusnorm
template<typename _T>
typename FullMatrixAccessor<_T>::T
FullMatrixAccessor<_T>::frobenius_norm() const
{
  return this->Base::frobenius_norm();
}




#endif //! FULLMATRIXACCESSOR_H
