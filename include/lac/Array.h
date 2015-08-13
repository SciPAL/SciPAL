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


#ifndef SCIPAL_ARRAY_H
#define SCIPAL_ARRAY_H

//! #include <lac/FullMatrixAccessor.h>

#include <lac/cublas_wrapper.hh>

#include <lac/blas_wrapper.hh>

//! The namespace SciPAL collects all data types needed for encapsulating
//! the blas libraries.
namespace SciPAL {

//! This class encapsulates the memory management provided by the blas wrappers.
//! Inheritance of the Data class of the blas wrapper is non-public to express the
//! "implemented with" relationship, cf. S. Meyers, <i>Effective C++</i>.
template <typename T, typename BW>
class Array : protected BW::template Data<T> {

    typedef typename BW::template Data<T> Base;

public:

    Array();

    //! Construct an array of length @p n.
    //! @param n: Number of elements to allocate.
    Array(size_t n);

    //! Resize an array to length @p n.
    //! @param n: New number of elements.
    void reinit(size_t n);

    //! Read-Write access to the data pointer. Use this function if you know what you are doing.
    T * val();

     //! Read-only access to the data pointer.
    const T * val() const;

    //! Return the number of elements this array currently has.
    int n_elements() const { return __n; }

    //! @param other: Array to copy.
    Array<T, BW> & operator = (const Array<T, BW>& other);

protected:
    //! The number of elements this array currently has. By having it in this class
    //! we do not have to put it again and again into the Data classes provided by the different blas wrappers.
    size_t __n;


    //! To avoid undesired side effects the copy constructor is made private.
    //! If gets accidentally used the compiler will throw an error message.
    Array(const Array<T, BW>& other) {}

    //! The swap function exchanges the attributes of two arrays.
    //! Internally, only the data pointers are exchanged.
    void swap(Array<T, BW>& other)
    {
        std::swap(__n, other.__n);
        this->Base::swap(other.base());
    }

    //! Read-Write access to the data of the base class.
    Base& base() { return *this; }

};

} //! namespace SciPAL END


template <typename T, typename BW>
SciPAL::Array<T, BW>::Array()
    :
      Base(),
      __n(0)
{}

//! @param n : Groesse des Arrays.
template <typename T, typename BW>
SciPAL::Array<T, BW>::Array(size_t n)
    :
      Base(n),
      __n(n)
{}


template <typename T, typename BW>
inline void SciPAL::Array<T, BW>::reinit(size_t n)
{
    if ( this->__n == n) return;

    AssertThrow(n > 0,
                dealii::ExcMessage("Reinitialization requires a positive number of elements.") );

    if (this->__n != n)
    {
        this->Base::resize(n);
        // Alocation errors are handled by the subclass. Thus we can go on.
        __n = n;
    }
}


template <typename T, typename BW>
T * SciPAL::Array<T, BW>::val()
{
#ifdef QT_NO_DEBUG
    AssertThrow(this->data() != 0,
                dealii::ExcMessage("Object not initialized"));
#else
    Assert(this->data() != 0,
           dealii::ExcMessage("Object not initialized") );
#endif

    return this->data();
}



template <typename T, typename BW>
const T * SciPAL::Array<T, BW>::val() const
{
#ifdef QT_NO_DEBUG
    AssertThrow(this->data() != 0,
                dealii::ExcMessage("Object not initialized"));
#else
    Assert(this->data() != 0,
           dealii::ExcMessage("Object not initialized") );
#endif

    return this->data();
}



template <typename T, typename BW>
SciPAL::Array<T, BW> &
SciPAL::Array<T, BW>::operator = (const Array<T, BW>& other)
{

    //Assert(other.__is_allocd,
      //     dealii::ExcMessage("You cannot copy from an uninitialized object!") );

    this->reinit(other.__n);

    // Element-wise copy of arrays.
    int inc_src  = 1;
    int inc_this = 1;

    BW::copy(this->__n, other.val(), inc_src,
             this->data(), inc_this);

    return *this;
}


#endif //! SCIPAL_ARRAY_H
